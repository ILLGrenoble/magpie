/**
 * hdf5 file handling
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18-oct-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/readdata_8cpp-example.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/writedata_8cpp-example.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1CommonFG.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1DataSet.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1DataSpace.html
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_H5__
#define __TLIBS2_H5__

#include <H5Cpp.h>

#include <type_traits>
#include <memory>
#include <string>
#include <vector>


namespace tl2 {


/**
 * get corresponding h5 primitive data type
 */
template<class T>
/*constexpr*/ H5::PredType get_h5_type()
{
	if constexpr(std::is_integral_v<T> && std::is_signed_v<T>)
	{
		if constexpr(sizeof(T) == 1)
			return H5::PredType::NATIVE_INT8;
		else if constexpr(sizeof(T) == 2)
			return H5::PredType::NATIVE_INT16;
		else if constexpr(sizeof(T) == 4)
			return H5::PredType::NATIVE_INT32;
		else if constexpr(sizeof(T) == 8)
			return H5::PredType::NATIVE_INT64;
		else
			return H5::PredType::NATIVE_INT;
	}
	else if constexpr(std::is_integral_v<T> && !std::is_signed_v<T>)
	{
		if constexpr(sizeof(T) == 1)
			return H5::PredType::NATIVE_UINT8;
		else if constexpr(sizeof(T) == 2)
			return H5::PredType::NATIVE_UINT16;
		else if constexpr(sizeof(T) == 4)
			return H5::PredType::NATIVE_UINT32;
		else if constexpr(sizeof(T) == 8)
			return H5::PredType::NATIVE_UINT64;
		else
			return H5::PredType::NATIVE_UINT;
	}
	else if constexpr(std::is_floating_point_v<T>)
	{
		if constexpr(sizeof(T) == 4)
			return H5::PredType::NATIVE_FLOAT;
		else if constexpr(sizeof(T) == 8)
			return H5::PredType::NATIVE_DOUBLE;
		else if constexpr(sizeof(T) == 16)
			return H5::PredType::NATIVE_LDOUBLE;
		else
			return H5::PredType::NATIVE_DOUBLE;
	}

	// default fallback type
	return H5::PredType::NATIVE_INT;
}


// --------------------------------------------------------------------------------


/**
 * get the entries of a given directory from an hdf5 file
 */
template<class t_strvec = std::vector<std::string>>
bool get_h5_entries(H5::H5File& file, const std::string& path, t_strvec& entries)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		// H5G_iterate_t callback, see: /usr/include/hdf5/serial/H5Dpublic.h
		H5G_iterate_t collect_elems = [](hid_t /*group*/, const char* name, void* _vec) -> herr_t
		{
			((t_strvec*)(_vec))->push_back(name);
			return 0;
		};

		file.iterateElems(path, 0, collect_elems, &entries);
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * get an entry's attribute
 */
template<class t_str = std::string>
t_str get_h5_attr(H5::H5File& file, const t_str& path, const t_str& attr_name, bool only_group = false)
{
	if(!file.nameExists(path))
		return "";

	H5::H5Object *obj = nullptr;
	H5::DataSet set;
	H5::Group grp;

	// get entry's type
	H5G_stat_t stat;
	file.getObjinfo(path, true, stat);
	if(stat.type == H5G_GROUP)
	{
		grp = file.openGroup(path);
		obj = &grp;
	}
	else if(!only_group && stat.type == H5G_DATASET)
	{
		set = file.openDataSet(path);
		obj = &set;
	}

	if(!obj || !obj->attrExists(attr_name))
		return "";

	// get entry's attribute
	t_str attr_val;
	try
	{
		H5::Attribute attr = obj->openAttribute(attr_name);

		hsize_t len = attr.getStorageSize();
		attr.read(attr.getStrType(), attr_val);
		attr.close();
	}
	catch(const H5::AttributeIException&)
	{}

	obj->close();
	return attr_val;
}


// --------------------------------------------------------------------------------


/**
 * get a scalar value from an hdf5 file
 */
template<class T>
bool get_h5_scalar(H5::H5File& file, const std::string& path, T& val)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		hsize_t rank = dspace.getSimpleExtentNdims();
		if(rank != 1 && rank != 0)
			return false;
		if(rank == 1)
		{
			hsize_t dim = 0;
			dspace.getSimpleExtentDims(&dim, 0);
			if(dim != 1)
				return false;
		}

		H5::PredType ty = get_h5_type<T>();
		dset.read(&val, ty);

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a scalar value to an hdf5 file
 */
template<class T>
bool set_h5_scalar(H5::H5File& file, const std::string& path, T val)
{
	try
	{
		H5::PredType ty = get_h5_type<T>();

		hsize_t dim = 0;
		H5::DataSpace dspace(0, &dim);
		H5::DataSet dset = file.createDataSet(path, ty, dspace);

		dset.write(&val, ty);

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


// --------------------------------------------------------------------------------


/**
 * get a string value from an hdf5 file
 */
template<class t_str = std::string>
bool get_h5_string(H5::H5File& file, const std::string& path, t_str& val)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		hsize_t rank = dspace.getSimpleExtentNdims();
		if(rank != 1)
			return false;

		hsize_t dim = 0;
		dspace.getSimpleExtentDims(&dim, 0);
		if(dim != 1)
			return false;

		dset.read(val, dset.getStrType());

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a string value to an hdf5 file
 */
template<class t_str = std::string>
bool set_h5_string(H5::H5File& file, const std::string& path, const t_str& val)
{
	try
	{
		hsize_t dim = 0;
		H5::DataSpace dspace(0, &dim);

		H5::DataSet dset = file.createDataSet(
			path, H5::StrType(H5::PredType::C_S1, H5T_VARIABLE), dspace);

		dset.write(val, dset.getStrType());

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


// --------------------------------------------------------------------------------


/**
 * get a vector of values from an hdf5 file
 */
template<class T, template<class...> class t_vec = std::vector>
bool get_h5_vector(H5::H5File& file, const std::string& path, t_vec<T>& vals)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		hsize_t rank = dspace.getSimpleExtentNdims();
		if(rank != 1)
			return false;

		hsize_t dim = 0;
		dspace.getSimpleExtentDims(&dim, 0);

		H5::PredType ty = get_h5_type<T>();

		vals.resize(dim);
		dset.read(vals.data(), ty);

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a vector of values to an hdf5 file
 */
template<class T, template<class...> class t_vec = std::vector>
bool set_h5_vector(H5::H5File& file, const std::string& path, const t_vec<T>& vals)
{
	try
	{
		H5::PredType ty = get_h5_type<T>();

		hsize_t dim = vals.size();
		H5::DataSpace dspace(1, &dim);
		H5::DataSet dset = file.createDataSet(path, ty, dspace);

		dset.write(vals.data(), ty);

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


// --------------------------------------------------------------------------------


/**
 * get a vector of string value from an hdf5 file
 */
template<class t_str_vec = std::vector<std::string>>
bool get_h5_string_vector(H5::H5File& file, const std::string& path, t_str_vec& vals)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		hsize_t rank = dspace.getSimpleExtentNdims();
		if(rank != 1)
			return false;

		hsize_t dim = 0;
		dspace.getSimpleExtentDims(&dim, 0);
		vals.resize(dim);

		const char **chs = new const char*[dim];
		dset.read(chs, H5::StrType(H5::PredType::C_S1, H5T_VARIABLE));
		for(hsize_t i=0; i<dim; ++i)
			vals[i] = chs[i];
		for(hsize_t i=0; i<dim; ++i)
			delete[] chs[i];
		delete[] chs;

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a vector of string values to an hdf5 file
 */
template<class t_str_vec = std::vector<std::string>>
bool set_h5_string_vector(H5::H5File& file, const std::string& path, const t_str_vec& val)
{
	try
	{
		hsize_t dim = val.size();
		H5::DataSpace dspace(1, &dim);

		H5::DataSet dset = file.createDataSet(
			path, H5::StrType(H5::PredType::C_S1, H5T_VARIABLE), dspace);

		const char **chs = new const char*[dim];
		for(std::size_t i=0; i<dim; ++i)
			chs[i] = val[i].c_str();
		dset.write(chs, dset.getStrType());
		delete[] chs;

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


// --------------------------------------------------------------------------------


/**
 * get a matrix of value from an hdf5 file
 */
template<class T, template<class...> class t_vec = std::vector>
bool get_h5_matrix(H5::H5File& file, const std::string& path, t_vec<t_vec<T>>& vals)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		hsize_t rank = dspace.getSimpleExtentNdims();
		if(rank != 2)
			return false;

		hsize_t dims[2] = { 0, 0 };
		dspace.getSimpleExtentDims(dims, 0);

		H5::PredType ty = get_h5_type<T>();

		T* buf = new T[dims[0]*dims[1]];
		dset.read(buf, ty);

		vals.resize(dims[0]);
		for(hsize_t row=0; row<dims[0]; ++row)
		{
			vals[row].resize(dims[1]);
			for(hsize_t col=0; col<dims[1]; ++col)
				vals[row][col] = buf[row*dims[1] + col];
		}

		delete[] buf;

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a matrix of values to an hdf5 file
 */
template<class T, template<class...> class t_vec = std::vector>
bool set_h5_matrix(H5::H5File& file, const std::string& path, const t_vec<t_vec<T>>& vals)
{
	try
	{
		H5::PredType ty = get_h5_type<T>();

		hsize_t dim1 = vals.size();
		hsize_t dim2 = dim1 ? vals[0].size() : 0;
		hsize_t dims[] = { dim1, dim2 };

		auto buf = std::make_unique<T[]>(dims[0] * dims[1]);
		for(hsize_t row=0; row<dims[0]; ++row)
			for(hsize_t col=0; col<dims[1]; ++col)
				buf[row*dims[1] + col] = vals[row][col];

		H5::DataSpace dspace(2, dims);
		H5::DataSet dset = file.createDataSet(path, ty, dspace);

		dset.write(buf.get(), ty);

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


// --------------------------------------------------------------------------------


/**
 * get a multidimensional array of value from an hdf5 file
 */
template<class T, template<class...> class t_vec = std::vector>
bool get_h5_multidim(H5::H5File& file, const std::string& path,
	hsize_t& rank, t_vec<hsize_t>& dims, t_vec<T>& vals)
{
	try
	{
		if(!file.nameExists(path))
			return false;

		H5::PredType ty = get_h5_type<T>();

		H5::DataSet dset = file.openDataSet(path);
		H5::DataSpace dspace = dset.getSpace();

		rank = dspace.getSimpleExtentNdims();

		dims.resize(rank);
		dspace.getSimpleExtentDims(dims.data(), 0);

		hsize_t linear_dim = rank ? dims[0] : 0;
		for(hsize_t i=1; i<rank; ++i)
			linear_dim *= dims[i];
		vals.resize(linear_dim);

		dset.read(vals.data(), ty);

		dspace.close();
		dset.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


/**
 * write a multidimensional array of values to an hdf5 file
 */
template<class T>
bool set_h5_multidim(H5::H5File& file, const std::string& path,
	hsize_t rank, const hsize_t* dims, const T* vals)
{
	try
	{
		H5::PredType ty = get_h5_type<T>();

		H5::DataSpace dspace(rank, dims);
		H5::DataSet dset = file.createDataSet(path, ty, dspace);

		dset.write(vals, ty);

		dset.close();
		dspace.close();
		return true;
	}
	catch(const H5::Exception& ex)
	{
		return false;
	}
}


}

#endif
