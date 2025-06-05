/**
 * extracts a camera detector image from an h5 file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 21/may/2025
 * @license GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// g++ -std=c++20 -O2 -I.. -I /usr/include/hdf5/serial -o h5camdet h5camdet.cpp -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_cpp -lhdf5 -lpng

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include <filesystem>
namespace fs = std::filesystem;

#include <boost/algorithm/string.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;

#include "libs/h5file.h"
#include "libs/file.h"


template<class t_img = gil::gray16_image_t, class t_real = double>
bool extract_image(const std::string& file, bool add_angle_infos = true,
	t_real scale = 1., t_real offs = 0.)
{
	using t_pix_val = std::uint32_t;
	using t_img_val = typename t_img::value_type;
	//const t_pix_val max_img_val = std::numeric_limits<t_img_val>::max();
	const t_pix_val max_img_val = (1 << (sizeof(t_img_val)*8)) - 1;

	// load the scan file
	H5::H5File h5file(file, H5F_ACC_RDONLY);

	std::vector<t_pix_val> data;
	hsize_t rank = 0;
	std::vector<hsize_t> dims;
	if(bool ok = tl2::get_h5_multidim(h5file, "/entry0/data/CameraDetector_data", rank, dims, data); !ok)
	{
		std::cerr << "Cannot read file \"" << file << "\"." << std::endl;
		return false;
	}

	if(rank < 2)
	{
		std::cerr << "Invalid rank in file \"" << file << "\"." << std::endl;
		return false;
	}

	double omega{}, tilt1{}, tilt2{};
	std::string instr{};
	tl2::get_h5_scalar(h5file, "entry0/instrument/omega/value", omega);
	tl2::get_h5_scalar(h5file, "entry0/instrument/tilt1/value", tilt1);
	tl2::get_h5_scalar(h5file, "entry0/instrument/tilt2/value", tilt2);
	tl2::get_h5_string(h5file, "entry0/instrument/name", instr);

	std::string user{}, local{}, title{}, prop{}, time{};
	double dur{};
	tl2::get_h5_string(h5file, "entry0/user/name", user);
	tl2::get_h5_string(h5file, "entry0/user/namelocalcontact", local);
	tl2::get_h5_string(h5file, "entry0/title", title);
	tl2::get_h5_string(h5file, "entry0/user/proposal", prop);
	tl2::get_h5_string(h5file, "entry0/start_time", time);
	tl2::get_h5_scalar(h5file, "entry0/duration", dur);


	// save it as png image
	t_img png(dims[0], dims[1]);
	auto png_view = gil::view(png);

	t_pix_val min_counts = std::numeric_limits<t_pix_val>::max();
	t_pix_val max_counts = 0;

	for(hsize_t y = 0; y < dims[1]; ++y)
	{
		auto png_row = png_view.row_begin(y);

		for(hsize_t x = 0; x < dims[0]; ++x)
		{
			t_pix_val pixel = data[x * dims[1] + y];
			min_counts = std::min(pixel, min_counts);
			max_counts = std::max(pixel, max_counts);

			// scale
			t_real scaled = (static_cast<t_real>(pixel) + offs) * scale;
			if(scaled < 0.)
				scaled = 0.;

			// clamp value to data range
			t_img_val val = static_cast<t_img_val>(
				std::min<t_real>(max_img_val, scaled));
			*(png_row + x) = static_cast<t_img_val>(val);
		}
	}

	h5file.close();


	std::string angle_infos;

	if(add_angle_infos)
	{
		// add angle infos and replace commas with underscores
		std::ostringstream ostr_angle_infos;
		ostr_angle_infos.precision(4);
		ostr_angle_infos << "_omega" << omega << "_tilta" << tilt1 << "_tiltb" << tilt2;
		angle_infos = ostr_angle_infos.str();
		boost::replace_all(angle_infos, ".", "_");
	}

	fs::path file_png = file;
	file_png.replace_filename(file_png.stem().string() + angle_infos);
	file_png.replace_extension("png");
	gil::write_view(file_png.filename(), png_view, gil::png_tag());


	// print some infos
	std::cout << file << " -> " << file_png.filename().string() << "\n";
	//std::cout << "\trank = " << rank << "\n";
	std::cout << "\timg_dims = ";
	for(hsize_t dim : dims)
		std::cout << dim << ", ";
	std::cout << "\n";
	std::cout << "\tω = " << omega << ", tilts = " << tilt1 << ", " << tilt2 << ",\n";
	std::cout << "\tuser = " << user << ", local_contact = " << local << ", instrument = " << instr << ",\n";
	std::cout << "\tproposal = " << prop << ", title = " << title << ",\n";
	std::cout << "\ttime = " << time << ", duration = " << dur << " s,\n";
	std::cout << "\tcount_range = [ " << min_counts << ", " << max_counts << " ].\n";
	std::cout.flush();

	return true;
}


template<class t_img = gil::gray16_image_t, class t_real = double>
void extract_images(const std::string& dir, bool add_angle_infos = true,
	t_real scale = 1., t_real offs = 0.)
{
	for(const std::string& file : tl2::get_all_files(dir))
	{
		std::string ext = tl2::get_fileext(file, true);
		if(ext != "nxs" && ext != "h5")
			continue;

		extract_image<t_img, t_real>(file, add_angle_infos, scale, offs);
		std::cout << std::endl;
	}
}


int main(int argc, char **argv)
{
	// types
	using t_img = gil::gray16_image_t;
	//using t_img = gil::gray8_image_t;
	using t_real = double;

	// options
	bool add_angle_infos = true;
	t_real scale = 1.;
	t_real offs = 0.;
	//t_real scale = 255. / 1000. * 2.5;
	//t_real offs = -75.;

	if(argc < 2)
	{
		std::cerr << "Please give a h5 file or a directory." << std::endl;
		return -1;
	}

	try
	{
		H5::Exception::dontPrint();

		std::string file(argv[1]);
		if(tl2::dir_exists(file))
			extract_images<t_img, t_real>(file, add_angle_infos, scale, offs);
		else if(tl2::file_exists(file))
			extract_image<t_img, t_real>(file, add_angle_infos, scale, offs);
		else
		{
			std::cerr << "File or directory \"" << file
				<< "\" does not exist." << std::endl;
			return -1;
		}
	}
	catch(const H5::Exception& ex)
	{
		std::cerr << ex.getDetailMsg() << std::endl;
	}

	return 0;
}
