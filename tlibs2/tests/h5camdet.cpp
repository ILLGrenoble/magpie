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

// g++ -std=c++20 -I.. -I /usr/include/hdf5/serial -o h5camdet h5camdet.cpp -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_cpp -lhdf5 -lpng

#include <iostream>
#include <sstream>
#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include <boost/algorithm/string.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;

#include "libs/h5file.h"
#include "libs/file.h"


bool extract_image(const std::string& file)
{
	// load the scan file
	H5::H5File h5file(file, H5F_ACC_RDONLY);

	std::vector<std::uint32_t> data;
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
	gil::gray16_image_t png(dims[0], dims[1]);
	auto png_view = gil::view(png);

	for(hsize_t y = 0; y < dims[1]; ++y)
	{
		auto png_row = png_view.row_begin(y);

		for(hsize_t x = 0; x < dims[0]; ++x)
			*(png_row + x) = data[x * dims[1] + y];
	}

	h5file.close();

	// add angle infos and replace commas with underscores
	std::ostringstream ostr_angle_infos;
	ostr_angle_infos.precision(4);
	ostr_angle_infos << "_omega" << omega << "_tilta" << tilt1 << "_tiltb" << tilt2;
	std::string angle_infos = ostr_angle_infos.str();
	boost::replace_all(angle_infos, ".", "_");

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
	std::cout << "\ttime = " << time << ", duration = " << dur << " s.\n";
	std::cout.flush();

	return true;
}


void extract_images(const std::string& dir)
{
	for(const std::string& file : tl2::get_all_files(dir))
	{
		std::string ext = tl2::get_fileext(file, true);
		if(ext != "nxs" && ext != "h5")
			continue;

		extract_image(file);
		std::cout << std::endl;
	}
}


int main(int argc, char **argv)
{
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
			extract_images(file);
		else if(tl2::file_exists(file))
			extract_image(file);
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
