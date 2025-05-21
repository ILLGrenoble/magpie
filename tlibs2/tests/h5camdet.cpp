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

// g++ -std=c++20 -I.. -I /usr/include/hdf5/serial -o h5camdet h5camdet.cpp -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_cpp -lhdf5

#include <iostream>
#include <string>
#include "libs/h5file.h"


bool extract_image(const char* file)
{
	H5::H5File h5file(file, H5F_ACC_RDONLY);

	std::vector<std::uint32_t> data;
	hsize_t rank = 0;
	std::vector<hsize_t> dims;
	if(bool ok = tl2::get_h5_multidim(h5file, "/entry0/data/CameraDetector_data", rank, dims, data); !ok)
	{
		std::cerr << "Cannot read file \"" << file << "\"." << std::endl;
		return false;
	}

	std::cout << "Rank: " << rank << std::endl;
	std::cout << "Dimensions: ";
	for(hsize_t dim : dims)
		std::cout << dim << " ";
	std::cout << std::endl;

	if(rank < 2)
	{
		std::cerr << "Invalid rank in file \"" << file << "\"." << std::endl;
		return false;
	}

	for(hsize_t y = 0; y < dims[1]; ++y)
	{
		for(hsize_t x = 0; x < dims[0]; ++x)
			std::cout << data[y * dims[0] + x] << " ";
		std::cout << std::endl;
	}

	h5file.close();
	return true;
}


int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cerr << "Please give a h5 file." << std::endl;
		return -1;
	}

	try
	{
		H5::Exception::dontPrint();
		extract_image(argv[1]);
	}
	catch(const H5::Exception& ex)
	{
		std::cerr << ex.getDetailMsg() << std::endl;
	}

	return 0;
}
