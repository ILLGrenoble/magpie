/**
 * calculates the rotation axis of a matrix
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8-jun-20, 18-sep-26
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// clang++ -std=c++20 -I /opt/homebrew/Cellar/boost/*/include -I.. -o rotation_axis rotation_axis.cpp 

#include <iostream>
#include <boost/math/quaternion.hpp>

#include "libs/maths.h"
using namespace tl2_ops;

using t_real = double;
using t_vec = tl2::vec<t_real>;
using t_mat = tl2::mat<t_real>;
using t_quat = boost::math::quaternion<t_real>;

static constexpr t_real eps = 1e-8;


int main(int argc, char** argv)
{
	if(argc < 3*3 + 1)
	{
		std::cerr << "Please enter a 3x3 rotation matrix, example:\n\t"
			<< argv[0] << " 1 0 0  0 \"cos(5/180*pi)\" \"sin(5/180*pi)\"  0 \"-sin(5/180*pi)\" \"cos(5/180*pi)\""
			<< std::endl;
		return -1;
	}

	// get matrix
	t_mat mat = tl2::unit<t_mat>(3, 3);
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			mat(i, j) = tl2::str_to_var_parse<t_real>(std::string(argv[1 + i*3 + j]));

	// get rotateion axis and angle
	t_quat quat = tl2::rot3_to_quat<t_mat, t_quat>(mat);
	auto [axis, angle] = tl2::rotation_axis<t_quat, t_vec>(quat);

	axis /= tl2::norm(axis);

	std::cout << "mat = " << mat << std::endl;
	std::cout << "quat = " << quat << std::endl;
	std::cout << "axis = " << axis << std::endl;
	std::cout << "angle = " << angle/tl2::pi<t_real>*180. << " deg" << std::endl;

	return 0;
}
