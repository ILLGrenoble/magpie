/**
 * rotation matrix unit tests
 * @author Tobias Weber <tweber@ill.fr>
 * @date march-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * @note Forked on 8-March-2022 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).

 *
 * References:
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -I .. -o rotation rotation.cpp
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "mathlibs" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#define BOOST_TEST_MODULE test_rotation

#include <iostream>
#include <tuple>

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace ty = boost::typeindex;

#include "libs/maths.h"


BOOST_AUTO_TEST_CASE_TEMPLATE(test_rotation,
	t_real, decltype(std::tuple<float, double/*, long double*/>{}))
{
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	using namespace tl2_ops;
	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 0.5);

	std::cout << "Testing with " << ty::type_id_with_cvr<t_real>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	// test 2d case
	const t_vec vecFrom2d = tl2::create<t_vec>({1, 1});
	const t_vec vecTo2d = tl2::create<t_vec>({1, 0});
	t_mat mat2d = tl2::rotation<t_mat, t_vec, t_real>(vecFrom2d, vecTo2d, nullptr, eps, false);
	t_vec vec2d = mat2d * vecFrom2d / tl2::norm<t_vec>(vecFrom2d) * tl2::norm<t_vec>(vecTo2d);
	std::cout << "2d rotation:\n";
	tl2::niceprint(std::cout, mat2d, 1e-4, 4);
	std::cout << std::endl;
	BOOST_TEST((tl2::equals<t_vec>(vecTo2d, vec2d, eps)));
	BOOST_TEST((tl2::equals<t_real>(tl2::det<t_mat>(mat2d), 1., eps)));

	// test 3d case
	const t_vec vecFrom3d = tl2::create<t_vec>({-1, -1, 0});
	const t_vec vecTo3d = tl2::create<t_vec>({0, 0, 1});
	t_mat mat3d = tl2::rotation<t_mat, t_vec, t_real>(vecFrom3d, vecTo3d, nullptr, eps, false);
	t_vec vec3d = mat3d * vecFrom3d / tl2::norm<t_vec>(vecFrom3d) * tl2::norm<t_vec>(vecTo3d);
	std::cout << "3d rotation:\n";
	tl2::niceprint(std::cout, mat3d, 1e-4, 4);
	std::cout << std::endl;
	//std::cout << "test rotation:\n";
	//tl2::niceprint(std::cout, mat3d*tl2::create<t_vec>({1., 2., 3}), 1e-4, 4);
	//std::cout << std::endl;
	BOOST_TEST((tl2::equals<t_vec>(vecTo3d, vec3d, eps)));
	BOOST_TEST((tl2::equals<t_real>(tl2::det<t_mat>(mat3d), 1., eps)));

	// test 4d case
	const t_vec vecFrom4d = tl2::create<t_vec>({1, 0, 1, 0});
	const t_vec vecTo4d = tl2::create<t_vec>({1, -1, 0, 0});
	t_mat mat4d = tl2::rotation<t_mat, t_vec, t_real>(vecFrom4d, vecTo4d, nullptr, eps, false);
	t_vec vec4d = mat4d * vecFrom4d / tl2::norm<t_vec>(vecFrom4d) * tl2::norm<t_vec>(vecTo4d);
	std::cout << "4d rotation:\n";
	tl2::niceprint(std::cout, mat4d, 1e-4, 4);
	std::cout << std::endl;
	BOOST_TEST((tl2::equals<t_vec>(vecTo4d, vec4d, eps)));
	BOOST_TEST((tl2::equals<t_real>(tl2::det<t_mat>(mat4d), 1., eps)));
}
