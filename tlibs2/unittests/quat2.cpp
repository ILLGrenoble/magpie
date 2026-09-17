/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8-jun-20
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

#define BOOST_TEST_MODULE Quat2
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>
#include <boost/math/quaternion.hpp>

#include "libs/maths.h"
using namespace tl2_ops;



using t_types = std::tuple</*long double,*/ double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_rot, t_real, t_types)
{
	#include "defs.h"
	using t_quat = boost::math::quaternion<t_real>;
	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-5 : 1e-8;

	t_vec axis = tl2::create<t_vec>({1, 2, 3});
	t_real angle = tl2::pi<t_real>/t_real{4};
	t_mat mat = tl2::rotation<t_mat, t_vec>(axis, angle, false);

	t_quat quat1 = tl2::rot3_to_quat<t_mat, t_quat>(mat);
	t_quat quat2 = tl2::rotation_quat<t_vec, t_quat>(axis, angle);

	auto [axis1, angle1] = tl2::rotation_axis<t_quat, t_vec>(quat1);
	auto [axis2, angle2] = tl2::rotation_axis<t_quat, t_vec>(quat2);

	t_mat mat2 = tl2::quat_to_rot3<t_quat, t_mat>(quat2);

	axis /= tl2::norm(axis);
	axis1 /= tl2::norm(axis1);
	axis2 /= tl2::norm(axis2);

	std::cout << "\nUsing " << sizeof(t_real)*8 << "-bit floats." << std::endl;
	std::cout << "q1 = " << quat1 << std::endl;
	std::cout << "q2 = " << quat2 << std::endl;
	std::cout << "axis  = " << axis  << ", angle  = " << angle << std::endl;
	std::cout << "axis1 = " << axis1 << ", angle1 = " << angle1 << std::endl;
	std::cout << "axis2 = " << axis2 << ", angle2 = " << angle2 << std::endl;

	BOOST_TEST(tl2::equals(axis, axis1, eps));
	BOOST_TEST(tl2::equals<t_real>(angle, angle1, eps));
	BOOST_TEST(tl2::equals(axis, axis2, eps));
	BOOST_TEST(tl2::equals<t_real>(angle, angle2, eps));
	BOOST_TEST(tl2::equals(quat1, quat2, eps));
	BOOST_TEST(tl2::equals(mat, mat2, eps));
}



BOOST_AUTO_TEST_CASE_TEMPLATE(test_rot_vec, t_real, t_types)
{
	#include "defs.h"
	using t_quat = boost::math::quaternion<t_real>;
	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-5 : 1e-8;

	t_vec vec1 = tl2::create<t_vec>({1, 0, 0});
	t_vec vec2 = tl2::create<t_vec>({1, 1, 0});

	t_mat mat = tl2::rotation<t_mat, t_vec>(vec1, vec2);
	t_quat quat = tl2::rotation_quat<t_vec, t_quat>(vec1, vec2);

	t_quat quat2 = tl2::rot3_to_quat<t_mat, t_quat>(mat);

	t_vec vec2b = tl2::quat_vec_prod<t_quat, t_vec>(quat, vec1);
	t_vec vec2c = tl2::prod_mv<t_mat, t_vec>(mat, vec1);

	vec2 /= tl2::norm(vec2);

	std::cout << "\nUsing " << sizeof(t_real)*8 << "-bit floats." << std::endl;
	std::cout << "q1 = " << quat << std::endl;
	std::cout << "q2 = " << quat2 << std::endl;
	std::cout << "vec2b = " << vec2b << std::endl;
	std::cout << "vec2c = " << vec2c << std::endl;

	BOOST_TEST(tl2::equals(quat, quat2, eps));
	BOOST_TEST(tl2::equals(vec2, vec2b, eps));
	BOOST_TEST(tl2::equals(vec2, vec2c, eps));
}
