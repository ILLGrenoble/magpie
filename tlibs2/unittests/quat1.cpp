/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8-jun-20
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#define BOOST_TEST_MODULE Quat1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>
#include <boost/math/quaternion.hpp>

#include "libs/maths.h"
using namespace tl2_ops;


using t_types = std::tuple<long double, double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_quat = boost::math::quaternion<t_real>;

	t_real angle = tl2::pi<t_real>/t_real{4};
	auto vec = tl2::create<t_vec>({1, 2, 3});
	auto quat = tl2::rotation_quat<t_vec, t_quat>(vec, angle);

	auto [vec2, angle2] = tl2::rotation_axis<t_quat, t_vec>(quat);

	vec /= tl2::norm(vec);
	vec2 /= tl2::norm(vec2);

	std::cout << "q = " << quat << std::endl;
	std::cout << "axis1 = " << vec << ", angle1 = " << angle << std::endl;
	std::cout << "axis2 = " << vec2 << ", angle2 = " << angle2 << std::endl;


	BOOST_TEST(tl2::equals(vec, vec2, 1e-5));
	BOOST_TEST(tl2::equals<t_real>(angle, angle2, 1e-5));
}
