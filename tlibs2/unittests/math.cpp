/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#define BOOST_TEST_MODULE Math0
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_math, t_real, t_types)
{
	const t_real eps = 1e-5;

	// coordinate trafos
	{
		t_real x = tl2::get_rand<t_real>(-10., 10.);
		t_real y = tl2::get_rand<t_real>(-10., 10.);
		t_real z = tl2::get_rand<t_real>(-10., 10.);

		auto [ rho, phi, theta ] = tl2::cart_to_sph(x, y, z);
		auto [ x2, y2, z2 ] = tl2::sph_to_cart(rho, phi, theta);

		BOOST_TEST(tl2::equals<t_real>(x, x2, eps));
		BOOST_TEST(tl2::equals<t_real>(y, y2, eps));
		BOOST_TEST(tl2::equals<t_real>(z, z2, eps));

		auto [ rho_c, phi_c, z_c ] = tl2::sph_to_cyl(rho, phi, theta);
		auto [ rho2, phi2, theta2 ] = tl2::cyl_to_sph(rho_c, phi_c, z_c);

		BOOST_TEST(tl2::equals<t_real>(rho, rho2, eps));
		BOOST_TEST(tl2::equals<t_real>(phi, phi2, eps));
		BOOST_TEST(tl2::equals<t_real>(theta, theta2, eps));

		auto [ x3, y3, z3 ] = tl2::cyl_to_cart(rho_c, phi_c, z_c);

		BOOST_TEST(tl2::equals<t_real>(x, x3, eps));
		BOOST_TEST(tl2::equals<t_real>(y, y3, eps));
		BOOST_TEST(tl2::equals<t_real>(z, z3, eps));

		auto [ rho_c2, phi_c2, z_c2 ] = tl2::cart_to_cyl(x3, y3, z3);

		BOOST_TEST(tl2::equals<t_real>(rho_c, rho_c2, eps));
		BOOST_TEST(tl2::equals<t_real>(phi_c, phi_c2, eps));
		BOOST_TEST(tl2::equals<t_real>(z_c, z_c2, eps));

		t_real u = tl2::get_rand<t_real>(0., 1.);
		t_real v = tl2::get_rand<t_real>(0., 1.);

		auto [ phi_uv, theta_uv ] = tl2::uv_to_sph(u, v);
		auto [ u2, v2 ] = tl2::sph_to_uv(phi_uv, theta_uv);

		BOOST_TEST(tl2::equals<t_real>(u, u2, eps));
		BOOST_TEST(tl2::equals<t_real>(v, v2, eps));
	}
}
