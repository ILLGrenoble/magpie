/**
 * crystal matrix text
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -o cryst cryst.cpp
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o cryst cryst.cpp -llapacke
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

#define BOOST_TEST_MODULE Xtal Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;


#include <iostream>
#include <vector>

#include "libs/maths.h"
#include "libs/phys.h"
using namespace tl2_ops;


#ifdef USE_LAPACK
	using t_types = std::tuple<double, float>;
#else
	using t_types = std::tuple<long double, double, float>;
#endif


BOOST_AUTO_TEST_CASE_TEMPLATE(test_xtal, t_real, t_types)
{
	std::cout << "Test 1" << std::endl;

	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	//using t_cplx = std::complex<t_real>;
	//using t_vec_cplx = std::vector<t_cplx>;
	//using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	auto A = tl2::A_matrix<t_mat, t_real>(
		4.56, 4.56, 4.56,
		90. / 180. * tl2::pi<t_real>,
		90. / 180. * tl2::pi<t_real>,
		90. / 180. * tl2::pi<t_real>);
	auto B = tl2::B_matrix<t_mat, t_real>(
		4.56, 4.56, 4.56,
		90. / 180. * tl2::pi<t_real>,
		90. / 180. * tl2::pi<t_real>,
		90. / 180. * tl2::pi<t_real>);
	auto [B2, ok] = tl2::inv<t_mat>(A);
	B2 = 2.*tl2::pi<t_real> * tl2::trans<t_mat>(B2);
	auto G = tl2::metric<t_mat>(B);

	t_vec vec1 = tl2::create<t_vec>({1, 1, 0});
	t_vec vec2 = tl2::create<t_vec>({1, -1, 0});
	t_vec vec3 = tl2::cross<t_mat, t_vec>(B, vec1, vec2);
	t_mat UB = tl2::UB_matrix<t_mat, t_vec>(B, vec1, vec2, vec3);

	std::cout << "A  = " << A << std::endl;
	std::cout << "B  = " << B << std::endl;
	std::cout << "B2 = " << B2 << std::endl;
	std::cout << "G  = " << G << std::endl;
	std::cout << "UB = " << UB << std::endl;
	//std::cout << tl2::levi<t_mat>(B, {0,1,2}) << std::endl;
	//std::cout << tl2::levi<t_mat>(B, {0,2,1}) << std::endl;
	//std::cout << tl2::levi<t_mat>(B, {1,2,1}) << std::endl;

	t_real ki = 1.5;
	t_real kf = 1.4;
	t_vec Q = tl2::create<t_vec>({1, -1, 0});
	auto [anglesok, a3, a4, dist] = tl2::calc_tas_a3a4<t_mat, t_vec>(B, ki, kf, Q, vec1, vec3);
	std::cout << "a3 = " << a3 / tl2::pi<t_real> * 180.
		<< ", a4 = " << a4 / tl2::pi<t_real> * 180.
		<< std::endl;
	std::cout << "distance of Q to scattering plane: " << dist << std::endl;


	// calculate back to Q
	t_real Qlen = tl2::calc_tas_Q_len<t_real>(ki, kf, a4);
	std::optional<t_vec> Qhkl = tl2::calc_tas_hkl<t_mat, t_vec>(B, ki, kf, Qlen, a3, vec1, vec3);
	std::cout << "Q = " << *Qhkl << std::endl;


	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals(B, B2, std::numeric_limits<t_real>::epsilon()*1e2));
	BOOST_TEST(tl2::equals<t_real>(dist, 0, std::numeric_limits<t_real>::epsilon()*1e2));
	BOOST_TEST(tl2::equals<t_real>(a3/tl2::pi<t_real>*180., 44.359, 1e-3));
	BOOST_TEST(tl2::equals<t_real>(a4/tl2::pi<t_real>*180., 84.359, 1e-3));
	BOOST_TEST((Qhkl.operator bool()));
	BOOST_TEST(tl2::equals<t_vec>(Q, *Qhkl, 1e-3));

	std::cout << std::endl;
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_xtal2, t_real, t_types)
{
	std::cout << "Test 2" << std::endl;

	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	auto B = tl2::B_matrix<t_mat, t_real>(
		5., 6., 7.,
		60. / 180. * tl2::pi<t_real>,
		60. / 180. * tl2::pi<t_real>,
		60. / 180. * tl2::pi<t_real>);

	t_vec vec1 = tl2::create<t_vec>({1, 0, 0});
	t_vec vec2 = tl2::create<t_vec>({0, 1, 0});
	t_vec vec3 = tl2::cross<t_mat, t_vec>(B, vec1, vec2);
	t_mat UB = tl2::UB_matrix<t_mat, t_vec>(B, vec1, vec2, vec3);
	std::cout << "vec1 = " << vec1 << std::endl;
	std::cout << "vec2 = " << vec2 << std::endl;
	std::cout << "vec3 = " << vec3 << std::endl;
	std::cout << "UB = " << UB << std::endl;

	t_real ki = 1.4;
	t_real kf = 1.5;
	t_vec Q = tl2::create<t_vec>({1, 1, 0});
	t_real a3_offs = tl2::pi<t_real>*0.5;
	auto [anglesok, a3, a4, dist] = tl2::calc_tas_a3a4<t_mat, t_vec, t_real>(
		B, ki, kf, Q, vec1, vec3, 1., a3_offs);
	std::cout << "a3 = " << a3 / tl2::pi<t_real> * 180.
		<< ", a4 = " << a4 / tl2::pi<t_real> * 180.
		<< std::endl;
	std::cout << "distance of Q to scattering plane: " << dist << std::endl;

	// calculate back to Q
	t_real Qlen = tl2::calc_tas_Q_len<t_real>(ki, kf, a4);
	std::optional<t_vec> Qhkl = tl2::calc_tas_hkl<t_mat, t_vec, t_real>(
		B, ki, kf, Qlen, a3, vec1, vec3, 1., a3_offs);
	std::cout << "Q = " << *Qhkl << std::endl;


	BOOST_TEST(tl2::equals<t_real>(dist, 0, std::numeric_limits<t_real>::epsilon()*1e2));
	BOOST_TEST(tl2::equals<t_real>(a3/tl2::pi<t_real>*180., -15.84, 1e-3));
	BOOST_TEST(tl2::equals<t_real>(a4/tl2::pi<t_real>*180., 68.895, 1e-3));
	BOOST_TEST((Qhkl.operator bool()));
	BOOST_TEST(tl2::equals<t_vec>(Q, *Qhkl, 1e-3));

	std::cout << std::endl;
}
