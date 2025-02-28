/**
 * regression test
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -o leastsq leastsq.cpp
 * g++ -std=c++20 -DUSE_LAPACK -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o leastsq leastsq.cpp -llapacke
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

#define BOOST_TEST_MODULE Least Squares Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"
using namespace tl2_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_leastsq, t_real, t_types)
{
	//using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	//using t_mat = tl2::mat<t_real, std::vector>;
	//using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	//using t_mat_cplx = tl2::mat<t_cplx, std::vector>;


	auto x = tl2::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
	auto y = tl2::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14});

	auto [params, ok] = tl2::leastsq<t_vec>(x, y, 1);
	std::cout << "ok: " << ok << ", params: " << params << std::endl;

	BOOST_TEST(ok);
	BOOST_TEST(params[0] == 3.9, testtools::tolerance(1e-3));
	BOOST_TEST(params[1] == 1.036, testtools::tolerance(1e-3));
}
