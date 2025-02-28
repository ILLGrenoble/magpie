/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-21
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

#define BOOST_TEST_MODULE Mat2
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_mat2, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	const t_real eps = 1e-4;


	{
		auto M = tl2::create<t_mat>({
			1., 2., 3.,
			3., 1., 4.,
			9., -4., 2.
		});

		// test determinant
		t_real det = tl2::det(M);
		std::cout << "M = " << M << std::endl;
		std::cout << "|M| = " << det << std::endl;

		BOOST_TEST(tl2::equals<t_real>(det, 15., eps));


		// test splitting
		auto [M_S, M_A] = tl2::split_symm(M);
		BOOST_TEST(tl2::equals<t_mat>(M, M_S + M_A, eps));
	}


	{
		auto M = tl2::create<t_mat>({
			1., 3., 2.,
			3., 4., 1.,
			9., 2., -4.
		});

		// test determinant
		t_real det = tl2::det(M);
		std::cout << "M = " << M << std::endl;
		std::cout << "|M| = " << det << std::endl;

		BOOST_TEST(tl2::equals<t_real>(det, -15., eps));


		// test splitting
		auto [M_S, M_A] = tl2::split_symm(M);
		BOOST_TEST(tl2::equals<t_mat>(M, M_S + M_A, eps));
	}


	{
		auto M = tl2::create<t_mat_cplx>({
			1., 2., 3.,
			3., 1., 4.,
			9., -4., 2
		});

		// test determinant
		t_cplx det = tl2::det(M);
		std::cout << "M = " << M << std::endl;
		std::cout << "|M| = " << det << std::endl;

		BOOST_TEST(tl2::equals<t_cplx>(det, 15., eps));


		// test splitting
		auto [M_S, M_A] = tl2::split_symm(M);
		BOOST_TEST(tl2::equals<t_mat_cplx>(M, M_S + M_A, eps));
	}
}
