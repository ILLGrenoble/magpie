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

	#include "defs.h"

	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-4 : 1e-8;

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

	std::cout << std::endl;

	{
		t_mat M = tl2::create<t_mat>({
			1., 2., 3., 4., 5.,
			6., 7., 8., 9., 10.,
			11., 12., 13., 14., 15.,
			16., 17., 18., 19., 20.,
			21., 22., 23., 24., 25.,
		});
		std::cout << "M [ " << (void*)(&M) << " ] = " << M << std::endl;

		t_mat M2{ M };
		std::cout << "M2 = " << M2 << std::endl;
		BOOST_TEST(tl2::equals<t_mat>(M, M2, eps));

		t_mat M3{ std::move(M2) };
		std::cout << "M3 = " << M3 << std::endl;

		t_real tr = tl2::trace<t_mat>(M3);
		BOOST_TEST(tl2::equals<t_mat>(M, M3, eps));
		BOOST_TEST(tl2::equals<t_real>(tr, 1. + 7. + 13. + 19. + 25., eps));
	}

	{
		t_mat_cplx M = tl2::create<t_mat_cplx>({
			1., 2., 3., 4., 5.,
			6., 7., 8., 9., 10.,
			11., 12., 13., 14., 15.,
			16., 17., 18., 19., 20.,
			21., 22., 23., 24., 25.,
		});


		t_mat_cplx M2{ M };
		BOOST_TEST(tl2::equals<t_mat_cplx>(M, M2, eps));

		t_mat_cplx M3{ std::move(M2) };
		t_cplx tr = tl2::trace<t_mat_cplx>(M3);
		BOOST_TEST(tl2::equals<t_mat_cplx>(M, M3, eps));
		BOOST_TEST(tl2::equals<t_cplx>(tr, 1. + 7. + 13. + 19. + 25., eps));

		t_mat_cplx M4{M3.size1(), M3.size2()};
		M4 = tl2::herm(M3);
		M4 = tl2::herm(M4);
		BOOST_TEST(tl2::equals<t_mat_cplx>(M3, M4, eps));
	}


	{
		t_mat M = tl2::rand<t_mat>(64, 64);

		t_mat M2{ M };
		BOOST_TEST(tl2::equals<t_mat>(M, M2, eps));

		t_mat M3{ std::move(M2) };
		BOOST_TEST(tl2::equals<t_mat>(M, M3, eps));
	}
}
