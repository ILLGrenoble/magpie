/**
 * test cholesky decomposition
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o chol chol.cpp -llapacke
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

#define BOOST_TEST_MODULE Cholesky Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/algos.h"
#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_chol, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	t_real eps = 1e-2;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};


	// real version
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		auto org_mat = tl2::zero<t_mat>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				org_mat(i, j) = rnddist(rndgen);

		t_mat mat = tl2::trans(org_mat) * org_mat;
		//std::cout << "C0  = " << org_mat << std::endl;
		std::cout << "M   = " << mat << std::endl;

		auto [ok, C] = tl2_la::chol<t_mat, t_vec>(mat);
		auto [ok2, C2, P] = tl2_la::chol_semi<t_mat, t_vec>(mat, eps);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;
		std::cout << "ok 2 = " << std::boolalpha << ok2 << std::endl;

		t_mat CtC = tl2::trans(C) * C;
		std::cout << "C    = " << C << std::endl;
		std::cout << "C^tC = " << CtC << std::endl;

		t_mat CtC2 = tl2::trans(C2) * C2;
		std::cout << "C2      = " << C2 << std::endl;
		std::cout << "C2^t C2 = " << CtC2 << std::endl;
		t_mat mat_pivot = tl2::trans(P) * mat * P;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals<t_mat>(mat, CtC, eps));
		BOOST_TEST(tl2::equals<t_mat>(mat_pivot, CtC2, eps));
		//BOOST_TEST(tl2::equals<t_mat>(org_mat, C, eps));
		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}


	// complex version
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		auto org_mat = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				org_mat(i,j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};

		t_mat_cplx mat = tl2::herm(org_mat) * org_mat;
		std::cout << "M   = " << mat << std::endl;

		auto [ok, C] = tl2_la::chol<t_mat_cplx, t_vec_cplx>(mat);
		auto [ok2, C2, P] = tl2_la::chol_semi<t_mat_cplx, t_vec_cplx>(mat, eps);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;
		std::cout << "ok 2 = " << std::boolalpha << ok2 << std::endl;

		t_mat_cplx ChC = tl2::herm(C) * C;
		std::cout << "C    = " << C << std::endl;
		std::cout << "C^hC = " << ChC << std::endl;

		t_mat_cplx ChC2 = tl2::herm(C2) * C2;
		std::cout << "C2      = " << C2 << std::endl;
		std::cout << "C2^h C2 = " << ChC2 << std::endl;
		t_mat_cplx mat_pivot = tl2::trans(P) * mat * P;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals<t_mat_cplx>(mat, ChC, eps));
		BOOST_TEST(tl2::equals<t_mat_cplx>(mat_pivot, ChC2, eps));

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}
}
