/**
 * test cholesky decomposition
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o chol2 chol2.cpp -llapacke
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
		auto mat = tl2::zero<t_mat>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=i; j<dim; ++j)
				mat(j, i) = mat(i, j) = rnddist(rndgen);

		std::cout << "M    = " << mat << std::endl;

		bool sym = 1;
		auto [ok, C, D] =
			tl2_la::chol_herm<t_mat, t_vec>(mat);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		t_mat CDCt = C * D * tl2::trans(C);
		std::cout << "C     = " << C << std::endl;
		std::cout << "CDC^t = " << CDCt << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_mat>(mat, CDCt, eps));
		//BOOST_TEST(tl2::equals<t_mat>(org_mat, C, eps));
		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}


	// complex version
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		auto mat = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
		{
			for(std::size_t j=i; j<dim; ++j)
			{
				mat(i,j) = rnddist(rndgen);
				if(i != j)
					mat(i,j) += rnddist(rndgen)*t_cplx{0, 1};
				mat(j,i) = std::conj(mat(i,j));
			}
		}

		std::cout << "M    = " << mat << std::endl;

		auto [ok, C, D] =
			tl2_la::chol_herm<t_mat_cplx, t_vec_cplx>(mat);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		t_mat_cplx CDCh = C * D * tl2::herm(C);
		std::cout << "C     = " << C << std::endl;
		std::cout << "CDC^h = " << CDCh << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_mat_cplx>(mat, CDCh, eps));

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}
}
