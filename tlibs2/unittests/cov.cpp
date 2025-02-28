/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -I.. -o cov cov.cpp
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

#define BOOST_TEST_MODULE Cov1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"
using namespace tl2_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	std::vector<t_vec> vecs{{
		tl2::create<t_vec>({1, 2, 3}),
		tl2::create<t_vec>({3, 2, 1}),
		tl2::create<t_vec>({5, -2, 0.5}),
	}};

	std::vector<t_real> probs{{ 1., 1., 0.5 }};

	auto [cov, corr] = tl2::covariance<t_mat, t_vec>(vecs, &probs);
	std::cout << "cov  = " << cov << std::endl;
	std::cout << "corr = " << corr << std::endl;

	// comparing with: np.cov([[1,2,3], [3,2,1], [5,-2,0.5]], rowvar=False, aweights=[1,1,0.5], ddof=0)
	const t_mat cov_cmp = tl2::create<t_mat>({
		 2.24, -1.92, -1.52,
		-1.92,  2.56,  0.96,
		-1.52,  0.96,  1.16
	});
	BOOST_TEST(tl2::equals(cov, cov_cmp, 1e-5));
}
