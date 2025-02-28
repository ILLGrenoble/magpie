/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 20-aug-20
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

#define BOOST_TEST_MODULE Stat1
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

	std::vector<t_vec> vecs = {{
		tl2::create<t_vec>({-1, 60}),
		tl2::create<t_vec>({20, 5}),
		tl2::create<t_vec>({-3, 40}),
		tl2::create<t_vec>({40, 3}),
		tl2::create<t_vec>({-5, 20}),
		tl2::create<t_vec>({60, 1}),
	}};

	auto [cov, cor] = tl2::covariance<t_mat, t_vec>(vecs);
	std::cout << "cov = " << cov << std::endl;
	std::cout << "cor = " << cor << std::endl;

	//BOOST_TEST(tl2::equals(vec, vec2, 1e-5));
}
