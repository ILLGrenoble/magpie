/**
 * minimisation test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 27-aug-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -I.. -o min min0.cpp ../libs/log.cpp -lMinuit2
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

#include <string>
#include <iostream>

#define BOOST_TEST_MODULE Min Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include "libs/fit.h"


using t_types_real = std::tuple<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_real, t_real, t_types_real)
{
	std::string func{"(x - 5.67)^2"};

	std::vector<std::string> params{{"x"}};
	std::vector<t_real> vals{{0.}};
	std::vector<t_real> errs{{0.1}};
	std::vector<bool> fixed{false};

	bool ok = tl2::minimise_expr(func, params, vals, errs, &fixed);

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals<t_real>(vals[0], 5.67, 1e-3));
}
