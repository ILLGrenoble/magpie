/**
 * fit test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-aug-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -I.. -o fit0 fit0.cpp ../libs/log.cpp -lMinuit2
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

#define BOOST_TEST_MODULE Fit Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;


#include "libs/fit.h"
#include "libs/maths.h"


using t_types_real = std::tuple<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_real, t_real, t_types_real)
{
	t_real amp = 1.;
	t_real freq = 2*tl2::pi<t_real>;
	t_real offs = 12.;

	std::vector<t_real> xs, ys, yerrs;
	for(t_real x=0.; x<1.; x+=0.05)
	{
		xs.push_back(x);
		ys.push_back(amp*std::sin(freq*x) + offs);
		yerrs.push_back(0.1);
	}

	std::string func{"amp*sin(freq*x) + offs"};

	std::vector<std::string> params{{"amp", "freq", "offs"}};
	std::vector<t_real> vals{{amp*t_real{1.2}, freq*t_real{0.8}, offs*t_real{1.}}};
	std::vector<t_real> errs{{0.5, 0.1, 1.}};
	std::vector<bool> fixed{false, false, true};

	bool ok = tl2::fit_expr(func, xs, ys, yerrs, "x", params, vals, errs, &fixed);

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals<t_real>(vals[0], amp, 1e-3));
	BOOST_TEST(tl2::equals<t_real>(vals[1], freq, 1e-3));
	BOOST_TEST(tl2::equals<t_real>(vals[2], offs, 1e-3));
}
