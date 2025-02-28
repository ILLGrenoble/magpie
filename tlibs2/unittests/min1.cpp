/**
 * minimisation test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 25-jul-24
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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



// minimise lambda function with fixed args
BOOST_AUTO_TEST_CASE_TEMPLATE(test_min1, t_real, t_types_real)
{
	auto func = [](t_real x)
	{
		return (x - 5.67) * (x - 5.67);
	};

	std::vector<std::string> params{{ "x", }};
	std::vector<t_real> vals{{ 0., }};
	std::vector<t_real> errs{{ 5, }};
	std::vector<bool> fixed{ false, };

	bool ok = tl2::minimise<t_real, 1>(func, params, vals, errs, &fixed);

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals<t_real>(vals[0], 5.67, 1e-3));
}



// minimise lambda function with non-fixed args
BOOST_AUTO_TEST_CASE_TEMPLATE(test_min1_dynarg, t_real, t_types_real)
{
	auto func = [](const std::vector<tl2::t_real_min>& args)
	{
		const tl2::t_real_min& x = args[0];
		return (x - 5.67) * (x - 5.67);
	};

	std::vector<std::string> params{{ "x", }};
	std::vector<t_real> vals{{ 0., }};
	std::vector<t_real> errs{{ 5, }};
	std::vector<bool> fixed{ false, };

	bool ok = tl2::minimise_dynargs<t_real>(1, func, params, vals, errs, &fixed);

	BOOST_TEST(ok);
	BOOST_TEST(tl2::equals<t_real>(vals[0], 5.67, 1e-3));
}
