/**
 * units test
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
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

#define BOOST_TEST_MODULE La1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/units.h"

using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	std::cout << "kB = " << tl2::kB<t_real> << std::endl;
	std::cout << "muB = " << tl2::mu_B<t_real> << std::endl;

	std::cout << "kB = "
		<< static_cast<t_real>(tl2::kB<t_real>/tl2::meV<t_real>*tl2::kelvin<t_real>)
		<< std::endl;
	std::cout << "muB = "
		<< static_cast<t_real>(tl2::mu_B<t_real>/tl2::meV<t_real>*tl2::tesla<t_real>)
		<< std::endl;
}
