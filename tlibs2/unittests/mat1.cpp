/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -o mat1 mat1.cpp
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

#define BOOST_TEST_MODULE Mat1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"
using namespace tl2_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_mat1, t_real, t_types)
{
	#include "defs.h"

	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-4 : 1e-8;


	std::cout << tl2::stoval<unsigned int>(std::string("123")) << std::endl;

	std::vector vec1{{
		tl2::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
		tl2::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14})
	}};
	std::vector vec2{{
		tl2::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10.5}),
		tl2::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14})
	}};
	std::vector vec3{{
		tl2::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10}),
		tl2::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14, 14})
	}};

	BOOST_TEST(tl2::equals_all(vec1, vec1, eps));
	BOOST_TEST(tl2::equals_all(vec3, vec3, eps));
	BOOST_TEST(!tl2::equals_all(vec1, vec2, eps));
	BOOST_TEST(!tl2::equals_all(vec1, vec3, eps));

	t_vec vec4{ vec1[0] };
	BOOST_TEST(tl2::equals<t_vec>(vec4, vec1[0], eps));

	t_vec vec5{ std::move(vec4) };
	BOOST_TEST(tl2::equals<t_vec>(vec5, vec1[0], eps));
	BOOST_TEST(vec5.size() = vec1[0].size());
}
