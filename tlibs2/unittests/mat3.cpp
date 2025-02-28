/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 30-jun-2024
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

#define BOOST_TEST_MODULE Mat3
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_mat3, t_real, t_types)
{
	using namespace tl2_ops;

	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	const t_real eps = 1e-4;

	// rotation trafo test
	{
		t_vec vec_from = tl2::create<t_vec>({
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
		});
		t_vec vec_to = tl2::create<t_vec>({
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
		});

		vec_from /= tl2::norm<t_vec>(vec_from);
		vec_to /= tl2::norm<t_vec>(vec_to);

		t_mat mat = tl2::rotation<t_mat, t_vec, t_real>(vec_from, vec_to);

		std::cout << "Matrix to rotate\n\t" << vec_from
			<< " into\n\t" << vec_to << ":"
			<< "\n\t" << mat << std::endl;

		t_vec vec_to_calc = mat * vec_from;
		BOOST_TEST(tl2::equals<t_vec>(vec_to, vec_to_calc, eps));
	}

	// arrow trafo test
	{
		t_vec vec_from = tl2::create<t_vec>({
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
		});
		t_vec vec_to = tl2::create<t_vec>({
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
			tl2::get_rand<t_real>(-1, 1),
		});

		vec_from /= tl2::norm<t_vec>(vec_from);
		vec_to /= tl2::norm<t_vec>(vec_to);

		t_vec vec_from_hom = tl2::create<t_vec>(
			{ vec_from[0], vec_from[1], vec_from[2], 1 });
		t_vec vec_to_hom = tl2::create<t_vec>(
			{ vec_to[0], vec_to[1], vec_to[2], 1 });

		t_mat mat = tl2::get_arrow_matrix<t_vec, t_mat, t_real>(
			vec_to,                           // to
			1,                                // post-scale
			tl2::create<t_vec>({ 0, 0, 0 }),  // post-translate
			vec_from,                         // from
			1,                                // pre-scale
			tl2::create<t_vec>({ 0, 0, 0 })); // pre-translate

		std::cout << "Matrix to rotate\n\t" << vec_from
			<< " into\n\t" << vec_to << ":"
			<< "\n\t" << mat << std::endl;

		t_vec vec_to_calc = mat * vec_from_hom;
		BOOST_TEST(tl2::equals<t_vec>(vec_to_hom, vec_to_calc, eps));
	}
}
