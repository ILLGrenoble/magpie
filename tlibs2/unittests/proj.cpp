/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 5-jul-2027
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#define BOOST_TEST_MODULE Proj
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

	const t_real eps = 1e-6;

	// projection axis
	t_vec axis = tl2::create<t_vec>({
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1),
	});

	// vector to project
	t_vec vec = tl2::create<t_vec>({
		tl2::get_rand<t_real>(-10, 10),
		tl2::get_rand<t_real>(-10, 10),
		tl2::get_rand<t_real>(-10, 10),
	});

	axis /= tl2::norm<t_vec>(axis);
	//vec /= tl2::norm<t_vec>(vec);

	t_mat matProj = tl2::projector<t_mat, t_vec>(axis, false);
	t_mat matOrthoProj = tl2::ortho_projector<t_mat, t_vec>(axis, false);
	BOOST_TEST(tl2::equals<t_real>(tl2::trace(matProj), 1, eps));
	BOOST_TEST(tl2::equals<t_real>(tl2::trace(matOrthoProj), 2, eps));

	t_vec vec_proj = matProj * vec;
	t_vec vec_ortho_proj = matOrthoProj * vec;
	BOOST_TEST(tl2::equals<t_vec>(vec_proj + vec_ortho_proj, vec, eps));
}
