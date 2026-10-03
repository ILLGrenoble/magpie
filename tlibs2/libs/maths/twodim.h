/**
 * tlibs2 -- maths library -- 2-dim algos
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_MATHS_2D_H__
#define __TLIBS2_MATHS_2D_H__

#include <cmath>
#include <tuple>
#include <vector>
#include <limits>

#include "decls.h"
#include "constants.h"
#include "projectors.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// 2-dim algos
// ----------------------------------------------------------------------------

/**
 * 2-dim cross product
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec, typename t_scalar = typename t_vec::value_type>
t_scalar cross_2d(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	t_scalar prod = t_scalar(0);

	// only valid for 2-vectors -> use first two components
	if(vec1.size() < 2 || vec2.size() < 2)
		return prod;

	return vec1[0]*vec2[1] - vec1[1]*vec2[0];
}


/**
 * SO(2) rotation matrix
 * @see https://en.wikipedia.org/wiki/Rotation_matrix
 */
template<class t_mat>
t_mat rotation_2d(const typename t_mat::value_type angle)
requires tl2::is_mat<t_mat>
{
	return givens<t_mat>(2, 0, 1, angle);
}


/**
 * polygon area
 * @see https://en.wikipedia.org/wiki/Shoelace_formula
 */
template<class t_vecs, class t_vec = typename t_vecs::value_type,
	class t_scalar = typename t_vec::value_type>
t_scalar area_2d(const t_vecs& vecs, bool is_loop = false /* last vertex == first vertex? */)
requires is_basic_vec<t_vec>
{
	t_scalar area = t_scalar(0);

	const std::size_t size = vecs.size();
	const std::size_t num_pts = is_loop ? size - 1 : size;
	for(std::size_t idx = 0; idx < num_pts; ++idx)
	{
		const t_vec& vec1 = vecs[idx];
		const t_vec& vec2 = vecs[(idx + 1) % size];

		area += cross_2d<t_vec, t_scalar>(vec1, vec2);
	}

	return area / t_scalar(2);
}
// ----------------------------------------------------------------------------

}
#endif
