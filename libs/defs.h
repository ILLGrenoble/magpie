/**
 * basic type definitions
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
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

#ifndef __MAGCORE_DEFS__
#define __MAGCORE_DEFS__

#include "vers.h"
#include "tlibs2/libs/maths.h"


#include <cstddef>
#include <complex>


// ----------------------------------------------------------------------------
// basic type definitions
// ----------------------------------------------------------------------------
using t_size = std::size_t;

//using t_real = float;
using t_real = double;

using t_cplx = std::complex<t_real>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// vector and matrix type definitions
// ----------------------------------------------------------------------------
using t_vec_real = tl2::vec<t_real, std::vector, __TLIBS2_DEFAULT_ALLOC__, 4>;

//using t_vec_real_for_mat = tl2::vec<t_real, std::vector, __TLIBS2_DEFAULT_ALLOC__, 8*8>;
//using t_vec_real_for_mat = tl2::vec<t_real, std::vector, std::allocator, 8*8>;
//using t_mat_real = tl2::mat<t_real, t_vec_real_for_mat>;
using t_mat_real = tl2::mat<t_real, std::vector<t_real, __TLIBS2_DEFAULT_ALLOC__<t_real>>>;
//using t_mat_real = tl2::mat<t_real, std::vector<t_real>>;


using t_vec = tl2::vec<t_cplx, std::vector, __TLIBS2_DEFAULT_ALLOC__, 4>;

//using t_vec_for_mat = tl2::vec<t_cplx, std::vector, __TLIBS2_DEFAULT_ALLOC__, 8*8>;
//using t_vec_for_mat = tl2::vec<t_cplx, std::vector, std::allocator, 8*8>;
//using t_mat = tl2::mat<t_cplx, t_vec_for_mat>;
using t_mat = tl2::mat<t_cplx, std::vector<t_cplx, __TLIBS2_DEFAULT_ALLOC__<t_cplx>>>;
//using t_mat = tl2::mat<t_cplx, std::vector<t_cplx>>;
// ----------------------------------------------------------------------------


#endif
