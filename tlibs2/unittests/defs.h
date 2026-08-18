/**
 * type definitions
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
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


//using t_real = double;
using t_vec = tl2::vec<t_real, std::vector, tl2::alloc_noinit, 4>;
using t_vec44 = tl2::vec<t_real, std::vector, tl2::alloc_noinit, 4*4>;
using t_mat = tl2::mat<t_real, t_vec44>;

using t_vec4 = tl2::vec_static<t_real, 4>;
using t_mat44 = tl2::mat_static<t_real, 4, 4>;

using t_cplx = std::complex<t_real>;
using t_vec_cplx = tl2::vec<t_cplx, std::vector, tl2::alloc_noinit, 4>;
using t_vec_cplx44 = tl2::vec<t_cplx, std::vector, tl2::alloc_noinit, 4*4>;
using t_mat_cplx = tl2::mat<t_cplx, t_vec_cplx44>;
