/**
 * basic type definitions
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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


#define MAGCORE_VER "2.8.5"


#endif
