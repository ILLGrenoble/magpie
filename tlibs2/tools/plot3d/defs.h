/**
 * type definitions and setting variables
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TL2_PLOTTER3D_DEFS__
#define __TL2_PLOTTER3D_DEFS__

#include <QtCore/QString>

#include <string>
#include <array>
#include <variant>
#include <optional>
#include <cstddef>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/qt/gl.h"



// ----------------------------------------------------------------------------
// type definitions
// ----------------------------------------------------------------------------
using t_size = std::size_t;

//using t_real = float;
using t_real = double;

using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;

using t_real_gl = tl2::t_real_gl;
using t_vec2_gl = tl2::t_vec2_gl;
using t_vec3_gl = tl2::t_vec3_gl;
using t_vec_gl = tl2::t_vec_gl;
using t_mat_gl = tl2::t_mat_gl;
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// number of threads for calculation
extern unsigned int g_num_threads;

// number precisions
extern int g_prec, g_prec_gui;

// epsilon
extern t_real g_eps;

// fraction of points to check if the stop button was pressed
extern unsigned int g_stop_check_fraction;

// camera fov
extern t_real g_fov;

// use native menubar and dialogs?
extern int g_use_native_menubar, g_use_native_dialogs;

// gui theme and font
extern QString g_theme, g_font;
// ----------------------------------------------------------------------------


#endif
