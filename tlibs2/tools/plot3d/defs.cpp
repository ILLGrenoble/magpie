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

#include <QtCore/QSettings>
#include <thread>

#include "defs.h"



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// number of threads for calculation
unsigned int g_num_threads = std::max<unsigned int>(
	1, std::thread::hardware_concurrency()/2);

// epsilons and precisions
int g_prec = 6;
int g_prec_gui = 3;
t_real g_eps = 1e-6;

// fraction of points to check if the stop button was pressed
unsigned int g_stop_check_fraction = 100;

t_real g_fov = 90.;

// gui theme
QString g_theme = "Fusion";

// gui font
QString g_font = "";

// use native menu bar?
int g_use_native_menubar = 0;

// use native dialogs?
int g_use_native_dialogs = 0;
// ----------------------------------------------------------------------------
