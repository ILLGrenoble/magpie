/**
 * brillouin zone tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date Maz-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include "globals.h"
#include <cmath>


/**
 * will give wrong results if epsilon is too high
 * (check with space group #146)
 */
t_real g_eps_bz = 1e-7;
int g_prec_bz = 7;
int g_prec_gui_bz = 4;


/**
 * sets new epsilon and precision values
 */
void set_eps_bz(t_real eps, int prec)
{
	// determine precision from epsilon
	if(prec < 0)
		prec = int(-std::log10(eps));

	g_eps_bz = eps;
	g_prec_bz = prec;
}
