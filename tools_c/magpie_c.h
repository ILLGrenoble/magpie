/**
 * magnetic dynamics c library interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15-january-2026
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


#ifndef __MAGDYN_C__
#define __MAGDYN_C__


typedef void* t_magpie;
typedef double t_magpie_real;


#ifdef __cplusplus
extern "C" {
#endif


// initialise the magpie library
t_magpie magpie_create();

// free the magpie library
void magpie_free(t_magpie _mag);

// load a magnetic model
int magpie_load(t_magpie _mag, const char* file);

// count the number of sites in the current magnetic model
unsigned int magpie_site_count(t_magpie _mag);

// get the maximum number of magnon branches in the current magnetic model
unsigned int magpie_branch_count(t_magpie _mag);

// calculate the energies and spin-spin correlation at the point Q = (hkl)
unsigned int magpie_calc_energies(t_magpie _mag,
	t_magpie_real h, t_magpie_real k, t_magpie_real l,
	t_magpie_real* Es, t_magpie_real* ws);

// calculate S(h, k, l, E)
t_magpie_real magpie_calc_S(t_magpie _mag,
	t_magpie_real h, t_magpie_real k, t_magpie_real l, t_magpie_real E,
	t_magpie_real sigma);

// save a dispersion branch going from (h0 k0 l0) to (h1 k1 l1)
int magpie_save_dispersion(t_magpie _mag,
	const char* file,
	t_magpie_real h0, t_magpie_real k0, t_magpie_real l0,
	t_magpie_real h1, t_magpie_real k1, t_magpie_real l1,
	unsigned int num_pts);


#ifdef __cplusplus
}
#endif

#endif
