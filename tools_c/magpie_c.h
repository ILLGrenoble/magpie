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


typedef double t_magpie_real;


#ifdef __cplusplus
extern "C" {
#endif


void* magpie_create();


void magpie_free(void *_mag);


int magpie_load(void *_mag, const char* file);


int magpie_save_dispersion(void *_mag,
 const char* file,
 t_magpie_real h0, t_magpie_real k0, t_magpie_real l0,
 t_magpie_real h1, t_magpie_real k1, t_magpie_real l1,
 unsigned int num_pts);


#ifdef __cplusplus
}
#endif

#endif
