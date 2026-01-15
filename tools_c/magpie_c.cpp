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

#include "libs/magdyn.h"
#include "magpie_c.h"


// types
using t_real = t_magpie_real;
using t_cplx = std::complex<t_real>;
using t_mat = tl2::mat<t_cplx>;
using t_vec = tl2::vec<t_cplx>;
using t_mat_real = tl2::mat<t_real>;
using t_vec_real = tl2::vec<t_real>;

using t_magdyn = tl2_mag::MagDyn<
	t_mat, t_vec,
	t_mat_real, t_vec_real,
	t_cplx, t_real,
	std::size_t>;



extern "C" void* magpie_create()
{
	return new t_magdyn{};
}



extern "C" void magpie_free(void *_mag)
{
	if(!_mag)
		return;

	t_magdyn *mag = reinterpret_cast<t_magdyn*>(_mag);
	delete mag;
}



extern "C" int magpie_load(void *_mag, const char* file)
{
	if(!_mag)
		return false;
	
	t_magdyn *mag = reinterpret_cast<t_magdyn*>(_mag);
	return mag->Load(file);
}



extern "C" int magpie_save_dispersion(void *_mag,
  const char* file,
	t_magpie_real h0, t_magpie_real k0, t_magpie_real l0,
	t_magpie_real h1, t_magpie_real k1, t_magpie_real l1,
	unsigned int num_pts)
{
	
	if(!_mag)
		return false;
	
	t_magdyn *mag = reinterpret_cast<t_magdyn*>(_mag);
	return mag->SaveDispersion(file, h0, k0, l0, h1, k1, l1, num_pts);
}
