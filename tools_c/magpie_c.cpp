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
#include "libs/maths.h"
#include "libs/defs.h"
#include "magpie_c.h"

#include <mutex>



// ----------------------------------------------------------------------------
// internals
// ----------------------------------------------------------------------------
// types

using t_E_and_S = std::vector<std::pair<t_real, t_real>>;

// cache S(Q, E) results
using t_rtree_vertex = boost::geometry::model::point<t_real, 3, boost::geometry::cs::cartesian>;
using t_rtree_leaf = std::tuple<t_rtree_vertex, t_E_and_S>;
using t_rtree = boost::geometry::index::rtree<t_rtree_leaf, boost::geometry::index::linear<16>>;

// magpie core
using t_magdyn = magdyn::MagDyn<
	t_mat, t_vec,
	t_mat_real, t_vec_real,
	t_cplx, t_real,
	std::size_t>;



/**
 * main magpie data structure
 */
struct MagpieData
{
	// calculation core
	t_magdyn mag;

	// cache S(Q, E) results
	t_rtree rtree_S{ };
	t_real rtree_rlu_eps{ 1e-4 };
	bool use_rtree{ true };
	std::mutex rtree_mutex;
};



/**
 * internal function to calculate E and S
 */
static t_E_and_S _calc_energies(t_magpie _mag,
	t_real h, t_real k, t_real l, bool ignore_weights = false)
{
	if(!_mag)
		return t_E_and_S{};

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);

	t_vec_real Q_rlu = tl2::create<t_vec_real>({ h, k, l });


	// look for cached S(Q, E) value
	t_rtree_vertex rtree_vert;
	if(dat->use_rtree)
	{
		tl2::to_geo_vertex<t_rtree_vertex, t_vec_real>(
			rtree_vert, Q_rlu, std::make_index_sequence<3>());

		std::vector<t_rtree_leaf> rtree_nearest;

		{
			std::lock_guard<std::mutex> _lck{dat->rtree_mutex};
			dat->rtree_S.query(boost::geometry::index::nearest(rtree_vert, 1), std::back_inserter(rtree_nearest));
		}

		if(rtree_nearest.size() > 0)
		{
			t_vec_real nearest_point =
				tl2::from_geo_vertex<t_vec_real, t_rtree_vertex>(
					std::get<0>(rtree_nearest[0]), std::make_index_sequence<3>());

			// use cached value if within epsilon distance
			if(tl2::norm(Q_rlu - nearest_point) < dat->rtree_rlu_eps)
			{
				//std::cout << "Using cached value at Q = (" << h << ", " << k << ", " << l << ")." << std::endl;
				return std::get<1>(rtree_nearest[0]);
			}
		}
	}


	// calculate energies and weights
	const auto _Es_and_S = dat->mag.CalcEnergies(Q_rlu, ignore_weights).E_and_S;

	t_E_and_S Es_and_S;
	Es_and_S.reserve(_Es_and_S.size());

	for(std::size_t branch_idx = 0; branch_idx < _Es_and_S.size(); ++branch_idx)
	{
		Es_and_S.emplace_back(std::make_pair(
			_Es_and_S[branch_idx].E, _Es_and_S[branch_idx].weight_perp));
	}


	// cache S(Q, E) values
	if(dat->use_rtree && !ignore_weights)
	{
		std::lock_guard<std::mutex> _lck{dat->rtree_mutex};
		dat->rtree_S.insert(std::make_tuple(rtree_vert, Es_and_S));
	}

	return Es_and_S;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// external API functions
// ----------------------------------------------------------------------------
/**
 * initialise the magpie library
 */
extern "C"
t_magpie magpie_create()
{
	MagpieData *dat = new MagpieData();

	dat->mag.SetPerformChecks(false);
	dat->mag.SetSilent(false);

	return dat;
}



/**
 * free the magpie library
 */
extern "C"
void magpie_free(t_magpie _mag)
{
	if(!_mag)
		return;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	delete dat;
}



/**
 * enable caching energies and weights
 */
extern "C"
void magpie_use_cache(t_magpie _mag, int enable)
{
	if(!_mag)
		return;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	dat->use_rtree = (enable != 0);
}



/**
 * load a magnetic model
 */
extern "C"
int magpie_load(t_magpie _mag, const char* file)
{
	if(!_mag)
		return false;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);

	dat->rtree_S.clear();
	return dat->mag.Load(file);
}



/**
 * count the number of sites in the current magnetic model
 */
extern "C"
unsigned int magpie_site_count(t_magpie _mag)
{
	if(!_mag)
		return 0;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);

	return static_cast<unsigned int>(dat->mag.GetMagneticSitesCount());
}



/**
 * get the maximum number of magnon branches in the current magnetic model
 */
extern "C"
unsigned int magpie_branch_count(t_magpie _mag)
{
	if(!_mag)
		return 0;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);

	// each magnetic site leads to a magnon branch
	std::size_t max_branch_count = dat->mag.GetMagneticSitesCount();

	// include magnon creation and annihilation
	max_branch_count *= 2;

	if(dat->mag.IsIncommensurate())
	{
		// include non-commensurate structures
		max_branch_count *= 3;
	}

	return static_cast<unsigned int>(max_branch_count);
}



/**
 * get the crystal lattice, length in angstrom, angles in radians
 */
extern "C"
void magpie_get_lattice(t_magpie _mag,
  t_magpie_real *a, t_magpie_real *b, t_magpie_real *c,
	t_magpie_real *alpha, t_magpie_real *beta, t_magpie_real *gamma)
{
	if(!_mag)
		return;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);

	t_vec_real latt = dat->mag.GetCrystalLattice();

	if(a) *a = latt[0];
	if(b) *b = latt[1];
	if(c) *c = latt[2];
	
	if(alpha) *alpha = latt[3];
	if(beta) *beta = latt[4];
	if(gamma) *gamma = latt[5];
}



/**
 * set temperature
 */
extern "C"
void magpie_set_temperature(t_magpie _mag, t_real T)
{
	if(!_mag)
		return;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	dat->mag.SetTemperature(T);
}



/**
 * set external magnetic field
 */
extern "C"
void magpie_set_field(t_magpie _mag, t_real B,
	t_real Bx, t_real By, t_real Bz,
	int align_spins, int keep_signs)
{
	if(!_mag)
		return;

	typename t_magdyn::ExternalField field{};
	field.mag = B;
	field.dir = tl2::create<t_vec_real>({ Bx, By, Bz });
	field.align_spins = (align_spins != 0);
	field.keep_signs = (keep_signs != 0);

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	dat->mag.SetExternalField(field);
}



/**
 * calculate the energies and spin-spin correlation at the point Q = (hkl)
 */
extern "C"
unsigned int magpie_calc_energies(t_magpie _mag,
	t_real h, t_real k, t_real l,
	t_real* Es, t_real* ws)
{
	if(!_mag)
		return 0;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	const t_E_and_S Es_and_S = _calc_energies(_mag, h, k, l, !ws);

	const std::size_t max_branches = 2*dat->mag.GetMagneticSitesCount();
	const std::size_t num_branches = std::min(Es_and_S.size(), max_branches);

	for(std::size_t branch_idx = 0; branch_idx < num_branches; ++branch_idx)
	{
		if(Es)
			Es[branch_idx] = Es_and_S[branch_idx].first;
		if(ws)
			ws[branch_idx] = Es_and_S[branch_idx].second;
	}

	return static_cast<unsigned int>(num_branches);
}



/**
 * calculate S(h, k, l, E)
 */
extern "C"
t_real magpie_calc_S(t_magpie _mag,
	t_real h, t_real k, t_real l, t_real E,
	t_real sigma)
{
	if(!_mag)
		return 0;

	const t_E_and_S Es_and_S = _calc_energies(_mag, h, k, l, false);

	t_real S = 0.;
	for(std::size_t branch_idx = 0; branch_idx < Es_and_S.size(); ++branch_idx)
	{
		const t_real E_branch = Es_and_S[branch_idx].first;
		const t_real S_branch = Es_and_S[branch_idx].second;

		S += S_branch * tl2::gauss_model(E, E_branch, sigma, 1., 0.);
	}

	return S;
}



/**
 * save a dispersion branch going from (h0 k0 l0) to (h1 k1 l1)
 */
extern "C"
int magpie_save_dispersion(t_magpie _mag,
	const char* file,
	t_real h0, t_real k0, t_real l0,
	t_real h1, t_real k1, t_real l1,
	unsigned int num_pts)
{
	if(!_mag)
		return false;

	MagpieData *dat = reinterpret_cast<MagpieData*>(_mag);
	return dat->mag.SaveDispersion(file, h0, k0, l0, h1, k1, l1, num_pts);
}
// ----------------------------------------------------------------------------
