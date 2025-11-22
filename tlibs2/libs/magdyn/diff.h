/**
 * tlibs2 -- magnetic dynamics -- differentiation / group velocity calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS2_MAGDYN_DIFF_H__
#define __TLIBS2_MAGDYN_DIFF_H__


#include "magdyn.h"



// --------------------------------------------------------------------
// differentiation / group velocity calculations
// --------------------------------------------------------------------

namespace tl2_mag {


/**
 * differentiate the dispersion branches (group velocities)
 */
MAGDYN_TEMPL
std::tuple<std::vector<t_real>, MAGDYN_TYPE::SofQE> MAGDYN_INST::CalcGroupVelocities(
	const t_vec_real& Q, const t_vec_real& deltaQ, const std::vector<t_size> *perm) const
{
	//SetUniteDegenerateEnergies(false);

	SofQE S1 = CalcEnergies(Q, false);
	SofQE S2 = CalcEnergies(Q + deltaQ, false);

	if(perm)
	{
		//M = tl2::reorder_cols<t_mat, t_vec>(M, *perm);
		S1.E_and_S = tl2::reorder(S1.E_and_S, *perm);
		S2.E_and_S = tl2::reorder(S2.E_and_S, *perm);
	}

	const t_size num_bands = std::min(S1.E_and_S.size(), S2.E_and_S.size());
	const t_real dQ = tl2::norm<t_vec_real>(deltaQ);

	std::vector<t_real> diffs;
	diffs.reserve(num_bands);
	for(t_size idx = 0; idx < num_bands; ++idx)
	{
		t_real v = (S2.E_and_S[idx].E - S1.E_and_S[idx].E) / dQ;
		diffs.push_back(v);
	}

	return std::make_tuple(diffs, S1);
}


}
// --------------------------------------------------------------------

#endif
