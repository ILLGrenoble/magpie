/**
 * magnetic dynamics -- landau-lifshitz approach
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8 November 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *
 * @note Forked on 5-November-2024 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 * @desc For references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAGDYN_LANDAU_H__
#define __MAGDYN_LANDAU_H__


#include "magdyn.h"


namespace magdyn {

/**
 * get fluctuation matrix
 */
MAGDYN_TEMPL
t_mat MAGDYN_INST::CalcFluctuationMatrix(const t_vec3_real& Qvec) const
{
	t_Jmap _J_Q0;
	if(!m_J_Q0_valid)
		_J_Q0 = CalcReciprocalJs();
	const t_Jmap& J_Q0 = m_J_Q0_valid ? m_J_Q0 : _J_Q0;
	const t_Jmap J_Q = CalcReciprocalJs(Qvec);
	const t_size N = GetMagneticSitesCount();

	t_mat F = tl2::zero<t_mat>(N*3, N*3);

	// Q-independent parts
	for(const auto& [sites, J] : J_Q0)
	{
		const t_size idx1 = sites.first;
		const t_size idx2 = sites.second;

		const MagneticSite& s_i = GetMagneticSite(idx1);
		const MagneticSite& s_j = GetMagneticSite(idx2);
		const t_vec3& v_i = s_i.trafo_z_calc;
		const t_vec3& v_j = s_j.trafo_z_calc;
		t_cplx JS = -tl2::inner_noconj(v_i, J * v_j);

		t_cplx field_energy = -CalcFieldEnergy(idx1);

		for(int i = 0; i < 3; ++i)
			F(idx1*3 + i, idx1*3 + i) += JS + field_energy;;
	}

	// Q-dependent parts
	for(const auto& [sites, J] : J_Q)
	{
		const t_size idx1 = sites.first;
		const t_size idx2 = sites.second;

		for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			F(idx1*3 + i, idx2*3 + j) += J(i, j);
	}

	//std::cout << "F =\n";
	//tl2::niceprint(std::cout, F, 1e-4, 4);
	return F;
}



/**
 * get Mx matrix
 */
MAGDYN_TEMPL
t_mat MAGDYN_INST::CalcMCrossMatrix() const
{
	const t_size N = GetMagneticSitesCount();
	t_mat Mx = tl2::zero<t_mat>(N*3, N*3);

	for(t_size site_idx = 0; site_idx < N; ++site_idx)
	{
		const MagneticSite& s = GetMagneticSite(site_idx);
		const t_vec3& v = s.trafo_z_calc;
		t_mat33 matx = tl2::skewsymmetric<t_mat33, t_vec3>(v);

		for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			Mx(site_idx*3 + i, site_idx*3 + j) = matx(i, j);
	}

	//std::cout << "Mx =\n";
	//tl2::niceprint(std::cout, Mx, 1e-4, 4);
	return Mx;
}

}

#endif
