/**
 * magnetic dynamics -- landau-lifshitz formalism
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
 * convert fluctuation matrix from holstein-primakoff formalism
 */
template<class t_magdyn,
	class t_mat = typename t_magdyn::t_mat_,
	class t_vec3_real = typename t_magdyn::t_vec3_real>
t_mat calc_fluctuation_matrix(const t_magdyn &mag, const t_vec3_real &Qvec)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec3_real>
#endif
{
	using t_cplx = typename t_magdyn::t_cplx_;
	using t_size = typename t_magdyn::t_size_;
	using t_vec3 = typename t_magdyn::t_vec3;
	using MagneticSite = typename t_magdyn::MagneticSite;
	using t_Jmap = typename t_magdyn::t_Jmap;

	const t_Jmap J_Q0 = mag.CalcReciprocalJs();
	const t_Jmap J_Q = mag.CalcReciprocalJs(Qvec);
	const t_size N = mag.GetMagneticSitesCount();

	t_mat F = tl2::zero<t_mat>(N*3, N*3);

	// Q-independent parts
	for(const auto& [sites, J] : J_Q0)
	{
		const t_size idx1 = sites.first;
		const t_size idx2 = sites.second;

		const MagneticSite& s_i = mag.GetMagneticSite(idx1);
		const MagneticSite& s_j = mag.GetMagneticSite(idx2);
		const t_vec3& v_i = s_i.trafo_z_calc;
		const t_vec3& v_j = s_j.trafo_z_calc;
		t_cplx JS = -tl2::inner_noconj(v_i, J * v_j);

		t_cplx field_energy = -mag.CalcFieldEnergy(idx1);

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
 * convert Mx matrix from holstein-primakoff formalism
 */
template<class t_magdyn, class t_mat = typename t_magdyn::t_mat_>
t_mat calc_mcross_matrix(const t_magdyn& mag)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
{
	using t_size = typename t_magdyn::t_size_;
	using t_vec3 = typename t_magdyn::t_vec3;
	using t_mat33 = typename t_magdyn::t_mat33;
	using MagneticSite = typename t_magdyn::MagneticSite;

	const t_size N = mag.GetMagneticSitesCount();
	t_mat Mx = tl2::zero<t_mat>(N*3, N*3);

	for(t_size site_idx = 0; site_idx < N; ++site_idx)
	{
		const MagneticSite& s = mag.GetMagneticSite(site_idx);
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



template<class t_mat, class t_vec,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_cplx::value_type>
std::vector<t_cplx> calc_energies_landau(const t_mat& F, const t_mat& Mx)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
{
	t_mat mat = -t_cplx(0, 1) * Mx * F;
	bool is_herm = false;
	//std::cout << "mat =\n";
	//tl2::niceprint(std::cout, mat, 1e-4, 4);

	const auto [evecs_ok, evals, evecs] =
		tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
			mat, true, is_herm, true);

	/*for(const t_cplx& eval : evals)
		std::cout << eval << ", ";
	std::cout << std::endl;*/
	return evals;
}


}

#endif
