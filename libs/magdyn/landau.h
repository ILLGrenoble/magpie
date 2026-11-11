/**
 * magnetic dynamics -- landau-lifshitz formalism
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8 November 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *	- https://en.wikipedia.org/wiki/Landau%E2%80%93Lifshitz%E2%80%93Gilbert_equation
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), https://kups.ub.uni-koeln.de/7937/
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

#include <limits>
#include <tuple>
#include <algorithm>


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
	using t_real = typename t_magdyn::t_real_;
	using t_size = typename t_magdyn::t_size_;
	using MagneticSite = typename t_magdyn::MagneticSite;
	using t_Jmap = typename t_magdyn::t_Jmap;

	const t_Jmap J_Q0 = mag.CalcReciprocalJs();
	const t_Jmap J_Q = mag.CalcReciprocalJs(Qvec);
	const t_size N = mag.GetMagneticSitesCount();

	constexpr const t_cplx imag(0., 1.);
	constexpr const t_real twopi = t_real(2) * tl2::pi<t_real>;
	constexpr const t_real phase_sign = -1;  // as in magdyn.h

	t_mat F = tl2::zero<t_mat>(N*3, N*3);

	// Q-independent parts
	for(const auto& [sites, J] : J_Q0)
	{
		const t_size idx1 = sites.first;
		const t_size idx2 = sites.second;

		const MagneticSite& s_i = mag.GetMagneticSite(idx1);
		const MagneticSite& s_j = mag.GetMagneticSite(idx2);
		t_cplx JS = tl2::inner_noconj(s_i.trafo_z_calc, J * s_j.trafo_z_calc);
		JS *= s_j.spin_mag_calc / s_i.spin_mag_calc;

		for(int i = 0; i < 3; ++i)
			F(idx1*3 + i, idx1*3 + i) -= JS;
	}

	// phases based on site positions
	std::vector<t_cplx> phases;
	phases.reserve(N);
	for(t_size site_idx = 0; site_idx < N; ++site_idx)
	{
		const MagneticSite& s = mag.GetMagneticSite(site_idx);
		phases.emplace_back(
			std::exp(-phase_sign * imag * twopi *
				tl2::inner<t_vec3_real>(s.pos_calc, Qvec)));
	}

	// Q-dependent parts
	for(const auto& [sites, J] : J_Q)
	{
		const t_size idx1 = sites.first;
		const t_size idx2 = sites.second;

		for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			F(idx1*3 + i, idx2*3 + j) += phases[idx1] * std::conj(phases[idx2]) * J(i, j);
	}

	// field
	for(t_size site_idx = 0; site_idx < N; ++site_idx)
	{
		const MagneticSite& s = mag.GetMagneticSite(site_idx);
		t_cplx field_energy = mag.CalcFieldEnergy(site_idx);
		for(int i = 0; i < 3; ++i)
			F(site_idx*3 + i, site_idx*3 + i) -= field_energy / s.spin_mag_calc;
	}

	//std::cout << "F =\n";
	//tl2::niceprint(std::cout, F, 1e-4, 4);
	return F;
}



/**
 * convert Mx matrix from holstein-primakoff formalism
 */
template<class t_magdyn, class t_mat = typename t_magdyn::t_mat_>
std::tuple<t_mat, t_mat, t_mat> calc_mcross_matrix(const t_magdyn& mag)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
{
	using t_size = typename t_magdyn::t_size_;
	using t_real = typename t_magdyn::t_real_;
	using t_cplx = typename t_magdyn::t_cplx_;
	using t_vec = typename t_magdyn::t_vec_;
	using t_vec3 = typename t_magdyn::t_vec3;
	using t_mat33 = typename t_magdyn::t_mat33;
	using MagneticSite = typename t_magdyn::MagneticSite;

	constexpr const t_cplx imag(0., 1.);
	constexpr const t_real eps = 1e-6;

	const t_size N = mag.GetMagneticSitesCount();
	t_mat Mx = tl2::zero<t_mat>(N*3, N*3);
	t_mat Nx = tl2::zero<t_mat>(N*3, N*3);

	for(t_size site_idx = 0; site_idx < N; ++site_idx)
	{
		const MagneticSite& s = mag.GetMagneticSite(site_idx);
		const t_vec3& v = s.trafo_z_calc;
		t_mat33 matx = tl2::skewsymmetric<t_mat33, t_vec3>(v);

		for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
		{
			Mx(site_idx*3 + i, site_idx*3 + j) = s.spin_mag_calc * matx(i, j);
			Nx(site_idx*3 + i, site_idx*3 + j) = -imag * matx(i, j) / s.spin_mag_calc;
		}
	}

	/*std::cout << "Mx =\n";
	tl2::niceprint(std::cout, Mx, 1e-4, 4);
	std::cout << "Nx =\n";
	tl2::niceprint(std::cout, Nx, 1e-4, 4);*/

	// basis removing longitudinal modes
	auto [evecs_ok, evals, evecs] =
		tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(-imag * Mx,
			false /*only_evals*/, true /*is_herm*/, false /*normalise*/);

	t_real eval_max = std::abs(std::max_element(evals.begin(), evals.end(),
		[](const t_cplx& val1, const t_cplx& val2)
		{
			return std::abs(val1.real()) < std::abs(val2.real());
		})->real());
	std::vector<t_vec> evecs_nonzero;
	for(t_size idx = 0; idx < evecs.size(); ++idx)
	{
		if(std::abs(evals[idx].real()) < eps * eval_max)
			continue;
		evecs_nonzero.push_back(evecs[idx]);
	}

	t_mat basis = tl2::create<t_mat, t_vec>(evecs_nonzero);
	//std::cout << "basis =\n";
	//tl2::niceprint(std::cout, basis, 1e-4, 4);

	return std::make_tuple(Mx, Nx, basis);
}



template<class t_mat, class t_vec,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_cplx::value_type>
std::tuple<std::vector<t_real>, std::vector<t_vec>, std::vector<t_real>>
calc_energies_landau(const t_mat& F, const t_mat& Mx, const t_mat& Nx, const t_mat& basis,
  bool only_energies = true, bool unite_degenerate = true)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
{
	using t_size = typename t_mat::size_type;
	constexpr const t_cplx imag(0., 1.);
	constexpr const t_real eps = 1e-6;

	t_size num_sites = Mx.size1() / 3;

	t_mat mat = -imag * Mx * F;
	t_mat mat_trans = tl2::herm(basis) * mat * basis;
	auto [evecs_ok, evals, evecs] =
		tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(mat_trans,
			only_energies, false /*is_herm*/, false /*normalise*/);
	for(t_vec& evec : evecs)
		evec = basis * evec;

	// sort states by energies in descending order
	const std::vector<t_size> sorting = tl2::get_perm(evals.size(),
		[&evals](t_size idx1, t_size idx2) -> bool
	{
		return evals[idx1].real() > evals[idx2].real();
	});
	evals = tl2::reorder(evals, sorting);
	if(!only_energies)
		evecs = tl2::reorder(evecs, sorting);

	//std::cout << "mat =\n";
	//tl2::niceprint(std::cout, mat, 1e-4, 4);
	/*for(std::size_t i = 0; i < evals.size(); ++i)
	{
		using namespace tl2_ops;
		std::cout << i << ", eval = " << evals[i] << ", evec = " << evecs[i] << std::endl;
	}*/

	std::vector<t_real> energies;
	std::vector<t_vec> states;
	energies.reserve(evals.size());
	states.reserve(evecs.size());

	if(!unite_degenerate)
	{
		for(const t_cplx& eval : evals)
			energies.push_back(eval.real());
		for(t_size idx = 0; idx < evecs.size(); ++idx)
		{
			// symplectic norm
			t_cplx norm = std::sqrt(tl2::inner(evecs[idx], Nx, evecs[idx]));
			if(tl2::equals_0<t_real>(std::abs(norm), eps))
				norm = 1;
			states.emplace_back(evecs[idx] / norm);
		}
	}
	else
	{
		// TODO
	}

	// spin-spin correlations
	std::vector<t_real> weights;
	if(!only_energies)
	{
		weights.reserve(states.size());

		for(t_size state_idx = 0; state_idx < states.size(); ++state_idx)
		{
			// TODO: use the sites' g factors instead
			t_vec state = std::abs(tl2::g_e<t_real>) * states[state_idx];

			t_vec state_sum = tl2::zero<t_vec>(3);
			for(t_size row_beg = 0; row_beg < state.size(); row_beg += 3)
				state_sum += tl2::subvec<t_vec>(state, row_beg, row_beg + 3);

			t_mat S = tl2::outer<t_mat, t_vec>(state_sum, state_sum);
			S /= t_real(num_sites);

			t_real weight = std::abs(tl2::trace(S).real());
			weights.push_back(weight);
		}
	}

	return std::make_tuple(energies, states, weights);
}


}

#endif
