/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS2_MAGDYN_GS_H__
#define __TLIBS2_MAGDYN_GS_H__

#include <vector>
#include <unordered_set>
#include <string>
#include <iostream>
#include <iomanip>

#include "../maths.h"
#include "../str.h"
#include "../fit.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * get the energy minimum
 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
 */
MAGDYN_TEMPL
t_real MAGDYN_INST::CalcMinimumEnergy() const
{
	// energies at (000)
	const auto E_and_S = CalcEnergies(0., 0., 0., true);

	// get minimum
	const auto min_iter = std::min_element(
		E_and_S.begin(), E_and_S.end(),
		[](const EnergyAndWeight& E_and_S_1, const EnergyAndWeight& E_and_S_2) -> bool
	{
		return std::abs(E_and_S_1.E) < std::abs(E_and_S_2.E);
	});

	if(min_iter == E_and_S.end())
		return 0.;

	return min_iter->E;
}



/**
 * get the ground-state energy
 * @note zero-operator term in expansion of equation (20) in (Toth 2015)
 */
MAGDYN_TEMPL
t_real MAGDYN_INST::CalcGroundStateEnergy() const
{
	t_real E = 0.;

	for(const ExchangeTerm& term : GetExchangeTerms())
	{
		// check if the site indices are valid
		if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
			continue;

		const MagneticSite& s_i = GetMagneticSite(term.site1_calc);
		const MagneticSite& s_j = GetMagneticSite(term.site2_calc);

		const t_mat J = CalcRealJ(term);  // Q = 0 -> no rotation needed

		const t_vec Si = s_i.spin_mag_calc * s_i.trafo_z_calc;
		const t_vec Sj = s_j.spin_mag_calc * s_j.trafo_z_calc;

		E += tl2::inner_noconj<t_vec>(Si, J * Sj).real();
	}

	return E;
}



/**
 * minimise energy to found ground state
 */
#if defined(__TLIBS2_USE_MINUIT__) && defined(__TLIBS2_MAGDYN_USE_MINUIT__)
MAGDYN_TEMPL
bool MAGDYN_INST::CalcGroundState(const std::unordered_set<std::string>* fixed_params,
	bool verbose, const bool *stop_request)
{
	// function to minimise the state's energy
	auto func = [this](const std::vector<tl2::t_real_min>& args)
	{
		// set new spin configuration
		auto dyn = *this;

		for(t_size site_idx = 0; site_idx < dyn.GetMagneticSitesCount(); ++site_idx)
		{
			MagneticSite& site = dyn.m_sites[site_idx];

			t_real u = args[site_idx * 2 + 0];
			t_real v = args[site_idx * 2 + 1];

			const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
			const auto [ x, y, z ] = tl2::sph_to_cart<t_real>(1., phi, theta);

			site.spin_dir[0] = tl2::var_to_str(x, m_prec);
			site.spin_dir[1] = tl2::var_to_str(y, m_prec);
			site.spin_dir[2] = tl2::var_to_str(z, m_prec);

			dyn.CalcMagneticSite(site);
		}

		// ground state energy with the new configuration
		return dyn.CalcGroundStateEnergy();
	};

	// add minimisation parameters and initial values
	t_size num_args = GetMagneticSitesCount() * 2;
	std::vector<std::string> params;
	std::vector<t_real> vals, errs;
	std::vector<t_real> lower_lims, upper_lims;
	std::vector<bool> fixed;
	params.reserve(num_args);
	vals.reserve(num_args);
	errs.reserve(num_args);
	lower_lims.reserve(num_args);
	upper_lims.reserve(num_args);
	fixed.reserve(num_args);

	for(const MagneticSite& site : GetMagneticSites())
	{
		const t_vec_real& S = site.spin_dir_calc;
		const auto [ rho, phi, theta ] =  tl2::cart_to_sph<t_real>(S[0], S[1], S[2]);
		const auto [ u, v ] = tl2::sph_to_uv<t_real>(phi, theta);

		std::string phi_name = site.name + "_phi";     // phi or u parameter
		std::string theta_name = site.name + "_theta"; // theta or v parameter
		params.push_back(phi_name);
		params.push_back(theta_name);

		bool fix_phi = false;
		bool fix_theta = false;
		if(fixed_params)
		{
			fix_phi = (fixed_params->find(phi_name) != fixed_params->end());
			fix_theta = (fixed_params->find(theta_name) != fixed_params->end());
		}
		fixed.push_back(fix_phi);
		fixed.push_back(fix_theta);

		vals.push_back(u);
		vals.push_back(v);

		lower_lims.push_back(0. -m_eps);
		lower_lims.push_back(0. -m_eps);

		upper_lims.push_back(1. + m_eps);
		upper_lims.push_back(1. + m_eps);

		errs.push_back(0.1);
		errs.push_back(0.1);
	}

	if(tl2::minimise_dynargs<t_real>(num_args, func,
		params, vals, errs, &fixed, &lower_lims, &upper_lims,
		verbose && !m_silent, stop_request))
	{
		// set the spins to the newly-found ground state
		for(t_size site_idx = 0; site_idx < GetMagneticSitesCount(); ++site_idx)
		{
			MagneticSite& site = m_sites[site_idx];

			t_real u = vals[site_idx * 2 + 0];
			t_real v = vals[site_idx * 2 + 1];
			tl2::set_eps_round<t_real>(u, m_eps);
			tl2::set_eps_round<t_real>(v, m_eps);

			const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
			auto [ x, y, z ] = tl2::sph_to_cart<t_real>(1., phi, theta);
			tl2::set_eps_round<t_real>(x, m_eps);
			tl2::set_eps_round<t_real>(y, m_eps);
			tl2::set_eps_round<t_real>(z, m_eps);

			site.spin_dir[0] = tl2::var_to_str(x, m_prec);
			site.spin_dir[1] = tl2::var_to_str(y, m_prec);
			site.spin_dir[2] = tl2::var_to_str(z, m_prec);

			CalcMagneticSite(site);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			using namespace tl2_ops;
			std::cout << site.name
				<< ": u = " << u << ", "
				<< "v = " << v << ", "
				<< "phi = " << phi << ", "
				<< "theta = " << theta << ", "
				<< "S = " << site.spin_dir_calc << std::endl;
#endif
		}

		return true;
	}
	else
	{
		CERR_OPT << "Magdyn error: Ground state minimisation did not converge."
			<< std::endl;

		return false;
	}
}
#else  // __TLIBS2_USE_MINUIT__
MAGDYN_TEMPL
bool MAGDYN_INST::CalcGroundState(const std::unordered_set<std::string>*, bool, const bool*)
{
	std::cerr << "Magdyn error: Ground state minimisation support disabled." << std::endl;
	return false;
}
#endif

// --------------------------------------------------------------------

#endif
