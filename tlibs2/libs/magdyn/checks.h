/**
 * tlibs2 -- magnetic dynamics -- sanity checks
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
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

#ifndef __TLIBS2_MAGDYN_CHECKS_H__
#define __TLIBS2_MAGDYN_CHECKS_H__


#include "magdyn.h"


/**
 * check if the site is valid
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckMagneticSite(t_size idx, bool print_error) const
{
	// always perform this check as the GUI deliberately sets invalid
	// indices when adding a new coupling to the table
	//if(!m_perform_checks)
	//	return true;

	const t_size N = m_sites.size();
	if(idx >= N)
	{
		if(print_error)
		{
			TL2_CERR_OPT << "Magdyn error: Site index " << idx
				<< " is out of bounds. Number of sites: " << N
				<< "." << std::endl;
		}

		return false;
	}
	
	if(!m_perform_checks)
		return true;

	const MagneticSite& site = m_sites[idx];
	for(t_size comp = 0; comp < site.pos_calc.size(); ++comp)
	{
		// position outside unit cell extents?
		if(site.pos_calc[comp] < m_uc_min - m_eps ||
			site.pos_calc[comp] > m_uc_max + m_eps)
		{
			if(print_error)
			{
				using namespace tl2_ops;
				TL2_CERR_OPT << "Magdyn error: Site " << idx
					<< " position is out of the unit cell extents."
					<< " Position: " << site.pos_calc
					<< ", unit cell: [" << m_uc_min << ", " << m_uc_max << "]"
					<< "." << std::endl;
			}
			break;
			//return false;  // TODO
		}
	}

	return true;
}



/**
 * check if the term is valid
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckExchangeTerm(t_size idx, bool print_error) const
{
	if(!m_perform_checks)
		return true;

	const t_size N = m_exchange_terms.size();
	if(idx >= N)
	{
		if(print_error)
		{
			TL2_CERR_OPT << "Magdyn error: Coupling index " << idx
				<< " is out of bounds. Number of couplings: " << N
				<< "." << std::endl;
		}

		return false;
	}

	return true;
}



/**
 * check if imaginary weights remain
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckImagWeights(const MAGDYN_TYPE::SofQE& S) const
{
	if(!m_perform_checks)
		return true;

	using namespace tl2_ops;
	bool ok = true;

	for(const EnergyAndWeight& EandS : S.E_and_S)
	{
		// imaginary parts should be gone after UniteEnergies()
		if(!tl2::equals_0(EandS.S_perp_sum.imag(), m_eps) ||
			!tl2::equals_0(EandS.S_sum.imag(), m_eps))
		{
			ok = false;

			TL2_CERR_OPT << "Magdyn warning: Remaining imaginary S(Q, E) component at Q = "
				<< S.Q_rlu << " and E = " << EandS.E
				<< ": imag(S) = " << EandS.S_sum.imag()
				<< ", imag(S_perp) = " << EandS.S_perp_sum.imag()
				<< "." << std::endl;
		}
	}

	return ok;
}

#endif
