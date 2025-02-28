/**
 * tlibs2 -- magnetic dynamics -- polarisation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - TODO
 *
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

#ifndef __TLIBS2_MAGDYN_POLARISATION_H__
#define __TLIBS2_MAGDYN_POLARISATION_H__

#include <vector>
#include <iostream>
#include <iomanip>

#include "../maths.h"
#include "../phys.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * TODO: implements the blume-maleev formalism
 */
MAGDYN_TEMPL
void MAGDYN_INST::CalcPolarisation(const t_vec_real& /*Q_rlu*/,
	MAGDYN_TYPE::EnergyAndWeight& /*E_and_S*/) const
{
	/*  TODO: polarisation via blume-maleev equation
	...
	...
	static const t_vec_real h_rlu = tl2::create<t_vec_real>({ 1., 0., 0. });
	static const t_vec_real l_rlu = tl2::create<t_vec_real>({ 0., 0., 1. });
	t_vec_real h_lab = m_xtalUB * h_rlu;
	t_vec_real l_lab = m_xtalUB * l_rlu;
	t_vec_real Q_lab = m_xtalUB * Q_rlu;
	t_mat_real rotQ = tl2::rotation<t_mat_real>(h_lab, Q_lab, &l_lab, m_eps, true);
	t_mat_real rotQ_hkl = m_xtalUBinv * rotQ * m_xtalUB;
	const auto [rotQ_hkl_inv, rotQ_hkl_inv_ok] = tl2::inv(rotQ_hkl);
	if(!rotQ_hkl_inv_ok)
		std::cerr << "Magdyn error: Cannot invert Q rotation matrix."
			<< std::endl;

	t_mat rotQ_hkl_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl);
	t_mat rotQ_hkl_inv_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl_inv);

	t_mat S_perp_pol = rotQ_hkl_cplx * E_and_S.S_perp * rotQ_hkl_inv_cplx;
	t_mat S_pol = rotQ_hkl_cplx * E_and_S.S * rotQ_hkl_inv_cplx;
	...
	...
	*/


	// test for helix
	/*for(std::uint8_t i = 0; i < 3; ++i)
	{
		t_mat pol = get_polarisation_incommensurate<t_mat>(i, false);
		t_mat S = pol * E_and_S.S;
		t_mat S_perp = pol * E_and_S.S_perp;

		// TODO: set the other components to 0
		E_and_S.S(i, i) = std::abs(tl2::trace<t_mat>(S).real());
		E_and_S.S_perp(i, i) = std::abs(tl2::trace<t_mat>(S_perp).real());
	}*/
}

// --------------------------------------------------------------------

#endif
