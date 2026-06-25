/**
 * tlibs2 -- magnetic dynamics -- polarisation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - TODO
 *
 * @desc For further references, see the 'LITERATURE' file.
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

#ifndef __TLIBS2_MAGDYN_POLARISATION_H__
#define __TLIBS2_MAGDYN_POLARISATION_H__


#include "../phys.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * implements the blume-maleev formalism
 */
MAGDYN_TEMPL
void MAGDYN_INST::CalcPolarisation(const t_vec_real& Q_rlu,
	MAGDYN_TYPE::EnergyAndWeight& E_and_S) const
{
	/*const t_mat_real& UB = m_xtalUB;
	const t_mat_real& UBinv = m_xtalUBinv;

	t_vec_real Q_lab = UB * Q_rlu;
	t_vec_real h_lab = UB * tl2::create<t_vec_real>({ 1., 0., 0. });
	t_vec_real l_lab = UB * tl2::create<t_vec_real>({ 0., 0., 1. });
	t_mat_real rotQ = tl2::rotation<t_mat_real>(h_lab, Q_lab, &l_lab, m_eps, true);
	t_mat_real rotQ_hkl = UBinv * rotQ * UB;*/


	t_vec_real Q_lab = m_xtalB * Q_rlu;
	Q_lab /= tl2::norm(Q_lab);

	t_vec_real up_lab = m_xtalB * m_scatteringplane[2];
	up_lab /= tl2::norm(up_lab);

	t_vec_real Qperp_lab = tl2::cross(up_lab, Q_lab);
	Qperp_lab /= norm(Qperp_lab);

	t_mat_real rotQ_hkl = tl2::create<t_mat_real>(3, 3);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ_hkl, Q_lab, 0);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ_hkl, Qperp_lab, 1);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ_hkl, up_lab, 2);


	//t_mat_real rotQ_hkl_inv = tl2::trans<t_mat_real>(rotQ_hkl);
	const auto [rotQ_hkl_inv, rotQ_hkl_inv_ok] = tl2::inv(rotQ_hkl);
	if(!rotQ_hkl_inv_ok)
	{
		TL2_CERR_OPT << "Magdyn error: Cannot invert Q rotation matrix."
			<< std::endl;
	}

	t_mat rotQ_hkl_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl);
	t_mat rotQ_hkl_inv_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_hkl_inv);

	E_and_S.S_pol_perp = rotQ_hkl_cplx * E_and_S.S_perp * rotQ_hkl_inv_cplx;
	E_and_S.S_pol = rotQ_hkl_cplx * E_and_S.S * rotQ_hkl_inv_cplx;


	/*using namespace tl2_ops;
	std::cout << "UB =\n";
	tl2::niceprint(std::cout, UB, 1e-4, 4);
	std::cout << "Q_rlu =\n";
	tl2::niceprint(std::cout, Q_rlu, 1e-4, 4);
	std::cout << "Q_lab =\n";
	tl2::niceprint(std::cout, Q_lab, 1e-4, 4);
	std::cout << "rotQ =\n";
	tl2::niceprint(std::cout, rotQ, 1e-4, 4);
	std::cout << "rotQ_hkl =\n";
	tl2::niceprint(std::cout, rotQ_hkl, 1e-4, 4);
	std::cout << std::endl;*/


	// test for helix
	/*for(std::uint8_t i = 0; i < 3; ++i)
	{
		t_mat pol = get_polarisation_incommensurate<t_mat>(i, false);
		t_mat S_pol = pol * E_and_S.S_pol;
		t_mat S_pol_perp = pol * E_and_S.S_pol_perp;

		// TODO: set the other components to 0
		E_and_S.S_pol(i, i) = std::abs(tl2::trace<t_mat>(S_pol).real());
		E_and_S.S_pol_perp(i, i) = std::abs(tl2::trace<t_mat>(S_pol_perp).real());
	}*/
}

// --------------------------------------------------------------------


#endif
