/**
 * tlibs2 -- magnetic dynamics -- polarisation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - results checked against: https://github.com/SunnySuite/Sunny.jl/blob/main/src/Measurements/MeasureSpec.jl
 *
 * @desc For further references, see the 'LITERATURE' file.
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
bool MAGDYN_INST::CalcPolarisation(const t_vec_real& Q_rlu,
	MAGDYN_TYPE::EnergyAndWeight& E_and_S) const
{
	if(tl2::equals_0<t_vec_real>(Q_rlu, m_eps))
	{
		TL2_CERR_OPT << "Magdyn error: Q coordinate system cannot be calculated for Q = 0."
			<< std::endl;
		return false;
	}


	/*t_vec_real Q_lab = m_xtalUB * Q_rlu;
	t_vec_real h_lab = m_xtalUB * tl2::create<t_vec_real>({ 1., 0., 0. });
	t_vec_real l_lab = m_xtalUB * tl2::create<t_vec_real>({ 0., 0., 1. });
	t_mat_real rotQ = tl2::rotation<t_mat_real>(h_lab, Q_lab, &l_lab, m_eps, true);
	t_mat_real rotQ_hkl = m_xtalUBinv * rotQ * m_xtalUB;*/


	t_vec_real Q_lab = m_xtalB * Q_rlu;
	Q_lab /= tl2::norm(Q_lab);

	t_vec_real up_lab = m_xtalB * m_scatteringplane[2];
	up_lab /= tl2::norm(up_lab);

	t_vec_real Qperp_lab = tl2::cross(up_lab, Q_lab);
	Qperp_lab /= norm(Qperp_lab);

	t_mat_real rotQ = tl2::create<t_mat_real>(3, 3);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ, Q_lab, 0);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ, Qperp_lab, 1);
	tl2::set_row<t_mat_real, t_vec_real>(rotQ, up_lab, 2);


	//t_mat_real rotQ_inv = tl2::trans<t_mat_real>(rotQ);
	const auto [rotQ_inv, rotQ_inv_ok] = tl2::inv(rotQ);
	if(!rotQ_inv_ok)
	{
		using namespace tl2_ops;
		TL2_CERR_OPT << "Magdyn error: Cannot invert Q rotation matrix for Q = "
			<< Q_rlu << "." << std::endl;
		return false;
	}

	t_mat rotQ_cplx = tl2::convert<t_mat, t_mat_real>(rotQ);
	t_mat rotQ_inv_cplx = tl2::convert<t_mat, t_mat_real>(rotQ_inv);

	E_and_S.S_pol_perp = rotQ_cplx * E_and_S.S_perp * rotQ_inv_cplx;
	E_and_S.S_pol = rotQ_cplx * E_and_S.S * rotQ_inv_cplx;


	/*using namespace tl2_ops;
	std::cout << "Q_rlu =\n";
	tl2::niceprint(std::cout, Q_rlu, 1e-4, 4);
	std::cout << "Q_lab =\n";
	tl2::niceprint(std::cout, Q_lab, 1e-4, 4);
	std::cout << "rotQ =\n";
	tl2::niceprint(std::cout, rotQ, 1e-4, 4);
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

	return true;
}

// --------------------------------------------------------------------


#endif
