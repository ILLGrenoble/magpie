/**
 * tlibs2 -- magnetic dynamics -- powder calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date july-2026
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

#ifndef __TLIBS2_MAGDYN_POWDER_H__
#define __TLIBS2_MAGDYN_POWDER_H__


#include "magdyn.h"



// --------------------------------------------------------------------
// powder calculation functions
// --------------------------------------------------------------------

/**
 * calculate powder energies for the given Q
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQEs
MAGDYN_INST::CalcPowder(t_real Q_invA,
	t_size num_points, t_size num_threads, bool calc_weights,
	std::function<bool(int, int)> *progress_fkt,
	std::function<void(const MAGDYN_TYPE::SofQE*)> *result_fkt) const
{
	auto [ Binv, inv_ok ] = tl2::inv(m_xtalB);

	if(!inv_ok)
	{
		TL2_CERR_OPT << "Magdyn error: B matrix is not invertible."
			<< std::endl;
	}

	std::vector<t_vec_real> Qvecs_rlu;
	Qvecs_rlu.reserve(num_points);

	// generate mc points
	for(t_size pt = 0; pt < num_points; ++pt)
	{
		// generate a random point on a sphere with radius Q_invA
		t_real u = tl2::get_rand<t_real>(0., 1.);
		t_real v = tl2::get_rand<t_real>(0., 1.);
		const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
		const auto [ x, y, z ] = tl2::sph_to_cart<t_real>(Q_invA, phi, theta);

		// generate Q vector
		t_vec_real Qvec_invA = tl2::create<t_vec_real>({ x, y, z });
		if(m_perform_checks)
		{
			if(!tl2::equals(tl2::norm(Qvec_invA), Q_invA, m_eps))
			{
				TL2_CERR_OPT << "Magdyn error: Invalid Q vector length for |Q| = "
					<< Q_invA << " / A." << std::endl;
			}
		}

		// transform Q vector to rlu
		t_vec_real Qvec_rlu = Binv * Qvec_invA;
		Qvecs_rlu.emplace_back(std::move(Qvec_rlu));
	}

	return CalcDispersion(Qvecs_rlu, num_threads, calc_weights, progress_fkt, result_fkt);
}



/**
 * calculate powder energies for the given Q, binning the energies
 */
MAGDYN_TEMPL
MAGDYN_TYPE::t_histo
MAGDYN_INST::CalcPowderBin(t_real Q_invA,
  t_real E_start, t_real E_end, t_size E_bins,
	t_size num_points, t_size num_threads, bool calc_weights,
	std::function<bool(int, int)> *progress_fkt,
	std::function<void(const MAGDYN_TYPE::SofQE*)> *result_fkt) const
{
	namespace histo = boost::histogram;
	t_histo histE = histo::make_histogram(histo::axis::regular<t_real>(E_bins, E_start, E_end));

	// calculate S(Q, E)
	auto SQEs = CalcPowder(Q_invA, num_points, num_threads, calc_weights,
		progress_fkt, result_fkt);

	// put S(Q, E) into energy bins
	for(const auto& SQE : SQEs)
		for(const auto& E_and_S : SQE.E_and_S)
			histE(E_and_S.E, histo::weight(E_and_S.weight_perp));

	return histE;
}

// --------------------------------------------------------------------


#endif
