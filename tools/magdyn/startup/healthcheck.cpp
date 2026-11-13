/**
 * magnetic dynamics health check
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-june-2024
 * @license GPLv3, see 'LICENSE' file
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

#include "libs/magdyn.h"
#include "defs.h"
#include "../../../unittests/analytic.h"


bool healthcheck()
{
	// magnon calculator
	t_magdyn magdyn{};


	// add variables
	t_real J = 0.456;
	typename t_magdyn::Variable varJ{};
	varJ.name = "J";
	varJ.value = J;
	magdyn.AddVariable(std::move(varJ));

	t_real S = 2.;
	typename t_magdyn::Variable varS{};
	varS.name = "S";
	varS.value = S;
	magdyn.AddVariable(std::move(varS));


	// add a site
	typename t_magdyn::MagneticSite site{};
	site.name = "site";

	site.pos[0] = "0";
	site.pos[1] = "0";
	site.pos[2] = "0";

	site.spin_dir[0] = "0";
	site.spin_dir[1] = "0";
	site.spin_dir[2] = "1";

	site.spin_mag = "S";

	//magdyn.CalcMagneticSite(site);
	magdyn.AddMagneticSite(std::move(site));


	// add a coupling
	typename t_magdyn::ExchangeTerm coupling{};

	coupling.name = "coupling";

	coupling.site1 = "site";
	coupling.site2 = "site";

	coupling.dist[0] = "1";
	coupling.dist[1] = "0";
	coupling.dist[2] = "0";

	coupling.J = "J";

	//magdyn.CalcExchangeTerm(coupling);
	magdyn.AddExchangeTerm(std::move(coupling));


	// set propagation vector
	t_vec3_real prop = tl2::create<t_vec3_real>({ 0.5, 0., 0. });
	magdyn.SetOrderingWavevector(prop);

	t_vec3_real rotax = tl2::create<t_vec3_real>({ 0., 1., 0. });
	magdyn.SetRotationAxis(rotax);


	// call this after setting the rotation axis to recalculate the Js
	magdyn.CalcMagneticSites();
	magdyn.CalcExchangeTerms();


	using SofQE = typename t_magdyn::SofQE;
	using EnergiesAndWeights = typename t_magdyn::EnergiesAndWeights;


	// compare dispersion and correlation with analytical formula
	auto check_E = [J, S](const SofQE& SQE) -> bool
	{
		const t_real eps = g_eps;

		auto [exact_E, exact_w] = antiferro<t_real>(SQE.Q_rlu[0], J, S);

		const EnergiesAndWeights& Es_and_S = SQE.E_and_S;
		if(Es_and_S.size() != 2)
			return false;
		if(!tl2::equals<t_real>(Es_and_S[0].E, -Es_and_S[1].E, eps))
			return false;
		if(!tl2::equals<t_real>(std::abs(Es_and_S[0].E), exact_E, eps))
			return false;
		if(!tl2::equals<t_real>(std::abs(Es_and_S[0].weight_perp), exact_w, eps))
			return false;

		return true;
	};


	// calculate a single point
	//EnergiesAndWeights Es_and_S = magdyn.CalcEnergies(0.1, 0., 0., false);

	// calculate a dispersion branch
	const t_real h_beg = 0.6;
	const t_real h_end = 0.9;
	const t_size num_Qs = 64;
	const t_size num_threads = 2;

	auto all_Es_and_S = magdyn.CalcDispersion(h_beg, 0., 0.,
		h_end, 0., 0., num_Qs, num_threads, true);

	if(all_Es_and_S.size() != num_Qs)
		return false;
	for(const SofQE& SQE : all_Es_and_S)
	{
		if(!check_E(SQE))
			return false;
	}

	return true;
}
