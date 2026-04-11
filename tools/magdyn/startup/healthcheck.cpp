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

#include "tlibs2/libs/magdyn.h"
#include "defs.h"


bool healthcheck()
{
	// magnon calculator
	t_magdyn magdyn{};


	// add a variable
	typename t_magdyn::Variable var{};
	var.name = "J";
	var.value = 0.5;

	magdyn.AddVariable(std::move(var));


	// add a site
	typename t_magdyn::MagneticSite site{};
	site.name = "site";

	site.pos[0] = "0";
	site.pos[1] = "0";
	site.pos[2] = "0";

	site.spin_dir[0] = "0";
	site.spin_dir[1] = "0";
	site.spin_dir[2] = "1";

	site.spin_mag = "1";

	magdyn.CalcMagneticSite(site);
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

	magdyn.CalcExchangeTerm(coupling);
	magdyn.AddExchangeTerm(std::move(coupling));


	// set propagation vector
	t_vec_real prop = tl2::create<t_vec_real>({ 0.5, 0., 0. });
	magdyn.SetOrderingWavevector(prop);

	t_vec_real rotax = tl2::create<t_vec_real>({ 0., 1., 0. });
	magdyn.SetRotationAxis(rotax);


	// calculate a point on the dispersion
	const t_real h = 0.1;
	const t_real w = 0.5878;
	const t_real w_perp = 0.6513;

	using EnergiesAndWeights = typename t_magdyn::EnergiesAndWeights;
	EnergiesAndWeights Es_and_S = magdyn.CalcEnergies(h, 0., 0., false);

	auto check_E = [](const EnergiesAndWeights& Es_and_S, t_real w, t_real w_perp) -> bool
	{
		const t_real eps = 1e-4; //g_eps;

		if(Es_and_S.size() != 2)
			return false;
		if(!tl2::equals<t_real>(Es_and_S[0].E, -Es_and_S[1].E, eps))
			return false;
		if(!tl2::equals<t_real>(std::abs(Es_and_S[0].E), w, eps))
			return false;
		if(!tl2::equals<t_real>(std::abs(Es_and_S[0].weight_perp), w_perp, eps))
			return false;
		return true;
	};

	if(!check_E(Es_and_S, w, w_perp))
		return false;


	// calculate a dispersion branch
	const t_real h_end = 0.9;
	const t_size num_Qs = 64;
	const t_size num_threads = 2;

	auto all_Es_and_S = magdyn.CalcDispersion(h, 0., 0.,
		h_end, 0., 0., num_Qs, num_threads, true);

	if(all_Es_and_S.size() != num_Qs)
		return false;
	if(!check_E(all_Es_and_S[0].E_and_S, w, w_perp))
		return false;

	return true;
}
