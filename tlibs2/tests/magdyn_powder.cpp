/**
 * magnetic dynamics test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-july-2026
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

// g++ -std=c++20 -I.. -I /opt/homebrew/Cellar/boost/1.90.0_1/include -L/opt/homebrew/Cellar/gcc/16.1.0/lib/gcc/16 -o magdyn_powder magdyn_powder.cpp -lgfortran -llapacke -llapack -lblas

#include "libs/magdyn.h"


// types
using t_real = double;
using t_cplx = std::complex<t_real>;
using t_size = std::size_t;

using t_mat = tl2::mat<t_cplx>;
using t_vec = tl2::vec<t_cplx>;
using t_mat_real = tl2::mat<t_real>;
using t_vec_real = tl2::vec<t_real>;

using t_magdyn = tl2_mag::MagDyn<
	t_mat, t_vec,
	t_mat_real, t_vec_real,
	t_cplx, t_real,
	std::size_t>;


int main()
{
	// magnon calculator
	t_magdyn magdyn{};
	magdyn.SetPerformChecks(true);
	magdyn.SetSilent(false);


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


	// calculate powder spectra
	t_size num_E = 64;
	t_size num_Q = 16384;
	auto E_histo = magdyn.CalcPowderBin(0.2, 0., 1., num_E, num_Q);

	int field_w = 15;
	std::cout << "# " << std::setw(field_w - 2) << "E" << std::setw(field_w) << "S" << std::endl;
	for(const auto& idx : boost::histogram::indexed(E_histo))
	{
		const auto& bin = idx.bin();
		t_real E = bin.lower() + (bin.upper() - bin.lower()) / 2.;
		t_real S = *idx / t_real(num_Q) / t_real(num_E);

		std::cout << std::setw(field_w) << E << std::setw(field_w) << S << std::endl;
	}


	// calculate a point on the dispersion
	//auto Es_and_S = magdyn.CalcEnergies(0.1, 0., 0., false);

	// calculate and save a dispersion branch
	//magdyn.SaveDispersion("disp_" + real_name + ".dat",
	//	-1. ,0., 0.,  // from
	//	+1., 0., 0.,  // to
	//	256);


	// save magnetic model
	//magdyn.Save("model_" + real_name + ".magdyn");
}
