/**
 * calculates the nonreciprocal energy difference
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-July-2025
 * @license see 'LICENSE' file
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

// g++ -std=c++20 -march=native -O2 -Wall -Wextra -Weffc++ -I .. -o nonrecip nonrecip.cpp -llapacke -llapack -lblas -lgfortran


#include "libs/magdyn.h"
using namespace tl2_ops;


// types
using t_real = double;
using t_cplx = std::complex<t_real>;
using t_mat = tl2::mat<t_cplx>;
using t_vec = tl2::vec<t_cplx>;
using t_mat_real = tl2::mat<t_real>;
using t_vec_real = tl2::vec<t_real>;
using t_size = std::size_t;
using t_magdyn = tl2_mag::MagDyn<
	t_mat, t_vec, t_mat_real, t_vec_real,
	t_cplx, t_real, t_size>;
using t_SofQE = typename t_magdyn::SofQE;
using t_field = typename t_magdyn::ExternalField;


static constexpr t_real eps = 1e-12;
static constexpr unsigned int print_prec = 5;



int main(int argc, char** argv)
{
	t_vec_real Q0 = tl2::create<t_vec_real>({ 0., 0., 0. });
	t_real q_eps = 0.001;

	unsigned int print_width = print_prec * 3;
	std::cout.precision(print_prec);

	if(argc < 2)
	{
		std::cerr << "Please specify a magdyn file." << std::endl;
		return -1;
	}

	/*t_field field001
	{
		.align_spins = true,
		.dir = tl2::create<t_vec_real>({ 0., 0., 1. }),
		.mag = 1.,
	};*/

	t_magdyn magdyn{};

	if(!magdyn.Load(argv[1]))
	{
		std::cerr << "Could not load model." << std::endl;
		return -1;
	}

	magdyn.SetEpsilon(eps);
	magdyn.SetUniteDegenerateEnergies(false);
	//magdyn.SetExternalField(field001);
	//magdyn.CalcExternalField();
	//magdyn.CalcMagneticSites();
	//magdyn.CalcExchangeTerms();


	std::cout << std::left << std::setw(print_width) << "# q" << " ";
	std::cout << std::left << std::setw(print_width) << "Ep_1" << " ";
	std::cout << std::left << std::setw(print_width) << "Em_1" << " ";
	std::cout << std::left << std::setw(print_width) << "dE_1" << " ";
	std::cout << std::left << std::setw(print_width) << "...";
	std::cout << std::endl;

	t_vec_real Q = Q0;
	for(t_real q = 0.; q < 1.; q += q_eps)
	{
		Q[0] = Q0[0] + q;
		std::cout << std::left << std::setw(print_width) << q << " ";

		t_SofQE S = magdyn.CalcEnergies(Q, true);
		t_size num_bands = S.E_and_S.size();

		for(t_size band = 0; band < num_bands / 2; ++band)
		{
			t_real Eneg = S.E_and_S[band].E;
			t_real Epos = S.E_and_S[num_bands - band - 1].E;

			std::cout
				<< std::left << std::setw(print_width) << Epos << " "
				<< std::left << std::setw(print_width) << Eneg << " "
				<< std::left << std::setw(print_width) << Epos + Eneg << " ";
		}

		std::cout << std::endl;
	}

	return 0;
}
