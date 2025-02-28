/**
 * calculates the berry curvatures
 * @author Tobias Weber <tweber@ill.fr>
 * @date 5-November-2024
 * @license see 'LICENSE' file
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

// g++ -std=c++20 -march=native -O2 -Wall -Wextra -Weffc++ -I .. -o berry berry.cpp -llapacke -llapack -lblas -lgfortran


#include "libs/magdyn.h"
using namespace tl2_ops;


#define CALC_CHERN_NUM 0


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
static constexpr t_real print_eps = 1e-5;
static constexpr unsigned int print_prec = 5;



void print_states(const t_SofQE& S)
{
	std::cout << "\nQ = " << S.Q_rlu << ", E = ";
	for(const auto& EandS : S.E_and_S)
		std::cout << EandS.E << ", ";

	std::cout << "states = \n";
	tl2::niceprint(std::cout, S.evec_mat, print_eps, print_prec);
	std::cout << std::endl;
}



int main(int argc, char** argv)
{
	t_real delta = eps;       // for differentiation
	t_real delta_int = 1e-3;  // for integration
	t_vec_real Q0 = tl2::create<t_vec_real>({ 0., 0., 0. });
	t_real max_curv = 100.;

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

#if CALC_CHERN_NUM != 0
	std::cout << "# Chern numbers: ";
	std::vector<t_cplx> cherns = magdyn.CalcChernNumbers(0.5, delta, delta_int);
	for(const t_cplx& chern : cherns)
		std::cout << chern << " ";
	std::cout << std::endl;
#endif

	std::cout << std::left << std::setw(print_width) << "# q" << " ";
	std::cout << std::left << std::setw(print_width) << "E_1" << " ";
	std::cout << std::left << std::setw(print_width) << "Re(b_1)" << " ";
	std::cout << std::left << std::setw(print_width) << "Im(b_1)" << " ";
	std::cout << std::left << std::setw(print_width) << "...";
	std::cout << std::endl;

	t_vec_real Q = Q0;
	for(t_real q = 0.; q < 1.; q += 0.001)
	{
		Q[0] = Q0[0] + q;
		std::cout << std::left << std::setw(print_width) << q << " ";

		t_SofQE S; // = magdyn.CalcEnergies(Q, false);
		//print_states(S);

		// band permutations
		/*t_size num_bands = S.E_and_S.size();
		std::vector<t_size> perm;
		perm.resize(num_bands);
		std::iota(perm.begin(), perm.end(), 0);
		if(q > 0.5)
		{
			std::swap(perm[0], perm[1]);
			std::swap(perm[2], perm[3]);
		}*/

		std::vector<t_cplx> curves;
		std::tie(curves, S) = magdyn.CalcBerryCurvatures(Q, delta/*, &perm*/);

		for(t_size band = 0; band < curves.size(); ++band)
		{
			const t_cplx& curve = curves[band];
			//t_real E = S.E_and_S[perm[band]].E;
			t_real E = S.E_and_S[band].E;
			std::cout << std::left << std::setw(print_width) << E << " ";

			if(std::abs(curve.real()) > max_curv)
				std::cout << std::left << std::setw(print_width) << "-" << " ";
			else
				std::cout << std::left << std::setw(print_width) << curve.real() << " ";

			if(std::abs(curve.imag()) > max_curv)
				std::cout << std::left << std::setw(print_width) << "-" << " ";
			else
				std::cout << std::left << std::setw(print_width) << curve.imag() << " ";
		}

		std::cout << std::endl;
	}

	return 0;
}
