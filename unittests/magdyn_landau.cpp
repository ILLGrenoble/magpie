/**
 * magdyn -- ll test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8-nov-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#define BOOST_TEST_MODULE Magdyn Test

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;
namespace ty = boost::typeindex;

#include "libs/magdyn.h"


using t_types_real = std::tuple<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_magdyn_ll, t_real, t_types_real)
{
	// real type name
	const std::string real_name = ty::type_id_with_cvr<t_real>().pretty_name();

	std::cout << "\n================================================================================" << std::endl;
	std::cout << "Running tests with " << real_name << " data type." << std::endl;
	std::cout << "================================================================================" << std::endl;

	// constants
	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-4 : 1e-6;

	// types
	using t_cplx = std::complex<t_real>;
	using t_mat = tl2::mat<t_cplx>;
	using t_vec = tl2::vec<t_cplx>;
	using t_size = std::size_t;

	using t_magdyn = magdyn::MagDyn<
		t_mat, t_vec,
		t_cplx, t_real,
		t_size>;

	using t_vec3_real = typename t_magdyn::t_vec3_real;
	using SofQE = typename t_magdyn::SofQE;

	// magnon calculator
	t_magdyn magdyn{};

	magdyn.Load("../../examples/ferromagnetic_chain.magpie");

	magdyn.SetEpsilon(eps);
	magdyn.SetUniteDegenerateEnergies(false);

	t_mat Mx = magdyn::calc_mcross_matrix<t_magdyn>(magdyn);

	t_vec3_real Q_start = tl2::create<t_vec3_real>({ -1., 0., 0. });
	t_vec3_real Q_end = tl2::create<t_vec3_real>({ +1., 0., 0. });

	t_size num_Qs = 16;
	for(t_size i = 0; i < num_Qs; ++i)
	{
		// landau-lifshitz method
		t_vec3_real Q = tl2::lerp(Q_start, Q_end, t_real(i) / t_real(num_Qs - 1));
		t_mat F = magdyn::calc_fluctuation_matrix<t_magdyn>(magdyn, Q);
		auto Es_landau = magdyn::calc_energies_landau<t_mat, t_vec, t_cplx, t_real>(F, Mx);

		// holdstein-primakoff method
		SofQE S = magdyn.CalcEnergies(Q, false);

		for(t_size mode = 0; mode < std::min(S.E_and_S.size(), Es_landau.size()); ++mode)
		{
			const t_cplx& E = S.E_and_S[mode].E;
			const t_cplx& E_landau = Es_landau[mode];

			BOOST_TEST(tl2::equals<t_cplx>(E, E_landau, eps));
		}
	}

	std::cout << "================================================================================" << std::endl;
}
