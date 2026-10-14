/**
 * magdyn -- topology test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-oct-2026
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


BOOST_AUTO_TEST_CASE_TEMPLATE(test_magdyn_fm, t_real, t_types_real)
{
	// real type name
	const std::string real_name = ty::type_id_with_cvr<t_real>().pretty_name();

	std::cout << "\n================================================================================" << std::endl;
	std::cout << "Running tests with " << real_name << " data type." << std::endl;
	std::cout << "================================================================================" << std::endl;

	// constants
	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-1 : 1e-4;
	static constexpr t_real diff_eps = std::is_same_v<t_real, float> ? 1e-3 : 1e-5;

	// types
	using t_cplx = std::complex<t_real>;
	using t_mat = tl2::mat<t_cplx>;
	using t_vec = tl2::vec<t_cplx>;

	using t_magdyn = magdyn::MagDyn<
		t_mat, t_vec,
		t_cplx, t_real,
		std::size_t>;

	// magnon calculator
	t_magdyn magdyn{};

	magdyn.Load("../../examples/honeycomb.magpie");
	//magdyn.Load("../../examples/ferromagnetic_chain.magpie");

	magdyn.SetEpsilon(eps);
	magdyn.SetUniteDegenerateEnergies(false);

	std::vector<t_cplx> cs = magdyn.CalcChernNumbers(
	  0.5,             // zone boundary
		diff_eps, 1e-1,  // diff and int deltas
		0, 1,            // dims
		false, true);
	for(t_cplx& c : cs)
	{
		BOOST_TEST(tl2::equals<t_real>(c.imag(), 0., eps));
		BOOST_TEST(tl2::equals<t_real>(std::abs(c.real()), 1., eps));
		tl2::set_eps_0(c, eps);
		std::cout << "chern number: " << c << std::endl;
	}

	std::cout << "================================================================================" << std::endl;
}
