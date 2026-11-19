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
	static constexpr t_real eps = std::is_same_v<t_real, float> ? 1e-2 : 1e-6;
	static constexpr t_real diff_eps = std::is_same_v<t_real, float> ? 1e-4 : 1e-5;

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

	magdyn.Load("../../examples/honeycomb.magpie");
	//magdyn.Load("../../examples/ferromagnetic_chain.magpie");

	magdyn.SetEpsilon(eps);
	magdyn.SetUniteDegenerateEnergies(false);

	t_vec3_real Q = tl2::create<t_vec3_real>({ 0.65, 0.65, 0. });
	t_vec3_real Q_start = tl2::create<t_vec3_real>({ -0.5, -0.5, 0. });
	t_vec3_real Q_dir1 = tl2::create<t_vec3_real>({ 1., 0., 0. });
	t_vec3_real Q_dir2 = tl2::create<t_vec3_real>({ 0., 1., 0. });


	// berry curvatures -- standard calculation
	std::vector<t_cplx> curvs = std::get<0>(
		magdyn.CalcBerryCurvatures(Q, Q_dir1/64., Q_dir2/64.));
	BOOST_TEST(curvs.size() == 4);
	for(t_cplx& curv : curvs)
	{
		BOOST_TEST(tl2::equals<t_real>(curv.imag(), 0., eps));
		BOOST_TEST(tl2::equals<t_real>(std::abs(curv.real()), std::abs(curvs[0].real()), eps));

		tl2::set_eps_0(curv, eps);
		std::cout << "berry curvature (stable method): " << curv << std::endl;
	}


	// berry curvatures -- alternative calculation
	const t_real delta_int = 2.5e-2;
	t_size dim1 = 0, dim2 = 1;
	std::vector<t_size> *perm = nullptr;
	auto evec_func = magdyn::get_evecmat_func<t_mat, t_vec, t_vec3_real, t_magdyn>(&magdyn, perm);

	std::vector<t_cplx> curvs2 = std::get<0>(
		magdyn::berry_curvatures<
			t_vec, t_vec3_real, t_mat, SofQE, t_cplx, t_real, t_size>(
				evec_func, Q, Q_dir1, Q_dir2,
				diff_eps, delta_int, eps,
				dim1, dim2, true));
	for(t_cplx& curv : curvs2)
	{
		BOOST_TEST(tl2::equals<t_real>(curv.imag(), 0., eps));
		BOOST_TEST(tl2::equals<t_real>(std::abs(curv.real()), std::abs(curvs2[0].real()), eps));

		tl2::set_eps_0(curv, eps);
		std::cout << "berry curvature (unstable method): " << curv << std::endl;
	}


	// chern numbers -- standard calculation
	t_vec cs = magdyn.CalcChernNumbers(
		Q_start, Q_dir1, Q_dir2,  // path
		64, 64);                  // grid

	BOOST_TEST(cs.size() == 4);
	if(cs.size() == 4)
	{
		BOOST_TEST(tl2::equals<t_real>(cs[0].real(), 1., eps));
		BOOST_TEST(tl2::equals<t_real>(cs[1].real(), -1., eps));
		BOOST_TEST(tl2::equals<t_real>(cs[2].real(), 1., eps));
		BOOST_TEST(tl2::equals<t_real>(cs[3].real(), -1., eps));
	}

	for(t_cplx& c : cs)
	{
		BOOST_TEST(tl2::equals<t_real>(c.imag(), 0., eps));
		BOOST_TEST(tl2::equals<t_real>(std::abs(c.real()), 1., eps));

		tl2::set_eps_0(c, eps);
		std::cout << "chern number (stable method): " << c << std::endl;
	}


	// chern number -- alternative calculation
	t_vec cs2 = magdyn::chern_numbers<t_vec, t_vec3_real, t_mat, SofQE,
		t_cplx, t_real, t_size>(
			evec_func, Q_start, Q_dir1, Q_dir2, eps,
			diff_eps, delta_int, dim1, dim2, true);

	for(t_cplx& c : cs2)
	{
		BOOST_TEST(tl2::equals<t_real>(c.imag(), 0., eps));
		BOOST_TEST(tl2::equals<t_real>(std::abs(c.real()), 1., /*eps*/ 5e-2));

		tl2::set_eps_0(c, eps);
		std::cout << "chern number (unstable method): " << c << std::endl;
	}


	std::cout << "================================================================================" << std::endl;
}
