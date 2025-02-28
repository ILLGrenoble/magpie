/**
 * magnetic dynamics test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-june-2024
 * @license GPLv3, see 'LICENSE' file
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

#define BOOST_TEST_MODULE Magdyn Test

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;
namespace ty = boost::typeindex;

#include "libs/magdyn.h"


using t_types_real = std::tuple<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_magdyn, t_real, t_types_real)
{
	// real type name
	const std::string real_name = ty::type_id_with_cvr<t_real>().pretty_name();

	std::cout << "\n================================================================================" << std::endl;
	std::cout << "Running tests with " << real_name << " data type." << std::endl;
	std::cout << "================================================================================" << std::endl;


	static constexpr t_real eps = 1e-4;

	// types
	using t_cplx = std::complex<t_real>;
	using t_mat = tl2::mat<t_cplx>;
	using t_vec = tl2::vec<t_cplx>;
	using t_mat_real = tl2::mat<t_real>;
	using t_vec_real = tl2::vec<t_real>;

	using t_magdyn = tl2_mag::MagDyn<
		t_mat, t_vec,
		t_mat_real, t_vec_real,
		t_cplx, t_real,
		std::size_t>;

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
	auto Es_and_S = magdyn.CalcEnergies(0.1, 0., 0., false);
	BOOST_TEST(Es_and_S.size() == 2);  // + and - energy branch
	BOOST_TEST(tl2::equals<t_real>(Es_and_S[0].E, -Es_and_S[1].E, eps));
	BOOST_TEST(tl2::equals<t_real>(std::abs(Es_and_S[0].E), 0.5878, eps));
	BOOST_TEST(tl2::equals<t_real>(std::abs(Es_and_S[0].weight), 0.6513, eps));


	// calculate and save a dispersion branch
	magdyn.SaveDispersion("disp_" + real_name + ".dat",
		-1. ,0., 0.,  // from
		+1., 0., 0.,  // to
		256);


	// save magnetic model
	magdyn.Save("model_" + real_name + ".magdyn");
	std::cout << "================================================================================" << std::endl;
}
