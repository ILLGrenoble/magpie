/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-21
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#define BOOST_TEST_MODULE Fft1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_fft1, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;

	std::vector<t_cplx> invec{{
		1., 2., 3., 4., 5., 6. ,7. ,8 }};

	auto outvec_fft = tl2::fft<t_real, t_cplx, std::vector>(invec, 0, 0);
	auto outvec_dft = tl2::dft<t_real, t_cplx, std::vector>(invec, 0, 0);

	BOOST_TEST((outvec_fft.size() == outvec_dft.size()));
	for(std::size_t i=0; i<std::min(outvec_dft.size(), outvec_fft.size()); ++i)
	{
		const t_cplx& c1 = outvec_fft[i];
		const t_cplx& c2 = outvec_dft[i];
		//std::cout << c1 << "  " << c2 << std::endl;
		BOOST_TEST((tl2::equals<t_cplx>(c1, c2, 1e-4)));
	}
}
