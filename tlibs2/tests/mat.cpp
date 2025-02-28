/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-24
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

#include <iostream>
#include <vector>

#include "libs/maths.h"


int main()
{
	using namespace tl2_ops;

	using t_real = double;
	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;


	{
		auto M = tl2::create<t_mat>({
			1., 2., 3.,
			4., 5., 6.,
			7., 8., 9.
		});

		auto [M_S, M_A] = tl2::split_symm(M);

		std::cout << "symmetric part:\n";
		tl2::niceprint(std::cout, M_S);
		std::cout << std::endl;

		std::cout << "skew-symmetric part:\n";
		tl2::niceprint(std::cout, M_A);
		std::cout << std::endl;
	}


	return 0;
}
