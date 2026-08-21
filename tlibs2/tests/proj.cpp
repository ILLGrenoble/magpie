/**
 * orthogonal projector test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 21-aug-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
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

// g++ -std=c++20 -O2 -march=native -I.. -o proj proj.cpp

#include <iostream>

#include "libs/maths.h"


using namespace tl2_ops;

using t_real = double;
using t_vec = tl2::vec_static<t_real, 3>;
using t_mat = tl2::mat_static<t_real, 3, 3>;


/**
 * calculate the spherical average of the projector and the orthogonal projector
 */
void proj_avg()
{
	const std::size_t num_iters = 1e6;
	const t_real Q_len = 0.1;
	const t_real eps = 1e-4;

	t_mat proj_avg = tl2::zero<t_mat>();
	t_mat ortho_proj_avg = tl2::zero<t_mat>();

	for(std::size_t i = 0; i < num_iters; ++i)
	{
		t_real u = tl2::get_rand<t_real>(0., 1.);
		t_real v = tl2::get_rand<t_real>(0., 1.);
		const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(u, v);
		const auto [ h, k, l ] = tl2::sph_to_cart<t_real>(Q_len, phi, theta);

		t_vec Q = tl2::create<t_vec>({ h, k, l });
		t_mat proj = tl2::projector<t_mat, t_vec>(Q, false);
		t_mat ortho_proj = tl2::ortho_projector<t_mat, t_vec>(Q, false);

		proj_avg += proj / t_real(num_iters);
		ortho_proj_avg += ortho_proj / t_real(num_iters);
	}

	std::cout << "Projector spherical average:\n";
	tl2::set_eps_0(proj_avg, eps);
	tl2::niceprint(std::cout, proj_avg);
	std::cout << std::endl;

	std::cout << "Orthogonal projector spherical average:\n";
	tl2::set_eps_0(ortho_proj_avg, eps);
	tl2::niceprint(std::cout, ortho_proj_avg);
	std::cout << std::endl;
}


int main()
{
	proj_avg();

	return 0;
}
