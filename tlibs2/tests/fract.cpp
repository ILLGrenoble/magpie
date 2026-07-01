/**
 * fractional coordinates test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 1-jul-2026
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
 *
 * clang++ -std=c++20 -I.. -Wall -Wextra -Weffc++ -o fract fract.cpp
 */

#include "libs/maths.h"
using namespace tl2_ops;


int main()
{
	using t_real = double;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	constexpr t_real eps = 1e-6;
	constexpr t_real pi = tl2::pi<t_real>;

	t_mat B = tl2::B_matrix<t_mat>(
		4., 5., 6.,
		90./180.*pi, 90./180.*pi, 120./180.*pi);
	tl2::set_eps_0(B, eps);
	std::cout << "B = \n";
	tl2::niceprint(std::cout, B, eps);
	std::cout << std::endl;

	t_mat G = tl2::metric<t_mat>(B);
	tl2::set_eps_0(G, eps);
	std::cout << "G = \n";
	tl2::niceprint(std::cout, G, eps);
	std::cout << std::endl;

	t_vec vec1 = tl2::create<t_vec>({ 1., -2., 0. });
	t_vec vec2 = tl2::create<t_vec>({ 1., -2., 0. });
	t_real angle = tl2::angle(G, vec1, vec2);
	std::cout << "Angle between " << vec1 << " and " << vec2 << ": ";
	std::cout << tl2::r2d(angle) << " deg." << std::endl;

	return 0;
}
