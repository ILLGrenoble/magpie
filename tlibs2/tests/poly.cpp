/**
 * 2d polygon test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 1-nov-2026
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

// g++ -std=c++20 -O2 -march=native -I.. -o poly poly.cpp

#include <iostream>

#include "libs/maths.h"


using namespace tl2_ops;

using t_int = int;
using t_vec = tl2::vec_static<t_int, 2>;


int main()
{
	std::vector<t_vec> poly
	{
		tl2::create<t_vec>({ 10, 1 }),
		tl2::create<t_vec>({ 30, 1 }),
		tl2::create<t_vec>({ 20, 8 }),
	};

	for(int y = 0; y < 10; ++y)
	{
		for(int x = 0; x < 50; ++x)
		{
			t_vec pt = tl2::create<t_vec>({ x, y });

			bool inside = tl2::in_poly_2d(pt, poly);
			std::cout << (inside ? "*" : " ");
		}

		std::cout << std::endl;
	}

	return 0;
}
