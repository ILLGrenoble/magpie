/**
 * special functions test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 25-nov-2026
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

// g++ -std=c++20 -O2 -march=native -I.. -o funcs funcs.cpp

#include <iostream>

#include "libs/maths.h"


using namespace tl2_ops;

using t_real = double;
using t_cplx = std::complex<t_real>;


void polylog()
{
	t_cplx z = -1;
	t_cplx dz = 0.01;

	for(std::size_t i = 0; i < 200; ++i)
	{
		t_cplx val = tl2::polylog<t_cplx>(z);
		std::cout << z.real() << " " << val.real() << " " << val.imag() << std::endl;

		z += dz;
	}
}


int main()
{
	polylog();

	return 0;
}
