/**
 * sorting test
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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
 * g++ -std=c++20 -I../ sort.cpp -o sort
 */

#include "libs/maths.h"
#include "libs/algos.h"

#include <iostream>


template<class t_vec>
void print(const t_vec& vec)
{
	for(std::size_t i=0; i<vec.size(); ++i)
		std::cout << vec[i] << " ";
	std::cout << std::endl;
}


int main()
{
	using t_vec = std::vector<double>;

	t_vec vec = tl2::rand<t_vec>(128);
	print(vec);

	std::vector<std::size_t> perm = tl2::get_perm(vec);
	vec = tl2::reorder(vec, perm);
	print(vec);

	return 0;
}
