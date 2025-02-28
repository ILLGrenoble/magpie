/**
 * magnetic space group library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 7-apr-18
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * magtools
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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

#include "libs/magsg.h"
#include "libs/maths.h"
using namespace m_ops;


int main()
{
	using t_real = double;
	using t_vec = std::vector<t_real>;
	using t_mat = tl2::mat<t_real, std::vector>;

	Spacegroups<t_mat, t_vec> sgs;
	sgs.Load("magsg.xml");

	return 0;
}
