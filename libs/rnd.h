/**
 * random values
 * @author Tobias Weber <tweber@ill.fr>
 * @date 14-mar-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#ifndef __MAGCORE_RND__
#define __MAGCORE_RND__


#include <sstream>
#include <string>
#include <random>


template<class t_str = std::string>
t_str get_random_colour()
{
	static std::mt19937 rndgen{tl2::epoch<unsigned int>()};

	std::ostringstream ostrcol;
	std::uniform_int_distribution<int> dist{0, 255};

	ostrcol << "#" << std::hex << std::setw(2) << std::setfill('0') << dist(rndgen)
		<< std::setw(2) << std::setfill('0') << dist(rndgen)
		<< std::setw(2) << std::setfill('0') << dist(rndgen);

	return ostrcol.str();
}


#endif
