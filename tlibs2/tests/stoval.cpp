/**
 * test locale
 * @author Tobias Weber <tweber@ill.fr>
 * @date 24-jul-2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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
 * g++ -std=c++20 -I.. -o stoval stoval.cpp
 */

#include <iostream>
#include <sstream>
#include <locale>

#include "libs/str.h"


int main()
{
	std::ios_base::sync_with_stdio(false);

	::setlocale(LC_ALL, "C");
	const std::locale& C_loc = std::locale::classic();
	std::locale::global(C_loc);

	std::string flt{"12.34"};
	std::cout << tl2::stoval<double>(flt) << std::endl;

	return 0;
}
