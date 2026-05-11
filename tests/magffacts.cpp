/**
 * tests loading tabulated magnetic form factors
 * @author Tobias Weber <tweber@ill.fr>
 * @date 11-may-2026
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

// clang++ -std=c++20 -I/opt/homebrew/include -I.. -o magffacts magffacts.cpp


#include "../libs/magffacts.h"

using t_real = double;


int main()
{
	const char* table = "../res/magffacts.xml";

	MagFormfactorTable<t_real> tab;
	if(!tab.LoadTable(table))
	{
		std::cerr << "Loading table \"" << table << "\" failed." << std::endl;
		return -1;
	}

	return 0;
}
