/**
 * generates symmetry-equivalent vectors
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "libs/loadcif.h"
#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/str.h"
#include "libs/defs.h"

#include <boost/algorithm/string.hpp>
#include <iostream>


using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;

constexpr t_real g_eps = 1e-6;


int main(int argc, char** argv)
{
	using namespace tl2_ops;

	std::string sgname;
	std::cout << "Enter space group name: ";
	std::getline(std::cin, sgname);

	std::vector<t_mat> ops = get_sg_ops<t_mat, t_real>(sgname);
	if(ops.size() == 0)
	{
		std::cerr << "Invalid space group." << std::endl;
		return -1;
	}

	for(std::size_t op_idx = 0; op_idx < ops.size(); ++op_idx)
	{
		std::cout << "Symmetry transformation " << op_idx + 1;
		std::cout << " (det_rot: ";
		std::cout << tl2::det<t_mat>(tl2::submat<t_mat>(ops[op_idx], 3, 3));
		std::cout  << "):\n";
		tl2::niceprint(std::cout, ops[op_idx], g_eps);
		std::cout << "\n";
	}


	while(true)
	{
		std::string strpos;
		std::cout << "Enter vector ['p' for pseudo-vector, 'e' or ENTER to end]: ";

		std::getline(std::cin, strpos);
		boost::trim(strpos);
		if(strpos == "" || strpos == "e")
			break;

		// tokenise
		std::vector<std::string> vecstr;
		boost::split(vecstr, strpos, [](char c) -> bool
		{
			return c == ' ' || c == '\t' || c == ',' || c == ';';
		}, boost::token_compress_on);


		// convert to real vector
		t_vec vec;
		bool pseudovec = false;
		for(const std::string& str : vecstr)
		{
			vec.emplace_back(tl2::stoval<t_real>(str));
			if(str == "p")
				pseudovec = true;
		}

		// fill up possibly missing coordinates
		while(vec.size() < 3)
			vec.push_back(0);
		vec.resize(3);

		std::cout << "Input " << (pseudovec ? "pseudo-vector: " : "vector: ") << vec << "\n";

		std::vector<t_vec> vecs = tl2::apply_ops_hom<t_vec, t_mat, t_real>(
			vec, ops, g_eps, false /*keep in uc*/, true /*ignore occupied*/,
			false /*return homogeneous*/, pseudovec);

		for(std::size_t sym_idx = 0; sym_idx < vecs.size(); ++sym_idx)
		{
			std::cout << "Symmetry-equivalent vector " << sym_idx + 1 << ": ";
			std::cout << vecs[sym_idx] << "\n";
		}
	}

	return 0;
}
