/**
 * brillouin zone tool, configuration file
 * @author Tobias Weber <tweber@ill.fr>
 * @date Maz-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include "bz_conf.h"
#include "libs/symops.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;

using namespace tl2_ops;


/**
 * loads a configuration xml file
 */
BZConfig load_bz_config(const std::string& filename, bool use_stdin)
{
	std::ifstream ifstr;
	std::istream *istr = use_stdin ? &std::cin : &ifstr;

	if(!use_stdin)
	{
		ifstr.open(filename);
		if(!ifstr)
			throw std::runtime_error("Cannot open file \"" + filename + "\".");
	}

	pt::ptree node;
	pt::read_xml(*istr, node);

	// check signature
	if(auto opt = node.get_optional<std::string>("bz.meta.info");
		!opt || *opt!=std::string{"bz_tool"})
	{
		throw std::runtime_error("Unrecognised file format.");
	}


	// load configuration settings
	BZConfig cfg;
	cfg.xtal_a = node.get_optional<t_real>("bz.xtal.a");
	cfg.xtal_b = node.get_optional<t_real>("bz.xtal.b");
	cfg.xtal_c = node.get_optional<t_real>("bz.xtal.c");
	cfg.xtal_alpha = node.get_optional<t_real>("bz.xtal.alpha");
	cfg.xtal_beta = node.get_optional<t_real>("bz.xtal.beta");
	cfg.xtal_gamma = node.get_optional<t_real>("bz.xtal.gamma");
	cfg.order = node.get_optional<int>("bz.order");
	cfg.cut_order = node.get_optional<int>("bz.cut.order");
	cfg.cut_x = node.get_optional<t_real>("bz.cut.x");
	cfg.cut_y = node.get_optional<t_real>("bz.cut.y");
	cfg.cut_z = node.get_optional<t_real>("bz.cut.z");
	cfg.cut_nx = node.get_optional<t_real>("bz.cut.nx");
	cfg.cut_ny = node.get_optional<t_real>("bz.cut.ny");
	cfg.cut_nz = node.get_optional<t_real>("bz.cut.nz");
	cfg.cut_d = node.get_optional<t_real>("bz.cut.d");
	cfg.sg_idx = node.get_optional<int>("bz.sg_idx");


	// symops
	if(auto symops = node.get_child_optional("bz.symops"); symops)
	{
		for(const auto &symop : *symops)
		{
			std::string op = symop.second.get<std::string>(
				"", "1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1");

			cfg.symops.emplace_back(str_to_op<t_mat>(op));
		}
	}

	// formulas
	if(auto formulas = node.get_child_optional("bz.formulas"); formulas)
	{
		for(const auto &formula : *formulas)
		{
			auto expr = formula.second.get_optional<std::string>("");
			if(expr && *expr != "")
				cfg.formulas.push_back(*expr);
		}
	}


	return cfg;
}
