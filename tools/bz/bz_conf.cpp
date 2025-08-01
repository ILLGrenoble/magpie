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
#include "tlibs2/libs/maths.h"

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

	// alternate plane definition with two in-plane vectors
	boost::optional<t_real> cut_x2, cut_y2, cut_z2;
	cut_x2 = node.get_optional<t_real>("bz.cut.x2");
	cut_y2 = node.get_optional<t_real>("bz.cut.y2");
	cut_z2 = node.get_optional<t_real>("bz.cut.z2");
	if(cut_x2 && cut_y2 && cut_z2 &&
		cfg.cut_x && cfg.cut_y && cfg.cut_z &&
		cfg.xtal_a && cfg.xtal_b && cfg.xtal_c &&
		cfg.xtal_alpha && cfg.xtal_beta && cfg.xtal_gamma)
	{
		t_vec_bz vec1_rlu = tl2::create<t_vec_bz>({
			*cfg.cut_x, *cfg.cut_y, *cfg.cut_z });
		t_vec_bz vec2_rlu = tl2::create<t_vec_bz>({
			*cut_x2, *cut_y2, *cut_z2 });

		t_mat_bz crystB = tl2::B_matrix<t_mat_bz>(
			*cfg.xtal_a, *cfg.xtal_b, *cfg.xtal_c,
			tl2::d2r<t_real>(*cfg.xtal_alpha),
			tl2::d2r<t_real>(*cfg.xtal_beta),
			tl2::d2r<t_real>(*cfg.xtal_gamma));

		// get plane normal
		if(auto [crystBinv, inv_ok] = tl2::inv(crystB); inv_ok)
		{
			t_vec_bz vec1_invA = crystB * vec1_rlu;
			t_vec_bz vec2_invA = crystB * vec2_rlu;

			t_vec_bz norm_invA = tl2::cross<t_vec_bz>(vec1_invA, vec2_invA);

			t_vec_bz norm_rlu = crystBinv * norm_invA;
			norm_rlu /= tl2::norm<t_vec_bz>(norm_rlu);

			cfg.cut_nx = norm_rlu[0];
			cfg.cut_ny = norm_rlu[1];
			cfg.cut_nz = norm_rlu[2];
		}
	}


	// symops
	if(auto symops = node.get_child_optional("bz.symops"); symops)
	{
		for(const auto &symop : *symops)
		{
			std::string op = symop.second.get<std::string>(
				"", "1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1");

			cfg.symops.emplace_back(str_to_op<t_mat_bz>(op));
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
