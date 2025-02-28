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

#ifndef __BZCONF_H__
#define __BZCONF_H__


#include <vector>
#include <string>
#include <boost/optional.hpp>

#include "globals.h"


/**
 * bz configuration for file loading
 */
struct BZConfig
{
	boost::optional<t_real> xtal_a, xtal_b, xtal_c;
	boost::optional<t_real> xtal_alpha, xtal_beta, xtal_gamma;

	boost::optional<int> order, cut_order;

	boost::optional<t_real> cut_x, cut_y, cut_z;
	boost::optional<t_real> cut_nx, cut_ny, cut_nz;
	boost::optional<t_real> cut_d;

	boost::optional<int> sg_idx;

	std::vector<t_mat> symops;
	std::vector<std::string> formulas;
};


extern BZConfig load_bz_config(const std::string& filename, bool use_stdin);


#endif
