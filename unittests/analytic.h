/**
 * analytical results for comparions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-nov-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAGPIE_ANALYTICS_H__
#define __MAGPIE_ANALYTICS_H__

#include <utility>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/phys.h"



/**
 * 1-d ferromagnetic chain
 */
template<class t_real = double>
std::pair<t_real, t_real> ferro(t_real h, t_real J = -1., t_real S = 1.)
{
	t_real E = -2.*S*J*(1. - std::cos(h * 2. * tl2::pi<t_real>));
	t_real w = std::pow(tl2::g_e<t_real>, 2.) * S / 2.;

	return std::make_pair(E, w);
}



/**
 * 1-d antiferromagnetic chain
 */
template<class t_real = double>
std::pair<t_real, t_real> antiferro(t_real h, t_real J = 1., t_real S = 1.)
{
	t_real E = 2.*S*J*std::abs(sin(h * 2. * tl2::pi<t_real>));
	t_real w = std::pow(tl2::g_e<t_real>, 2.) * std::abs(std::tan(h * tl2::pi<t_real>)) * S / 2.;

	return std::make_pair(E, w);
}



#endif
