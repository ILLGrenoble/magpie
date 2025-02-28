/**
 * symmetry operator helpers
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
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

#ifndef __MAGCORE_SYMOPS_H__
#define __MAGCORE_SYMOPS_H__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/str.h"



/**
 * converts a symmetry operation matrix to a string
 */
template<class t_mat, class t_real = typename t_mat::value_type>
std::string op_to_str(const t_mat& op, int prec = -1, t_real eps = 1e-6)
{
	using namespace tl2_ops;

	std::ostringstream ostr;
	if(prec >= 0)
		ostr.precision(prec);

	for(std::size_t row = 0; row < op.size1(); ++row)
	{
		for(std::size_t col = 0; col < op.size2(); ++col)
		{
			t_real elem = op(row, col);
			tl2::set_eps_0(elem);

			// write periodic values as fractions to avoid
			// calculation errors involving epsilon
			if(tl2::equals<t_real>(elem, 1./3., eps))
				ostr << "1/3";
			else if(tl2::equals<t_real>(elem, 2./3., eps))
				ostr << "2/3";
			else if(tl2::equals<t_real>(elem, 1./6., eps))
				ostr << "1/6";
			else if(tl2::equals<t_real>(elem, 5./6., eps))
				ostr << "5/6";
			else if(tl2::equals<t_real>(elem, -1./3., eps))
				ostr << "-1/3";
			else if(tl2::equals<t_real>(elem, -2./3., eps))
				ostr << "-2/3";
			else if(tl2::equals<t_real>(elem, -1./6., eps))
				ostr << "-1/6";
			else if(tl2::equals<t_real>(elem, -5./6., eps))
				ostr << "-5/6";
			else
				ostr << elem;

			if(col != op.size2() - 1)
				ostr << " ";
		}

		if(row != op.size1()-1)
			ostr << " \n";
	}

	return ostr.str();
}



/**
 * convert a string to a symmetry operation matrix
 */
template<class t_mat>
t_mat str_to_op(const std::string& str)
{
	using namespace tl2_ops;

	using t_real = typename t_mat::value_type;
	t_mat op = tl2::unit<t_mat>(4);

	std::istringstream istr(str);
	for(std::size_t row = 0; row < op.size1(); ++row)
	{
		for(std::size_t col = 0; col < op.size2(); ++col)
		{
			std::string op_str;
			istr >> op_str;

			auto [ok, val] = tl2::eval_expr<std::string, t_real>(op_str);
			op(row, col) = val;

			if(!ok)
			{
				std::cerr << "Error evaluating symop component \""
					<< op_str << "\"" << std::endl;
			}
		}
	}

	return op;
}



/**
 * get the properties of a symmetry operation
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::string get_op_properties(const t_mat& op, t_real eps = 1e-6)
{
	using namespace tl2_ops;

	const t_mat rot = tl2::submat<t_mat>(op, 3, 3);
	const t_vec trans = tl2::subvec<t_vec>(tl2::col<t_mat, t_vec>(op, 3), 3);
	const t_mat sym = tl2::unit<t_mat>(3) - rot;

	std::string prop;

	if(tl2::is_unit<t_mat>(op, eps))
	{
		if(prop.size())
			prop += ", ";
		prop += "identity";
	}

	if(tl2::hom_is_centring<t_mat>(op, eps))
	{
		if(prop.size())
			prop += ", ";
		prop += "centring";
	}

	if(tl2::det<t_mat>(rot) < 0.)
	{
		if(prop.size())
			prop += ", ";
		prop += "reflecting";
	}

	// check if the translation cannot be removed by a trafo
	// @see https://sunwoo-kim.github.io/en/projects/nonsymmorphic/, p. 6
	if(!tl2::equals_0<t_vec>(trans, eps) &&
		tl2::equals_0<t_real>(tl2::det<t_mat>(sym), eps))
	{
		if(prop.size())
			prop += ", ";
		prop += "non-symmorphic";
	}

	return prop;
}



/**
 * checks for allowed Bragg reflections
 *
 * algorithm based on Clipper's HKL_class
 *   constructor in clipper/core/coords.cpp by K. Cowtan, 2013
 * @see http://www.ysbl.york.ac.uk/~cowtan/clipper/
 *
 * symmetry operation S on position r: R*r + t
 * F = sum<S>( exp(2*pi*i * (R*r + t)*G) )
 *   = sum<S>( exp(2*pi*i * ((R*r)*G + t*G)) )
 *   = sum<S>( exp(2*pi*i * (r*(G*R) + t*G)) )
 *   = sum<S>( exp(2*pi*i * (r*(G*R)))  *  exp(2*pi*i * (G*t)) )
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type,
template<class...> class t_cont = std::vector>
std::pair<bool, std::size_t>
is_reflection_allowed(const t_vec& Q, const t_cont<t_mat>& symops, t_real eps)
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
{
	for(std::size_t opidx = 0; opidx < symops.size(); ++opidx)
	{
		const t_mat& mat = symops[opidx];
		t_mat rot = tl2::submat<t_mat>(mat, 0,0, 3,3); // rotation part of the symop
		rot = tl2::trans(rot);                         // recip -> transpose

		// does Q not transform into itself?
		if(!tl2::equals<t_vec>(Q, rot*Q, eps))
			continue;

		t_vec trans = tl2::create<t_vec>({ mat(0,3), mat(1,3), mat(2,3) });

		// does Q translate to multiples of the lattice vector?
		if(tl2::is_integer<t_real>(tl2::inner<t_vec>(trans, Q), eps))
			continue;

		return std::make_pair(false, opidx);
	}

	return std::make_pair(true, symops.size());
}



/**
 * get an "x, y, z" form description of a symop
 */
template<class t_mat, class t_real = typename t_mat::value_type>
std::string symop_to_xyz(const t_mat& symop, int prec = -1, t_real eps = 1e-6)
requires tl2::is_mat<t_mat>
{
	using t_size = decltype(symop.size1());

	static const std::string varnames[] = { "x", "y", "z" };
	std::string symop_xyz;

	for(t_size row = 0; row < symop.size1() - 1; ++row)
	{
		std::ostringstream ostr;
		if(prec >= 0)
			ostr.precision(prec);

		for(t_size col = 0; col < symop.size2() - 1; ++col)
		{
			// rotation part
			t_real rval = symop(row, col);
			tl2::set_eps_0(rval, eps);

			if(!tl2::equals_0(rval, eps))
			{
				if(rval < t_real(0))
					ostr << "-";
				else
					ostr << "+";

				if(!tl2::equals(std::abs(rval), t_real(1), eps))
					ostr << std::abs(rval);

				ostr << varnames[col];
			}
		}

		// translation part
		t_real tval = symop(row, symop.size2() - 1);
		tl2::set_eps_0(tval, eps);

		if(!tl2::equals_0(tval, eps))
		{
			if(tval < t_real(0))
				ostr << "-";
			else
				ostr << "+";

			tval = std::abs(tval);

			if(tl2::equals<t_real>(tval, 1./2., eps))
				ostr << "1/2";
			else if(tl2::equals<t_real>(tval, 1./3., eps))
				ostr << "1/3";
			else if(tl2::equals<t_real>(tval, 2./3., eps))
				ostr << "2/3";
			else if(tl2::equals<t_real>(tval, 1./6., eps))
				ostr << "1/6";
			else if(tl2::equals<t_real>(tval, 5./6., eps))
				ostr << "5/6";
			else
				ostr << tval;
		}

		if(ostr.str().length() == 0)
			ostr << "0";

		if(row < symop.size1() - 2)
			ostr << ", ";

		symop_xyz += ostr.str();
	}

	return symop_xyz;
}


#endif
