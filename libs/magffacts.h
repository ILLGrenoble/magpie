/**
 * loads tabulated magnetic form factors
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

#ifndef __MAGCORE_MAGFFTABLE__
#define __MAGCORE_MAGFFTABLE__


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <string>
#include <vector>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>

#include "tlibs2/libs/file.h"
#include "tlibs2/libs/phys.h"



/**
 * magnetic form factor for a given ion
 */
template<class t_real = double>
struct MagFormfactor
{
	std::string name{};   // ion name
	std::string terms{};  // term symbol
	std::tuple<t_real, t_real, t_real> SLJ;

	// j_0, j_2, j_4, ... coefficients
	std::vector<std::vector<t_real>> coefficients{};


	/**
	 * evaluate the form factor at the given Q
	 * @see https://mcphase.github.io/webpage/manual/node164.html
	 * @see https://github.com/SunnySuite/Sunny.jl/blob/main/src/FormFactor.jl
	 */
	t_real eval(t_real Q, t_real g = -1., bool is_Q = true) const
	{
		if(is_Q)
			Q /= 4. * tl2::pi<t_real>;

		if(g < 0.)  // calculate g
		{
			t_real S = std::get<0>(SLJ);
			t_real L = std::get<1>(SLJ);
			t_real J = std::get<2>(SLJ);

			g = tl2::eff_gJ<t_real>(S, L, J);
		}

		t_real ffact = 0.;

		for(std::size_t coeff_idx = 0; coeff_idx < coefficients.size(); ++coeff_idx)
		{
			const auto& coeffs = coefficients[coeff_idx];
			t_real prefactor = (coeff_idx == 0 ? 1. : Q*Q * (2./g - 1.));

			std::size_t i = 0;
			for(i = 0; i + 1 < coeffs.size(); i += 2)
				ffact += prefactor * coeffs[i] * std::exp(-Q*Q * coeffs[i + 1]);
			if(i < coeffs.size())
				ffact += prefactor * coeffs[i];
		}

		return ffact;
	}


	/**
	 * evaluate the form factor at the given Q
	 * @see https://mcphase.github.io/webpage/manual/node164.html
	 * @see https://github.com/SunnySuite/Sunny.jl/blob/main/src/FormFactor.jl
	 */
	std::string to_string(t_real g = -1.) const
	{
		// add brackets around negative values
		auto add_brackets = [](t_real val) -> std::string
		{
			std::ostringstream ostr;
			if(val >= 0.)
				ostr << val;
			else
				ostr << "(" << val << ")";
			return ostr.str();
		};

		if(g < 0.)  // calculate g
		{
			t_real S = std::get<0>(SLJ);
			t_real L = std::get<1>(SLJ);
			t_real J = std::get<2>(SLJ);

			g = tl2::eff_gJ<t_real>(S, L, J);
		}

		std::ostringstream expr;
		expr << "g = " << g << ";\n\n";

		for(std::size_t coeff_idx = 0; coeff_idx < coefficients.size(); ++coeff_idx)
		{
			const auto& coeffs = coefficients[coeff_idx];
			std::string prefactor_start = (coeff_idx == 0 ? "" : "s*s * (");
			std::string prefactor_end = (coeff_idx == 0 ? "" : ")");

			std::size_t i = 0;
			if(coeff_idx != 0)
				expr << "(2/g - 1)*";
			expr << prefactor_start;
			for(i = 0; i + 1 < coeffs.size(); i += 2)
			{
				expr << add_brackets(coeffs[i])
					<< " * exp(-s*s * " << add_brackets(coeffs[i + 1]) << ")";
				if(i < coeffs.size() - 3)
					expr << " +\n";
			}
			if(i < coeffs.size())
				expr << "\n+ " << add_brackets(coeffs[i]);
			expr << prefactor_end;

			if(coeff_idx < coefficients.size() - 1)
				expr << "\n+\n";
		}

		return expr.str();
	}
};



/**
 * table of magnetic form factors
 */
template<class t_real = double>
class MagFormfactorTable
{
public:
	using t_magffact = MagFormfactor<t_real>;


public:
	MagFormfactorTable() = default;
	~MagFormfactorTable() = default;


	/**
	 * loads form factors from an xml table
	 */
	bool LoadTable(const std::string& file)
	{
		try
		{
			if(!tl2::file_exists(file))
				return false;

			// read xml file
			std::ifstream ifstr{file};
			boost::property_tree::ptree node;
			boost::property_tree::read_xml(ifstr, node);

			// iterate ion list
			const auto& ffacts = node.get_child("magnetic_form_factors");
			for(const auto &ion : ffacts)
			{
				if(ion.first != "ion")
					continue;

				t_magffact ffact;
				ffact.name = ion.second.get<std::string>("<xmlattr>.name", "");
				if(ffact.name == "")
					continue;
				ffact.terms = ion.second.get<std::string>("<xmlattr>.terms", "");
				ffact.SLJ = tl2::from_termsymbol<t_real, std::string>(ffact.terms);

				for(std::size_t coeff_idx = 0; coeff_idx < 10; coeff_idx += 2)
				{
					std::vector<t_real> coeffs;

					tl2::get_tokens<t_real, std::string>(
						ion.second.get<std::string>(
						  "coefficients_" + tl2::var_to_str(coeff_idx), ""), " ", coeffs);

					// remove zero terms
					for(std::size_t i = 2; i + 1 < coeffs.size();)
					{
						if(tl2::equals_0(coeffs[i]))
							coeffs.erase(coeffs.begin() + i, coeffs.begin() + i + 2);
						else
							i += 2;
					}

					if(!coeffs.size())
						break;

					ffact.coefficients.emplace_back(std::move(coeffs));
				}

				m_formfactors.emplace_back(std::move(ffact));
			}
		}
		catch(const std::exception& ex)
		{
			std::cerr << "Form factor table: " << ex.what() << std::endl;
			return false;
		}

		return true;
	}


	const std::vector<t_magffact>& GetFormfactors() const
	{
		return m_formfactors;
	}


	const t_magffact* GetFormfactor(const std::string& name) const
	{
		for(const t_magffact& ffact : m_formfactors)
		{
			if(ffact.name == name)
				return &ffact;
		}

		return nullptr;
	}


private:
	std::vector<t_magffact> m_formfactors{};
};



#endif
