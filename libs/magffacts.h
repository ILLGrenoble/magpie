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
#include <fstream>
#include <iostream>

#include "tlibs2/libs/file.h"



/**
 * magnetic form factor for a given ion
 */
template<class t_real = double>
struct MagFormfactor
{
	std::string name{};   // ion name
	std::string terms{};  // term symbol

	std::vector<t_real> coefficients_0{};
	std::vector<t_real> coefficients_2{};
	std::vector<t_real> coefficients_4{};
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

				tl2::get_tokens<t_real, std::string>(
					ion.second.get<std::string>("coefficients_0", ""), " ", ffact.coefficients_0);
				tl2::get_tokens<t_real, std::string>(
					ion.second.get<std::string>("coefficients_2", ""), " ", ffact.coefficients_2);
				tl2::get_tokens<t_real, std::string>(
					ion.second.get<std::string>("coefficients_4", ""), " ", ffact.coefficients_4);

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


private:
	std::vector<t_magffact> m_formfactors{};
};



#endif
