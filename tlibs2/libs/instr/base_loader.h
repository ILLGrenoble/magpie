/**
 * tlibs2 -- instrument-file library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015-2021
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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
 */

#ifndef __TLIBS2_LOADINSTR_BASE_LOADER_H__
#define __TLIBS2_LOADINSTR_BASE_LOADER_H__


#include "base.h"
#include "frm.h"
#include "psi.h"
#include "tax.h"
#include "macs.h"
#include "trisp.h"
#include "raw.h"


namespace tl2 {


// automatically choose correct instrument
template<class t_real>
std::shared_ptr<FileInstrBase<t_real>> FileInstrBase<t_real>::LoadInstr(const char* pcFile)
{
	std::shared_ptr<FileInstrBase<t_real>> pDat;

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return nullptr;

	std::string strLine, strLine2, strLine3;
	std::getline(ifstr, strLine);
	std::getline(ifstr, strLine2);
	std::getline(ifstr, strLine3);

	trim(strLine);
	trim(strLine2);
	trim(strLine3);
	strLine = str_to_lower(strLine);
	strLine2 = str_to_lower(strLine2);
	strLine3 = str_to_lower(strLine3);

	if(strLine == "")
		return nullptr;

	const std::string strNicos("nicos data file");
	const std::string strMacs("ice");
	const std::string strPsi("tas data");
	const std::string strPsiOld("instr:");
	const std::string strTax("scan =");

	if(strLine.find(strNicos) != std::string::npos)
	{ // frm file
		//std::cout << pcFile << " is an frm file." << std::endl;
		pDat = std::make_shared<FileFrm<t_real>>();
	}
	else if(strLine.find('#') != std::string::npos &&
		strLine.find(strMacs) != std::string::npos &&
		strLine2.find('#') != std::string::npos)
	{ // macs file
		//std::cout << pcFile << " is a macs file." << std::endl;
		pDat = std::make_shared<FileMacs<t_real>>();
	}
	else if(strLine2.find("scan start") != std::string::npos)
	{ // trisp file
		//std::cout << pcFile << " is a trisp file." << std::endl;
		pDat = std::make_shared<FileTrisp<t_real>>();
	}
	else if(strLine.find('#') == std::string::npos &&
		strLine2.find('#') == std::string::npos &&
		(strLine3.find(strPsi) != std::string::npos ||
		strLine.find(strPsiOld) != std::string::npos))
	{ // psi or ill file
		//std::cout << pcFile << " is an ill or psi file." << std::endl;
		pDat = std::make_shared<FilePsi<t_real>>();
	}
	else if(strLine.find('#') != std::string::npos &&
		strLine.find(strTax) != std::string::npos &&
		strLine2.find('#') != std::string::npos &&
		strLine3.find('#') != std::string::npos)
	{ // tax file
		//std::cout << pcFile << " is a tax file." << std::endl;
		pDat = std::make_shared<FileTax<t_real>>();
	}
	else
	{ // raw file
		std::cerr << "Data file loader: "
			<< "\"" << pcFile << "\" is of unknown type, falling back to raw loader."
			<< std::endl;
		pDat = std::make_shared<FileRaw<t_real>>();
	}

	if(pDat && !pDat->Load(pcFile))
		pDat = nullptr;

	return pDat;
}

}

#endif
