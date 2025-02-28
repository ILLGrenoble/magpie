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

#ifndef __TLIBS2_LOADINSTR_BASE_H__
#define __TLIBS2_LOADINSTR_BASE_H__

#include <unordered_map>
#include <vector>
#include <array>
#include <memory>
#include <regex>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../file.h"
#include "../phys.h"
#include "../str.h"


namespace tl2 {


// interface for instrument-specific data files
template<class _t_real = double>
class FileInstrBase
{
	public:
		using t_real = _t_real;

		using t_mapParams = std::unordered_map<std::string, std::string>;
		using t_vecColNames = std::vector<std::string>;
		using t_vecVals = std::vector<t_real>;
		using t_vecDat = std::vector<t_vecVals>;

	protected:
		void RenameDuplicateCols();

		std::array<t_real, 5> GetScanHKLKiKf(const char* pcH, const char* pcK,
			const char* pcL, const char* pcE, std::size_t i) const;

	public:
		FileInstrBase() = default;
		virtual ~FileInstrBase() = default;

		virtual bool Load(const char* pcFile) = 0;

		virtual std::array<t_real, 3> GetSampleLattice() const = 0;
		virtual std::array<t_real, 3> GetSampleAngles() const = 0;
		virtual std::array<t_real, 2> GetMonoAnaD() const = 0;

		virtual std::array<bool, 3> GetScatterSenses() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane0() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane1() const = 0;

		virtual t_real GetKFix() const = 0;
		virtual bool IsKiFixed() const = 0;

		virtual bool HasCol(const std::string& strName) const;
		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const = 0;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) = 0;

		virtual std::array<t_real, 4> GetPosHKLE() const = 0;	// zero pos.

		virtual std::size_t GetScanCount() const = 0;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const = 0;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat);
		virtual void SmoothData(const std::string& strCol, t_real dEps, bool bIterate=true);

		virtual std::string GetTitle() const = 0;
		virtual std::string GetUser() const = 0;
		virtual std::string GetLocalContact() const = 0;
		virtual std::string GetScanNumber() const = 0;
		virtual std::string GetSampleName() const = 0;
		virtual std::string GetSpacegroup() const = 0;
		virtual std::string GetTimestamp() const = 0;

		virtual const t_vecDat& GetData() const = 0;
		virtual t_vecDat& GetData() = 0;
		virtual const t_vecColNames& GetColNames() const = 0;
		virtual const t_mapParams& GetAllParams() const = 0;

		virtual std::vector<std::string> GetScannedVars() const = 0;
		virtual std::string GetCountVar() const = 0;
		virtual std::string GetMonVar() const = 0;
		virtual std::string GetCountErr() const;
		virtual std::string GetMonErr() const;

		virtual std::string GetScanCommand() const = 0;

		// polarisation stuff
		virtual void ParsePolData();
		virtual std::size_t NumPolChannels() const;
		virtual const std::vector<std::array<t_real, 6>>& GetPolStates() const;
		virtual void SetPolNames(const char* pVec1, const char* pVec2,
			const char* pCur1, const char* pCur2);

	public:
		virtual bool MatchColumn(const std::string& strRegex,
			std::string& strColName, bool bSortByCounts=false, bool bFilterEmpty=true) const;

		static std::shared_ptr<FileInstrBase<t_real>> LoadInstr(const char* pcFile);
};


// ----------------------------------------------------------------------------
// implementations


template<class t_real>
void FileInstrBase<t_real>::RenameDuplicateCols()
{
	using t_mapCols = std::unordered_map<std::string, std::size_t>;
	t_mapCols mapCols;

	t_vecColNames& vecCols = const_cast<t_vecColNames&>(this->GetColNames());
	for(std::string& strCol : vecCols)
	{
		t_mapCols::iterator iter = mapCols.find(strCol);
		if(iter == mapCols.end())
		{
			mapCols.insert(std::make_pair(strCol, 0));
		}
		else
		{
			std::cerr << "Data file loader: "
				<< "Column \"" << strCol << "\" is duplicate, renaming it."
				<< std::endl;

			++iter->second;
			strCol += "_" + var_to_str(iter->second);
		}
	}
}


template<class t_real>
std::array<t_real, 5> FileInstrBase<t_real>::GetScanHKLKiKf(const char* pcH, const char* pcK,
	const char* pcL, const char* pcE, std::size_t i) const
{
	// zero position to fallback if no position is given in scan rows
	const std::array<t_real, 4> arrZeroPos = GetPosHKLE();

	const t_vecVals& vecH = GetCol(pcH);
	const t_vecVals& vecK = GetCol(pcK);
	const t_vecVals& vecL = GetCol(pcL);
	const t_vecVals& vecE = GetCol(pcE);

	t_real h = i < vecH.size() ? vecH[i] : std::get<0>(arrZeroPos);
	t_real k = i < vecK.size() ? vecK[i] : std::get<1>(arrZeroPos);
	t_real l = i < vecL.size() ? vecL[i] : std::get<2>(arrZeroPos);
	t_real E = i < vecE.size() ? vecE[i] : std::get<3>(arrZeroPos);

	bool bKiFix = IsKiFixed();
	t_real kfix = GetKFix();
	t_real kother = get_other_k<units::si::system, t_real>
		(E*meV<t_real>, kfix/angstrom<t_real>, bKiFix) * angstrom<t_real>;

	return std::array<t_real,5>{{h,k,l, bKiFix?kfix:kother, bKiFix?kother:kfix}};
}


template<class t_real>
bool FileInstrBase<t_real>::MatchColumn(const std::string& strRegex,
	std::string& strColName, bool bSortByCounts, bool bFilterEmpty) const
{
	const FileInstrBase<t_real>::t_vecColNames& vecColNames = GetColNames();
	std::regex rx(strRegex, std::regex::ECMAScript | std::regex_constants::icase);

	using t_pairCol = std::pair<std::string, t_real>;
	std::vector<t_pairCol> vecMatchedCols;

	for(const std::string& strCurColName : vecColNames)
	{
		std::smatch m;
		if(std::regex_match(strCurColName, m, rx))
		{
			const typename FileInstrBase<t_real>::t_vecVals& vecVals = GetCol(strCurColName);

			t_real dSum = std::accumulate(vecVals.begin(), vecVals.end(), 0.,
				[](t_real t1, t_real t2) -> t_real
				{
					return t1+t2;
				});

			if(!bFilterEmpty || !equals<t_real>(dSum, 0.))
				vecMatchedCols.push_back(t_pairCol{strCurColName, dSum});
		}
	}

	if(bSortByCounts)
	{
		std::sort(vecMatchedCols.begin(), vecMatchedCols.end(),
		[](const t_pairCol& pair1, const t_pairCol& pair2) -> bool
			{
				return pair1.second > pair2.second;
			});
	}

	if(vecMatchedCols.size())
	{
		strColName = vecMatchedCols[0].first;
		return true;
	}
	return false;
}


template<class t_real>
bool FileInstrBase<t_real>::MergeWith(const FileInstrBase<t_real>* pDat)
{
	if(this->GetColNames().size() != pDat->GetColNames().size())
	{
		std::cerr << "Data file loader: Cannot merge: Mismatching number of columns." << std::endl;
		return false;
	}

	for(const std::string& strCol : GetColNames())
	{
		t_vecVals& col1 = this->GetCol(strCol);
		const t_vecVals& col2 = pDat->GetCol(strCol);

		if(col1.size() == 0 || col2.size() == 0)
		{
			std::cerr << "Data file loader: "
				<< "Cannot merge: Column \"" << strCol << "\" is empty."
				<< std::endl;
			return false;
		}

		col1.insert(col1.end(), col2.begin(), col2.end());
	}

	return true;
}


template<class t_real>
void FileInstrBase<t_real>::SmoothData(const std::string& strCol, t_real dEps, bool bIterate)
{
	std::size_t iIdxCol;
	this->GetCol(strCol, &iIdxCol);		// get column index
	if(iIdxCol == GetColNames().size())	// no such column?
	{
		std::cerr << "Data file loader: "
			<< "No such data column: \"" << strCol << "\"." << std::endl;
		return;
	}

	while(true)
	{
		t_vecDat& vecDatOld = this->GetData();
		const std::size_t iNumCols = vecDatOld.size();
		const std::size_t iNumRows = vecDatOld[0].size();
		t_vecDat vecDatNew(iNumCols);
		std::vector<bool> vecValidRows(iNumRows, 1);

		for(std::size_t iPt1=0; iPt1<iNumRows; ++iPt1)
		{
			if(!vecValidRows[iPt1]) continue;

			t_vecVals vecVals(iNumCols, 0);
			std::size_t iNumUnited = 0;
			for(std::size_t iPt2=iPt1; iPt2<iNumRows; ++iPt2)
			{
				if(!vecValidRows[iPt2]) continue;
				if(std::abs(vecDatOld[iIdxCol][iPt1]-vecDatOld[iIdxCol][iPt2]) <= dEps)
				{
					for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
						vecVals[iCol] += vecDatOld[iCol][iPt2];
					++iNumUnited;
					vecValidRows[iPt2] = 0;
				}
			}

			for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
			{
				vecVals[iCol] /= t_real(iNumUnited);
				vecDatNew[iCol].push_back(vecVals[iCol]);
			}
		}

		vecDatOld = vecDatNew;
		// if no iteration requested or no more change -> break
		if(!bIterate || iNumRows==vecDatNew[0].size())
			break;
	}
}


template<class t_real>
bool FileInstrBase<t_real>::HasCol(const std::string& strName) const
{
	const t_vecColNames& cols = GetColNames();

	for(std::size_t i=0; i<cols.size(); ++i)
	{
		if(str_to_lower(cols[i]) == str_to_lower(strName))
			return true;
	}

	return false;
}


template<class t_real>
void FileInstrBase<t_real>::ParsePolData()
{}


template<class t_real>
void FileInstrBase<t_real>::SetPolNames(const char* /*pVec1*/, const char* /*pVec2*/,
	const char* /*pCur1*/, const char* /*pCur2*/)
{}


template<class t_real>
std::size_t FileInstrBase<t_real>::NumPolChannels() const
{
	return 0;
}


template<class t_real>
const std::vector<std::array<t_real, 6>>& FileInstrBase<t_real>::GetPolStates() const
{
	static const std::vector<std::array<t_real, 6>> vecNull;
	return vecNull;
}


template<class t_real>
std::string FileInstrBase<t_real>::GetCountErr() const
{
	return "";
}


template<class t_real>
std::string FileInstrBase<t_real>::GetMonErr() const
{
	return "";
}


// -----------------------------------------------------------------------------


}

#endif
