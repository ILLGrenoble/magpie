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

#ifndef __TLIBS2_LOADINSTR_MACS_H__
#define __TLIBS2_LOADINSTR_MACS_H__


#include "base.h"


namespace tl2 {


// macs files
template<class _t_real = double>
class FileMacs : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams{};
		t_vecColNames m_vecQuantities{};
		t_vecDat m_vecData{};

	public:
		FileMacs() = default;
		virtual ~FileMacs() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// ----------------------------------------------------------------------------
// implementation


template<class t_real>
void FileMacs<t_real>::ReadHeader(std::istream& istr)
{
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);

		trim(strLine);
		if(strLine.length()==0 || strLine[0]!='#')
			continue;

		strLine = strLine.substr(1);

		std::pair<std::string, std::string> pairLine =
			split_first<std::string>(strLine, " \t", 1);

		if(pairLine.first == "")
			continue;
		else if(pairLine.first == "Columns")
		{
			get_tokens<std::string, std::string>(pairLine.second, " \t", m_vecQuantities);
			FileInstrBase<t_real>::RenameDuplicateCols();
			continue;
		}
		else
		{
			typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
	}
}


template<class t_real>
void FileMacs<t_real>::ReadData(std::istream& istr)
{
	m_vecData.resize(m_vecQuantities.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length() == 0 || strLine[0] == '#')
			continue;

		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecQuantities.size())
		{
			std::cerr << "Data file loader: Line size mismatch." << std::endl;

			// add zeros
			while(m_vecQuantities.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}
}


template<class t_real>
bool FileMacs<t_real>::Load(const char* pcFile)
{
	for(int iStep : { 0, 1 })
	{
		std::ifstream ifstr(pcFile);
		if(!ifstr.is_open())
			return false;

		if(iStep == 0)
			ReadHeader(ifstr);
		else if(iStep == 1)
			ReadData(ifstr);
	}

	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileMacs<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileMacs*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileMacs<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i = 0; i < m_vecQuantities.size(); ++i)
	{
		if(m_vecQuantities[i] == strName)
		{
			if(pIdx)
				*pIdx = i;
			return m_vecData[i];
		}
	}

	if(pIdx)
		*pIdx = m_vecQuantities.size();
	return vecNull;
}


template<class t_real>
std::array<t_real, 3> FileMacs<t_real>::GetSampleLattice() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Lattice");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vecToks;
	get_tokens<t_real, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		std::cerr << "Data file loader: Invalid sample lattice array size." << std::endl;
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{ vecToks[0], vecToks[1], vecToks[2] }};
}


template<class t_real>
std::array<t_real, 3> FileMacs<t_real>::GetSampleAngles() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Lattice");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vecToks;
	get_tokens<t_real, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		std::cerr << "Data file loader: Invalid sample lattice array size." << std::endl;
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{ d2r(vecToks[3]), d2r(vecToks[4]), d2r(vecToks[5]) }};
}


template<class t_real>
std::array<t_real, 2> FileMacs<t_real>::GetMonoAnaD() const
{
	typename t_mapParams::const_iterator iterM = m_mapParams.find("MonoSpacing");
	typename t_mapParams::const_iterator iterA = m_mapParams.find("AnaSpacing");

	t_real m = (iterM!=m_mapParams.end() ? str_to_var<t_real>(iterM->second) : 3.355);
	t_real a = (iterA!=m_mapParams.end() ? str_to_var<t_real>(iterA->second) : 3.355);

	return std::array<t_real,2>{{ m, a }};
}


template<class t_real>
std::array<bool, 3> FileMacs<t_real>::GetScatterSenses() const
{
	return std::array<bool,3>{{ false, true, false }};
}


template<class t_real>
std::array<t_real, 3> FileMacs<t_real>::GetScatterPlane0() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Orient");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vecToks;
	get_tokens<t_real, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		std::cerr << "Data file loader: Invalid sample orientation array size." << std::endl;
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{ vecToks[0], vecToks[1], vecToks[2] }};
}


template<class t_real>
std::array<t_real, 3> FileMacs<t_real>::GetScatterPlane1() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Orient");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vecToks;
	get_tokens<t_real, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		std::cerr << "Data file loader: Invalid sample orientation array size." << std::endl;
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{ vecToks[3], vecToks[4], vecToks[5]} };
}


template<class t_real>
std::array<t_real, 4> FileMacs<t_real>::GetPosHKLE() const
{
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol("QX");
	const t_vecVals& vecK = GetCol("QY");
	const t_vecVals& vecL = GetCol("QZ");
	const t_vecVals& vecE = GetCol("E");

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileMacs<t_real>::GetKFix() const
{
	// 1) look in data columns
	const std::string strKey = (IsKiFixed() ? "Ei" : "Ef");
	const t_vecVals& vecVals = GetCol(strKey);
	if(vecVals.size() != 0)
	{
		bool bImag;
		t_real k = E2k<units::si::system, t_real>
			(vecVals[0] * meV<t_real>, bImag) * angstrom<t_real>;
		return k;
	}


	// 2) look in header
	typename t_mapParams::const_iterator iter = m_mapParams.find("FixedE");
	if(iter == m_mapParams.end())
	{
		std::cerr << "Data file loader: Cannot determine kfix." << std::endl;
		return 0.;
	}

	std::vector<std::string> vecToks;
	get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

	if(vecToks.size() < 2)
	{
		std::cerr << "Data file loader: Cannot determine kfix." << std::endl;
		return 0.;
	}

	t_real dEfix = str_to_var<t_real>(vecToks[1]);
	bool bImag;
	t_real k = E2k<units::si::system, t_real>(dEfix * meV<t_real>, bImag) * angstrom<t_real>;
	return k;
}


template<class t_real>
bool FileMacs<t_real>::IsKiFixed() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("FixedE");
	if(iter==m_mapParams.end())
		return false;	// assume ckf

	std::vector<std::string> vecToks;
	get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

	if(vecToks.size() == 0)
		return false;	// assume ckf

	std::string strFixedE = vecToks[0];
	trim(strFixedE);

	if(strFixedE == "Ef")
		return false;
	else if(strFixedE == "Ei")
		return true;

	return false;		// assume ckf
}


template<class t_real>
std::size_t FileMacs<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FileMacs<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("QX", "QY", "QZ", "E", i);
}


template<class t_real>
bool FileMacs<t_real>::MergeWith(const FileInstrBase<t_real>* pDat)
{
	if(!FileInstrBase<t_real>::MergeWith(pDat))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		typename t_mapParams::iterator iter = m_mapParams.find("Filename");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


template<class t_real>
std::string FileMacs<t_real>::GetTitle() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("ExptID");
	if(iter != m_mapParams.end())
		strTitle = iter->second;

	iter = m_mapParams.find("ExptName");
	if(iter != m_mapParams.end() && iter->second != "")
		strTitle += " - " + iter->second;

	return strTitle;
}


template<class t_real>
std::string FileMacs<t_real>::GetUser() const
{
	std::string str;
	typename t_mapParams::const_iterator iter = m_mapParams.find("User");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}


template<class t_real>
std::string FileMacs<t_real>::GetLocalContact() const
{
	// TODO
	return "";
}


template<class t_real>
std::string FileMacs<t_real>::GetScanNumber() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Filename");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileMacs<t_real>::GetSampleName() const
{
	return "";
}


template<class t_real>
std::string FileMacs<t_real>::GetSpacegroup() const
{
	return "";
}


template<class t_real>
std::vector<std::string> FileMacs<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecScan;

	typename t_mapParams::const_iterator iter = m_mapParams.find("Scan");
	if(iter != m_mapParams.end())
	{
		std::vector<std::string> vecToks;
		get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

		if(vecToks.size() >= 2)
			vecScan.push_back(vecToks[1]);
	}

	if(!vecScan.size())
	{
		std::cerr << "Data file loader: Could not determine scan variable." << std::endl;
		if(m_vecQuantities.size() >= 1)
		{
			std::cerr << "Data file loader: "
				<< "Using first column: \"" << m_vecQuantities[0]
				<< "\"." << std::endl;
			vecScan.push_back(m_vecQuantities[0]);
		}
	}

	return vecScan;
}


template<class t_real>
std::string FileMacs<t_real>::GetCountVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn("spec[a-z0-9]*", strRet))
		return strRet;
	return "";
}


template<class t_real>
std::string FileMacs<t_real>::GetMonVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn("mon[a-z0-9]*", strRet))
		return strRet;
	return "";
}


template<class t_real>
std::string FileMacs<t_real>::GetScanCommand() const
{
	// TODO
	return "";
}


template<class t_real>
std::string FileMacs<t_real>::GetTimestamp() const
{
	std::string str;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Date");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}
// -----------------------------------------------------------------------------

}

#endif
