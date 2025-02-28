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

#ifndef __TLIBS2_LOADINSTR_PSI_H__
#define __TLIBS2_LOADINSTR_PSI_H__


#include "base.h"


namespace tl2{


// psi & ill files
template<class _t_real = double>
class FilePsi : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

		// internal parameters in m_mapParams
		using t_mapIParams = std::map<std::string, t_real>;

	protected:
		t_mapParams m_mapParams{};
		t_mapIParams m_mapParameters{};
		t_mapIParams m_mapZeros{};
		t_mapIParams m_mapVariables{};
		t_mapIParams m_mapPosHkl{};
		t_mapIParams m_mapScanSteps{};
		t_vecColNames m_vecColNames{};
		t_vecDat m_vecData{};

		// automatically look and parse for polarisation data
		bool m_bAutoParsePol = false;

		// incoming and outgoing polarisation states
		std::vector<std::array<t_real, 6>> m_vecPolStates{};

		// instrument-specific device names
		std::string m_strPolVec1 = "p1";
		std::string m_strPolVec2 = "p2";
		std::string m_strPolCur1 = "i1";
		std::string m_strPolCur2 = "i2";

	protected:
		void ReadData(std::istream& istr);
		void GetInternalParams(const std::string& strAll, t_mapIParams& mapPara);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		virtual bool Load(const char* pcFile) override;

		void PrintParams(std::ostream& ostr) const;
		const t_mapParams& GetParams() const { return m_mapParams; }

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }

		const t_vecVals& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

	public:
		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos
		std::array<t_real, 4> GetDeltaHKLE() const;	// scan steps

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecColNames; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;

		virtual std::size_t NumPolChannels() const override;

	public:
		void SetAutoParsePolData(bool b) { m_bAutoParsePol = b; }
		virtual void ParsePolData() override;

		virtual const std::vector<std::array<t_real, 6>>& GetPolStates() const override
		{
			return m_vecPolStates;
		}

		virtual void SetPolNames(const char* pVec1, const char* pVec2,
			const char* pCur1, const char* pCur2) override
		{
			m_strPolVec1 = pVec1; m_strPolVec2 = pVec2;
			m_strPolCur1 = pCur1; m_strPolCur2 = pCur2;
		}
};



// ----------------------------------------------------------------------------
// implementation

template<class t_real>
void FilePsi<t_real>::ReadData(std::istream& istr)
{
	std::size_t iLine = 0;

	// header
	std::string strHdr;
	std::getline(istr, strHdr);
	trim(strHdr);
	++iLine;

	get_tokens<std::string, std::string, t_vecColNames>(strHdr, " \t", m_vecColNames);
	for(std::string& _str : m_vecColNames)
	{
		boost::replace_all<std::string, std::string, std::string>(_str, "\n", "");
		boost::replace_all<std::string, std::string, std::string>(_str, "\r", "");
	}

	m_vecData.resize(m_vecColNames.size());
	FileInstrBase<t_real>::RenameDuplicateCols();


	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		++iLine;
		trim(strLine);

		if(strLine.length() == 0)
			continue;
		if(strLine[0] == '#')
			continue;

		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecColNames.size())
		{
			std::cerr << "Data file loader: "
				<< "Column size mismatch in data line " << iLine
				<< ": Expected " << m_vecColNames.size()
				<< ", got " << vecToks.size()
				<< "." << std::endl;

			// add zeros
			while(m_vecColNames.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}
}


template<class t_real>
void FilePsi<t_real>::GetInternalParams(const std::string& strAll, FilePsi<t_real>::t_mapIParams& mapPara)
{
	std::vector<std::string> vecToks;
	get_tokens<std::string, std::string>(strAll, ",\n", vecToks);

	for(const std::string& strTok : vecToks)
	{
		std::pair<std::string, std::string> pair =
				split_first<std::string>(strTok, "=", 1);

		if(pair.first == "")
			continue;

		t_real dVal = str_to_var<t_real>(pair.second);
		mapPara.insert(typename t_mapIParams::value_type(pair.first, dVal));
	}
}


template<class t_real>
void FilePsi<t_real>::ParsePolData()
{
	m_vecPolStates.clear();
	typename t_mapParams::const_iterator iter = m_mapParams.find("POLAN");
	if(iter == m_mapParams.end())
		return;

	std::vector<std::string> vecLines;
	get_tokens<std::string, std::string>(iter->second, ",", vecLines);

	// initial and final polarisation states
	t_real Pix = t_real(0), Piy = t_real(0), Piz = t_real(0);
	t_real Pfx = t_real(0), Pfy = t_real(0), Pfz = t_real(0);
	t_real Pi_sign = t_real(1);
	t_real Pf_sign = t_real(1);


	// check if a string describes a number
	auto is_num = [](const std::string& str) -> bool
	{
		for(typename std::string::value_type c : str)
		{
			if(c!='0' && c!='1' && c!='2' && c!='3' && c!='4'
				&& c!='5' && c!='6' && c!='7' && c!='8' && c!='9'
				&& c!='+'&& c!='-'&& c!='.')
				return false;
		}
		return true;
	};


	// iterate command lines
	for(std::string& strLine : vecLines)
	{
		trim(strLine);
		strLine = str_to_lower(strLine);

		std::vector<std::string> vecLine;
		get_tokens<std::string, std::string>(strLine, " \t", vecLine);

		if(vecLine.size() == 0)
			continue;

		if(vecLine[0] == "dr")	// polarisation vector or current driven
		{
			std::string strCurDev = "";

			std::size_t iCurComp = 0;
			for(std::size_t iDr=1; iDr<vecLine.size(); ++iDr)
			{
				const std::string& strWord = vecLine[iDr];

				if(is_num(strWord))	// value to drive to
				{
					t_real dNum = str_to_var<t_real>(strWord);

					if(strCurDev == m_strPolVec1)
					{	// incoming polarisation vector changed
						switch(iCurComp)
						{
							case 0: Pix = dNum; break;
							case 1: Piy = dNum; break;
							case 2: Piz = dNum; break;
						}
					}
					else if(strCurDev == m_strPolVec2)
					{	// outgoing polarisation vector changed
						switch(iCurComp)
						{
							case 0: Pfx = dNum; break;
							case 1: Pfy = dNum; break;
							case 2: Pfz = dNum; break;
						}
					}
					else if(strCurDev == m_strPolCur1)
					{	// sign of polarisation vector 1 changed
						if(iCurComp == 0)
						{
							if(dNum >= t_real(0))
								Pi_sign = t_real(1);
							else
								Pi_sign = t_real(-1);
						}
					}
					else if(strCurDev == m_strPolCur2)
					{	// sign of polarisation vector 2 changed
						if(iCurComp == 0)
						{
							if(dNum >= t_real(0))
								Pf_sign = t_real(1);
							else
								Pf_sign = t_real(-1);
						}
					}

					++iCurComp;
				}
				else	// (next) device to drive
				{
					strCurDev = strWord;
					iCurComp = 0;
				}
			}
		}
		else if(vecLine[0] == "co")	// count command issued -> save current spin states
		{
			m_vecPolStates.push_back(std::array<t_real,6>({{
				Pi_sign*Pix, Pi_sign*Piy, Pi_sign*Piz,
				Pf_sign*Pfx, Pf_sign*Pfy, Pf_sign*Pfz }}));

			//std::cout << Pi_sign*Pix << " " << Pi_sign*Piy << " " << Pi_sign*Piz << " -> "
			//	<< Pf_sign*Pfx << " " << Pf_sign*Pfy << " " << Pf_sign*Pfz << std::endl;
		}
	}


	// cleanup
	for(std::size_t iPol=0; iPol<m_vecPolStates.size(); ++iPol)
	{
		for(unsigned iComp=0; iComp<6; ++iComp)
			set_eps_0(m_vecPolStates[iPol][iComp]);
	}
}


template<class t_real>
bool FilePsi<t_real>::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);

		if(strLine.substr(0,4) == "RRRR")
			skip_after_line<char>(ifstr, "VVVV", true);

		std::pair<std::string, std::string> pairLine =
				split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "DATA_")
			ReadData(ifstr);
		else if(pairLine.first == "")
			continue;
		else
		{
			typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
	}

	typename t_mapParams::const_iterator iterParams = m_mapParams.find("PARAM"),
		iterZeros = m_mapParams.find("ZEROS"),
		iterVars = m_mapParams.find("VARIA"),
		iterPos = m_mapParams.find("POSQE"),
		iterSteps = m_mapParams.find("STEPS");

	if(iterParams!=m_mapParams.end()) GetInternalParams(iterParams->second, m_mapParameters);
	if(iterZeros!=m_mapParams.end()) GetInternalParams(iterZeros->second, m_mapZeros);
	if(iterVars!=m_mapParams.end()) GetInternalParams(iterVars->second, m_mapVariables);
	if(iterPos!=m_mapParams.end()) GetInternalParams(iterPos->second, m_mapPosHkl);
	if(iterSteps!=m_mapParams.end()) GetInternalParams(iterSteps->second, m_mapScanSteps);

	if(m_bAutoParsePol)
		ParsePolData();
	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FilePsi<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FilePsi*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FilePsi<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i = 0; i < m_vecColNames.size(); ++i)
	{
		if(str_to_lower(m_vecColNames[i]) == str_to_lower(strName))
		{
			if(pIdx)
				*pIdx = i;
			return m_vecData[i];
		}
	}

	if(pIdx)
		*pIdx = m_vecColNames.size();
	return vecNull;
}


template<class t_real>
void FilePsi<t_real>::PrintParams(std::ostream& ostr) const
{
	for(const typename t_mapParams::value_type& val : m_mapParams)
	{
		ostr << "Param: " << val.first
			<< ", Val: " << val.second << "\n";
	}
}


template<class t_real>
std::array<t_real,3> FilePsi<t_real>::GetSampleLattice() const
{
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("AS");
	typename t_mapIParams::const_iterator iterB = m_mapParameters.find("BS");
	typename t_mapIParams::const_iterator iterC = m_mapParameters.find("CS");

	t_real a = (iterA!=m_mapParameters.end() ? iterA->second : 0.);
	t_real b = (iterB!=m_mapParameters.end() ? iterB->second : 0.);
	t_real c = (iterC!=m_mapParameters.end() ? iterC->second : 0.);

	return std::array<t_real,3>{{a,b,c}};
}


template<class t_real>
std::array<t_real,3> FilePsi<t_real>::GetSampleAngles() const
{
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("AA");
	typename t_mapIParams::const_iterator iterB = m_mapParameters.find("BB");
	typename t_mapIParams::const_iterator iterC = m_mapParameters.find("CC");

	t_real alpha = (iterA!=m_mapParameters.end() ? d2r(iterA->second) : pi<t_real>/2.);
	t_real beta = (iterB!=m_mapParameters.end() ? d2r(iterB->second) : pi<t_real>/2.);
	t_real gamma = (iterC!=m_mapParameters.end() ? d2r(iterC->second) : pi<t_real>/2.);

	return std::array<t_real,3>{{alpha, beta, gamma}};
}

template<class t_real>
std::array<t_real,2> FilePsi<t_real>::GetMonoAnaD() const
{
	typename t_mapIParams::const_iterator iterM = m_mapParameters.find("DM");
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("DA");

	t_real m = (iterM!=m_mapParameters.end() ? iterM->second : 3.355);
	t_real a = (iterA!=m_mapParameters.end() ? iterA->second : 3.355);

	return std::array<t_real,2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FilePsi<t_real>::GetScatterSenses() const
{
	typename t_mapIParams::const_iterator iterM = m_mapParameters.find("SM");
	typename t_mapIParams::const_iterator iterS = m_mapParameters.find("SS");
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("SA");

	bool m = (iterM!=m_mapParameters.end() ? (iterM->second>0) : false);
	bool s = (iterS!=m_mapParameters.end() ? (iterS->second>0) : true);
	bool a = (iterA!=m_mapParameters.end() ? (iterA->second>0) : false);

	return std::array<bool,3>{{ m, s, a }};
}


template<class t_real>
std::array<t_real, 3> FilePsi<t_real>::GetScatterPlane0() const
{
	typename t_mapIParams::const_iterator iterX = m_mapParameters.find("AX");
	typename t_mapIParams::const_iterator iterY = m_mapParameters.find("AY");
	typename t_mapIParams::const_iterator iterZ = m_mapParameters.find("AZ");

	t_real x = (iterX!=m_mapParameters.end() ? iterX->second : 1.);
	t_real y = (iterY!=m_mapParameters.end() ? iterY->second : 0.);
	t_real z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
std::array<t_real, 3> FilePsi<t_real>::GetScatterPlane1() const
{
	typename t_mapIParams::const_iterator iterX = m_mapParameters.find("BX");
	typename t_mapIParams::const_iterator iterY = m_mapParameters.find("BY");
	typename t_mapIParams::const_iterator iterZ = m_mapParameters.find("BZ");

	t_real x = (iterX!=m_mapParameters.end() ? iterX->second : 0.);
	t_real y = (iterY!=m_mapParameters.end() ? iterY->second : 1.);
	t_real z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
t_real FilePsi<t_real>::GetKFix() const
{
	typename t_mapIParams::const_iterator iterK = m_mapParameters.find("KFIX");
	t_real k = (iterK!=m_mapParameters.end() ? iterK->second : 0.);

	return k;
}


template<class t_real>
std::array<t_real, 4> FilePsi<t_real>::GetPosHKLE() const
{
	typename t_mapIParams::const_iterator iterH = m_mapPosHkl.find("QH");
	typename t_mapIParams::const_iterator iterK = m_mapPosHkl.find("QK");
	typename t_mapIParams::const_iterator iterL = m_mapPosHkl.find("QL");
	typename t_mapIParams::const_iterator iterE = m_mapPosHkl.find("EN");

	t_real h = (iterH!=m_mapPosHkl.end() ? iterH->second : 0.);
	t_real k = (iterK!=m_mapPosHkl.end() ? iterK->second : 0.);
	t_real l = (iterL!=m_mapPosHkl.end() ? iterL->second : 0.);
	t_real E = (iterE!=m_mapPosHkl.end() ? iterE->second : 0.);

	return std::array<t_real,4>{{h,k,l,E}};
}


template<class t_real>
std::array<t_real, 4> FilePsi<t_real>::GetDeltaHKLE() const
{
	typename t_mapIParams::const_iterator iterH = m_mapScanSteps.find("DQH");
	if(iterH==m_mapScanSteps.end()) iterH = m_mapScanSteps.find("QH");

	typename t_mapIParams::const_iterator iterK = m_mapScanSteps.find("DQK");
	if(iterK==m_mapScanSteps.end()) iterK = m_mapScanSteps.find("QK");

	typename t_mapIParams::const_iterator iterL = m_mapScanSteps.find("DQL");
	if(iterL==m_mapScanSteps.end()) iterL = m_mapScanSteps.find("QL");

	typename t_mapIParams::const_iterator iterE = m_mapScanSteps.find("DEN");
	if(iterE==m_mapScanSteps.end()) iterE = m_mapScanSteps.find("EN");


        t_real h = (iterH!=m_mapScanSteps.end() ? iterH->second : 0.);
        t_real k = (iterK!=m_mapScanSteps.end() ? iterK->second : 0.);
        t_real l = (iterL!=m_mapScanSteps.end() ? iterL->second : 0.);
        t_real E = (iterE!=m_mapScanSteps.end() ? iterE->second : 0.);

        return std::array<t_real,4>{{h,k,l,E}};
}


template<class t_real>
bool FilePsi<t_real>::MergeWith(const FileInstrBase<t_real>* pDat)
{
	if(!FileInstrBase<t_real>::MergeWith(pDat))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		typename t_mapParams::iterator iter = m_mapParams.find("FILE_");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


template<class t_real>
bool FilePsi<t_real>::IsKiFixed() const
{
	typename t_mapIParams::const_iterator iter = m_mapParameters.find("FX");
	t_real val = (iter!=m_mapParameters.end() ? iter->second : 2.);

	return equals<t_real>(val, 1., 0.25);
}


template<class t_real>
std::size_t FilePsi<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FilePsi<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	// default column names
	const char *h = "QH";
	const char *k = "QK";
	const char *l = "QL";
	const char *E = "EN";

	// alternate column names
	if(!FileInstrBase<t_real>::HasCol("QH") && FileInstrBase<t_real>::HasCol("H"))
		h = "H";
	if(!FileInstrBase<t_real>::HasCol("QK") && FileInstrBase<t_real>::HasCol("K"))
		k = "K";
	if(!FileInstrBase<t_real>::HasCol("QL") && FileInstrBase<t_real>::HasCol("L"))
		l = "L";
	if(!FileInstrBase<t_real>::HasCol("EN") && FileInstrBase<t_real>::HasCol("E"))
		E = "E";

	return FileInstrBase<t_real>::GetScanHKLKiKf(h, k, l, E, i);
}


template<class t_real>
std::string FilePsi<t_real>::GetTitle() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("TITLE");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FilePsi<t_real>::GetUser() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("USER_");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FilePsi<t_real>::GetLocalContact() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("LOCAL");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FilePsi<t_real>::GetScanNumber() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("FILE_");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FilePsi<t_real>::GetSampleName() const
{
	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetSpacegroup() const
{
	return "";
}


template<class t_real>
std::vector<std::string> FilePsi<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// steps parameter
	for(const typename t_mapIParams::value_type& pair : m_mapScanSteps)
	{
		if(!equals<t_real>(pair.second, 0.) && pair.first.length())
		{
			if(std::tolower(pair.first[0]) == 'd')
				vecVars.push_back(pair.first.substr(1));
			else
				vecVars.push_back(pair.first);
		}
	}


	// nothing found yet -> try scan command instead
	if(!vecVars.size())
	{
		typename t_mapParams::const_iterator iter = m_mapParams.find("COMND");
		if(iter != m_mapParams.end())
		{
			std::vector<std::string> vecToks;
			get_tokens<std::string, std::string>(iter->second, " \t", vecToks);
			for(std::string& strTok : vecToks)
				trim(strTok);

			std::transform(vecToks.begin(), vecToks.end(), vecToks.begin(), str_to_lower<std::string>);
			typename std::vector<std::string>::iterator iterTok
				= std::find(vecToks.begin(), vecToks.end(), "dqh");

			if(iterTok != vecToks.end())
			{
				t_real dh = str_to_var<t_real>(*(++iterTok));
				t_real dk = str_to_var<t_real>(*(++iterTok));
				t_real dl = str_to_var<t_real>(*(++iterTok));
				t_real dE = str_to_var<t_real>(*(++iterTok));

				if(!equals<t_real>(dh, 0.)) vecVars.push_back("QH");
				if(!equals<t_real>(dk, 0.)) vecVars.push_back("QK");
				if(!equals<t_real>(dl, 0.)) vecVars.push_back("QL");
				if(!equals<t_real>(dE, 0.)) vecVars.push_back("EN");
			}


			// still nothing found, try regex
			if(!vecVars.size())
			{
				const std::string strRegex = "(sc|scan)[ \\t]+([a-z0-9]+)[ \\t]+[0-9\\.-]+[ \\t]+[d|D]([a-z0-9]+).*";
				std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);
				std::smatch m;
				if(std::regex_search(iter->second, m, rx) && m.size()>3)
				{
					const std::string& strSteps = m[3];
					vecVars.push_back(str_to_upper(strSteps));
				}
			}
		}
	}

	if(!vecVars.size())
	{
		std::cerr << "Data file loader: Could not determine scan variable." << std::endl;
		if(m_vecColNames.size() >= 1)
		{
			std::cerr << "Data file loader: "
				<< "Using first column: \"" << m_vecColNames[0] <<
				"\"." << std::endl;
			vecVars.push_back(m_vecColNames[0]);
		}
	}

	return vecVars;
}


template<class t_real>
std::string FilePsi<t_real>::GetCountVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn("cnts", strRet))
		return strRet;
	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetMonVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn("m[0-9]", strRet))
		return strRet;
	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetScanCommand() const
{
	std::string strCmd;
	typename t_mapParams::const_iterator iter = m_mapParams.find("COMND");
	if(iter != m_mapParams.end())
		strCmd = iter->second;
	return strCmd;
}


template<class t_real>
std::string FilePsi<t_real>::GetTimestamp() const
{
	std::string strDate;
	typename t_mapParams::const_iterator iter = m_mapParams.find("DATE_");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}


template<class t_real>
std::size_t FilePsi<t_real>::NumPolChannels() const
{
	return m_vecPolStates.size();
}
// ----------------------------------------------------------------------------

}

#endif
