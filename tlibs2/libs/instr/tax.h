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

#ifndef __TLIBS2_LOADINSTR_TAX_H__
#define __TLIBS2_LOADINSTR_TAX_H__


#include "base.h"
#include "raw.h"


namespace tl2{

/**
 * tax data files
 */
template<class _t_real = double>
class FileTax : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		DatFile<t_real, char> m_dat;
		t_vecColNames m_vecCols;

	public:
		FileTax() = default;
		virtual ~FileTax() = default;

	protected:
		std::array<t_real, 3> GetScatterPlaneVector(int i) const;

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

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

		virtual const t_vecDat& GetData() const override;
		virtual t_vecDat& GetData() override;
		virtual const t_vecColNames& GetColNames() const override;
		virtual const t_mapParams& GetAllParams() const override;

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};



// ----------------------------------------------------------------------------
// implementation

template<class t_real>
bool FileTax<t_real>::Load(const char* pcFile)
{
	// load data columns
	m_dat.SetCommentChar('#');
	m_dat.SetSeparatorChars("=");
	bool bOk = m_dat.Load(pcFile);

	// get column header names
	m_vecCols.clear();

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);
		trim<std::string>(strLine);
		if(strLine.length() == 0)
			continue;

		if(strLine == "# col_headers =")
		{
			// column headers in next line
			std::getline(ifstr, strLine);
			trim<std::string>(strLine);
			strLine = strLine.substr(1);

			get_tokens<std::string, std::string>(strLine, " \t", m_vecCols);
			break;
		}
	}

	if(m_vecCols.size() != m_dat.GetColumnCount())
	{
		std::cerr << "Data file loader: "
			<< "Mismatch between the number of data columns ("
			<< m_dat.GetColumnCount() << ") and column headers ("
			<< m_vecCols.size() << ")."
			<< std::endl;
	}

	// fill the rest with dummy column names
	for(std::size_t iCol=m_vecCols.size(); iCol<m_dat.GetColumnCount(); ++iCol)
		m_vecCols.emplace_back(var_to_str(iCol+1));

	return bOk;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileTax<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileTax*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileTax<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	auto iter = std::find(m_vecCols.begin(), m_vecCols.end(), strName);
	if(iter == m_vecCols.end())
		return vecNull;

	std::size_t iCol = iter - m_vecCols.begin();
	if(iCol < m_dat.GetColumnCount())
	{
		if(pIdx)
			*pIdx = iCol;
		return m_dat.GetColumn(iCol);
	}

	if(pIdx)
		*pIdx = m_dat.GetColumnCount();

	return vecNull;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecDat&
FileTax<t_real>::GetData() const
{
	return m_dat.GetData();
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecDat&
FileTax<t_real>::GetData()
{
	return m_dat.GetData();
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecColNames&
FileTax<t_real>::GetColNames() const
{
	return m_vecCols;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_mapParams&
FileTax<t_real>::GetAllParams() const
{
	return m_dat.GetHeader();
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetSampleLattice() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	typename t_map::const_iterator iter = params.find("latticeconstants");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);

		if(vecToks.size() > 0)
			a = vecToks[0];
		if(vecToks.size() > 1)
			b = vecToks[1];
		if(vecToks.size() > 2)
			c = vecToks[2];
	}

	return std::array<t_real, 3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetSampleAngles() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	typename t_map::const_iterator iter = params.find("latticeconstants");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);

		if(vecToks.size() > 3)
			a = d2r(vecToks[3]);
		if(vecToks.size() > 4)
			b = d2r(vecToks[4]);
		if(vecToks.size() > 5)
			c = d2r(vecToks[5]);
	}

	return std::array<t_real, 3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 2> FileTax<t_real>::GetMonoAnaD() const
{
	t_real m{0}, a{0};

	// TODO

	return std::array<t_real, 2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FileTax<t_real>::GetScatterSenses() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool m{0}, s{1}, a{0};

	typename t_map::const_iterator iter = params.find("sense");
	if(iter != params.end())
	{
		if(iter->second.length() > 0)
			m = (iter->second[0] == '+');
		if(iter->second.length() > 1)
			s = (iter->second[1] == '+');
		if(iter->second.length() > 2)
			a = (iter->second[2] == '+');
	}

	return std::array<bool, 3>{{m, s, a}};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlaneVector(int /*i*/) const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	//const t_vecVals& colGl = GetCol("sgl");
	//const t_vecVals& colGu = GetCol("sgu");
	//t_real gl = d2r(mean(colGl));
	//t_real gu = d2r(mean(colGu));

	t_real x{0}, y{0}, z{0};
	typename t_map::const_iterator iter = params.find("ubmatrix");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);

		//std::array<t_real, 3> lattice = GetSampleLattice();
		//std::array<t_real, 3> angles = GetSampleAngles();

		// TODO
		std::cerr << "Error: Scattering plane determination is not yet implemented for this file type." << std::endl;
	}

	return std::array<t_real, 3>{{ x, y, z }};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlane0() const
{
	return GetScatterPlaneVector(0);
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlane1() const
{
	return GetScatterPlaneVector(1);
}


template<class t_real>
std::array<t_real, 4> FileTax<t_real>::GetPosHKLE() const
{
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol("h");
	const t_vecVals& vecK = GetCol("k");
	const t_vecVals& vecL = GetCol("l");
	const t_vecVals& vecE = GetCol("e");

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileTax<t_real>::GetKFix() const
{
	bool ki_fixed = IsKiFixed();
	const t_vecVals& colEfix = GetCol(ki_fixed ? "ei" : "ef");

	t_real E = mean(colEfix);
	bool imag;
	t_real k = E2k<units::si::system, t_real>(E * meV<t_real>, imag) * angstrom<t_real>;
	return k;
}


template<class t_real>
bool FileTax<t_real>::IsKiFixed() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool b{0};

	typename t_map::const_iterator iter = params.find("mode");
	if(iter != params.end())
		b = (str_to_var<int>(iter->second) != 0);

	return b;
}


template<class t_real>
std::size_t FileTax<t_real>::GetScanCount() const
{
	if(m_dat.GetColumnCount() != 0)
		return m_dat.GetRowCount();
	return 0;
}


template<class t_real>
std::array<t_real, 5> FileTax<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("h", "k", "l", "e", i);
}


template<class t_real> std::vector<std::string> FileTax<t_real>::GetScannedVars() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColVars;
	typename t_map::const_iterator iter = params.find("def_x");
	if(iter != params.end())
		strColVars = iter->second;

	std::vector<std::string> vecVars;
	get_tokens<std::string, std::string>(strColVars, ",;", vecVars);

	// if nothing is given, default to E
	if(!vecVars.size())
		vecVars.push_back("e");

	return vecVars;
}


template<class t_real> std::string FileTax<t_real>::GetCountVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "detector";
	typename t_map::const_iterator iter = params.find("def_y");
	if(iter != params.end())
		strColCtr = iter->second;

	return strColCtr;
}


template<class t_real> std::string FileTax<t_real>::GetMonVar() const
{
	return "monitor";
}


template<class t_real>
bool FileTax<t_real>::MergeWith(const FileInstrBase<t_real>* pDat)
{
	return FileInstrBase<t_real>::MergeWith(pDat);
}


template<class t_real> std::string FileTax<t_real>::GetTitle() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("experiment");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetUser() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("users");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetLocalContact() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("local_contact");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetScanNumber() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("scan");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetSampleName() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("samplename");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetSpacegroup() const
{
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetScanCommand() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("command");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetTimestamp() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string timestamp;

	typename t_map::const_iterator iterDate = params.find("date");
	if(iterDate != params.end())
		timestamp = iterDate->second;

	typename t_map::const_iterator iterTime = params.find("time");
	if(iterTime != params.end())
	{
		if(timestamp.length() > 0)
			timestamp += ", ";
		timestamp += iterTime->second;
	}

	return timestamp;
}

}
// -----------------------------------------------------------------------------


#endif
