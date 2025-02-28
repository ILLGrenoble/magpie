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

#ifndef __TLIBS2_LOADINSTR_RAW_H__
#define __TLIBS2_LOADINSTR_RAW_H__


#include "base.h"


namespace tl2 {


template<class t_real = double, class t_char = char>
class DatFile
{
	public:
		using t_str = std::basic_string<t_char>;
		using t_col = std::vector<t_real>;
		using t_dat = std::vector<t_col>;
		using t_map = std::unordered_map<t_str, t_str>;

	protected:
		bool m_bOk = false;

		t_char m_chComm = '#';            // comment
		t_str m_strSeps = {'=', ':'};     // key-value separator
		t_str m_strDatSep = {' ', '\t'};  // data separators

		t_dat m_vecCols{};
		std::size_t m_iCurLine = 0;

		t_map m_mapHdr{};
		std::vector<t_str> m_vecRawHdr{};

	protected:
		void ReadHeaderLine(const t_str& _strLine)
		{
			std::size_t iBeg = _strLine.find_first_not_of(m_chComm);
			if(iBeg == t_str::npos) return;

			t_str strLine(_strLine, iBeg, t_str::npos);
			trim(strLine);
			if(strLine.length() == 0)
				return;
			//std::cout << "Header line: " << strLine << std::endl;

			bool bInsert = true;
			std::pair<t_str, t_str> pair =
				split_first<t_str>(strLine, t_str({m_strSeps}), 1);
			if(pair.first.length()==0 && pair.second.length()==0)
				bInsert=0;

			if(bInsert)
				m_mapHdr.insert(std::move(pair));
			m_vecRawHdr.emplace_back(std::move(strLine));
		}

		void ReadDataLine(const t_str& strLine)
		{
			std::vector<t_real> vecToks;
			get_tokens<t_real, t_str>(strLine, m_strDatSep, vecToks);
			if(vecToks.size() == 0)
				return;

			if(vecToks.size() > m_vecCols.size())
			{
				if(m_vecCols.size() == 0)	// first line with data
				{
					m_vecCols.resize(vecToks.size());
				}
				else				// add zero columns
				{
					std::size_t iRest = vecToks.size()-m_vecCols.size();
					for(std::size_t iCol=0; iCol<iRest; ++iCol)
					{
						t_col vecNew(m_vecCols[0].size(), t_real(0));
						m_vecCols.emplace_back(std::move(vecNew));
					}

					std::cerr << "Data file loader: "
						<< "Too many elements in line "
						<< m_iCurLine << "." << std::endl;
				}
			}
			if(m_vecCols.size() < vecToks.size())
			{
				std::cerr << "Data file loader: "
					<< "Too few elements in line "
					<< m_iCurLine << "." << std::endl;
			}

			for(std::size_t iCol=0; iCol<vecToks.size(); ++iCol)
				m_vecCols[iCol].push_back(vecToks[iCol]);
			// fill with 0 if less tokens than columns
			for(std::size_t iCol=vecToks.size(); iCol<m_vecCols.size(); ++iCol)
				m_vecCols[iCol].push_back(t_real(0));
		}

	public:
		void Unload()
		{
			m_bOk = false;
			m_iCurLine = 0;
			m_vecRawHdr.clear();
			m_mapHdr.clear();
			m_vecCols.clear();
		}

		bool Save(std::basic_ostream<t_char>& ostr)
		{
			if(!m_bOk)
				return false;

			for(const typename t_map::value_type& val : m_mapHdr)
				ostr << m_chComm << " " <<
					val.first << " " << m_strSeps[0] << " " << val.second << "\n";

			std::size_t iRows = GetRowCount();
			std::size_t iCols = GetColumnCount();

			for(std::size_t iRow=0; iRow<iRows; ++iRow)
			{
				for(std::size_t iCol=0; iCol<iCols; ++iCol)
					ostr << std::setw(16) << m_vecCols[iCol][iRow] << m_strDatSep[0];
				ostr << "\n";
			}
			return true;
		}

		bool Save(const t_str& strFile)
		{
			std::ofstream ofstr(strFile);
			if(!ofstr)
				return false;
			return Save(ofstr);
		}

		bool Load(std::basic_istream<t_char>& istr)
		{
			Unload();
			while(!istr.eof())
			{
				t_str strLine;
				std::getline(istr, strLine);
				++m_iCurLine;
				trim<t_str>(strLine);
				if(strLine.length() == 0)
					continue;
				//std::cout << "Line: " << strLine << std::endl;

				if(strLine[0] == m_chComm)
					ReadHeaderLine(strLine);
				else
					ReadDataLine(strLine);
			}

			return true;
		}

		bool Load(const t_str& strFile)
		{
			//std::cout << "File: " << strFile << std::endl;
			std::basic_ifstream<t_char> ifstr(wstr_to_str(strFile));
			if(!ifstr)
				return false;
			std::basic_istream<t_char> *pIstr = &ifstr;
			return Load(*pIstr);
		}

		bool IsOk() const { return m_bOk; }
		operator bool() const { return IsOk(); }

		void SetCommentChar(t_char ch) { m_chComm = ch; }
		void SetSeparatorChars(const t_str& str) { m_strSeps = str; }
		void SetDataSeparators(const t_str& str) { m_strDatSep = str; };

		const t_dat& GetData() const { return m_vecCols; }
		t_dat& GetData() { return m_vecCols; }
		const t_col& GetColumn(std::size_t iCol) const { return m_vecCols[iCol]; }
		t_col& GetColumn(std::size_t iCol) { return m_vecCols[iCol]; }
		std::size_t GetColumnCount() const { return m_vecCols.size(); }
		std::size_t GetRowCount() const { return m_vecCols[0].size(); }

		const std::vector<t_str>& GetRawHeader() const { return m_vecRawHdr; }
		const t_map& GetHeader() const { return m_mapHdr; }

	public:
		DatFile() = default;
		virtual ~DatFile() = default;

		DatFile(const t_str& strFile)
		{
			Load(strFile);
		}
};


// ----------------------------------------------------------------------------


// raw data files
template<class _t_real = double>
class FileRaw : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		DatFile<t_real, char> m_dat{};
		t_vecColNames m_vecCols{};

	public:
		FileRaw() = default;
		virtual ~FileRaw() = default;

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

		virtual const t_vecDat& GetData() const override;
		virtual t_vecDat& GetData() override;
		virtual const t_vecColNames& GetColNames() const override;
		virtual const t_mapParams& GetAllParams() const override;

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;
		virtual std::string GetCountErr() const override;
		virtual std::string GetMonErr() const override;

		virtual std::string GetScanCommand() const override;

		std::string GetColNameFromParam(const std::string& paramName,
			const std::string& defaultValue) const;
};



// ----------------------------------------------------------------------------
// implementation


template<class t_real>
bool FileRaw<t_real>::Load(const char* pcFile)
{
	bool bOk = m_dat.Load(pcFile);
	m_vecCols.clear();
	for(std::size_t iCol=0; iCol<m_dat.GetColumnCount(); ++iCol)
		m_vecCols.emplace_back(var_to_str(iCol+1));
	return bOk;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileRaw<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileRaw*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileRaw<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	std::size_t iCol = str_to_var<std::size_t>(strName) - 1;
	if(iCol < m_dat.GetColumnCount())
	{
		if(pIdx)
			*pIdx = iCol;
		return m_dat.GetColumn(iCol);
	}

	static std::vector<t_real> vecNull;
	if(pIdx)
		*pIdx = m_dat.GetColumnCount();
	return vecNull;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecDat&
FileRaw<t_real>::GetData() const
{
	return m_dat.GetData();
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecDat&
FileRaw<t_real>::GetData()
{
	return m_dat.GetData();
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecColNames&
FileRaw<t_real>::GetColNames() const
{
	return m_vecCols;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_mapParams&
FileRaw<t_real>::GetAllParams() const
{
	return m_dat.GetHeader();
}


template<class t_real>
std::array<t_real,3> FileRaw<t_real>::GetSampleLattice() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	{
		typename t_map::const_iterator iter = params.find("sample_a");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sample_b");
		if(iter != params.end())
			b = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sample_c");
		if(iter != params.end())
			c = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{ a, b, c }};
}


template<class t_real>
std::array<t_real,3> FileRaw<t_real>::GetSampleAngles() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	{
		typename t_map::const_iterator iter = params.find("sample_alpha");
		if(iter != params.end())
			a = d2r(str_to_var<t_real>(iter->second));
	}
	{
		typename t_map::const_iterator iter = params.find("sample_beta");
		if(iter != params.end())
			b = d2r(str_to_var<t_real>(iter->second));
	}
	{
		typename t_map::const_iterator iter = params.find("sample_gamma");
		if(iter != params.end())
			c = d2r(str_to_var<t_real>(iter->second));
	}

	return std::array<t_real,3>{{ a, b, c }};
}


template<class t_real>
std::array<t_real,2> FileRaw<t_real>::GetMonoAnaD() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real m{0}, a{0};

	{
		typename t_map::const_iterator iter = params.find("mono_d");
		if(iter != params.end())
			m = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("ana_d");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,2>{{ m, a }};
}


template<class t_real>
std::array<bool, 3> FileRaw<t_real>::GetScatterSenses() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real m{0}, s{1}, a{0};

	{
		typename t_map::const_iterator iter = params.find("sense_m");
		if(iter != params.end())
			m = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sense_s");
		if(iter != params.end())
			s = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sense_a");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}

	return std::array<bool,3>{{ m > 0., s > 0., a > 0. }};
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetScatterPlane0() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real x{0}, y{0}, z{0};

	{
		typename t_map::const_iterator iter = params.find("orient1_x");
		if(iter != params.end())
			x = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient1_y");
		if(iter != params.end())
			y = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient1_z");
		if(iter != params.end())
			z = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{ x, y, z }};
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetScatterPlane1() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real x{0}, y{0}, z{0};

	{
		typename t_map::const_iterator iter = params.find("orient2_x");
		if(iter != params.end())
			x = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient2_y");
		if(iter != params.end())
			y = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient2_z");
		if(iter != params.end())
			z = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{ x, y, z }};
}


template<class t_real>
std::string FileRaw<t_real>::GetColNameFromParam(
	const std::string& paramName, const std::string& defaultVal) const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find(paramName);
	if(iter != params.end())
		return iter->second;

	return defaultVal;
}


template<class t_real>
std::array<t_real, 4> FileRaw<t_real>::GetPosHKLE() const
{
	// get column names
	std::string strColH = GetColNameFromParam("col_h", "1");
	std::string strColK = GetColNameFromParam("col_k", "2");
	std::string strColL = GetColNameFromParam("col_l", "3");
	std::string strColE = GetColNameFromParam("col_E", "4");

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol(strColH);
	const t_vecVals& vecK = GetCol(strColK);
	const t_vecVals& vecL = GetCol(strColL);
	const t_vecVals& vecE = GetCol(strColE);

	// get values
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileRaw<t_real>::GetKFix() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real k{0};

	typename t_map::const_iterator iter = params.find("k_fix");
	if(iter != params.end())
		k = str_to_var<t_real>(iter->second);

	return k;
}


template<class t_real>
bool FileRaw<t_real>::IsKiFixed() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool b{false};

	typename t_map::const_iterator iter = params.find("is_ki_fixed");
	if(iter != params.end())
		b = (str_to_var<int>(iter->second) != 0);

	return b;
}


template<class t_real>
std::size_t FileRaw<t_real>::GetScanCount() const
{
	if(m_dat.GetColumnCount() != 0)
		return m_dat.GetRowCount();
	return 0;
}


template<class t_real>
std::array<t_real, 5> FileRaw<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	// get column names
	std::string strColH = GetColNameFromParam("col_h", "1");
	std::string strColK = GetColNameFromParam("col_k", "2");
	std::string strColL = GetColNameFromParam("col_l", "3");
	std::string strColE = GetColNameFromParam("col_E", "4");

	return FileInstrBase<t_real>::GetScanHKLKiKf(
		strColH.c_str(), strColK.c_str(), strColL.c_str(),
		strColE.c_str(), i);
}


template<class t_real> std::vector<std::string> FileRaw<t_real>::GetScannedVars() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColVars;

	{
		typename t_map::const_iterator iter = params.find("cols_scanned");
		if(iter != params.end())
			strColVars = iter->second;
	}

	std::vector<std::string> vecVars;
	get_tokens<std::string, std::string>(strColVars, ",;", vecVars);

	// if nothing is given, default to E
	if(!vecVars.size())
		vecVars.push_back("4");

	return vecVars;
}


template<class t_real> std::string FileRaw<t_real>::GetCountVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "5";

	{
		typename t_map::const_iterator iter = params.find("col_ctr");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetMonVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "6";

	{
		typename t_map::const_iterator iter = params.find("col_mon");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetCountErr() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr;

	{
		typename t_map::const_iterator iter = params.find("col_ctr_err");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetMonErr() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr;

	{
		typename t_map::const_iterator iter = params.find("col_mon_err");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real>
bool FileRaw<t_real>::MergeWith(const FileInstrBase<t_real>* pDat)
{
	return FileInstrBase<t_real>::MergeWith(pDat);
}


template<class t_real> std::string FileRaw<t_real>::GetTitle() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetUser() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetLocalContact() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetScanNumber() const { return "0"; }
template<class t_real> std::string FileRaw<t_real>::GetSampleName() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetSpacegroup() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetScanCommand() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetTimestamp() const { return ""; }

// -----------------------------------------------------------------------------

}

#endif
