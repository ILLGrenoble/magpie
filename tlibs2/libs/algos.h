/**
 * tlibs2 -- algorithm library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2018-2021
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Forked 2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __TLIBS2_ALGOS_H__
#define __TLIBS2_ALGOS_H__

#include <algorithm>
#include <numeric>
#include <string>
#include <chrono>
#include <vector>

#include <boost/date_time/c_time.hpp>


namespace tl2 {

/**
 * copy algorithm with interleave
 */
template<class T1, class T2>
void copy_interleave(T1 inIter, T1 inEnd, T2 outIter,
	std::size_t interleave, std::size_t startskip)
{
	std::advance(inIter, startskip);

	while(std::distance(inIter, inEnd) > 0)
	{
		*outIter = *inIter;

		++outIter;
		std::advance(inIter, interleave);
	}
}



/**
 * count number of ocurrences of a sub-string in a string
 */
template<class t_str=std::string>
std::size_t count_occurrences(const t_str &str, const t_str &tok)
{
	std::size_t num = 0;
	std::size_t start = 0;
	const std::size_t len_tok = tok.length();

	while(true)
	{
		std::size_t idx = str.find(tok, start);
		if(idx == t_str::npos)
			break;

		++num;
		start += idx+len_tok;
	}

	return num;
}



/**
 * merge containers
 */
template<class t_cont>
t_cont arrayunion(const std::initializer_list<t_cont>& lst)
{
	t_cont contRet;
	for(const t_cont& cont : lst)
		contRet.insert(contRet.end(), cont.begin(), cont.end());
	return contRet;
}



/**
 * minimum of four numbers
 */
template<typename T=double>
T min4(T t1, T t2, T t3, T t4)
{
	T tmin = t1;
	tmin = std::min(tmin, t2);
	tmin = std::min(tmin, t3);
	tmin = std::min(tmin, t4);
	return tmin;
}



/**
 * like std::chrono::seconds/minutes/hours, but with variable type
 */
template<typename T = long >
using t_dur_secs = std::chrono::duration<T, std::ratio<1, 1>>;
template<typename T = long >
using t_dur_mins = std::chrono::duration<T, std::ratio<60, 1>>;
template<typename T = long >
using t_dur_hours = std::chrono::duration<T, std::ratio<60*60, 1>>;

template<typename T = long >
using t_dur_days = std::chrono::duration<T, std::ratio<60*60*24, 1>>;

template<typename T = long >
using t_dur_weeks = std::chrono::duration<T, std::ratio<60*60*24*7, 1>>;



/**
 * duration since epoch
 */
template<typename t_dur = std::chrono::seconds>
t_dur epoch_dur()
{
	namespace ch = std::chrono;
	return ch::duration_cast<t_dur>(ch::system_clock::now().time_since_epoch());
}



/**
 * seconds since epoch
 */
template<typename T=double>
T epoch()
{
	return epoch_dur<t_dur_secs<T>>().count();
}



/**
 * create a string representation of epoch
 */
template<typename T=double>
std::string epoch_to_str(T tSeconds, const char *pcFmt="%a %Y-%b-%d %H:%M:%S %Z")
{
	namespace ch = std::chrono;
	using boost::date_time::c_time;

	t_dur_secs<T> secs(tSeconds);
	ch::system_clock::time_point tp(ch::duration_cast<ch::seconds>(secs));

	std::time_t t = ch::system_clock::to_time_t(tp);
	std::tm tm;
	c_time::localtime(&t, &tm);

	char cTime[256];
	std::strftime(cTime, sizeof cTime, pcFmt, &tm);
	return std::string(cTime);
}



/**
 * get the permutation of indices to sort a container
 */
template<class Comp, class t_cont = std::vector<std::size_t>>
t_cont get_perm(std::size_t num_elems, Comp comp)
{
	t_cont perm(num_elems);
	std::iota(perm.begin(), perm.end(), 0);

	std::stable_sort(perm.begin(), perm.end(), comp);
	return perm;
}



/**
 * get the permutation of indices to sort a container
 */
template<class t_cont, class t_cont_perm = std::vector<std::size_t>>
t_cont_perm get_perm(const t_cont& cont)
{
	t_cont_perm perm = get_perm(
		cont.size(),
		[&cont](std::size_t idx1, std::size_t idx2) -> bool
		{
			return cont[idx1] < cont[idx2];
		});

	return perm;
}



/**
 * reorder a vector according to a permutation
 */
template<class t_vec, class t_perm = std::vector<std::size_t>>
t_vec reorder(const t_vec& vec, const t_perm& perm)
{
	t_vec vec_new;
	vec_new.reserve(vec.size());

	for(decltype(vec.size()) i = 0; i < vec.size(); ++i)
		vec_new.push_back(vec[perm[i]]);

	return vec_new;
}



/**
 * reorder a vector according to a permutation
 */
template<class t_vec, class t_func, class t_perm = std::vector<std::size_t>>
t_vec reorder(const t_func& get_elem_ptr, std::size_t N, const t_perm& perm)
{
	using t_elem = typename t_vec::value_type;

	t_vec vec_new;
	vec_new.reserve(N);

	for(std::size_t i = 0; i < N; ++i)
	{
		const t_elem* elem = get_elem_ptr(perm[i]);
		if(elem)
			vec_new.push_back(*elem);
	}

	return vec_new;
}



/*
 * count how many bits are needed for the given number
 */
template<typename T = std::size_t>
T count_needed_bits(T imax)
{
	T inum = 0;

	for(; imax != 0; imax >>= 1)
		++inum;

	return inum;
}



/**
 * reverse the bits of the given number
 */
template<typename T = std::size_t>
T bit_reverse(T imax, T inum)
{
	if(imax<2)
		return inum;

	T irev = 0;
	T ibitcnt = count_needed_bits<T>(imax)-2;

	for(T i = 1; i < imax; i <<= 1)
	{
		if(inum & i)
			irev |= (1 << ibitcnt);
		--ibitcnt;
	}

	return irev;
}



template<typename T = std::size_t,
	template<class...> class t_cont = std::vector>
t_cont<T> bit_reverse_indices(T imax)
{
	t_cont<T> vec;
	vec.reserve(imax);

	for(T i = 0; i < imax; ++i)
		vec.push_back(bit_reverse<T>(imax, i));

	return vec;
}



template<class T = double>
class Stopwatch
{
	public:
		typedef std::chrono::system_clock::time_point t_tp_sys;
		typedef std::chrono::steady_clock::time_point t_tp_st;
		typedef std::chrono::duration<T> t_dur;
		typedef std::chrono::system_clock::duration t_dur_sys;

	protected:
		t_tp_sys m_timeStart{};
		t_tp_st m_timeStart_st{}, m_timeStop_st{};

		t_dur m_dur{};
		t_dur_sys m_dur_sys{};

		T m_dDur{};

	public:
		Stopwatch() = default;
		~Stopwatch() = default;


		/**
		 * start the timer
		 */
		void start()
		{
			m_timeStart = std::chrono::system_clock::now();
			m_timeStart_st = std::chrono::steady_clock::now();
		}


		/**
		 * stop the timer
		 */
		void stop()
		{
			m_timeStop_st = std::chrono::steady_clock::now();

			m_dur = std::chrono::duration_cast<t_dur>(m_timeStop_st-m_timeStart_st);
			m_dur_sys = std::chrono::duration_cast<t_dur_sys>(m_dur);

			m_dDur = T(t_dur::period::num)/T(t_dur::period::den) * T(m_dur.count());
		}


		/**
		 * get the amount of time that passed
		 */
		T GetDur() const
		{
			return m_dDur;
		}


		/**
		 * convert a time point to a string representation
		 */
		static std::string to_str(const t_tp_sys& t)
		{
			using boost::date_time::c_time;

			std::time_t tStart = std::chrono::system_clock::to_time_t(t);
			std::tm tmStart;
			c_time::localtime(&tStart, &tmStart);

			char cTime[256];
			std::strftime(cTime, sizeof cTime, "%a %Y-%b-%d %H:%M:%S %Z", &tmStart);
			return std::string(cTime);
		}


		/**
		 * get the start time as string
		 */
		std::string GetStartTimeStr() const
		{
			return to_str(m_timeStart);
		}


		/**
		 * get the stop time as string
		 */
		std::string GetStopTimeStr() const
		{
			return to_str(m_timeStart + m_dur_sys);
		}


		/**
		 * get the estimated stop time
		 */
		t_tp_sys GetEstStopTime(T dProg) const
		{
			t_tp_st timeStop_st = std::chrono::steady_clock::now();
			t_dur dur = std::chrono::duration_cast<t_dur>(timeStop_st - m_timeStart_st);
			dur *= (T(1)/dProg);

			t_dur_sys dur_sys = std::chrono::duration_cast<t_dur_sys>(dur);
			t_tp_sys tpEnd = m_timeStart + dur_sys;
			return tpEnd;
		}


		/**
		 * get the estimated stop time as string
		 */
		std::string GetEstStopTimeStr(T dProg) const
		{
			return to_str(GetEstStopTime(dProg));
		}
};

}

#endif
