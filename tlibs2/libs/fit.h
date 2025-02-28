/**
 * tlibs2 -- fitting and minimisation library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012-2024
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS2_FITTER_H__
#define __TLIBS2_FITTER_H__

#if __has_include(<Minuit2/FCNBase.h>) && __has_include(<Minuit2/MnTraceObject.h>)
	#include <Minuit2/FCNBase.h>
	#include <Minuit2/MnFcn.h>
	#include <Minuit2/FunctionMinimum.h>
	#include <Minuit2/MnMigrad.h>
	#include <Minuit2/MnPrint.h>

	#define __TLIBS2_USE_MINUIT__
#else
	//#pragma message("tlibs2: Disabling Minuit library (not found).")
#endif

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <exception>

#include "expr.h"
#include "maths.h"


namespace tl2 {


// ----------------------------------------------------------------------------
// stop request handling
// ----------------------------------------------------------------------------
struct StopRequestException : public std::runtime_error
{
	StopRequestException(const char* msg = "") : runtime_error{msg}
	{}
};


class StopRequest
{
public:
	StopRequest() = default;
	virtual ~StopRequest() = default;


	StopRequest(const StopRequest& req)
		: m_stop_requested{req.m_stop_requested}
	{
	}


	const StopRequest& operator=(const StopRequest& req)
	{
		m_stop_requested = req.m_stop_requested;
		return *this;
	}


	void SetStopRequest(const bool *b)
	{
		m_stop_requested = b;
	}


	void HandleStopRequest() const
	{
		if(!m_stop_requested || !*m_stop_requested)
			return;

		// the only way to get out of an ongoing minuit operation
		throw StopRequestException("Stop requested.");
	}


private:
	const bool *m_stop_requested{};
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Minuit interface
// @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
// ----------------------------------------------------------------------------

#ifdef __TLIBS2_USE_MINUIT__
using t_real_min = std::invoke_result_t<
	decltype(&ROOT::Minuit2::MnFcn::Up), ROOT::Minuit2::MnFcn>;



template<class t_real>
class FitterFuncModel
{
public:
	virtual ~FitterFuncModel() = default;

	virtual bool SetParams(const std::vector<t_real>& params) = 0;
	virtual t_real operator()(t_real x) const = 0;
	virtual FitterFuncModel<t_real>* copy() const = 0;
};



/**
 * interface using supplied functions with a fixed number of parameters
 * num_args also includes the "x" parameter to the function, m_vals does not
 */
template<class t_real, std::size_t num_args, typename t_func>
class FitterLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	t_func m_func{};
	std::vector<t_real> m_vals{};
	bool m_separate_free_param{true};  // separate "x" from parameters (for fitter)

public:
	FitterLamFuncModel(t_func func, bool bSeparateX = true)
		: m_func{func}, m_vals{}, m_separate_free_param{bSeparateX}
	{
		m_vals.resize(m_separate_free_param ? num_args - 1 : num_args);
	}


	virtual bool SetParams(const std::vector<t_real>& params) override
	{
		for(std::size_t i = 0; i < std::min(params.size(), m_vals.size()); ++i)
			m_vals[i] = params[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		std::vector<t_real> valsWithX;
		valsWithX.reserve(num_args);

		if(m_separate_free_param)
		{
			valsWithX.push_back(x);
			for(t_real d : m_vals)
				valsWithX.push_back(d);
		}

		const std::vector<t_real> *pvals = m_separate_free_param ? &valsWithX : &m_vals;
		t_real funcval = call<num_args, t_func, t_real, std::vector>(m_func, *pvals);
		return funcval;
	}


	virtual FitterLamFuncModel* copy() const override
	{
		FitterLamFuncModel<t_real, num_args, t_func>* model =
			new FitterLamFuncModel<t_real, num_args, t_func>(m_func);

		model->m_vals = this->m_vals;
		model->m_separate_free_param = this->m_separate_free_param;

		return model;
	}
};



/**
 * interface using supplied functions with a dynamic number of parameters
 */
template<class t_real, typename t_func>
class FitterDynLamFuncModel : public FitterFuncModel<t_real>
{
protected:
	std::size_t m_num_args{1};
	t_func m_func{};
	std::vector<t_real> m_vals{};
	bool m_separate_free_param{true};  // separate "x" from parameters (for fitter)

public:
	/**
	 * num_args also includes the "x" parameter to the function, m_vals does not
	 */
	FitterDynLamFuncModel(std::size_t num_args, t_func func, bool bSeparateX = true)
	: m_num_args{num_args}, m_func{func}, m_vals{}, m_separate_free_param{bSeparateX}
	{
		m_vals.resize(m_separate_free_param ? m_num_args - 1 : m_num_args);
	}


	virtual bool SetParams(const std::vector<t_real>& params) override
	{
		for(std::size_t i = 0; i < std::min(params.size(), m_vals.size()); ++i)
			m_vals[i] = params[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		std::vector<t_real> valsWithX;
		valsWithX.reserve(m_num_args);

		if(m_separate_free_param)
		{
			valsWithX.push_back(x);
			for(t_real d : m_vals)
				valsWithX.push_back(d);
		}

		return m_func(m_separate_free_param ? valsWithX : m_vals);
	}


	virtual FitterDynLamFuncModel* copy() const override
	{
		FitterDynLamFuncModel<t_real, t_func>* model =
		new FitterDynLamFuncModel<t_real, t_func>(m_num_args, m_func);

		model->m_vals = this->m_vals;
		model->m_separate_free_param = this->m_separate_free_param;

		return model;
	}
};



/**
 * interface using supplied functions
 * num_args also includes the "x" parameter to the function, m_vals does not
 */
template<class t_real>
class FitterParsedFuncModel : public FitterFuncModel<t_real>
{
protected:
	std::string m_func;

	std::string m_xName = "x";
	const std::vector<std::string>& m_names;
	std::vector<t_real> m_vals;

	ExprParser<t_real> m_expr{};


public:
	FitterParsedFuncModel(const std::string& func, const std::string& xName,
		const std::vector<std::string>& vecNames)
		: m_func{func}, m_xName{xName}, m_names{vecNames}
	{
		if(!m_expr.parse(m_func))
			throw std::runtime_error("Could not parse function.");
	}


	virtual bool SetParams(const std::vector<t_real>& params) override
	{
		m_vals.resize(params.size());
		for(std::size_t i = 0; i < std::min(params.size(), m_vals.size()); ++i)
			m_vals[i] = params[i];
		return true;
	}


	virtual t_real operator()(t_real x = t_real(0)) const override
	{
		// copy the parsed expression to be thread safe
		ExprParser<t_real> expr = m_expr;

		// x is not used for minimiser
		if(m_xName != "")
			expr.register_var(m_xName, x);

		for(std::size_t i = 0; i < m_vals.size(); ++i)
			expr.register_var(m_names[i], m_vals[i]);

		t_real val = expr.eval();
		return val;
	}


	virtual FitterParsedFuncModel* copy() const override
	{
		return new FitterParsedFuncModel<t_real>(m_func, m_xName, m_names);
	}
};

// ----------------------------------------------------------------------------



/**
 * generic chi^2 calculation for fitting
 * @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */
template<class t_real = t_real_min>
class Chi2Function : public ROOT::Minuit2::FCNBase, public StopRequest
{
protected:
	const FitterFuncModel<t_real_min> *m_fkt = nullptr;

	std::size_t m_num_pts = 0;
	const t_real *m_x = nullptr;
	const t_real *m_y = nullptr;
	const t_real *m_dy = nullptr;

	t_real_min m_sigma = 1.;
	bool m_debug = false;


public:
	Chi2Function(const FitterFuncModel<t_real_min> *fkt = nullptr,
		std::size_t num_pts = 0, const t_real *px = nullptr,
		const t_real *py = nullptr, const t_real *pdy = nullptr)
		: m_fkt{fkt}, m_num_pts{num_pts}, m_x{px}, m_y{py}, m_dy{pdy}
	{}

	virtual ~Chi2Function() = default;


	const Chi2Function<t_real>& operator=(const Chi2Function<t_real>& other)
	{
		StopRequest::operator=(*this);

		this->m_fkt = other.m_fkt;
		this->m_x = other.m_x;
		this->m_y = other.m_y;
		this->m_dy = other.m_dy;
		this->m_sigma = other.m_sigma;
		this->m_debug = other.m_debug;

		return *this;
	}

	Chi2Function(const Chi2Function<t_real>& other)
	{
		operator=(other);
	}


	/*
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& params) const
	{
		// cannot operate on m_fkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_fkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(params);
		return tl2::chi2<t_real_min, decltype(*pfkt), const t_real*>(
			*pfkt, m_num_pts, m_x, m_y, m_dy);
	}

	virtual t_real_min Up() const override
	{
		return m_sigma*m_sigma;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& params) const override
	{
		HandleStopRequest();

		t_real_min dChi2 = chi2(params);
		if(m_debug)
			std::cerr << "Fitter: chi2 = " << dChi2 << "." << std::endl;
		return dChi2;
	}

	void SetSigma(t_real_min dSig)
	{
		m_sigma = dSig;
	}

	t_real_min GetSigma() const
	{
		return m_sigma;
	}

	void SetDebug(bool b)
	{
		m_debug = b;
	}
};



/**
 * function adaptor for minimisation
 * @see http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 */
template<class t_real = t_real_min>
class MiniFunction : public ROOT::Minuit2::FCNBase, public StopRequest
{
protected:
	const FitterFuncModel<t_real_min> *m_fkt = nullptr;
	t_real_min m_sigma = 1.;

public:
	MiniFunction(const FitterFuncModel<t_real_min>* fkt = nullptr)
		: m_fkt(fkt)
	{}

	virtual ~MiniFunction() = default;

	MiniFunction(const MiniFunction<t_real>& other)
		: m_fkt(other.m_fkt), m_sigma(other.m_sigma)
	{}

	const MiniFunction<t_real>& operator=(const MiniFunction<t_real>& other)
	{
		StopRequest::operator=(*this);

		this->m_fkt = other.m_fkt;
		this->m_sigma = other.m_sigma;

		return *this;
	}

	virtual t_real_min Up() const override
	{
		return m_sigma*m_sigma;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& params) const override
	{
		HandleStopRequest();

		// cannot operate on m_fkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<FitterFuncModel<t_real_min>> uptrFkt(m_fkt->copy());
		FitterFuncModel<t_real_min>* pfkt = uptrFkt.get();

		pfkt->SetParams(params);
		return (*pfkt)(t_real_min(0));	// "0" is an ignored dummy value here
	}

	void SetSigma(t_real_min dSig)
	{
		m_sigma = dSig;
	}

	t_real_min GetSigma() const
	{
		return m_sigma;
	}
};



// ----------------------------------------------------------------------------



/**
 * fit function to x,y,dy data points
 */
template<class t_real = t_real_min, std::size_t num_args, typename t_func>
bool fit(t_func&& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::vector<std::string>& param_names,	// size: num_args-1
	std::vector<t_real>& vals,
	std::vector<t_real>& errs,
	const std::vector<bool>* fixed = nullptr,

	bool debug = true, const bool *stop_request = nullptr)
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			std::cerr << "Fitter: No data given." << std::endl;
			return false;
		}

		// check if all params are fixed
		if(fixed && std::all_of(fixed->begin(), fixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX)
				vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY)
				vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr)
				vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterLamFuncModel<t_real_min, num_args, t_func> mod(func);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		}
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecXConverted.size(), vecXConverted.data(),
				vecYConverted.data(), vecYErrConverted.data());
		}
		chi2->SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			params.Add(param_names[param_idx],
				static_cast<t_real_min>(vals[param_idx]),
				static_cast<t_real_min>(errs[param_idx]));
			if(fixed && (*fixed)[param_idx])
				params.Fix(param_names[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool fit_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			vals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(param_names[param_idx]));
			errs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(param_names[param_idx])));
		}

		if(debug)
			std::cerr << mini << std::endl;

		return fit_valid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * fit expression to x,y,dy data points
 */
template<class t_real = t_real_min>
bool fit_expr(const std::string& func,

	const std::vector<t_real>& vecX,
	const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr,

	const std::string& x_name,
	const std::vector<std::string>& param_names,	// size: num_args-1
	std::vector<t_real>& vals,
	std::vector<t_real>& errs,
	const std::vector<bool>* fixed = nullptr,

	bool debug = true, const bool *stop_request = nullptr)
{
	try
	{
		if(!vecX.size() || !vecY.size() || !vecYErr.size())
		{
			std::cerr << "Fitter: No data given." << std::endl;
			return false;
		}

		// check if all params are fixed
		if(fixed && std::all_of(fixed->begin(), fixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		// convert vectors if value types don't match with minuit's type
		std::vector<t_real_min> vecXConverted, vecYConverted, vecYErrConverted;
		if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			vecXConverted.reserve(vecX.size());
			vecYConverted.reserve(vecY.size());
			vecYErrConverted.reserve(vecYErr.size());

			for(t_real d : vecX)
				vecXConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecY)
				vecYConverted.push_back(static_cast<t_real_min>(d));
			for(t_real d : vecYErr)
				vecYErrConverted.push_back(static_cast<t_real_min>(d));
		}

		FitterParsedFuncModel<t_real_min> mod(func, x_name, param_names);

		std::unique_ptr<Chi2Function<t_real_min>> chi2;
		if constexpr(std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecX.size(), vecX.data(), vecY.data(), vecYErr.data());
		}
		else if constexpr(!std::is_same_v<t_real, t_real_min>)
		{
			chi2 = std::make_unique<Chi2Function<t_real_min>>(
				&mod, vecXConverted.size(), vecXConverted.data(),
				vecYConverted.data(), vecYErrConverted.data());
		}
		chi2->SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			params.Add(param_names[param_idx], static_cast<t_real_min>(vals[param_idx]), static_cast<t_real_min>(errs[param_idx]));
			if(fixed && (*fixed)[param_idx])
				params.Fix(param_names[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(*chi2, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool fit_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			vals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(param_names[param_idx]));
			errs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(param_names[param_idx])));
		}

		if(debug)
			std::cerr << mini << std::endl;

		return fit_valid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum using a lambda function with fixed args
 */
template<class t_real = t_real_min, std::size_t num_args, typename t_func>
bool minimise(t_func&& func, const std::vector<std::string>& param_names,
	std::vector<t_real>& vals, std::vector<t_real>& errs,
	const std::vector<bool>* fixed = nullptr,
	const std::vector<t_real>* lower_limits = nullptr,
	const std::vector<t_real>* upper_limits = nullptr,
	bool debug = true, const bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(fixed && std::all_of(fixed->begin(), fixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterLamFuncModel<t_real_min, num_args, t_func> mod(func, false);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			params.Add(param_names[param_idx],
				static_cast<t_real_min>(vals[param_idx]),
				static_cast<t_real_min>(errs[param_idx]));
			if(lower_limits && upper_limits)
				params.SetLimits(param_names[param_idx], (*lower_limits)[param_idx], (*upper_limits)[param_idx]);
			else if(lower_limits && !upper_limits)
				params.SetLowerLimit(param_names[param_idx], (*lower_limits)[param_idx]);
			else if(upper_limits && !lower_limits)
				params.SetUpperLimit(param_names[param_idx], (*upper_limits)[param_idx]);
			if(fixed && (*fixed)[param_idx])
				params.Fix(param_names[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool minimum_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			vals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(param_names[param_idx]));
			errs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(param_names[param_idx])));
		}

		if(debug)
			std::cerr << mini << std::endl;

		return minimum_valid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum using a lambda function with variable args
 */
template<class t_real = t_real_min, typename t_func>
bool minimise_dynargs(std::size_t num_args, t_func&& func,
	const std::vector<std::string>& param_names,
	std::vector<t_real>& vals, std::vector<t_real>& errs,
	const std::vector<bool>* fixed = nullptr,
	const std::vector<t_real>* lower_limits = nullptr,
	const std::vector<t_real>* upper_limits = nullptr,
	bool debug = true, const bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(fixed && std::all_of(fixed->begin(), fixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterDynLamFuncModel<t_real_min, t_func> mod(num_args, func, false);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			params.Add(param_names[param_idx],
				static_cast<t_real_min>(vals[param_idx]),
				static_cast<t_real_min>(errs[param_idx]));
			if(lower_limits && upper_limits)
				params.SetLimits(param_names[param_idx], (*lower_limits)[param_idx], (*upper_limits)[param_idx]);
			else if(lower_limits && !upper_limits)
				params.SetLowerLimit(param_names[param_idx], (*lower_limits)[param_idx]);
			else if(upper_limits && !lower_limits)
				params.SetUpperLimit(param_names[param_idx], (*upper_limits)[param_idx]);
			if(fixed && (*fixed)[param_idx])
				params.Fix(param_names[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool minimum_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			vals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(param_names[param_idx]));
			errs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(param_names[param_idx])));
		}

		if(debug)
			std::cerr << mini << std::endl;

		return minimum_valid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}



/**
 * find function minimum for an expression
 */
template<class t_real = t_real_min>
bool minimise_expr(const std::string& func, const std::vector<std::string>& param_names,
	std::vector<t_real>& vals, std::vector<t_real>& errs,
	const std::vector<bool>* fixed = nullptr,
	bool debug = true, const bool *stop_request = nullptr)
{
	try
	{
		// check if all params are fixed
		if(fixed && std::all_of(fixed->begin(), fixed->end(),
			[](bool b) -> bool { return b; }))
			{
				std::cerr << "Fitter: All parameters are fixed." << std::endl;
				return false;
			}

		FitterParsedFuncModel<t_real_min> mod(func, "", param_names);
		MiniFunction<t_real_min> minfunc(&mod);
		minfunc.SetStopRequest(stop_request);

		ROOT::Minuit2::MnUserParameters params;
		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			params.Add(param_names[param_idx],
				static_cast<t_real_min>(vals[param_idx]),
				static_cast<t_real_min>(errs[param_idx]));
			if(fixed && (*fixed)[param_idx])
				params.Fix(param_names[param_idx]);
		}

		ROOT::Minuit2::MnMigrad migrad(minfunc, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bool minimum_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();

		for(std::size_t param_idx = 0; param_idx < param_names.size(); ++param_idx)
		{
			vals[param_idx] = static_cast<t_real>(
				mini.UserState().Value(param_names[param_idx]));
			errs[param_idx] = static_cast<t_real>(
				std::fabs(mini.UserState().Error(param_names[param_idx])));
		}

		if(debug)
			std::cerr << mini << std::endl;

		return minimum_valid;
	}
	catch(const tl2::StopRequestException&)
	{
		throw;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Fitter: " << ex.what() << std::endl;
	}

	return false;
}

#endif	// __TLIBS2_USE_MINUIT__
// ----------------------------------------------------------------------------


}

#endif
