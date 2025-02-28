/**
 * expression test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-mar-20
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

#define BOOST_TEST_MODULE Expr Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <complex>

#include "libs/expr.h"
#include "libs/str.h"
#include "libs/maths.h"


using t_types_real = std::tuple<double, float>;
using t_types_int = std::tuple<int, long>;
using t_types_cplx = std::tuple<std::complex<double>, std::complex<float>>;



BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_real, t_real, t_types_real)
{
	static constexpr t_real eps = 1e-6;

	for(bool codegen : { false, true })
	{
		tl2::ExprParser<t_real> parser;
		//parser.SetDebug(true);

		bool ok = parser.parse("1 + 2*3", codegen);
		auto result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 7, eps));

		ok = parser.parse("4 + 5*6", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 34, eps));

		ok = parser.parse(" - (sqrt(4)-5)^3 - 5/2 ", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 24.5, eps));

		ok = parser.parse("-cos(sin(1.23*pi))^(-1.2 + 3.2)", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, -0.6228, 1e-3));

		ok = parser.parse("-1.23e1 + 5e-4", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, -12.2995, 1e-3));

		ok = parser.parse("a = 4 + 5*6", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 34, eps));

		ok = parser.parse("a", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 34, eps));

		ok = parser.parse("b = 12 + 34*56 - 20; b + 4;", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 1900, eps));

		ok = parser.parse("x = 12 + 5; y = 4 * 2; x * y", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(result, 136, eps));
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_cplx, t_cplx, t_types_cplx)
{
	using t_real = typename t_cplx::value_type;
	static constexpr t_real eps = 1e-5;

	for(bool codegen : { false, true })
	{
		tl2::ExprParser<t_cplx> parser;

		bool ok = parser.parse("imag * imag", codegen);
		auto result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -1, eps));

		ok = parser.parse("1 + 2*3", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, 7, eps));

		ok = parser.parse("4 + 5*6", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, 34, eps));

		ok = parser.parse(" - (sqrt(4)-5)^3 - 5/2 ", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, 24.5, eps));

		ok = parser.parse("-cos(sin(1.23*pi))^(-1.2 + 3.2)", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -0.6228, 1e-3));

		ok = parser.parse("-1.23e1 + 5e-4", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -12.2995, 1e-3));

		ok = parser.parse("a = imag * imag", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -1, eps));

		ok = parser.parse("a", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -1, eps));

		ok = parser.parse("b = 5 * imag * imag; b + 1;", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(result, -4, eps));
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_func, t_real, t_types_real)
{
	static constexpr t_real eps = 1e-6;
	auto tupres = tl2::eval_expr<std::string, t_real>("\t2 + \t2*3*4\n");

	BOOST_TEST(
		std::get<0>(tupres),
		tl2::equals<t_real>(std::get<1>(tupres), 14., eps));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_int, t_int, t_types_int)
{
	for(bool codegen : { false, true })
	{
		tl2::ExprParser<t_int> parser;

		bool ok = parser.parse("1 + 2*3", codegen);
		t_int result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(result == 7);

		ok = parser.parse("4 + 5*6", codegen);
		result = parser.eval();
		BOOST_TEST(ok);
		BOOST_TEST(result == 34);
	}
}
