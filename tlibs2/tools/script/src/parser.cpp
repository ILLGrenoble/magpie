/**
 * parser entry point
 * @author Tobias Weber <tweber@ill.fr>
 * @date 20-dec-19
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * matrix_calc
 * Copyright (C) 2020       Tobias WEBER (privately developed).
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

#include "ast.h"
#include "parser.h"
#include "llasm.h"
#include "printast.h"
#include "str.h"
#include "file.h"

#include <iostream>
#include <fstream>
#include <locale>

#include <boost/predef/os.h>
#include <boost/program_options.hpp>
namespace args = boost::program_options;


/**
 * Lexer message
 */
void yy::Lexer::LexerOutput(const char* str, int /*len*/)
{
	std::cerr << "Lexer output (line " << GetCurLine() << "): " << str
		<< "." << std::endl;
}


/**
 * Lexer error output
 */
void yy::Lexer::LexerError(const char* err)
{
	std::cerr << "Lexer error in line " << GetCurLine() << ": " << err
		<< "." << std::endl;
}


/**
 * Parser error output
 */
void yy::Parser::error(const std::string& err)
{
	std::cerr << "Parser error in line " << context.GetCurLine() << ": " << err
		<< "." << std::endl;
}


/**
 * call lexer from parser
 */
extern yy::Parser::symbol_type yylex(yy::ParserContext &context)
{
	return context.GetLexer().lex();
}


int main(int argc, char** argv)
{
	try
	{
		std::ios_base::sync_with_stdio(0);
		std::locale loc{};
		std::locale::global(loc);

		std::cout << "--------------------------------------------------------------------------------\n";
		std::cout << "This is the tlibs2 script compiler.\n";
		std::cout << "Author: Tobias Weber <tweber@ill.fr>, 2020.\n";
		std::cout << "Licensed under GPLv3.\n";
		std::cout << "--------------------------------------------------------------------------------\n";


		// llvm toolchain
		std::string tool_opt = "opt";
		std::string tool_bc = "llvm-as";
		std::string tool_s = "llc";
		std::string tool_o = "clang";
		std::string tool_exec = "clang";
		std::string tool_strip = "llvm-strip";


		// --------------------------------------------------------------------
		// get program arguments
		// --------------------------------------------------------------------
		std::vector<std::string> vecProgs;
		bool optimise = false;
		bool show_symbols = false;
		bool show_ast = false;
		std::string outprog;

		args::options_description arg_descr("Compiler arguments");
		arg_descr.add_options()
			("out,o", args::value(&outprog), "compiled program output")
			("optimise,O", args::bool_switch(&optimise), "optimise program")
			("symbols,s", args::bool_switch(&show_symbols), "output symbol table")
			("ast,a", args::bool_switch(&show_ast), "output syntax tree")
			("program", args::value<decltype(vecProgs)>(&vecProgs), "input program to compile");

		args::positional_options_description posarg_descr;
		posarg_descr.add("program", -1);

		args::options_description arg_descr_toolchain("Toolchain programs");
		arg_descr_toolchain.add_options()
			("tool_opt", args::value(&tool_opt), "llvm optimiser")
			("tool_bc", args::value(&tool_bc), "llvm bitcode assembler")
			("tool_bccomp", args::value(&tool_s), "llvm bitcode compiler")
			("tool_asm", args::value(&tool_o), "native assembler")
			("tool_link", args::value(&tool_exec), "native linker")
			("tool_strip", args::value(&tool_strip), "strip tool");
		arg_descr.add(arg_descr_toolchain);

		auto argparser = args::command_line_parser{argc, argv};
		argparser.style(args::command_line_style::default_style);
		argparser.options(arg_descr);
		argparser.positional(posarg_descr);

		args::variables_map mapArgs;
		auto parsedArgs = argparser.run();
		args::store(parsedArgs, mapArgs);
		args::notify(mapArgs);

		if(vecProgs.size() == 0)
		{
			std::cerr << "Please specify an input program.\n";
			std::cerr << arg_descr << std::endl;
			return 0;
		}

		const std::string& inprog = vecProgs[0];
		if(outprog == "")
			outprog = tl2::get_file_noext(tl2::get_file_nodir(inprog));
		if(outprog == "")
			outprog = "out";

		std::string outprog_ast = outprog + "_ast.xml";
		std::string outprog_syms = outprog + "_syms.txt";

		std::string outprog_3ac = outprog + ".asm";
		std::string outprog_3ac_opt = outprog + "_opt.asm";
		std::string outprog_bc = outprog + ".bc";
		std::string outprog_s = outprog + ".s";
		std::string outprog_o = outprog + ".o";
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// parse input
		// --------------------------------------------------------------------
		std::cout << "Parsing \"" << inprog << "\"..." << std::endl;

		std::ifstream ifstr{inprog};
		if(!ifstr)
		{
			std::cerr << "Cannot open \"" << inprog << "\"." << std::endl;
			return -1;
		}
		yy::ParserContext ctx{ifstr};

		// register runtime functions
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "pow", SymbolType::SCALAR, {SymbolType::SCALAR, SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "sqrt", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "cbrt", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "exp", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "exp2", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "exp10", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "log", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "log2", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "log10", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "sin", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "cos", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "tan", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "asin", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "acos", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "atan", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "atan2", SymbolType::SCALAR, {SymbolType::SCALAR, SymbolType::SCALAR}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "sinh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "cosh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "tanh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "asinh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "acosh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "atanh", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "round", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "ceil", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "floor", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "fabs", SymbolType::SCALAR, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "labs", SymbolType::INT, {SymbolType::INT}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "strlen", SymbolType::INT, {SymbolType::STRING}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "set_eps", SymbolType::VOID, {SymbolType::SCALAR}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "get_eps", SymbolType::SCALAR, {}, nullptr, nullptr, true);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "set_debug", SymbolType::VOID, {SymbolType::INT}, nullptr, nullptr, true);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putstr", SymbolType::VOID, {SymbolType::STRING}, nullptr, nullptr, false);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putflt", SymbolType::VOID, {SymbolType::SCALAR}, nullptr, nullptr, false);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "putint", SymbolType::VOID, {SymbolType::INT}, nullptr, nullptr, false);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "getflt", SymbolType::SCALAR, {SymbolType::STRING}, nullptr, nullptr, false);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "getint", SymbolType::INT, {SymbolType::STRING}, nullptr, nullptr, false);

		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "flt_to_str", SymbolType::VOID, {SymbolType::SCALAR, SymbolType::STRING, SymbolType::INT}, nullptr, nullptr, false);
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "int_to_str", SymbolType::VOID, {SymbolType::INT, SymbolType::STRING, SymbolType::INT}, nullptr, nullptr, false);

		std::vector<SymbolType> qr_rettypes{{ SymbolType::MATRIX, SymbolType::MATRIX }};
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "qr", SymbolType::COMP, {SymbolType::MATRIX}, nullptr, &qr_rettypes, true);

		std::vector<SymbolType> eigenvals_rettypes{{ SymbolType::VECTOR, SymbolType::VECTOR }};
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "eigenvals", SymbolType::COMP, {SymbolType::MATRIX}, nullptr, &eigenvals_rettypes, true);

		std::vector<SymbolType> eigenvecs_rettypes{{ SymbolType::VECTOR, SymbolType::VECTOR, SymbolType::MATRIX, SymbolType::MATRIX }};
		ctx.GetSymbols().AddFunc(ctx.GetScopeName(), "eigenvecs", SymbolType::COMP, {SymbolType::MATRIX}, nullptr, &eigenvecs_rettypes, true);


		yy::Parser parser(ctx);
		ctx.SetParser(&parser);

		int res = parser.parse();
		if(res != 0)
		{
			std::cerr << "Parser reports failure." << std::endl;
			return res;
		}

		if(show_symbols)
		{
			std::cout << "Writing symbol table to \"" << outprog_syms << "\"..." << std::endl;

			std::ofstream ostrSyms{outprog_syms};
			ostrSyms << ctx.GetSymbols() << std::endl;
		}

		if(show_ast)
		{
			std::cout << "Writing AST to \"" << outprog_ast << "\"..." << std::endl;

			std::ofstream ostrAST{outprog_ast};
			ASTPrinter printer{&ostrAST};

			ostrAST << "<ast>\n";
			auto stmts = ctx.GetStatements()->GetStatementList();
			for(auto iter=stmts.rbegin(); iter!=stmts.rend(); ++iter)
			{
				(*iter)->accept(&printer);
				ostrAST << "\n";
			}
			ostrAST << "</ast>" << std::endl;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// 3AC generation
		// --------------------------------------------------------------------
		std::cout << "Generating intermediate code: \""
			<< inprog << "\" -> \"" << outprog_3ac
			<< "\"..." << std::endl;

		std::ofstream ofstr{outprog_3ac};
		std::ostream* ostr = &ofstr;
		LLAsm llasm{&ctx.GetSymbols(), ostr};
		auto stmts = ctx.GetStatements()->GetStatementList();
		for(auto iter=stmts.rbegin(); iter!=stmts.rend(); ++iter)
		{
			(*iter)->accept(&llasm);
			(*ostr) << std::endl;
		}


		(*ostr) << "; -----------------------------------------------------------------------------\n";
		(*ostr) << "; imported external functions\n";
		(*ostr) << LLAsm::get_function_declarations(ctx.GetSymbols()) << std::endl;
		(*ostr) << "; -----------------------------------------------------------------------------\n";


		// additional runtime/startup code
		(*ostr) << "\n" << R"START(
; -----------------------------------------------------------------------------
; further imported external functions
declare i8* @llvm.stacksave()
declare void @llvm.stackrestore(i8*)

declare i8* @strncpy(i8*, i8*, i64)
declare i8* @strncat(i8*, i8*, i64)
declare i32 @strncmp(i8*, i8*, i64)

declare i32 @puts(i8*)
declare i32 @snprintf(i8*, i64, i8*, ...)
declare i32 @printf(i8*, ...)
declare i32 @scanf(i8*, ...)

declare i8* @memcpy(i8*, i8*, i64)

declare i8* @ext_heap_alloc(i64, i64)
declare void @ext_heap_free(i8*)

declare void @ext_init()
declare void @ext_deinit()
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; functions from runtime.cpp
declare double @ext_determinant(double*, i64)
declare i64 @ext_power(double*, double*, i64, i64)
declare i64 @ext_transpose(double*, double*, i64, i64)
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; constants
@__strfmt_s = constant [3 x i8] c"%s\00"
@__strfmt_lg = constant [4 x i8] c"%lg\00"
@__strfmt_ld = constant [4 x i8] c"%ld\00"
@__str_vecbegin = constant [3 x i8] c"[ \00"
@__str_vecend = constant [3 x i8] c" ]\00"
@__str_vecsep = constant [3 x i8] c", \00"
@__str_matsep = constant [3 x i8] c"; \00"
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; runtime functions

; returns 0 if flt <= eps
define double @zero_eps(double %flt)
{
	%eps = call double @get_eps()
	%fltabs = call double (double) @fabs(double %flt)

	%cond = fcmp ole double %fltabs, %eps
	br i1 %cond, label %labelIf, label %labelEnd
labelIf:
	ret double 0.
labelEnd:
	ret double %flt
}

; double -> string
define void @flt_to_str(double %flt, i8* %strptr, i64 %len)
{
	%fmtptr = bitcast [4 x i8]* @__strfmt_lg to i8*
	%theflt = call double (double) @zero_eps(double %flt)
	call i32 (i8*, i64, i8*, ...) @snprintf(i8* %strptr, i64 %len, i8* %fmtptr, double %theflt)
	ret void
}

; int -> string
define void @int_to_str(i64 %i, i8* %strptr, i64 %len)
{
	%fmtptr = bitcast [4 x i8]* @__strfmt_ld to i8*
	call i32 (i8*, i64, i8*, ...) @snprintf(i8* %strptr, i64 %len, i8* %fmtptr, i64 %i)
	ret void
}

; output a string
define void @putstr(i8* %val)
{
	call i32 (i8*) @puts(i8* %val)
	ret void
}

; output a float
define void @putflt(double %val)
{
	; convert to string
	%strval = alloca [64 x i8]
	%strvalptr = bitcast [64 x i8]* %strval to i8*
	call void @flt_to_str(double %val, i8* %strvalptr, i64 64)

	; output string
	call void (i8*) @putstr(i8* %strvalptr)
	ret void
}

; output an int
define void @putint(i64 %val)
{
	; convert to string
	%strval = alloca [64 x i8]
	%strvalptr = bitcast [64 x i8]* %strval to i8*
	call void @int_to_str(i64 %val, i8* %strvalptr, i64 64)

	; output string
	call void (i8*) @putstr(i8* %strvalptr)
	ret void
}

; input a float
define double @getflt(i8* %str)
{
	; output given string
	%fmtptr_s = bitcast [3 x i8]* @__strfmt_s to i8*
	call i32 (i8*, ...) @printf(i8* %fmtptr_s, i8* %str)

	; alloc double
	%d_ptr = alloca double

	; read double from stdin
	%fmtptr_g = bitcast [4 x i8]* @__strfmt_lg to i8*
	call i32 (i8*, ...) @scanf(i8* %fmtptr_g, double* %d_ptr)

	%d = load double, double* %d_ptr
	ret double %d
}

; input an int
define i64 @getint(i8* %str)
{
	; output given string
	%fmtptr_s = bitcast [3 x i8]* @__strfmt_s to i8*
	call i32 (i8*, ...) @printf(i8* %fmtptr_s, i8* %str)

	; alloc int
	%i_ptr = alloca i64

	; read int from stdin
	%fmtptr_ld = bitcast [4 x i8]* @__strfmt_ld to i8*
	call i32 (i8*, ...) @scanf(i8* %fmtptr_ld, i64* %i_ptr)

	%i = load i64, i64* %i_ptr
	ret i64 %i
}

; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
; main entry point for llvm
define i32 @main()
{
	call void @ext_init()

	; call entry function
	call void @start()

	call void @ext_deinit()
	ret i32 0
}
; -----------------------------------------------------------------------------
)START";

		(*ostr) << std::endl;
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// 3AC optimisation
		// --------------------------------------------------------------------
		if(optimise)
		{
			std::cout << "Optimising intermediate code: \""
				<< outprog_3ac << "\" -> \"" << outprog_3ac_opt
				<< "\"..." << std::endl;

			std::string cmd_opt = tool_opt + " -stats -S --strip-debug -o "
				+ outprog_3ac_opt + " " + outprog_3ac;
			if(std::system(cmd_opt.c_str()) != 0)
			{
				std::cerr << "Failed." << std::endl;
				return -1;
			}

			outprog_3ac = outprog_3ac_opt;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// Bitcode generation
		// --------------------------------------------------------------------
		std::cout << "Assembling bitcode: \""
			<< outprog_3ac << "\" -> \"" << outprog_bc
			<< "\"..." << std::endl;

		std::string cmd_bc = tool_bc + " -o " + outprog_bc + " " + outprog_3ac;
		if(std::system(cmd_bc.c_str()) != 0)
		{
			std::cerr << "Failed." << std::endl;
			return -1;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// Native compilation
		// --------------------------------------------------------------------
		std::cout << "Generating native assembly \""
			<< outprog_bc << "\" -> \"" << outprog_s
			<< "\"..." << std::endl;

		std::string opt_flag_s = optimise ? "-O2" : "";
		std::string cmd_s = tool_s + " " + opt_flag_s + " -o " + outprog_s + " " + outprog_bc;
		if(std::system(cmd_s.c_str()) != 0)
		{
			std::cerr << "Failed." << std::endl;
			return -1;
		}


		std::cout << "Assembling native code \"" << outprog_s << "\" -> \"" << outprog_o
			<< "\"..." << std::endl;

		std::string opt_flag_o = optimise ? "-O2" : "";
		std::string cmd_o = tool_o + " " + opt_flag_o + " -c -o " + outprog_o + " " + outprog_s;
		if(std::system(cmd_o.c_str()) != 0)
		{
			std::cerr << "Failed." << std::endl;
			return -1;
		}
		// --------------------------------------------------------------------



		// --------------------------------------------------------------------
		// Linking
		// --------------------------------------------------------------------
		std::cout << "Generating native executable \""
			<< outprog_o << "\" -> \"" << outprog
			<< "\"..." << std::endl;

		std::string opt_flag_exec = optimise ? "-O2" : "";
		std::string exec_libs = "-L. -lmcalc_rt -lm";
		std::string cmd_exec = tool_exec + " " + opt_flag_exec + " -o " + outprog + " " + outprog_o + " " + exec_libs;
		if(std::system(cmd_exec.c_str()) != 0)
		{
			std::cerr << "Failed." << std::endl;
			return -1;
		}
		// --------------------------------------------------------------------



#if defined(BOOST_OS_MACOS_AVAILABLE)
		// --------------------------------------------------------------------
		// Change linkage of runtime
		// --------------------------------------------------------------------
		std::cout << "Adjusting linkage of runtime library..." << std::endl;

		std::string cmd_rtlink = "install_name_tool -change @rpath/libmcalc_rt.dylib ./libmcalc_rt.dylib " + outprog;
		if(std::system(cmd_rtlink.c_str()) != 0)
		{
			std::cerr << "Failed." << std::endl;
			return -1;
		}
		// --------------------------------------------------------------------
#endif


		if(optimise)
		{
			std::cout << "Stripping debug symbols from \"" << outprog 
				<< "\"..." << std::endl;

			std::string cmd_strip = tool_strip + " " + outprog;
			if(std::system(cmd_strip.c_str()) != 0)
			{
				std::cerr << "Failed." << std::endl;
				return -1;
			}
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
