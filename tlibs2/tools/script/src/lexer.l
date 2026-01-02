/**
 * lexer
 * @author Tobias Weber <tweber@ill.fr>
 * @date 09-dec-19
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * References:
 *   - https://github.com/westes/flex/tree/master/examples/manual
 *   - http://www.gnu.org/software/bison/manual/html_node/index.html
 *   - http://git.savannah.gnu.org/cgit/bison.git/tree/examples
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

%option c++
%option noyywrap
%option yyclass = "yy::Lexer"


%{
	#include <string>
	#include <sstream>
	#include <iomanip>

	#include "parser.h"
	#include "str.h"
%}


white	[ \t\r]

num_bin	[01]
num_oct	[0-7]
num_dec	[0-9]
num_hex	[0-9a-fA-F]

int_bin "0b"{num_bin}+
int_dec {num_dec}+
int_hex "0x"{num_hex}+

real	{num_dec}+(\.{num_dec}*)?
realeng	{real}([Ee][+-]?{num_dec}+)?

ident	[A-Za-z_]+[A-Za-z0-9_]*


%%


"\n"		{ IncCurLine(); }
{white}+	/* whitespace */
"#".*		/* comment */

"+="		{ return yy::Parser::make_PLUSASSIGN(); }
"-="		{ return yy::Parser::make_MINUSASSIGN(); }
"*="		{ return yy::Parser::make_MULTASSIGN(); }
"/="		{ return yy::Parser::make_DIVASSIGN(); }
"%="		{ return yy::Parser::make_MODASSIGN(); }
"^="		{ return yy::Parser::make_POWASSIGN(); }

";"		{ return yytext[0]; }
"+"|"-"		{ return yytext[0]; }
"*"|"/"|"%"	{ return yytext[0]; }
"^"|"'"		{ return yytext[0]; }
"("|")"		{ return yytext[0]; }
"{"|"}"		{ return yytext[0]; }
"["|"]"		{ return yytext[0]; }
"|"|"|"		{ return yytext[0]; }
","		{ return yytext[0]; }

"~"		{ return yy::Parser::make_RANGE(); }
"=="		{ return yy::Parser::make_EQU(); }
"<>"		{ return yy::Parser::make_NEQ(); }
">="		{ return yy::Parser::make_GEQ(); }
"<="		{ return yy::Parser::make_LEQ(); }
">"		{ return yy::Parser::make_GT(); }
"<"		{ return yy::Parser::make_LT(); }
"and"		{ return yy::Parser::make_AND(); }
"xor"		{ return yy::Parser::make_XOR(); }
"or"		{ return yy::Parser::make_OR(); }
"not"		{ return yy::Parser::make_NOT(); }

"assign"	{ return yy::Parser::make_ASSIGN(); }
"="		{ return yytext[0]; }

"scalar"|"var"	{ return yy::Parser::make_SCALARDECL(); }
"vec"		{ return yy::Parser::make_VECTORDECL(); }
"mat"		{ return yy::Parser::make_MATRIXDECL(); }
"str"		{ return yy::Parser::make_STRINGDECL(); }
"int"		{ return yy::Parser::make_INTDECL(); }

"func"		{ return yy::Parser::make_FUNC(); }
"ret"		{ return yy::Parser::make_RET(); }

"if"		{ return yy::Parser::make_IF(); }
"then"		{ return yy::Parser::make_THEN(); }
"else"		{ return yy::Parser::make_ELSE(); }

"loop"		{ return yy::Parser::make_LOOP(); }
"do"		{ return yy::Parser::make_DO(); }
"endloop"	{ return yy::Parser::make_ENDLOOP(); }
"endloops"	{ return yy::Parser::make_ENDLOOPS(); }
"skiploop"	{ return yy::Parser::make_SKIPLOOP(); }


{int_bin}	{
	long dec = 0;

	std::size_t pos = 0;
	for(const char* pc=yytext+yyleng-1; pc>=yytext+2; --pc)
	{
		if(*pc != '0')
			dec += (1<<pos);
		++pos;
	}

	return yy::Parser::make_INT(dec);
}

{int_dec}	{
	return yy::Parser::make_INT(
		tl2::stoval<long, std::string>(yytext));
}

{int_hex}	{
	long l = 0;
	std::istringstream{yytext} >> std::hex >> l;
	return yy::Parser::make_INT(l);
}

{realeng}	{
	return yy::Parser::make_REAL(
		tl2::stoval<double, std::string>(yytext));
}


{ident}		{
	return yy::Parser::make_IDENT(yytext);
}


"\""[^\"]*"\""	{
	// string
	std::string str{yytext+1, yytext+yyleng-1};
	return yy::Parser::make_STRING(str); }


. {
	std::string err = "Unknown token: \"";
	err += yytext;
	err += "\".";

	LexerError(err.c_str());
}


%%
