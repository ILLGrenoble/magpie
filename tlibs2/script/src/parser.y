/**
 * parser
 * @author Tobias Weber <tweber@ill.fr>
 * @date 20-dec-19
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * References:
 *   - https://github.com/westes/flex/tree/master/examples/manual
 *   - http://www.gnu.org/software/bison/manual/html_node/index.html
 *   - http://git.savannah.gnu.org/cgit/bison.git/tree/examples
 *   - https://de.wikipedia.org/wiki/LL(k)-Grammatik
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

// parser options
%skeleton "lalr1.cc"
//%glr-parser
%require "3.2"

%define api.parser.class { Parser }
%define api.namespace { yy }
%define api.value.type variant	// instead of union
%define api.token.constructor	// symbol constructors


// code for parser_impl.cpp
// before inclusion of definitions header
%code requires
{
	#include "ast.h"
	#include "sym.h"
}

// after inclusion of definitions header
%code
{
	#include "parser.h"
	#include <cmath>
	#include <cstdint>

	#define DEFAULT_STRING_SIZE 128


	std::shared_ptr<ASTVar> get_var(yy::ParserContext& context, const std::string& ident)
	{
		const Symbol* sym = context.FindScopedSymbol(ident);
		if(sym)
			++sym->refcnt;
		else
			context.GetParser()->error("Cannot find symbol \"" + ident + "\".");

		return std::make_shared<ASTVar>(ident);
	}
}


// (forward) declarations for parser_defs.h
%code requires
{
	namespace yy
	{
		class ParserContext;
	}
}


// parameter to use for parser and yylex
%param
{
	yy::ParserContext &context
}


// terminal symbols
%token<std::string> IDENT
%token<double> REAL
%token<std::int64_t> INT
%token<std::string> STRING
%token FUNC RET ASSIGN
%token SCALARDECL VECTORDECL MATRIXDECL STRINGDECL INTDECL
%token IF THEN ELSE
%token LOOP DO ENDLOOP ENDLOOPS SKIPLOOP
%token EQU NEQ GT LT GEQ LEQ
%token AND XOR OR NOT
%token PLUSASSIGN MINUSASSIGN MULTASSIGN DIVASSIGN MODASSIGN POWASSIGN
%token RANGE


// nonterminals
%type<std::shared_ptr<AST>> expr
%type<std::shared_ptr<AST>> statement
%type<std::shared_ptr<ASTStmts>> statements
%type<std::shared_ptr<ASTVarDecl>> variables
%type<std::shared_ptr<ASTArgNames>> full_argumentlist
%type<std::shared_ptr<ASTArgNames>> argumentlist
%type<std::shared_ptr<ASTArgNames>> identlist
%type<std::shared_ptr<ASTArgNames>> typelist
%type<std::shared_ptr<ASTStmts>> block
%type<std::shared_ptr<ASTFunc>> function
%type<std::shared_ptr<ASTTypeDecl>> typedecl
%type<std::shared_ptr<ASTExprList>> exprlist
%type<std::shared_ptr<AST>> opt_assign


// precedences and left/right-associativity
// see: https://en.wikipedia.org/wiki/Order_of_operations
%nonassoc RET
%left ','
%right '='
%right PLUSASSIGN MINUSASSIGN MULTASSIGN DIVASSIGN MODASSIGN POWASSIGN
%left XOR
%left OR
%left AND
%left GT LT GEQ LEQ
%left EQU NEQ
%left '+' '-'
%left '*' '/' '%'
%right NOT
%right UNARY_OP
%right '^' '\''
%left '(' '[' '{' '|'

// for the if/else r/s conflict shift "else"
// see: https://www.gnu.org/software/bison/manual/html_node/Non-Operators.html
%precedence IF THEN
%precedence ELSE

%precedence IDENT


%%
// non-terminals / grammar

/**
 * program start symbol
 */
program
	: statements[stmts]		{ context.SetStatements($stmts); }
	;


/**
 * a list of statements
 */
statements[res]
	: statement[stmt] statements[lst]	{ $lst->AddStatement($stmt, true); $res = $lst; }
	| /* epsilon */			{ $res = std::make_shared<ASTStmts>(); }
	;


/**
 * variables
 */
variables[res]
	// several variables
	: IDENT[name] ',' variables[lst] {
			std::string symName = context.AddScopedSymbol($name)->scoped_name;
			$lst->AddVariable(symName);
			$res = $lst;
		}

	// a variable
	| IDENT[name] {
			std::string symName = context.AddScopedSymbol($name)->scoped_name;
			$res = std::make_shared<ASTVarDecl>();
			$res->AddVariable(symName);
		}

	// a variable with an assignment
	| IDENT[name] '=' expr[term] {
			std::string symName = context.AddScopedSymbol($name)->scoped_name;
			$res = std::make_shared<ASTVarDecl>(std::make_shared<ASTAssign>($name, $term));
			$res->AddVariable(symName);
		}
	;


/**
 * statement
 */
statement[res]
	: expr[term] ';'	{ $res = $term; }
	| block[blk]		{ $res = $blk; }

	// function
	| function[func]	{ $res = $func;  }

	// (multiple) return(s)
	| RET exprlist[terms] ';' {
			$res = std::make_shared<ASTReturn>($terms);
		}

	// simple return
	| RET ';' {
			$res = std::make_shared<ASTReturn>(nullptr);
		}

	// variable declarations
	// scalar / double
	| SCALARDECL {
			context.PushSymType(SymbolType::SCALAR);
		}
		variables[vars] ';'	{ 
			$res = $vars;
			context.PopSymType();
		}

	// vector
	| VECTORDECL INT[dim] {
			context.PushSymType(SymbolType::VECTOR);
			context.PushSymDims(std::size_t($dim));
		}
		variables[vars] ';' {
			$res = $vars;
			context.PopSymType();
			context.PopSymDims();
		}

	// matrix
	| MATRIXDECL INT[dim1] INT[dim2] {
			context.PushSymType(SymbolType::MATRIX);
			context.PushSymDims(std::size_t($dim1), std::size_t($dim2));
		}
		variables[vars] ';' {
			$res = $vars;
			context.PopSymType();
			context.PopSymDims();
		}

	// string with default size
	| STRINGDECL {
			context.PushSymType(SymbolType::STRING);
			context.PushSymDims(std::size_t(DEFAULT_STRING_SIZE));
		}
		variables[vars] ';'	{
			$res = $vars;
			context.PopSymType();
			context.PopSymDims();
		}

	// string with a given (static) size
	| STRINGDECL INT[dim] {
			context.PushSymType(SymbolType::STRING);
			context.PushSymDims(std::size_t(std::size_t($dim)));
		}
		variables[vars] ';'	{ 
			$res = $vars; 
			context.PopSymType();
			context.PopSymDims();
		}

	// int
	| INTDECL {
			context.PushSymType(SymbolType::INT);
		}
		variables[vars] ';'	{ 
			$res = $vars; 
			context.PopSymType();
		}

	// conditional
	| IF expr[cond] THEN statement[if_stmt] {
		$res = std::make_shared<ASTCond>($cond, $if_stmt); }
	| IF expr[cond] THEN statement[if_stmt] ELSE statement[else_stmt] {
		$res = std::make_shared<ASTCond>($cond, $if_stmt, $else_stmt); }

	// loop
	| LOOP expr[cond] DO statement[stmt] {
		$res = std::make_shared<ASTLoop>($cond, $stmt); }

	// simple loop with automatic counter variable
	| LOOP IDENT[name] '=' expr[idx1] RANGE {
			auto stmts = std::make_shared<ASTStmts>();
			context.PushCachedAST(stmts);

			// int idx = idx1;
			std::string ctrname = context.AddScopedSymbol($name, SymbolType::INT)->scoped_name;
			auto vardecl = std::make_shared<ASTVarDecl>(std::make_shared<ASTAssign>($name, $idx1));
			vardecl->AddVariable(ctrname);
			stmts->AddStatement(vardecl);
		} 
		expr[idx2] DO statement[stmt] {
			auto stmts = std::static_pointer_cast<ASTStmts>(context.PopCachedAST());

			// idx <= idx2
			auto ctrvar = get_var(context, $name);
			auto cond = std::make_shared<ASTComp>(ctrvar, $idx2, ASTComp::LEQ);

			// ++idx
			auto one = std::make_shared<ASTNumConst<std::int64_t>>(1);
			auto opres = std::make_shared<ASTPlus>(ctrvar, one, 0);
			auto inc_ctr = std::make_shared<ASTAssign>($name, opres);

			auto loop = std::make_shared<ASTLoop>(cond, $stmt, inc_ctr);
			stmts->AddStatement(loop);

			$res = stmts;
		}

	| SKIPLOOP ';' {
		$res = std::make_shared<ASTLoopJump>(ASTLoopJump::SKIP); }
	| SKIPLOOP INT[num] ';' {
		$res = std::make_shared<ASTLoopJump>(ASTLoopJump::SKIP, std::size_t($num)); }
	| ENDLOOP ';' {
		$res = std::make_shared<ASTLoopJump>(ASTLoopJump::END); }
	| ENDLOOP INT[num] ';' {
		$res = std::make_shared<ASTLoopJump>(ASTLoopJump::END, std::size_t($num)); }
	| ENDLOOPS ';' {
		$res = std::make_shared<ASTLoopJump>(ASTLoopJump::ENDALL); }
	;


/**
 * function
 */
function[res]
	// single return value
	: FUNC typedecl[rettype] IDENT[ident] {
			context.EnterScope($ident);
		}
		'(' full_argumentlist[args] ')' {
			// register argument variables
			for(const auto& arg : $args->GetArgs())
			{
				Symbol* sym = context.AddScopedSymbol(std::get<0>(arg));
				sym->ty = std::get<1>(arg);
				std::get<0>(sym->dims) = std::get<2>(arg);
				std::get<1>(sym->dims) = std::get<3>(arg);
			}

			std::array<std::size_t, 2> retdims{{$rettype->GetDim(0), $rettype->GetDim(1)}};
			context.GetSymbols().AddFunc(
				context.GetScopeName(1), $ident,
				$rettype->GetType(), $args->GetArgTypes(), &retdims);
		}
			block[blk] {
			$res = std::make_shared<ASTFunc>($ident, $rettype, $args, $blk);
			context.LeaveScope($ident);
		}

	// no return value
	| FUNC IDENT[ident] {
			context.EnterScope($ident);
		}
		'(' full_argumentlist[args] ')' {
			// register argument variables
			for(const auto& arg : $args->GetArgs())
			{
				Symbol* sym = context.AddScopedSymbol(std::get<0>(arg));
				sym->ty = std::get<1>(arg);
				std::get<0>(sym->dims) = std::get<2>(arg);
				std::get<1>(sym->dims) = std::get<3>(arg);
			}

			context.GetSymbols().AddFunc(
				context.GetScopeName(1), $ident,
				SymbolType::VOID, $args->GetArgTypes());
		}
		block[blk] {
			auto rettype = std::make_shared<ASTTypeDecl>(SymbolType::VOID);
			$res = std::make_shared<ASTFunc>($ident, rettype, $args, $blk);
			context.LeaveScope($ident);
		}

	// multiple return values
	| FUNC '(' typelist[retargs] ')' IDENT[ident] {
			context.EnterScope($ident);
		}
		'(' full_argumentlist[args] ')' {
			// register argument variables
			for(const auto& arg : $args->GetArgs())
			{
				Symbol* sym = context.AddScopedSymbol(std::get<0>(arg));
				sym->ty = std::get<1>(arg);
				std::get<0>(sym->dims) = std::get<2>(arg);
				std::get<1>(sym->dims) = std::get<3>(arg);
			}

			std::vector<SymbolType> multirettypes = $retargs->GetArgTypes();
			context.GetSymbols().AddFunc(
				context.GetScopeName(1), $ident,
				SymbolType::COMP, $args->GetArgTypes(),
				nullptr, &multirettypes);
		}
		block[blk] {
			auto rettype = std::make_shared<ASTTypeDecl>(SymbolType::COMP);
			$res = std::make_shared<ASTFunc>($ident, rettype, $args, $blk, $retargs);
			context.LeaveScope($ident);
		}
	;


/**
 * declaration of variables
 */
typedecl[res]
	// scalars
	: SCALARDECL {
		$res = std::make_shared<ASTTypeDecl>(SymbolType::SCALAR);
	}

	// vectors
	| VECTORDECL INT[dim] {
		$res = std::make_shared<ASTTypeDecl>(SymbolType::VECTOR, $dim);
	}

	// matrices
	| MATRIXDECL INT[dim1] INT[dim2]
	{
		$res = std::make_shared<ASTTypeDecl>(SymbolType::MATRIX, $dim1, $dim2);
	}

	// strings with default size
	| STRINGDECL {
		$res = std::make_shared<ASTTypeDecl>(SymbolType::STRING, DEFAULT_STRING_SIZE);
	}

	// strings with given size
	| STRINGDECL INT[dim] {
		$res = std::make_shared<ASTTypeDecl>(SymbolType::STRING, $dim);
	}

	// ints
	| INTDECL {
		$res = std::make_shared<ASTTypeDecl>(SymbolType::INT);
	}
	;


full_argumentlist[res]
	: argumentlist[args]	{ $res = $args; }
	| /*epsilon*/		{ $res = std::make_shared<ASTArgNames>(); }
	;


/**
 * a comma-separated list of type names and variable identifiers
 */
argumentlist[res]
	: typedecl[ty] IDENT[argname] ',' argumentlist[lst] {
			$lst->AddArg($argname, $ty->GetType(), $ty->GetDim(0), $ty->GetDim(1));
			$res = $lst;
		}

	| typedecl[ty] IDENT[argname] {
			$res = std::make_shared<ASTArgNames>();
			$res->AddArg($argname, $ty->GetType(), $ty->GetDim(0), $ty->GetDim(1));
		}
	;


/**
 * a comma-separated list of variable identifiers
 */
identlist[res]
	: IDENT[argname] ',' identlist[lst] {
			$lst->AddArg($argname);
			$res = $lst;
		}

	| IDENT[argname] {
			$res = std::make_shared<ASTArgNames>();
			$res->AddArg($argname);
		}
	;


/**
 * a comma-separated list of type names
 */
typelist[res]
	: typedecl[ty] ',' typelist[lst] {
			$lst->AddArg("ret", $ty->GetType(), $ty->GetDim(0), $ty->GetDim(1));
			$res = $lst;
		}

	| typedecl[ty] {
			$res = std::make_shared<ASTArgNames>();
			$res->AddArg("ret", $ty->GetType(), $ty->GetDim(0), $ty->GetDim(1));
		}
	;


/**
 * a comma-separated list of expressions
 */
exprlist[res]
	: expr[num] ',' exprlist[lst] {
			$lst->AddExpr($num);
			$res = $lst;
		}
	| expr[num] {
			$res = std::make_shared<ASTExprList>();
			$res->AddExpr($num);
		}
	;


/**
 * a block of statements
 */
block[res]
	: '{' statements[stmts] '}'		{ $res = $stmts; }
	;


/**
 * expression
 */
expr[res]
	: '(' expr[term] ')'	{ $res = $term; }

	// unary expressions
	| '+' expr[term] %prec UNARY_OP		{ $res = $term; }
	| '-' expr[term] %prec UNARY_OP		{ $res = std::make_shared<ASTUMinus>($term); }
	| '|' expr[term] '|'	{ $res = std::make_shared<ASTNorm>($term); }
	| expr[term] '\''	{ $res = std::make_shared<ASTTransp>($term); }

	// unary boolean expression
	| NOT expr[term]	{ $res = std::make_shared<ASTBool>($term, ASTBool::NOT); }

	// binary expressions
	| expr[term1] '+' expr[term2]	{ $res = std::make_shared<ASTPlus>($term1, $term2, 0); }
	| expr[term1] '-' expr[term2]	{ $res = std::make_shared<ASTPlus>($term1, $term2, 1); }
	| expr[term1] '*' expr[term2]	{ $res = std::make_shared<ASTMult>($term1, $term2, 0); }
	| expr[term1] '/' expr[term2]	{ $res = std::make_shared<ASTMult>($term1, $term2, 1); }
	| expr[term1] '%' expr[term2]	{ $res = std::make_shared<ASTMod>($term1, $term2); }
	| expr[term1] '^' expr[term2]	{ $res = std::make_shared<ASTPow>($term1, $term2); }

	// binary boolean expressions
	| expr[term1] AND expr[term2]	{ $res = std::make_shared<ASTBool>($term1, $term2, ASTBool::AND); }
	| expr[term1] OR expr[term2]	{ $res = std::make_shared<ASTBool>($term1, $term2, ASTBool::OR); }
	| expr[term1] XOR expr[term2]	{ $res = std::make_shared<ASTBool>($term1, $term2, ASTBool::XOR); }

	// comparison expressions
	| expr[term1] EQU expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::EQU); }
	| expr[term1] NEQ expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::NEQ); }
	| expr[term1] GT expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::GT); }
	| expr[term1] LT expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::LT); }
	| expr[term1] GEQ expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::GEQ); }
	| expr[term1] LEQ expr[term2]	{ $res = std::make_shared<ASTComp>($term1, $term2, ASTComp::LEQ); }

	// binary expressions followed by assignment
	| IDENT[ident] PLUSASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTPlus>(var, $term, 0);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	| IDENT[ident] MINUSASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTPlus>(var, $term, 1);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	| IDENT[ident] MULTASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTMult>(var, $term, 0);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	| IDENT[ident] DIVASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTMult>(var, $term, 1);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	| IDENT[ident] MODASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTMod>(var, $term);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	| IDENT[ident] POWASSIGN expr[term] {
			auto var = get_var(context, $ident);
			auto opres = std::make_shared<ASTPow>(var, $term);
			$res = std::make_shared<ASTAssign>($ident, opres);
		}

	// constants
	| REAL[num]		{ $res = std::make_shared<ASTNumConst<double>>($num); }
	| INT[num]		{ $res = std::make_shared<ASTNumConst<std::int64_t>>($num); }
	| STRING[str]	{ $res = std::make_shared<ASTStrConst>($str); }
	| '[' exprlist[arr] ']'	{	// scalar array
			$arr->SetScalarArray(true);
			$res = $arr;
		}

	// variable
	| IDENT[ident] %prec IDENT	{
			// does the identifier name a constant?
			auto pair = context.GetConst($ident);
			if(std::get<0>(pair))
			{
				auto variant = std::get<1>(pair);
				if(std::holds_alternative<double>(variant))
					$res = std::make_shared<ASTNumConst<double>>(std::get<double>(variant));
				else if(std::holds_alternative<std::int64_t>(variant))
					$res = std::make_shared<ASTNumConst<std::int64_t>>(std::get<std::int64_t>(variant));
				else if(std::holds_alternative<std::string>(variant))
					$res = std::make_shared<ASTStrConst>(std::get<std::string>(variant));
			}

			// identifier names a variable
			else
			{
				auto var = get_var(context, $ident);
				$res = std::make_shared<ASTVar>($ident);
			}
		}

	// vector access and assignment
	| expr[term] '[' expr[idx] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any vector expression
				$res = std::make_shared<ASTArrayAccess>($term, $idx);
			}
			else
			{	// assignment of a vector element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx);
				}
			}
		}

	// vector ranged access and assignment
	| expr[term] '[' expr[idx1] RANGE expr[idx2] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any vector expression
				$res = std::make_shared<ASTArrayAccess>(
					$term, $idx1, $idx2, nullptr, nullptr, true);
			}
			else
			{	// assignment of a vector element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx1, $idx2, nullptr, nullptr, true);
				}
			}
		}

	// matrix access and assignment
	| expr[term] '[' expr[idx1] ',' expr[idx2] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any matrix expression
				$res = std::make_shared<ASTArrayAccess>($term, $idx1, $idx2);
			}
			else
			{	// assignment of a matrix element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx1, $idx2);
				}
			}
		}

	// matrix ranged access and assignment
	| expr[term] '[' expr[idx1] RANGE expr[idx2] ',' expr[idx3] RANGE expr[idx4] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any matrix expression
				$res = std::make_shared<ASTArrayAccess>(
					$term, $idx1, $idx2, $idx3, $idx4, true, true);
			}
			else
			{	// assignment of a matrix element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx1, $idx2, $idx3, $idx4, true, true);
				}
			}
		}

	// matrix mixed ranged/non-ranged access and assignment
	| expr[term] '[' expr[idx1] RANGE expr[idx2] ',' expr[idx3] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any matrix expression
				$res = std::make_shared<ASTArrayAccess>(
					$term, $idx1, $idx2, $idx3, $idx3, true, true);
			}
			else
			{	// assignment of a matrix element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx1, $idx2, $idx3, $idx3, true, true);
				}
			}
		}

	// matrix mixed non-ranged/ranged access and assignment
	| expr[term] '[' expr[idx1] ',' expr[idx3] RANGE expr[idx4] ']' opt_assign[opt_term] {
			if(!$opt_term)
			{	// array access into any matrix expression
				$res = std::make_shared<ASTArrayAccess>(
					$term, $idx1, $idx1, $idx3, $idx4, true, true);
			}
			else
			{	// assignment of a matrix element
				if($term->type() != ASTType::Var)
				{
					error("Can only assign to an l-value symbol.");
					$res = nullptr;
				}
				else
				{
					auto var = std::static_pointer_cast<ASTVar>($term);
					$res = std::make_shared<ASTArrayAssign>(
						var->GetIdent(), $opt_term, $idx1, $idx1, $idx3, $idx4, true, true);
				}
			}
		}


	// function calls
	| IDENT[ident] '(' ')' {
			const Symbol* sym = context.GetSymbols().FindSymbol($ident);
			if(sym && sym->ty == SymbolType::FUNC)
				++sym->refcnt;
			else
				error("Cannot find function \"" + $ident + "\".");

			$res = std::make_shared<ASTCall>($ident);
		}

	| IDENT[ident] '(' exprlist[args] ')' {
			const Symbol* sym = context.GetSymbols().FindSymbol($ident);
			if(sym && sym->ty == SymbolType::FUNC)
				++sym->refcnt;
			else
				error("Cannot find function \"" + $ident + "\".");

			$res = std::make_shared<ASTCall>($ident, $args);
		}


	// single assignment
	| IDENT[ident] '=' expr[term] %prec '=' {
			$res = std::make_shared<ASTAssign>($ident, $term);
		}

	// multiple assignments
	| ASSIGN identlist[idents] '=' expr[term] %prec '=' {
			$res = std::make_shared<ASTAssign>($idents->GetArgIdents(), $term);
		}
	;


/**
 * optional assignment
 */
opt_assign[res]
	: '=' expr[term]	{ $res = $term; }
	| /*epsilon*/		{ $res = nullptr; }
	;

%%
