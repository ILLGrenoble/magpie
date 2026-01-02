/**
 * parser
 * @author Tobias Weber <tweber@ill.fr>
 * @date 27-may-18
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

#ifndef __PARSER_H__
#define __PARSER_H__

#undef yyFlexLexer

#include <iostream>
#include <string>
#include <array>
#include <stack>
#include <variant>
#include <unordered_map>
#include <cmath>

#include <FlexLexer.h>

#include "ast.h"
#include "sym.h"
#include "parser_defs.h"


namespace yy
{
	class ParserContext;


	/**
	 * lexer
	 */
	class Lexer : public yyFlexLexer
	{
	private:
		std::size_t m_curline = 1;

	protected:
		ParserContext* m_context = nullptr;

	public:
		Lexer() : yyFlexLexer{std::cin, std::cerr} {}
		Lexer(ParserContext* context, std::istream& istr)
			: yyFlexLexer{istr, std::cerr}, m_context(context) {}
		virtual ~Lexer() = default;

		virtual yy::Parser::symbol_type lex();

		virtual void LexerOutput(const char* str, int len) override;
		virtual void LexerError(const char* err) override;

		void IncCurLine() { ++m_curline; }
		std::size_t GetCurLine() const { return m_curline; }

	private:
		virtual int yylex() final { return -1; }
	};



	/**
 	* holds parser state
 	*/
	class ParserContext
	{
	private:
		yy::Lexer m_lex;
		yy::Parser* m_parser = nullptr;

		std::shared_ptr<ASTStmts> m_statements;
		std::stack<std::shared_ptr<AST>> m_cachedast;

		SymTab m_symbols;
		std::unordered_map<std::string, std::variant<double, std::int64_t, std::string>> m_consts
		{{
			{"pi", double(M_PI)},
		}};

		// information about currently parsed symbol
		std::vector<std::string> m_curscope;
		std::stack<SymbolType> m_symtype;
		std::stack<std::array<std::size_t, 2>> m_symdims;

	public:
		ParserContext(std::istream& istr = std::cin) : m_lex{this, istr}, m_statements{}
		{
			// defaults
			m_symtype.push(SymbolType::SCALAR);
			m_symdims.push(std::array<std::size_t, 2>{1, 1});
		}


		yy::Lexer& GetLexer() { return m_lex; }
		yy::Parser* GetParser() { return m_parser; }

		void SetParser(yy::Parser* parser) { m_parser = parser; }


		// --------------------------------------------------------------------
		void SetStatements(std::shared_ptr<ASTStmts> stmts) { m_statements = stmts; }
		const std::shared_ptr<ASTStmts> GetStatements() const { return m_statements; }
		// --------------------------------------------------------------------


		// --------------------------------------------------------------------
		void PushCachedAST(std::shared_ptr<AST> stmts) { m_cachedast.push(stmts); }
		std::shared_ptr<AST> PopCachedAST() { auto elem = m_cachedast.top(); m_cachedast.pop(); return elem; }
		// --------------------------------------------------------------------


		// --------------------------------------------------------------------
		// current function scope
		const std::vector<std::string>& GetScope() const
		{
			return m_curscope;
		}

		std::string GetScopeName(std::size_t up=0) const
		{
			std::string name;
			for(std::size_t i=0; i<m_curscope.size()-up; ++i)
				name += m_curscope[i] + "::";	// scope name separator
			return name;
		}

		void EnterScope(const std::string& name)
		{
			m_curscope.push_back(name);
		}

		void LeaveScope(const std::string& name)
		{
			const std::string& curscope = *m_curscope.rbegin();

			if(curscope != name)
			{
				std::cerr << "Error in line " << GetCurLine()
					<< ": Trying to leave scope " << name
					<< ", but the top scope is " << curscope
					<< "." << std::endl;
			}

			m_curscope.pop_back();
		}
		// --------------------------------------------------------------------


		// --------------------------------------------------------------------
		/**
		 * add a symbol with the type and dimensions given by the respective members
		 */
		Symbol* AddScopedSymbol(const std::string& name)
		{
			const std::string& scope = GetScopeName();
			return m_symbols.AddSymbol(scope, name, m_symtype.top(), m_symdims.top());
		}

		/**
		 * add a symbol with the type and dimensions given by the parameters
		 */
		Symbol* AddScopedSymbol(const std::string& name, SymbolType symtype, const std::array<std::size_t, 2>& symdims = {1, 1})
		{
			const std::string& scope = GetScopeName();
			return m_symbols.AddSymbol(scope, name, symtype, symdims);
		}

		const Symbol* FindScopedSymbol(const std::string& name) const
		{
			const std::string& scope = GetScopeName();
			return m_symbols.FindSymbol(scope + name);
		}

		const SymTab& GetSymbols() const { return m_symbols; }
		SymTab& GetSymbols() { return m_symbols; }

		/**
		  * type of current symbol
		  */
		void PushSymType(SymbolType ty) { m_symtype.push(ty); }
		SymbolType PopSymType() { SymbolType ty = m_symtype.top(); m_symtype.pop(); return ty; }


		/**
		  * dimensions of vector and matrix symbols
		  */
		void PushSymDims(std::size_t dim1, std::size_t dim2=1)
		{
			m_symdims.push(std::array<std::size_t, 2>{dim1, dim2});
		}

		std::array<std::size_t, 2> PopSymDims()
		{
			auto dims = m_symdims.top();
			m_symdims.pop();
			return dims;
		}


		std::pair<bool, std::variant<double, std::int64_t, std::string>>
		GetConst(const std::string& name) const
		{
			auto iter = m_consts.find(name);
			if(iter == m_consts.end())
				return std::make_pair(0, 0.);

			return std::make_pair(1, iter->second);
		}

		// --------------------------------------------------------------------

		std::size_t GetCurLine() const { return m_lex.GetCurLine(); }
	};
}



// yylex definition for lexer
#undef YY_DECL
#define YY_DECL yy::Parser::symbol_type yy::Lexer::lex()

// yylex function which the parser calls
extern yy::Parser::symbol_type yylex(yy::ParserContext &context);


// stop parsing
#define yyterminate() { return yy::Parser::by_type::kind_type(0); }


#endif
