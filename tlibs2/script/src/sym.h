/**
 * symbol table
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-apr-20
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

#ifndef __SYMTAB_H__
#define __SYMTAB_H__

#include <memory>
#include <string>
#include <unordered_map>
#include <array>
#include <iostream>
#include <iomanip>


struct Symbol;
using SymbolPtr = std::shared_ptr<Symbol>;


enum class SymbolType
{
	SCALAR,
	VECTOR,
	MATRIX,
	STRING,
	INT,
	VOID,

	COMP,	// compound
	FUNC,	// function pointer

	UNKNOWN,
};


struct Symbol
{
	std::string name{};
	std::string scoped_name{};

	SymbolType ty = SymbolType::VOID;
	std::array<std::size_t, 2> dims{{1,1}};

	// for functions
	std::vector<SymbolType> argty{{}};
	SymbolType retty = SymbolType::VOID;
	std::array<std::size_t, 2> retdims{{1,1}};

	// for compound type
	std::vector<SymbolPtr> elems{};

	bool tmp = false;		// temporary or declared variable?
	bool on_heap = false;		// heap or stack variable?
	bool is_external = false;

	mutable std::size_t refcnt = 0;	// number of reference to this symbol


	/**
	 * get the corresponding data type name
	 */
	static std::string get_type_name(SymbolType ty)
	{
		switch(ty)
		{
			case SymbolType::SCALAR: return "scalar";
			case SymbolType::VECTOR: return "vec";
			case SymbolType::MATRIX: return "mat";
			case SymbolType::STRING: return "str";
			case SymbolType::INT: return "int";
			case SymbolType::VOID: return "void";
			case SymbolType::COMP: return "comp";
			case SymbolType::FUNC: return "func";

			case SymbolType::UNKNOWN: return "unknown";
		}

		return "invalid";
	}
};


class SymTab
{
public:
	Symbol* AddSymbol(const std::string& scope,
		const std::string& name, SymbolType ty,
		const std::array<std::size_t, 2>& dims,
		bool is_temp = false, bool on_heap = false)
	{
		Symbol sym{.name = name, .scoped_name = scope+name,
			.ty = ty, .dims = dims, .elems = {},
			.tmp = is_temp, .on_heap = on_heap, .is_external = false,
			.refcnt = 0};
		auto pair = m_syms.insert_or_assign(scope+name, sym);
		return &pair.first->second;
	}


	Symbol* AddFunc(const std::string& scope,
		const std::string& name, SymbolType retty,
		const std::vector<SymbolType>& argtypes,
		const std::array<std::size_t, 2>* retdims = nullptr,
		const std::vector<SymbolType>* multirettypes = nullptr,
		bool externalfunc = false)
	{
		Symbol sym{.name = name, .scoped_name = scope+name,
			.ty = SymbolType::FUNC,
			.argty = argtypes, .retty = retty,
			.is_external = externalfunc, .refcnt = 0};
		if(retdims)
			sym.retdims = *retdims;
		if(multirettypes)
		{
			for(SymbolType ty : *multirettypes)
			{
				auto retsym = std::make_shared<Symbol>();
				retsym->ty = ty;
				sym.elems.emplace_back(retsym);
			}
		}
		auto pair = m_syms.insert_or_assign(scope+name, sym);
		return &pair.first->second;
	}


	const Symbol* FindSymbol(const std::string& name) const
	{
		auto iter = m_syms.find(name);
		if(iter == m_syms.end())
			return nullptr;
		return &iter->second;
	}


    const std::unordered_map<std::string, Symbol>& GetSymbols() const
    {
            return m_syms;
    }


	friend std::ostream& operator<<(std::ostream& ostr, const SymTab& tab)
	{
		ostr << std::left << std::setw(32) << "full name"
			<< std::left << std::setw(16) << "type"
			<< std::left << std::setw(8) << "refs"
			<< std::left << std::setw(8) << "dim1"
			<< std::left << std::setw(8) << "dim2"
			<< "\n";
		ostr << "--------------------------------------------------------------------------------\n";
		for(const auto& pair : tab.m_syms)
			ostr << std::left << std::setw(32) << pair.first
				<< std::left << std::setw(16) << Symbol::get_type_name(pair.second.ty)
				<< std::left << std::setw(8) << pair.second.refcnt
				<< std::left << std::setw(8) << std::get<0>(pair.second.dims)
				<< std::left << std::setw(8) << std::get<1>(pair.second.dims)
				<< "\n";

		return ostr;
	}


private:
	std::unordered_map<std::string, Symbol> m_syms;
};


#endif
