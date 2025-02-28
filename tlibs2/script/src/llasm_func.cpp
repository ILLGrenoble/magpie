/**
 * llvm three-address code generator -- functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date apr/may-2020
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 *
 * References:
 *   - https://llvm.org/docs/LangRef.html
 *   - https://llvm.org/docs/tutorial/MyFirstLanguageFrontend/LangImpl03.html
 *   - https://llvm.org/docs/GettingStarted.html
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

#include "llasm.h"
#include <sstream>


t_astret LLAsm::visit(const ASTCall* ast)
{
	const std::string& funcname = ast->GetIdent();
	t_astret func = get_sym(funcname);

	if(func == nullptr)
		throw std::runtime_error("ASTCall: Function \"" + funcname + "\" not in symbol table.");
	if(ast->GetArgumentList().size() != func->argty.size())
		throw std::runtime_error("ASTCall: Invalid number of function parameters for \"" + funcname + "\".");


	// prepare arguments
	std::vector<t_astret> args;
	std::size_t _idx=0;
	for(const auto& curarg : ast->GetArgumentList())
	{
		t_astret arg = curarg->accept(this);

		// cast if needed
		t_astret arg_casted = convert_sym(arg, func->argty[_idx]);
		if(arg_casted->ty == SymbolType::STRING)
		{
			// string arguments are of type i8*, so use a pointer to the string's array
			t_astret strptr = get_tmp_var(arg_casted->ty, &arg_casted->dims);

			(*m_ostr) << "%" << strptr->name << " = getelementptr ["
				<< std::get<0>(arg_casted->dims) << " x i8], ["
				<< std::get<0>(arg_casted->dims) << " x i8]* %"
				<< arg_casted->name << ", i64 0, i64 0\n";

			args.push_back(strptr);
		}

		else if(arg_casted->ty == SymbolType::VECTOR || arg_casted->ty == SymbolType::MATRIX)
		{
			// array arguments are of type double*, so use a pointer to the array
			t_astret arrptr = get_tmp_var(arg_casted->ty, &arg_casted->dims);
			std::size_t dim = get_arraydim(arg_casted);

			(*m_ostr) << "%" << arrptr->name << " = getelementptr ["
				<< dim << " x double], ["
				<< dim << " x double]* %"
				<< arg_casted->name << ", i64 0, i64 0\n";

			args.push_back(arrptr);

			// add vector and matrix size arguments for external calls
			if(func->is_external)
			{
				t_astret dim1 = get_tmp_var(SymbolType::INT);
				t_astret dim1val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << dim1->name << " = alloca i64\n";
				(*m_ostr) << "store i64 " << std::get<0>(arg_casted->dims) << ", i64* %" << dim1->name << "\n";
				(*m_ostr) << "%" << dim1val->name << " = load i64, i64* %" << dim1->name << "\n";

				args.push_back(dim1val);

				if(arg_casted->ty == SymbolType::MATRIX)
				{
					t_astret dim2 = get_tmp_var(SymbolType::INT);
					t_astret dim2val = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << dim2->name << " = alloca i64\n";
					(*m_ostr) << "store i64 " << std::get<1>(arg_casted->dims) << ", i64* %" << dim2->name << "\n";
					(*m_ostr) << "%" << dim2val->name << " = load i64, i64* %" << dim2->name << "\n";

					args.push_back(dim2val);
				}
			}
		}

		else
		{
			args.push_back(arg_casted);
		}

		++_idx;
	}


	t_astret retvar = get_tmp_var(func->retty, &func->retdims);
	std::string retty = LLAsm::get_type_name(func->retty);

	// call function
	if(func->retty != SymbolType::VOID)
		(*m_ostr) << "%" << retvar->name << " = ";
	(*m_ostr) << "call " << retty << " @" << funcname << "(";
	for(std::size_t idx=0; idx<args.size(); ++idx)
	{
		(*m_ostr) << LLAsm::get_type_name(args[idx]->ty) << " %" << args[idx]->name;
		if(idx < args.size()-1)
			(*m_ostr) << ", ";
	}
	(*m_ostr) << ")\n";


	// prepare return values
	// allocate memory for local string copy
	if(func->retty == SymbolType::STRING)
	{
		t_astret symcpy = get_tmp_var(func->retty, &func->retdims);
		retvar = cp_mem_str(retvar, symcpy, true);

		// free heap return value (TODO: check if it really is on the heap)
		(*m_ostr) << "call void @ext_heap_free(i8* %" << retvar->name << ")\n";
	}

	// allocate memory for local array copy
	else if(func->retty == SymbolType::VECTOR || func->retty == SymbolType::MATRIX)
	{
		// cast mem to i8*
		t_astret memcast = get_tmp_var();
		(*m_ostr) << "%" << memcast->name << " = bitcast double* %" << retvar->name << " to i8*\n";

		t_astret symcpy = get_tmp_var(func->retty, &func->retdims);
		retvar = cp_mem_vec(memcast, symcpy, true);

		// free heap return value (TODO: check if it really is on the heap)
		(*m_ostr) << "call void @ext_heap_free(i8* %" << memcast->name << ")\n";
	}

	// multiple return values
	else if(func->retty == SymbolType::COMP)
	{
		// copy multi-return types
		const_cast<Symbol*>(retvar)->elems = func->elems;
	}

	return retvar;
}


t_astret LLAsm::visit(const ASTFunc* ast)
{
	m_funcstack.push(Func{.func = ast});
	m_curscope.push_back(ast->GetIdent());

	std::string rettype = LLAsm::get_type_name(std::get<0>(ast->GetRetType()));
	(*m_ostr) << "define " << rettype << " @" << ast->GetIdent() << "(";

	auto argnames = ast->GetArgs();
	std::size_t idx=0;
	for(const auto& [argname, argtype, dim1, dim2] : argnames)
	{
		const std::string arg = std::string{"__arg_"} + argname;
		(*m_ostr) << LLAsm::get_type_name(argtype) << " %" << arg;
		if(idx < argnames.size()-1)
			(*m_ostr) << ", ";
		++idx;
	}
	(*m_ostr) << ")\n{\n";


	// create local copies of the arguments
	for(const auto& [argname, argtype, dim1, dim2] : argnames)
	{
		const std::string arg = std::string{"__arg_"} + argname;
		std::array<std::size_t, 2> argdims{{dim1, dim2}};

		// create a local variable for each function parameter
		//t_astret symparam = get_tmp_var(argtype, &argdims, &argname);
		t_astret symparam = get_sym(argname);

		if(argtype == SymbolType::SCALAR || argtype == SymbolType::INT)
		{
			std::string ty = LLAsm::get_type_name(argtype);
			(*m_ostr) << "%" << symparam->name << " = alloca " << ty << "\n";
			(*m_ostr) << "store " << ty << " %" << arg << ", " << ty << "* %" << symparam->name << "\n";
		}
		else if(argtype == SymbolType::STRING)
		{
			// allocate memory for local string copy
			(*m_ostr) << "%" << symparam->name << " = alloca [" << std::get<0>(argdims) << " x i8]\n";

			t_astret strptr = get_tmp_var();
			(*m_ostr) << "%" << strptr->name << " = getelementptr [" << std::get<0>(argdims) << " x i8], ["
				<< std::get<0>(argdims) << " x i8]* %" << symparam->name << ", i64 0, i64 0\n";

			// copy string
			(*m_ostr) << "call i8* @strncpy(i8* %" << strptr->name << ", i8* %" << arg
				<< ", i64 " << std::get<0>(argdims) << ")\n";
		}
		else if(argtype == SymbolType::VECTOR || argtype == SymbolType::MATRIX)
		{
			std::size_t argdim = ::get_arraydim(argdims);

			// allocate memory for local array copy
			(*m_ostr) << "%" << symparam->name << " = alloca [" << argdim << " x double]\n";

			t_astret arrptr = get_tmp_var();
			(*m_ostr) << "%" << arrptr->name << " = getelementptr [" << argdim << " x double], ["
				<< argdim << " x double]* %" << symparam->name << ", i64 0, i64 0\n";

			// copy array
			t_astret arrptr_cast = get_tmp_var();
			t_astret arg_cast = get_tmp_var();

			// cast to memcpy argument pointer type
			(*m_ostr) << "%" << arrptr_cast->name << " = bitcast double* %" << arrptr->name << " to i8*\n";
			(*m_ostr) << "%" << arg_cast->name << " = bitcast double* %" << arg << " to i8*\n";

			(*m_ostr) << "call i8* @memcpy(i8* %" << arrptr_cast->name << ", i8* %" << arg_cast->name
				<< ", i64 " << argdim*sizeof(double) << ")\n";
		}
		else
		{
			throw std::runtime_error("ASTFunc: Argument \"" + argname + "\" has invalid type.");
		}
	}


	t_astret lastres = ast->GetStatements()->accept(this);


	if(std::get<0>(ast->GetRetType()) == SymbolType::VOID)
	{
		(*m_ostr) << "ret void\n";
	}
	else
	{
		// TODO: string and array pointers cannot be returned as they refer to the local stack

		// return result of last expression
		if(lastres)
			(*m_ostr) << "ret " << rettype << " %" << lastres->name << "\n";
		else
			(*m_ostr) << "ret " << rettype << " 0" << "\n";
	}

	(*m_ostr) << "}\n";
	m_curscope.pop_back();
	m_funcstack.pop();

	return nullptr;
}


t_astret LLAsm::visit(const ASTReturn* ast)
{
	const Func& actfunc = m_funcstack.top();
	const ASTFunc* thisfunc = actfunc.func;

	const auto& rets = ast->GetRets();
	if(!rets)
	{
		// no defined return value
		(*m_ostr) << "ret void\n";
		return nullptr;
	}

	const auto& retvals = rets->GetList();
	std::size_t numRets = retvals.size();

	// multiple return values
	if(numRets > 1)
	{
		const auto& thisfuncrets = thisfunc->GetRets();

		// evaluate all symbols and get total size
		std::size_t retsize = 0;
		std::vector<t_astret> retsyms;
		std::vector<std::size_t> elemindices;
		auto funcretiter = thisfuncrets.begin();
		for(const auto& retast : retvals)
		{
			t_astret retsym = retast->accept(this);

			if(!check_sym_compat(
				retsym->ty, std::get<0>(retsym->dims), std::get<1>(retsym->dims),
				std::get<1>(*funcretiter), std::get<2>(*funcretiter), std::get<3>(*funcretiter)))
			{
				std::ostringstream ostrErr;
				ostrErr << "ASTReturn: Multi-return type or dimension is incompatible with declared function return type: ";
				ostrErr << Symbol::get_type_name(retsym->ty) << "["
					<< std::get<0>(retsym->dims) << ", " << std::get<1>(retsym->dims)
					<< "] != " << Symbol::get_type_name(std::get<1>(*funcretiter)) << "["
					<< std::get<2>(*funcretiter) << ", " << std::get<3>(*funcretiter)
					<< "].";
				throw std::runtime_error(ostrErr.str());
			}
			// implicitly convert to declared function return type
			retsym = convert_sym(retsym, std::get<1>(*funcretiter));
			retsyms.push_back(retsym);

			elemindices.push_back(retsize);
			retsize += get_bytesize(retsym);
			std::advance(funcretiter, 1);
		}

		std::array<std::size_t, 2> memsize{{retsize, 1}};
		t_astret memblock = get_tmp_var(SymbolType::STRING, &memsize, nullptr, true);
		(*m_ostr) << "%" << memblock->name << " = call i8* @ext_heap_alloc(i64 "
			<< retsize << ", i64 1" << ")\n";

		// write memory block
		for(std::size_t retsymidx=0; retsymidx<retsyms.size(); ++retsymidx)
		{
			const t_astret retsym = retsyms[retsymidx];
			std::size_t elemidx = elemindices[retsymidx];

			t_astret varmemptr = get_tmp_var(SymbolType::STRING);
			(*m_ostr) << "%" << varmemptr->name << " = getelementptr i8, i8* %"
				<< memblock->name << ", i64 " << elemidx << "\n";

			// directly copy scalar value to memory
			if(retsym->ty == SymbolType::SCALAR || retsym->ty == SymbolType::INT)
			{
				t_astret varptr = get_tmp_var(retsym->ty);
				(*m_ostr) << "%" << varptr->name << " = bitcast i8* %" << varmemptr->name
					<< " to " << LLAsm::get_type_name(retsym->ty) << "*\n";
				(*m_ostr) << "store " << LLAsm::get_type_name(retsym->ty)
					<< " %" << retsym->name << ", "
					<< LLAsm::get_type_name(retsym->ty) << "* %"
					<< varptr->name << "\n";
			}

			// write string array to memory block
			else if(retsym->ty == SymbolType::STRING)
			{
				cp_str_mem(retsym, varmemptr);
			}

			// write double array to memory block
			else if(retsym->ty == SymbolType::VECTOR || retsym->ty == SymbolType::MATRIX)
			{
				cp_vec_mem(retsym, varmemptr);
			}

			// nested compound symbols
			else if(retsym->ty == SymbolType::COMP)
			{
				cp_comp_mem(retsym, varmemptr);
			}

			else
			{
				throw std::runtime_error("ASTReturn: Unhandled multi-return type \""
					+ LLAsm::get_type_name(retsym->ty) + "\".");
			}
		}

		(*m_ostr) << "ret i8* %" << memblock->name << "\n";
		return memblock;
	}

	// one return value
	else if(numRets == 1)
	{
		t_astret retsym = (*retvals.begin())->accept(this);

		if(!check_sym_compat(
			retsym->ty, std::get<0>(retsym->dims), std::get<1>(retsym->dims),
			std::get<0>(thisfunc->GetRetType()), std::get<1>(thisfunc->GetRetType()), std::get<2>(thisfunc->GetRetType())))
		{
			std::ostringstream ostrErr;
			ostrErr << "ASTReturn: Return type or dimension is incompatible with declared function return type: ";
			ostrErr << Symbol::get_type_name(retsym->ty) << "["
				<< std::get<0>(retsym->dims) << ", " << std::get<1>(retsym->dims)
				<< "] != " << Symbol::get_type_name(std::get<0>(thisfunc->GetRetType())) << "["
				<< std::get<1>(thisfunc->GetRetType()) << ", " << std::get<2>(thisfunc->GetRetType())
				<< "].";
			throw std::runtime_error(ostrErr.str());
		}

		// implicitly convert to declared function return type
		retsym = convert_sym(retsym, std::get<0>(thisfunc->GetRetType()));

		if(retsym->ty == SymbolType::SCALAR || retsym->ty == SymbolType::INT)
		{
			(*m_ostr) << "ret " << LLAsm::get_type_name(retsym->ty) << " %" << retsym->name << "\n";
			return retsym;
		}

		// string and array pointers cannot be returned directly as they refer to the local stack
		// returning a pointer to a copy on the heap instead
		else if(retsym->ty == SymbolType::STRING)
		{
			t_astret strret = cp_str_mem(retsym);
			(*m_ostr) << "ret i8* %" << strret->name << "\n";
			return strret;
		}

		// string and array pointers cannot be returned directly as they refer to the local stack
		// returning a pointer to a copy on the heap instead
		else if(retsym->ty == SymbolType::VECTOR || retsym->ty == SymbolType::MATRIX)
		{
			t_astret arrret_double = cp_vec_mem(retsym);
			(*m_ostr) << "ret double* %" << arrret_double->name << "\n";
			return arrret_double;
		}

		else
		{
			throw std::runtime_error("ASTReturn: Unhandled return type \""
				+ LLAsm::get_type_name(retsym->ty) + "\".");
		}
	}

	// no return values
	else
	{
		(*m_ostr) << "ret void\n";
	}

	return nullptr;
}
