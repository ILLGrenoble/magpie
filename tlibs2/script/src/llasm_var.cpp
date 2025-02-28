/**
 * llvm three-address code generator -- variables
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


t_astret LLAsm::visit(const ASTVar* ast)
{
	t_astret sym = get_sym(ast->GetIdent());
	if(sym == nullptr)
		throw std::runtime_error("ASTVar: Symbol \"" + ast->GetIdent() + "\" not in symbol table.");

	std::string var = std::string{"%"} + sym->name;

	if(sym->ty == SymbolType::SCALAR || sym->ty == SymbolType::INT)
	{
		t_astret retvar = get_tmp_var(sym->ty, &sym->dims);
		std::string ty = LLAsm::get_type_name(sym->ty);
		(*m_ostr) << "%" << retvar->name << " = load "
			<< ty  << ", " << ty << "* " << var << "\n";
		return retvar;
	}
	else if(sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::MATRIX)
	{
		return sym;
	}
	else if(sym->ty == SymbolType::STRING)
	{
		return sym;
	}
	else
	{
		throw std::runtime_error("ASTVar: Invalid type for visited variable: \"" + sym->name + "\".");
	}

	return nullptr;
}


t_astret LLAsm::visit(const ASTVarDecl* ast)
{
	for(const auto& _var : ast->GetVariables())
	{
		t_astret sym = get_sym(_var);
		std::string ty = LLAsm::get_type_name(sym->ty);

		if(sym->ty == SymbolType::SCALAR || sym->ty == SymbolType::INT)
		{
			(*m_ostr) << "%" << sym->name << " = alloca " << ty << "\n";
		}
		else if(sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::MATRIX)
		{
			std::size_t dim = get_arraydim(sym);

			// allocate the array's memory
			(*m_ostr) << "%" << sym->name << " = alloca [" << dim << " x double]\n";
		}
		else if(sym->ty == SymbolType::STRING)
		{
			std::size_t dim = std::get<0>(sym->dims);

			// allocate the string's memory
			(*m_ostr) << "%" << sym->name << " = alloca [" << dim << " x i8]\n";

			// get a pointer to the string
			t_astret strptr = get_tmp_var(sym->ty);
			(*m_ostr) << "%" << strptr->name << " = getelementptr [" << dim << " x i8], ["
				<< dim << " x i8]* %" << sym->name << ", i64 0, i64 0\n";

			// set first element to zero
			(*m_ostr) << "store i8 0, i8* %"  << strptr->name << "\n";
		}
		else
		{
			throw std::runtime_error("ASTVarDecl: Invalid type in declaration: \"" + sym->name + "\".");
		}


		// optional assignment
		if(ast->GetAssignment())
			ast->GetAssignment()->accept(this);
	}

	return nullptr;
}


t_astret LLAsm::visit(const ASTAssign* ast)
{
	t_astret expr = ast->GetExpr()->accept(this);

	// multiple assignments
	if(ast->IsMultiAssign())
	{
		if(expr->ty != SymbolType::COMP)
		{
			throw std::runtime_error("ASTAssign: Need a compound symbol for multi-assignment.");
		}

		const auto& vars = ast->GetIdents();

		if(expr->elems.size() != vars.size())
		{
			std::ostringstream ostrerr;
			ostrerr << "ASTAssign: Mismatch in multi-assign size, "
				<< "expected " << vars.size() << " symbols, received "
				<< expr->elems.size() << " symbols.";
			throw std::runtime_error(ostrerr.str());
		}

		std::size_t elemidx = 0;
		for(std::size_t idx=0; idx<vars.size(); ++idx)
		{
			const SymbolPtr retsym = expr->elems[idx];

			const std::string& var = vars[idx];
			t_astret sym = get_sym(var);

			// get current memory block pointer
			t_astret varmemptr = get_tmp_var(SymbolType::STRING);
			(*m_ostr) << "%" << varmemptr->name << " = getelementptr i8, i8* %"
				<< expr->name << ", i64 " << elemidx << "\n";

			// directly read scalar value from memory block
			if(sym->ty == SymbolType::SCALAR || sym->ty == SymbolType::INT)
			{
				std::string symty = LLAsm::get_type_name(sym->ty);
				std::string retty = LLAsm::get_type_name(retsym->ty);

				t_astret varptr = get_tmp_var(sym->ty);
				(*m_ostr) << "%" << varptr->name << " = bitcast i8* %" << varmemptr->name
					<< " to " << retty << "*\n";

				t_astret _varval = get_tmp_var(sym->ty);
				(*m_ostr) << "%" << _varval->name << " = load " <<
					retty << ", " << retty << "* %" << varptr->name << "\n";

				if(!check_sym_compat(
					sym->ty, std::get<0>(sym->dims), std::get<1>(sym->dims),
					retsym->ty, std::get<0>(retsym->dims), std::get<1>(retsym->dims)))
				{
					std::ostringstream ostrErr;
					ostrErr << "ASTAssign: Multi-assignment type or dimension mismatch: ";
					ostrErr << Symbol::get_type_name(sym->ty) << "["
						<< std::get<0>(sym->dims) << ", " << std::get<1>(sym->dims)
						<< "] != " << Symbol::get_type_name(retsym->ty) << "["
						<< std::get<0>(retsym->dims) << ", " << std::get<1>(retsym->dims)
						<< "].";
					throw std::runtime_error(ostrErr.str());
				}

				// cast if needed
				_varval = convert_sym(_varval, sym->ty);

				(*m_ostr) << "store " << symty << " %" << _varval->name << ", "<< symty << "* %" << var << "\n";
			}

			// read double array from memory block
			else if(sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::MATRIX)
			{
				cp_mem_vec(varmemptr, sym, false);
			}

			// read char array from memory block
			else if(sym->ty == SymbolType::STRING)
			{
				cp_mem_str(varmemptr, sym, false);
			}

			// nested compound symbols
			else if(sym->ty == SymbolType::COMP)
			{
				cp_mem_comp(varmemptr, sym);
			}

			elemidx += get_bytesize(sym);
		}

		// free heap return value (TODO: check if it really is on the heap)
		(*m_ostr) << "call void @ext_heap_free(i8* %" << expr->name << ")\n";
	}

	// single assignment
	else
	{
		std::string var = ast->GetIdent();
		t_astret sym = get_sym(var);

		if(!check_sym_compat(
			sym->ty, std::get<0>(sym->dims), std::get<1>(sym->dims),
			expr->ty, std::get<0>(expr->dims), std::get<1>(expr->dims)))
		{
			std::ostringstream ostrErr;
			ostrErr << "ASTAssign: Assignment type or dimension mismatch: ";
			ostrErr << Symbol::get_type_name(sym->ty) << "["
				<< std::get<0>(sym->dims) << ", " << std::get<1>(sym->dims)
				<< "] != " << Symbol::get_type_name(expr->ty) << "["
				<< std::get<0>(expr->dims) << ", " << std::get<1>(expr->dims)
				<< "].";
			throw std::runtime_error(ostrErr.str());
		}

		// cast if needed
		expr = convert_sym(expr, sym->ty);

		if(expr->ty == SymbolType::SCALAR || expr->ty == SymbolType::INT)
		{
			std::string ty = LLAsm::get_type_name(expr->ty);
			(*m_ostr) << "store " << ty << " %" << expr->name << ", "<< ty << "* %" << var << "\n";
		}

		else if(sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::MATRIX)
		{
			std::size_t dimDst = get_arraydim(sym);
			std::size_t dimSrc = get_arraydim(expr);

			// copy elements in a loop
			generate_loop(0, dimDst, [this, expr, sym, dimDst, dimSrc](t_astret ctrval)
			{
				// loop statements
				t_astret elemptr_src = get_tmp_var(SymbolType::SCALAR);
				(*m_ostr) << "%" << elemptr_src->name << " = getelementptr [" << dimSrc << " x double], ["
					<< dimSrc << " x double]* %" << expr->name << ", i64 0, i64 %" << ctrval->name << "\n";
				t_astret elemptr_dst = get_tmp_var(SymbolType::SCALAR);
				(*m_ostr) << "%" << elemptr_dst->name << " = getelementptr [" << dimDst << " x double], ["
					<< dimDst << " x double]* %" << sym->name << ", i64 0, i64 %" << ctrval->name << "\n";
				t_astret elem_src = get_tmp_var(SymbolType::SCALAR);
				(*m_ostr) << "%" << elem_src->name << " = load double, double* %" << elemptr_src->name << "\n";

				(*m_ostr) << "store double %" << elem_src->name << ", double* %" << elemptr_dst->name << "\n";
			});
		}

		else if(sym->ty == SymbolType::STRING)
		{
			std::size_t src_dim = std::get<0>(expr->dims);
			std::size_t dst_dim = std::get<0>(sym->dims);
			//if(src_dim > dst_dim)	// TODO
			//	throw std::runtime_error("ASTAssign: Buffer of string \"" + sym->name + "\" is not large enough.");
			std::size_t dim = std::min(src_dim, dst_dim);


			// copy elements in a loop
			generate_loop(0, dim, [this, expr, sym, src_dim, dst_dim](t_astret ctrval)
			{
				// loop statements
				t_astret elemptr_src = get_tmp_var();
				(*m_ostr) << "%" << elemptr_src->name << " = getelementptr [" << src_dim << " x i8], ["
					<< src_dim << " x i8]* %" << expr->name << ", i64 0, i64 %" << ctrval->name << "\n";
				t_astret elemptr_dst = get_tmp_var();
				(*m_ostr) << "%" << elemptr_dst->name << " = getelementptr [" << dst_dim << " x i8], ["
					<< dst_dim << " x i8]* %" << sym->name << ", i64 0, i64 %" << ctrval->name << "\n";
				t_astret elem_src = get_tmp_var();
				(*m_ostr) << "%" << elem_src->name << " = load i8, i8* %" << elemptr_src->name << "\n";

				(*m_ostr) << "store i8 %" << elem_src->name << ", i8* %" << elemptr_dst->name << "\n";
			});
		}

		return expr;
	}

	return nullptr;
}


t_astret LLAsm::visit(const ASTNumConst<double>* ast)
{
	double val = ast->GetVal();

	t_astret retvar = get_tmp_var(SymbolType::SCALAR);
	t_astret retval = get_tmp_var(SymbolType::SCALAR);
	(*m_ostr) << "%" << retvar->name << " = alloca double\n";
	(*m_ostr) << "store double " << std::scientific << val << ", double* %" << retvar->name << "\n";
	(*m_ostr) << "%" << retval->name << " = load double, double* %" << retvar->name << "\n";

	return retval;
}


t_astret LLAsm::visit(const ASTNumConst<std::int64_t>* ast)
{
	std::int64_t val = ast->GetVal();

	t_astret retvar = get_tmp_var(SymbolType::INT);
	t_astret retval = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << retvar->name << " = alloca i64\n";
	(*m_ostr) << "store i64 " << val << ", i64* %" << retvar->name << "\n";
	(*m_ostr) << "%" << retval->name << " = load i64, i64* %" << retvar->name << "\n";

	return retval;
}


t_astret LLAsm::visit(const ASTStrConst* ast)
{
	const std::string& str = ast->GetVal();
	std::size_t dim = str.length()+1;

	std::array<std::size_t, 2> dims{{dim, 1}};
	t_astret str_mem = get_tmp_var(SymbolType::STRING, &dims);

	// allocate the string's memory
	(*m_ostr) << "%" << str_mem->name << " = alloca [" << dim << " x i8]\n";

	// set the individual chars
	for(std::size_t idx=0; idx<dim; ++idx)
	{
		t_astret ptr = get_tmp_var();
		(*m_ostr) << "%" << ptr->name << " = getelementptr [" << dim << " x i8], ["
			<< dim << " x i8]* %" << str_mem->name << ", i64 0, i64 " << idx << "\n";

		int val = (idx<dim-1) ? str[idx] : 0;
		(*m_ostr) << "store i8 " << val << ", i8* %"  << ptr->name << "\n";
	}

	return str_mem;
}


t_astret LLAsm::visit(const ASTExprList* ast)
{
	// only double arrays are handled here
	if(!ast->IsScalarArray())
	{
		throw std::runtime_error("ASTExprList: General expression list should not be directly evaluated.");
	}

	// array values and size
	const auto& lst = ast->GetList();
	std::size_t len = lst.size();
	std::array<std::size_t, 2> dims{{len, 1}};

	// allocate double array
	t_astret vec_mem = get_tmp_var(SymbolType::VECTOR, &dims);
	(*m_ostr) << "%" << vec_mem->name << " = alloca [" << len << " x double]\n";

	// set the individual array elements
	auto iter = lst.begin();
	for(std::size_t idx=0; idx<len; ++idx)
	{
		t_astret ptr = get_tmp_var();
		(*m_ostr) << "%" << ptr->name << " = getelementptr [" << len << " x double], ["
			<< len << " x double]* %" << vec_mem->name << ", i64 0, i64 " << idx << "\n";

		t_astret val = (*iter)->accept(this);
		val = convert_sym(val, SymbolType::SCALAR);

		(*m_ostr) << "store double %" << val->name << ", double* %"  << ptr->name << "\n";
		++iter;
	}

	return vec_mem;
}
