/**
 * llvm three-address code generator
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


LLAsm::LLAsm(SymTab* syms, std::ostream* ostr) : m_syms{syms}, m_ostr{ostr}
{}


t_astret LLAsm::visit(const ASTStmts* ast)
{
	t_astret lastres = nullptr;

	for(const auto& stmt : ast->GetStatementList())
		lastres = stmt->accept(this);

	return lastres;
}


t_astret LLAsm::visit(const ASTCond* ast)
{
	generate_cond([this, ast]() -> t_astret
	{
		return ast->GetCond()->accept(this);
	}, [this, ast]
	{
		ast->GetIf()->accept(this);
	}, [this, ast]
	{
		ast->GetElse()->accept(this);
	}, ast->HasElse());

	return nullptr;
}


t_astret LLAsm::visit(const ASTLoop* ast)
{
	auto footer = std::make_optional([this, ast]
	{
		if(ast->GetLoopFooter())
			ast->GetLoopFooter()->accept(this);
	});

	if(!ast->GetLoopFooter())
		footer = std::nullopt;


	generate_loop([this, ast]() -> t_astret
	{
		t_astret cond = ast->GetCond()->accept(this);
		return cond;
	}, [this, ast]
	{
		ast->GetLoopStmt()->accept(this);
	}, footer);

	return nullptr;
}


t_astret LLAsm::visit(const ASTLoopJump* ast)
{
	Func& actfunc = m_funcstack.top();
	std::size_t levels = ast->GetLevels();

	const auto& beforectr = actfunc.loopFooterLabels;
	const std::size_t beforectrsize = beforectr.size();

	const auto& endloop = actfunc.loopEndLabels;
	const std::size_t endloopsize = endloop.size();
	if(actfunc.loopEndLabels.size() == 0)
		throw std::runtime_error("No active loop.");

	switch(ast->GetKind())
	{
		case ASTLoopJump::SKIP:
		{
			if(beforectrsize == 0)
				throw std::runtime_error("No active managed loop.");

			levels = std::min(levels, beforectrsize);
			std::string labelBeforeCtr = beforectr[beforectrsize-levels];

			if(labelBeforeCtr == "")
				throw std::runtime_error("No loop label found before counter.");
			(*m_ostr) << "br label %" << labelBeforeCtr << "\n";
			break;
		}

		case ASTLoopJump::ENDALL:
		{
			levels = endloopsize;
			// no break, execute next block
		}
		case ASTLoopJump::END:
		{
			levels = std::min(levels, endloopsize);
			std::string endLabel = endloop[endloopsize-levels];

			if(endLabel == "")
				throw std::runtime_error("No loop end label found.");
			(*m_ostr) << "br label %" << endLabel << "\n";
			break;
		}
	}

	return nullptr;
}


// ----------------------------------------------------------------------------


t_astret LLAsm::get_tmp_var(SymbolType ty,
	const std::array<std::size_t, 2>* dims,
	const std::string* name, bool on_heap)
{
	std::string var;
	if(name)
		var = *name;

	// create a unique temporary name
	if(!name || var == "")
	{
		var = "__t";
		var += std::to_string(m_varCount);
		++m_varCount;
	}

	// if the symbol is already known (e.g. for parameters), update and use it
	if(name)
	{
		Symbol *sym = const_cast<Symbol*>(get_sym(*name));
		if(sym)
		{
			sym->ty = ty;
			sym->on_heap = on_heap;
			if(dims)
				sym->dims = *dims;
			return sym;
		}
	}

	if(dims)
		return m_syms->AddSymbol("", var, ty, *dims, true, on_heap);
	else
		return m_syms->AddSymbol("", var, ty, {1,1}, true, on_heap);
}


std::string LLAsm::get_label()
{
	std::string lab{"__l"};
	lab += std::to_string(m_labelCount);
	++m_labelCount;
	return lab;
}


/**
 * get the corresponding data type name
 */
std::string LLAsm::get_type_name(SymbolType ty)
{
	switch(ty)
	{
		case SymbolType::SCALAR: return "double";
		case SymbolType::VECTOR: return "double*";
		case SymbolType::MATRIX: return "double*";
		case SymbolType::STRING: return "i8*";
		case SymbolType::INT: return "i64";
		case SymbolType::VOID: return "void";
		case SymbolType::COMP: return "i8*";	// pointer to memory block
		//case SymbolType::FUNC: return "func";	// TODO: function pointers
	}

	return "invalid";
}


/**
 * get the element type for an array type
 */
SymbolType LLAsm::get_element_type(SymbolType ty)
{
	if(ty == SymbolType::VECTOR)
		return SymbolType::SCALAR;
	else if(ty == SymbolType::MATRIX)
		return SymbolType::SCALAR;
	else if(ty == SymbolType::STRING)
		return SymbolType::STRING;

	return SymbolType::UNKNOWN;
}


/**
 * output declarations for registered functions
 */
std::string LLAsm::get_function_declarations(const SymTab& symtab, bool only_externals)
{
	std::ostringstream ostr;

	for(const auto& _sym : symtab.GetSymbols())
	{
		const Symbol& sym = std::get<1>(_sym);

		if(sym.ty != SymbolType::FUNC)
			continue;

		// consider only external runtime functions
		if(only_externals && !sym.is_external)
			continue;

		ostr << "declare " << get_type_name(sym.retty)
			<< " @" << sym.name << "(";

		for(std::size_t arg=0; arg<sym.argty.size(); ++arg)
		{
 			SymbolType ty = sym.argty[arg];
			ostr << get_type_name(ty);

			// add size arguments for array types
			if(ty == SymbolType::VECTOR || ty == SymbolType::MATRIX)
				ostr << ", i64";	// first dim
			if(ty == SymbolType::MATRIX)
				ostr << ", i64";	// second dim

			if(arg < sym.argty.size()-1)
				ostr << ", ";
		}

		ostr << ")\n";
	}

	return ostr.str();
}


/**
 * get the (static) byte size of a symbol
 */
std::size_t LLAsm::get_bytesize(t_astret sym)
{
	SymbolType ty = sym->ty;

	switch(ty)
	{
		case SymbolType::SCALAR:
			return sizeof(double);
		case SymbolType::VECTOR:
			return sizeof(double)*std::get<0>(sym->dims);
		case SymbolType::MATRIX:
			return sizeof(double)*std::get<0>(sym->dims)*std::get<1>(sym->dims);
		case SymbolType::STRING:
			return sizeof(char)*std::get<0>(sym->dims);
		case SymbolType::INT:
			return sizeof(std::int64_t);
		case SymbolType::VOID:
			return 0;
		case SymbolType::COMP:
		{
			std::size_t size = 0;
			for(const SymbolPtr& thesym : sym->elems)
				size += get_bytesize(thesym.get());
			return size;
		}
		case SymbolType::FUNC:
			return sizeof(void*);
	}

	return 0;
}


/**
 * get the dimensions of an array type
 */
std::size_t LLAsm::get_arraydim(t_astret sym)
{
	return ::get_arraydim<2>(sym->dims);
}


/**
 * find the symbol with a specific name in the symbol table
 */
t_astret LLAsm::get_sym(const std::string& name) const
{
	std::string scoped_name;
	for(const std::string& scope : m_curscope)
		scoped_name += scope + "::";	// scope name separator
	scoped_name += name;

	const Symbol* sym = nullptr;
	if(m_syms)
	{
		sym = m_syms->FindSymbol(scoped_name);

		// try global scope instead
		if(!sym)
			sym = m_syms->FindSymbol(name);
	}

	if(sym==nullptr)
		throw std::runtime_error("get_sym: \"" + scoped_name + "\" does not have an associated symbol.");

	//++sym->refcnt;	// increment symbol's reference counter
	return sym;
}


/**
 * convert symbol to another type
 */
t_astret LLAsm::convert_sym(t_astret sym, SymbolType ty_to)
{
	// already the correct type
	if(sym->ty == ty_to)
		return sym;

	// re-interpret vector as matrix
	else if(ty_to == SymbolType::MATRIX && sym->ty == SymbolType::VECTOR)
		return sym;


	// scalar conversions
	if(ty_to == SymbolType::SCALAR || ty_to == SymbolType::INT)
	{
		std::string op;
		if(sym->ty == SymbolType::INT && ty_to == SymbolType::SCALAR)
			op = "sitofp";
		else if(sym->ty == SymbolType::SCALAR && ty_to == SymbolType::INT)
			op = "fptosi";

		if(op == "")
			throw std::runtime_error("Invalid scalar type conversion.");

		std::string from = LLAsm::get_type_name(sym->ty);
		std::string to = LLAsm::get_type_name(ty_to);

		t_astret var = get_tmp_var(ty_to, &sym->dims);
		(*m_ostr) << "%" << var->name << " = " << op << " " << from << "%" << sym->name << " to " << to << "\n";
		return var;
	}

	// conversions to string
	else if(ty_to == SymbolType::STRING)
	{
		// ... from int
		if(sym->ty == SymbolType::INT)
		{
			std::size_t len = 32;
			std::array<std::size_t, 2> dims{{len, 1}};
			t_astret str_mem = get_tmp_var(SymbolType::STRING, &dims);
			t_astret strptr = get_tmp_var(SymbolType::STRING, &dims);

			(*m_ostr) << "%" << str_mem->name << " = alloca [" << len << " x i8]\n";
			(*m_ostr) << "%" << strptr->name << " = getelementptr ["
				<< len << " x i8], [" << len << " x i8]* %"
				<< str_mem->name << ", i64 0, i64 0\n";

			(*m_ostr) << "call void @int_to_str(i64 %"  << sym->name
				<< ", i8* %" << strptr->name << ", i64 " << len << ")\n";
			return str_mem;
		}

		// ... from (double) scalar
		else if(sym->ty == SymbolType::SCALAR)
		{
			std::size_t len = 32;
			std::array<std::size_t, 2> dims{{len, 1}};
			t_astret str_mem = get_tmp_var(SymbolType::STRING, &dims);
			t_astret strptr = get_tmp_var(SymbolType::STRING, &dims);

			(*m_ostr) << "%" << str_mem->name << " = alloca [" << len << " x i8]\n";
			(*m_ostr) << "%" << strptr->name << " = getelementptr ["
				<< len << " x i8], [" << len << " x i8]* %"
				<< str_mem->name << ", i64 0, i64 0\n";

			(*m_ostr) << "call void @flt_to_str(double %"  << sym->name
				<< ", i8* %" << strptr->name << ", i64 " << len << ")\n";
			return str_mem;
		}

		// ... from vector or matrix
		else if(sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::MATRIX)
		{
			std::size_t num_floats = get_arraydim(sym);
			std::size_t len = 32 * num_floats;
			std::array<std::size_t, 2> dims{{len, 1}};

			t_astret str_mem = get_tmp_var(SymbolType::STRING, &dims);
			t_astret strptr = get_tmp_var(SymbolType::STRING, &dims);

			(*m_ostr) << "%" << str_mem->name << " = alloca [" << len << " x i8]\n";
			(*m_ostr) << "%" << strptr->name << " = getelementptr ["
				<< len << " x i8], [" << len << " x i8]* %"
				<< str_mem->name << ", i64 0, i64 0\n";


			// prepare "[ ", "] ", ", ", and "; " strings
			t_astret vecbegin = get_tmp_var(SymbolType::STRING);
			t_astret vecend = get_tmp_var(SymbolType::STRING);
			t_astret vecsep = get_tmp_var(SymbolType::STRING);
			t_astret matsep = nullptr;

			(*m_ostr) << "%" << vecbegin->name << " = bitcast [3 x i8]* @__str_vecbegin to i8*\n";
			(*m_ostr) << "%" << vecend->name << " = bitcast [3 x i8]* @__str_vecend to i8*\n";
			(*m_ostr) << "%" << vecsep->name << " = bitcast [3 x i8]* @__str_vecsep to i8*\n";

			if(sym->ty == SymbolType::MATRIX)
			{
				matsep = get_tmp_var(SymbolType::STRING);
				(*m_ostr) << "%" << matsep->name << " = bitcast [3 x i8]* @__str_matsep to i8*\n";
			}


			// vector start: "[ "
			(*m_ostr) << "call i8* @strncpy(i8* %" << strptr->name << ", i8* %" << vecbegin->name << ", i64 3)\n";


			// output array elements in a loop
			generate_loop(0, num_floats, [this, num_floats, sym, strptr, matsep, vecsep](t_astret ctrval)
			{	// loop statements
				// get vector/matrix element
				t_astret elemptr = get_tmp_var();
				t_astret elem = get_tmp_var();

				(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << num_floats << " x double], ["
					<< num_floats << " x double]* %" << sym->name << ", i64 0, i64 %" << ctrval->name << "\n";
				(*m_ostr) << "%" << elem->name << " = load double, double* %" << elemptr->name << "\n";


				// convert vector/matrix component to string
				std::size_t lenComp = 32;
				std::array<std::size_t, 2> dimsComp{{lenComp, 1}};
				t_astret strComp_mem = get_tmp_var(SymbolType::STRING, &dimsComp);
				t_astret strCompptr = get_tmp_var(SymbolType::STRING, &dimsComp);

				(*m_ostr) << "%" << strComp_mem->name << " = alloca [" << lenComp << " x i8]\n";
				(*m_ostr) << "%" << strCompptr->name << " = getelementptr ["
					<< lenComp << " x i8], [" << lenComp << " x i8]* %"
					<< strComp_mem->name << ", i64 0, i64 0\n";

				(*m_ostr) << "call void @flt_to_str(double %"  << elem->name
					<< ", i8* %" << strCompptr->name << ", i64 " << lenComp << ")\n";

				(*m_ostr) << "call i8* @strncat(i8* %" << strptr->name << ", i8* %"
					<< strCompptr->name << ", i64 " << lenComp << ")\n";


				// separator ", " or "; "
				if(sym->ty == SymbolType::MATRIX)
				{
					generate_cond([this, ctrval, num_floats]() -> t_astret
					{
						// don't output last "; " or ", ": i < num_floats-1
						t_astret _cond = get_tmp_var(SymbolType::INT);
						(*m_ostr) << "%" << _cond->name << " = icmp slt i64 %"
							<< ctrval->name << ", " << num_floats-1 << "\n";
						return _cond;
					}, [this, sym, ctrval, strptr, matsep, vecsep]
					{
						t_astret ctr_1 = get_tmp_var(SymbolType::INT);
						t_astret ctr_mod_cols = get_tmp_var(SymbolType::INT);

						generate_cond([this, sym, ctrval, ctr_1, ctr_mod_cols]() -> t_astret
						{
							//output ';' or ','?: ((i+1) % std::get<1>(sym->dims)) == 0
							t_astret cond_M_or_v = get_tmp_var(SymbolType::INT);
							(*m_ostr) << "%" << ctr_1->name << " = add i64 %" << ctrval->name << ", 1\n";
							(*m_ostr) << "%" << ctr_mod_cols->name << " = srem i64 %" << ctr_1->name
								<< ", " << std::get<0>(sym->dims) << "\n";
							(*m_ostr) << "%" << cond_M_or_v->name << " = icmp eq i64 %"
								<< ctr_mod_cols->name << ", 0\n";
							return cond_M_or_v;
						}, [this, strptr, matsep]
						{
							// if: matrix separator case
							(*m_ostr) << "call i8* @strncat(i8* %" << strptr->name
								<< ", i8* %" << matsep->name << ", i64 3)\n";
						}, [this, strptr, vecsep]
						{
							// else: vector separator case
							(*m_ostr) << "call i8* @strncat(i8* %" << strptr->name
								<< ", i8* %" << vecsep->name << ", i64 3)\n";

						}, true);
					}, []{}, false);
				}

				else if(sym->ty == SymbolType::VECTOR)
				{
					generate_cond([this, ctrval, num_floats]() -> t_astret
					{
						// don't output last ", ": i < num_floats-1
						t_astret _cond = get_tmp_var(SymbolType::INT);
						(*m_ostr) << "%" << _cond->name << " = icmp slt i64 %"
							<< ctrval->name << ", " << num_floats-1 << "\n";
						return _cond;
					}, [this, strptr, vecsep]
					{
						(*m_ostr) << "call i8* @strncat(i8* %" << strptr->name
							<< ", i8* %" << vecsep->name << ", i64 3)\n";
					}, []{}, false);
				}
			});

			// vector end: " ]"
			(*m_ostr) << "call i8* @strncat(i8* %" << strptr->name << ", i8* %" << vecend->name << ", i64 3)\n";

			return str_mem;
		}
	}

	// unknown conversion
	else
	{
		throw std::runtime_error("Invalid type conversion.");
	}

	return nullptr;
}


bool LLAsm::check_sym_compat(
	SymbolType ty1, std::size_t dim1_1, std::size_t dim1_2,
	SymbolType ty2, std::size_t dim2_1, std::size_t dim2_2)
{
	// allow truncating assignments of the form: vec 2 v = [1, 2, 3, 4];
	// needed for ranged array access
	// TODO: remove as soon as ranged access uses dynamic array size
	if(ty1==SymbolType::VECTOR && ty2==SymbolType::VECTOR)
		return dim1_1 <= dim2_1 && (dim2_2==0 || dim2_2==1);

	// allow assignments of the form: mat 2 2 A = [1, 2, 3, 4];
	// also allow truncated assignments ("<="), TODO: remove in the future
	if(ty1==SymbolType::MATRIX && ty2==SymbolType::VECTOR)
		return dim1_1*dim1_2 <= dim2_1 && (dim2_2==0 || dim2_2==1);

	if(ty1==SymbolType::MATRIX && ty2==SymbolType::MATRIX)
		return dim1_1==dim2_1 && dim1_2==dim2_2;
	if(ty1==SymbolType::VECTOR && ty2==SymbolType::VECTOR)
		return dim1_1==dim2_1;
	if(ty1==SymbolType::MATRIX && ty2==SymbolType::VECTOR)
		return false;
	if(ty1==SymbolType::VECTOR && ty2==SymbolType::MATRIX)
		return false;

	return true;
}


/**
 * product between a scalar and a matrix (or vector)
 */
t_astret LLAsm::scalar_matrix_prod(t_astret scalar, t_astret matrix, bool mul_or_div)
{
	scalar = convert_sym(scalar, SymbolType::SCALAR);

	std::size_t dim = std::get<0>(matrix->dims);
	if(matrix->ty == SymbolType::MATRIX)
		dim *= std::get<1>(matrix->dims);

	// allocate double array for result
	t_astret vec_mem = get_tmp_var(matrix->ty, &matrix->dims);
	(*m_ostr) << "%" << vec_mem->name << " = alloca [" << dim << " x double]\n";

	// copy matrix elements in a loop
	generate_loop(0, dim, [this, dim, matrix, scalar, vec_mem, mul_or_div](t_astret ctrval)
	{
		// loop statements
		t_astret elemptr_src = get_tmp_var();
		t_astret elem_src = get_tmp_var();
		(*m_ostr) << "%" << elemptr_src->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << matrix->name << ", i64 0, i64 %" << ctrval->name << "\n";
		(*m_ostr) << "%" << elem_src->name << " = load double, double* %" << elemptr_src->name << "\n";

		t_astret elemptr_dst = get_tmp_var();
		(*m_ostr) << "%" << elemptr_dst->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << vec_mem->name << ", i64 0, i64 %" << ctrval->name << "\n";

		// multiply element value with scalar
		t_astret elem_dst = get_tmp_var(SymbolType::SCALAR);

		if(mul_or_div)
		{
			(*m_ostr) << "%" << elem_dst->name << " = fmul double %"
				<< elem_src->name << ", %" << scalar->name << "\n";
		}
		else
		{
			(*m_ostr) << "%" << elem_dst->name << " = fdiv double %"
				<< elem_src->name << ", %" << scalar->name << "\n";
		}

		// save result in array
		(*m_ostr) << "store double %" << elem_dst->name << ", double* %" << elemptr_dst->name << "\n";
	});

	return vec_mem;
}


/**
 * copy the memory of a compound symbol
 */
t_astret LLAsm::cp_comp_mem(t_astret sym, t_astret mem)
{
	std::size_t len = get_bytesize(sym);

	(*m_ostr) << "call i8* @memcpy(i8* %" << mem->name << ", i8* %" << sym->name
		<< ", i64 " << len << ")\n";

	return mem;
}


/**
 * copy a vector into a memory block
 */
t_astret LLAsm::cp_vec_mem(t_astret sym, t_astret mem)
{
	std::size_t dim = std::get<0>(sym->dims);
	if(sym->ty == SymbolType::MATRIX)
		dim *= std::get<1>(sym->dims);

	t_astret termptr = get_tmp_var(sym->ty);
	(*m_ostr) << "%" << termptr->name << " = bitcast [" << dim << " x double]* %" << sym->name << " to i8*\n";

	// if no memory block is given, allocate one
	if(!mem)
	{
		mem = get_tmp_var(sym->ty, &sym->dims, nullptr, true);
		(*m_ostr) << "%" << mem->name << " = call i8* @ext_heap_alloc(i64 "
			<< dim << ", i64 " << sizeof(double) << ")\n";
	}

	(*m_ostr) << "call i8* @memcpy(i8* %" << mem->name << ", i8* %" << termptr->name
		<< ", i64 " << dim*sizeof(double) << ")\n";

	// cast to the actual double*
	t_astret mem_double = get_tmp_var(sym->ty, &sym->dims);
	(*m_ostr) << "%" << mem_double->name << " = bitcast i8* %" << mem->name << " to double*\n";

	return mem_double;
}


/**
 * copy a string into a memory block
 */
t_astret LLAsm::cp_str_mem(t_astret sym, t_astret mem)
{
	std::size_t dim = std::get<0>(sym->dims);

	t_astret termptr = get_tmp_var(SymbolType::STRING);
	(*m_ostr) << "%" << termptr->name << " = bitcast [" << dim << " x i8]* %" << sym->name << " to i8*\n";

	t_astret strretlen = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << strretlen->name << " = call i64 @strlen(i8* %" << termptr->name << ")\n";

	t_astret strretlen_z = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << strretlen_z->name << " = add i64 %" << strretlen->name << ", 1\n";

	// if no memory block is given, allocate one
	if(!mem)
	{
		mem = get_tmp_var(SymbolType::STRING, &sym->dims, nullptr, true);
		(*m_ostr) << "%" << mem->name << " = call i8* @ext_heap_alloc(i64 %"
			<< strretlen_z->name << ", i64 1" << ")\n";
	}

	(*m_ostr) << "call i8* @strncpy(i8* %" << mem->name << ", i8* %" << termptr->name
		<< ", i64 " << dim << ")\n";

	return mem;
}


/**
 * copy a a memory block to a composite symbol
 */
t_astret LLAsm::cp_mem_comp(t_astret mem, t_astret sym)
{
	std::size_t len = get_bytesize(sym);

	(*m_ostr) << "call i8* @memcpy(i8* %" << sym->name << ", i8* %" << mem->name
		<< ", i64 " << len << ")\n";

	return sym;
}


/**
 * copy a a memory block to a vector
 */
t_astret LLAsm::cp_mem_vec(t_astret mem, t_astret sym, bool alloc_sym)
{
	std::size_t argdim = std::get<0>(sym->dims);
	if(sym->ty == SymbolType::MATRIX)
		argdim *= std::get<1>(sym->dims);

	// allocate memory for local array copy
	if(alloc_sym)
		(*m_ostr) << "%" << sym->name << " = alloca [" << argdim << " x double]\n";

	t_astret arrptr = get_tmp_var();
	(*m_ostr) << "%" << arrptr->name << " = getelementptr [" << argdim << " x double], ["
		<< argdim << " x double]* %" << sym->name << ", i64 0, i64 0\n";

	// copy array
	t_astret arrptr_cast = get_tmp_var();

	// cast to memcpy argument pointer type
	(*m_ostr) << "%" << arrptr_cast->name << " = bitcast double* %" << arrptr->name << " to i8*\n";

	(*m_ostr) << "call i8* @memcpy(i8* %" << arrptr_cast->name << ", i8* %" << mem->name
		<< ", i64 " << argdim*sizeof(double) << ")\n";

	return sym;
}


/**
 * copy a a memory block to a string
 */
t_astret LLAsm::cp_mem_str(t_astret mem, t_astret sym, bool alloc_sym)
{
	if(alloc_sym)
		(*m_ostr) << "%" << sym->name << " = alloca [" << std::get<0>(sym->dims) << " x i8]\n";

	t_astret strptr = get_tmp_var();
	(*m_ostr) << "%" << strptr->name << " = getelementptr [" << std::get<0>(sym->dims) << " x i8], ["
		<< std::get<0>(sym->dims) << " x i8]* %" << sym->name << ", i64 0, i64 0\n";

	// copy string
	(*m_ostr) << "call i8* @strncpy(i8* %" << strptr->name << ", i8* %" << mem->name
		<< ", i64 " << std::get<0>(sym->dims) << ")\n";

	return sym;
}


t_astret LLAsm::safe_array_index(t_astret idx, std::size_t size)
{
	// modidx = idx % size
	t_astret modidx = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << modidx->name << " = srem i64 %"
		<< idx->name << ", " << size << "\n";

	// if(modidx < 0) modidx2 = size - modidx
	t_astret modidx2 = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << modidx2->name << " = alloca i64\n";

	generate_cond(
	[this, modidx]() -> t_astret
	{
		t_astret _cond = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << _cond->name << " = icmp slt i64 %"
			<< modidx->name << ", 0\n";
		return _cond;
	}, [this, modidx, modidx2, size]
	{
		t_astret diff = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << diff->name << " = sub i64 " << size
			<< ", %" << modidx->name << "\n";

		(*m_ostr) << "store i64 %" << diff->name
			<< ", i64* %" << modidx2->name << "\n";
	}, [this, modidx, modidx2]
	{
		(*m_ostr) << "store i64 %" << modidx->name
			<< ", i64* %" << modidx2->name << "\n";
	}, true);

	// write result back to an i64 variable
	t_astret modidx3 = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << modidx3->name
		<< " = load i64, i64* %" << modidx2->name << "\n";
	return modidx3;
}
