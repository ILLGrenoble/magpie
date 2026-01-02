/**
 * llvm three-address code generator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 11-apr-20
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

#ifndef __LLASM_H__
#define __LLASM_H__

#include "ast.h"
#include "sym.h"

#include <stack>
#include <optional>
#include <functional>


/**
 * multiply the elements of a container
 */
template<class t_cont, std::size_t ...seq>
constexpr typename t_cont::value_type multiply_elements(
	const t_cont& cont, const std::index_sequence<seq...>&)
{
	return (std::get<seq>(cont) * ...);
}


/**
 * multiply all dimensions of an array type
 */
template<std::size_t NUM_DIMS=2>
std::size_t get_arraydim(const std::array<std::size_t, NUM_DIMS>& dims)
{
	return multiply_elements(dims, std::make_index_sequence<NUM_DIMS>());
}


struct Func
{
	// currently active function
	const ASTFunc* func = nullptr;

	// label after loop
	std::vector<std::string> loopEndLabels{};

	// label skipping loop body, but before incrementing counter
	std::vector<std::string> loopFooterLabels{};
};



class LLAsm : public ASTVisitor
{
public:
	LLAsm(SymTab* syms, std::ostream* ostr=&std::cout);
	virtual ~LLAsm() = default;

	virtual t_astret visit(const ASTUMinus* ast) override;
	virtual t_astret visit(const ASTPlus* ast) override;
	virtual t_astret visit(const ASTMult* ast) override;
	virtual t_astret visit(const ASTMod* ast) override;
	virtual t_astret visit(const ASTPow* ast) override;
	virtual t_astret visit(const ASTTransp* ast) override;
	virtual t_astret visit(const ASTNorm* ast) override;
	virtual t_astret visit(const ASTVar* ast) override;
	virtual t_astret visit(const ASTCall* ast) override;
	virtual t_astret visit(const ASTStmts* ast) override;
	virtual t_astret visit(const ASTVarDecl* ast) override;
	virtual t_astret visit(const ASTFunc* ast) override;
	virtual t_astret visit(const ASTReturn* ast) override;
	virtual t_astret visit(const ASTAssign* ast) override;
	virtual t_astret visit(const ASTArrayAccess* ast) override;
	virtual t_astret visit(const ASTArrayAssign* ast) override;
	virtual t_astret visit(const ASTComp* ast) override;
	virtual t_astret visit(const ASTCond* ast) override;
	virtual t_astret visit(const ASTBool* ast) override;
	virtual t_astret visit(const ASTLoop* ast) override;
	virtual t_astret visit(const ASTLoopJump* ast) override;
	virtual t_astret visit(const ASTStrConst* ast) override;
	virtual t_astret visit(const ASTExprList* ast) override;
	virtual t_astret visit(const ASTNumConst<double>* ast) override;
	virtual t_astret visit(const ASTNumConst<std::int64_t>* ast) override;

	// ------------------------------------------------------------------------
	// internally handled dummy nodes
	// ------------------------------------------------------------------------
	virtual t_astret visit(const ASTArgNames*) override { return nullptr; }
	virtual t_astret visit(const ASTTypeDecl*) override { return nullptr; }
	// ------------------------------------------------------------------------


	/**
	 * output declarations for registered functions
	 */
	static std::string get_function_declarations(const SymTab& symtab, bool only_externals=1);


protected:
	t_astret get_tmp_var(SymbolType ty = SymbolType::SCALAR,
		const std::array<std::size_t, 2>* dims = nullptr,
		const std::string* name = nullptr, bool on_heap = false);

	/**
	 * create a label for branch instructions
	 */
	std::string get_label();

	/**
	 * get the type name for a symbol
	 */
	static std::string get_type_name(SymbolType ty);

	/**
	 * get the element type for an array type
	 */
	static SymbolType get_element_type(SymbolType ty);

	/**
	 * get the (static) byte size of a symbol
	 */
	static std::size_t get_bytesize(t_astret sym);

	/**
	 * get the dimensions of an array type
	 */
	static std::size_t get_arraydim(t_astret sym);

	/**
	 * find the symbol with a specific name in the symbol table
	 */
	t_astret get_sym(const std::string& name) const;

	/**
	 * convert symbol to another type
	 */
	t_astret convert_sym(t_astret sym, SymbolType ty_to);

	/**
	 * check if two symbols can be converted to one another
	 */
	static bool check_sym_compat(
		SymbolType ty1, std::size_t dim1_1, std::size_t dim1_2,
		SymbolType ty2, std::size_t dim2_1, std::size_t dim2_2);

	/**
	 * limits the array index to the [0, size[ range
	 */
	t_astret safe_array_index(t_astret idx, std::size_t size);

	/**
	 * generates if-then-else code
	 */
	template<class t_funcCond, class t_funcBody, class t_funcElseBody>
	void generate_cond(t_funcCond funcCond, t_funcBody funcBody, t_funcElseBody funcElseBody, bool hasElse=0);

	/**
	 * generates loop code
	 */
	template<class t_funcCond, class t_funcBody, class t_funcFooter=std::function<void()>>
	void generate_loop(t_funcCond funcCond, t_funcBody funcBody, std::optional<t_funcFooter> funcFooter=std::nullopt);

	/**
	 * generates loop code with managed counter
	 */
	template<class t_funcBody>
	void generate_loop(std::int64_t start, std::int64_t end, t_funcBody funcBody);


	t_astret cp_comp_mem(t_astret sym, t_astret mem);
	t_astret cp_vec_mem(t_astret sym, t_astret mem=nullptr);
	t_astret cp_str_mem(t_astret sym, t_astret mem=nullptr);

	t_astret cp_mem_comp(t_astret mem, t_astret sym);
	t_astret cp_mem_vec(t_astret mem, t_astret sym, bool alloc_sym=true);
	t_astret cp_mem_str(t_astret mem, t_astret sym, bool alloc_sym=true);


private:
	std::size_t m_varCount = 0; 	// # of tmp vars
	std::size_t m_labelCount = 0;	// # of labels

	std::vector<std::string> m_curscope;

	SymTab* m_syms = nullptr;

	std::ostream* m_ostr = &std::cout;

	// helper functions to reduce code redundancy
	t_astret scalar_matrix_prod(t_astret scalar, t_astret matrix, bool mul_or_div=1);

	// TODO: nested functions
	std::stack<Func> m_funcstack;
};



/**
 * generates if-then-else code
 */
template<class t_funcCond, class t_funcBody, class t_funcElseBody>
void LLAsm::generate_cond(t_funcCond funcCond, t_funcBody funcBody, t_funcElseBody funcElseBody, bool hasElse)
{
	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; condition head\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	t_astret cond = funcCond();
	(*m_ostr) << ";-------------------------------------------------------------\n";

	std::string labelIf = get_label();
	std::string labelElse = hasElse ? get_label() : "";
	std::string labelEnd = get_label();

	if(hasElse)
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelIf << ", label %" << labelElse << "\n";
	else
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelIf << ", label %" << labelEnd << "\n";

	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "; condition body\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << labelIf << ":\n";
	funcBody();
	(*m_ostr) << ";-------------------------------------------------------------\n";

	(*m_ostr) << "br label %" << labelEnd << "\n";

	if(hasElse)
	{
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "; condition \"else\" body\n";
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << labelElse << ":\n";
		funcElseBody();
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "br label %" << labelEnd << "\n";
	}

	(*m_ostr) << labelEnd << ":\n";
	(*m_ostr) << ";-------------------------------------------------------------\n\n";
}



/**
 * generates loop code
 */
template<class t_funcCond, class t_funcBody, class t_funcFooter>
void LLAsm::generate_loop(t_funcCond funcCond, t_funcBody funcBody, std::optional<t_funcFooter> funcFooter)
{
	std::string labelStart = get_label();
	std::string labelBegin = get_label();
	std::string labelEnd = get_label();
	std::string block = get_label();
	std::string labelFooter{};

	Func& actfunc = m_funcstack.top();
	actfunc.loopEndLabels.push_back(labelEnd);
	if(funcFooter)
	{
		labelFooter = get_label();
		actfunc.loopFooterLabels.push_back(labelFooter);
	}

	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; loop head\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "br label %" << labelStart << "\n";
	(*m_ostr) << labelStart << ":\n";
	(*m_ostr) << "%" << block << " = call i8* @llvm.stacksave()\n";
	t_astret cond = funcCond();
	(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelBegin << ", label %" << labelEnd << "\n";

	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << "; loop body\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	(*m_ostr) << labelBegin << ":\n";
	funcBody();
	(*m_ostr) << "call void @llvm.stackrestore(i8* %" << block << ")\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";

	if(funcFooter)
	{
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "; loop foot\n";
		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "br label %" << labelFooter << "\n";	// dummy jump to prevent assembler error
		(*m_ostr) << labelFooter << ":\n";
		(*funcFooter)();
		(*m_ostr) << ";-------------------------------------------------------------\n";
	}

	(*m_ostr) << "br label %" << labelStart << "\n";
	(*m_ostr) << labelEnd << ":\n";
	(*m_ostr) << "call void @llvm.stackrestore(i8* %" << block << ")\n";
	(*m_ostr) << ";-------------------------------------------------------------\n\n";

	if(funcFooter)
		actfunc.loopFooterLabels.pop_back();
	actfunc.loopEndLabels.pop_back();
}



/**
 * generates a loop code with managed counter
 */
template<class t_funcBody>
void LLAsm::generate_loop(std::int64_t start, std::int64_t end, t_funcBody funcBody)
{
	(*m_ostr) << "\n;-------------------------------------------------------------\n";
	(*m_ostr) << "; loop counter\n";
	(*m_ostr) << ";-------------------------------------------------------------\n";
	t_astret ctr = get_tmp_var(SymbolType::INT);
	t_astret ctrval = get_tmp_var(SymbolType::INT);
	(*m_ostr) << "%" << ctr->name << " = alloca i64\n";
	(*m_ostr) << "store i64 " << start << ", i64* %" << ctr->name << "\n";

	generate_loop([this, ctr, ctrval, end]() -> t_astret
	{
		// loop condition: ctr < end
		(*m_ostr) << "%" << ctrval->name << " = load i64, i64* %" << ctr->name << "\n";

		t_astret _cond = get_tmp_var();
		(*m_ostr) << "%" << _cond->name << " = icmp slt i64 %" << ctrval->name <<  ", " << end << "\n";

		return _cond;
	}, [this, &funcBody, ctr, ctrval]
	{
		funcBody(ctrval);

		(*m_ostr) << ";-------------------------------------------------------------\n";
		(*m_ostr) << "; increment loop counter\n";
		(*m_ostr) << ";-------------------------------------------------------------\n";
		t_astret newctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << newctrval->name << " = add i64 %" << ctrval->name << ", 1\n";
		(*m_ostr) << "store i64 %" << newctrval->name << ", i64* %" << ctr->name << "\n";
	});
}

#endif
