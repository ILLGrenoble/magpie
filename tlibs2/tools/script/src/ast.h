/**
 * syntax tree
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

#ifndef __AST_H__
#define __AST_H__

#include <memory>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdint>

#include "sym.h"


class AST;
class ASTUMinus;
class ASTPlus;
class ASTMult;
class ASTMod;
class ASTPow;
class ASTTransp;
class ASTNorm;
class ASTStrConst;
class ASTVar;
class ASTStmts;
class ASTVarDecl;
class ASTArgNames;
class ASTTypeDecl;
class ASTFunc;
class ASTReturn;
class ASTCall;
class ASTAssign;
class ASTArrayAssign;
class ASTArrayAccess;
class ASTComp;
class ASTBool;
class ASTCond;
class ASTLoop;
class ASTLoopJump;
class ASTExprList;
template<class> class ASTNumConst;


enum class ASTType
{
	UMinus,
	Plus,
	Mult,
	Mod,
	Pow,
	Transp,
	Norm,
	StrConst,
	Var,
	Stmts,
	VarDecl,
	ArgNames,
	TypeDecl,
	Func,
	Return,
	Call,
	Assign,
	ArrayAssign,
	ArrayAccess,
	Comp,
	Bool,
	Cond,
	Loop,
	ExprList,
	NumConst,
};


using ASTPtr = std::shared_ptr<AST>;
using t_astret = const Symbol*;


/**
 * ast visitor
 */
class ASTVisitor
{
public:
	virtual ~ASTVisitor() {}

	virtual t_astret visit(const ASTUMinus* ast) = 0;
	virtual t_astret visit(const ASTPlus* ast) = 0;
	virtual t_astret visit(const ASTMult* ast) = 0;
	virtual t_astret visit(const ASTMod* ast) = 0;
	virtual t_astret visit(const ASTPow* ast) = 0;
	virtual t_astret visit(const ASTTransp* ast) = 0;
	virtual t_astret visit(const ASTNorm* ast) = 0;
	virtual t_astret visit(const ASTVar* ast) = 0;
	virtual t_astret visit(const ASTStmts* ast) = 0;
	virtual t_astret visit(const ASTVarDecl* ast) = 0;
	virtual t_astret visit(const ASTArgNames* ast) = 0;
	virtual t_astret visit(const ASTTypeDecl* ast) = 0;
	virtual t_astret visit(const ASTFunc* ast) = 0;
	virtual t_astret visit(const ASTReturn* ast) = 0;
	virtual t_astret visit(const ASTCall* ast) = 0;
	virtual t_astret visit(const ASTAssign* ast) = 0;
	virtual t_astret visit(const ASTArrayAssign* ast) = 0;
	virtual t_astret visit(const ASTArrayAccess* ast) = 0;
	virtual t_astret visit(const ASTComp* ast) = 0;
	virtual t_astret visit(const ASTBool* ast) = 0;
	virtual t_astret visit(const ASTCond* ast) = 0;
	virtual t_astret visit(const ASTLoop* ast) = 0;
	virtual t_astret visit(const ASTLoopJump* ast) = 0;
	virtual t_astret visit(const ASTStrConst* ast) = 0;
	virtual t_astret visit(const ASTExprList* ast) = 0;
	virtual t_astret visit(const ASTNumConst<double>* ast) = 0;
	virtual t_astret visit(const ASTNumConst<std::int64_t>* ast) = 0;
};


#define ASTVISITOR_ACCEPT virtual t_astret accept(ASTVisitor* visitor) const override { return visitor->visit(this); }


/**
 * ast node base
 */
class AST
{
public:
	virtual t_astret accept(ASTVisitor* visitor) const = 0;
	virtual ASTType type() = 0;

	virtual ~AST() {}
};


class ASTUMinus : public AST
{
public:
	ASTUMinus(ASTPtr term)
	: term{term}
	{}

	const ASTPtr GetTerm() const { return term; }

	virtual ASTType type() override { return ASTType::UMinus; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term;
};


class ASTPlus : public AST
{
public:
	ASTPlus(ASTPtr term1, ASTPtr term2,
		bool invert = 0)
		: term1{term1}, term2{term2}, inverted{invert}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }
	bool IsInverted() const { return inverted; }

	virtual ASTType type() override { return ASTType::Plus; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
	bool inverted = 0;
};


class ASTMult : public AST
{
public:
	ASTMult(ASTPtr term1, ASTPtr term2,
		bool invert = 0)
		: term1{term1}, term2{term2}, inverted{invert}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }
	bool IsInverted() const { return inverted; }

	virtual ASTType type() override { return ASTType::Mult; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
	bool inverted = 0;
};


class ASTMod : public AST
{
public:
	ASTMod(ASTPtr term1, ASTPtr term2)
		: term1{term1}, term2{term2}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }

	virtual ASTType type() override { return ASTType::Mod; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
};


class ASTPow : public AST
{
public:
	ASTPow(ASTPtr term1, ASTPtr term2)
		: term1{term1}, term2{term2}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }

	virtual ASTType type() override { return ASTType::Pow; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
};


class ASTTransp : public AST
{
public:
	ASTTransp(ASTPtr term) : term{term}
	{}

	const ASTPtr GetTerm() const { return term; }

	virtual ASTType type() override { return ASTType::Transp; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term;
};


class ASTNorm : public AST
{
public:
	ASTNorm(ASTPtr term) : term{term}
	{}

	const ASTPtr GetTerm() const { return term; }

	virtual ASTType type() override { return ASTType::Norm; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term;
};


class ASTVar : public AST
{
public:
	ASTVar(const std::string& ident)
		: ident{ident}
	{}

	const std::string& GetIdent() const { return ident; }

	virtual ASTType type() override { return ASTType::Var; }
	ASTVISITOR_ACCEPT

private:
	std::string ident;
};


class ASTStmts : public AST
{
public:
	ASTStmts() : stmts{}
	{}

	void AddStatement(ASTPtr stmt, bool front=false)
	{
		if(front)
			stmts.push_front(stmt);
		else
			stmts.push_back(stmt);
	}

	const std::list<ASTPtr>& GetStatementList() const
	{
		return stmts;
	}

	virtual ASTType type() override { return ASTType::Stmts; }
	ASTVISITOR_ACCEPT

private:
	std::list<ASTPtr> stmts;
};


class ASTVarDecl : public AST
{
public:
	ASTVarDecl()
		: vars{}
	{}

	ASTVarDecl(std::shared_ptr<ASTAssign> optAssign)
		: vars{}, optAssign{optAssign}
	{}

	void AddVariable(const std::string& var) { vars.push_front(var); }
	const std::list<std::string>& GetVariables() const { return vars; }

	const std::shared_ptr<ASTAssign> GetAssignment() const { return optAssign; }

	virtual ASTType type() override { return ASTType::VarDecl; }
	ASTVISITOR_ACCEPT

private:
	std::list<std::string> vars;

	// optional assignment
	std::shared_ptr<ASTAssign> optAssign;
};


class ASTArgNames : public AST
{
public:
	ASTArgNames() : argnames{}
	{}

	void AddArg(const std::string& argname,
		SymbolType ty=SymbolType::UNKNOWN,
		std::size_t dim1=1, std::size_t dim2=1)
	{
		argnames.push_front(std::make_tuple(argname, ty, dim1, dim2));
	}

	const std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>>& GetArgs() const
	{
		return argnames;
	}

	std::vector<std::string> GetArgIdents() const
	{
		std::vector<std::string> idents;
		for(const auto& arg : argnames)
			idents.push_back(std::get<0>(arg));
		return idents;
	}

	std::vector<SymbolType> GetArgTypes() const
	{
		std::vector<SymbolType> ty;
		for(const auto& arg : argnames)
			ty.push_back(std::get<1>(arg));
		return ty;
	}

	virtual ASTType type() override { return ASTType::ArgNames; }
	ASTVISITOR_ACCEPT

private:
	std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>> argnames;
};


class ASTTypeDecl : public AST
{
public:
	ASTTypeDecl(SymbolType ty, std::size_t dim1=1, std::size_t dim2=1)
		: ty{ty}, dim1{dim1}, dim2{dim2}
	{}

	SymbolType GetType() const { return ty; }

	std::size_t GetDim(int i=0) const
	{
		if(i==0) return dim1;
		else if(i==1) return dim2;
		return 0;
	}

	std::tuple<SymbolType, std::size_t, std::size_t> GetRet() const
	{
		return std::make_tuple(ty, dim1, dim2);
	}

	virtual ASTType type() override { return ASTType::TypeDecl; }
	ASTVISITOR_ACCEPT

private:
	SymbolType ty;
	std::size_t dim1=1, dim2=1;
};


class ASTFunc : public AST
{
public:
	ASTFunc(const std::string& ident, std::shared_ptr<ASTTypeDecl>& rettype,
		std::shared_ptr<ASTArgNames> args, std::shared_ptr<ASTStmts> stmts,
		std::shared_ptr<ASTArgNames> rets = nullptr)
		: ident{ident}, rettype{rettype->GetRet()}, args{args->GetArgs()},
			stmts{stmts}, rets{}
	{
		if(rets)
			this->rets = rets->GetArgs();
	}

	const std::string& GetIdent() const { return ident; }
	std::tuple<SymbolType, std::size_t, std::size_t> GetRetType() const { return rettype; }

	const std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>>&
	GetArgs() const { return args; }

	const std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>>&
	GetRets() const { return rets; }

	const std::shared_ptr<ASTStmts> GetStatements() const { return stmts; }

	virtual ASTType type() override { return ASTType::Func; }
	ASTVISITOR_ACCEPT

private:
	std::string ident;
	std::tuple<SymbolType, std::size_t, std::size_t> rettype;
	std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>> args;
	std::shared_ptr<ASTStmts> stmts;
	std::list<std::tuple<std::string, SymbolType, std::size_t, std::size_t>> rets;
};


class ASTReturn : public AST
{
public:
	ASTReturn(const std::shared_ptr<ASTExprList>& rets = nullptr,
		bool multi_return=false)
		: rets{rets}, multi_return{multi_return}
	{}

	const std::shared_ptr<ASTExprList> GetRets() const { return rets; }

	bool IsMultiReturn() const { return multi_return; }

	virtual ASTType type() override { return ASTType::Return; }
	ASTVISITOR_ACCEPT

private:
	std::shared_ptr<ASTExprList> rets;
	bool multi_return = false;
};


class ASTExprList : public AST
{
public:
	ASTExprList()
	{}

	void AddExpr(ASTPtr expr)
	{
		exprs.push_front(expr);
	}

	const std::list<ASTPtr>& GetList() const
	{
		return exprs;
	}

	/**
	 * specialised use as an array of scalars
	 */
	void SetScalarArray(bool b)
	{
		m_scalararray = b;
	}

	bool IsScalarArray() const
	{
		return m_scalararray;
	}

	virtual ASTType type() override { return ASTType::ExprList; }
	ASTVISITOR_ACCEPT

private:
	std::list<ASTPtr> exprs;
	bool m_scalararray = false;
};


class ASTCall : public AST
{
public:
	ASTCall(const std::string& ident)
		: ident{ident}, args{std::make_shared<ASTExprList>()}
	{}

	ASTCall(const std::string& ident, std::shared_ptr<ASTExprList> args)
		: ident{ident}, args{args}
	{}

	const std::string& GetIdent() const { return ident; }
	const std::list<ASTPtr>& GetArgumentList() const { return args->GetList(); }

	virtual ASTType type() override { return ASTType::Call; }
	ASTVISITOR_ACCEPT

private:
	std::string ident;
	std::shared_ptr<ASTExprList> args;
};


class ASTAssign : public AST
{
public:
	ASTAssign(const std::string& ident, ASTPtr expr)
		: idents{{ident}}, expr{expr}
	{}

	ASTAssign(const std::vector<std::string>& idents, ASTPtr expr)
		: idents{idents}, expr{expr}
	{}

	const std::vector<std::string>& GetIdents() const { return idents; }
	const std::string& GetIdent() const { return GetIdents()[0]; }
	const ASTPtr GetExpr() const { return expr; }

	bool IsMultiAssign() const { return idents.size() > 1; }

	virtual ASTType type() override { return ASTType::Assign; }
	ASTVISITOR_ACCEPT

private:
	std::vector<std::string> idents;
	ASTPtr expr;
};


class ASTComp : public AST
{
public:
	enum CompOp
	{
		EQU, NEQ,
		GT, LT, GEQ, LEQ
	};

public:
	ASTComp(ASTPtr term1, ASTPtr term2, CompOp op)
		: term1{term1}, term2{term2}, op{op}
	{}

	ASTComp(ASTPtr term1, CompOp op)
		: term1{term1}, term2{nullptr}, op{op}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }
	CompOp GetOp() const { return op; }

	virtual ASTType type() override { return ASTType::Comp; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
	CompOp op;
};


class ASTBool : public AST
{
public:
	enum BoolOp
	{
		NOT,
		AND, OR, XOR
	};

public:
	ASTBool(ASTPtr term1, ASTPtr term2, BoolOp op)
		: term1{term1}, term2{term2}, op{op}
	{}

	ASTBool(ASTPtr term1, BoolOp op)
		: term1{term1}, term2{nullptr}, op{op}
	{}

	const ASTPtr GetTerm1() const { return term1; }
	const ASTPtr GetTerm2() const { return term2; }
	BoolOp GetOp() const { return op; }

	virtual ASTType type() override { return ASTType::Bool; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term1, term2;
	BoolOp op;
};


class ASTCond : public AST
{
public:
	ASTCond(const ASTPtr cond, ASTPtr if_stmt)
		: cond{cond}, if_stmt{if_stmt}
	{}
	ASTCond(const ASTPtr cond, ASTPtr if_stmt, ASTPtr else_stmt)
		: cond{cond}, if_stmt{if_stmt}, else_stmt{else_stmt}
	{}

	const ASTPtr GetCond() const { return cond; }
	const ASTPtr GetIf() const { return if_stmt; }
	const ASTPtr GetElse() const { return else_stmt; }
	bool HasElse() const { return else_stmt != nullptr; }

	virtual ASTType type() override { return ASTType::Cond; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr cond;
	ASTPtr if_stmt, else_stmt;
};


class ASTLoop : public AST
{
public:
	ASTLoop(const ASTPtr cond, ASTPtr stmt, ASTPtr footer=nullptr)
		: cond{cond}, stmt{stmt}, footer{footer}
	{}

	const ASTPtr GetCond() const { return cond; }
	const ASTPtr GetLoopStmt() const { return stmt; }
	const ASTPtr GetLoopFooter() const { return footer; }

	virtual ASTType type() override { return ASTType::Loop; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr cond, stmt, footer;
};


class ASTLoopJump : public AST
{
public:
	enum LoopJumpType
	{
		SKIP,
		END,
		ENDALL
	};

public:
	ASTLoopJump(LoopJumpType kind, std::size_t levels=1)
		: m_kind{kind}, m_levels{levels}
	{}

	LoopJumpType GetKind() const { return m_kind; }
	std::size_t GetLevels() const { return m_levels; }

	virtual ASTType type() override { return ASTType::Loop; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr cond, stmt;
	LoopJumpType m_kind{LoopJumpType::SKIP};
	std::size_t m_levels{1};
};


class ASTArrayAccess : public AST
{
public:
	ASTArrayAccess(ASTPtr term,
		ASTPtr num1, ASTPtr num2 = nullptr,
		ASTPtr num3 = nullptr, ASTPtr num4 = nullptr,
		bool ranged12 = false, bool ranged34 = false)
		: term{term}, num1{num1}, num2{num2}, num3{num3}, num4{num4},
			ranged12{ranged12}, ranged34{ranged34}
	{}

	const ASTPtr GetTerm() const { return term; }

	const ASTPtr GetNum1() const { return num1; }
	const ASTPtr GetNum2() const { return num2; }
	const ASTPtr GetNum3() const { return num3; }
	const ASTPtr GetNum4() const { return num4; }

	bool IsRanged12() const { return ranged12; }
	bool IsRanged34() const { return ranged34; }

	virtual ASTType type() override { return ASTType::ArrayAccess; }
	ASTVISITOR_ACCEPT

private:
	ASTPtr term;

	ASTPtr num1, num2;
	ASTPtr num3, num4;
	bool ranged12 = false;
	bool ranged34 = false;
};


class ASTArrayAssign : public AST
{
public:
	ASTArrayAssign(const std::string& ident, ASTPtr expr,
		ASTPtr num1, ASTPtr num2 = nullptr,
		ASTPtr num3 = nullptr, ASTPtr num4 = nullptr,
		bool ranged12 = false, bool ranged34 = false)
		: ident{ident}, expr{expr}, num1{num1}, num2{num2},
			num3{num3}, num4{num4}, ranged12{ranged12}, ranged34{ranged34}
	{}

	const std::string& GetIdent() const { return ident; }
	const ASTPtr GetExpr() const { return expr; }

	const ASTPtr GetNum1() const { return num1; }
	const ASTPtr GetNum2() const { return num2; }
	const ASTPtr GetNum3() const { return num3; }
	const ASTPtr GetNum4() const { return num4; }

	bool IsRanged12() const { return ranged12; }
	bool IsRanged34() const { return ranged34; }

	virtual ASTType type() override { return ASTType::ArrayAssign; }
	ASTVISITOR_ACCEPT

private:
	std::string ident;
	ASTPtr expr;

	ASTPtr num1, num2;
	ASTPtr num3, num4;
	bool ranged12 = false;
	bool ranged34 = false;
};


template<class t_num>
class ASTNumConst : public AST
{
public:
	ASTNumConst(t_num val) : val{val}
	{}

	t_num GetVal() const { return val; }

	virtual ASTType type() override { return ASTType::NumConst; }
	ASTVISITOR_ACCEPT

private:
	t_num val{};
};


class ASTStrConst : public AST
{
public:
	ASTStrConst(const std::string& str) : val{str}
	{}

	const std::string& GetVal() const { return val; }

	virtual ASTType type() override { return ASTType::StrConst; }
	ASTVISITOR_ACCEPT

private:
	std::string val;
};


#endif
