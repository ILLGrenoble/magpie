/**
 * outputs the syntax tree
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
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

#include "printast.h"
#include "sym.h"


ASTPrinter::ASTPrinter(std::ostream* ostr)
	: m_ostr{ostr}
{
}


t_astret ASTPrinter::visit(const ASTUMinus* ast)
{
	(*m_ostr) << "<UMinus>\n";
	ast->GetTerm()->accept(this);
	(*m_ostr) << "</UMinus>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTPlus* ast)
{
	(*m_ostr) << "<Plus>\n";
	ast->GetTerm1()->accept(this);
	ast->GetTerm2()->accept(this);
	(*m_ostr) << "</Plus>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTMult* ast)
{
	(*m_ostr) << "<Mult>\n";
	ast->GetTerm1()->accept(this);
	ast->GetTerm2()->accept(this);
	(*m_ostr) << "</Mult>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTMod* ast)
{
	(*m_ostr) << "<Mod>\n";
	ast->GetTerm1()->accept(this);
	ast->GetTerm2()->accept(this);
	(*m_ostr) << "</Mod>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTPow* ast)
{
	(*m_ostr) << "<Pow>\n";
	ast->GetTerm1()->accept(this);
	ast->GetTerm2()->accept(this);
	(*m_ostr) << "</Pow>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTTransp* ast)
{
	(*m_ostr) << "<Transp>\n";
	ast->GetTerm()->accept(this);
	(*m_ostr) << "</Transp>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTNorm* ast)
{
	(*m_ostr) << "<Norm>\n";
	ast->GetTerm()->accept(this);
	(*m_ostr) << "</Norm>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTVar* ast)
{
	(*m_ostr) << "<Var ident=\"" << ast->GetIdent() << "\" />\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTCall* ast)
{
	(*m_ostr) << "<Call ident=\"" << ast->GetIdent() << "\">\n";

	std::size_t argidx = 0;
	for(const auto& arg : ast->GetArgumentList())
	{
		(*m_ostr) << "<arg_" << argidx << ">\n";
		arg->accept(this);
		(*m_ostr) << "</arg_" << argidx << ">\n";
		++argidx;
	}

	(*m_ostr) << "</Call>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTStmts* ast)
{
	(*m_ostr) << "<Stmts>\n";

	std::size_t stmtidx = 0;
	for(const auto& stmt : ast->GetStatementList())
	{
		(*m_ostr) << "<stmt_" << stmtidx << ">\n";
		stmt->accept(this);
		(*m_ostr) << "</stmt_" << stmtidx << ">\n";
		++stmtidx;
	}

	(*m_ostr) << "</Stmts>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTVarDecl* ast)
{
	(*m_ostr) << "<VarDecl>\n";

	std::size_t varidx = 0;
	for(const auto& var : ast->GetVariables())
	{
		(*m_ostr) << "<var_" << varidx << " ident=\"" << var << "\" />\n";
		++varidx;
	}

	if(ast->GetAssignment())
		ast->GetAssignment()->accept(this);

	(*m_ostr) << "</VarDecl>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTFunc* ast)
{
	(*m_ostr) << "<Func ident=\"" << ast->GetIdent() << "\"" << ">\n";

	auto ret = ast->GetRetType();
	std::string retTypeName = Symbol::get_type_name(std::get<0>(ret));
	(*m_ostr) << "<ret type=\"" << retTypeName << "\" dim1=\"" << std::get<1>(ret) << "\" dim2=\"" << std::get<2>(ret) << "\"";
	(*m_ostr) << " />\n";

	std::size_t argidx = 0;
	for(const auto& arg : ast->GetArgs())
	{
		std::string argTypeName = Symbol::get_type_name(std::get<1>(arg));

		(*m_ostr) << "<arg_" << argidx << " name=\"" << std::get<0>(arg) << "\"";
		(*m_ostr) << " type=\"" << argTypeName << "\" dim1=\"" << std::get<2>(arg) << "\" dim2=\"" << std::get<3>(arg) << "\"";
		(*m_ostr) << " />\n";
		++argidx;
	}

	ast->GetStatements()->accept(this);
	(*m_ostr) << "</Func>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTReturn* ast)
{
	const auto& rets = ast->GetRets();
	if(!rets)
	{
		(*m_ostr) << "<Return />\n";
		return nullptr;
	}

	const auto& retvals = rets->GetList();
	std::size_t numRets = retvals.size();

	if(numRets == 0)
	{
		(*m_ostr) << "<Return />\n";
	}
	else if(numRets == 1)
	{
		(*m_ostr) << "<Return>\n";
		(*retvals.begin())->accept(this);
		(*m_ostr) << "</Return>\n";
	}
	else if(numRets > 1)
	{
		(*m_ostr) << "<MultiReturn>\n";
		std::size_t elemnr = 0;
		for(const auto& elem : retvals)
		{
			(*m_ostr) << "<val_" << elemnr << ">\n";
			elem->accept(this);
			(*m_ostr) << "</val_" << elemnr << ">\n";

			++elemnr;
		}
		(*m_ostr) << "</MultiReturn>\n";
	}

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTAssign* ast)
{
	// multiple assignments
	if(ast->IsMultiAssign())
	{
		(*m_ostr) << "<MultiAssign>\n";

		const auto& idents = ast->GetIdents();
		std::size_t identidx = 0;
		for(const auto& ident : idents)
		{
			(*m_ostr) << "<ident_" << identidx << ">";
			(*m_ostr) << ident;
			(*m_ostr) << "</ident_" << identidx << ">\n";
			++identidx;
		}

		(*m_ostr) << "<rhs>\n";
		ast->GetExpr()->accept(this);
		(*m_ostr) << "</rhs>\n";

		(*m_ostr) << "</MultiAssign>\n";
	}

	// single assignment
	else
	{
		(*m_ostr) << "<Assign ident=\"" << ast->GetIdent() << "\">\n";
		ast->GetExpr()->accept(this);
		(*m_ostr) << "</Assign>\n";
	}

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTArrayAccess* ast)
{
	(*m_ostr) << "<ArrayAccess"
		<< " is_range_12=\"" << ast->IsRanged12() << "\""
		<< " is_range_34=\"" << ast->IsRanged34() << "\""
		<< ">\n";

	(*m_ostr) << "<idx1>\n";
	ast->GetNum1()->accept(this);
	(*m_ostr) << "</idx1>\n";

	if(ast->GetNum2())
	{
		(*m_ostr) << "<idx2>\n";
		ast->GetNum2()->accept(this);
		(*m_ostr) << "</idx2>\n";
	}

	if(ast->GetNum3())
	{
		(*m_ostr) << "<idx3>\n";
		ast->GetNum3()->accept(this);
		(*m_ostr) << "</idx3>\n";
	}

	if(ast->GetNum4())
	{
		(*m_ostr) << "<idx4>\n";
		ast->GetNum4()->accept(this);
		(*m_ostr) << "</idx4>\n";
	}

	(*m_ostr) << "<term>\n";
	ast->GetTerm()->accept(this);
	(*m_ostr) << "</term>\n";

	(*m_ostr) << "</ArrayAccess>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTArrayAssign* ast)
{
	(*m_ostr) << "<ArrayAssign ident=\"" << ast->GetIdent() << "\""
		<< " is_range_12=\"" << ast->IsRanged12() << "\""
		<< " is_range_34=\"" << ast->IsRanged34() << "\""
		<< ">\n";

	(*m_ostr) << "<idx1>\n";
	ast->GetNum1()->accept(this);
	(*m_ostr) << "</idx1>\n";

	if(ast->GetNum2())
	{
		(*m_ostr) << "<idx2>\n";
		ast->GetNum2()->accept(this);
		(*m_ostr) << "</idx2>\n";
	}

	if(ast->GetNum3())
	{
		(*m_ostr) << "<idx3>\n";
		ast->GetNum3()->accept(this);
		(*m_ostr) << "</idx3>\n";
	}

	if(ast->GetNum4())
	{
		(*m_ostr) << "<idx4>\n";
		ast->GetNum2()->accept(this);
		(*m_ostr) << "</idx4>\n";
	}

	(*m_ostr) << "<expr>\n";
	ast->GetExpr()->accept(this);
	(*m_ostr) << "</expr>\n";

	(*m_ostr) << "</ArrayAssign>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTComp* ast)
{
	(*m_ostr) << "<Comp op=\"";
	switch(ast->GetOp())
	{
		case ASTComp::EQU: (*m_ostr) << "equ"; break;
		case ASTComp::NEQ: (*m_ostr) << "neq"; break;
		case ASTComp::GT: (*m_ostr) << "gt"; break;
		case ASTComp::LT: (*m_ostr) << "lt"; break;
		case ASTComp::GEQ: (*m_ostr) << "geq"; break;
		case ASTComp::LEQ: (*m_ostr) << "leq"; break;
		default: (*m_ostr) << "unknown"; break;
	}
	(*m_ostr) << "\">\n";

	ast->GetTerm1()->accept(this);
	ast->GetTerm2()->accept(this);
	(*m_ostr) << "</Comp>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTCond* ast)
{
	(*m_ostr) << "<Cond>\n";

	(*m_ostr) << "<cond>\n";
	ast->GetCond()->accept(this);
	(*m_ostr) << "</cond>\n";

	(*m_ostr) << "<if>\n";
	ast->GetIf()->accept(this);
	(*m_ostr) << "</if>\n";

	if(ast->GetElse())
	{
		(*m_ostr) << "<else>\n";
		ast->GetElse()->accept(this);
		(*m_ostr) << "</else>\n";
	}

	(*m_ostr) << "</Cond>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTBool* ast)
{
	(*m_ostr) << "<Bool op=\"";
	switch(ast->GetOp())
	{
		case ASTBool::NOT: (*m_ostr) << "not"; break;
		case ASTBool::AND: (*m_ostr) << "and"; break;
		case ASTBool::OR: (*m_ostr) << "or"; break;
		case ASTBool::XOR: (*m_ostr) << "xor"; break;
		default: (*m_ostr) << "unknown"; break;
	}
	(*m_ostr) << "\">\n";

	ast->GetTerm1()->accept(this);
	if(ast->GetTerm2())
		ast->GetTerm2()->accept(this);

	(*m_ostr) << "</Bool>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTLoop* ast)
{
	(*m_ostr) << "<Loop>\n";

	(*m_ostr) << "<cond>\n";
	ast->GetCond()->accept(this);
	(*m_ostr) << "</cond>\n";

	(*m_ostr) << "<stmt>\n";
	ast->GetLoopStmt()->accept(this);
	(*m_ostr) << "</stmt>\n";

	(*m_ostr) << "</Loop>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTLoopJump* ast)
{
	(*m_ostr) << "<LoopJump kind=\"";

	switch(ast->GetKind())
	{
		case ASTLoopJump::SKIP: (*m_ostr) << "skip"; break;
		case ASTLoopJump::END: (*m_ostr) << "end"; break;
		case ASTLoopJump::ENDALL: (*m_ostr) << "endall"; break;
		default: (*m_ostr) << "unknown"; break;
	}

	(*m_ostr) << "\", levels=\"" << ast->GetLevels() << "\"/>\n";
	return nullptr;
}


t_astret ASTPrinter::visit(const ASTStrConst* ast)
{
	(*m_ostr) << "<Const type=\"str\" val=\"" << ast->GetVal() << "\" />\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTExprList* ast)
{
	(*m_ostr) << "<ExprList>\n";

	std::size_t expridx = 0;
	for(const auto& expr : ast->GetList())
	{
		(*m_ostr) << "<expr_" << expridx << ">\n";
		expr->accept(this);
		(*m_ostr) << "</expr_" << expridx << ">\n";
		++expridx;
	}

	(*m_ostr) << "</ExprList>\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTNumConst<double>* ast)
{
	(*m_ostr) << "<Const type=\"scalar\" val=\"" << ast->GetVal() << "\" />\n";

	return nullptr;
}


t_astret ASTPrinter::visit(const ASTNumConst<std::int64_t>* ast)
{
	(*m_ostr) << "<Const type=\"int\" val=\"" << ast->GetVal() << "\" />\n";

	return nullptr;
}
