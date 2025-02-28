/**
 * llvm three-address code generator -- operators
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


t_astret LLAsm::visit(const ASTUMinus* ast)
{
	t_astret term = ast->GetTerm()->accept(this);
	t_astret var = get_tmp_var(term->ty, &term->dims);

	// array types
	if(term->ty == SymbolType::VECTOR || term->ty == SymbolType::MATRIX)
	{
		std::size_t dim = std::get<0>(term->dims);
		if(term->ty == SymbolType::MATRIX)
			dim *= std::get<1>(term->dims);

		// allocate double array for result
		t_astret vec_mem = get_tmp_var(term->ty, &term->dims);
		(*m_ostr) << "%" << vec_mem->name << " = alloca [" << dim << " x double]\n";

		// copy array elements in a loop
		generate_loop(0, dim, [this, term, vec_mem, dim](t_astret ctrval)
		{
			t_astret elemptr_src = get_tmp_var();
			t_astret elem_src = get_tmp_var();

			(*m_ostr) << "%" << elemptr_src->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term->name << ", i64 0, i64 %" << ctrval->name << "\n";
			(*m_ostr) << "%" << elem_src->name << " = load double, double* %" << elemptr_src->name << "\n";

			t_astret elemptr_dst = get_tmp_var();
			(*m_ostr) << "%" << elemptr_dst->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << vec_mem->name << ", i64 0, i64 %" << ctrval->name << "\n";

			// unary minus
			t_astret elem_dst = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << elem_dst->name << " = fsub double 0." << ", %" << elem_src->name << "\n";

			// save result in array
			(*m_ostr) << "store double %" << elem_dst->name << ", double* %" << elemptr_dst->name << "\n";
		});

		return vec_mem;
	}

	else if(term->ty == SymbolType::SCALAR)
	{
		(*m_ostr) << "%" << var->name << " = fneg " << LLAsm::get_type_name(term->ty)
			<< " %" << term->name << "\n";
		return var;
	}

	else if(term->ty == SymbolType::INT)
	{
		(*m_ostr) << "%" << var->name << " = sub " << LLAsm::get_type_name(term->ty) << " "
			<< "0, %" << term->name << "\n";
		return var;
	}

	throw std::runtime_error("ASTUMinus: Invalid unary subtraction operation of \"" + term->name + "\".");
	return nullptr;
}


t_astret LLAsm::visit(const ASTPlus* ast)
{
	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = ast->GetTerm2()->accept(this);

	// array types
	if(term1->ty == SymbolType::VECTOR || term1->ty == SymbolType::MATRIX)
	{
		if(term2->ty != term1->ty)
		{
			throw std::runtime_error("ASTPlus: Type mismatch in addition/subtraction of \""
				+ term1->name + "\" and \"" + term2->name + "\".");
		}

		if(std::get<0>(term1->dims) != std::get<0>(term2->dims))
		{
			throw std::runtime_error("ASTPlus: Dimension mismatch in addition/subtraction of \""
				+ term1->name + "\" and \"" + term2->name + "\".");
		}

		std::size_t dim = std::get<0>(term1->dims);
		if(term1->ty == SymbolType::MATRIX)
		{
			throw std::runtime_error("ASTPlus: Dimension mismatch in addition/subtraction of \""
				+ term1->name + "\" and \"" + term2->name + "\".");

			dim *= std::get<1>(term1->dims);
		}

		// allocate double array for result
		t_astret vec_mem = get_tmp_var(term1->ty, &term1->dims);
		(*m_ostr) << "%" << vec_mem->name << " = alloca [" << dim << " x double]\n";

		std::string op = ast->IsInverted() ? "fsub" : "fadd";

		// copy elements in a loop
		generate_loop(0, dim, [this, term1, term2, dim, vec_mem, &op](t_astret ctrval)
		{
			// loop statements
			t_astret elemptr_src1 = get_tmp_var();
			(*m_ostr) << "%" << elemptr_src1->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term1->name << ", i64 0, i64 %" << ctrval->name << "\n";
			t_astret elemptr_src2 = get_tmp_var();
			(*m_ostr) << "%" << elemptr_src2->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term2->name << ", i64 0, i64 %" << ctrval->name << "\n";

			t_astret elem_src1 = get_tmp_var();
			(*m_ostr) << "%" << elem_src1->name << " = load double, double* %" << elemptr_src1->name << "\n";
			t_astret elem_src2 = get_tmp_var();
			(*m_ostr) << "%" << elem_src2->name << " = load double, double* %" << elemptr_src2->name << "\n";

			t_astret elemptr_dst = get_tmp_var();
			(*m_ostr) << "%" << elemptr_dst->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << vec_mem->name << ", i64 0, i64 %" << ctrval->name << "\n";

			// add/subtract
			t_astret elem_dst = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << elem_dst->name << " = " << op << " "
				<< "double %" << elem_src1->name << ", %" << elem_src2->name << "\n";

			// save result in array
			(*m_ostr) << "store double %" << elem_dst->name << ", double* %" << elemptr_dst->name << "\n";
		});

		return vec_mem;
	}

	// concatenate strings
	else if(term1->ty == SymbolType::STRING || term2->ty == SymbolType::STRING)
	{
		// cast to string if needed
		term1 = convert_sym(term1, SymbolType::STRING);
		term2 = convert_sym(term2, SymbolType::STRING);

		// TODO: get individual and total string lengths

		// get string pointers
		t_astret strptr1 = get_tmp_var();
		t_astret strptr2 = get_tmp_var();
		(*m_ostr) << "%" << strptr1->name << " = getelementptr [" << std::get<0>(term1->dims) << " x i8], ["
			<< std::get<0>(term1->dims) << " x i8]* %" << term1->name << ", i64 0, i64 0\n";
		(*m_ostr) << "%" << strptr2->name << " = getelementptr [" << std::get<0>(term2->dims) << " x i8], ["
			<< std::get<0>(term2->dims) << " x i8]* %" << term2->name << ", i64 0, i64 0\n";

		// allocate memory for concatenated string
		std::array<std::size_t, 2> dim{{ std::get<0>(term1->dims)+std::get<0>(term2->dims)-1, 1 }};
		t_astret res = get_tmp_var(SymbolType::STRING, &dim);

		// allocate the concatenated string's memory
		// TODO: use actual new string size, not the (maximum) allocated sizes of the terms in "dim"
		(*m_ostr) << "%" << res->name << " = alloca [" << std::get<0>(dim) << " x i8]\n";

		// get a pointer to the concatenated string
		t_astret resptr = get_tmp_var();
		(*m_ostr) << "%" << resptr->name << " = getelementptr [" << std::get<0>(dim) << " x i8], ["
			<< std::get<0>(dim) << " x i8]* %" << res->name << ", i64 0, i64 0\n";

		// copy first string
		(*m_ostr) << "call i8* @strncpy(i8* %" << resptr->name << ", i8* %" << strptr1->name
			<< ", i64 " << std::get<0>(term1->dims) << ")\n";

		// concatenate second string
		(*m_ostr) << "call i8* @strncat(i8* %" << resptr->name << ", i8* %" << strptr2->name
			<< ", i64 " << std::get<0>(term2->dims) << ")\n";

		return res;
	}

	// scalar types
	else
	{
		// cast if needed
		SymbolType ty = term1->ty;
		if(term1->ty==SymbolType::SCALAR || term2->ty==SymbolType::SCALAR)
			ty = SymbolType::SCALAR;
		t_astret var = get_tmp_var(ty, &term1->dims);

		term1 = convert_sym(term1, ty);
		term2 = convert_sym(term2, ty);

		std::string op = ast->IsInverted() ? "sub" : "add";
		if(ty == SymbolType::SCALAR)
			op = "f" + op;

		(*m_ostr) << "%" << var->name << " = " << op << " "
			<< LLAsm::get_type_name(ty) << " %" << term1->name << ", %" << term2->name << "\n";

		return var;
	}

	throw std::runtime_error("ASTPlus: Invalid addition operation between \""
		+ term1->name + "\" and \"" + term2->name + "\".");
	return nullptr;
}


t_astret LLAsm::visit(const ASTMult* ast)
{
	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = ast->GetTerm2()->accept(this);

	// inner product of vectors: s = v^i v_i
	if(term1->ty == SymbolType::VECTOR && term2->ty == SymbolType::VECTOR)
	{
		if(ast->IsInverted())
		{
			throw std::runtime_error("ASTMult: Cannot divide vector \"" + term1->name
				+ "\" by vector \"" + term2->name + "\".\n");
		}

		if(std::get<0>(term1->dims) != std::get<0>(term2->dims))
		{
			throw std::runtime_error("ASTMult: Dimension mismatch in inner product of \""
				+ term1->name + "\" and \"" + term2->name + "\".");
		}

		std::size_t dim = std::get<0>(term1->dims);


		// result variable
		t_astret dotptr = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << dotptr->name << " = alloca double\n";
		(*m_ostr) << "store double 0., double* %" << dotptr->name << "\n";


		// add elements in a loop
		generate_loop(0, dim, [this, term1, term2, dotptr, dim](t_astret ctrval)
		{
			// loop statements
			t_astret elemptr_src1 = get_tmp_var();
			(*m_ostr) << "%" << elemptr_src1->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term1->name << ", i64 0, i64 %" << ctrval->name << "\n";
			t_astret elemptr_src2 = get_tmp_var();
			(*m_ostr) << "%" << elemptr_src2->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term2->name << ", i64 0, i64 %" << ctrval->name << "\n";

			t_astret elem_src1 = get_tmp_var();
			(*m_ostr) << "%" << elem_src1->name << " = load double, double* %" << elemptr_src1->name << "\n";
			t_astret elem_src2 = get_tmp_var();
			(*m_ostr) << "%" << elem_src2->name << " = load double, double* %" << elemptr_src2->name << "\n";

			// multiply elements
			t_astret mul_elems = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << mul_elems->name << " = fmul "
				<< "double %" << elem_src1->name << ", %" << elem_src2->name << "\n";

			// increment dot product variable
			t_astret cur_dot = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << cur_dot->name << " = load double, double* %" << dotptr->name << "\n";
			t_astret sum_dot = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << sum_dot->name << " = fadd "
				<< "double %" << cur_dot->name << ", %" << mul_elems->name << "\n";
			(*m_ostr) << "store double %" << sum_dot->name << ", double* %" << dotptr->name << "\n";
		});


		t_astret dot = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << dot->name << " = load double, double* %" << dotptr->name << "\n";
		return dot;
	}

	// matrix-vector product: w^i = M^i_j v^j
	else if(term1->ty == SymbolType::MATRIX && term2->ty == SymbolType::VECTOR)
	{
		if(ast->IsInverted())
		{
			throw std::runtime_error("ASTMult: Cannot divide matrix \"" + term1->name
				+ "\" by vector \"" + term2->name + "\".\n");
		}

		if(std::get<1>(term1->dims) != std::get<0>(term2->dims))
		{
			throw std::runtime_error("ASTMult: Dimension mismatch in matrix-vector product of \""
				+ term1->name + "\" and \"" + term2->name + "\".");
		}

		std::size_t dim_i = std::get<0>(term1->dims);
		std::size_t dim_j = std::get<1>(term1->dims);

		// result vector w
		std::array<std::size_t, 2> w_dims{{dim_i, 1}};
		t_astret w_mem = get_tmp_var(SymbolType::VECTOR, &w_dims);
		(*m_ostr) << "%" << w_mem->name << " = alloca [" << dim_i << " x double]\n";

		// loop i
		generate_loop(0, dim_i, [this, term1, term2, dim_i, dim_j, w_mem](t_astret ctr_i_val)
		{
			// w[i] pointer
			t_astret elemptr_w_i = get_tmp_var();
			(*m_ostr) << "%" << elemptr_w_i->name << " = getelementptr [" << dim_i << " x double], ["
				<< dim_i << " x double]* %" << w_mem->name << ", i64 0, i64 %" << ctr_i_val->name << "\n";

			// w[i] = 0
			(*m_ostr) << "store double 0., double* %" << elemptr_w_i->name << "\n";

			// loop i index of M: i*dim_j
			t_astret M_idx_i = get_tmp_var();
			(*m_ostr) << "%" << M_idx_i->name << " = mul i64 %" << ctr_i_val->name << ", " << dim_j << "\n";

			// loop j
			generate_loop(0, dim_j, [this, term1, term2, dim_i, dim_j, M_idx_i, elemptr_w_i](t_astret ctr_j_val)
			{
				// v[j] pointer
				t_astret elemptr_v_j = get_tmp_var();
				(*m_ostr) << "%" << elemptr_v_j->name << " = getelementptr [" << dim_j << " x double], ["
					<< dim_j << " x double]* %" << term2->name << ", i64 0, i64 %" << ctr_j_val->name << "\n";

				// M index: i*dim_j + j
				t_astret M_idx = get_tmp_var();
				(*m_ostr) << "%" << M_idx->name << " = add i64 %" << M_idx_i->name << ", %" << ctr_j_val->name << "\n";

				// M[i,j] pointer
				t_astret elemptr_M_ij = get_tmp_var();
				(*m_ostr) << "%" << elemptr_M_ij->name << " = getelementptr [" << dim_i*dim_j << " x double], ["
					<< dim_i*dim_j << " x double]* %" << term1->name << ", i64 0, i64 %"
					<< M_idx->name << "\n";

				// v[j] value
				t_astret elem_v_j = get_tmp_var();
				(*m_ostr) << "%" << elem_v_j->name << " = load double, double* %" << elemptr_v_j->name << "\n";

				// M[i,j] value
				t_astret elem_M_ij = get_tmp_var();
				(*m_ostr) << "%" << elem_M_ij->name << " = load double, double* %" << elemptr_M_ij->name << "\n";

				// m = M[i,j]*v[j]
				t_astret m = get_tmp_var();
				(*m_ostr) << "%" << m->name << " = fmul double %" << elem_M_ij->name << ", %" << elem_v_j->name << "\n";

				// current w[i] value
				t_astret elem_w_i = get_tmp_var();
				(*m_ostr) << "%" << elem_w_i->name << " = load double, double* %" << elemptr_w_i->name << "\n";

				// d = w[i] + m
				t_astret d = get_tmp_var();
				(*m_ostr) << "%" << d->name << " = fadd double %" << elem_w_i->name << ", %" << m->name << "\n";

				// w[i] = d
				(*m_ostr) << "store double %" << d->name << ", double* %" << elemptr_w_i->name << "\n";
			});
		});

		return w_mem;
	}

	// matrix-matrix product: L^i_j = M^i_k N^k_j
	else if(term1->ty == SymbolType::MATRIX && term2->ty == SymbolType::MATRIX)
	{
		if(ast->IsInverted())
		{
			throw std::runtime_error("ASTMult: Cannot divide matrix \"" + term1->name
				+ "\" by matrix \"" + term2->name + "\".\n");
		}

		if(std::get<1>(term1->dims) != std::get<0>(term2->dims))
		{
			throw std::runtime_error("ASTMult: Dimension mismatch in matrix-matrix product of \""
				+ term1->name + "\" and \"" + term2->name + "\".");
		}

		std::size_t dim_i = std::get<0>(term1->dims);
		std::size_t dim_k = std::get<1>(term1->dims);
		std::size_t dim_j = std::get<1>(term2->dims);

		// result Matrix L
		std::array<std::size_t, 2> L_dims{{dim_i, dim_j}};
		t_astret L_mem = get_tmp_var(SymbolType::MATRIX, &L_dims);
		(*m_ostr) << "%" << L_mem->name << " = alloca [" << dim_i*dim_j << " x double]\n";

		// loop i
		generate_loop(0, dim_i, [this, term1, term2, dim_i, dim_j, dim_k, L_mem](t_astret ctr_i_val)
		{
			// loop i index of L: i*dim_j
			t_astret L_idx_i = get_tmp_var();
			(*m_ostr) << "%" << L_idx_i->name << " = mul i64 %" << ctr_i_val->name << ", " << dim_j << "\n";

			// loop i index of M: i*dim_k
			t_astret M_idx_i = get_tmp_var();
			(*m_ostr) << "%" << M_idx_i->name << " = mul i64 %" << ctr_i_val->name << ", " << dim_k << "\n";

			// loop j
			generate_loop(0, dim_j, [this, term1, term2, dim_i, dim_j, dim_k, L_idx_i, M_idx_i, L_mem](t_astret ctr_j_val)
			{
				// d = 0
				t_astret d = get_tmp_var(SymbolType::SCALAR);
				(*m_ostr) << "%" << d->name << " = alloca double\n";
				(*m_ostr) << "store double 0., double* %" << d->name << "\n";

				// loop k
				generate_loop(0, dim_k, [this, term1, term2, dim_i, dim_j, dim_k, M_idx_i, ctr_j_val, d](t_astret ctr_k_val)
				{
					// M index: i*dim_k + k
					t_astret M_idx = get_tmp_var();
					(*m_ostr) << "%" << M_idx->name << " = add i64 %" << M_idx_i->name << ", %" << ctr_k_val->name << "\n";

					// M[i,k] pointer
					t_astret elemptr_M_ik = get_tmp_var();
					(*m_ostr) << "%" << elemptr_M_ik->name << " = getelementptr [" << dim_i*dim_k << " x double], ["
						<< dim_i*dim_k << " x double]* %" << term1->name << ", i64 0, i64 %"
						<< M_idx->name << "\n";

					// M[i,k] value
					t_astret elem_M_ik = get_tmp_var();
					(*m_ostr) << "%" << elem_M_ik->name << " = load double, double* %" << elemptr_M_ik->name << "\n";


					// loop k index of N: k*dim_j
					t_astret N_idx_k = get_tmp_var();
					(*m_ostr) << "%" << N_idx_k->name << " = mul i64 %" << ctr_k_val->name << ", " << dim_j << "\n";

					// N index: k*dim_j + j
					t_astret N_idx = get_tmp_var();
					(*m_ostr) << "%" << N_idx->name << " = add i64 %" << N_idx_k->name << ", %" << ctr_j_val->name << "\n";

					// N[k,j] pointer
					t_astret elemptr_N_kj = get_tmp_var();
					(*m_ostr) << "%" << elemptr_N_kj->name << " = getelementptr [" << dim_k*dim_j << " x double], ["
						<< dim_k*dim_j << " x double]* %" << term2->name << ", i64 0, i64 %"
						<< N_idx->name << "\n";

					// N[k,l] value
					t_astret elem_N_kj = get_tmp_var();
					(*m_ostr) << "%" << elem_N_kj->name << " = load double, double* %" << elemptr_N_kj->name << "\n";


					// val = M[i,k] * N[k,j]
					t_astret val = get_tmp_var();
					(*m_ostr) << "%" << val->name << " = fmul double %" << elem_M_ik->name << ", %" << elem_N_kj->name << "\n";

					// d += val
					t_astret d_val_old = get_tmp_var();
					(*m_ostr) << "%" << d_val_old->name << " = load double, double* %" << d->name << "\n";
					t_astret d_val_new = get_tmp_var();
					(*m_ostr) << "%" << d_val_new->name << " = fadd double %" << d_val_old->name << ", %" << val->name << "\n";
					(*m_ostr) << "store double %" << d_val_new->name <<  ", double* %" << d->name << "\n";
				});

				// L index: i*dim_j + j
				t_astret L_idx = get_tmp_var();
				(*m_ostr) << "%" << L_idx->name << " = add i64 %" << L_idx_i->name << ", %" << ctr_j_val->name << "\n";

				// L[i,j] pointer
				t_astret elemptr_L_ij = get_tmp_var();
				(*m_ostr) << "%" << elemptr_L_ij->name << " = getelementptr [" << dim_i*dim_j << " x double], ["
					<< dim_i*dim_j << " x double]* %" << L_mem->name << ", i64 0, i64 %"
					<< L_idx->name << "\n";

				// L[i,j] = d
				t_astret d_val = get_tmp_var();
				(*m_ostr) << "%" << d_val->name << " = load double, double* %" << d->name << "\n";
				(*m_ostr) << "store double %" << d_val->name << ", double* %" << elemptr_L_ij->name << "\n";
			});
		});

		return L_mem;
	}

	// scalar-matrix/vector product
	else if((term1->ty == SymbolType::SCALAR || term1->ty == SymbolType::INT)
		&& (term2->ty == SymbolType::MATRIX || term2->ty == SymbolType::VECTOR))
	{
		if(ast->IsInverted())
		{
			throw std::runtime_error("ASTMult: Cannot divide scalar \"" + term1->name
				+ "\" by vector or matrix \"" + term2->name + "\".\n");
		}

		return scalar_matrix_prod(term1, term2, 1);
	}

	// scalar-matrix/vector product
	else if((term2->ty == SymbolType::SCALAR || term2->ty == SymbolType::INT)
		&& (term1->ty == SymbolType::MATRIX || term1->ty == SymbolType::VECTOR))
	{
		bool bMul = !ast->IsInverted();
		return scalar_matrix_prod(term2, term1, bMul);
	}

	// scalar types
	else
	{
		// cast if needed
		SymbolType ty = term1->ty;
		if(term1->ty==SymbolType::SCALAR || term2->ty==SymbolType::SCALAR)
			ty = SymbolType::SCALAR;
		t_astret var = get_tmp_var(ty, &term1->dims);

		term1 = convert_sym(term1, ty);
		term2 = convert_sym(term2, ty);


		std::string op = ast->IsInverted() ? "div" : "mul";
		if(ty == SymbolType::SCALAR)
			op = "f" + op;
		else if(ty == SymbolType::INT && ast->IsInverted())
			op = "s" + op;

		(*m_ostr) << "%" << var->name << " = " << op << " "
			<< LLAsm::get_type_name(ty) << " %" << term1->name << ", %" << term2->name << "\n";

		return var;
	}

	throw std::runtime_error("ASTMult: Invalid multiplication operation between \""
		+ term1->name + "\" and \"" + term2->name + "\".");
}


t_astret LLAsm::visit(const ASTMod* ast)
{
	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = ast->GetTerm2()->accept(this);

	// cast if needed
	SymbolType ty = term1->ty;
	if(term1->ty==SymbolType::SCALAR || term2->ty==SymbolType::SCALAR)
		ty = SymbolType::SCALAR;
	t_astret var = get_tmp_var(ty, &term1->dims);

	term1 = convert_sym(term1, ty);
	term2 = convert_sym(term2, ty);


	std::string op = "rem";
	if(ty == SymbolType::SCALAR)
		op = "f" + op;
	else if(ty == SymbolType::INT)
		op = "s" + op;

	(*m_ostr) << "%" << var->name << " = " << op << " "
		<< LLAsm::get_type_name(ty) << " %" << term1->name << ", %" << term2->name << "\n";

	return var;
}


t_astret LLAsm::visit(const ASTPow* ast)
{
	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = ast->GetTerm2()->accept(this);

	if(term1->ty == SymbolType::MATRIX)
	{
		// only integer powers are supported for matrix powers
		term2 = convert_sym(term2, SymbolType::INT);

		std::size_t dim1 = std::get<0>(term1->dims);
		std::size_t dim2 = std::get<1>(term1->dims);
		std::size_t dim = dim1*dim2;

		if(dim1 != dim2)
			throw std::runtime_error("ASTPow: Matrix power needs square matrix\n");

		// cast input matrix pointer to element pointer
		t_astret termptr = get_tmp_var();
		(*m_ostr) << "%" << termptr->name << " = bitcast ["
			<< dim << " x double]* %" << term1->name << " to double*\n";

		// allocate result matrix
		t_astret result_mem = get_tmp_var(SymbolType::MATRIX, &term1->dims);
		(*m_ostr) << "%" << result_mem->name << " = alloca [" << dim << " x double]\n";

		// cast result matrix pointer to element pointer
		t_astret result_ptr = get_tmp_var();
		(*m_ostr) << "%" << result_ptr->name << " = bitcast ["
			<< dim << " x double]* %" << result_mem->name << " to double*\n";

		// call runtime function
		t_astret result_status = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << result_status->name << " = call i64 @ext_power(double* %"
			<< termptr->name << ", double* %" << result_ptr->name
			<< ", i64 " << dim1 << ", i64 %" << term2->name << ")\n";

		// TODO: check result_status
		return result_mem;
	}

	// scalar types
	else
	{
		// cast if needed
		SymbolType ty = term1->ty;
		if(term1->ty==SymbolType::SCALAR || term2->ty==SymbolType::SCALAR)
			ty = SymbolType::SCALAR;
		t_astret var = get_tmp_var(ty, &term1->dims);

		term1 = convert_sym(term1, ty);
		term2 = convert_sym(term2, ty);

		(*m_ostr) << "%" << var->name << " = call double @pow("
			<< LLAsm::get_type_name(ty) << " %" << term1->name << ", "
			<< LLAsm::get_type_name(ty) << " %" << term2->name << ")\n";

		return var;
	}

	throw std::runtime_error("ASTPow: Invalid power operation between \""
		+ term1->name + "\" and \"" + term2->name + "\".");
}


t_astret LLAsm::visit(const ASTTransp* ast)
{
	t_astret term = ast->GetTerm()->accept(this);

	if(term->ty == SymbolType::MATRIX)
	{
		std::size_t dim1 = std::get<0>(term->dims);
		std::size_t dim2 = std::get<1>(term->dims);
		std::size_t dim = dim1*dim2;

		// allocate result matrix
		std::array<std::size_t, 2> dimtrans{{dim2, dim1}};
		t_astret result_mem = get_tmp_var(SymbolType::MATRIX, &dimtrans);
		(*m_ostr) << "%" << result_mem->name << " = alloca [" << dim << " x double]\n";


		// copy elements in a loop
		generate_loop(0, dim, [this, term, dim, dim1, dim2, result_mem](t_astret ctrval)
		{
			// loop statements
			// ctr_result = (ctr%cols)*rows + ctr/cols
			t_astret ctr_mod_cols = get_tmp_var(SymbolType::INT);
			t_astret ctr_div_cols = get_tmp_var(SymbolType::INT);
			t_astret ctr_mul_rows = get_tmp_var(SymbolType::INT);
			t_astret ctr_result = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << ctr_mod_cols->name << " = srem i64 %" << ctrval->name << ", " << dim2 << "\n";
			(*m_ostr) << "%" << ctr_div_cols->name << " = sdiv i64 %" << ctrval->name << ", " << dim2 << "\n";
			(*m_ostr) << "%" << ctr_mul_rows->name << " = mul i64 %" << ctr_mod_cols->name << ", " << dim1 << "\n";
			(*m_ostr) << "%" << ctr_result->name << " = add i64 %" << ctr_mul_rows->name
				<< ", %" << ctr_div_cols->name << "\n";

			// get matrix element pointers
			t_astret elemptr_term = get_tmp_var();
			(*m_ostr) << "%" << elemptr_term->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term->name << ", i64 0, i64 %" << ctrval->name << "\n";
			t_astret elemptr_result = get_tmp_var();
			(*m_ostr) << "%" << elemptr_result->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << result_mem->name << ", i64 0, i64 %" << ctr_result->name << "\n";

			t_astret elem_term = get_tmp_var();
			(*m_ostr) << "%" << elem_term->name << " = load double, double* %" << elemptr_term->name << "\n";

			// save result
			(*m_ostr) << "store double %" << elem_term->name << ", double* %" << elemptr_result->name << "\n";
		});

		return result_mem;
	}

	throw std::runtime_error("ASTTrans: Transposing is not possible for \"" + term->name + "\".");
}


t_astret LLAsm::visit(const ASTNorm* ast)
{
	t_astret term = ast->GetTerm()->accept(this);

	if(term->ty == SymbolType::SCALAR)
	{
		t_astret var = get_tmp_var(term->ty);
		(*m_ostr) << "%" << var->name << " = call double @fabs(double %" << term->name << ")\n";
		return var;
	}
	else if(term->ty == SymbolType::INT)
	{
		t_astret var = get_tmp_var(term->ty);
		(*m_ostr) << "%" << var->name << " = call i64 @labs(i64 %" << term->name << ")\n";
		return var;
	}
	else if(term->ty == SymbolType::VECTOR)
	{	// 2-norm
		std::size_t dim = std::get<0>(term->dims);

		// result variable
		t_astret dotptr = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << dotptr->name << " = alloca double\n";
		(*m_ostr) << "store double 0., double* %" << dotptr->name << "\n";

		// add elements in a loop
		generate_loop(0, dim, [this, term, dotptr, dim](t_astret ctrval)
		{
			// loop statements
			t_astret elemptr_src1 = get_tmp_var();
			(*m_ostr) << "%" << elemptr_src1->name << " = getelementptr [" << dim << " x double], ["
				<< dim << " x double]* %" << term->name << ", i64 0, i64 %" << ctrval->name << "\n";

			t_astret elem_src1 = get_tmp_var();
			(*m_ostr) << "%" << elem_src1->name << " = load double, double* %" << elemptr_src1->name << "\n";

			// square element
			t_astret mul_elems = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << mul_elems->name << " = fmul "
				<< "double %" << elem_src1->name << ", %" << elem_src1->name << "\n";

			// increment dot product variable
			t_astret cur_dot = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << cur_dot->name << " = load double, double* %" << dotptr->name << "\n";
			t_astret sum_dot = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << sum_dot->name << " = fadd "
				<< "double %" << cur_dot->name << ", %" << mul_elems->name << "\n";
			(*m_ostr) << "store double %" << sum_dot->name << ", double* %" << dotptr->name << "\n";
		});


		t_astret dot = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << dot->name << " = load double, double* %" << dotptr->name << "\n";

		t_astret dot_sqrt = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << dot_sqrt->name << " = call double @sqrt(double %" << dot->name << ")\n";

		return dot_sqrt;
	}
	else if(term->ty == SymbolType::MATRIX)
	{
		std::size_t dim1 = std::get<0>(term->dims);
		std::size_t dim2 = std::get<1>(term->dims);
		std::size_t dim = dim1*dim2;

		if(dim1 != dim2)
			throw std::runtime_error("ASTNorm: Determinant needs square matrix\n");

		// cast array pointer to element pointer
		t_astret termptr = get_tmp_var();
		(*m_ostr) << "%" << termptr->name << " = bitcast [" << dim << " x double]* %" << term->name << " to double*\n";

		t_astret det = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << det->name << " = call double @ext_determinant(double* %"
			<< termptr->name << ", i64 " << dim1 << ")\n";

		return det;
	}

	throw std::runtime_error("ASTNorm: Invalid symbol type for \"" + term->name + "\".");
}


t_astret LLAsm::visit(const ASTComp* ast)
{
	// code generation for equality test of two scalars
	auto scalars_equal = [this](const t_astret& term1, const t_astret& term2, ASTComp::CompOp op=ASTComp::EQU) -> t_astret
	{
		// get epsilon
		t_astret eps = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << eps->name << " = call double @get_eps()\n";

		// diff = term1 - term2
		t_astret diff = get_tmp_var(SymbolType::SCALAR);
		t_astret diff_abs_ptr = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << diff_abs_ptr->name << " = alloca double\n";

		(*m_ostr) << "%" << diff->name << " = fsub double %" << term1->name << ", %" << term2->name << "\n";

		// diff_neg = |diff|
		generate_cond([this, diff]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = fcmp olt double %" << diff->name << ", 0.\n";
			return cond;
		}, [this, diff, diff_abs_ptr]
		{
			t_astret diff_neg = get_tmp_var(SymbolType::SCALAR);
			(*m_ostr) << "%" << diff_neg->name << " = fneg double %" << diff->name << "\n";
			(*m_ostr) << "store double %" << diff_neg->name << ", double* %" << diff_abs_ptr->name << "\n";
		}, [this, diff, diff_abs_ptr]
		{
			(*m_ostr) << "store double %" << diff->name << ", double* %" << diff_abs_ptr->name << "\n";
		}, true);


		t_astret diff_abs = get_tmp_var(SymbolType::SCALAR);
		t_astret varbool = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << diff_abs->name << " = load double, double* %" << diff_abs_ptr->name << "\n";
		// diff_abs <= eps
		if(op == ASTComp::EQU)
		{
			(*m_ostr) << "%" << varbool->name << " = fcmp ole double %"
				<< diff_abs->name << ", %" << eps->name << "\n";
		}
		// diff_abs > eps
		else if(op == ASTComp::NEQ)
		{
			(*m_ostr) << "%" << varbool->name << " = fcmp ogt double %"
				<< diff_abs->name << ", %" << eps->name << "\n";
		}

		return varbool;
	};


	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = ast->GetTerm2()->accept(this);

	// string comparison
	if(term1->ty == SymbolType::STRING || term2->ty == SymbolType::STRING)
	{
		term1 = convert_sym(term1, SymbolType::STRING);
		term2 = convert_sym(term2, SymbolType::STRING);

		std::size_t maxdim = std::max(std::get<0>(term1->dims), std::get<0>(term2->dims));

		// get string pointer
		t_astret term1ptr = get_tmp_var(SymbolType::STRING);
		t_astret term2ptr = get_tmp_var(SymbolType::STRING);
		(*m_ostr) << "%" << term1ptr->name << " = bitcast [" << std::get<0>(term1->dims)
			<< " x i8]* %" << term1->name << " to i8*\n";
		(*m_ostr) << "%" << term2ptr->name << " = bitcast [" << std::get<0>(term2->dims)
			<< " x i8]* %" << term2->name << " to i8*\n";

		// compare strings
		t_astret strcmp = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << strcmp->name << " = call i32 @strncmp(i8* %"
			<< term1ptr->name << ", i8* %" << term2ptr->name << ", i64 " << maxdim << ")\n";

		std::string op;
		switch(ast->GetOp())
		{
			case ASTComp::EQU: op = "eq"; break;
			case ASTComp::NEQ: op = "ne"; break;
			case ASTComp::GT: op = "sgt"; break;
			case ASTComp::LT: op = "slt"; break;
			case ASTComp::GEQ: op = "sge"; break;
			case ASTComp::LEQ: op = "sle"; break;
		}

		t_astret cmp = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << cmp->name << " = icmp " << op << " i32 %" << strcmp->name << ", 0\n";
		return cmp;
	}

	// vector comparison
	else if(term1->ty == SymbolType::VECTOR && term2->ty == SymbolType::VECTOR &&
		(ast->GetOp() == ASTComp::EQU || ast->GetOp() == ASTComp::NEQ))
	{
		bool bNeq = (ast->GetOp() == ASTComp::NEQ);

		t_astret varbool = get_tmp_var(SymbolType::INT);
		t_astret varbool_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << varbool_ptr->name << " = alloca i1\n";
		(*m_ostr) << "store i1 1, i1* %" << varbool_ptr->name << "\n";

		std::size_t dim1 = std::get<0>(term1->dims);
		std::size_t dim2 = std::get<0>(term2->dims);

		// vectors are unequal if they are of different size
		if(dim1 != dim2)
		{
			if(bNeq)
				(*m_ostr) << "store i1 true, i1* %" << varbool_ptr->name << "\n";
			else
				(*m_ostr) << "store i1 false, i1* %" << varbool_ptr->name << "\n";

			(*m_ostr) << "%" << varbool->name << " = load i1, i1* %" << varbool_ptr->name << "\n";
			return varbool;
		}


		// compare elements in a loop
		// TODO: use generate_loop
		std::string labelStart = get_label();
		std::string labelBegin = get_label();
		std::string labelEnd = get_label();
		std::string labelCont = get_label();

		// loop counter
		t_astret ctr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctr->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctr->name << "\n";

		(*m_ostr) << "br label %" << labelStart << "\n";
		(*m_ostr) << labelStart << ":\n";

		// loop condition: ctr < dim
		t_astret ctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrval->name << " = load i64, i64* %" << ctr->name << "\n";

		t_astret cond = get_tmp_var();
		(*m_ostr) << "%" << cond->name << " = icmp slt i64 %" << ctrval->name <<  ", " << dim1 << "\n";
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelBegin << ", label %" << labelEnd << "\n";

		(*m_ostr) << labelBegin << ":\n";

		// ---------------
		// loop statements
		t_astret elemptr1 = get_tmp_var(SymbolType::SCALAR);
		t_astret elemptr2 = get_tmp_var(SymbolType::SCALAR);

		(*m_ostr) << "%" << elemptr1->name << " = getelementptr [" << dim1 << " x double], ["
			<< dim1 << " x double]* %" << term1->name << ", i64 0, i64 %" << ctrval->name << "\n";
		(*m_ostr) << "%" << elemptr2->name << " = getelementptr [" << dim2 << " x double], ["
			<< dim2 << " x double]* %" << term2->name << ", i64 0, i64 %" << ctrval->name << "\n";

		t_astret elem1 = get_tmp_var(SymbolType::SCALAR);
		t_astret elem2 = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elem1->name << " = load double, double* %" << elemptr1->name << "\n";
		(*m_ostr) << "%" << elem2->name << " = load double, double* %" << elemptr2->name << "\n";

		t_astret elems_equal = scalars_equal(elem1, elem2, ASTComp::EQU);
		(*m_ostr) << "store i1 %" << elems_equal->name << ", i1* %" << varbool_ptr->name << "\n";

		// break loop if not equal
		(*m_ostr) << "br i1 %" << elems_equal->name << ", label %" << labelCont << ", label %" << labelEnd << "\n";
		(*m_ostr) << labelCont << ":\n";


		// increment counter
		t_astret newctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << newctrval->name << " = add i64 %" << ctrval->name << ", 1\n";
		(*m_ostr) << "store i64 %" << newctrval->name << ", i64* %" << ctr->name << "\n";
		// ---------------

		(*m_ostr) << "br label %" << labelStart << "\n";
		(*m_ostr) << labelEnd << ":\n";


		t_astret elems_equal2 = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elems_equal2->name << " = load i1, i1* %" << varbool_ptr->name << "\n";

		if(bNeq)
		{
			// elems_equal = !elems_notequal
			t_astret elems_notequal = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << elems_notequal->name << " = xor i1 1, %" << elems_equal2->name << "\n";

			(*m_ostr) << "store i1 %" << elems_notequal->name << ", i1* %" << varbool_ptr->name << "\n";
		}

		(*m_ostr) << "%" << varbool->name << " = load i1, i1* %" << varbool_ptr->name << "\n";
		return varbool;
	}

	// matrix comparison
	else if(term1->ty == SymbolType::MATRIX && term2->ty == SymbolType::MATRIX && (ast->GetOp() == ASTComp::EQU || ast->GetOp() == ASTComp::NEQ))
	{
		bool bNeq = (ast->GetOp() == ASTComp::NEQ);

		t_astret varbool = get_tmp_var(SymbolType::INT);
		t_astret varbool_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << varbool_ptr->name << " = alloca i1\n";
		(*m_ostr) << "store i1 1, i1* %" << varbool_ptr->name << "\n";

		std::size_t dim1_1 = std::get<0>(term1->dims);
		std::size_t dim1_2 = std::get<1>(term1->dims);
		std::size_t dim2_1 = std::get<0>(term2->dims);
		std::size_t dim2_2 = std::get<1>(term2->dims);

		// matrices are unequal if they are of different size
		if(dim1_1 != dim2_1 || dim1_2 != dim2_2)
		{
			if(bNeq)
				(*m_ostr) << "store i1 true, i1* %" << varbool_ptr->name << "\n";
			else
				(*m_ostr) << "store i1 false, i1* %" << varbool_ptr->name << "\n";

			(*m_ostr) << "%" << varbool->name << " = load i1, i1* %" << varbool_ptr->name << "\n";
			return varbool;
		}


		// compare elements in a loop
		// TODO: use generate_loop
		std::size_t dim = dim1_1*dim1_2;
		std::string labelStart = get_label();
		std::string labelBegin = get_label();
		std::string labelEnd = get_label();
		std::string labelCont = get_label();

		// loop counter
		t_astret ctr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctr->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctr->name << "\n";

		(*m_ostr) << "br label %" << labelStart << "\n";
		(*m_ostr) << labelStart << ":\n";

		// loop condition: ctr < dim
		t_astret ctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrval->name << " = load i64, i64* %" << ctr->name << "\n";

		t_astret cond = get_tmp_var();
		(*m_ostr) << "%" << cond->name << " = icmp slt i64 %" << ctrval->name <<  ", " << dim << "\n";
		(*m_ostr) << "br i1 %" << cond->name << ", label %" << labelBegin << ", label %" << labelEnd << "\n";

		(*m_ostr) << labelBegin << ":\n";

		// ---------------
		// loop statements
		t_astret elemptr1 = get_tmp_var(SymbolType::SCALAR);
		t_astret elemptr2 = get_tmp_var(SymbolType::SCALAR);

		(*m_ostr) << "%" << elemptr1->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << term1->name << ", i64 0, i64 %" << ctrval->name << "\n";
		(*m_ostr) << "%" << elemptr2->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << term2->name << ", i64 0, i64 %" << ctrval->name << "\n";

		t_astret elem1 = get_tmp_var(SymbolType::SCALAR);
		t_astret elem2 = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elem1->name << " = load double, double* %" << elemptr1->name << "\n";
		(*m_ostr) << "%" << elem2->name << " = load double, double* %" << elemptr2->name << "\n";

		t_astret elems_equal = scalars_equal(elem1, elem2, ASTComp::EQU);
		(*m_ostr) << "store i1 %" << elems_equal->name << ", i1* %" << varbool_ptr->name << "\n";

		// break loop if not equal
		(*m_ostr) << "br i1 %" << elems_equal->name << ", label %" << labelCont << ", label %" << labelEnd << "\n";
		(*m_ostr) << labelCont << ":\n";


		// increment counter
		t_astret newctrval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << newctrval->name << " = add i64 %" << ctrval->name << ", 1\n";
		(*m_ostr) << "store i64 %" << newctrval->name << ", i64* %" << ctr->name << "\n";
		// ---------------

		(*m_ostr) << "br label %" << labelStart << "\n";
		(*m_ostr) << labelEnd << ":\n";


		t_astret elems_equal2 = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elems_equal2->name << " = load i1, i1* %" << varbool_ptr->name << "\n";

		if(bNeq)
		{
			// elems_equal = !elems_notequal
			t_astret elems_notequal = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << elems_notequal->name << " = xor i1 1, %" << elems_equal2->name << "\n";

			(*m_ostr) << "store i1 %" << elems_notequal->name << ", i1* %" << varbool_ptr->name << "\n";
		}

		(*m_ostr) << "%" << varbool->name << " = load i1, i1* %" << varbool_ptr->name << "\n";
		return varbool;
	}

	// scalar types
	else if((term1->ty == SymbolType::SCALAR || term1->ty == SymbolType::INT)
		&& (term2->ty == SymbolType::SCALAR || term2->ty == SymbolType::INT))
	{
		// cast if needed
		SymbolType ty = term1->ty;
		if(term1->ty==SymbolType::SCALAR || term2->ty==SymbolType::SCALAR)
			ty = SymbolType::SCALAR;

		term1 = convert_sym(term1, ty);
		term2 = convert_sym(term2, ty);


		// comparison within an epsilon range
		if(ty == SymbolType::SCALAR && (ast->GetOp() == ASTComp::EQU || ast->GetOp() == ASTComp::NEQ))
		{
			return scalars_equal(term1, term2, ast->GetOp());
		}

		// exact comparison
		else
		{
			std::string op;
			switch(ast->GetOp())
			{
				case ASTComp::EQU: op = "eq"; break;
				case ASTComp::NEQ: op = "ne"; break;
				case ASTComp::GT: op = "gt"; break;
				case ASTComp::LT: op = "lt"; break;
				case ASTComp::GEQ: op = "ge"; break;
				case ASTComp::LEQ: op = "le"; break;
			}

			std::string cmpop;
			switch(ty)
			{
				case SymbolType::SCALAR:
				{
					cmpop = "fcmp";
					op = "o" + op;
					break;
				}
				case SymbolType::INT:
				{
					cmpop = "icmp";
					if(op != "eq" && op != "ne")
						op = "s" + op;	// signed
					break;
				}
				default:
				{
					throw std::runtime_error("ASTComp: Invalid comparison between variables of type "
						+ LLAsm::get_type_name(ty) + ".");
					break;
				}
			}

			t_astret var = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << var->name << " = " << cmpop << " " << op << " "
				<< LLAsm::get_type_name(ty) << " %" << term1->name << ", %" << term2->name << "\n";
			return var;
		}
	}

	throw std::runtime_error("ASTComp: Invalid comparison of \"" + term1->name + "\" and \"" + term2->name + "\".");
}


t_astret LLAsm::visit(const ASTBool* ast)
{
	t_astret term1 = ast->GetTerm1()->accept(this);
	t_astret term2 = nullptr;
	if(ast->GetTerm2())
		term2 = ast->GetTerm2()->accept(this);

	t_astret ret = get_tmp_var(SymbolType::INT);
	switch(ast->GetOp())
	{
		case ASTBool::XOR:
		{
			(*m_ostr) << "%" << ret->name << " = xor i1 %" << term1->name
				<< ", %" << term2->name << "\n";
			break;
		}
		case ASTBool::OR:
		{
			(*m_ostr) << "%" << ret->name << " = or i1 %" << term1->name
				<< ", %" << term2->name << "\n";
			break;
		}
		case ASTBool::AND:
		{
			(*m_ostr) << "%" << ret->name << " = and i1 %" << term1->name
				<< ", %" << term2->name << "\n";
			break;
		}
		case ASTBool::NOT:
		{
			(*m_ostr) << "%" << ret->name << " = xor i1 1, %"
				<< term1->name << "\n";
			break;
		}
		default:
		{
			throw std::runtime_error("ASTBool: Invalid operation.");
		}
	}

	return ret;
}
