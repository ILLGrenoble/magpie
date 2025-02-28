/**
 * llvm three-address code generator -- arrays
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


t_astret LLAsm::visit(const ASTArrayAccess* ast)
{
	// get array indices
	t_astret num1 = ast->GetNum1()->accept(this);
	num1 = convert_sym(num1, SymbolType::INT);

	t_astret num2 = nullptr;
	t_astret num3 = nullptr;
	t_astret num4 = nullptr;

	if(ast->GetNum2())
	{
		num2 = ast->GetNum2()->accept(this);
		num2 = convert_sym(num2, SymbolType::INT);
	}

	if(ast->GetNum3())
	{
		num3 = ast->GetNum3()->accept(this);
		num3 = convert_sym(num3, SymbolType::INT);
	}

	if(ast->GetNum4())
	{
		num4 = ast->GetNum4()->accept(this);
		num4 = convert_sym(num4, SymbolType::INT);
	}

	t_astret term = ast->GetTerm()->accept(this);

	// single-element vector access
	if(term->ty == SymbolType::VECTOR && !ast->IsRanged12())
	{
		if(num2 || num3 || num4)	// no further arguments needed
			throw std::runtime_error("ASTArrayAccess: Invalid access operator for vector \"" + term->name + "\".");

		std::size_t dim = std::get<0>(term->dims);
		num1 = safe_array_index(num1, dim);

		// vector element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << term->name << ", i64 0, i64 %" << num1->name << "\n";

		// vector element
		t_astret elem = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elem->name << " = load double, double* %" << elemptr->name << "\n";
		return elem;
	}

	// single-element string access
	else if(term->ty == SymbolType::STRING && !ast->IsRanged12())
	{
		if(num2 || num3 || num4)	// no further arguments needed
			throw std::runtime_error("ASTArrayAccess: Invalid access operator for string \"" + term->name + "\".");

		std::size_t dim = std::get<0>(term->dims);
		num1 = safe_array_index(num1, dim);

		// string element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim << " x i8], ["
			<< dim << " x i8]* %" << term->name << ", i64 0, i64 %" << num1->name << "\n";

		// char at that pointer position
		t_astret elem = get_tmp_var();
		(*m_ostr) << "%" << elem->name << " = load i8, i8* %" << elemptr->name << "\n";

		// return string
		std::array<std::size_t, 2> retdims{{2, 1}};	// one char and a 0
		t_astret str_mem = get_tmp_var(SymbolType::STRING, &retdims);

		// allocate the return string's memory
		(*m_ostr) << "%" << str_mem->name << " = alloca [" << std::get<0>(retdims) << " x i8]\n";

		// store the char
		t_astret ptr0 = get_tmp_var();
		(*m_ostr) << "%" << ptr0->name << " = getelementptr [" << std::get<0>(retdims) << " x i8], ["
			<< std::get<0>(retdims) << " x i8]* %" << str_mem->name << ", i64 0, i64 0\n";

		(*m_ostr) << "store i8 %" << elem->name << ", i8* %" << ptr0->name << "\n";

		// add zero termination
		t_astret ptr1 = get_tmp_var();
		(*m_ostr) << "%" << ptr1->name << " = getelementptr [" << std::get<0>(retdims) << " x i8], ["
			<< std::get<0>(retdims) << " x i8]* %" << str_mem->name << ", i64 0, i64 1\n";

		(*m_ostr) << "store i8 0, i8* %" << ptr1->name << "\n";

		return str_mem;
	}

	// ranged vector and string access
	else if((term->ty == SymbolType::VECTOR || term->ty == SymbolType::STRING) && ast->IsRanged12())
	{
		if(num3 || num4)	// no further arguments needed
			throw std::runtime_error("ASTArrayAccess: Invalid access operator for \"" + term->name + "\".");

		std::string strty, strptrty;
		if(term->ty == SymbolType::VECTOR)
			strty = "double";
		else if(term->ty == SymbolType::STRING)
			strty = "i8";
		strptrty = strty + "*";

		std::size_t dim = std::get<0>(term->dims);
		num1 = safe_array_index(num1, dim);
		num2 = safe_array_index(num2, dim);

		// if num1 > num2, the loop has to be in the reversed order
		t_astret num2larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num1, num2]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num1->name
				<< ", %" << num2->name << "\n";
			return cond;
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num2larger_ptr->name << "\n";
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num2larger_ptr->name << "\n";
		}, true);

		t_astret num2larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger->name << " = load i1, i1* %"
			<< num2larger_ptr->name << "\n";

		// create a new vector for the results
		// for the moment with full static size, which will be truncated on assignment
		// TODO
		t_astret vec_mem = get_tmp_var(term->ty, &term->dims);
		(*m_ostr) << "%" << vec_mem->name << " = alloca [" << dim << " x " << strty << "]\n";

		// assign elements in a loop
		// source vector index counter
		t_astret ctrsrc = get_tmp_var(SymbolType::INT);
		t_astret ctrsrcval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsrc->name << " = alloca i64\n";
		(*m_ostr) << "store i64 %" << num1->name << ", i64* %" << ctrsrc->name << "\n";

		// destination vector index counter
		t_astret ctrdst = get_tmp_var(SymbolType::INT);
		t_astret ctrdstval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrdst->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctrdst->name << "\n";

		generate_loop([this, ctrsrc, ctrsrcval, ctrdst, ctrdstval, num2, num2larger]() -> t_astret
		{
			(*m_ostr) << "%" << ctrsrcval->name << " = load i64, i64* %"
				<< ctrsrc->name << "\n";

			(*m_ostr) << "%" << ctrdstval->name << " = load i64, i64* %"
				<< ctrdst->name << "\n";

			t_astret condptr = get_tmp_var();
			(*m_ostr) << "%" << condptr->name << " = alloca i1\n";

			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, condptr, num2, ctrsrcval]
			{
				// ctrsrcval <= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
					<< ctrsrcval->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << condptr->name << "\n";
			}, [this, condptr, num2, ctrsrcval]
			{
				// ctrsrcval >= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
					<< ctrsrcval->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << condptr->name << "\n";
			}, true);

			t_astret condval = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << condval->name << " = load i1, i1* %" << condptr->name << "\n";
			return condval;
		}, [this, &strty, &strptrty, ctrsrc, ctrsrcval, ctrdst, ctrdstval,
			term, dim, vec_mem, num2larger]
		{
			// source vector element pointer
			t_astret srcelemptr = get_tmp_var();
			(*m_ostr) << "%" << srcelemptr->name << " = getelementptr ["
				<< dim << " x "<< strty << "], [" << dim << " x " << strty << "]* %"
				<< term->name << ", i64 0, i64 %" << ctrsrcval->name << "\n";

			// destination vector element pointer
			t_astret dstelemptr = get_tmp_var();
			(*m_ostr) << "%" << dstelemptr->name << " = getelementptr ["
				<< dim << " x " << strty << "], [" << dim << " x " << strty << "]* %"
				<< vec_mem->name << ", i64 0, i64 %" << ctrdstval->name << "\n";

			// source vector element
			t_astret srcelem = get_tmp_var();
			(*m_ostr) << "%" << srcelem->name << " = load " << strty
				<< ", " << strptrty << " %" << srcelemptr->name << "\n";

			// store to destination vector element pointer
			(*m_ostr) << "store " << strty << " %" << srcelem->name
				<< ", " << strptrty << " %" << dstelemptr->name << "\n";

			// increment counters
			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, ctrsrcval, ctrsrc]
			{
				// ++ctrsrcval
				t_astret newctrsrcval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsrcval->name << " = add i64 %"
					<< ctrsrcval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsrcval->name << ", i64* %"
					<< ctrsrc->name << "\n";
			}, [this, ctrsrcval, ctrsrc]
			{
				// --ctrsrcval
				t_astret newctrsrcval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsrcval->name << " = sub i64 %"
					<< ctrsrcval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsrcval->name << ", i64* %"
					<< ctrsrc->name << "\n";
			}, true);

			t_astret newctrdstval = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << newctrdstval->name << " = add i64 %"
				<< ctrdstval->name << ", 1\n";
			(*m_ostr) << "store i64 %" << newctrdstval->name << ", i64* %"
				<< ctrdst->name << "\n";
		});

		if(term->ty == SymbolType::STRING)
		{
			// add zero termination
			t_astret _ctrdstval = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << _ctrdstval->name << " = load i64, i64* %"
				<< ctrdst->name << "\n";
			t_astret dstelemptr = get_tmp_var();
			(*m_ostr) << "%" << dstelemptr->name << " = getelementptr ["
				<< dim << " x " << strty << "], [" << dim << " x " << strty << "]* %"
				<< vec_mem->name << ", i64 0, i64 %" << _ctrdstval->name << "\n";

			(*m_ostr) << "store i8 0, i8* %" << dstelemptr->name << "\n";
		}

		return vec_mem;
	}

	// single-element matrix access
	else if(term->ty == SymbolType::MATRIX && !ast->IsRanged12() && !ast->IsRanged34())
	{
		if(!num2 || num3 || num4)	// second argument needed
			throw std::runtime_error("ASTArrayAccess: Invalid access operator for matrix \"" + term->name + "\".");

		std::size_t dim1 = std::get<0>(term->dims);
		std::size_t dim2 = std::get<1>(term->dims);
		num1 = safe_array_index(num1, dim1);
		num2 = safe_array_index(num2, dim2);

		// idx = num1*dim2 + num2
		t_astret idx1 = get_tmp_var();
		t_astret idx = get_tmp_var();
		(*m_ostr) << "%" << idx1->name << " = mul i64 %" << num1->name << ", " << dim2 << "\n";
		(*m_ostr) << "%" << idx->name << " = add i64 %" << idx1->name << ", %" << num2->name << "\n";

		// matrix element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim1*dim2 << " x double], ["
			<< dim1*dim2 << " x double]* %" << term->name << ", i64 0, i64 %" << idx->name << "\n";

		// matrix element
		t_astret elem = get_tmp_var(SymbolType::SCALAR);
		(*m_ostr) << "%" << elem->name << " = load double, double* %" << elemptr->name << "\n";
		return elem;
	}

	// ranged matrix access
	// TODO: allow mixed forms of the type mat[1, 2~3]
	else if(term->ty == SymbolType::MATRIX && (ast->IsRanged12() && ast->IsRanged34()))
	{
		std::size_t dim1 = std::get<0>(term->dims);
		std::size_t dim2 = std::get<1>(term->dims);

		num1 = safe_array_index(num1, dim1);
		num2 = safe_array_index(num2, dim1);
		num3 = safe_array_index(num3, dim2);
		num4 = safe_array_index(num4, dim2);

		// if num1 > num2, the loop has to be in the reversed order
		t_astret num2larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num1, num2]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num1->name
				<< ", %" << num2->name << "\n";
			return cond;
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num2larger_ptr->name << "\n";
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num2larger_ptr->name << "\n";
		}, true);

		t_astret num2larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger->name << " = load i1, i1* %"
			<< num2larger_ptr->name << "\n";

		// if num3 > num4, the loop has to be in the reversed order
		t_astret num4larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num4larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num3, num4]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num3->name
				<< ", %" << num4->name << "\n";
			return cond;
		}, [this, num4larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num4larger_ptr->name << "\n";
		}, [this, num4larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num4larger_ptr->name << "\n";
		}, true);

		t_astret num4larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num4larger->name << " = load i1, i1* %"
			<< num4larger_ptr->name << "\n";

		// create a new vector for the results
		// for the moment with full static size, which will be truncated on assignment
		// TODO
		std::array<std::size_t, 2> vecdims{{dim1*dim2, 1}};
		t_astret vec_mem = get_tmp_var(SymbolType::VECTOR, &vecdims);
		(*m_ostr) << "%" << vec_mem->name << " = alloca [" << dim1*dim2 << " x double]\n";

		// assign elements in a loop
		// source vector index counter
		t_astret ctrsrc1 = get_tmp_var(SymbolType::INT);
		t_astret ctrsrc1val = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsrc1->name << " = alloca i64\n";
		(*m_ostr) << "store i64 %" << num1->name << ", i64* %" << ctrsrc1->name << "\n";

		t_astret ctrsrc3 = get_tmp_var(SymbolType::INT);
		t_astret ctrsrc3val = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsrc3->name << " = alloca i64\n";

		// destination vector index counter
		t_astret ctrdst = get_tmp_var(SymbolType::INT);
		t_astret ctrdstval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrdst->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctrdst->name << "\n";

		generate_loop([this, ctrsrc1, ctrsrc1val, num2, num2larger]() -> t_astret
		{
			(*m_ostr) << "%" << ctrsrc1val->name << " = load i64, i64* %"
				<< ctrsrc1->name << "\n";

			t_astret cond1ptr = get_tmp_var();
			(*m_ostr) << "%" << cond1ptr->name << " = alloca i1\n";

			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, cond1ptr, num2, ctrsrc1val]
			{
				// ctrsrc1val <= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
					<< ctrsrc1val->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << cond1ptr->name << "\n";
			}, [this, cond1ptr, num2, ctrsrc1val]
			{
				// ctrsrc1val >= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
					<< ctrsrc1val->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << cond1ptr->name << "\n";
			}, true);

			t_astret cond1val = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond1val->name << " = load i1, i1* %"
				<< cond1ptr->name << "\n";
			return cond1val;
		}, [this, ctrsrc1, ctrsrc1val, ctrdst, ctrdstval, term,
			dim1, dim2, vec_mem, num2larger, num3, ctrsrc3, ctrsrc3val,
			num4larger, num4]
		{
			// -----------------------------------------------------------------
			// inner loop
			// -----------------------------------------------------------------
			(*m_ostr) << "store i64 %" << num3->name << ", i64* %" << ctrsrc3->name << "\n";

			generate_loop([this, ctrsrc3, ctrsrc3val, ctrdst, ctrdstval,
				num4, num4larger]() -> t_astret
			{
				(*m_ostr) << "%" << ctrsrc3val->name << " = load i64, i64* %"
					<< ctrsrc3->name << "\n";

				(*m_ostr) << "%" << ctrdstval->name << " = load i64, i64* %"
					<< ctrdst->name << "\n";

				t_astret cond3ptr = get_tmp_var();
				(*m_ostr) << "%" << cond3ptr->name << " = alloca i1\n";

				// if num3 > num4, the loop has to be in the reversed order
				generate_cond([num4larger]() -> t_astret
				{
					return num4larger;
				}, [this, cond3ptr, num4, ctrsrc3val]
				{
					// ctrsrc3val <= num4
					t_astret _cond = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
						<< ctrsrc3val->name <<  ", %" << num4->name << "\n";
					(*m_ostr) << "store i1 %" << _cond->name
						<< ", i1* %" << cond3ptr->name << "\n";
				}, [this, cond3ptr, num4, ctrsrc3val]
				{
					// ctrsrc3val >= num4
					t_astret _cond = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
						<< ctrsrc3val->name <<  ", %" << num4->name << "\n";
					(*m_ostr) << "store i1 %" << _cond->name
						<< ", i1* %" << cond3ptr->name << "\n";
				}, true);

				t_astret cond3val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << cond3val->name << " = load i1, i1* %"
					<< cond3ptr->name << "\n";
				return cond3val;
			}, [this, ctrsrc1val, ctrsrc3, ctrsrc3val, ctrdst, ctrdstval,
				term, dim1, dim2, vec_mem, num4larger]
			{
				// idx = num1*dim2 + num2
				t_astret idx1 = get_tmp_var();
				t_astret idx = get_tmp_var();
				(*m_ostr) << "%" << idx1->name << " = mul i64 %"
					<< ctrsrc1val->name << ", " << dim2 << "\n";
				(*m_ostr) << "%" << idx->name << " = add i64 %"
					<< idx1->name << ", %" << ctrsrc3val->name << "\n";

				// source matrix element pointer
				t_astret srcelemptr = get_tmp_var();
				(*m_ostr) << "%" << srcelemptr->name << " = getelementptr ["
					<< dim1*dim2 << " x double], [" << dim1*dim2 << " x double]* %"
					<< term->name << ", i64 0, i64 %" << idx->name << "\n";

				// destination matrix element pointer
				t_astret dstelemptr = get_tmp_var();
				(*m_ostr) << "%" << dstelemptr->name << " = getelementptr ["
					<< dim1*dim2 << " x double], [" << dim1*dim2 << " x double]* %"
					<< vec_mem->name << ", i64 0, i64 %" << ctrdstval->name << "\n";

				// source matrix element
				t_astret srcelem = get_tmp_var(SymbolType::SCALAR);
				(*m_ostr) << "%" << srcelem->name << " = load double, double* %"
					<< srcelemptr->name << "\n";

				// store to destination matrix element pointer
				(*m_ostr) << "store double %" << srcelem->name << ", double* %"
					<< dstelemptr->name << "\n";

				// increment counters
				// if num3 > num4, the loop has to be in the reversed order
				generate_cond([num4larger]() -> t_astret
				{
					return num4larger;
				}, [this, ctrsrc3val, ctrsrc3]
				{
					// ++ctrsrcval
					t_astret newctrsrc3val = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << newctrsrc3val->name << " = add i64 %"
						<< ctrsrc3val->name << ", 1\n";
					(*m_ostr) << "store i64 %" << newctrsrc3val->name << ", i64* %"
						<< ctrsrc3->name << "\n";
				}, [this, ctrsrc3val, ctrsrc3]
				{
					// --ctrsrcval
					t_astret newctrsrc3val = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << newctrsrc3val->name << " = sub i64 %"
						<< ctrsrc3val->name << ", 1\n";
					(*m_ostr) << "store i64 %" << newctrsrc3val->name << ", i64* %"
						<< ctrsrc3->name << "\n";
				}, true);

				t_astret newctrdstval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrdstval->name << " = add i64 %"
					<< ctrdstval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrdstval->name << ", i64* %"
					<< ctrdst->name << "\n";
			});
			// -----------------------------------------------------------------

			// increment counters
			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, ctrsrc1val, ctrsrc1]
			{
				// ++ctrsrcval
				t_astret newctrsrc1val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsrc1val->name << " = add i64 %"
					<< ctrsrc1val->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsrc1val->name << ", i64* %"
					<< ctrsrc1->name << "\n";
			}, [this, ctrsrc1val, ctrsrc1]
			{
				// --ctrsrcval
				t_astret newctrsrc1val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsrc1val->name << " = sub i64 %"
					<< ctrsrc1val->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsrc1val->name << ", i64* %"
					<< ctrsrc1->name << "\n";
			}, true);
		});

		return vec_mem;
	}

	throw std::runtime_error("ASTArrayAccess: Invalid array access to \"" + term->name + "\".");
}


t_astret LLAsm::visit(const ASTArrayAssign* ast)
{
	std::string var = ast->GetIdent();
	t_astret sym = get_sym(var);

	t_astret expr = ast->GetExpr()->accept(this);
	bool expr_is_array = expr->ty==SymbolType::VECTOR || expr->ty==SymbolType::MATRIX
		|| expr->ty==SymbolType::STRING;

	t_astret num1 = ast->GetNum1()->accept(this);
	num1 = convert_sym(num1, SymbolType::INT);

	t_astret num2 = nullptr;
	t_astret num3 = nullptr;
	t_astret num4 = nullptr;

	if(ast->GetNum2())
	{
		num2 = ast->GetNum2()->accept(this);
		num2 = convert_sym(num2, SymbolType::INT);
	}

	if(ast->GetNum3())
	{
		num3 = ast->GetNum3()->accept(this);
		num3 = convert_sym(num3, SymbolType::INT);
	}

	if(ast->GetNum4())
	{
		num4 = ast->GetNum4()->accept(this);
		num4 = convert_sym(num4, SymbolType::INT);
	}


	// single-element vector assignment
	if(sym->ty == SymbolType::VECTOR && !ast->IsRanged12())
	{
		// cast if needed
		if(expr->ty != SymbolType::SCALAR)
			expr = convert_sym(expr, SymbolType::SCALAR);

		if(num2)	// no second argument needed
			throw std::runtime_error("ASTArrayAssign: Invalid element assignment for vector \"" + sym->name + "\".");

		std::size_t dim = std::get<0>(sym->dims);
		num1 = safe_array_index(num1, dim);

		// vector element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim << " x double], ["
			<< dim << " x double]* %" << sym->name << ", i64 0, i64 %" << num1->name << "\n";

		// assign vector element
		(*m_ostr) << "store double %" << expr->name << ", double* %" << elemptr->name << "\n";
	}

	// single-element string assignment
	else if(sym->ty == SymbolType::STRING && !ast->IsRanged12())
	{
		// cast if needed
		expr = convert_sym(expr, SymbolType::STRING);

		if(num2)	// no second argument needed
			throw std::runtime_error("ASTArrayAssign: Invalid element assignment for string \"" + sym->name + "\".");

		// target string dimensions
		std::size_t dim = std::get<0>(sym->dims);
		// source string dimensions
		std::size_t dim_src = std::get<0>(expr->dims);

		num1 = safe_array_index(num1, dim);

		// source string element pointer
		t_astret elemptr_src = get_tmp_var();
		(*m_ostr) << "%" << elemptr_src->name << " = getelementptr [" << dim_src << " x i8], ["
			<< dim_src << " x i8]* %" << expr->name << ", i64 0, i64 0\n";

		// char at that pointer position
		t_astret elem_src = get_tmp_var();
		(*m_ostr) << "%" << elem_src->name << " = load i8, i8* %" << elemptr_src->name << "\n";

		// target string element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim << " x i8], ["
			<< dim << " x i8]* %" << sym->name << ", i64 0, i64 %" << num1->name << "\n";

		// assign string element
		(*m_ostr) << "store i8 %" << elem_src->name << ", i8* %" << elemptr->name << "\n";
	}

	// single-element matrix assignment
	else if(sym->ty == SymbolType::MATRIX && !ast->IsRanged12() && !ast->IsRanged34())
	{
		// cast if needed
		expr = convert_sym(expr, SymbolType::SCALAR);

		if(!num2)	// second argument needed
			throw std::runtime_error("ASTArrayAssign: Invalid element assignment for matrix \"" + sym->name + "\".");

		std::size_t dim1 = std::get<0>(sym->dims);
		std::size_t dim2 = std::get<1>(sym->dims);
		num1 = safe_array_index(num1, dim1);
		num2 = safe_array_index(num2, dim2);

		// idx = num1*dim2 + num2
		t_astret idx1 = get_tmp_var();
		t_astret idx = get_tmp_var();
		(*m_ostr) << "%" << idx1->name << " = mul i64 %" << num1->name << ", " << dim2 << "\n";
		(*m_ostr) << "%" << idx->name << " = add i64 %" << idx1->name << ", %" << num2->name << "\n";

		// matrix element pointer
		t_astret elemptr = get_tmp_var();
		(*m_ostr) << "%" << elemptr->name << " = getelementptr [" << dim1*dim2 << " x double], ["
			<< dim1*dim2 << " x double]* %" << sym->name << ", i64 0, i64 %" << idx->name << "\n";

		// assign matrix element
		(*m_ostr) << "store double %" << expr->name << ", double* %" << elemptr->name << "\n";
	}

	// ranged vector and string assignment
	else if((sym->ty == SymbolType::VECTOR || sym->ty == SymbolType::STRING) && ast->IsRanged12())
	{
		if(num3 || num4)	// no further arguments needed
			throw std::runtime_error("ASTArrayAssign: Invalid access operator for \"" + sym->name + "\".");

		std::string symty, symptrty;
		if(sym->ty == SymbolType::VECTOR)
			symty = "double";
		else if(sym->ty == SymbolType::STRING)
			symty = "i8";
		symptrty = symty + "*";

		std::size_t dim = std::get<0>(sym->dims);
		std::size_t exprdim = std::get<0>(expr->dims);
		num1 = safe_array_index(num1, dim);
		num2 = safe_array_index(num2, dim);

		// if num1 > num2, the loop has to be in the reversed order
		t_astret num2larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num1, num2]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num1->name
				<< ", %" << num2->name << "\n";
			return cond;
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num2larger_ptr->name << "\n";
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num2larger_ptr->name << "\n";
		}, true);

		t_astret num2larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger->name << " = load i1, i1* %"
			<< num2larger_ptr->name << "\n";

		// assign elements in a loop
		// source vector index counter
		t_astret ctrsym = get_tmp_var(SymbolType::INT);
		t_astret ctrsymval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsym->name << " = alloca i64\n";
		(*m_ostr) << "store i64 %" << num1->name << ", i64* %" << ctrsym->name << "\n";

		// destination vector index counter
		t_astret ctrexpr = get_tmp_var(SymbolType::INT);
		t_astret ctrexprval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrexpr->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctrexpr->name << "\n";

		generate_loop([this, ctrsym, ctrsymval, ctrexpr, ctrexprval, num2, num2larger]() -> t_astret
		{
			(*m_ostr) << "%" << ctrsymval->name << " = load i64, i64* %"
				<< ctrsym->name << "\n";

			(*m_ostr) << "%" << ctrexprval->name << " = load i64, i64* %"
				<< ctrexpr->name << "\n";

			t_astret condptr = get_tmp_var();
			(*m_ostr) << "%" << condptr->name << " = alloca i1\n";

			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, condptr, num2, ctrsymval]
			{
				// ctrsymval <= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
					<< ctrsymval->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << condptr->name << "\n";
			}, [this, condptr, num2, ctrsymval]
			{
				// ctrsymval >= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
					<< ctrsymval->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << condptr->name << "\n";
			}, true);

			t_astret condval = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << condval->name << " = load i1, i1* %" << condptr->name << "\n";
			return condval;
		}, [this, &symty, &symptrty, ctrsym, ctrsymval, ctrexpr, ctrexprval,
			expr, dim, exprdim, expr_is_array, sym, num2larger]
		{
			t_astret srcelem = nullptr;
			if(expr_is_array)
			{
				// source vector element pointer
				t_astret srcelemptr = get_tmp_var();
				(*m_ostr) << "%" << srcelemptr->name << " = getelementptr ["
					<< exprdim << " x "<< symty << "], [" << exprdim << " x " << symty << "]* %"
					<< expr->name << ", i64 0, i64 %" << ctrexprval->name << "\n";

				// source vector element
				srcelem = get_tmp_var(get_element_type(expr->ty));
				(*m_ostr) << "%" << srcelem->name << " = load " << symty
					<< ", " << symptrty << " %" << srcelemptr->name << "\n";

				// conversion
				srcelem = convert_sym(srcelem, get_element_type(sym->ty));
			}
			else
			{
				// single source value
				srcelem = convert_sym(expr, get_element_type(sym->ty));
			}

			// destination vector element pointer
			t_astret dstelemptr = get_tmp_var();
			(*m_ostr) << "%" << dstelemptr->name << " = getelementptr ["
				<< dim << " x " << symty << "], [" << dim << " x " << symty << "]* %"
				<< sym->name << ", i64 0, i64 %" << ctrsymval->name << "\n";

			// store to destination vector element pointer
			(*m_ostr) << "store " << symty << " %" << srcelem->name
				<< ", " << symptrty << " %" << dstelemptr->name << "\n";

			// increment counters
			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, ctrsymval, ctrsym]
			{
				// ++ctrsymval
				t_astret newctrsymval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsymval->name << " = add i64 %"
					<< ctrsymval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsymval->name << ", i64* %"
					<< ctrsym->name << "\n";
			}, [this, ctrsymval, ctrsym]
			{
				// --ctrsymval
				t_astret newctrsymval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsymval->name << " = sub i64 %"
					<< ctrsymval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsymval->name << ", i64* %"
					<< ctrsym->name << "\n";
			}, true);

			t_astret newctrexprval = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << newctrexprval->name << " = add i64 %"
				<< ctrexprval->name << ", 1\n";
			(*m_ostr) << "store i64 %" << newctrexprval->name << ", i64* %"
				<< ctrexpr->name << "\n";
		});
	}

	// ranged matrix assignment
	else if(sym->ty == SymbolType::MATRIX && (ast->IsRanged12() && ast->IsRanged34()))
	{
		std::size_t dim1 = std::get<0>(sym->dims);
		std::size_t dim2 = std::get<1>(sym->dims);
		std::size_t exprdim1 = std::get<0>(expr->dims);
		std::size_t exprdim2 = std::get<1>(expr->dims);

		num1 = safe_array_index(num1, dim1);
		num2 = safe_array_index(num2, dim1);
		num3 = safe_array_index(num3, dim2);
		num4 = safe_array_index(num4, dim2);

		// if num1 > num2, the loop has to be in the reversed order
		t_astret num2larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num1, num2]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num1->name
				<< ", %" << num2->name << "\n";
			return cond;
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num2larger_ptr->name << "\n";
		}, [this, num2larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num2larger_ptr->name << "\n";
		}, true);

		t_astret num2larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num2larger->name << " = load i1, i1* %"
			<< num2larger_ptr->name << "\n";

		// if num3 > num4, the loop has to be in the reversed order
		t_astret num4larger_ptr = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num4larger_ptr->name << " = alloca i1\n";

		generate_cond([this, num3, num4]() -> t_astret
		{
			t_astret cond = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond->name << " = icmp sle i64 %" << num3->name
				<< ", %" << num4->name << "\n";
			return cond;
		}, [this, num4larger_ptr]
		{
			(*m_ostr) << "store i1 1, i1* %" << num4larger_ptr->name << "\n";
		}, [this, num4larger_ptr]
		{
			(*m_ostr) << "store i1 0, i1* %" << num4larger_ptr->name << "\n";
		}, true);

		t_astret num4larger = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << num4larger->name << " = load i1, i1* %"
			<< num4larger_ptr->name << "\n";

		// assign elements in a loop
		// symbol vector index counter
		t_astret ctrsym1 = get_tmp_var(SymbolType::INT);
		t_astret ctrsym1val = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsym1->name << " = alloca i64\n";
		(*m_ostr) << "store i64 %" << num1->name << ", i64* %" << ctrsym1->name << "\n";

		t_astret ctrsym3 = get_tmp_var(SymbolType::INT);
		t_astret ctrsym3val = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrsym3->name << " = alloca i64\n";

		// expression vector index counter
		t_astret ctrexpr = get_tmp_var(SymbolType::INT);
		t_astret ctrexprval = get_tmp_var(SymbolType::INT);
		(*m_ostr) << "%" << ctrexpr->name << " = alloca i64\n";
		(*m_ostr) << "store i64 0, i64* %" << ctrexpr->name << "\n";

		generate_loop([this, ctrsym1, ctrsym1val, num2, num2larger]() -> t_astret
		{
			(*m_ostr) << "%" << ctrsym1val->name << " = load i64, i64* %"
				<< ctrsym1->name << "\n";

			t_astret cond1ptr = get_tmp_var();
			(*m_ostr) << "%" << cond1ptr->name << " = alloca i1\n";

			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, cond1ptr, num2, ctrsym1val]
			{
				// ctrsym1val <= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
					<< ctrsym1val->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << cond1ptr->name << "\n";
			}, [this, cond1ptr, num2, ctrsym1val]
			{
				// ctrsym1val >= num2
				t_astret _cond = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
					<< ctrsym1val->name <<  ", %" << num2->name << "\n";
				(*m_ostr) << "store i1 %" << _cond->name
					<< ", i1* %" << cond1ptr->name << "\n";
			}, true);

			t_astret cond1val = get_tmp_var(SymbolType::INT);
			(*m_ostr) << "%" << cond1val->name << " = load i1, i1* %"
				<< cond1ptr->name << "\n";
			return cond1val;
		}, [this, ctrsym1, ctrsym1val, ctrexpr, ctrexprval, expr,
			dim1, dim2, exprdim1, exprdim2, expr_is_array, sym, num2larger,
			num3, ctrsym3, ctrsym3val, num4larger, num4]
		{
			// -----------------------------------------------------------------
			// inner loop
			// -----------------------------------------------------------------
			(*m_ostr) << "store i64 %" << num3->name << ", i64* %" << ctrsym3->name << "\n";

			generate_loop([this, ctrsym3, ctrsym3val, ctrexpr, ctrexprval,
				num4, num4larger]() -> t_astret
			{
				(*m_ostr) << "%" << ctrsym3val->name << " = load i64, i64* %"
					<< ctrsym3->name << "\n";

				(*m_ostr) << "%" << ctrexprval->name << " = load i64, i64* %"
					<< ctrexpr->name << "\n";

				t_astret cond3ptr = get_tmp_var();
				(*m_ostr) << "%" << cond3ptr->name << " = alloca i1\n";

				// if num3 > num4, the loop has to be in the reversed order
				generate_cond([num4larger]() -> t_astret
				{
					return num4larger;
				}, [this, cond3ptr, num4, ctrsym3val]
				{
					// ctrsym3val <= num4
					t_astret _cond = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << _cond->name << " = icmp sle i64 %"
						<< ctrsym3val->name <<  ", %" << num4->name << "\n";
					(*m_ostr) << "store i1 %" << _cond->name
						<< ", i1* %" << cond3ptr->name << "\n";
				}, [this, cond3ptr, num4, ctrsym3val]
				{
					// ctrsym3val >= num4
					t_astret _cond = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << _cond->name << " = icmp sge i64 %"
						<< ctrsym3val->name <<  ", %" << num4->name << "\n";
					(*m_ostr) << "store i1 %" << _cond->name
						<< ", i1* %" << cond3ptr->name << "\n";
				}, true);

				t_astret cond3val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << cond3val->name << " = load i1, i1* %"
					<< cond3ptr->name << "\n";
				return cond3val;
			}, [this, ctrsym1val, ctrsym3, ctrsym3val, ctrexpr, ctrexprval,
				expr, dim1, dim2, exprdim1, exprdim2, expr_is_array, sym, num4larger]
			{
				// idx = num1*dim2 + num2
				t_astret idx1 = get_tmp_var();
				t_astret idx = get_tmp_var();
				(*m_ostr) << "%" << idx1->name << " = mul i64 %"
					<< ctrsym1val->name << ", " << dim2 << "\n";
				(*m_ostr) << "%" << idx->name << " = add i64 %"
					<< idx1->name << ", %" << ctrsym3val->name << "\n";

				t_astret srcelem = nullptr;
				if(expr_is_array)
				{
					// source matrix element pointer
					t_astret srcelemptr = get_tmp_var();
					(*m_ostr) << "%" << srcelemptr->name << " = getelementptr ["
						<< exprdim1*exprdim2 << " x double], [" << exprdim1*exprdim2 << " x double]* %"
						<< expr->name << ", i64 0, i64 %" << ctrexprval->name << "\n";

					// source matrix element
					srcelem = get_tmp_var(get_element_type(expr->ty));
					(*m_ostr) << "%" << srcelem->name << " = load double, double* %"
						<< srcelemptr->name << "\n";

					// conversion
					srcelem = convert_sym(srcelem, get_element_type(sym->ty));
				}
				else
				{
					// single source value
					srcelem = convert_sym(expr, SymbolType::SCALAR);
				}

				// destination matrix element pointer
				t_astret dstelemptr = get_tmp_var();
				(*m_ostr) << "%" << dstelemptr->name << " = getelementptr ["
					<< dim1*dim2 << " x double], [" << dim1*dim2 << " x double]* %"
					<< sym->name << ", i64 0, i64 %" << idx->name << "\n";

				// store to destination matrix element pointer
				(*m_ostr) << "store double %" << srcelem->name << ", double* %"
					<< dstelemptr->name << "\n";

				// increment counters
				// if num3 > num4, the loop has to be in the reversed order
				generate_cond([num4larger]() -> t_astret
				{
					return num4larger;
				}, [this, ctrsym3val, ctrsym3]
				{
					// ++ctrsym3
					t_astret newctrsym3val = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << newctrsym3val->name << " = add i64 %"
						<< ctrsym3val->name << ", 1\n";
					(*m_ostr) << "store i64 %" << newctrsym3val->name << ", i64* %"
						<< ctrsym3->name << "\n";
				}, [this, ctrsym3val, ctrsym3]
				{
					// --ctrsym3
					t_astret newctrsym3val = get_tmp_var(SymbolType::INT);
					(*m_ostr) << "%" << newctrsym3val->name << " = sub i64 %"
						<< ctrsym3val->name << ", 1\n";
					(*m_ostr) << "store i64 %" << newctrsym3val->name << ", i64* %"
						<< ctrsym3->name << "\n";
				}, true);

				t_astret newctrexprval = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrexprval->name << " = add i64 %"
					<< ctrexprval->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrexprval->name << ", i64* %"
					<< ctrexpr->name << "\n";
			});
			// -----------------------------------------------------------------

			// increment counters
			// if num1 > num2, the loop has to be in the reversed order
			generate_cond([num2larger]() -> t_astret
			{
				return num2larger;
			}, [this, ctrsym1val, ctrsym1]
			{
				// ++ctrsym1
				t_astret newctrsym1val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsym1val->name << " = add i64 %"
					<< ctrsym1val->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsym1val->name << ", i64* %"
					<< ctrsym1->name << "\n";
			}, [this, ctrsym1val, ctrsym1]
			{
				// --ctrsym1
				t_astret newctrsym1val = get_tmp_var(SymbolType::INT);
				(*m_ostr) << "%" << newctrsym1val->name << " = sub i64 %"
					<< ctrsym1val->name << ", 1\n";
				(*m_ostr) << "store i64 %" << newctrsym1val->name << ", i64* %"
					<< ctrsym1->name << "\n";
			}, true);
		});
	}

	return expr;
}
