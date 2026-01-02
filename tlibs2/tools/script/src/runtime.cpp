/**
 * runtime library using the tlibs2 math library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 17-apr-20
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

#include <cstdint>
#include <cfloat>
#include <vector>
#include <unordered_set>
#include <iostream>

#include "maths.h"
using namespace tl2_ops;


using t_real = double;
using t_int = int64_t;
using t_byte = int8_t;

using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;

static t_real g_eps = DBL_EPSILON;
static int8_t g_debug = 0;


extern "C" {

// ----------------------------------------------------------------------------
// heap management
// ----------------------------------------------------------------------------
static std::unordered_set<void*> heap;


void* ext_heap_alloc(uint64_t num, uint64_t elemsize)
{
	void *mem = new uint8_t[num * elemsize];
	heap.insert(mem);

	if(g_debug)
	{
		std::cerr << __func__
			<< ": count=" << num
			<< ", elem_size=" << elemsize
			<< ", mem=" << std::setw(8) << std::setfill('0') << std::hex << mem
			<< std::endl;
	}

	return mem;
}


void ext_heap_free(void* mem)
{
	if(!mem)
		return;

	delete[] reinterpret_cast<uint8_t*>(mem);
	heap.erase(mem);

	if(g_debug)
	{
		std::cerr << __func__
			<< ": mem=" << std::setw(8) << std::setfill('0') << std::hex << mem
			<< std::endl;
	}
}


void ext_init()
{
	heap.clear();
}


void ext_deinit()
{
	if(g_debug)
	{
		std::cerr << __func__ << ": " << heap.size() << " memory leaks detected." << std::endl;

		auto iter = heap.begin();
		for(std::size_t i=0; i<heap.size(); ++i)
		{
			std::cerr << "Leak " << i << ": "
				<< ", mem=" << std::setw(8) << std::setfill('0') << std::hex << (void*)*iter
				<< std::endl;
			++iter;
		}
	}
}


void set_debug(t_int dbg)
{
	g_debug = (dbg!=0);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// mathematical functions
// ----------------------------------------------------------------------------

/**
 * set float epsilon
 */
void set_eps(t_real eps)
{
	g_eps = eps;
}


/**
 * get float epsilon
 */
t_real get_eps()
{
	return g_eps;
}


/**
 * tests equality of floating point numbers
 */
int ext_equals(t_real x, t_real y, t_real eps)
{
	return tl2::equals<t_real>(x, y, eps);
}


/**
 * removes a given row and column of a square matrix
 */
void ext_submat(const t_real* M, t_int N, t_real* M_new, t_int iremove, t_int jremove)
{
	t_vec mat(N*N, M);
	t_vec submat = tl2::flat_submat<t_vec>(mat, N, N, iremove, jremove);
	submat.to_array(M_new);
}


/**
 * calculates the determinant
 */
t_real ext_determinant(const t_real* M, t_int N)
{
	t_mat mat(N, N, M);
	return tl2::det<t_mat>(mat);
}


/**
 * inverted matrix
 */
t_int ext_inverse(const t_real* M, t_real* I, t_int N)
{
	t_real fullDet = ext_determinant(M, N);

	// fail if determinant is zero
	if(ext_equals(fullDet, 0., g_eps))
		return 0;

	t_mat mat(N, N, M);
	auto [inv, ok] = tl2::inv<t_mat>(mat);
	inv.to_array(I);

	return ok==true;
}


/**
 * matrix power
 */
t_int ext_power(const t_real* M, t_real* P, t_int N, t_int POW)
{
	t_mat mat(N, N, M);
	auto [result, ok] = tl2::pow(mat, POW);
	result.to_array(P);

	return ok==true;
}


/**
 * transposed matrix
 */
void ext_transpose(const t_real* M, t_real* T, t_int rows, t_int cols)
{
	for(t_int ctr=0; ctr<rows*cols; ++ctr)
	{
		t_int i = ctr/cols;
		t_int j = ctr%cols;
		T[j*rows + i] = M[i*cols + j];
	}
}


/**
 * eigenvalues
 */
t_byte* eigenvals(const t_real* M, t_int rows, t_int cols)
{
	t_mat mat(rows, cols, M);
	set_eps_0(mat, g_eps);
	bool is_symm = tl2::is_symm_or_herm<t_mat>(mat, g_eps);

	auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
		tl2_la::eigenvec<t_mat, t_vec, t_real>(mat, true, is_symm, true);

	t_int N = rows;
	t_real* mem = reinterpret_cast<t_real*>(ext_heap_alloc(2*N, sizeof(t_real)));

	// store eigenvalues
	for(t_int i=0; i<N; ++i)
	{
		mem[i] = evals_re[i];
		mem[i + N] = evals_im[i];
	}

	return reinterpret_cast<t_byte*>(mem);
}


/**
 * eigenvalues and -vectors
 */
t_byte* eigenvecs(const t_real* M, t_int rows, t_int cols)
{
	t_mat mat(rows, cols, M);
	set_eps_0(mat, g_eps);

	bool is_symm = tl2::is_symm_or_herm<t_mat>(mat, g_eps);
	//is_symm = false;
	bool normalise = true;

	//using namespace tl2_ops;
	//std::cout << std::scientific << "mat = " << mat << ", symm: " << std::boolalpha << is_symm << std::endl;
	auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
		tl2_la::eigenvec<t_mat, t_vec, t_real>(mat, false, is_symm, normalise);

	t_int N = rows;
	t_real* mem = reinterpret_cast<t_real*>(ext_heap_alloc(2*N + 2*N*N, sizeof(t_real)));

	// store eigenvalues
	for(t_int i=0; i<N; ++i)
	{
		mem[i] = evals_re[i];
		mem[i + N] = evals_im[i];
	}

	// store eigenvectors
	t_int memidx = 2*N;
	for(t_int i=0; i<N; ++i)
	{
		for(t_int j=0; j<N; ++j)
		{
			mem[memidx] = evecs_re[i][j];
			mem[memidx + N*N] = evecs_im[i][j];

			++memidx;
		}
	}

	return reinterpret_cast<t_byte*>(mem);
}


/**
 * qr decomposition
 */
t_byte* qr(const t_real* M, t_int rows, t_int cols)
{
	t_mat mat(rows, cols, M);
	auto [ok, matQ, matR] = tl2::qr<t_mat, t_vec>(mat);

	t_real* mem = reinterpret_cast<t_real*>(ext_heap_alloc(matQ.size1()*matQ.size2() + matR.size1()*matR.size2(), sizeof(t_real)));
	matQ.to_array(mem);
	matR.to_array(mem + matQ.size1()*matQ.size2());

	return reinterpret_cast<t_byte*>(mem);
}
// ----------------------------------------------------------------------------


}
