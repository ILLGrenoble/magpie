/**
 * tlibs2 maths library -- complex algorithms
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_MATHS_CPLX_H__
#define __TLIBS2_MATHS_CPLX_H__

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

#include "decls.h"
#include "ndim.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// complex algorithms
// ----------------------------------------------------------------------------

/**
 * split a complex vector into two vectors with the real and imag parts
 */
template<class t_vec_cplx, class t_vec_real>
std::tuple<t_vec_real, t_vec_real> split_cplx(const t_vec_cplx& vec)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	using t_real = typename t_vec_real::value_type;
	const auto vec_size = vec.size();

	t_vec_real vecRe = zero<t_vec_real>(vec_size);
	t_vec_real vecIm = zero<t_vec_real>(vec_size);

	std::size_t iter = 0;
	std::size_t iterRe = 0;
	std::size_t iterIm = 0;

	for(; iter < vec_size; )
	{
		vecRe[iterRe] = t_real{vec[iter].real()};
		vecIm[iterIm] = t_real{vec[iter].imag()};

		++iter;
		++iterRe;
		++iterIm;
	}

	return std::make_tuple(vecRe, vecIm);
}


/**
 * split a complex matrix into two matrices with the real and imag parts
 */
template<class t_mat_cplx, class t_mat_real>
std::tuple<t_mat_real, t_mat_real> split_cplx(const t_mat_cplx& mat)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	const auto mat_size1 = mat.size1();
	const auto mat_size2 = mat.size2();

	t_mat_real matRe = zero<t_mat_real>(mat_size1, mat_size2);
	t_mat_real matIm = zero<t_mat_real>(mat_size1, mat_size2);

	for(std::size_t i = 0; i < mat_size1; ++i)
	{
		for(std::size_t j = 0; j < mat_size2; ++j)
		{
			matRe(i,j) = mat(i,j).real();
			matIm(i,j) = mat(i,j).imag();
		}
	}

	return std::make_tuple(matRe, matIm);
}


/**
 * split a complex vector into a real vector with the real and imag components
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_real split_cplx_real(const t_vec_cplx& vec)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	const auto vec_size = vec.size();

	using t_real = typename t_vec_real::value_type;
	t_vec_real vecReal = zero<t_vec_real>(vec_size * 2);

	std::size_t iter = 0;
	std::size_t iterRe = 0;
	std::size_t iterIm = vec_size;

	for(; iter < vec_size; )
	{
		vecReal[iterRe] = t_real{vec[iter].real()};
		vecReal[iterIm] = t_real{vec[iter].imag()};

		++iter;
		++iterRe;
		++iterIm;
	}

	return vecReal;
}


/**
 * split a complex matrix into a real matrix with the real and imag components
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_real split_cplx_real(const t_mat_cplx& mat)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	const auto mat_size1 = mat.size1();
	const auto mat_size2 = mat.size2();

	t_mat_real matRe = zero<t_mat_real>(mat_size1*2, mat_size2*2);

	for(std::size_t i = 0; i < mat_size1; ++i)
	{
		for(std::size_t j = 0; j < mat_size2; ++j)
		{
			matRe(i, j) = mat(i,j).real();
			matRe(i + mat_size1, j + mat_size2) = mat(i,j).imag();
		}
	}

	return matRe;
}


/**
 * create a complex vector from two real vectors
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_cplx unite_cplx(const t_vec_real& vecRe, const t_vec_real& vecIm)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	const auto vecRe_size = vecRe.size();

	assert(vecRe.size() == vecIm.size());
	t_vec_cplx vec = zero<t_vec_cplx>(vecRe_size);

	std::size_t iter = 0;
	std::size_t iterRe = 0;
	std::size_t iterIm = 0;

	for(; iterRe < vecRe_size; )
	{
		vec[iter].real(vecRe[iterRe]);
		vec[iter].imag(vecIm[iterIm]);

		++iter;
		++iterRe;
		++iterIm;
	}

	return vec;
}


/**
 * create a complex matrix from two real matrices
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_cplx unite_cplx(const t_mat_real& matRe, const t_mat_real& matIm)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	const auto matRe_size1 = matRe.size1();
	const auto matRe_size2 = matRe.size2();

	assert(matRe_size1 == matIm.size1());
	assert(matRe_size2 == matIm.size2());

	t_mat_cplx mat = zero<t_mat_cplx>(matRe_size1, matRe_size2);

	for(std::size_t i = 0; i < matRe_size1; ++i)
	{
		for(std::size_t j = 0; j < matRe_size2; ++j)
		{
			mat(i, j).real(matRe(i, j));
			mat(i, j).imag(matIm(i, j));
		}
	}

	return mat;
}


/**
 * create a complex vector from a real vectors having real and imag components
 */
template<class t_vec_cplx, class t_vec_real>
t_vec_cplx unite_cplx_real(const t_vec_real& vecReIm)
requires is_complex<typename t_vec_cplx::value_type> && is_vec<t_vec_cplx> && is_vec<t_vec_real>
{
	const auto vec_size = vecReIm.size() / 2;
	t_vec_cplx vec = zero<t_vec_cplx>(vec_size);

	std::size_t iter = 0;
	std::size_t iterRe = 0;
	std::size_t iterIm = vec_size;

	for(; iter < vec_size; )
	{
		vec[iter].real(vecReIm[iterRe]);
		vec[iter].imag(vecReIm[iterIm]);

		++iter;
		++iterRe;
		++iterIm;
	}

	return vec;
}


/**
 * create a complex matrix from a real matrix having real and imag components
 */
template<class t_mat_cplx, class t_mat_real>
t_mat_cplx unite_cplx_real(const t_mat_real& matReIm)
requires is_complex<typename t_mat_cplx::value_type> && is_mat<t_mat_cplx> && is_mat<t_mat_real>
{
	const auto mat_size1 = matReIm.size1() / 2;
	const auto mat_size2 = matReIm.size2() / 2;
	t_mat_cplx mat = zero<t_mat_cplx>(mat_size1, mat_size2);

	for(std::size_t i = 0; i < mat_size1; ++i)
	{
		for(std::size_t j = 0; j < mat_size2; ++j)
		{
			mat(i, j).real(matReIm(i, j));
			mat(i, j).imag(matReIm(i + mat_size1, j + mat_size2));
		}
	}

	return mat;
}


/**
 * SU(2) generators, pauli matrices sig_i = 2*S_i
 * @see (Arfken 2013), p. 110
 */
template<class t_mat>
const t_mat& su2_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx c0(0, 0);
	constexpr t_cplx c1(1, 0);
	constexpr t_cplx cI(0, 1);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0, c1}, { c1,  c0}}),  // x
		create<t_mat>({{c0, cI}, {-cI,  c0}}),  // y
		create<t_mat>({{c1, c0}, { c0, -c1}}),  // z
	};

	return mat[which];
}


/**
 * get a vector of pauli matrices
 * @see (Arfken 2013), p. 110
 */
template<class t_vec>
t_vec su2_matrices(bool bIncludeUnit = false)
requires is_basic_vec<t_vec> && is_mat<typename t_vec::value_type>
	&& is_complex<typename t_vec::value_type::value_type>
{
	using t_mat = typename t_vec::value_type;

	t_vec vec;
	vec.reserve(4);

	if(bIncludeUnit)
		vec.emplace_back(unit<t_mat>(2));
	for(std::size_t i=0; i<3; ++i)
		vec.emplace_back(su2_matrix<t_mat>(i));

	return vec;
}


/**
 * project the vector of SU(2) matrices onto a vector
 * proj = <sigma|vec>
 */
template<class t_vec, class t_mat>
t_mat proj_su2(const t_vec& vec, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	typename t_vec::value_type len = 1;
	if(!is_normalised)
		len = norm<t_vec>(vec);

	const auto sigma = su2_matrices<std::vector<t_mat>>(false);
	return inner<std::vector<t_mat>, t_vec>(sigma, vec);
}


/**
 * SU(2) ladders
 * @see https://en.wikipedia.org/wiki/Ladder_operator
 */
template<class t_mat>
const t_mat& su2_ladder(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx cI(0, 1);
	constexpr t_cplx c05(0.5, 0);

	static const t_mat mat[] =
	{
		c05*su2_matrix<t_mat>(0) + c05*cI*su2_matrix<t_mat>(1),  // up
		c05*su2_matrix<t_mat>(0) - c05*cI*su2_matrix<t_mat>(1),  // down
	};

	return mat[which];
}


/**
 * SU(3) generators, Gell-Mann matrices
 * @see https://de.wikipedia.org/wiki/Gell-Mann-Matrizen
 */
template<class t_mat>
const t_mat& su3_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx c0(0, 0);
	constexpr t_cplx c1(1, 0);
	constexpr t_cplx c2(2, 0);
	constexpr t_cplx cI(0, 1);
	constexpr t_real s3 = std::sqrt(3.);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0,c1,c0}, {c1,c0,c0}, {c0,c0,c0}}),               // 1
		create<t_mat>({{c0,cI,c0}, {-cI,c0,c0}, {c0,c0,c0}}),              // 2
		create<t_mat>({{c1,c0,c0}, {c0,-c1,c0}, {c0,c0,c0}}),              // 3
		create<t_mat>({{c0,c0,c1}, {c0,c0,c0}, {c1,c0,c0}}),               // 4
		create<t_mat>({{c0,c0,cI}, {c0,c0,c0}, {-cI,c0,c0}}),              // 5
		create<t_mat>({{c0,c0,c0}, {c0,c0,c1}, {c0,c1,c0}}),               // 6
		create<t_mat>({{c0,c0,c0}, {c0,c0,cI}, {c0,-cI,c0}}),              // 7
		create<t_mat>({{c1/s3,c0,c0}, {c0,c1/s3,c0}, {c0,c0,-c2/s3*c1}}),  // 8
	};

	return mat[which];
}


/**
 * conjugate complex vector
 */
template<class t_vec>
t_vec conj(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	const std::size_t N = vec.size();
	t_vec vecConj = zero<t_vec>(N);

	for(std::size_t comp = 0; comp < N; ++comp)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			vecConj[comp] = std::conj(vec[comp]);
		else	// simply copy non-complex vector
			vecConj[comp] = vec[comp];
	}

	return vecConj;
}


/**
 * conjugate complex matrix
 */
template<class t_mat>
t_mat conj(const t_mat& mat)
requires is_basic_mat<t_mat>
{
	const auto mat_size1 = mat.size1();
	const auto mat_size2 = mat.size2();

	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat_size2, mat_size1);

	for(std::size_t i = 0; i < mat_size1; ++i)
	{
		for(std::size_t j = 0; j < mat_size2; ++j)
		{
			if constexpr(is_complex<typename t_mat::value_type>)
				mat2(i, j) = std::conj(mat(i, j));
			else	// simply transpose non-complex matrix
				mat2(i, j) = mat(i, j);
		}
	}

	return mat2;
}


/**
 * hermitian conjugate complex matrix
 */
template<class t_mat>
t_mat herm(const t_mat& mat)
requires is_basic_mat<t_mat>
{
	const auto mat_size1 = mat.size1();
	const auto mat_size2 = mat.size2();

	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat_size2, mat_size1);

	for(std::size_t i = 0; i < mat_size1; ++i)
	{
		for(std::size_t j = 0; j < mat_size2; ++j)
		{
			if constexpr(is_complex<typename t_mat::value_type>)
				mat2(j, i) = std::conj(mat(i, j));
			else	// simply transpose non-complex matrix
				mat2(j, i) = mat(i, j);
		}
	}

	return mat2;
}

}

#endif
