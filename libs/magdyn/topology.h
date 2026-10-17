/**
 * magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (McClarty 2022) https://doi.org/10.1146/annurev-conmatphys-031620-104715
 *   - (Nakamura 2024) https://doi.org/10.1103/PhysRevB.110.075144
 *
 * @note Forked on 5-November-2024 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 * @desc For references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAGDYN_TOPO_H__
#define __MAGDYN_TOPO_H__


#include "tlibs2/libs/algos.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// topological calculations
// --------------------------------------------------------------------

namespace magdyn {

/**
 * calculates the berry connections
 * @see equ. 7 in (McClarty 2022)
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_vec, class t_vec_real, class t_mat, class t_S,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type,
	typename t_size = std::size_t>
std::tuple<std::vector<t_vec>, t_S> berry_connections(
	const std::function<std::tuple<t_mat, t_S>(const t_vec_real& Q, bool)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::sqrt(std::numeric_limits<t_real>::epsilon()),
	bool commutator = true, const t_mat* evecs_ref = nullptr)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
{
	constexpr const t_cplx imag{0, 1};

	// get the eigenvector phases to the same as the reference eigenvectors' phases
	// see (McClarty 2022)
	auto align_phases = [evecs_ref, &imag](t_mat& evecs)
	{
		if(!evecs_ref)
			return;

		// iterate bands
		for(t_size col = 0; col < evecs.size2(); ++col)
		{
			t_vec evec_ref = tl2::col<t_mat, t_vec>(*evecs_ref, col);
			t_vec evec = tl2::col<t_mat, t_vec>(evecs, col);

			t_cplx ph = std::exp(-imag * std::arg(tl2::inner(evec_ref, evec)));

			for(t_size row = 0; row < evecs.size1(); ++row)
				evecs(row, col) *= ph;
		}
	};

	auto [ evecs, S ] = get_evecs(Q, commutator);
	align_phases(evecs);

	const t_size BANDS = evecs.size2();
	const t_size DIM = Q.size();

	std::vector<t_vec> connections{};
	connections.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
		connections.emplace_back(tl2::create<t_vec>(DIM));

	for(t_size dim = 0; dim < DIM; ++dim)
	{
		t_vec_real Q1 = Q;
		Q1[dim] += delta;

		t_mat evecs_Q1 = std::get<0>(get_evecs(Q1, commutator));
		align_phases(evecs_Q1);

		// differentiate eigenvector matrix
		t_mat evecs_diff = (evecs - evecs_Q1) / delta;

		t_mat evecs_H = tl2::herm(evecs);
		t_mat C;
		if(commutator)
			C = mult_comm<t_mat>(evecs_H) * mult_comm<t_mat>(evecs_diff);
		else
			C = evecs_H * evecs_diff;

		for(t_size band = 0; band < BANDS; ++band)
			connections[band][dim] = C(band, band) * imag;
	}

	return std::make_tuple(connections, S);
}



/**
 * see equ. 2 in (Nakamura 2024)
 */
template<class t_vec, class t_cplx = typename t_vec::value_type,
	class t_real = typename t_cplx::value_type>
t_cplx berry_dot(const t_vec& v1, const t_vec& v2, bool commutator, t_real eps)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec>
#endif
{
	t_cplx s = commutator
		? tl2::inner(v1, mult_comm<t_vec>(v2))
		: tl2::inner(v1, v2);

	t_real len = std::abs(s);
	if(len > eps)
		s /= len;

	return s;
}



/**
 * calculates the berry curvatures
 * @see equ. 8 in (McClarty 2022)
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 *
 * t_mat: either matrix or vector container
 */
template<class t_vec, class t_vec_real, class t_mat, class t_S,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type,
	typename t_size = std::size_t, int calc_method = 1>
std::tuple<std::vector<t_cplx>, t_S> berry_curvatures(
	const std::function<std::tuple<t_mat, t_S>(const t_vec_real& Q, bool)>& get_evecs,
	const t_vec_real& Q,
	t_real delta_diff = 1e-6, // std::sqrt(std::numeric_limits<t_real>::epsilon()),
	t_real delta_int = 1e-2,
	t_size dim1 = 0, t_size dim2 = 1,
	bool commutator = true,
	t_real eps = std::numeric_limits<t_real>::epsilon())
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
{
	// only valid in three dimensions
	assert(Q.size() == 3);

	const auto [ evecs, S ] = get_evecs(Q, commutator);
	t_size BANDS = evecs.size2();

	std::vector<t_cplx> curvatures{};
	curvatures.reserve(BANDS);

	// direct implementation of equ. 8 in (McClarty 2022)
	if constexpr(calc_method == 0)
	{
		t_vec_real h = Q, k = Q;
		h[dim1] += delta_diff;
		k[dim2] += delta_diff;

		std::vector<t_vec> connections =
			std::get<0>(berry_connections<t_vec, t_vec_real, t_mat, t_S, t_cplx, t_real, t_size>(
				get_evecs, Q, delta_diff, commutator, &evecs));
		std::vector<t_vec> connections_h =
			std::get<0>(berry_connections<t_vec, t_vec_real, t_mat, t_S, t_cplx, t_real, t_size>(
				get_evecs, h, delta_diff, commutator, &evecs));
		std::vector<t_vec> connections_k =
			std::get<0>(berry_connections<t_vec, t_vec_real, t_mat, t_S, t_cplx, t_real, t_size>(
				get_evecs, k, delta_diff, commutator, &evecs));

		for(t_size band = 0; band < BANDS; ++band)
		{
			// differentiate connection's y component with respect to h
			t_cplx curv1 = (connections_h[band][dim2] - connections[band][dim2]) / delta_diff;
			// differentiate connection's x component with respect to k
			t_cplx curv2 = (connections_k[band][dim1] - connections[band][dim1]) / delta_diff;

			curvatures.emplace_back(curv1 - curv2);
		}
	}

	// FHS method, see (Nakamura 2024)
	else
	{
		t_vec_real Q01 = Q;  Q01[dim2] += delta_int;
		t_vec_real Q10 = Q;  Q10[dim1] += delta_int;
		t_vec_real Q11 = Q;  Q11[dim1] += delta_int;  Q11[dim2] += delta_int;

		const t_mat& evecs00 = evecs;
		const auto [ evecs01, S01 ] = get_evecs(Q01, commutator);
		const auto [ evecs10, S10 ] = get_evecs(Q10, commutator);
		const auto [ evecs11, S11 ] = get_evecs(Q11, commutator);

		// see equs. 3a & 3b in (Nakamura 2024)
		for(t_size band = 0; band < BANDS; ++band)
		{
				t_vec v00 = tl2::col<t_mat, t_vec>(evecs00, band);
				t_vec v01 = tl2::col<t_mat, t_vec>(evecs01, band);
				t_vec v10 = tl2::col<t_mat, t_vec>(evecs10, band);
				t_vec v11 = tl2::col<t_mat, t_vec>(evecs11, band);

				t_cplx curv = std::arg(
					berry_dot(v00, v10, commutator, eps) * berry_dot(v01, v00, commutator, eps) *
					berry_dot(v11, v01, commutator, eps) * berry_dot(v10, v11, commutator, eps));
				curv /= delta_int*delta_int;

				curvatures.emplace_back(std::move(curv));
		}
	}

	return std::make_tuple(curvatures, S);
}



/**
 * calculates the chern numbers
 * @see equ. 9 in (McClarty 2022)
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_vec, class t_vec_real, class t_mat, class t_S,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type,
	typename t_size = std::size_t, int calc_method = 2>
std::vector<t_cplx> chern_numbers(
	const std::function<std::tuple<t_mat, t_S>(const t_vec_real& Q, bool)>& get_evecs,
	t_real bz = 0.5,  // brillouin zone boundary
	t_real eps = std::numeric_limits<t_real>::epsilon(),
	t_real delta_diff = 1e-6, //std::sqrt(std::numeric_limits<t_real>::epsilon()),
	t_real delta_int = 1e-2, //std::cbrt(std::numeric_limits<t_real>::epsilon()),
	t_size dim1 = 0, t_size dim2 = 1, bool commutator = true)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
{
	std::vector<t_cplx> chern_nums;

	// calculate via area over berry curvatures
	// direct implementation of equ. 9 in (McClarty 2022)
	if constexpr(calc_method == 0)
	{
		t_vec_real Q = tl2::zero<t_vec_real>(3);

		for(Q[dim1] = -bz; Q[dim1] < bz; Q[dim1] += delta_int)
		for(Q[dim2] = -bz; Q[dim2] < bz; Q[dim2] += delta_int)
		{
			std::vector<t_cplx> curvs = std::get<0>(berry_curvatures<
				t_vec, t_vec_real, t_mat, t_S, t_cplx, t_real, t_size>(
					get_evecs, Q, delta_diff, delta_int, dim1, dim2, commutator));

			// initialise by resetting chern numbers to zeros
			if(!chern_nums.size())
				chern_nums.resize(curvs.size(), t_cplx{});

			// numerically integrate the brillouin zone area for each band
			for(t_size band = 0; band < curvs.size(); ++band)
				chern_nums[band] += curvs[band] * delta_int * delta_int;
		}

		for(t_cplx& num : chern_nums)
			num /= t_real(2) * tl2::pi<t_real>;
	}

	// calculate via boundary over berry connections
	// direct implementation of equ. 9 in (McClarty 2022)
	else if constexpr(calc_method == 1)
	{
		auto int_boundary = [get_evecs, delta_diff, delta_int, bz, dim1, dim2,
			&chern_nums, commutator](t_size dim, t_vec_real& Q, t_real sign)
		{
			for(Q[dim] = -bz; Q[dim] < bz; Q[dim] += delta_int)
			{
				std::vector<t_vec> conns =
					std::get<0>(berry_connections<t_vec, t_vec_real, t_mat, t_S, t_cplx, t_real, t_size>(
						get_evecs, Q, delta_diff, commutator));
	
				// initialise by resetting chern numbers to zeros
				if(!chern_nums.size())
					chern_nums.resize(conns.size(), t_cplx{});
	
				// numerically integrate along boundary segment
				for(t_size band = 0; band < conns.size(); ++band)
				{
					chern_nums[band] += conns[band][dim1] * delta_int * sign;
					chern_nums[band] -= conns[band][dim2] * delta_int * sign;
				}
			}
		};

		// bottom part of boundary
		t_vec_real Q = tl2::zero<t_vec_real>(3);
		Q[dim2] = -bz;
		int_boundary(dim1, Q, 1.);

		// top part of boundary
		Q[dim2] = +bz;
		int_boundary(dim1, Q, -1.);

		// left part of boundary
		Q = tl2::zero<t_vec_real>(3);
		Q[dim1] = -bz;
		int_boundary(dim2, Q, 1.);

		// right part of boundary
		Q[dim1] = +bz;
		int_boundary(dim2, Q, -1.);
	}

	// FHS method, see (Nakamura 2024)
	else
	{		
		// create grid
		const t_size N = t_size(2.*bz / delta_int);
		std::vector<t_mat> grid(N*N);

		t_vec_real Q = tl2::zero<t_vec_real>(3);
		Q[dim1] = -bz;

		for(t_size i = 0; i < N; ++i)
		{
			Q[dim2] = -bz;

			for(t_size j = 0; j < N; ++j)
			{
				auto [ evecs, S ] = get_evecs(Q, commutator);
				grid[i*N + j] = std::move(evecs);

				Q[dim2] += delta_int;
			}

			Q[dim1] += delta_int;
		}

		// initialise by resetting chern numbers to zeros
		const t_size num_bands = grid[0].size2();
		if(!chern_nums.size())
			chern_nums.resize(num_bands, t_cplx{});

		// calculate chern numbers per magnon band
		for(t_size band = 0; band < num_bands; ++band)
		{
			t_cplx chern_num{};

			// see equs. 3a & 3b in (Nakamura 2024)
			for(t_size i = 0; i < N; ++i)
			for(t_size j = 0; j < N; ++j)
			{
				t_vec v00 = tl2::col<t_mat, t_vec>(grid[i*N           + j          ], band);
				t_vec v01 = tl2::col<t_mat, t_vec>(grid[i*N           + ((j + 1)%N)], band);
				t_vec v10 = tl2::col<t_mat, t_vec>(grid[((i + 1)%N)*N + j          ], band);
				t_vec v11 = tl2::col<t_mat, t_vec>(grid[((i + 1)%N)*N + ((j + 1)%N)], band);

				chern_num += std::arg(
					berry_dot(v00, v10, commutator, eps) * berry_dot(v01, v00, commutator, eps) *
					berry_dot(v11, v01, commutator, eps) * berry_dot(v10, v11, commutator, eps));
			}

			chern_num /= t_real(2) * tl2::pi<t_real>;
			chern_nums[band] = std::move(chern_num);
		}
	}

	// should be integers
	return chern_nums;
}



/**
 * get a function returning the matrix of eigenstates at specific Q
 */
template<class t_mat, class t_vec, class t_vec_real, class t_magdyn, class t_perm>
static std::function<std::tuple<t_mat, typename t_magdyn::SofQE>(const t_vec_real& Q, bool)>
get_evecmat_func(const t_magdyn *magdyn, const t_perm *perm = nullptr)
{
	return [magdyn, perm](const t_vec_real& Q, bool commutator)
		-> std::tuple<t_mat, typename t_magdyn::SofQE>
	{
		typename t_magdyn::SofQE S = magdyn->CalcEnergies(Q, false);

		const t_mat& evec_mat = commutator ? S.evec_mat_comm : S.evec_mat;
		//t_mat M = evec_mat;  // incommensurate parts not considered in evec_mat

		t_mat M = tl2::create<t_mat>(evec_mat.size1(), S.E_and_S.size());
		for(std::size_t E_idx = 0; E_idx < S.E_and_S.size(); ++E_idx)
		{
			const auto& E_and_S = S.E_and_S[E_idx];
			set_col<t_mat, t_vec>(M,
				commutator ? E_and_S.state_comm : E_and_S.state, E_idx);
		}

		if(perm)
		{
			M = tl2::reorder_cols<t_mat, t_vec>(M, *perm);
			S.E_and_S = tl2::reorder(S.E_and_S, *perm);
		}

		return std::make_tuple(M, S);
	};
}
}   // namespace magdyn



/**
 * get the berry connection for each magnon band
 */
MAGDYN_TEMPL
std::tuple<std::vector<t_vec>, MAGDYN_TYPE::SofQE>
MAGDYN_INST::CalcBerryConnections(const t_vec3_real& Q, t_real delta,
	const std::vector<t_size>* perm, bool commutator) const
{
	//SetUniteDegenerateEnergies(false);
	if(m_unite_degenerate_energies)
	{
		MAGDYN_CERR_OPT << "Magdyn warning: "
			<< "Merging of degeneracies should be disabled for topological calculations."
			<< std::endl;
	}

	auto evec_func = get_evecmat_func<t_mat, t_vec, t_vec3_real>(this, perm);

	return berry_connections<t_vec, t_vec3_real, t_mat, SofQE, t_cplx, t_real, t_size>(
		evec_func, Q, delta, commutator);
}



/**
 * get the berry curvature for each magnon band
 */
MAGDYN_TEMPL
std::tuple<std::vector<t_cplx>, MAGDYN_TYPE::SofQE> MAGDYN_INST::CalcBerryCurvatures(
	const t_vec3_real& Q, t_real delta_diff, t_real delta_int,
	const std::vector<t_size>* perm, t_size dim1, t_size dim2,
	bool commutator) const
{
	//SetUniteDegenerateEnergies(false);
	if(m_unite_degenerate_energies)
	{
		MAGDYN_CERR_OPT << "Magdyn warning: "
			<< "Merging of degeneracies should be disabled for topological calculations."
			<< std::endl;
	}

	auto evec_func = get_evecmat_func<t_mat, t_vec, t_vec3_real>(this, perm);
	constexpr const int calc_method = 1;

	return berry_curvatures<
		t_vec, t_vec3_real, t_mat, SofQE, t_cplx, t_real, t_size, calc_method>(
			evec_func, Q, delta_diff, delta_int, dim1, dim2, commutator, m_eps);
}



/**
 * get the chern number for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_cplx> MAGDYN_INST::CalcChernNumbers(
	t_real bz, t_real delta_diff, t_real delta_int, t_size dim1, t_size dim2,
	bool commutator) const
{
	//SetUniteDegenerateEnergies(false);
	if(m_unite_degenerate_energies)
	{
		MAGDYN_CERR_OPT << "Magdyn warning: "
			<< "Merging of degeneracies should be disabled for topological calculations."
			<< std::endl;
	}


	constexpr const int calc_method = 2;
	std::vector<t_size> *perm = nullptr;

	auto evec_func = get_evecmat_func<t_mat, t_vec, t_vec3_real>(this, perm);

	return chern_numbers<t_vec, t_vec3_real, t_mat, SofQE,
		t_cplx, t_real, t_size, calc_method>(
			evec_func, bz, m_eps, delta_diff, delta_int, dim1, dim2,
			commutator);
}

// --------------------------------------------------------------------

#endif
