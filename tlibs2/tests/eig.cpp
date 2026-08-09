/**
 * benchmark eigensystem calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 9-aug-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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
 *
 * reference: clang++ -std=c++20 -O3 -march=native -DUSE_LAPACK -I.. -I/usr/local/include -I/opt/homebrew/include/ -L/usr/local/lib -L/opt/homebrew/Cellar/gcc/16.1.0/lib/gcc/16/ -o eig eig.cpp -llapacke -llapack -lblas -lgfortran
 * openblas : clang++ -std=c++20 -O3 -march=native -DUSE_LAPACK -I.. -I/usr/local/include/openblas -I/opt/homebrew/include/ -L/usr/local/lib -L/opt/homebrew/Cellar/gcc/16.1.0/lib/gcc/16/ -o eig eig.cpp -lopenblas -lgfortran 
 */

#define BOOST_TEST_MODULE Eigenvector Benchmark
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>
#include <random>

#include "libs/algos.h"
#include "libs/maths.h"

//#define DO_CHECKS


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_eig, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	t_real eps = 1e-7;
	if constexpr(std::is_same_v<t_real, float>)
		eps = 1e-1;

	std::cout.precision(8);
	const int w = 16;
	std::size_t maxdim = 1024;

	std::mt19937 rndgen{/*tl2::epoch<unsigned int>()*/ 1234};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	// real version
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "real, " << sizeof(t_real)*8 << "-bit" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(w) << std::left << "size" << " "
		<< std::setw(w) << std::left << "time (s)" << std::endl;
	for(std::size_t dim = 2; dim <= maxdim; dim *= 2)
	{
		tl2::Stopwatch<t_real> stopwatch;
		stopwatch.start();

		auto mat = tl2::zero<t_mat>(dim,dim);
		for(std::size_t i = 0; i < dim; ++i)
			for(std::size_t j = 0; j < dim; ++j)
				mat(i, j) = rnddist(rndgen);

		bool sym = false;
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(
				mat, false, sym, true);

		stopwatch.stop();
		t_real dur = stopwatch.GetDur();
		std::cout << std::setw(w) << std::left << dim << " "
			<< std::setw(w) << std::left << dur << std::endl;

		BOOST_TEST(ok);

#ifdef DO_CHECKS
		// checks
		t_mat_cplx mat_cplx = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i = 0; i < dim; ++i)
			for(std::size_t j = 0; j < dim; ++j)
				mat_cplx(i, j) = mat(i, j);

		for(std::size_t i = 0; i < dim; ++i)
		{
			t_cplx eval = t_cplx{evals_re[i], evals_im[i]};
			t_vec_cplx evec = tl2::zero<t_vec_cplx>(dim);

			for(std::size_t j = 0; j < dim; ++j)
				evec[j] = t_cplx{evecs_re[i][j], evecs_im[i][j]};

			t_vec_cplx tstvec1 = mat_cplx * evec;
			t_vec_cplx tstvec2 = eval * evec;
			bool is_equal = tl2::equals<t_vec_cplx, t_cplx>(tstvec1, tstvec2, eps);
			BOOST_TEST(is_equal);
		}
#endif
	}
	std::cout << "--------------------------------------------------------------------------------\n" << std::endl;


	// real version, symmetric
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "real, symmetric, " << sizeof(t_real)*8 << "-bit" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(w) << std::left << "size" << " "
		<< std::setw(w) << std::left << "time, sym. (s)" << " "
		<< std::setw(w) << std::left << "time, gen. (s)" << std::endl;
	for(std::size_t dim = 2; dim <= maxdim; dim*=2)
	{
		tl2::Stopwatch<t_real> stopwatch;
		stopwatch.start();

		auto mat = tl2::zero<t_mat>(dim,dim);
		for(std::size_t i = 0; i < dim; ++i)
			for(std::size_t j = i; j < dim; ++j)
				mat(j, i) = mat(i, j) = rnddist(rndgen);

		bool sym = true;
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(
				mat, false, sym, true);

		stopwatch.stop();
		t_real dur1 = stopwatch.GetDur();

		// compare results with non-symmetric calculation
		sym = false;
		auto [ok_gen, evals_re_gen, evals_im_gen, evecs_re_gen, evecs_im_gen] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(
				mat, false, sym, true);

		stopwatch.stop();
		t_real dur2 = stopwatch.GetDur();
		std::cout << std::setw(w) << std::left << dim << " "
			<< std::setw(w) << std::left << dur1 << " " << dur2 << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok_gen);

#ifdef DO_CHECKS
		t_mat_cplx mat_cplx = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i = 0; i < dim; ++i)
			for(std::size_t j = 0; j < dim; ++j)
				mat_cplx(i, j) = mat(i, j);

		for(std::size_t i = 0; i < dim; ++i)
		{
			t_cplx eval = t_cplx{evals_re[i], evals_im[i]};
			t_vec_cplx evec = tl2::zero<t_vec_cplx>(dim);
			for(std::size_t j = 0; j < dim; ++j)
				evec[j] = t_cplx{evecs_re[i][j], evecs_im[i][j]};

			t_cplx eval_gen = t_cplx{evals_re_gen[i], evals_im_gen[i]};
			t_vec_cplx evec_gen = tl2::zero<t_vec_cplx>(dim);
			for(std::size_t j = 0; j < dim; ++j)
				evec_gen[j] = t_cplx{evecs_re_gen[i][j], evecs_im_gen[i][j]};

			t_vec_cplx tstvec1 = mat_cplx * evec;
			t_vec_cplx tstvec2 = eval * evec;

			t_vec_cplx tstvec1_gen = mat_cplx * evec_gen;
			t_vec_cplx tstvec2_gen = eval_gen * evec_gen;

			bool is_equal = tl2::equals<t_vec_cplx, t_cplx>(tstvec1, tstvec2, eps);
			bool is_equal_gen = tl2::equals<t_vec_cplx, t_cplx>(tstvec1_gen, tstvec2_gen, eps);

			BOOST_TEST(is_equal);
			BOOST_TEST(is_equal_gen);
		}
#endif
	}
	std::cout << "--------------------------------------------------------------------------------\n" << std::endl;


	// complex version
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "complex, " << sizeof(t_real)*8 << "-bit" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(w) << std::left << "size" << " "
		<< std::setw(w) << std::left << "time (s)" << std::endl;
	for(std::size_t dim = 2; dim <= maxdim; dim*=2)
	{
		tl2::Stopwatch<t_real> stopwatch;
		stopwatch.start();

		auto mat = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i = 0; i < dim; ++i)
			for(std::size_t j = 0; j < dim; ++j)
				mat(i, j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};

		bool herm = false;
		auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(
				mat, false, herm, true);

		stopwatch.stop();
		t_real dur = stopwatch.GetDur();
		std::cout << std::setw(w) << std::left << dim << " "
			<< std::setw(w) << std::left << dur << std::endl;

		BOOST_TEST(ok);

#ifdef DO_CHECKS
		for(std::size_t i = 0; i < dim; ++i)
		{
			t_vec_cplx tstvec1 = mat*evecs[i];
			t_vec_cplx tstvec2 = evals[i]*evecs[i];
			bool is_equal = tl2::equals<t_vec_cplx, t_cplx>(tstvec1, tstvec2, eps);
			BOOST_TEST(is_equal);
		}
#endif
	}
	std::cout << "--------------------------------------------------------------------------------\n" << std::endl;


	// complex version, hermitian
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "complex, hermitian, " << sizeof(t_real)*8 << "-bit" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(w) << std::left << "size" << " "
		<< std::setw(w) << std::left << "time, herm. (s)" << " "
		<< std::setw(w) << std::left << "time, gen. (s)" << std::endl;
	for(std::size_t dim = 2; dim <= maxdim; dim*=2)
	{
		tl2::Stopwatch<t_real> stopwatch;
		stopwatch.start();

		auto mat = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i = 0; i < dim; ++i)
		{
			for(std::size_t j = i + 1; j < dim; ++j)
			{
				mat(i, j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};
				mat(j, i) = std::conj(mat(i, j));
			}
		}
		for(std::size_t i = 0; i < dim; ++i)
			mat(i, i) = rnddist(rndgen);

		bool norm = true;
		bool herm = true;
		auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(
				mat, false, herm, norm);

		stopwatch.stop();
		t_real dur1 = stopwatch.GetDur();

		// comparison with non-hermitian calulation
		stopwatch.start();
		herm = false;
		auto [ok_gen, evals_gen, evecs_gen] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(
				mat, false, herm, norm);

		stopwatch.stop();
		t_real dur2 = stopwatch.GetDur();
		std::cout << std::setw(w) << std::left << dim << " "
			<< std::setw(w) << std::left << dur1 << " " << dur2 << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok_gen);

#ifdef DO_CHECKS
		// test diagonalisation
		t_mat_cplx eval_mat = tl2::diag<t_mat_cplx>(evals);
		t_mat_cplx evec_mat = tl2::create<t_mat_cplx>(evecs);
		t_mat_cplx diag = tl2::herm(evec_mat) * mat * evec_mat;
		bool diag_equ = tl2::equals<t_mat_cplx, t_cplx>(eval_mat, diag, eps);
		BOOST_TEST(diag_equ);
		tl2::set_eps_0<t_mat_cplx>(eval_mat, eps);
		tl2::set_eps_0<t_mat_cplx>(diag, eps);

		t_mat_cplx eval_mat_gen = tl2::diag<t_mat_cplx>(evals_gen);
		t_mat_cplx evec_mat_gen = tl2::create<t_mat_cplx>(evecs_gen);
		t_mat_cplx diag_gen = tl2::herm(evec_mat_gen) * mat * evec_mat_gen;
		bool diag_equ_gen = tl2::equals<t_mat_cplx, t_cplx>(eval_mat_gen, diag_gen, eps);
		BOOST_TEST(diag_equ_gen);
		tl2::set_eps_0<t_mat_cplx>(eval_mat_gen, eps);
		tl2::set_eps_0<t_mat_cplx>(diag_gen, eps);

		for(std::size_t i = 0; i < dim; ++i)
		{
			t_vec_cplx tstvec1 = mat*evecs[i];
			t_vec_cplx tstvec2 = evals[i]*evecs[i];

			t_vec_cplx tstvec1_gen = mat*evecs_gen[i];
			t_vec_cplx tstvec2_gen = evals_gen[i]*evecs_gen[i];

			bool is_equal = tl2::equals<t_vec_cplx, t_cplx>(tstvec1, tstvec2, eps);
			bool is_equal_gen = tl2::equals<t_vec_cplx, t_cplx>(tstvec1_gen, tstvec2_gen, eps);

			BOOST_TEST(is_equal);
			BOOST_TEST(is_equal_gen);
		}
#endif
	}
	std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
}
