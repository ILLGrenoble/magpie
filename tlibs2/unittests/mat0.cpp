/**
 * math lib and lapack test
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o mat0 mat0.cpp -llapacke
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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
 */

#define BOOST_TEST_MODULE Mat0
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_mat0, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 1./2.);
	std::cout << "eps = " << eps << std::endl;


	auto M = tl2::create<t_mat>({1, 2, 3, 3, 2, 6, 4, 2, 4});
	auto Z = tl2::create<t_mat_cplx>({1, 2, 3, 3, 2, 6, 4, 2, 4});
	std::cout << "M = " << M << std::endl;
	std::cout << "Z = " << Z << std::endl;


	{
		auto [ok, Q, R] = tl2::qr<t_mat, t_vec>(M);
		t_real d = tl2::det<t_mat>(Q);
		auto QR = Q*R;

		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "det(Q) = " << d << std::endl;
		std::cout << "QR = " << QR << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_real>(d, 1, eps));
		BOOST_TEST(tl2::equals(QR, M, eps));

#ifdef USE_LAPACK
		auto [ok2, P, L, U] = tl2_la::lu<t_mat>(M);
		auto PLU = P*L*U;

		std::cout << "\nok2 = " << std::boolalpha << ok2 << std::endl;
		std::cout << "P = " << P << std::endl;
		std::cout << "L = " << L << std::endl;
		std::cout << "U = " << U << std::endl;
		std::cout << "PLU = " << PLU << std::endl;

		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals(PLU, M, eps));
#endif
	}

	{
		auto [ok, Q, R] = tl2::qr<t_mat_cplx, t_vec_cplx>(Z);
		t_cplx d = tl2::det<t_mat_cplx>(Q);
		auto QR = Q*R;

		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "det(Q) = " << d << std::endl;
		std::cout << "QR = " << QR << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals<t_cplx>(d, 1, eps));
		BOOST_TEST(tl2::equals(QR, Z, eps));

#ifdef USE_LAPACK
		auto [ok2, P, L, U] = tl2_la::lu<t_mat_cplx>(Z);
		auto PLU = P*L*U;

		std::cout << "\nok2 = " << std::boolalpha << ok2 << std::endl;
		std::cout << "P = " << P << std::endl;
		std::cout << "L = " << L << std::endl;
		std::cout << "U = " << U << std::endl;
		std::cout << "PLU = " << PLU << std::endl;

		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals(PLU, Z, eps));
#endif
	}

#ifdef USE_LAPACK
	{
		auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(Z, 0, 0, 1);
		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		for(std::size_t i=0; i<evals.size(); ++i)
			std::cout << "eval: " << evals[i] << ", evec: " << evecs[i] << std::endl;


		auto [ok2, U, Vh, vals] = tl2_la::singval<t_mat_cplx>(Z);
		std::cout << "\nok = " << std::boolalpha << ok2 << std::endl;
		std::cout << "singvals: ";
		for(std::size_t i=0; i<vals.size(); ++i)
			std::cout << vals[i] << " ";
		std::cout << std::endl;
		std::cout << "U = " << U << "\nVh = " << Vh << std::endl;

		std::cout << "diag{vals} * UVh = " << U*tl2::diag<t_mat_cplx>(vals)*Vh << std::endl;


		auto [inva, ok3a] = tl2_la::pseudoinv<t_mat_cplx>(Z);
		auto [invb, ok3b] = tl2::inv<t_mat_cplx>(Z);
		std::cout << "\nok = " << std::boolalpha << ok3a << ", " << ok3b << std::endl;
		std::cout << "pseudoinv = " << inva << std::endl;
		std::cout << "      inv  = " << invb << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(ok3a);
		BOOST_TEST(ok3b);

		auto ident = tl2::unit<t_mat_cplx>(Z.size1(), Z.size2());
		BOOST_TEST((tl2::is_unit<t_mat_cplx>(ident, t_cplx(eps))));

		auto mata1 = inva*Z;
		auto mata2 = Z*inva;
		auto matb1 = invb*Z;
		auto matb2 = Z*invb;
		BOOST_TEST(tl2::equals(mata1, ident, eps));
		BOOST_TEST(tl2::equals(matb1, ident, eps));
		BOOST_TEST(tl2::equals(mata2, ident, eps));
		BOOST_TEST(tl2::equals(matb2, ident, eps));
	}

	{
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(M, 0, 0, 1);
		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		for(std::size_t i=0; i<evals_re.size(); ++i)
			std::cout << "eval: " << evals_re[i] << " + i*" << evals_im[i]
			<< ", evec: " << evecs_re[i] << " +i*" << evecs_im[i] << std::endl;


		auto [ok2, U, Vt, vals] = tl2_la::singval<t_mat>(M);
		std::cout << "\nok = " << std::boolalpha << ok2 << std::endl;
		std::cout << "singvals: ";
		for(std::size_t i=0; i<vals.size(); ++i)
			std::cout << vals[i] << " ";
		std::cout << std::endl;
		std::cout << "U = " << U << "\nVt = " << Vt << std::endl;

		std::cout << "diag{vals} * UVt = " << U*tl2::diag<t_mat>(vals)*Vt << std::endl;


		auto [inva, ok3a] = tl2_la::pseudoinv<t_mat>(M);
		auto [invb, ok3b] = tl2::inv<t_mat>(M);
		std::cout << "\nok = " << std::boolalpha << ok3a << ", " << ok3b << std::endl;
		std::cout << "pseudoinv = " << inva << std::endl;
		std::cout << "      inv  = " << invb << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(ok3a);
		BOOST_TEST(ok3b);

		auto ident = tl2::unit<t_mat>(M.size1(), M.size2());
		auto mata = inva*M;
		auto matb = invb*M;
		BOOST_TEST(tl2::equals(mata, ident, eps));
		BOOST_TEST(tl2::equals(matb, ident, eps));
	}

	// test rotation
	for(std::size_t iteration=0; iteration<1000; ++iteration)
	{
		t_vec axis = tl2::create<t_vec>({
			tl2::get_rand<t_real>(-10., 10.),
			tl2::get_rand<t_real>(-10., 10.),
			tl2::get_rand<t_real>(-10., 10.) });
		axis /= tl2::norm<t_vec>(axis);
		t_real angle = tl2::get_rand<t_real>(-tl2::pi<t_real>, tl2::pi<t_real>);
		t_mat rot = tl2::rotation<t_mat, t_vec>(axis, angle, 1);

		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(rot, 0, 0, 1);
		//std::cout << "\nok = " << std::boolalpha << ok << std::endl;

		bool axis_found = false;
		for(std::size_t i=0; i<evals_re.size(); ++i)
		{
			if(tl2::equals<t_real>(evals_re[i], 1., eps) && tl2::equals_0<t_real>(evals_im[i], eps))
			{
				bool evec_ok = (tl2::equals<t_vec>(axis, evecs_re[i], eps) ||
					tl2::equals<t_vec>(-axis, evecs_re[i], eps)) &&
					tl2::equals_0<t_vec>(evecs_im[i], eps);
				axis_found = true;
				BOOST_TEST((evec_ok));

				if(!evec_ok)
				{
					std::cerr << "rotation axis: " << axis << ", angle: " << angle << std::endl;
					for(std::size_t j=0; j<evals_re.size(); ++j)
					{
						std::cerr << "eval: " << evals_re[j] << " + i*" << evals_im[j]
							<< ", evec: " << evecs_re[j] << " +i*" << evecs_im[j] 
							<< std::endl;
					}
				}
			}
		}

		BOOST_TEST((axis_found));
		if(!axis_found)
		{
			std::cerr << "Error: axis not found!" << std::endl;
			std::cerr << "rotation axis: " << axis << ", angle: " << angle << std::endl;
			for(std::size_t i=0; i<evals_re.size(); ++i)
			{
				std::cerr << "eval: " << evals_re[i] << " + i*" << evals_im[i]
					<< ", evec: " << evecs_re[i] << " +i*" << evecs_im[i] 
					<< std::endl;
			}
		}
	}
#endif

	std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
}
