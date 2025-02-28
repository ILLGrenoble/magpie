/**
 * calculates Q depending on detector shape
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-2021
 * @license GPLv3, see 'LICENSE' file
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
 *
 * g++-10 -std=c++20 -I.. -o detshape detshape.cpp
 */

#include <vector>

#include "libs/maths.h"


using t_real = double;
using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;


t_vec get_spherical_pos(t_real theta, t_real phi)
{
	theta = -theta;
	theta += tl2::pi<t_real> * t_real(0.5);
	phi += tl2::pi<t_real> * t_real(0.5);

	t_real x = std::cos(phi)*std::sin(theta);
	t_real y = std::sin(phi)*std::sin(theta);
	t_real z = std::cos(theta);

	t_vec vec = tl2::create<t_vec>({x, y, z});
	//t_real len = tl2::norm<t_vec>(vec);
	//std::cout << len << std::endl;

	return vec;
}


t_vec get_cylindrical_pos(t_real z, t_real phi, t_real rad)
{
	phi += tl2::pi<t_real> * t_real(0.5);

	t_real x = rad*std::cos(phi);
	t_real y = rad*std::sin(phi);

	t_vec vec = tl2::create<t_vec>({x, y, z});
	t_real len = tl2::norm<t_vec>(vec);
	//std::cout << len << std::endl;
	vec /= len;

	return vec;
}


t_vec get_flat_pos(t_real z, t_real phi, t_real rad)
{
	phi = -phi;

	t_real x = rad*std::tan(phi);
	t_real y = rad;
	z = z;

	t_vec vec = tl2::create<t_vec>({x, y, z});
	t_real len = tl2::norm<t_vec>(vec);
	//std::cout << len << std::endl;
	vec /= len;

	return vec;
}


int main()
{
	using namespace tl2_ops;

	// y is along beam, z is up
	t_real k = 1.;
	t_vec ki = tl2::create<t_vec>({0, 1, 0}) * k;

	// select a pixel on the detector
	t_real det_dist = 5.;
	t_real z = 2.5;
	t_real phi = t_real(45) / t_real(180) * tl2::pi<t_real>;

	// get the Q vector corresponding to this pixel
	{
		//t_real theta = t_real(10) / t_real(180) * tl2::pi<t_real>;
		t_real theta = std::atan2(z, det_dist);
		t_vec kf = get_spherical_pos(theta, phi) * k;

		t_vec Q = ki - kf;
		std::cout << "spherical detector\n";
		std::cout << "ki = " << ki << "\n";
		std::cout << "kf = " << kf << "\n";
		std::cout << "Q = " << Q << std::endl;
	}

	{
		t_vec kf = get_cylindrical_pos(z, phi, det_dist) * k;

		t_vec Q = ki - kf;
		std::cout << "\ncylindrical detector\n";
		std::cout << "ki = " << ki << "\n";
		std::cout << "kf = " << kf << "\n";
		std::cout << "Q = " << Q << std::endl;
	}

	{
		t_vec kf = get_flat_pos(z, phi, det_dist) * k;

		t_vec Q = ki - kf;
		std::cout << "\nflat detector\n";
		std::cout << "ki = " << ki << "\n";
		std::cout << "kf = " << kf << "\n";
		std::cout << "Q = " << Q << std::endl;
	}

	return 0;
}
