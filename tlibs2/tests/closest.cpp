/**
 * closest point test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 3-aug-2025
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is from my following project:
 *   - "mathlibs" (https://github.com/t-weber/mathlibs).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * "mathlibs" project
 * Copyright (C) 2021-2025  Tobias WEBER (privately developed).
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
 * clang++ -std=c++20 -I.. -Wall -Wextra -Weffc++ -o closest closest.cpp -lgmp
 */

#include <boost/geometry/index/rtree.hpp>
namespace geo = boost::geometry;
namespace geoidx = geo::index;

#include "libs/maths.h"
using namespace tl2_ops;


int main()
{
	using t_real = double;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;

	using t_vertex = geo::model::point<t_real, 3, geo::cs::cartesian>;
	using t_rtree_leaf = std::tuple<t_vertex, std::size_t>;
	using t_rtree = geoidx::rtree<t_rtree_leaf, geoidx::dynamic_linear>;

	constexpr t_real eps = 1e-6;
	constexpr t_real pi = tl2::pi<t_real>;

	t_mat B = tl2::B_matrix<t_mat>(
		4., 5., 6.,
		90./180.*pi, 90./180.*pi, 60./180.*pi);
	tl2::set_eps_0(B, eps);
	std::cout << "B = " << B << ".\n" << std::endl;

	// (hkl) reference points
	std::vector<t_vec> points;
	points.reserve(8 * 8 * 8);
	for(t_real h = -4.; h < 4.; h += 1.)
	for(t_real k = -4.; k < 4.; k += 1.)
	for(t_real l = -4.; l < 4.; l += 1.)
		points.emplace_back(tl2::create<t_vec>({ h, k, l }));

	// rtree with B
	t_rtree rt = tl2::make_rtree<
		t_real, t_vec, t_mat, 3, t_vertex, t_rtree_leaf, t_rtree>(
			points, &B);

	while(true)
	{
		// query point
		t_vec query_pt = tl2::create<t_vec>({ 0., 0., 0. });
		std::cout << "Enter (hkl) point to query: ";
		std::cin >> query_pt[0] >> query_pt[1] >> query_pt[2];

		std::vector<std::pair<std::size_t, t_real>> pt_indices =
			tl2::closest_point_rtree<
				t_real, t_vec, t_mat, 3, 16, t_vertex, t_rtree_leaf, t_rtree, std::vector>(
					rt, query_pt, &B, eps);

		// results
		std::cout << "Closest point(s) to " << query_pt << ":\n";
		for(auto [pt_idx, dist] : pt_indices)
			std::cout << "\t" << points[pt_idx] << " (distance: " << dist << " 1/A)\n";
		std::cout << std::endl;
	}

	return 0;
}
