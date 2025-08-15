/**
 * tlibs2 maths library -- r-tree algorithms
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * mathlibs" and "geo" projects
 * Copyright (C) 2020-2025  Tobias WEBER (privately developed).
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

#ifndef __TL2_RTREE_H__
#define __TL2_RTREE_H__


#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <vector>
#include <algorithm>

#include "ndim.h"


namespace tl2 {

/**
 * converts a vector to a boost::geometry vertex
 */
template<class t_vertex, class t_vec, std::size_t ...indices>
constexpr void to_geo_vertex(t_vertex& vert,
  const t_vec& vec, const std::index_sequence<indices...>&)
{
	(boost::geometry::set<indices>(vert, vec[indices]), ...);
}


/**
 * creates a vector from a boost::geometry vertex
 */
template<class t_vec, class t_vertex, std::size_t ...indices>
constexpr t_vec from_geo_vertex(const t_vertex& vert,
	const std::index_sequence<indices...>&)
{
	return tl2::create<t_vec>({ boost::geometry::get<indices>(vert) ... });
}


/**
 * creates an r-tree out of a collection of points
 * the points are transformed using B if given
 */
template<class t_real = double,
	class t_vec = tl2::vec<t_real, std::vector>,
	class t_mat = tl2::mat<t_real, std::vector>,
	std::size_t dim = 3,
	class t_vertex = boost::geometry::model::point<t_real, dim, boost::geometry::cs::cartesian>,
	class t_rtree_leaf = std::tuple<t_vertex, std::size_t>,
	class t_rtree = boost::geometry::index::rtree<t_rtree_leaf, boost::geometry::index::dynamic_linear>,
	template<class...> class t_cont = std::vector>
requires tl2::is_vec<t_vec> && tl2::is_mat<t_mat>
t_rtree make_rtree(const t_cont<t_vec>& points, const t_mat* B = nullptr)
{
	using namespace tl2_ops;

	t_rtree rt(typename t_rtree::parameters_type(points.size()));

	for(std::size_t ptidx = 0; ptidx < points.size(); ++ptidx)
	{
		t_vertex vert;
		tl2::to_geo_vertex<t_vertex, t_vec>(vert,
			B ? (*B) * points[ptidx] : points[ptidx],
			std::make_index_sequence<dim>());

		rt.insert(std::make_tuple(vert, ptidx));
	}

	return rt;
}


/**
 * finds the nearest neighbours in an r-tree
 * the query point is transformed using B if given
 * returns indices into the array used to create the rtree and distances to the query point
 */
template<class t_real = double,
	class t_vec = tl2::vec<t_real, std::vector>,
	class t_mat = tl2::mat<t_real, std::vector>,
	std::size_t dim = 3, std::size_t max_num = 16,
	class t_vertex = boost::geometry::model::point<t_real, dim, boost::geometry::cs::cartesian>,
	class t_rtree_leaf = std::tuple<t_vertex, std::size_t>,
	class t_rtree = boost::geometry::index::rtree<t_rtree_leaf, boost::geometry::index::dynamic_linear>,
	template<class...> class t_cont = std::vector>
requires tl2::is_vec<t_vec> && tl2::is_mat<t_mat>
t_cont<std::pair<std::size_t, t_real>> closest_point_rtree(
	const t_rtree& rt, const t_vec& query_point,
	const t_mat* B = nullptr, t_real eps = 1e-6)
{
	using namespace tl2_ops;
	
	using t_elem = std::pair<std::size_t, t_real>;
	t_cont<t_rtree_leaf> query_answers;
	t_cont<t_elem> found_points;
	query_answers.reserve(max_num);
	found_points.reserve(max_num);

	t_vertex query_pt;
	to_geo_vertex<t_vertex, t_vec>(query_pt,
		B ? (*B) * query_point : query_point,
		std::make_index_sequence<dim>());

	rt.query(boost::geometry::index::nearest(query_pt, max_num),
		std::back_inserter(query_answers));
	for(const t_rtree_leaf& ans : query_answers)
	{
		const t_vertex& pt = std::get<0>(ans);
		std::size_t idx = std::get<1>(ans);
		t_real dist = tl2::norm<t_vec>(
		    (B ? (*B) * query_point : query_point) -
		    from_geo_vertex<t_vec, t_vertex>(pt, std::make_index_sequence<dim>()));

		found_points.emplace_back(std::make_pair(idx, dist));
	}

	if(found_points.size() == 0)
		return found_points;

	// sort points by distance
	std::stable_sort(found_points.begin(), found_points.end(),
		[](const t_elem& elem1, const t_elem& elem2) -> bool
	{
			return std::get<1>(elem1) < std::get<1>(elem2);
	});

	// remove all points that are furhter away than the nearest neighbour
	t_real closest = std::get<1>(found_points[0]);
	for(auto iter = found_points.begin() + 1; iter != found_points.end();)
	{
		if(std::get<1>(*iter) - closest > eps)
			iter = found_points.erase(iter);
		else
			++iter;
	}

	return found_points;
}


} // namespace tl2
#endif
