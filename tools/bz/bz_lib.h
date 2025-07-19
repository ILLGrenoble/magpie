/**
 * brillouin zone calculation library
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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


#ifndef __TAKIN_BZLIB_H__
#define __TAKIN_BZLIB_H__

#include <vector>
#include <optional>

#include "tlibs2/libs/maths.h"
#include "libs/symops.h"
#include "libs/loadcif.h"

#if __has_include("pathslib/libs/voronoi.h")
	#include "pathslib/libs/voronoi.h"
#else
	#include "libs/voronoi.h"
#endif


/**
 * brillouin zone calculation
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
//requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
class BZCalc
{
public:
	BZCalc() = default;
	~BZCalc() = default;


	/**
	 * clear old BZ results
	 */
	void ClearBZ()
	{
		m_vertices.clear();

		m_triags.clear();
		m_triags_voroidx.clear();

		m_all_triags.clear();
		m_all_triags_voroidx.clear();
		m_all_triags_faceidx.clear();

		m_face_polygons.clear();
		m_face_norms.clear();
		m_face_dists.clear();
	}


	/**
	 * clear old BZ cut results
	 */
	void ClearBZCut()
	{
		m_cut_lines.clear();
		m_cut_lines000.clear();

		m_peaks_cut_in_plane.clear();
		m_peaks_cut_in_plane_invA.clear();
	}


	/**
	 * clear old symops
	 */
	void ClearSymOps()
	{
		m_symops.clear();
	}


	/**
	 * clear old peaks
	 */
	void ClearPeaks()
	{
		m_peaks.clear();
		m_peaks_invA.clear();
	}


	/**
	 * clear old peaks
	 */
	void ClearPeaksCut()
	{
		m_peaks_cut.clear();
		m_peaks_cut_invA.clear();
	}


	/**
	 * clear all
	 */
	void Clear()
	{
		ClearBZ();
		ClearBZCut();
		ClearSymOps();
		ClearPeaks();
		ClearPeaksCut();
	}


	// --------------------------------------------------------------------------------
	// getter and setter
	// --------------------------------------------------------------------------------
	static std::size_t GetErrIdx()
	{
		return s_erridx;
	}


	void SetEps(t_real eps)
	{
		m_eps = eps;
	}


	void SetCrystalA(const t_mat& A)
	{
		m_crystA = A;
	}


	void SetCrystalB(const t_mat& B)
	{
		m_crystB = B;
		m_crystB_ortho = orthonorm_sys<t_mat, t_vec>(m_crystB);
	}


	/**
	 * sets up a crystal lattice and angles
	 */
	void SetCrystal(t_real a, t_real b, t_real c,
		t_real alpha = 90., t_real beta = 90., t_real gamma = 90.)
	{
		t_mat crystB = tl2::B_matrix<t_mat>(a, b, c,
			tl2::d2r<t_real>(alpha), tl2::d2r<t_real>(beta), tl2::d2r<t_real>(gamma));

		SetCrystalB(crystB);
	}


	void SetPeaks(const std::vector<t_vec>& peaks, bool inv_A = false, bool for_cut = false)
	{
		if(inv_A)
		{
			if(for_cut)
				m_peaks_cut_invA = peaks;
			else
				m_peaks_invA = peaks;
		}
		else
		{
			if(for_cut)
				m_peaks_cut = peaks;
			else
				m_peaks = peaks;
		}
	}


	const std::vector<t_vec>& GetPeaks(bool inv_A = false, bool for_cut = false) const
	{
		if(inv_A)
			return for_cut ? m_peaks_cut_invA : m_peaks_invA;
		else
			return for_cut ? m_peaks_cut : m_peaks;
	}


	const std::vector<t_vec>& GetPeaksOnPlane(bool inv_A = false) const
	{
		return inv_A ? m_peaks_cut_in_plane_invA : m_peaks_cut_in_plane;
	}


	const std::vector<t_vec>& GetVertices() const
	{
		return m_vertices;
	}


	const std::vector<std::vector<t_vec>>& GetTriangles() const
	{
		return m_triags;
	}


	const std::vector<std::vector<std::size_t>>& GetTrianglesVoronoiIndices() const
	{
		return m_triags_voroidx;
	}


	const std::vector<std::vector<std::size_t>>& GetFacesIndices() const
	{
		return m_face_polygons;
	}


	const std::vector<t_vec>& GetFaceNormals() const
	{
		return m_face_norms;
	}


	const std::vector<t_real>& GetFaceDistances() const
	{
		return m_face_dists;
	}


	const std::vector<t_vec>& GetAllTriangles() const
	{
		return m_all_triags;
	}


	const std::vector<std::size_t>& GetAllTrianglesVoronoiIndices() const
	{
		return m_all_triags_voroidx;
	}


	const std::vector<std::size_t>& GetAllTrianglesFaceIndices() const
	{
		return m_all_triags_faceidx;
	}


	const t_mat& GetCrystalA() const
	{
		return m_crystA;
	}


	const t_mat& GetCrystalB(bool ortho = false) const
	{
		return ortho ? m_crystB_ortho : m_crystB;
	}


	t_real GetCutNormScale() const
	{
		return m_cut_norm_scale;
	}


	std::tuple<t_real, t_real, t_real, t_real> GetCutMinMax() const
	{
		return std::make_tuple(m_min_x, m_min_y, m_max_x, m_max_y);
	}


	/**
	 * cutting plane (or its inverse) in inverse angstroms
	 */
	const t_mat& GetCutPlane(bool inv = false) const
	{
		return inv ? m_cut_plane_inv : m_cut_plane;
	}


	/**
	 * cutting plane in inverse angstroms (or rlu)
	 */
	t_real GetCutPlaneD(bool rlu = false) const
	{
		return rlu ? m_d_rlu : m_d_invA;
	}


        const std::vector<std::tuple<t_vec, t_vec, std::array<t_real, 3>>>&
	GetCutLines(bool b000 = false) const
	{
	        // [x, y, Q]
		return b000 ? m_cut_lines000 : m_cut_lines;
	}


	/**
	 * get the index of the (000) bragg peak
	 */
	std::size_t Get000Peak(bool for_cut = false) const
	{
		const std::optional<std::size_t>& idx000 =
			for_cut ? m_idx000_cut : m_idx000;

		if(idx000)
			return *idx000;
		return s_erridx;
	}


	void AddSymOp(
		t_real r00, t_real r01, t_real r02,
		t_real r10, t_real r11, t_real r12,
		t_real r20, t_real r21, t_real r22,
		t_real t0, t_real t1, t_real t2)
	{
		t_mat symop = tl2::unit<t_mat>(4);
		// rotation
		symop(0, 0) = r00; symop(0, 1) = r01; symop(0, 2) = r02;
		symop(1, 0) = r10; symop(1, 1) = r11; symop(1, 2) = r12;
		symop(2, 0) = r20; symop(2, 1) = r21; symop(2, 2) = r22;
		// translation
		symop(0, 3) = t0;  symop(1, 3) = t1;  symop(2, 3) = t2;

		m_symops.emplace_back(std::move(symop));
	}


	/**
	 * set up a list of symmetry operations (given by the space group)
	 * @returns number of actually set symops
	 */
	std::size_t SetSymOps(const std::vector<t_mat>& ops, bool are_centring = false)
	{
		ClearSymOps();

		if(are_centring)
		{
			// symops are already purely centring, add all
			m_symops = ops;
		}
		else
		{
			// only add centring ops
			m_symops.clear();
			m_symops.reserve(ops.size());

			for(const t_mat& op : ops)
			{
				if(!tl2::hom_is_centring<t_mat>(op, m_eps))
					continue;

				m_symops.push_back(op);
			}
		}

		return m_symops.size();
	}


	std::size_t SetSymOpsFromSpaceGroup(const std::string& sgname)
	{
		std::vector<t_mat> ops = get_sg_ops<t_mat, t_real>(sgname);
		return SetSymOps(ops, false);
	}
	// --------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------
	// calculations
	// --------------------------------------------------------------------------------
	/**
	 * calculate the nuclear bragg peaks in lab coordinates
	 * @returns number of created peaks
	 */
	std::size_t CalcPeaksInvA(bool for_cut = false)
	{
		const std::vector<t_vec>& peaks = for_cut ? m_peaks_cut : m_peaks;
		std::vector<t_vec>& peaks_invA = for_cut ? m_peaks_cut_invA : m_peaks_invA;
		std::optional<std::size_t>& idx000 = for_cut ? m_idx000_cut : m_idx000;

		// calculate the peaks in lab coordinates
		peaks_invA.clear();
		peaks_invA.reserve(peaks.size());

		for(const t_vec& Q : peaks)
		{
			if(!is_reflection_allowed<t_mat, t_vec, t_real>(Q, m_symops, m_eps).first)
				continue;

			// also get the index of the (000) peak
			if(tl2::equals_0(Q, m_eps))
				idx000 = peaks_invA.size();

			peaks_invA.emplace_back(m_crystB * Q);
		}

		return peaks_invA.size();
	}


	/**
	 * create nuclear bragg peaks up to the given order
	 * @returns number of created peaks
	 */
	std::size_t CalcPeaks(int order, bool calc_invA = false, bool for_cut = false)
	{
		std::vector<t_vec>& peaks = for_cut ? m_peaks_cut : m_peaks;

		peaks.clear();
		peaks.reserve((2*order+1)*(2*order+1)*(2*order+1));

		for(int h = -order; h <= order; ++h)
		for(int k = -order; k <= order; ++k)
		for(int l = -order; l <= order; ++l)
		{
			peaks.emplace_back(
				tl2::create<t_vec>(
					{ t_real(h), t_real(k), t_real(l) }));
		}

		if(calc_invA)
			CalcPeaksInvA(for_cut);

		return peaks.size();
	}


	/**
	 * calculate the index of the nuclear (000) peak
	 */
	void Calc000Peak(bool for_cut = false)
	{
		const std::vector<t_vec>& peaks = for_cut ? m_peaks_cut : m_peaks;
		std::optional<std::size_t>& idx000 = for_cut ? m_idx000_cut : m_idx000;

		for(std::size_t idx = 0; idx < peaks.size(); ++idx)
		{
			const t_vec& Q = peaks[idx];
			if(!tl2::equals_0(Q, m_eps))
				continue;
			idx000 = idx;
			break;
		}
	}


	/**
	 * calculate the brillouin zone
	 */
	bool CalcBZ()
	{
		ClearBZ();

		if(!m_peaks.size())
			return false;
		if(!m_idx000)
			Calc000Peak();

		// calculate the voronoi diagram's vertices
		std::tie(m_vertices, std::ignore, std::ignore) =
			geo::calc_delaunay(3, m_peaks_invA, false, false, m_idx000);
		m_vertices = tl2::remove_duplicates(m_vertices, m_eps);
		if(!m_vertices.size())
			return false;

		for(t_vec& vertex : m_vertices)
			tl2::set_eps_0(vertex, m_eps);

		// calculate the faces of the BZ
		std::tie(std::ignore, m_triags, std::ignore) =
			geo::calc_delaunay(3, m_vertices, true, false);
		if(!m_triags.size())
			return false;

		// calculate all BZ triangles
		for(std::size_t triag_idx = 0; triag_idx < m_triags.size(); ++triag_idx)
		{
			std::vector<t_vec>& bz_triag = m_triags[triag_idx];
			if(bz_triag.size() < 3)
				continue;

			// calculate face plane
			t_vec norm = tl2::cross<t_vec>(
				{ bz_triag[1] - bz_triag[0], bz_triag[2] - bz_triag[0] });
			t_real norm_len = tl2::norm<t_vec>(norm);
			if(!tl2::equals_0<t_real>(norm_len, m_eps))
				norm /= tl2::norm<t_vec>(norm);

			t_real dist = tl2::inner<t_vec>(bz_triag[0], norm);
			if(dist < 0.)
			{
				norm = -norm;
				dist = -dist;

				std::reverse(bz_triag.begin(), bz_triag.end());
			}

			// find out if we've already got this face
			bool face_found = false;
			std::size_t face_idx = 0;
			for(face_idx = 0; face_idx < m_face_polygons.size(); ++face_idx)
			{
				if(tl2::equals<t_vec>(m_face_norms[face_idx], norm, m_eps)
					&& tl2::equals<t_real>(m_face_dists[face_idx], dist, m_eps))
				{
					face_found = true;
					break;
				}
			}

			// add the triangle to the face
			if(face_found)
			{
				m_face_polygons[face_idx].push_back(triag_idx);
			}
			else
			{
				m_face_polygons.emplace_back(std::vector<std::size_t>({triag_idx}));

				tl2::set_eps_0(norm, m_eps);
				tl2::set_eps_0(dist, m_eps);

				m_face_norms.emplace_back(std::move(norm));
				m_face_dists.push_back(dist);
			}

			// iterate triangle vertices
			std::vector<std::size_t> triagindices;
			for(t_vec& vert : bz_triag)
			{
				tl2::set_eps_0(vert, m_eps);

				// find index of vertex among voronoi vertices
				std::ptrdiff_t voroidx = -1;
				if(auto voro_iter = std::find_if(m_vertices.begin(), m_vertices.end(),
					[&vert, this](const t_vec& vec) -> bool
					{
						return tl2::equals<t_vec>(vec, vert, m_eps);
					}); voro_iter != m_vertices.end())
				{
					voroidx = voro_iter - m_vertices.begin();
				}

				std::size_t idx = (voroidx >= 0 ? voroidx : s_erridx);
				triagindices.push_back(idx);

				m_all_triags.push_back(vert);
				m_all_triags_voroidx.push_back(idx);
			}  // vertices

			m_all_triags_faceidx.push_back(face_idx);
			m_triags_voroidx.emplace_back(std::move(triagindices));
		}  // triangles

		return true;
	}


	/**
	 * calculate the 2d brillouin zone cut
	 */
	bool CalcBZCut(const t_vec& vec1_rlu, const t_vec& norm_rlu, t_real d_rlu,
		bool calc_bzcut_hull = true)
	{
		ClearBZCut();

		// calculate the 3d bz if not yet done beforehand
		if(!GetTriangles().size())
			CalcBZ();

		if(!m_peaks_cut.size())
			return false;
		if(!m_idx000_cut)
			Calc000Peak(true);

		m_vec1_rlu = vec1_rlu;
		m_norm_rlu = norm_rlu;
		m_d_rlu = d_rlu;

		m_vec1_rlu /= tl2::norm<t_vec>(m_vec1_rlu);
		m_norm_rlu /= tl2::norm<t_vec>(m_norm_rlu);

		t_vec vec1_invA = m_crystB * m_vec1_rlu;
		t_vec norm_invA = m_crystB * m_norm_rlu;
		m_cut_norm_scale = tl2::norm<t_vec>(norm_invA);
		norm_invA /= m_cut_norm_scale;
		m_d_invA = d_rlu*m_cut_norm_scale;

		t_vec vec2_invA = tl2::cross<t_vec>(norm_invA, vec1_invA);
		vec1_invA = tl2::cross<t_vec>(vec2_invA, norm_invA);
		vec1_invA /= tl2::norm<t_vec>(vec1_invA);
		vec2_invA /= tl2::norm<t_vec>(vec2_invA);

		t_mat B_inv = m_crystA / (t_real(2)*tl2::pi<t_real>);
		m_vec2_rlu = B_inv * vec2_invA;
		m_vec2_rlu /= tl2::norm<t_vec>(m_vec2_rlu);

		m_cut_plane = tl2::create<t_mat, t_vec>({ vec1_invA, vec2_invA, norm_invA }, false);
		m_cut_plane_inv = tl2::trans<t_mat>(m_cut_plane);

		for(std::size_t idx = 0; idx < m_peaks_cut.size(); ++idx)
		{
			const t_vec& Q = m_peaks_cut[idx];

			if(!is_reflection_allowed<t_mat, t_vec, t_real>(
				Q, m_symops, m_eps).first)
				continue;

			// (000) peak?
			bool is_000 = false;
			if(m_idx000_cut)
				is_000 = (idx == *m_idx000_cut);
			else
				is_000 = tl2::equals_0(Q, m_eps);

			t_vec Q_invA = m_crystB * Q;

			// Q on cutting plane?
			if(tl2::equals(tl2::inner(norm_invA, Q_invA), m_d_invA, m_eps))
			{
				m_peaks_cut_in_plane.push_back(Q);
				m_peaks_cut_in_plane_invA.push_back(m_cut_plane_inv * Q_invA);
			}

			std::vector<t_vec> cut_verts;
			std::optional<t_real> z_comp;

			for(const auto& _bz_poly : GetTriangles())
			{
				// centre bz around bragg peak
				auto bz_poly = _bz_poly;
				for(t_vec& vec : bz_poly)
					vec += Q_invA;

				auto vecs = tl2::intersect_plane_poly<t_vec>(
					norm_invA, m_d_invA, bz_poly, m_eps);
				vecs = tl2::remove_duplicates(vecs, m_eps);

				// calculate the hull of the bz cut
				if(calc_bzcut_hull)
				{
					for(const t_vec& vec : vecs)
					{
						t_vec vec_rot = m_cut_plane_inv * vec;
						tl2::set_eps_0(vec_rot, m_eps);

						cut_verts.emplace_back(
							tl2::create<t_vec>({
								vec_rot[0],
								vec_rot[1] }));

						// z component is the same for every vector
						if(!z_comp)
							z_comp = vec_rot[2];
					}
				}
				// alternatively use the lines directly
				else if(vecs.size() >= 2)
				{
					t_vec pt1 = m_cut_plane_inv * vecs[0];
					t_vec pt2 = m_cut_plane_inv * vecs[1];
					tl2::set_eps_0(pt1, m_eps);
					tl2::set_eps_0(pt2, m_eps);

					m_cut_lines.emplace_back(std::make_tuple(
						pt1, pt2,
						std::array<t_real,3>{Q[0], Q[1], Q[2]}));
					if(is_000)
					{
						m_cut_lines000.emplace_back(std::make_tuple(
							pt1, pt2,
							std::array<t_real,3>{Q[0], Q[1], Q[2]}));
					}
				}
			}

			// calculate the hull of the bz cut
			if(calc_bzcut_hull)
			{
				cut_verts = tl2::remove_duplicates(cut_verts, m_eps);
				if(cut_verts.size() < 3)
					continue;

				// calculate the faces of the BZ
				auto [bz_verts, bz_triags, bz_neighbours] =
					geo::calc_delaunay(2, cut_verts, true, false);

				for(std::size_t bz_idx = 0; bz_idx < bz_verts.size(); ++bz_idx)
				{
					std::size_t bz_idx2 = (bz_idx+1) % bz_verts.size();
					t_vec pt1 = tl2::create<t_vec>({
						bz_verts[bz_idx][0],
						bz_verts[bz_idx][1],
						z_comp ? *z_comp : t_real(0.) });
					t_vec pt2 = tl2::create<t_vec>({
						bz_verts[bz_idx2][0],
						bz_verts[bz_idx2][1],
						z_comp ? *z_comp : t_real(0.) });
					tl2::set_eps_0(pt1, m_eps);
					tl2::set_eps_0(pt2, m_eps);

					m_cut_lines.emplace_back(std::make_tuple(
						pt1, pt2,
						std::array<t_real,3>{Q[0], Q[1], Q[2]}));
					if(is_000)
					{
						m_cut_lines000.emplace_back(std::make_tuple(
							pt1, pt2,
							std::array<t_real,3>{Q[0], Q[1], Q[2]}));
					}
				}
			}
		}

		// get ranges
		m_min_x = std::numeric_limits<t_real>::max();
		m_max_x = -m_min_x;
		m_min_y = std::numeric_limits<t_real>::max();
		m_max_y = -m_min_y;

		for(const auto& tup : m_cut_lines)
		{
			const auto& pt1 = std::get<0>(tup);
			const auto& pt2 = std::get<1>(tup);

			m_min_x = std::min(m_min_x, pt1[0]);
			m_min_x = std::min(m_min_x, pt2[0]);
			m_max_x = std::max(m_max_x, pt1[0]);
			m_max_x = std::max(m_max_x, pt2[0]);

			m_min_y = std::min(m_min_y, pt1[1]);
			m_min_y = std::min(m_min_y, pt2[1]);
			m_max_y = std::max(m_max_y, pt1[1]);
			m_max_y = std::max(m_max_y, pt2[1]);
		}

		return true;
	}
	// --------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------
	// output
	// --------------------------------------------------------------------------------
	/**
	 * print a description of the bz
	 */
	std::string Print(int prec = 6) const
	{
		using namespace tl2_ops;
		const t_mat& B = GetCrystalB(false);
		t_mat BorthoT = tl2::trans(GetCrystalB(true));

		std::ostringstream ostr;
		ostr.precision(prec);

#ifdef DEBUG
		ostr << "# centring symmetry operations\n";
		for(const t_mat& op : m_symops)
			ostr << op << "\n";
#endif

		// voronoi vertices forming the vertices of the bz
		// (rotated by orthonormal part of B^(-1))
		const std::vector<t_vec>& voronoiverts = GetVertices();
		ostr << "# Brillouin zone vertices (Å⁻¹, B⁻¹-rotated)\n";
		for(std::size_t idx = 0; idx < voronoiverts.size(); ++idx)
		{
			t_vec voro = BorthoT * voronoiverts[idx];
			tl2::set_eps_0(voro, m_eps);

			ostr << "vertex " << idx << ": (" << voro << ")\n";
		}

		ostr << "\n# Brillouin zone vertices (Å⁻¹, raw)\n";
		for(std::size_t idx = 0; idx < voronoiverts.size(); ++idx)
		{
			const t_vec& voro = voronoiverts[idx];
			ostr << "raw vertex " << idx << ": (" << voro << ")\n";
		}

		// voronoi bisectors / polygons
		const auto& bz_polys = GetTriangles();
		const auto& bz_polys_idx = GetTrianglesVoronoiIndices();

		ostr << "\n# Brillouin zone polygons (Å⁻¹)\n";
		for(std::size_t idx_triag = 0; idx_triag < bz_polys.size(); ++idx_triag)
		{
			const auto& triag = bz_polys[idx_triag];
			const auto& triag_idx = bz_polys_idx[idx_triag];

			ostr << "polygon " << idx_triag << ": \n";
			for(std::size_t idx_vert = 0; idx_vert < triag.size(); ++idx_vert)
			{
				const t_vec& vert = triag[idx_vert];
				std::size_t voroidx = triag_idx[idx_vert];

				ostr << "\tvertex " << voroidx << ": (" << vert << ")\n";
			}
		}

		// faces and their planes
		const auto& bz_faces_idx = GetFacesIndices();
		const std::vector<t_vec>& norms = GetFaceNormals();
		const std::vector<t_real>& dists = GetFaceDistances();

		ostr << "\n# Brillouin zone faces (Å⁻¹):\n";
		for(std::size_t idx_face = 0; idx_face < bz_faces_idx.size(); ++idx_face)
		{
			const auto& triag_idx = bz_faces_idx[idx_face];
			const t_vec& norm = norms[idx_face];
			t_real dist = dists[idx_face];

			ostr << "face " << idx_face << ": \n";
			ostr << "\tplane normal: " << "(" << norm[0] << ", " << norm[1] << ", " << norm[2] << ")\n";
			ostr << "\tplane distance: " << dist << "\n";
			ostr << "\tpolygons: ";
			for(std::size_t idx_poly = 0; idx_poly < triag_idx.size(); ++idx_poly)
			{
				std::size_t poly_idx = triag_idx[idx_poly];

				ostr << poly_idx;
				if(idx_poly < triag_idx.size() - 1)
					ostr << ", ";
			}
			ostr << "\n";
		}

		// reciprocal crystal B matrix
		ostr << "\n# B matrix:\n";
		ostr << "# The B matrix transforms from Q in rlu to the orthogonal lab system in Å⁻¹: Q_Å⁻¹ = B · Q_rlu\n";
		ostr << "# The columns of the B matrix are the reciprocal basis vectors.\n";
		for(std::size_t i = 0; i < B.size1(); ++i)
		{
			ostr << "\t";
			for(std::size_t j = 0; j < B.size2(); ++j)
			{
				t_real elem = B(i, j);
				tl2::set_eps_0(elem, m_eps);
				ostr << elem;
				if(j < B.size2() - 1)
					ostr << " ";
			}
			ostr << "\n";
		}

		// real crystal A matrix
		if(auto [A, A_ok] = tl2::inv(B); A_ok)
		{
			ostr << "\n# A matrix:\n";
			ostr << "# The columns of the A matrix are the real basis vectors.\n";
			A = tl2::trans(A);
			A *= t_real(2) * tl2::pi<t_real>;

			for(std::size_t i = 0; i < A.size1(); ++i)
			{
				ostr << "\t";
				for(std::size_t j = 0; j < A.size2(); ++j)
				{
					t_real elem = A(i, j);
					tl2::set_eps_0(elem, m_eps);
					ostr << elem;
					if(j < A.size2() - 1)
						ostr << " ";
				}
				ostr << "\n";
			}
		}

		return ostr.str();
	}


	/**
	 * export a description of the bz in json format
	 */
	std::string PrintJSON(int prec = 6) const
	{
		const t_mat& B = GetCrystalB(false);
		const t_mat& Bortho = GetCrystalB(true);
		t_mat BorthoT = tl2::trans(Bortho);

		std::ostringstream ostr;
		ostr.precision(prec);

		ostr << "{\n";

		// voronoi vertices forming the vertices of the bz
		// (rotated by orthonormal part of B^(-1))
		const std::vector<t_vec>& voronoiverts = GetVertices();
		ostr << "\"vertices\" : [\n";
		for(std::size_t idx = 0; idx < voronoiverts.size(); ++idx)
		{
			t_vec voro = BorthoT * voronoiverts[idx];
			tl2::set_eps_0(voro, m_eps);

			ostr << "\t[ " << voro[0] << ", " << voro[1] << ", " << voro[2] << " ]";
			if(idx < voronoiverts.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// raw voronoi vertices
		ostr << "\"vertices_nonrot\" : [\n";
		for(std::size_t idx = 0; idx < voronoiverts.size(); ++idx)
		{
			const t_vec& voro = voronoiverts[idx];
			ostr << "\t[ " << voro[0] << ", " << voro[1] << ", " << voro[2] << " ]";
			if(idx < voronoiverts.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// voronoi bisectors / polygons
		const auto& bz_polys_idx = GetTrianglesVoronoiIndices();
		ostr << "\"polygons\" : [\n";
		for(std::size_t idx_triag = 0; idx_triag < bz_polys_idx.size(); ++idx_triag)
		{
			const auto& triag_idx = bz_polys_idx[idx_triag];

			ostr << "\t[ ";
			for(std::size_t idx_vert = 0; idx_vert < triag_idx.size(); ++idx_vert)
			{
				ostr << triag_idx[idx_vert];
				if(idx_vert < triag_idx.size() - 1)
					ostr << ", ";
			}
			ostr << " ]";
			if(idx_triag < bz_polys_idx.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// faces
		const auto& bz_faces_idx = GetFacesIndices();
		ostr << "\"faces\" : [\n";
		for(std::size_t idx_face = 0; idx_face < bz_faces_idx.size(); ++idx_face)
		{
			const auto& triag_idx = bz_faces_idx[idx_face];

			ostr << "\t[ ";
			for(std::size_t idx_poly = 0; idx_poly < triag_idx.size(); ++idx_poly)
			{
				ostr << triag_idx[idx_poly];
				if(idx_poly < triag_idx.size() - 1)
					ostr << ", ";
			}
			ostr << " ]";
			if(idx_face < bz_faces_idx.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// face planes
		const std::vector<t_vec>& norms = GetFaceNormals();
		const std::vector<t_real>& dists = GetFaceDistances();
		ostr << "\"face_planes\" : [\n";
		for(std::size_t idx = 0; idx < norms.size(); ++idx)
		{
			const t_vec& norm = norms[idx];
			t_real dist = dists[idx];

			ostr << "\t[ " << norm[0] << ", " << norm[1] << ", " << norm[2] << ", " << dist << " ]";
			if(idx < norms.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// crystal B matrix
		ostr << "\"crystal_B\" : [ ";
		for(std::size_t i = 0; i < B.size1(); ++i)
		{
			for(std::size_t j = 0; j < B.size2(); ++j)
			{
				t_real elem = B(i, j);
				tl2::set_eps_0(elem, m_eps);
				ostr << elem;
				if(i < B.size1() - 1 || j < B.size2() - 1)
					ostr << ",";
				ostr << " ";
			}
		}
		ostr << "],\n\n";

		// orthonormal (rotation) part of crystal B matrix
		ostr << "\"crystal_B_ortho\" : [ ";
		for(std::size_t i = 0; i < Bortho.size1(); ++i)
		{
			for(std::size_t j = 0; j < Bortho.size2(); ++j)
			{
				t_real elem = Bortho(i, j);
				tl2::set_eps_0(elem, m_eps);
				ostr << elem;
				if(i < Bortho.size1() - 1 || j < Bortho.size2() - 1)
					ostr << ",";
				ostr << " ";
			}
		}
		ostr << "]\n";

		ostr << "}\n";
		return ostr.str();
	}


	/**
	 * print a description of the bz cut
	 */
	std::string PrintCut(int prec = 6) const
	{
		using namespace tl2_ops;

		t_vec vec1_invA = tl2::col<t_mat, t_vec>(GetCutPlane(), 0);
		t_vec vec2_invA = tl2::col<t_mat, t_vec>(GetCutPlane(), 1);
		t_vec norm_invA = tl2::col<t_mat, t_vec>(GetCutPlane(), 2);
		t_vec vec1_rlu = m_vec1_rlu, vec2_rlu = m_vec2_rlu, norm_rlu = m_norm_rlu;
		t_real d_rlu = m_d_rlu, d_invA = m_d_invA;

		tl2::set_eps_0(norm_invA, m_eps); tl2::set_eps_0(norm_rlu, m_eps);
		tl2::set_eps_0(vec1_invA, m_eps); tl2::set_eps_0(vec1_rlu, m_eps);
		tl2::set_eps_0(vec2_invA, m_eps); tl2::set_eps_0(vec2_rlu, m_eps);

		std::ostringstream ostr;
		ostr.precision(prec);

		// description of the cutting plane
		ostr << "# Cutting plane";
		ostr << "\nin relative lattice units:";
		ostr << "\n\tnormal: [" << norm_rlu << "] rlu";
		ostr << "\n\tin-plane vector 1: [" << vec1_rlu << "] rlu";
		ostr << "\n\tin-plane vector 2: [" << vec2_rlu << "] rlu";
		ostr << "\n\tplane offset: " << d_rlu << " rlu";

		ostr << "\nin lab units:";
		ostr << "\n\tnormal: [" << norm_invA << "] Å⁻¹";
		ostr << "\n\tin-plane vector 1: [" << vec1_invA << "] Å⁻¹";
		ostr << "\n\tin-plane vector 2: [" << vec2_invA << "] Å⁻¹";
		ostr << "\n\tplane offset: " << d_invA << " Å⁻¹";
		ostr << "\n" << std::endl;

		// description of the bz cut
		ostr << "# Brillouin zone cut (Å⁻¹)" << std::endl;
		for(std::size_t i = 0; i < GetCutLines(true).size(); ++i)
		{
			const auto& line = GetCutLines(true)[i];

			ostr << "line " << i
				<< ":\n\tvertex 0: (" << std::get<0>(line) << ")"
				<< "\n\tvertex 1: (" << std::get<1>(line) << ")"
				<< std::endl;
		}

		return ostr.str();
	}
	// --------------------------------------------------------------------------------


private:
	// --------------------------------------------------------------------------------
	// inputs
	// --------------------------------------------------------------------------------
	t_real m_eps{ 1e-7 };                          // calculation epsilon

	t_mat m_crystA{ tl2::unit<t_mat>(3) };         // crystal A matrix
	t_mat m_crystB{ tl2::unit<t_mat>(3) };         // crystal B matrix
	t_mat m_crystB_ortho{ tl2::unit<t_mat>(3) };   // orthonormal part of crystal B matrix

	std::vector<t_mat> m_symops{ };                // space group centring symmetry operations
	std::vector<t_vec> m_peaks{ };                 // nuclear bragg peaks
	std::vector<t_vec> m_peaks_cut{ };             // nuclear bragg peaks for zone cut
	std::vector<t_vec> m_peaks_invA { };           // nuclear bragg peaks in lab coordinates
	std::vector<t_vec> m_peaks_cut_invA{ };        // nuclear bragg peaks for zone cut in lab coordinates
	std::optional<std::size_t> m_idx000{};         // index of the (000) peak
	std::optional<std::size_t> m_idx000_cut{};     // index of the (000) peak for zone cut
	// --------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------
	// calculated brillouin zone infos
	// --------------------------------------------------------------------------------
	std::vector<t_vec> m_vertices{ };              // voronoi vertices

	std::vector<std::vector<t_vec>> m_triags{};    // bz triangles
	std::vector<std::vector<std::size_t>> m_triags_voroidx{}; // voronoi vertex indices

	std::vector<t_vec> m_all_triags{ };                // all brillouin zone triangles
	std::vector<std::size_t> m_all_triags_voroidx{ };  // voronoi vertex indices
	std::vector<std::size_t> m_all_triags_faceidx{ };  // face indices

	std::vector<std::vector<std::size_t>> m_face_polygons{ };
	std::vector<t_vec> m_face_norms{ };
	std::vector<t_real> m_face_dists{ };
	// --------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------
	// calculated brillouin zone cut infos
	// --------------------------------------------------------------------------------
	t_vec m_vec1_rlu{}, m_vec2_rlu{}, m_norm_rlu{};  // cutting plane in rlu
	t_real m_cut_norm_scale { 1. };              // convert 1/A to rlu lengths along the normal
	t_mat m_cut_plane = tl2::unit<t_mat>(3);     // cutting plane in inverse angstroms
	t_mat m_cut_plane_inv = tl2::unit<t_mat>(3); // ...and its inverse
	t_real m_d_rlu{ 0. }, m_d_invA{ 0. };        // cutting plane distance

	// [x, y, Q]
	std::vector<std::tuple<t_vec, t_vec, std::array<t_real, 3>>> m_cut_lines{}, m_cut_lines000{};

	t_real m_min_x{ 1. }, m_max_x{ -1. };          // plot ranges for curves
	t_real m_min_y{ 1. }, m_max_y{ -1. };          // plot ranges for curves

	std::vector<t_vec> m_peaks_cut_in_plane{ };      // bragg peaks on the cutting plane in rlu
	std::vector<t_vec> m_peaks_cut_in_plane_invA{ }; // ... and in 1/A and transformed into the plane system
	// --------------------------------------------------------------------------------

	static const std::size_t s_erridx{0xffffffff}; // index for reporting errors
};


#endif
