/**
 * magnetic dynamics -- 3d dispersion plotting
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#include <cstdlib>
#include <sstream>

#include "dispersion3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * get a unique colour for a given magnon band
 */
std::array<int, 3> Dispersion3DDlg::GetBranchColour(t_size branch_idx, t_size num_branches) const
{
	if(num_branches <= 1)
		return std::array<int, 3>{{ 0xff, 0x00, 0x00 }};

	return std::array<int, 3>{{
		int(std::lerp(1., 0., t_real(branch_idx) / t_real(num_branches - 1)) * 255.),
		0x00,
		int(std::lerp(0., 1., t_real(branch_idx) / t_real(num_branches - 1)) * 255.),
	}};
}



/**
 * plot the calculated dispersion
 */
void Dispersion3DDlg::Plot(bool clear_settings)
{
	if(!m_dispplot || !m_dispplot->GetRenderer())
		return;

	// keep some settings from previous plot, e.g. the band visibility flags
	std::vector<bool> active_bands;
	if(!clear_settings)
	{
		active_bands.reserve(m_table_bands->rowCount());
		for(int row = 0; row < m_table_bands->rowCount(); ++row)
			active_bands.push_back(IsBandEnabled(t_size(row)));
	}

	bool use_E_min = m_enable_E_range[0]->isChecked();
	bool use_E_max = m_enable_E_range[1]->isChecked();
	t_real E_min_sel = m_E_range[0]->value();
	t_real E_max_sel = m_E_range[1]->value();

	t_real E_scale = m_E_scale->value();
	t_real Q_scale1 = m_Q_scale1->value();
	t_real Q_scale2 = m_Q_scale2->value();

	// reset plot
	ClearBands();
	m_cam_centre = tl2::zero<t_vec_gl>(3);
	m_dispplot->GetRenderer()->RemoveObjects();
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * Q_scale1), static_cast<t_real_gl>(1.5 * Q_scale2),
			static_cast<t_real_gl>(2.5 * (m_minmax_E[1] + 4.) * E_scale) }));
	m_dispplot->GetRenderer()->SetLight(1, -tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * Q_scale1), static_cast<t_real_gl>(1.5 * Q_scale2),
			static_cast<t_real_gl>(2.5 * (m_minmax_E[1] + 4.) * E_scale) }));


	// set coordinate cube size
	if(m_dispplot->GetRenderer()->GetCoordCube().size())
	{
		t_real E_min = use_E_min ? E_min_sel : m_minmax_E[0];
		t_real E_max = use_E_max ? E_max_sel : m_minmax_E[1];
		t_real E_range = E_max - E_min;
		t_real E_mean = E_min + E_range*0.5;

		t_mat_gl obj_scale = tl2::hom_scaling<t_mat_gl>(
			0.5 * Q_scale1, 0.5 * Q_scale2, 0.5 * E_scale * E_range);
		t_mat_gl obj_shift = tl2::hom_translation<t_mat_gl>(t_real(0.), t_real(0.), E_scale * E_mean);

		for(auto obj : m_dispplot->GetRenderer()->GetCoordCube())
		{
			using namespace tl2_ops;
			m_dispplot->GetRenderer()->SetObjectMatrix(obj, obj_shift * obj_scale);
		}

		auto [Q_idx_1, Q_idx_2] = GetQIndices();
		/*static const std::array<std::string, 3> labels
		{{
			"h (rlu)", "k (rlu)", "l (rlu)"
		}};*/

		if(m_dispplot && m_minmax_Q1[0].size() == 3)
		{
			//using namespace tl2_ops;
			//std::cout << m_minmax_Q1[0] << " -> " << m_minmax_Q1[1] << ", " << Q_idx_1 << std::endl;
			//std::cout << m_minmax_Q2[0] << " -> " << m_minmax_Q2[1] << ", " << Q_idx_2 << std::endl;
			m_dispplot->GetRenderer()->UpdateCoordCubeTextures(
				m_minmax_Q1[0][Q_idx_1], m_minmax_Q1[1][Q_idx_1], -1.,
				m_minmax_Q2[0][Q_idx_2], m_minmax_Q2[1][Q_idx_2], -1.,
				E_min, E_max, -1.);
		}
	}


	// plot the magnon bands
	m_first_band = 0;
	t_size bands_end = m_data.size();
	if(use_E_min || use_E_max)
		std::tie(m_first_band, bands_end) = BandIndicesInRange();

	t_size num_active_bands = 0;
	for(t_size band_idx = m_first_band; band_idx < bands_end; ++band_idx)
	{
		bool band_active = band_idx - m_first_band < active_bands.size()
			? active_bands[band_idx - m_first_band] : true;

		// colour for this magnon band
		std::array<int, 3> col = GetBranchColour(band_idx - m_first_band, bands_end - m_first_band);
		const QColor colFull(col[0], col[1], col[2]);

		if(band_active)
		{
			t_data_Qs& data = m_data[band_idx];

			// dispersion surface function
			auto patch_fkt = [this, &data, E_scale](
				t_real_gl /*x2*/, t_real_gl /*x1*/, t_size idx_1, t_size idx_2)
					-> std::pair<t_real_gl, bool>
			{
				t_size idx = idx_1 * m_Q_count_2 + idx_2;
				if(idx >= data.size())
					return std::make_pair(0., false);

				bool valid = std::get<6>(data[idx]);
				if(std::get<3>(data[idx]) != idx_1 || std::get<4>(data[idx]) != idx_2)
				{
					std::cerr << "Error: Patch index mismatch: "
						<< "Expected " << std::get<3>(data[idx]) << " for x, but got " << idx_1
						<< "; expected " << std::get<4>(data[idx]) << " for y, but got " << idx_2
						<< "." << std::endl;
					valid = false;
				}

				t_real_gl E = std::get<1>(data[idx]);
				return std::make_pair(E * E_scale, valid);
			};

			// dispersion surface border line 1
			auto line1_fkt = [&patch_fkt](
				t_real_gl x, t_size idx) -> std::pair<t_real_gl, bool>
			{
				return patch_fkt(x, 0., idx, 0);
			};

			// dispersion surface border line 2
			auto line2_fkt = [this, &patch_fkt](
				t_real_gl x, t_size idx) -> std::pair<t_real_gl, bool>
			{
				return patch_fkt(x, 0., idx, m_Q_count_2 - 1);
			};

			// dispersion surface border line 3
			auto line3_fkt = [&patch_fkt](
				t_real_gl y, t_size idx) -> std::pair<t_real_gl, bool>
			{
				return patch_fkt(0., y, 0, idx);
			};

			// dispersion surface border line 4
			auto line4_fkt = [this, &patch_fkt](
				t_real_gl y, t_size idx) -> std::pair<t_real_gl, bool>
			{
				return patch_fkt(0., y, m_Q_count_1 - 1, idx);
			};

			std::ostringstream objLabel;
			objLabel << "Band #" << (band_idx + 1);

			// colours
			t_real_gl r = t_real_gl(col[0]) / t_real_gl(255.);
			t_real_gl g = t_real_gl(col[1]) / t_real_gl(255.);
			t_real_gl b = t_real_gl(col[2]) / t_real_gl(255.);
			t_real_gl r_line = r * t_real_gl(0.2);
			t_real_gl g_line = g * t_real_gl(0.2);
			t_real_gl b_line = b * t_real_gl(0.2);

			// surface
			std::size_t obj = m_dispplot->GetRenderer()->AddPatch(patch_fkt, 0., 0., 0.,
				Q_scale1, Q_scale2, m_Q_count_1, m_Q_count_2, r, g, b, 1.);
			m_dispplot->GetRenderer()->SetObjectLabel(obj, objLabel.str());
			m_band_objs.insert(std::make_pair(obj, band_idx));

			// border line 1
			t_real_gl y1 = -Q_scale2*0.5;
			/*std::size_t objLine1 =*/ m_dispplot->GetRenderer()->AddLine(line1_fkt, 0., 0., 0.,
				Q_scale1, m_Q_count_1, r_line, g_line, b_line, 1., y1, false);

			// border line 2
			t_real_gl y2 = -Q_scale2*0.5 + Q_scale2;
			/*std::size_t objLine2 =*/ m_dispplot->GetRenderer()->AddLine(line2_fkt, 0., 0., 0.,
				Q_scale1, m_Q_count_1, r_line, g_line, b_line, 1., y2, false);

			// border line 3
			t_real_gl x1 = -Q_scale1*0.5;
			/*std::size_t objLine3 =*/ m_dispplot->GetRenderer()->AddLine(line3_fkt, 0., 0., 0.,
				Q_scale2, m_Q_count_2, r_line, g_line, b_line, 1., x1, true);

			// border line 4
			t_real_gl x2 = -Q_scale1*0.5 + Q_scale1;
			/*std::size_t objLine4 =*/ m_dispplot->GetRenderer()->AddLine(line4_fkt, 0., 0., 0.,
				Q_scale2, m_Q_count_2, r_line, g_line, b_line, 1., x2, true);

			m_cam_centre[2] += GetMeanEnergy(band_idx);
			++num_active_bands;
		}  // band active

		AddBand("#" + tl2::var_to_str(band_idx + 1 - m_first_band), colFull, band_active);
	}  // bands


	if(m_show_main_Q && m_Qstart.size() == 3 && m_Qend.size() == 3)
	{
		// indicate dispersion from main dialog
		t_vec_real mainDir = m_Qend - m_Qstart;
		mainDir /= tl2::norm(mainDir);

		//const t_mat_real& xtalB = m_dyn->GetCrystalBTrafo();
		//const t_vec_real* plane = m_dyn->GetScatteringPlane();

		auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();
		Q_dir_1 /= tl2::norm(Q_dir_1);
		Q_dir_2 /= tl2::norm(Q_dir_2);

		t_real proj1 = tl2::inner(mainDir, Q_dir_1);
		t_real proj2 = tl2::inner(mainDir, Q_dir_2);
		t_real posproj1 = tl2::inner(m_Qstart - Q_origin, Q_dir_1);
		t_real posproj2 = tl2::inner(m_Qstart - Q_origin, Q_dir_2);
		t_real Q_scale = proj1*Q_scale1 + proj2*Q_scale2;

		t_real E_min = use_E_min ? E_min_sel : m_minmax_E[0];
		t_real E_max = use_E_max ? E_max_sel : m_minmax_E[1];
		t_real E_range = E_max - E_min;

		std::size_t plane = m_dispplot->GetRenderer()->AddPlane(
		  proj2, -proj1, 0.,  // flipping components for normal
			(posproj1 - 0.5*proj2) * Q_scale1,  // x
			(posproj2 - 0.5*proj1) * Q_scale2,  // y
			(E_min + 0.5*E_range) * E_scale,    // z
			Q_scale/2., E_scale*E_range/2., 0.75, 0.75, 0.75, 0.5,
			true);

		m_dispplot->GetRenderer()->SetObjectVisible(plane, true);
		m_dispplot->GetRenderer()->SetObjectPriority(plane, 0);
	}


	if(num_active_bands)
		m_cam_centre[2] /= static_cast<t_real>(num_active_bands);

	if(clear_settings)
		CentrePlotCamera();    // centre camera and update
	else
		m_dispplot->update();  // only update renderer
}



/**
 * conversion from (Q, E) to plot position
 */
t_vec_real Dispersion3DDlg::QEToPlotXYZ(const t_vec_real& Q, t_real E) const
{
	// coordinate scaling
	t_real E_scale = m_E_scale->value();
	t_real Q_scale1 = m_Q_scale1->value();
	t_real Q_scale2 = m_Q_scale2->value();

	// get plot coordinate system
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	// Q = Q_origin + Q_dir_1*(0.5*Q_scale1 + x) / Q_scale1 + Q_dir_2*(0.5*Q_scale2 + y) / Q_scale2;
	// Q - Q_origin = Q_dir_1*(0.5*Q_scale1 + x) / Q_scale1 + Q_dir_2*(0.5*Q_scale2 + y) / Q_scale2;
	// Q - Q_origin = Q_dir_1*(0.5 + x/Q_scale1) + Q_dir_2*(0.5 + y/Q_scale2);
	// Q - Q_origin = Q_dir_1*0.5 + x/Q_scale1*Q_dir1 + Q_dir_2*0.5 + y/Q_scale2*Q_dir2;
	// Q - Q_origin - Q_dir_1*0.5 - Q_dir_2*0.5 = x*Q_dir1/Q_scale1 + y*Q_dir2/Q_scale2;
	// Q - Q_origin - Q_dir_1*0.5 - Q_dir_2*0.5 = (Q_dir1/Q_scale1, Q_dir2/Q_scale2) * (x y)^T
	// (Q_dir1/Q_scale1, Q_dir2/Q_scale2)^(-1) (Q - Q_origin - Q_dir_1*0.5 - Q_dir_2*0.5) = (x y)^T

	t_mat_real mat = tl2::create<t_mat_real, t_vec_real>({ Q_dir_1/Q_scale1, Q_dir_2/Q_scale2 });
	t_vec_real vec = Q - Q_origin - Q_dir_1*0.5 - Q_dir_2*0.5;

	auto [xy, ok] = tl2::leastsq<t_vec_real, t_mat_real>(mat, vec);
	if(!ok)
		return t_vec_real{};

	return tl2::create<t_vec_real>({ xy[0], xy[1], E * E_scale });
}



/**
 * conversion from plot position to (Q, E)
 */
std::pair<t_vec_real, t_real> Dispersion3DDlg::PlotXYZToQE(t_real x, t_real y, t_real z) const
{
	// coordinate scaling
	t_real E_scale = m_E_scale->value();
	t_real Q_scale1 = m_Q_scale1->value();
	t_real Q_scale2 = m_Q_scale2->value();

	// get plot coordinate system
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	// reconstruct Q position
	t_real Q1param = (0.5*Q_scale1 + x) / Q_scale1;
	t_real Q2param = (0.5*Q_scale2 + y) / Q_scale2;

	const t_vec_real Q = Q_origin + Q_dir_1*Q1param + Q_dir_2*Q2param;
	const t_real E = z / E_scale;

	return std::make_pair(Q, E);
}



/**
 * mouse intersection with dispersion band
 */
void Dispersion3DDlg::PlotPickerIntersection(
	[[maybe_unused]] const t_vec3_gl* pos,
	[[maybe_unused]] std::size_t objIdx,
	[[maybe_unused]] std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl* posSphere)
{
	m_status->setText("");
	m_cur_obj = std::nullopt;

	m_dispplot->GetRenderer()->SetObjectsHighlight(false);

	if(!pos || !m_dyn)
		return;

	m_cur_obj = objIdx;

	auto [Q, E] = PlotXYZToQE((*pos)[0], (*pos)[1], (*pos)[2]);

	// test
	//using namespace tl2_ops;
	//t_vec_real xyz = QEToPlotXYZ(Q, E);
	//std::cout << xyz << "; " << (*pos)[0] << " " << (*pos)[1] << " " << (*pos)[2] << std::endl;

	const t_mat_real& B = m_dyn->GetCrystalBTrafo();
	const t_vec_real Qvec_invA = B * Q;
	const t_real Q_invA = tl2::norm(Qvec_invA);

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr
		<< "Q = (" << Q[0] << ", " << Q[1] << ", " << Q[2] << ") rlu, "
		<< "|Q| = " << Q_invA << " Å⁻¹, "
		<< "E = " << E << " meV";

	const std::string& label = m_dispplot->GetRenderer()->GetObjectLabel(objIdx);
	if(label != "")
		ostr << ", " << label;
	ostr << ".";

	m_status->setText(ostr.str().c_str());
}
