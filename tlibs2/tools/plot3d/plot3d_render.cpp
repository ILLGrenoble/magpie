/**
 * 3d plotting
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "plot3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * get a unique colour for a given surface
 */
std::array<int, 3> Plot3DDlg::GetSurfColour(t_size surf_idx, t_size num_surfs) const
{
	if(num_surfs <= 1)
		return std::array<int, 3>{{ 0xff, 0x00, 0x00 }};

	return std::array<int, 3>{{
		int(std::lerp(1., 0., t_real(surf_idx) / t_real(num_surfs - 1)) * 255.),
		0x00,
		int(std::lerp(0., 1., t_real(surf_idx) / t_real(num_surfs - 1)) * 255.),
	}};
}



/**
 * plot the calculated surfaces
 */
void Plot3DDlg::Plot(bool clear_settings)
{
	if(!m_dispplot)
		return;

	// keep some settings from previous plot, e.g. the surface visibility flags
	std::vector<bool> active_surfs;
	if(!clear_settings)
	{
		active_surfs.reserve(m_table_surfs->rowCount());
		for(int row = 0; row < m_table_surfs->rowCount(); ++row)
			active_surfs.push_back(IsSurfaceEnabled(t_size(row)));
	}

	t_real x_scale = m_x_scale->value();
	t_real y_scale = m_y_scale->value();
	t_real z_scale = m_z_scale->value();

	// reset plot
	ClearSurfaces();
	m_cam_centre = tl2::zero<t_vec_gl>(3);
	m_dispplot->GetRenderer()->RemoveObjects();
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * x_scale), static_cast<t_real_gl>(1.5 * y_scale),
			static_cast<t_real_gl>(2.5 * (m_minmax_z[1] + 4.) * z_scale) }));
	m_dispplot->GetRenderer()->SetLight(1, -tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * x_scale), static_cast<t_real_gl>(1.5 * y_scale),
			static_cast<t_real_gl>(2.5 * (m_minmax_z[1] + 4.) * z_scale) }));

	// set coordinate cube size
	if(m_dispplot->GetRenderer()->GetCoordCube().size())
	{
		t_real E_range = m_minmax_z[1];
		t_real E_min = m_minmax_z[0];
		E_range -= E_min;

		t_real E_mean = m_minmax_z[0] + 0.5*E_range;

		t_mat_gl obj_scale = tl2::hom_scaling<t_mat_gl>(
			0.5 * x_scale, 0.5 * y_scale, 0.5 * z_scale * E_range);
		t_mat_gl obj_shift = tl2::hom_translation<t_mat_gl>(t_real(0.), t_real(0.), z_scale * E_mean);

		for(auto obj : m_dispplot->GetRenderer()->GetCoordCube())
		{
			using namespace tl2_ops;
			m_dispplot->GetRenderer()->SetObjectMatrix(obj, obj_shift * obj_scale);
		}

		m_dispplot->GetRenderer()->UpdateCoordCubeTextures(
			m_minmax_x[0], m_minmax_x[1], -1.,
			m_minmax_y[0], m_minmax_y[1], -1.,
			E_min, m_minmax_z[1], -1.);
	}

	// plot the surfaces
	t_size num_surfs = m_data.size();
	t_size num_active_surfs = 0;
	for(t_size surf_idx = 0; surf_idx < num_surfs; ++surf_idx)
	{
		bool surf_active = surf_idx < active_surfs.size() ? active_surfs[surf_idx] : true;

		// colour for this surface
		std::array<int, 3> col = GetSurfColour(surf_idx, num_surfs);
		const QColor colFull(col[0], col[1], col[2]);

		if(surf_active)
		{
			t_data_vec& data = m_data[surf_idx];

			auto patch_fkt = [this, &data, z_scale](
				t_real_gl /*x2*/, t_real_gl /*x1*/, t_size idx_1, t_size idx_2)
					-> std::pair<t_real_gl, bool>
			{
				t_size idx = idx_1 * m_y_count + idx_2;
				if(idx >= data.size())
					return std::make_pair(0., false);

				bool valid = std::get<4>(data[idx]);
				if(std::get<2>(data[idx]) != idx_1 || std::get<3>(data[idx]) != idx_2)
				{
					std::cerr << "Error: Patch index mismatch: "
						<< "Expected " << std::get<2>(data[idx]) << " for x, but got " << idx_1
						<< "; expected " << std::get<3>(data[idx]) << " for y, but got " << idx_2
						<< "." << std::endl;
					valid = false;
				}

				t_real_gl E = std::get<1>(data[idx]);
				return std::make_pair(E * z_scale, valid);
			};

			std::ostringstream objLabel;
			objLabel << "Surface #" << (surf_idx + 1);

			t_real_gl r = t_real_gl(col[0]) / t_real_gl(255.);
			t_real_gl g = t_real_gl(col[1]) / t_real_gl(255.);
			t_real_gl b = t_real_gl(col[2]) / t_real_gl(255.);

			std::size_t obj = m_dispplot->GetRenderer()->AddPatch(patch_fkt, 0., 0., 0.,
				x_scale, y_scale, m_x_count, m_y_count, r, g, b);
			m_dispplot->GetRenderer()->SetObjectLabel(obj, objLabel.str());
			m_surf_objs.insert(std::make_pair(obj, surf_idx));

			m_cam_centre[2] += GetMeanZ(surf_idx);
			++num_active_surfs;
		}

		AddSurface("#" + tl2::var_to_str(surf_idx + 1), colFull, surf_active);
	}

	if(num_active_surfs)
		m_cam_centre[2] /= static_cast<t_real>(num_active_surfs);

	if(clear_settings)
		CentrePlotCamera();    // centre camera and update
	else
		m_dispplot->update();  // only update renderer
}



/**
 * mouse intersection with surface
 */
void Plot3DDlg::PlotPickerIntersection(
	[[maybe_unused]] const t_vec3_gl* pos,
	[[maybe_unused]] std::size_t objIdx,
	[[maybe_unused]] std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl* posSphere)
{
	m_status->setText("");
	m_cur_obj = std::nullopt;

	m_dispplot->GetRenderer()->SetObjectsHighlight(false);

	if(!pos)
		return;

	m_cur_obj = objIdx;

	auto [ gx, gy, gz ] = PlotToGraphCoords((*pos)[0], (*pos)[1], (*pos)[2]);

	/*auto [ x, y, z ] = GraphToPlotCoords(gx, gy, gz);
	std::cout << x << " " << y << " " << z << "; "
		<< (*pos)[0] << " " << (*pos)[1] << " " << (*pos)[2]
		<< std::endl;*/

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr << "(" << gx << ", " << gy << ", " << gz << ")";

	const std::string& label = m_dispplot->GetRenderer()->GetObjectLabel(objIdx);
	if(label != "")
		ostr << ", " << label;
	ostr << ".";

	m_status->setText(ostr.str().c_str());
}



/**
 * conversion from plot position to graph coordinates
 */
std::tuple<t_real, t_real, t_real> Plot3DDlg::PlotToGraphCoords(t_real x, t_real y, t_real z) const
{
	// coordinate scaling
	t_real x_scale = m_x_scale->value();
	t_real y_scale = m_y_scale->value();
	t_real z_scale = m_z_scale->value();

	// get plot coordinate system
	t_real x_start = m_xrange[0]->value();
	t_real y_start = m_yrange[0]->value();
	t_real x_end = m_xrange[1]->value();
	t_real y_end = m_yrange[1]->value();

	// coordinate trafo
	t_real gx = std::lerp(x_start, x_end, x / x_scale + 0.5);
	t_real gy = std::lerp(y_start, y_end, y / y_scale + 0.5);
	t_real gz = z / z_scale;

	return std::make_tuple(gx, gy, gz);
}



/**
 * conversion from graph coordinates to plot position
 */
std::tuple<t_real, t_real, t_real> Plot3DDlg::GraphToPlotCoords(t_real gx, t_real gy, t_real gz) const
{
	// coordinate scaling
	t_real z_scale = m_z_scale->value();
	t_real x_scale = m_x_scale->value();
	t_real y_scale = m_y_scale->value();

	// get plot coordinate system
	t_real x_start = m_xrange[0]->value();
	t_real y_start = m_yrange[0]->value();
	t_real x_end = m_xrange[1]->value();
	t_real y_end = m_yrange[1]->value();

	// gx = x_start + (x_end - x_start) * (x / x_scale + 0.5)
	// gx - x_start = (x_end - x_start) * (x / x_scale + 0.5)
	// (gx - x_start) / (x_end - x_start) = x / x_scale + 0.5
	// ((gx - x_start) / (x_end - x_start)) * x_scale - 0.5 = x
	t_real x = ((gx - x_start) / (x_end - x_start)) * x_scale - 0.5;
	t_real y = ((gy - y_start) / (y_end - y_start)) * y_scale - 0.5;
	t_real z = gz * z_scale;

	return std::make_tuple(x, y, z);
}
