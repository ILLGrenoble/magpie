/**
 * 3d plotter calculation
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

#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <cstdlib>
#include <mutex>
#include <sstream>

#include "plot3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * calculate the dispersion
 */
void Plot3DDlg::Calculate()
{
	m_minmax_z[0] = +std::numeric_limits<t_real>::max();
	m_minmax_z[1] = -std::numeric_limits<t_real>::max();
	m_minmax_x[0] = m_minmax_x[1] = tl2::zero<t_vec>(3);
	m_minmax_y[0] = m_minmax_y[1] = tl2::zero<t_vec>(3);
	m_data.clear();

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	m_x_count = m_num_points[0]->value();
	m_y_count = m_num_points[1]->value();

	// TODO
	t_vec Q_origin = tl2::zero<t_vec>(3);
	t_vec Q_dir_1 = tl2::zero<t_vec>(3);
	t_vec Q_dir_2 = tl2::zero<t_vec>(3);

	t_vec Q_step_1 = Q_dir_1 / t_real(m_x_count);
	t_vec Q_step_2 = Q_dir_2 / t_real(m_y_count);

	m_minmax_x[0] = Q_origin;
	m_minmax_x[1] = Q_origin + Q_step_1*t_real(m_x_count - 1);
	m_minmax_y[0] = Q_origin;
	m_minmax_y[1] = Q_origin + Q_step_2*t_real(m_y_count - 1);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stop_requested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(m_x_count * m_y_count);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(m_x_count * m_y_count);

	t_size expected_bands = 1;
	m_data.resize(expected_bands);

	for(t_size Q_idx_1 = 0; Q_idx_1 < m_x_count; ++Q_idx_1)
	for(t_size Q_idx_2 = 0; Q_idx_2 < m_y_count; ++Q_idx_2)
	{
		auto task = [this, &mtx, Q_idx_1, Q_idx_2,
			&Q_origin, &Q_step_1, &Q_step_2, expected_bands]()
		{
			// calculate the dispersion at the given Q point
			t_vec Q = Q_origin + Q_step_1*t_real(Q_idx_1) + Q_step_2*t_real(Q_idx_2);

			// iterate the energies for this Q point
			t_size data_band_idx = 0;
			for(t_size band_idx = 0; band_idx < /*Es_and_S.size()*/ 1 && data_band_idx < expected_bands; ++band_idx, ++data_band_idx)
			{
				bool valid = true;
				t_real E = 0.;
				if(std::isnan(E) || std::isinf(E))
					valid = false;

				t_real weight = -1;

				if(valid)
				{
					m_minmax_z[0] = std::min(m_minmax_z[0], E);
					m_minmax_z[1] = std::max(m_minmax_z[1], E);
				}

				// count energy degeneracy
				t_size degeneracy = 0.;

				// generate and add data point
				t_data_Q dat{std::make_tuple(Q, E, weight, Q_idx_1, Q_idx_2, degeneracy, valid)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[data_band_idx].emplace_back(std::move(dat));
			}

			// fill up band data in case some indices were skipped due to invalid hamiltonians
			for(; data_band_idx < expected_bands; ++data_band_idx)
			{
				t_data_Q dat{std::make_tuple(Q, 0., 0., Q_idx_1, Q_idx_2, 1, false)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[data_band_idx].emplace_back(std::move(dat));
			}
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating dispersion in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		// process events to see if the stop button was clicked
		// only do this for a fraction of the points to avoid gui overhead
		if(task_idx % std::max<t_size>(tasks.size() / std::sqrt(g_stop_check_fraction), 1) == 0)
			qApp->processEvents();

		if(m_stop_requested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	// finish parallel calculations
	pool.join();

	// get sorting of data by Q
	for(t_size band_idx = 0; band_idx < m_data.size(); ++band_idx)
	{
		std::vector<std::size_t> perm = tl2::get_perm(m_data[band_idx].size(),
			[this, band_idx](std::size_t idx1, std::size_t idx2) -> bool
		{
			// sorting by Q indices
			t_size Q1_idx_1 = std::get<3>(m_data[band_idx][idx1]);
			t_size Q1_idx_2 = std::get<4>(m_data[band_idx][idx1]);
			t_size Q2_idx_1 = std::get<3>(m_data[band_idx][idx2]);
			t_size Q2_idx_2 = std::get<4>(m_data[band_idx][idx2]);

			if(Q1_idx_1 != Q2_idx_1)
				return Q1_idx_1 < Q2_idx_1;
			return Q1_idx_2 < Q2_idx_2;
		});

		m_data[band_idx] = tl2::reorder(m_data[band_idx], perm);
	}

	// remove fully invalid bands
	for(auto iter = m_data.begin(); iter != m_data.end();)
	{
		if(!IsValid(*iter))
			iter = m_data.erase(iter);
		else
			++iter;
	}

	// sort band energies in descending order
	std::stable_sort(m_data.begin(), m_data.end(), [this](const t_data_Qs& dat1, const t_data_Qs& dat2)
	{
		return GetMeanZ(dat1) >= GetMeanZ(dat2);
	});

	// show elapsed time
	stopwatch.stop();
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Dispersion calculation";
	if(m_stop_requested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	Plot(true);
}



/**
 * determine if a band is valid or only contains invalid points
 */
bool Plot3DDlg::IsValid(const t_data_Qs& data) const
{
	for(const t_data_Q& data : data)
	{
		// data point valid?
		if(std::get<6>(data))
			return true;
	}

	// all points invalid
	return false;
}



/**
 * count the number of valid and invalid points in a band
 */
std::pair<t_size, t_size> Plot3DDlg::NumValid(const t_data_Qs& data) const
{
	t_size valid = 0, invalid = 0;

	for(const t_data_Q& data : data)
	{
		// data point valid?
		if(std::get<6>(data))
			++valid;
		else
			++invalid;
	}

	return std::make_pair(valid, invalid);
}



/**
 * calculate the mean surface z value
 */
t_real Plot3DDlg::GetMeanZ(const Plot3DDlg::t_data_Qs& data) const
{
	t_real E_mean = 0.;
	t_size num_pts = 0;

	for(const t_data_Q& data : data)
	{
		// data point invalid?
		if(!std::get<6>(data))
			continue;

		E_mean += std::get<1>(data);
		++num_pts;
	}

	if(num_pts)
		E_mean /= static_cast<t_real>(num_pts);
	return E_mean;
}



/**
 * calculate the mean surface z value
 */
t_real Plot3DDlg::GetMeanZ(t_size band_idx) const
{
	if(band_idx >= m_data.size())
		return 0.;

	return GetMeanZ(m_data[band_idx]);
}



/**
 * get a unique colour for a given surface
 */
std::array<int, 3> Plot3DDlg::GetBranchColour(t_size branch_idx, t_size num_branches) const
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
void Plot3DDlg::Plot(bool clear_settings)
{
	if(!m_dispplot)
		return;

	// keep some settings from previous plot, e.g. the band visibility flags
	std::vector<bool> active_bands;
	if(!clear_settings)
	{
		active_bands.reserve(m_table_bands->rowCount());
		for(int row = 0; row < m_table_bands->rowCount(); ++row)
			active_bands.push_back(IsSurfaceEnabled(t_size(row)));
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
		t_mat_gl obj_shift = tl2::hom_translation<t_mat_gl>(0., 0., z_scale * E_mean);

		for(auto obj : m_dispplot->GetRenderer()->GetCoordCube())
		{
			using namespace tl2_ops;
			m_dispplot->GetRenderer()->SetObjectMatrix(obj, obj_shift * obj_scale);
		}

		if(m_dispplot || m_minmax_x[0].size() == 3)
		{
			m_dispplot->GetRenderer()->UpdateCoordCubeTextures(
				m_minmax_x[0][0], m_minmax_x[1][0], -1.,
				m_minmax_y[0][1], m_minmax_y[1][1], -1.,
				E_min, m_minmax_z[1], -1.);
		}
	}

	// plot the surfaces
	t_size num_bands = m_data.size();
	t_size num_active_bands = 0;
	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		bool band_active = band_idx < active_bands.size() ? active_bands[band_idx] : true;

		// colour for this surface
		std::array<int, 3> col = GetBranchColour(band_idx, num_bands);
		const QColor colFull(col[0], col[1], col[2]);

		if(band_active)
		{
			t_data_Qs& data = m_data[band_idx];

			auto patch_fkt = [this, &data, z_scale](
				t_real_gl /*x2*/, t_real_gl /*x1*/, t_size idx_1, t_size idx_2)
					-> std::pair<t_real_gl, bool>
			{
				t_size idx = idx_1 * m_y_count + idx_2;
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
				return std::make_pair(E * z_scale, valid);
			};

			std::ostringstream objLabel;
			objLabel << "Surface #" << (band_idx + 1);

			t_real_gl r = t_real_gl(col[0]) / t_real_gl(255.);
			t_real_gl g = t_real_gl(col[1]) / t_real_gl(255.);
			t_real_gl b = t_real_gl(col[2]) / t_real_gl(255.);

			std::size_t obj = m_dispplot->GetRenderer()->AddPatch(patch_fkt, 0., 0., 0.,
				x_scale, y_scale, m_x_count, m_y_count, r, g, b);
			m_dispplot->GetRenderer()->SetObjectLabel(obj, objLabel.str());
			m_band_objs.insert(std::make_pair(obj, band_idx));

			m_cam_centre[2] += GetMeanZ(band_idx);
			++num_active_bands;
		}

		AddSurface("#" + tl2::var_to_str(band_idx + 1), colFull, band_active);
	}

	if(num_active_bands)
		m_cam_centre[2] /= static_cast<t_real>(num_active_bands);

	if(clear_settings)
		CentrePlotCamera();    // centre camera and update
	else
		m_dispplot->update();  // only update renderer
}



/**
 * clears the table of surfaces
 */
void Plot3DDlg::ClearSurfaces()
{
	m_band_objs.clear();
	m_cur_obj = std::nullopt;

	m_table_bands->clearContents();
	m_table_bands->setRowCount(0);
}



/**
 * adds a surface to the table
 */
void Plot3DDlg::AddSurface(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_bands)
		return;

	int row = m_table_bands->rowCount();
	m_table_bands->insertRow(row);

	QTableWidgetItem *item = new QTableWidgetItem{name.c_str()};
	item->setFlags(item->flags() & ~Qt::ItemIsEditable);

	QBrush bg = item->background();
	bg.setColor(colour);
	bg.setStyle(Qt::SolidPattern);
	item->setBackground(bg);

	QBrush fg = item->foreground();
	fg.setColor(QColor{0xff, 0xff, 0xff});
	fg.setStyle(Qt::SolidPattern);
	item->setForeground(fg);

	QCheckBox *checkBand = new QCheckBox(m_table_bands);
	checkBand->setChecked(enabled);
	connect(checkBand, &QCheckBox::toggled, [this]() { Plot(false); });

	m_table_bands->setItem(row, COL_BC_BAND, item);
	m_table_bands->setCellWidget(row, COL_BC_ACTIVE, checkBand);
}



/**
 * verifies if the surface's checkbox is checked
 */
bool Plot3DDlg::IsSurfaceEnabled(t_size idx) const
{
	if(!m_table_bands || int(idx) >= m_table_bands->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(
		m_table_bands->cellWidget(int(idx), COL_BC_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * mouse intersection with dispersion band
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

	// coordinate scaling
	t_real x_scale = m_x_scale->value();
	t_real y_scale = m_y_scale->value();
	t_real z_scale = m_z_scale->value();

	t_real_gl x = (*pos)[0] / x_scale;
	t_real_gl y = (*pos)[1] / y_scale;
	t_real_gl z = (*pos)[2] / z_scale;

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr << "(" << x << ", " << y << ", " << z << ")";

	const std::string& label = m_dispplot->GetRenderer()->GetObjectLabel(objIdx);
	if(label != "")
		ostr << ", " << label;
	ostr << ".";

	m_status->setText(ostr.str().c_str());
}
