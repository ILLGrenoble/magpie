/**
 * magnetic dynamics -- 3d dispersion calculation
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

#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <cstdlib>
#include <mutex>
#include <sstream>

#include "dispersion3d.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"


#define COORD_PADDING    1.2     // space between magnon bands and coordinate planes



/**
 * set a pointer to the main magdyn kernel
 */
void Dispersion3DDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * save the Q start and end points from the main window's dispersion
 */
void Dispersion3DDlg::SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend)
{
	m_Qstart = Qstart;
	m_Qend = Qend;
}



/**
 * set the Q position and directions from the main window's Q start and end points
 */
void Dispersion3DDlg::FromMainQ()
{
	if(m_Qstart.size() < 3 || m_Qend.size() < 3 || !m_dyn)
		return;

	const t_mat_real& xtalB = m_dyn->GetCrystalBTrafo();
	const t_vec_real* plane = m_dyn->GetScatteringPlane();
	if(!plane)
		return;

	// direction 1 is from the start to the end point
	t_vec_real Qdir1 = m_Qend - m_Qstart;
	// direction 2 is perpendicular to direction 1 inside the scattering plane
	t_vec_real Qdir2 = tl2::cross(xtalB, plane[2], Qdir1);

	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin[i]->setValue(m_Qstart[i]);
		m_Q_dir1[i]->setValue(Qdir1[i]);
		m_Q_dir2[i]->setValue(Qdir2[i]);
	}
}



/**
 * get the Q origin and direction vectors
 */
std::tuple<t_vec_real, t_vec_real, t_vec_real> Dispersion3DDlg::GetQVectors() const
{
	t_vec_real Q_origin = tl2::create<t_vec_real>(
	{
		m_Q_origin[0]->value(),
		m_Q_origin[1]->value(),
		m_Q_origin[2]->value(),
	});

	t_vec_real Q_dir_1 = tl2::create<t_vec_real>(
	{
		m_Q_dir1[0]->value(),
		m_Q_dir1[1]->value(),
		m_Q_dir1[2]->value(),
	});

	t_vec_real Q_dir_2 = tl2::create<t_vec_real>(
	{
		m_Q_dir2[0]->value(),
		m_Q_dir2[1]->value(),
		m_Q_dir2[2]->value(),
	});

	return std::make_tuple(std::move(Q_origin), std::move(Q_dir_1), std::move(Q_dir_2));
}



/**
 * calculate the dispersion
 */
void Dispersion3DDlg::Calculate()
{
	if(!m_dyn)
		return;

	m_minmax_E[0] = +std::numeric_limits<t_real>::max();
	m_minmax_E[1] = -std::numeric_limits<t_real>::max();
	m_data.clear();

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	// get coordinates
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	m_Q_count_1 = m_num_Q_points[0]->value();
	m_Q_count_2 = m_num_Q_points[1]->value();

	t_vec_real Q_step_1 = Q_dir_1 / t_real(m_Q_count_1);
	t_vec_real Q_step_2 = Q_dir_2 / t_real(m_Q_count_2);

	t_real min_S = m_S_filter->value();
	bool use_weights = m_S_filter_enable->isChecked();
	bool use_projector = true;
	bool unite_degen = m_unite_degeneracies->isChecked();

	// calculate the dispersion
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(unite_degen);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stop_requested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(m_Q_count_1 * m_Q_count_2);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(m_Q_count_1 * m_Q_count_2);

	t_size expected_bands = dyn.GetMagneticSitesCount() * 2;
	if(dyn.IsIncommensurate())
		expected_bands *= 3;
	m_data.resize(expected_bands);

	for(t_size Q_idx_1 = 0; Q_idx_1 < m_Q_count_1; ++Q_idx_1)
	for(t_size Q_idx_2 = 0; Q_idx_2 < m_Q_count_2; ++Q_idx_2)
	{
		auto task = [this, &mtx, &dyn, Q_idx_1, Q_idx_2,
			&Q_origin, &Q_step_1, &Q_step_2, expected_bands,
			unite_degen, use_weights, use_projector, min_S]()
		{
			// calculate the dispersion at the given Q point
			t_vec_real Q = Q_origin + Q_step_1*t_real(Q_idx_1) + Q_step_2*t_real(Q_idx_2);
			auto Es_and_S = dyn.CalcEnergies(Q, !use_weights).E_and_S;

			// iterate the energies for this Q point
			t_size data_band_idx = 0;
			for(t_size band_idx = 0; band_idx < Es_and_S.size() && data_band_idx < expected_bands; ++band_idx, ++data_band_idx)
			{
				const auto& E_and_S = Es_and_S[band_idx];

				bool valid = true;
				t_real E = E_and_S.E;
				if(std::isnan(E) || std::isinf(E))
					valid = false;

				t_real weight = -1;
				if(use_weights)
				{
					weight = E_and_S.weight;

					if(!use_projector)
					{
						const t_mat& S = E_and_S.S;
						weight = tl2::trace<t_mat>(S).real();
					}

					// filter invalid S(Q, E)
					if(std::isnan(weight) || std::isinf(weight))
						weight = 0.;

					// filter minimum S(Q, E)
					if(min_S >= 0. && std::abs(weight) <= min_S)
						valid = false;
				}

				if(valid)
				{
					m_minmax_E[0] = std::min(m_minmax_E[0], E);
					m_minmax_E[1] = std::max(m_minmax_E[1], E);
				}

				// count energy degeneracy
				t_size degeneracy = E_and_S.degeneracy;
				for(t_size band_idx2 = 0; band_idx2 < Es_and_S.size(); ++band_idx2)
				{
					if(band_idx2 == band_idx)
						continue;

					if(tl2::equals(E, Es_and_S[band_idx2].E, g_eps))
						degeneracy += Es_and_S[band_idx2].degeneracy;
				}

				/*if(degeneracy > 1)
				{
					std::cout << "degenerate point: Q indices: " << Q_idx_1 << " " << Q_idx_2
						<< ", band index: " << band_idx << " (" << data_band_idx
						<< "): " << E << " meV (" << degeneracy << "x)" << std::endl;
				}*/

				// generate and add data point
				t_data_Q dat{std::make_tuple(Q, E, weight, Q_idx_1, Q_idx_2, degeneracy, valid)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[data_band_idx].emplace_back(std::move(dat));

				if(unite_degen && degeneracy > 1)
				{
					// skip degeneracies to keep matching bands together
					for(t_size band_idx2 = data_band_idx + 1; band_idx2 < data_band_idx + degeneracy; ++band_idx2)
						m_data[band_idx2].emplace_back(std::make_tuple(Q, 0., 0., Q_idx_1, Q_idx_2, 1, false));

					data_band_idx += degeneracy - 1;
				}
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
			/*
			// sorting by Q components
			t_real h1 = std::get<0>(m_data[0][idx1])[0];
			t_real k1 = std::get<0>(m_data[0][idx1])[1];
			t_real l1 = std::get<0>(m_data[0][idx1])[2];

			t_real h2 = std::get<0>(m_data[0][idx2])[0];
			t_real k2 = std::get<0>(m_data[0][idx2])[1];
			t_real l2 = std::get<0>(m_data[0][idx2])[2];

			if(!tl2::equals(h1, h2, g_eps))
				return h1 < h2;
			if(!tl2::equals(k1, k2, g_eps))
				return k1 < k2;

			return l1 < l2;*/

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

	// move degenerate points to the bands where most of the other points are
	if(unite_degen && m_data.size())
	{
		for(t_size band_idx = 0; band_idx < m_data.size() - 1; ++band_idx)
		{
			t_data_Qs& band1 = m_data[band_idx];
			t_data_Qs& band2 = m_data[band_idx + 1];
			if(band1.size() != band2.size())
				continue;

			auto [_valid1, invalid1] = NumValid(band1);
			auto [_valid2, invalid2] = NumValid(band2);

			for(t_size Q_idx = 0; Q_idx < band1.size(); ++Q_idx)
			{
				t_size degen1 = std::get<5>(band1[Q_idx]);
				bool valid1 = std::get<6>(band1[Q_idx]);
				bool valid2 = std::get<6>(band2[Q_idx]);

				if(degen1 > 1 && valid1 && !valid2 && invalid1 > invalid2)
					std::swap(band1[Q_idx], band2[Q_idx]);
			}
		}
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
		return GetMeanEnergy(dat1) >= GetMeanEnergy(dat2);
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
 * count number of bands with positive energies (should be half of all bands)
 */
t_size Dispersion3DDlg::NumPositive() const
{
	t_size num_pos = 0;

	for(const t_data_Qs& band : m_data)
	{
		if(IsPositive(band))
			++num_pos;
	}

	return num_pos;
}



/**
 * determine if a band only contains positive energies
 */
bool Dispersion3DDlg::IsPositive(const t_data_Qs& data) const
{
	for(const t_data_Q& data : data)
	{
		t_real E = std::get<1>(data);
		if(E < 0.)
			return false;
	}

	return true;
}



/**
 * determine if a band is valid or only contains invalid points
 */
bool Dispersion3DDlg::IsValid(const t_data_Qs& data) const
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
std::pair<t_size, t_size> Dispersion3DDlg::NumValid(const t_data_Qs& data) const
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
 * calculate the mean band energy
 */
t_real Dispersion3DDlg::GetMeanEnergy(const Dispersion3DDlg::t_data_Qs& data) const
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
 * calculate the mean band energy
 */
t_real Dispersion3DDlg::GetMeanEnergy(t_size band_idx) const
{
	if(band_idx >= m_data.size())
		return 0.;

	return GetMeanEnergy(m_data[band_idx]);
}



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
	if(!m_dispplot)
		return;

	// keep some settings from previous plot, e.g. the band visibility flags
	std::vector<bool> active_bands;
	if(!clear_settings)
	{
		active_bands.reserve(m_table_bands->rowCount());
		for(int row = 0; row < m_table_bands->rowCount(); ++row)
			active_bands.push_back(IsBandEnabled(t_size(row)));
	}

	bool only_pos_E = m_only_pos_E->isChecked();
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
		t_real E_range = m_minmax_E[1];

		if(!only_pos_E)
			E_range -= m_minmax_E[0];

		t_real E_mean = only_pos_E ? 0.5*E_range : m_minmax_E[0] + 0.5*E_range;

		t_mat_gl obj_scale = tl2::hom_scaling<t_mat_gl>(
			0.5 * COORD_PADDING * Q_scale1, 0.5 * COORD_PADDING * Q_scale2,
			0.5 * COORD_PADDING * E_scale * E_range);
		t_mat_gl obj_shift = tl2::hom_translation<t_mat_gl>(0., 0., E_scale * E_mean);

		for(auto obj : m_dispplot->GetRenderer()->GetCoordCube())
		{
			using namespace tl2_ops;
			m_dispplot->GetRenderer()->SetObjectMatrix(obj, obj_shift * obj_scale);
		}

		// TODO
		m_dispplot->GetRenderer()->UpdateCoordCubeTextures();
	}

	// plot the magnon bands
	t_size num_bands = m_data.size();
	if(only_pos_E)
	{
		num_bands = NumPositive();
		//num_bands /= 2;
	}
	t_size num_active_bands = 0;
	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		bool band_active = band_idx < active_bands.size() ? active_bands[band_idx] : true;

		// colour for this magnon band
		std::array<int, 3> col = GetBranchColour(band_idx, num_bands);
		const QColor colFull(col[0], col[1], col[2]);

		if(band_active)
		{
			t_data_Qs& data = m_data[band_idx];

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

			std::ostringstream objLabel;
			objLabel << "Band #" << (band_idx + 1);

			t_real_gl r = t_real_gl(col[0]) / t_real_gl(255.);
			t_real_gl g = t_real_gl(col[1]) / t_real_gl(255.);
			t_real_gl b = t_real_gl(col[2]) / t_real_gl(255.);

			std::size_t obj = m_dispplot->GetRenderer()->AddPatch(patch_fkt, 0., 0., 0.,
				Q_scale1, Q_scale2, m_Q_count_1, m_Q_count_2, r, g, b);
			m_dispplot->GetRenderer()->SetObjectLabel(obj, objLabel.str());
			m_band_objs.insert(std::make_pair(obj, band_idx));

			m_cam_centre[2] += GetMeanEnergy(band_idx);
			++num_active_bands;
		}

		AddBand("#" + tl2::var_to_str(band_idx + 1), colFull, band_active);
	}

	if(num_active_bands)
		m_cam_centre[2] /= static_cast<t_real>(num_active_bands);

	if(clear_settings)
		CentrePlotCamera();    // centre camera and update
	else
		m_dispplot->update();  // only update renderer
}



/**
 * clears the table of magnon bands
 */
void Dispersion3DDlg::ClearBands()
{
	m_band_objs.clear();
	m_cur_obj = std::nullopt;

	m_table_bands->clearContents();
	m_table_bands->setRowCount(0);
}



/**
 * adds a magnon band to the table
 */
void Dispersion3DDlg::AddBand(const std::string& name, const QColor& colour, bool enabled)
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
 * verifies if the band's checkbox is checked
 */
bool Dispersion3DDlg::IsBandEnabled(t_size idx) const
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
void Dispersion3DDlg::PlotPickerIntersection(
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
	t_real E_scale = m_E_scale->value();
	t_real Q_scale1 = m_Q_scale1->value();
	t_real Q_scale2 = m_Q_scale2->value();

	// get Q and E position at cursor intersection
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	// reconstruct Q position
	t_real_gl Q1param = (0.5*Q_scale1 + (*pos)[0]) / Q_scale1;
	t_real_gl Q2param = (0.5*Q_scale2 + (*pos)[1]) / Q_scale2;
	t_vec_real Q = Q_origin + Q_dir_1*Q1param + Q_dir_2*Q2param;

	t_real_gl E = (*pos)[2] / E_scale;

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr
		<< "Q = (" << Q[0] << ", " << Q[1] << ", " << Q[2] << ") rlu, "
		<< "E = " << E << " meV";

	const std::string& label = m_dispplot->GetRenderer()->GetObjectLabel(objIdx);
	if(label != "")
		ostr << ", " << label;
	ostr << ".";

	m_status->setText(ostr.str().c_str());
}
