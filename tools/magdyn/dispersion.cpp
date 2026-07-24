/**
 * magnetic dynamics -- calculations for dispersion plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * magpie & mag-core
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

// these need to be included before all other things on mingw
#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include "magdyn.h"

#include <QtWidgets/QApplication>

#include <future>
#include <mutex>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"

using namespace tl2_ops;



/**
 * clears the dispersion graph
 */
void MagDynDlg::ClearDispersion(bool replot)
{
	m_graphs.clear();

	if(m_plot)
	{
		m_plot->clearPlottables();
		if(replot)
			m_plot->replot();
	}

	m_qs_data.clear();
	m_Es_data.clear();
	m_ws_data.clear();
	m_degen_data.clear();

	for(int i = 0; i < 2*3*3; ++i)
	{
		m_qs_data_channel[i].clear();
		m_Es_data_channel[i].clear();
		m_ws_data_channel[i].clear();
		m_ws_total_channel[i] = 0.;
	}

	m_Q_idx = 0;
}



/**
 * get the dispersion's start and end points in Q
 */
std::pair<t_vec_real, t_vec_real> MagDynDlg::GetDispersionQ() const
{
	t_vec_real Q_start = tl2::create<t_vec_real>(
	{
		(t_real)m_Q_start[0]->value(),
		(t_real)m_Q_start[1]->value(),
		(t_real)m_Q_start[2]->value(),
	});

	t_vec_real Q_end = tl2::create<t_vec_real>(
	{
		(t_real)m_Q_end[0]->value(),
		(t_real)m_Q_end[1]->value(),
		(t_real)m_Q_end[2]->value(),
	});

	return std::make_pair(std::move(Q_start), std::move(Q_end));
}



/**
 * get the dispersion's start and end points in E
 */
std::pair<t_real, t_real> MagDynDlg::GetDispersionE() const
{
	t_real E_min = m_E_min;
	if(E_min < 0. && m_ignore_annihilation->isChecked())
		E_min = 0.;

	return std::make_pair(E_min, m_E_max);
}



/**
 * a new start or end Q coordinate has been entered
 */
void MagDynDlg::DispersionQChanged(bool calc_dyn)
{
	if(this->m_autocalc->isChecked() && calc_dyn)
		this->CalcDispersion();

	t_vec_real Q_start, Q_end;
	t_real E_start{0}, E_end{1};
	if(m_topo_dlg || m_diff_dlg || m_powder_dlg || m_disp3d_dlg || m_bzscene || m_bz_dlg)
	{
		std::tie(Q_start, Q_end) = GetDispersionQ();
		std::tie(E_start, E_end) = GetDispersionE();
	}

	if(m_topo_dlg)
		m_topo_dlg->SetDispersionQ(Q_start, Q_end);
	if(m_diff_dlg)
		m_diff_dlg->SetDispersionQ(Q_start, Q_end);
	if(m_powder_dlg)
		m_powder_dlg->SetDispersionQE(Q_start, Q_end, E_start, E_end);
	if(m_disp3d_dlg)
		m_disp3d_dlg->SetDispersionQ(Q_start, Q_end);

	// show current scan in scattering plane plot
	if(m_bzscene)
	{
		// draw scan line if inside scattering plane
		m_bzscene->ClearLines();

		t_mat_real UB = m_dyn.GetCrystalUBTrafo();
		t_vec_real pt_start = UB*Q_start;
		t_vec_real pt_end = UB*Q_end;

		// add line if both points are in the scattering plane
		if(tl2::equals_0<t_real>(pt_start[2], g_eps) &&
			tl2::equals_0<t_real>(pt_end[2], g_eps))
			m_bzscene->AddLine(pt_start, pt_end, true,
			"Scan Direction", "Start Q", "End Q");
	}

	// show current scan in bz plot
	if(m_bz_dlg)
	{
		// draw scan line
		m_bz_dlg->ClearLines(false);

		t_mat_real B = m_dyn.GetCrystalBTrafo();
		t_vec_real pt_start = B*Q_start;
		t_vec_real pt_end = B*Q_end;

		m_bz_dlg->AddLine(pt_start, pt_end, false);
	}
}



/**
 * calculate the dispersion branches
 */
void MagDynDlg::CalcDispersion()
{
	if(m_ignoreCalc)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableInput(true);
	} BOOST_SCOPE_EXIT_END
	EnableInput(false);

	// nothing to calculate?
	if(m_dyn.GetMagneticSitesCount() == 0 || m_dyn.GetExchangeTermsCount() == 0)
	{
		ClearDispersion(true);
		return;
	}
	else
	{
		ClearDispersion(false);
	}


	// get Qs
	auto [Q_start, Q_end] = GetDispersionQ();
	t_vec_real Q_range = Q_end - Q_start;
	for(int i = 0; i < 3; ++i)
		Q_range[i] =  std::abs(Q_range[i]);

	// Q component with maximum range
	m_Q_idx = 0;
	if(Q_range[1] > Q_range[m_Q_idx])
		m_Q_idx = 1;
	if(Q_range[2] > Q_range[m_Q_idx])
		m_Q_idx = 2;

	m_Q_min = Q_start[m_Q_idx];
	m_Q_max = Q_end[m_Q_idx];

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx] > Q_end[m_Q_idx])
		std::swap(Q_start, Q_end);


	// reserve vector memory
	t_size num_pts = m_num_points->value();

	m_qs_data.reserve(num_pts*10);
	m_Es_data.reserve(num_pts*10);
	m_ws_data.reserve(num_pts*10);
	m_degen_data.reserve(num_pts*10);

	for(t_size i = 0; i < 2*3*3; ++i)
	{
		m_qs_data_channel[i].reserve(num_pts*10);
		m_Es_data_channel[i].reserve(num_pts*10);
		m_ws_data_channel[i].reserve(num_pts*10);
	}


	// options
	const bool is_comm = !m_dyn.IsIncommensurate();
	const bool use_min_E = false;
	const bool unite_degeneracies = m_unite_degeneracies->isChecked();
	const bool ignore_annihilation = m_ignore_annihilation->isChecked();
	const bool use_weights = m_use_weights->isChecked();
	const bool use_projector = m_use_projector->isChecked();
	const bool use_polcoords = m_use_polcoords->isChecked();;
	const bool force_incommensurate = m_force_incommensurate->isChecked();

	t_real E0 = use_min_E ? m_dyn.CalcMinimumEnergy() : 0.;
	m_dyn.SetUniteDegenerateEnergies(unite_degeneracies);
	m_dyn.SetForceIncommensurate(force_incommensurate);
	m_dyn.SetCalcHamiltonian(
		m_hamiltonian_comp[0]->isChecked() || is_comm,  // always calculate commensurate case
		m_hamiltonian_comp[1]->isChecked(),
		m_hamiltonian_comp[2]->isChecked());

	// thread pool and mutex to protect m_qs_data, m_Es_data, and m_ws_data
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(num_pts);
	m_progress->setValue(0);
	m_status->setText(QString("Starting dispersion calculation using %1 threads.").arg(g_num_threads));
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(num_pts);

	for(t_size i = 0; i < num_pts; ++i)
	{
		auto task = [this, &mtx, i, num_pts, E0, &Q_start, &Q_end,
			use_projector, use_weights, use_polcoords, ignore_annihilation]()
		{
			const t_vec_real Q = num_pts > 1
				? tl2::create<t_vec_real>(
				{
					std::lerp(Q_start[0], Q_end[0], t_real(i) / t_real(num_pts - 1)),
					std::lerp(Q_start[1], Q_end[1], t_real(i) / t_real(num_pts - 1)),
					std::lerp(Q_start[2], Q_end[2], t_real(i) / t_real(num_pts - 1)),
				})
				: tl2::create<t_vec_real>({ Q_start[0], Q_start[1], Q_start[2] });

			auto S = m_dyn.CalcEnergies(Q, !use_weights);

			for(const auto& E_and_S : S.E_and_S)
			{
				if(m_stopRequested)
					break;

				t_real E = E_and_S.E - E0;
				if(std::isnan(E) || std::isinf(E))
					continue;
				if(ignore_annihilation && E < t_real(0))
					continue;

				std::lock_guard<std::mutex> _lck{mtx};

				// weights
				if(use_weights)
				{
					const t_mat *S = &E_and_S.S_perp;
					t_real weight = 0.;

					if(use_polcoords)
					{
						if(use_projector)
						{
							S = &E_and_S.S_pol_perp;
							weight = E_and_S.weight_pol_perp;
						}
						else
						{
							S = &E_and_S.S_pol;
							weight = E_and_S.weight_pol_full;
						}
					}
					else
					{
						if(use_projector)
						{
							S = &E_and_S.S_perp;
							weight = E_and_S.weight_perp;
						}
						else
						{
							S = &E_and_S.S;
							weight = E_and_S.weight_full;
						}
					}

					if(std::isnan(weight) || std::isinf(weight))
						continue;

					m_ws_data.push_back(weight);

					for(t_size channel_i = 0; channel_i < 3; ++channel_i)
					for(t_size channel_j = 0; channel_j < 3; ++channel_j)
					{
						t_real weight_channel_re = 0., weight_channel_im = 0.;
						if(S->size1() > channel_i && S->size2() > channel_j)
						{
							weight_channel_re = std::abs((*S)(channel_i, channel_j).real());
							weight_channel_im = std::abs((*S)(channel_i, channel_j).imag());
						}

						if(!tl2::equals_0<t_real>(weight_channel_re, g_eps))
						{
							m_qs_data_channel[channel_i*3 + channel_j].push_back(Q[m_Q_idx]);
							m_Es_data_channel[channel_i*3 + channel_j].push_back(E);
							m_ws_data_channel[channel_i*3 + channel_j].push_back(weight_channel_re);
						}
						if(!tl2::equals_0<t_real>(weight_channel_im, g_eps))
						{
							m_qs_data_channel[channel_i*3 + channel_j + 3*3].push_back(Q[m_Q_idx]);
							m_Es_data_channel[channel_i*3 + channel_j + 3*3].push_back(E);
							m_ws_data_channel[channel_i*3 + channel_j + 3*3].push_back(weight_channel_im);
						}
					}  // channel
				}  // weights

				m_qs_data.push_back(Q[m_Q_idx]);
				m_Es_data.push_back(E);
				m_degen_data.push_back(static_cast<int>(E_and_S.degeneracy - 1));
			}  // E_and_S
		};  // task

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}  // num_pts

	m_status->setText(QString("Calculating dispersion in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		// process events to see if the stop button was clicked
		// only do this for a fraction of the points to avoid gui overhead
		bool process_evts = (task_idx % std::max<t_size>(tasks.size() / g_stop_check_fraction, 1) == 0);
		if(process_evts)
			qApp->processEvents();

		if(m_stopRequested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();

		if(process_evts || task_idx + 1 == tasks.size())
			m_progress->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	auto sort_data = [](QVector<qreal>& qvec,
		QVector<qreal>& Evec,
		QVector<qreal>& wvec,
		QVector<int>* dvec = nullptr) -> std::vector<std::size_t>
	{
		// sort vectors by q component
		std::vector<std::size_t> perm = tl2::get_perm(qvec.size(),
			[&qvec, &Evec](std::size_t idx1, std::size_t idx2) -> bool
		{
			// if q components are equal, sort by E
			if(tl2::equals(qvec[idx1], qvec[idx2], (qreal)g_eps))
				return Evec[idx1] < Evec[idx2];
			return qvec[idx1] < qvec[idx2];
		});

		qvec = tl2::reorder(qvec, perm);
		Evec = tl2::reorder(Evec, perm);
		wvec = tl2::reorder(wvec, perm);

		if(dvec)
			*dvec = tl2::reorder(*dvec, perm);

		return perm;
	};

	//std::vector<std::size_t> perm =
	sort_data(m_qs_data, m_Es_data, m_ws_data, &m_degen_data);
	for(t_size i = 0; i < 2*3*3; ++i)
	{
		//tl2::reorder(m_ws_data_channel[i], perm);
		sort_data(m_qs_data_channel[i], m_Es_data_channel[i], m_ws_data_channel[i]);

		// total weight in channel
		m_ws_total_channel[i] = std::accumulate(
			m_ws_data_channel[i].begin(), m_ws_data_channel[i].end(), 0.);
	}

	PlotDispersion();
}



/**
 * set the number of Q points on the dispersion to calculate
 */
void MagDynDlg::SetNumQPoints(t_size num_Q_pts)
{
	m_num_points->setValue((int)num_Q_pts);
}



/**
 * set the current dispersion path and the hamiltonian to the given one
 */
void MagDynDlg::SetCoordinates(const t_vec_real& Qi, const t_vec_real& Qf, bool calc_dynamics)
{
	m_ignoreCalc = true;

	BOOST_SCOPE_EXIT(this_, calc_dynamics)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked() && calc_dynamics)
		{
			this_->CalcDispersion();
			this_->CalcHamiltonian();
		}
	} BOOST_SCOPE_EXIT_END

	const bool set_Qi = (Qi.size() >= 3);
	const bool set_Qf = (Qf.size() >= 3);

	// calculate the dispersion from Qi to Qf
	if(set_Qi)
	{
		m_Q_start[0]->setValue(Qi[0]);
		m_Q_start[1]->setValue(Qi[1]);
		m_Q_start[2]->setValue(Qi[2]);

		// calculate the hamiltonian for Qi
		m_Q[0]->setValue(Qi[0]);
		m_Q[1]->setValue(Qi[1]);
		m_Q[2]->setValue(Qi[2]);
	}

	if(set_Qf)
	{
		m_Q_end[0]->setValue(Qf[0]);
		m_Q_end[1]->setValue(Qf[1]);
		m_Q_end[2]->setValue(Qf[2]);
	}
}



/**
 * set the selected coordinate path as the current one
 */
void MagDynDlg::SetCurrentCoordinate(int which)
{
	using t_item = tl2::NumericTableWidgetItem<t_real>;

	int idx_i = m_coordinates_cursor_row;
	if(idx_i < 0 || idx_i >= m_coordinatestab->rowCount())
		return;

	const auto* hi = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_H));
	const auto* ki = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_K));
	const auto* li = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_L));

	// set dispersion start and end coordinates
	if(which == 0)
	{
		int idx_f = idx_i + 1;

		// wrap around
		if(idx_f >= m_coordinatestab->rowCount())
			idx_f = 0;

		if(idx_f == idx_i)
			return;
		if(idx_f < 0 || idx_f >= m_coordinatestab->rowCount())
			return;

		const auto* hf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_H));
		const auto* kf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_K));
		const auto* lf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_L));

		if(!hi || !ki || !li || !hf || !kf || !lf)
			return;

		m_ignoreCalc = true;

		BOOST_SCOPE_EXIT(this_)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
				this_->CalcDispersion();
		} BOOST_SCOPE_EXIT_END

		m_Q_start[0]->setValue(hi->GetValue());
		m_Q_start[1]->setValue(ki->GetValue());
		m_Q_start[2]->setValue(li->GetValue());
		m_Q_end[0]->setValue(hf->GetValue());
		m_Q_end[1]->setValue(kf->GetValue());
		m_Q_end[2]->setValue(lf->GetValue());
	}

	// send initial Q coordinates to hamiltonian calculation
	else if(which == 1)
	{
		if(!hi || !ki || !li)
			return;

		m_ignoreCalc = true;

		BOOST_SCOPE_EXIT(this_)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
				this_->CalcHamiltonian();
		} BOOST_SCOPE_EXIT_END

		m_Q[0]->setValue(hi->GetValue());
		m_Q[1]->setValue(ki->GetValue());
		m_Q[2]->setValue(li->GetValue());
	}
}
