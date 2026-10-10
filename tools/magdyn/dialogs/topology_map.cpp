/**
 * magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
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

#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <limits>
#include <mutex>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QFileDialog>

#include "topology.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"
#include "libs/helpers.h"


// ============================================================================
// calculate berry curvature map
// ============================================================================

/**
 * create the panel for the berry curvature tab
 */
QWidget* TopologyDlg::CreateBerryCurvatureMapPanel()
{
	QWidget *panel = new QWidget(this);

	// plotter
	m_plot_bcm = new QCustomPlot(panel);
	m_plot_bcm->setFont(font());
	m_plot_bcm->xAxis->setLabel("Momentum Transfer Q1 (rlu)");
	m_plot_bcm->yAxis->setLabel("Momentum Transfer Q1 (rlu)");
	m_plot_bcm->setInteraction(QCP::iRangeDrag, true);
	m_plot_bcm->setInteraction(QCP::iRangeZoom, true);
	m_plot_bcm->setSelectionRectMode(QCP::srmZoom);
	m_plot_bcm->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	m_plot_colour_bcm = new QCPColorScale(m_plot_bcm);
	m_plot_colour_bcm->axis()->setScaleType(QCPAxis::stLinear);
	m_plot_colour_bcm->axis()->setLabel("Berry Curvature B");
	m_plot_colour_bcm->setRangeDrag(true);
	m_plot_colour_bcm->setRangeZoom(true);
	m_plot_colour_bcm->setType(QCPAxis::atRight);
	m_plot_bcm->plotLayout()->addElement(0, 1, m_plot_colour_bcm);

	m_plot_map_bcm = new QCPColorMap(m_plot_bcm->xAxis, m_plot_bcm->yAxis);
	m_plot_map_bcm->setTightBoundary(true);
	m_plot_map_bcm->setColorScale(m_plot_colour_bcm);
	m_plot_map_bcm->setGradient(QCPColorGradient::gpHot);
	m_plot_map_bcm->setInterpolate(true);
	m_plot_map_bcm->setDataScaleType(QCPAxis::stLinear);

	// context menu for plotter
	m_menuPlot_bcm = new QMenu("Plotter", panel);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot_bcm);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot_bcm);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot_bcm);

	if(g_use_icons)
	{
		acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
		acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	}

	m_menuPlot_bcm->addAction(acRescalePlot);
	m_menuPlot_bcm->addSeparator();
	m_menuPlot_bcm->addAction(acSaveFigure);
	m_menuPlot_bcm->addAction(acSaveData);

	m_Q_origin_bcm[0] = new QDoubleSpinBox(panel);
	m_Q_origin_bcm[1] = new QDoubleSpinBox(panel);
	m_Q_origin_bcm[2] = new QDoubleSpinBox(panel);
	m_Q_dir1_bcm[0] = new QDoubleSpinBox(panel);
	m_Q_dir1_bcm[1] = new QDoubleSpinBox(panel);
	m_Q_dir1_bcm[2] = new QDoubleSpinBox(panel);
	m_Q_dir2_bcm[0] = new QDoubleSpinBox(panel);
	m_Q_dir2_bcm[1] = new QDoubleSpinBox(panel);
	m_Q_dir2_bcm[2] = new QDoubleSpinBox(panel);
	m_num_Q_bcm[0] = new QSpinBox(panel);
	m_num_Q_bcm[1] = new QSpinBox(panel);
	m_num_Q_bcm[0]->setToolTip("Number of grid points along the first momentum axis.");
	m_num_Q_bcm[1]->setToolTip("Number of grid points along the second momentum axis.");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin_bcm[i]->setDecimals(4);
		m_Q_origin_bcm[i]->setMinimum(-99.9999);
		m_Q_origin_bcm[i]->setMaximum(+99.9999);
		m_Q_origin_bcm[i]->setSingleStep(0.01);
		m_Q_origin_bcm[i]->setValue(0.);
		m_Q_origin_bcm[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_origin_bcm[i]->setPrefix(hklPrefix[i]);
		m_Q_origin_bcm[i]->setToolTip("Starting momentum transfer.");

		m_Q_dir1_bcm[i]->setDecimals(4);
		m_Q_dir1_bcm[i]->setMinimum(-99.9999);
		m_Q_dir1_bcm[i]->setMaximum(+99.9999);
		m_Q_dir1_bcm[i]->setSingleStep(0.01);
		m_Q_dir1_bcm[i]->setValue(i == 0 ? 1. : 0.);
		m_Q_dir1_bcm[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir1_bcm[i]->setPrefix(hklPrefix[i]);
		m_Q_dir1_bcm[i]->setToolTip("Direction of momentum transfer along the first axis.");

		m_Q_dir2_bcm[i]->setDecimals(4);
		m_Q_dir2_bcm[i]->setMinimum(-99.9999);
		m_Q_dir2_bcm[i]->setMaximum(+99.9999);
		m_Q_dir2_bcm[i]->setSingleStep(0.01);
		m_Q_dir2_bcm[i]->setValue(i == 1 ? 1. : 0.);
		m_Q_dir2_bcm[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir2_bcm[i]->setPrefix(hklPrefix[i]);
		m_Q_dir2_bcm[i]->setToolTip("Direction of momentum transfer along the second axis.");
	}

	for(int i = 0; i < 2; ++i)
	{
		m_num_Q_bcm[i]->setMinimum(2);
		m_num_Q_bcm[i]->setMaximum(9999);
		m_num_Q_bcm[i]->setSingleStep(1);
		m_num_Q_bcm[i]->setValue(64);
	}

	// main dispersion button
	QPushButton *btnMainQ = new QPushButton("From Main Q", panel);
	btnMainQ->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	btnMainQ->setToolTip("Set the Q origin and directions from the dispersion in the main window.");

	// coordinate components
	m_coords_bcm[0] = new QSpinBox(panel);
	m_coords_bcm[0]->setMinimum(0);
	m_coords_bcm[0]->setMaximum(2);
	m_coords_bcm[0]->setValue(0);
	m_coords_bcm[0]->setPrefix("i = ");
	m_coords_bcm[0]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_coords_bcm[0]->setToolTip("First component index of B_ij matrix.");

	m_coords_bcm[1] = new QSpinBox(panel);
	m_coords_bcm[1]->setMinimum(0);
	m_coords_bcm[1]->setMaximum(2);
	m_coords_bcm[1]->setValue(1);
	m_coords_bcm[1]->setPrefix("j = ");
	m_coords_bcm[1]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_coords_bcm[1]->setToolTip("Second component index of B_ij matrix.");

	m_band_bcm = new QSpinBox(panel);
	m_band_bcm->setMinimum(0);
	m_band_bcm->setMaximum(999);
	m_band_bcm->setValue(0);
	m_band_bcm->setPrefix("band = ");
	m_band_bcm->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_band_bcm->setToolTip("Magnon band index.");

	// maximum cutoff for filtering numerical artefacts in berry curvature
	m_B_filter_enable_bcm = new QCheckBox("Maximum B:", panel);
	m_B_filter_enable_bcm->setChecked(true);
	m_B_filter_enable_bcm->setToolTip("Enable maximum cutoff Berry curvature for filtering numerical artefacts.");

	m_B_filter_bcm = new QDoubleSpinBox(panel);
	m_B_filter_bcm->setDecimals(2);
	m_B_filter_bcm->setMinimum(0.);
	m_B_filter_bcm->setMaximum(999999.99);
	m_B_filter_bcm->setSingleStep(1.);
	m_B_filter_bcm->setValue(m_B_filter_bc->maximum());
	m_B_filter_bcm->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_B_filter_bcm->setToolTip("Maximum cutoff Berry curvature for filtering numerical artefacts.");

	// minimum cutoff for filtering S(Q, E)
	m_S_filter_enable_bcm = new QCheckBox("Minimum S(Q, E):", panel);
	m_S_filter_enable_bcm->setChecked(false);
	m_S_filter_enable_bcm->setToolTip("Enable minimum S(Q, E).");

	m_S_filter_bcm = new QDoubleSpinBox(panel);
	m_S_filter_bcm->setDecimals(5);
	m_S_filter_bcm->setMinimum(0.);
	m_S_filter_bcm->setMaximum(9999.99999);
	m_S_filter_bcm->setSingleStep(0.01);
	m_S_filter_bcm->setValue(0.01);
	m_S_filter_bcm->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_S_filter_bcm->setToolTip("Minimum S(Q, E) to keep.");

	// progress bar
	m_progress_bcm = new QProgressBar(panel);
	m_progress_bcm->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btnStartStop_bcm = new QPushButton("Calculate", panel);
	m_btnStartStop_bcm->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	int y = 0;
	QGridLayout *grid = new QGridLayout(panel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);
	grid->addWidget(m_plot_bcm, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Q Origin (rlu):", this), y, 0, 1, 1);
	grid->addWidget(m_Q_origin_bcm[0], y, 1, 1, 1);
	grid->addWidget(m_Q_origin_bcm[1], y, 2, 1, 1);
	grid->addWidget(m_Q_origin_bcm[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Direction 1 (rlu):", this), y, 0, 1, 1);
	grid->addWidget(m_Q_dir1_bcm[0], y, 1, 1, 1);
	grid->addWidget(m_Q_dir1_bcm[1], y, 2, 1, 1);
	grid->addWidget(m_Q_dir1_bcm[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Direction 2 (rlu):", this), y, 0, 1, 1);
	grid->addWidget(m_Q_dir2_bcm[0], y, 1, 1, 1);
	grid->addWidget(m_Q_dir2_bcm[1], y, 2, 1, 1);
	grid->addWidget(m_Q_dir2_bcm[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Grid Points:", this), y, 0, 1, 1);
	grid->addWidget(m_num_Q_bcm[0], y, 1, 1, 1);
	grid->addWidget(m_num_Q_bcm[1], y, 2, 1, 1);
	grid->addWidget(btnMainQ, y++, 3, 1, 1);
	grid->addWidget(new QLabel("B Component:", panel), y, 0, 1, 1);
	grid->addWidget(m_coords_bcm[0], y, 1, 1, 1);
	grid->addWidget(m_coords_bcm[1], y, 2, 1, 1);
	grid->addWidget(m_band_bcm, y++, 3, 1, 1);
	grid->addWidget(m_B_filter_enable_bcm, y, 0, 1, 1);
	grid->addWidget(m_B_filter_bcm, y, 1, 1, 1);
	grid->addWidget(m_S_filter_enable_bcm, y, 2, 1, 1);
	grid->addWidget(m_S_filter_bcm, y++, 3, 1, 1);
	grid->addWidget(m_progress_bcm, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop_bcm, y++, 3, 1, 1);

	// calculation
	connect(m_btnStartStop_bcm, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled_bcm)
			CalculateBerryCurvatureMap();
		else
			m_stopRequested_bcm = true;
	});
	
	connect(btnMainQ, &QAbstractButton::clicked, this, &TopologyDlg::MapQFromMainQ);
	connect(m_B_filter_enable_bcm, &QCheckBox::toggled, m_B_filter_bcm, &QDoubleSpinBox::setEnabled);
	connect(m_S_filter_enable_bcm, &QCheckBox::toggled, m_S_filter_bcm, &QDoubleSpinBox::setEnabled);
	connect(m_plot_bcm, &QCustomPlot::mouseMove, this, &TopologyDlg::BerryCurvatureMapPlotMouseMove);
	connect(m_plot_bcm, &QCustomPlot::mousePress, this, &TopologyDlg::BerryCurvatureMapPlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &TopologyDlg::RescaleBerryCurvatureMapPlot);
	connect(acSaveFigure, &QAction::triggered, this, &TopologyDlg::SaveBerryCurvatureMapPlotFigure);
	connect(acSaveData, &QAction::triggered, this, &TopologyDlg::SaveBerryCurvatureMapData);

	m_B_filter_bcm->setEnabled(m_B_filter_enable_bcm->isChecked());
	m_S_filter_bcm->setEnabled(m_S_filter_enable_bcm->isChecked());
	EnableBerryCurvatureMapCalculation();

	return panel;
}



/**
 * set the Q position and directions from the main window's Q start and end points
 */
void TopologyDlg::MapQFromMainQ()
{
	if(m_Qstart.size() < 3 || m_Qend.size() < 3 || !m_dyn)
		return;

	const t_mat33_real& xtalB = m_dyn->GetCrystalBTrafo();
	const t_vec3_real* plane = m_dyn->GetScatteringPlane();
	if(!plane)
		return;

	// direction 1 is from the start to the end point
	t_vec3_real Qdir1 = m_Qend - m_Qstart;
	// direction 2 is perpendicular to direction 1 inside the scattering plane
	t_vec3_real Qdir2 = tl2::cross(xtalB, plane[2], Qdir1);

	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin_bcm[i]->setValue(m_Qstart[i]);
		m_Q_dir1_bcm[i]->setValue(Qdir1[i]);
		m_Q_dir2_bcm[i]->setValue(Qdir2[i]);
	}
}



/**
 * get the Q origin and direction vectors
 */
std::tuple<t_vec3_real, t_vec3_real, t_vec3_real> TopologyDlg::GetMapQVectors() const
{
	t_vec3_real Q_origin = tl2::create<t_vec3_real>(
	{
		(t_real)m_Q_origin_bcm[0]->value(),
		(t_real)m_Q_origin_bcm[1]->value(),
		(t_real)m_Q_origin_bcm[2]->value(),
	});

	t_vec3_real Q_dir_1 = tl2::create<t_vec3_real>(
	{
		(t_real)m_Q_dir1_bcm[0]->value(),
		(t_real)m_Q_dir1_bcm[1]->value(),
		(t_real)m_Q_dir1_bcm[2]->value(),
	});

	t_vec3_real Q_dir_2 = tl2::create<t_vec3_real>(
	{
		(t_real)m_Q_dir2_bcm[0]->value(),
		(t_real)m_Q_dir2_bcm[1]->value(),
		(t_real)m_Q_dir2_bcm[2]->value(),
	});

	return std::make_tuple(std::move(Q_origin), std::move(Q_dir_1), std::move(Q_dir_2));
}



/**
 * converts array indices to Q position
 */
t_vec3_real TopologyDlg::GetMapQFromIndices(std::size_t idx1, std::size_t idx2) const
{
	t_size Q_count_1 = m_num_Q_bcm[0]->value();
	t_size Q_count_2 = m_num_Q_bcm[1]->value();

	// get coordinates
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetMapQVectors();

	t_vec3_real Q_step_1 = Q_dir_1 / t_real(Q_count_1 - 1);
	t_vec3_real Q_step_2 = Q_dir_2 / t_real(Q_count_2 - 1);

	return Q_origin + Q_step_1*t_real(idx1) + Q_step_2*t_real(idx2);
}



/**
 * toggle between "calculate" and "stop" button
 */
void TopologyDlg::EnableBerryCurvatureMapCalculation(bool enable)
{
	m_calcEnabled_bcm = enable;

	if(enable)
	{
		m_btnStartStop_bcm->setText("Calculate");
		m_btnStartStop_bcm->setToolTip("Start calculation.");
		if(g_use_icons)
			m_btnStartStop_bcm->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop_bcm->setText("Stop");
		m_btnStartStop_bcm->setToolTip("Stop running calculation.");
		if(g_use_icons)
			m_btnStartStop_bcm->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * calculate the berry curvature
 */
void TopologyDlg::CalculateBerryCurvatureMap()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableBerryCurvatureMapCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableBerryCurvatureMapCalculation(false);

	ClearBerryCurvatureMapPlot(false);

	// get coordinates
	//auto [Q_origin, Q_dir1, Q_dir2] = GetMapQVectors();

/*
	// get Q component with maximum range
	t_vec3_real Q_range = Q_end - Q_start;
	m_Q_idx_bc = 0;
	if(std::abs(Q_range[1]) > std::abs(Q_range[m_Q_idx_bc]))
		m_Q_idx_bc = 1;
	if(std::abs(Q_range[2]) > std::abs(Q_range[m_Q_idx_bc]))
		m_Q_idx_bc = 2;

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx_bc] > Q_end[m_Q_idx_bc])
		std::swap(Q_start, Q_end);

	// Q range
	m_Q_min_bc = Q_start[m_Q_idx_bc];
	m_Q_max_bc = Q_end[m_Q_idx_bc];
*/
	// get settings
	t_size Q_count_1 = m_num_Q_bcm[0]->value();
	t_size Q_count_2 = m_num_Q_bcm[1]->value();

	std::vector<t_size> *perm = nullptr;
	t_size dim1 = m_coords_bcm[0]->value();
	t_size dim2 = m_coords_bcm[1]->value();

	// calculate berry curvature
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested_bcm = false;
	m_progress_bcm->setMinimum(0);
	m_progress_bcm->setMaximum(Q_count_1*Q_count_2);
	m_progress_bcm->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count_1*Q_count_2);

	//m_data_bcm.clear();
	//m_data_bcm.reserve(Q_count);

	const bool enforce_commutator = g_enforce_commutator;

	for(t_size Q_idx_1 = 0; Q_idx_1 < Q_count_1; ++Q_idx_1)
	for(t_size Q_idx_2 = 0; Q_idx_2 < Q_count_2; ++Q_idx_2)
	{
		auto task = [this, &mtx, &dyn, Q_idx_1, Q_idx_2, Q_count_1, Q_count_2, perm,
			dim1, dim2, enforce_commutator]()
		{
			// calculate the berry curvature at the given Q point
			const t_vec3_real Q = GetMapQFromIndices(Q_idx_1, Q_idx_2);
/*
			BerryCurvatureData data_bc;
			data_bc.momentum = Q;
			typename t_magdyn::SofQE S;
			std::tie(data_bc.curvatures, S) = dyn.CalcBerryCurvatures(
				Q, g_delta_diff, perm, dim1, dim2, g_evecs_ortho != 0,
				enforce_commutator);
			t_size num_bands = data_bc.curvatures.size();
			data_bc.energies.reserve(num_bands);
			data_bc.weights.reserve(num_bands);

			// calculate energies per band
			assert(S.E_and_S.size() == num_bands);
			for(t_size band = 0; band < num_bands; ++band)
			{
				data_bc.energies.push_back(S.E_and_S[band].E);
				data_bc.weights.push_back(S.E_and_S[band].weight_perp);
			}

			std::lock_guard<std::mutex> _lck{mtx};
			m_data_bc.emplace_back(std::move(data_bc));
*/
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		// process events to see if the stop button was clicked
		// only do this for a fraction of the points to avoid gui overhead
		bool process_evts = (task_idx % std::max<t_size>(tasks.size() / g_stop_check_fraction, 1) == 0);
		if(process_evts)
			qApp->processEvents();

		if(m_stopRequested_bcm)
		{
			pool.stop();
			break;
		}

		task->get_future().get();

		if(process_evts || task_idx + 1 == tasks.size())
			m_progress_bcm->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested_bc)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort raw unfiltered data by Q
	/*std::vector<std::size_t> perm_all = tl2::get_perm(m_data_bc.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data_bc[idx1].momentum[m_Q_idx_bc]
			< m_data_bc[idx2].momentum[m_Q_idx_bc];
	});

	m_data_bcm = tl2::reorder(m_data_bc, perm_all);*/

	PlotBerryCurvatureMap();
}



/**
 * clears the berry curature map
 */
void TopologyDlg::ClearBerryCurvatureMapPlot(bool replot)
{
	if(!m_plot_bcm || !m_plot_map_bcm)
		return;

	m_plot_map_bcm->data()->setSize(0, 0);
	//m_plot_bcm->clearPlottables();

	if(replot)
		m_plot_bcm->replot();
}



/**
 * rescale plot axes to fit the content
 */
void TopologyDlg::RescaleBerryCurvatureMapPlot()
{
	if(!m_plot_bcm || !m_plot_map_bcm)
		return;

	m_plot_map_bcm->rescaleDataRange();
	m_plot_bcm->rescaleAxes();
	m_plot_bcm->replot();
}



void TopologyDlg::PlotBerryCurvatureMap()
{
	
}



/**
 * save plot as image file
 */
void TopologyDlg::SaveBerryCurvatureMapPlotFigure()
{
	if(!m_plot_bcm)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("topology/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("topology/dir", QFileInfo(filename).path());

	if(!m_plot_bcm->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void TopologyDlg::SaveBerryCurvatureMapData()
{
	
}



/**
 * show current cursor coordinates
 */
void TopologyDlg::BerryCurvatureMapPlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q1 = m_plot_bcm->xAxis->pixelToCoord(evt->pos().x());
	t_real Q2 = m_plot_bcm->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 rlu, Q2 = %2 rlu.");
	status = status.arg(Q1, 0, 'g', g_prec_gui).arg(Q2, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void TopologyDlg::BerryCurvatureMapPlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot_bcm)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot_bcm->popup(pos);
		evt->accept();
	}
}
// ============================================================================
