/**
 * magnetic dynamics -- powder excitations
 * @author Tobias Weber <tweber@ill.fr>
 * @date July 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
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

#include "powder.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



// ============================================================================
// powder dialog
// ============================================================================

/**
 * sets up the differentiation dialog
 */
PowderDlg::PowderDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Powder Spectrum");
	setSizeGripEnabled(true);

	// status bar
	m_status = new QLabel(this);
	m_status->setFrameShape(QFrame::Panel);
	m_status->setFrameShadow(QFrame::Sunken);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// close button
	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// save data button
	QPushButton *btnSaveData = btnbox->addButton("Save Data...", QDialogButtonBox::ActionRole);
	btnSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	btnSaveData->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// main grid
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(CreatePowderPanel(), 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnbox, 1, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("diff/geo"))
			restoreGeometry(m_sett->value("powder/geo").toByteArray());
		else
			resize(640, 800);
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &PowderDlg::accept);
	connect(btnSaveData, &QAbstractButton::clicked, this, &PowderDlg::SavePowderData);
}



PowderDlg::~PowderDlg()
{
}



/**
 * set a pointer to the main magdyn kernel
 */
void PowderDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * set the Q start and end points from the main window's dispersion
 */
void PowderDlg::SetDispersionQE(
	const t_vec_real& Qstart, const t_vec_real& Qend,
	const t_real Estart, const t_real Eend)
{
	m_Qstart = Qstart;
	m_Qend = Qend;

	m_Estart = Estart;
	m_Eend = Eend;
}



void PowderDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * dialog is closing
 */
void PowderDlg::accept()
{
	if(m_sett)
		m_sett->setValue("powder/geo", saveGeometry());

	QDialog::accept();
}
// ============================================================================



// ============================================================================
// calculate powder spectrum
// ============================================================================

/**
 * create the panel for the powder tab
 */
QWidget* PowderDlg::CreatePowderPanel()
{
	QWidget *panelPowder = new QWidget(this);

	// plotter
	m_plot_powder = new QCustomPlot(panelPowder);
	m_plot_powder->setFont(font());
	m_plot_powder->xAxis->setScaleType(QCPAxis::stLinear);
	m_plot_powder->yAxis->setScaleType(QCPAxis::stLinear);
	m_plot_powder->xAxis->setLabel("Momentum Transfer Q (Å⁻¹)");
	m_plot_powder->yAxis->setLabel("Energy Transfer E (meV)");
	m_plot_powder->setInteraction(QCP::iRangeDrag, true);
	m_plot_powder->setInteraction(QCP::iRangeZoom, true);
	m_plot_powder->setSelectionRectMode(QCP::srmZoom);
	m_plot_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	m_plot_colour = new QCPColorScale(m_plot_powder);
	m_plot_colour->axis()->setScaleType(QCPAxis::stLinear);
	m_plot_colour->axis()->setLabel("S(Q, E)");
	m_plot_colour->setRangeDrag(true);
	m_plot_colour->setRangeZoom(true);
	m_plot_colour->setType(QCPAxis::atRight);
	m_plot_powder->plotLayout()->addElement(0, 1, m_plot_colour);

	m_plot_map = new QCPColorMap(m_plot_powder->xAxis, m_plot_powder->yAxis);
	m_plot_map->setTightBoundary(true);
	m_plot_map->setColorScale(m_plot_colour);
	m_plot_map->setGradient(QCPColorGradient::gpHot);
	m_plot_map->setInterpolate(true);
	m_plot_map->setDataScaleType(QCPAxis::stLinear);

	// context menu for plotter
	m_menuPlot_powder = new QMenu("Plotter", panelPowder);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot_powder);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot_powder);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot_powder);

	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));

	m_menuPlot_powder->addAction(acRescalePlot);
	m_menuPlot_powder->addSeparator();
	m_menuPlot_powder->addAction(acSaveFigure);
	m_menuPlot_powder->addAction(acSaveData);

	// start and stop Q
	m_Q_start_powder = new QDoubleSpinBox(panelPowder);
	m_Q_start_powder->setDecimals(4);
	m_Q_start_powder->setMinimum(0.);
	m_Q_start_powder->setMaximum(99.9999);
	m_Q_start_powder->setSingleStep(0.01);
	m_Q_start_powder->setValue(0.);
	m_Q_start_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_Q_start_powder->setToolTip("Dispersion initial momentum transfer (Å⁻¹).");

	m_Q_end_powder = new QDoubleSpinBox(panelPowder);
	m_Q_end_powder->setDecimals(4);
	m_Q_end_powder->setMinimum(0.);
	m_Q_end_powder->setMaximum(99.9999);
	m_Q_end_powder->setSingleStep(0.01);
	m_Q_end_powder->setValue(1.);
	m_Q_end_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_Q_end_powder->setToolTip("Dispersion final momentum transfer (Å⁻¹).");

	// number of Q points in the plot
	m_num_Q_powder = new QSpinBox(panelPowder);
	m_num_Q_powder->setMinimum(1);
	m_num_Q_powder->setMaximum(99999);
	m_num_Q_powder->setValue(64);
	m_num_Q_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_Q_powder->setToolTip("Number of Q points to calculate.");

	// number of Q vectors per point
	m_num_Qvecs_powder = new QSpinBox(panelPowder);
	m_num_Qvecs_powder->setMinimum(1);
	m_num_Qvecs_powder->setMaximum(99999);
	m_num_Qvecs_powder->setValue(512);
	m_num_Qvecs_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_Qvecs_powder->setToolTip("Number of Q vectors to calculate per point.");

	// start and stop E
	m_E_start_powder = new QDoubleSpinBox(panelPowder);
	m_E_start_powder->setDecimals(4);
	m_E_start_powder->setMinimum(0.);
	m_E_start_powder->setMaximum(999.9999);
	m_E_start_powder->setSingleStep(0.01);
	m_E_start_powder->setValue(0.);
	m_E_start_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_E_start_powder->setToolTip("Dispersion initial energy transfer (meV).");

	m_E_end_powder = new QDoubleSpinBox(panelPowder);
	m_E_end_powder->setDecimals(4);
	m_E_end_powder->setMinimum(0.);
	m_E_end_powder->setMaximum(999.9999);
	m_E_end_powder->setSingleStep(0.01);
	m_E_end_powder->setValue(1.);
	m_E_end_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_E_end_powder->setToolTip("Dispersion final energy transfer (meV).");

	// number of E points in the plot
	m_num_E_powder = new QSpinBox(panelPowder);
	m_num_E_powder->setMinimum(1);
	m_num_E_powder->setMaximum(99999);
	m_num_E_powder->setValue(64);
	m_num_E_powder->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_E_powder->setToolTip("Number of E points to calculate.");

	// orthogonal projector
	m_use_proj = new QCheckBox(panelPowder);
	m_use_proj->setText("Use Projector");
	m_use_proj->setToolTip("Use the projector orthogonal to Q in neutron scattering.");
	m_use_proj->setChecked(true);

	QCheckBox *use_log = new QCheckBox(panelPowder);
	use_log->setText("Logarithmic");
	use_log->setToolTip("Use logarithmic S(Q, E) scale.");
	use_log->setChecked(false);

	// dispersion Q button
	QPushButton *btnQE = new QPushButton("Set Main Q && E", panelPowder);
	btnQE->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	btnQE->setToolTip("Set the Q and E start and end points from the dispersion in the main window.");

	// progress bar
	m_progress_powder = new QProgressBar(panelPowder);
	m_progress_powder->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btnStartStop_powder = new QPushButton("Calculate", panelPowder);
	m_btnStartStop_powder->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	// component grid
	auto grid = new QGridLayout(panelPowder);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_plot_powder, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Start Q (Å⁻¹):", panelPowder), y, 0, 1, 1);
	grid->addWidget(m_Q_start_powder, y, 1, 1, 1);
	grid->addWidget(new QLabel("End Q (Å⁻¹):", panelPowder), y, 2, 1, 1);
	grid->addWidget(m_Q_end_powder, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Count:", panelPowder), y, 0, 1, 1);
	grid->addWidget(m_num_Q_powder, y, 1, 1, 1);
	grid->addWidget(new QLabel("Q Vec. Count:", panelPowder), y, 2, 1, 1);
	grid->addWidget(m_num_Qvecs_powder, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Start E (meV):", panelPowder), y, 0, 1, 1);
	grid->addWidget(m_E_start_powder, y, 1, 1, 1);
	grid->addWidget(new QLabel("End E (meV):", panelPowder), y, 2, 1, 1);
	grid->addWidget(m_E_end_powder, y++, 3, 1, 1);
	grid->addWidget(new QLabel("E Count:", panelPowder), y, 0, 1, 1);
	grid->addWidget(m_num_E_powder, y, 1, 1, 1);
	grid->addWidget(m_use_proj, y++, 2, 1, 2);
	grid->addWidget(use_log, y, 0, 1, 2);
	grid->addWidget(btnQE, y++, 3, 1, 1);
	grid->addWidget(m_progress_powder, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop_powder, y++, 3, 1, 1);

	// connections
	connect(m_plot_powder, &QCustomPlot::mouseMove, this, &PowderDlg::PowderPlotMouseMove);
	connect(m_plot_powder, &QCustomPlot::mousePress, this, &PowderDlg::PowderPlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &PowderDlg::RescalePowderPlot);
	connect(acSaveFigure, &QAction::triggered, this, &PowderDlg::SavePowderPlotFigure);
	connect(acSaveData, &QAction::triggered, this, &PowderDlg::SavePowderData);
	connect(btnQE, &QAbstractButton::clicked, this, &PowderDlg::SetPowderQE);

	connect(use_log, &QCheckBox::toggled, [this](bool checked)
	{
		if(!m_plot_colour || !m_plot_powder)
			return;

		m_plot_colour->axis()->setScaleType(
			checked ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
		m_plot_powder->replot();
	});

	// calculation
	connect(m_btnStartStop_powder, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled_powder)
			CalculatePowder();
		else
			m_stopRequested_powder = true;
	});

	EnablePowderCalculation();
	return panelPowder;
}



/**
 * plots the powder spectra
 */
void PowderDlg::PlotPowder()
{
	if(!m_plot_powder)
		return;

	ClearPowderPlot(false);

	const t_size Qvecs_count = m_num_Qvecs_powder->value();
	const t_size Q_count = m_data_powder.size(); //m_num_Q_powder->value();
	const t_size E_count = m_num_E_powder->value();

	m_plot_map->data()->setSize((int)Q_count, (int)E_count);
	m_plot_map->data()->setRange(
		QCPRange{ m_Q_min_powder, m_Q_max_powder },
		QCPRange{ m_E_min_powder, m_E_max_powder });

	for(t_size Qidx = 0; Qidx < Q_count; ++Qidx)
	{
		const auto& E_histo = m_data_powder[Qidx].histogram;

		t_size Eidx = 0;
		//for(t_size Eidx = 0; Eidx < E_count; ++Eidx)
		for(const auto& Ebinidx : boost::histogram::indexed(E_histo))
		{
			t_real S = *Ebinidx / t_real(Qvecs_count) / t_real(E_count);
			m_plot_map->data()->setCell((int)Qidx, (int)Eidx, S);
			++Eidx;
		}
	}

	RescalePowderPlot();
}



/**
 * calculate the powder spectra
 */
void PowderDlg::CalculatePowder()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnablePowderCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnablePowderCalculation(false);

	ClearPowderPlot(false);

	// get coordinates
	t_real Q_start = m_Q_start_powder->value();
	t_real Q_end = m_Q_end_powder->value();
	t_real E_start = m_E_start_powder->value();
	t_real E_end = m_E_end_powder->value();

	// keep Q and E in ascending order
	if(Q_start > Q_end)
		std::swap(Q_start, Q_end);
	if(E_start > E_end)
		std::swap(E_start, E_end);

	// final ranges
	m_Q_min_powder = Q_start;
	m_Q_max_powder = Q_end;
	m_E_min_powder = E_start;
	m_E_max_powder = E_end;

	// get settings
	t_size Q_count = m_num_Q_powder->value();
	t_size Qvec_count = m_num_Qvecs_powder->value();
	t_size E_count = m_num_E_powder->value();
	bool ortho_proj = m_use_proj->isChecked();

	// calculate powder spectra
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested_powder = false;
	m_progress_powder->setMinimum(0);
	m_progress_powder->setMaximum(Q_count);
	m_progress_powder->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count);

	m_data_powder.clear();
	m_data_powder.reserve(Q_count);

	for(t_size Q_idx = 0; Q_idx < Q_count; ++Q_idx)
	{
		auto task = [this, &mtx, &dyn,
			Q_start, Q_end, Q_idx, Q_count, Qvec_count,
			E_start, E_end, E_count, ortho_proj]()
		{
			const t_real Q = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(Q_idx) / t_real(Q_count - 1))
				: Q_start;

			PowderData data_powder;
			data_powder.momentum = Q;
			data_powder.histogram = dyn.CalcPowderBin(Q, E_start, E_end, E_count,
				Qvec_count, 1, true, ortho_proj);

			std::lock_guard<std::mutex> _lck{mtx};
			m_data_powder.emplace_back(std::move(data_powder));
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

		if(m_stopRequested_powder)
		{
			pool.stop();
			break;
		}

		task->get_future().get();

		if(process_evts || task_idx + 1 == tasks.size())
			m_progress_powder->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested_powder)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort raw unfiltered data by Q
	std::vector<std::size_t> perm_all = tl2::get_perm(m_data_powder.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data_powder[idx1].momentum
			< m_data_powder[idx2].momentum;
	});

	m_data_powder = tl2::reorder(m_data_powder, perm_all);

	PlotPowder();
}



/**
 * clears the dispersion graph
 */
void PowderDlg::ClearPowderPlot(bool replot)
{
	if(!m_plot_powder || !m_plot_map)
		return;

	m_plot_map->data()->setSize(0, 0);
	//m_plot_powder->clearPlottables();

	if(replot)
		m_plot_powder->replot();
}



/**
 * show current cursor coordinates
 */
void PowderDlg::PowderPlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot_powder->xAxis->pixelToCoord(evt->pos().x());
	t_real E = m_plot_powder->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 Å⁻¹, E = %2 meV.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(E, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void PowderDlg::PowderPlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot_powder)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot_powder->popup(pos);
		evt->accept();
	}
}



/**
 * rescale plot axes to fit the content
 */
void PowderDlg::RescalePowderPlot()
{
	if(!m_plot_powder || !m_plot_map)
		return;

	m_plot_map->rescaleDataRange();
	m_plot_powder->rescaleAxes();
	m_plot_powder->replot();
}



/**
 * save plot as image file
 */
void PowderDlg::SavePowderPlotFigure()
{
	if(!m_plot_powder)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("powder/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("powder/dir", QFileInfo(filename).path());

	if(!m_plot_powder->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void PowderDlg::SavePowderData()
{
	if(m_data_powder.size() == 0)
	{
		ShowError("No data to save.");
		return;
	}

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("powder/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("powder/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename).toStdString().c_str());
		return;
	}

	ofstr.precision(g_prec);
	int field_len = g_prec * 2.5;

	// write meta header
	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr << "#\n"
		<< "# Created by Magpie " << MAGPIE_VER << "\n"
		<< "# Author: Tobias Weber\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n\n";

	const t_size Qvecs_count = m_num_Qvecs_powder->value();
	const t_size Q_count = m_data_powder.size();
	const t_size E_count = m_num_E_powder->value();
	ofstr << "# row_label: Q (1/A)" << "\n";
	ofstr << "# col_label: E (meV)" << "\n";
	ofstr << "# row_count: " << Q_count << "\n";
	ofstr << "# col_count: " << E_count << "\n";

	std::vector<t_real> Qs;
	Qs.reserve(Q_count);
	ofstr << "# row_values: ";
	for(const PowderData& data : m_data_powder)
	{
		ofstr << data.momentum << " ";
		Qs.push_back(data.momentum);
	}
	ofstr << "\n";

	std::vector<t_real> Es;
	Es.reserve(E_count);
	ofstr << "# col_values: ";
	for(const auto& Ebinidx : boost::histogram::indexed(m_data_powder[0].histogram))
	{
		const auto& Ebin = Ebinidx.bin();
		const t_real E = Ebin.lower() + (Ebin.upper() - Ebin.lower()) / 2.;
		ofstr << E << " ";
		Es.push_back(E);
	}
	ofstr << "\n";

	// plot command
	t_real E_range = *Es.rbegin() - *Es.begin();
	t_real Q_range = *Qs.rbegin() - *Qs.begin();
	ofstr << "# plot_with: plot \"" << QFileInfo(filename).fileName().toStdString() << "\" "
		<< "matrix using "
		<< "($2*" << Q_range << "/" << Q_count << "+" << Qs[0] << "):"
		<< "($1*" << E_range << "/" << E_count << "+" << Es[0] << "):"
		<< "($3) with image";
	ofstr << "\n\n";

	// write data
	for(const PowderData& data : m_data_powder)
	{
		const auto& E_histo = data.histogram;

		for(const auto& Ebinidx : boost::histogram::indexed(E_histo))
		{
			t_real S = *Ebinidx / t_real(Qvecs_count) / t_real(E_count);
			ofstr << std::setw(field_len) << std::left << S << " ";
		}
	
		ofstr << "\n";
	}

	ofstr << std::endl;
}



/**
 * toggle between "calculate" and "stop" button
 */
void PowderDlg::EnablePowderCalculation(bool enable)
{
	m_calcEnabled_powder = enable;

	if(enable)
	{
		m_btnStartStop_powder->setText("Calculate");
		m_btnStartStop_powder->setToolTip("Start calculation.");
		m_btnStartStop_powder->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop_powder->setText("Stop");
		m_btnStartStop_powder->setToolTip("Stop running calculation.");
		m_btnStartStop_powder->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * sets the Q and E ranges to the main window's dispersion
 */
void PowderDlg::SetPowderQE()
{
	// E
	m_E_start_powder->setValue(m_Estart);
	m_E_end_powder->setValue(m_Eend);

	// Q
	if(m_Qstart.size() < 3 || m_Qend.size() < 3 || !m_dyn)
		return;

	const t_mat_real& B = m_dyn->GetCrystalBTrafo();

	t_vec_real Qvec_start_invA = B * m_Qstart;
	t_real Q_start_invA = tl2::norm(Qvec_start_invA);
	m_Q_start_powder->setValue(Q_start_invA);

	t_vec_real Qvec_end_invA = B * m_Qend;
	t_real Q_end_invA = tl2::norm(Qvec_end_invA);
	m_Q_end_powder->setValue(Q_end_invA);
}
// ============================================================================
