/**
 * magnetic dynamics -- form factor plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date May 2026
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

#include "ffact.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



// ============================================================================
// form factors dialog
// ============================================================================

/**
 * sets up the form factor dialog
 */
FormFactorDlg::FormFactorDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Magnetic Form Factors");
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

	// main grid
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(CreateFormFactorPanel(), 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnbox, 1, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("ffact/geo"))
			restoreGeometry(m_sett->value("ffact/geo").toByteArray());
		else
			resize(640, 800);

		if(m_sett->contains("ffact/splitter"))
			m_split_plot_ff->restoreState(m_sett->value("ffact/splitter").toByteArray());
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &FormFactorDlg::accept);
}



FormFactorDlg::~FormFactorDlg()
{
}



void FormFactorDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * dialog is closing
 */
void FormFactorDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("ffact/geo", saveGeometry());
		m_sett->setValue("ffact/splitter", m_split_plot_ff->saveState());
	}

	QDialog::accept();
}
// ============================================================================



// ============================================================================
// calculate form factors
// ============================================================================

/**
 * parse the form factor formulas and plot them
 */
void FormFactorDlg::PlotFormFactors(const std::vector<std::string>& ffacts)
{
	m_parsers_ff.clear();

	for(t_size idx = 0; idx < ffacts.size(); ++idx)
	{
		const std::string& ffact = ffacts[idx];
		tl2::ExprParser<t_cplx> parser;

		parser.SetAutoregisterVariables(false);
		parser.SetInvalid0(false);

		parser.register_var("Q", 0.);
		parser.register_var("Q2", 0.);
		parser.register_var("s", 0.);
		parser.register_var("s2", 0.);

		if(!parser.parse_noexcept(ffact))
		{
			std::cerr << "Magdyn error: Magnetic form factor formula #" << idx
				<< ": \"" << ffact << "\" could not be parsed."
				<< std::endl;
		}

		m_parsers_ff.emplace_back(std::move(parser));
	}

	CalculateFormFactors();
}



/**
 * column indices in table for the form factor indices
 */
enum : int
{
	COL_FF_INDEX = 0,
	COL_FF_ACTIVE,
	NUM_COLS_FF,
};



/**
 * create the panel for the form factor tab
 */
QWidget* FormFactorDlg::CreateFormFactorPanel()
{
	QWidget *panelFFact = new QWidget(this);

	// plotter
	m_plot_ff = new QCustomPlot(panelFFact);
	m_plot_ff->setFont(font());
	m_plot_ff->xAxis->setLabel("Momentum Transfer Q (Å⁻¹)");
	m_plot_ff->yAxis->setLabel("Magnetic Form Factor F_M(Q)");
	m_plot_ff->setInteraction(QCP::iRangeDrag, true);
	m_plot_ff->setInteraction(QCP::iRangeZoom, true);
	m_plot_ff->setSelectionRectMode(QCP::srmZoom);
	m_plot_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// form factor curves
	QWidget *panel_ffsel = new QWidget(panelFFact);
	m_table_ff = new QTableWidget(panel_ffsel);
	m_table_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_ff->setShowGrid(true);
	m_table_ff->setSortingEnabled(false);
	m_table_ff->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_ff->setSelectionMode(QTableWidget::SingleSelection);
	m_table_ff->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_ff->verticalHeader()->setVisible(false);
	m_table_ff->setColumnCount(NUM_COLS_FF);
	m_table_ff->setHorizontalHeaderItem(COL_FF_INDEX, new QTableWidgetItem{"Idx."});
	m_table_ff->setHorizontalHeaderItem(COL_FF_ACTIVE, new QTableWidgetItem{"Act."});
	m_table_ff->setColumnWidth(COL_FF_INDEX, 40);
	m_table_ff->setColumnWidth(COL_FF_ACTIVE, 25);
	m_table_ff->resizeColumnsToContents();

	// splitter for plot and indices list
	m_split_plot_ff = new QSplitter(panelFFact);
	m_split_plot_ff->setOrientation(Qt::Horizontal);
	m_split_plot_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot_ff->addWidget(m_plot_ff);
	m_split_plot_ff->addWidget(panel_ffsel);
	m_split_plot_ff->setCollapsible(0, false);
	m_split_plot_ff->setCollapsible(1, true);
	m_split_plot_ff->setStretchFactor(m_split_plot_ff->indexOf(m_plot_ff), 24);
	m_split_plot_ff->setStretchFactor(m_split_plot_ff->indexOf(panel_ffsel), 1);

	// context menu for plotter
	m_menuPlot_ff = new QMenu("Plotter", panelFFact);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot_ff);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot_ff);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot_ff);

	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));

	m_menuPlot_ff->addAction(acRescalePlot);
	m_menuPlot_ff->addSeparator();
	m_menuPlot_ff->addAction(acSaveFigure);
	m_menuPlot_ff->addAction(acSaveData);

	// indices panel grid
	int y_indices = 0;
	QGridLayout *grid_indices = new QGridLayout(panel_ffsel);
	grid_indices->setSpacing(4);
	grid_indices->setContentsMargins(6, 6, 6, 6);
	grid_indices->addWidget(m_table_ff, y_indices++, 0, 1, 1);

	// start and stop coordinates
	m_Q_start_ff = new QDoubleSpinBox(panelFFact);
	m_Q_end_ff = new QDoubleSpinBox(panelFFact);

	m_Q_start_ff->setToolTip("Dispersion initial momentum transfer, Q (Å⁻¹).");
	m_Q_end_ff->setToolTip("Dispersion final momentum transfer, Q (Å⁻¹).");

	m_Q_start_ff->setDecimals(4);
	m_Q_start_ff->setMinimum(0.);
	m_Q_start_ff->setMaximum(+99.9999);
	m_Q_start_ff->setSingleStep(0.01);
	m_Q_start_ff->setValue(0.);
	m_Q_start_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_Q_start_ff->setPrefix("Q = ");

	m_Q_end_ff->setDecimals(4);
	m_Q_end_ff->setMinimum(0.);
	m_Q_end_ff->setMaximum(+99.9999);
	m_Q_end_ff->setSingleStep(0.01);
	m_Q_end_ff->setValue(5.);
	m_Q_end_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_Q_end_ff->setPrefix("Q = ");

	// number of Q points in the plot
	m_num_Q_ff = new QSpinBox(panelFFact);
	m_num_Q_ff->setMinimum(1);
	m_num_Q_ff->setMaximum(99999);
	m_num_Q_ff->setValue(128);
	m_num_Q_ff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_Q_ff->setToolTip("Number of Q points to calculate.");

	// plot |F|^2
	m_ff_squared = new QCheckBox(panelFFact);
	m_ff_squared->setText("Plot |F(Q)|^2");
	m_ff_squared->setChecked(false);

	// progress bar
	m_progress_ff = new QProgressBar(panelFFact);
	m_progress_ff->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btnStartStop_ff = new QPushButton("Calculate", panelFFact);

	// component grid
	auto grid = new QGridLayout(panelFFact);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_split_plot_ff, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Start Q (Å⁻¹):", panelFFact), y, 0, 1, 1);
	grid->addWidget(m_Q_start_ff, y, 1, 1, 1);
	grid->addWidget(new QLabel("End Q (Å⁻¹):", panelFFact), y, 2, 1, 1);
	grid->addWidget(m_Q_end_ff, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Count:", panelFFact), y, 0, 1, 1);
	grid->addWidget(m_num_Q_ff, y, 1, 1, 1);
	grid->addWidget(m_ff_squared, y++, 2, 1, 2);
	grid->addWidget(m_progress_ff, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop_ff, y++, 3, 1, 1);

	// connections
	connect(m_ff_squared, &QCheckBox::toggled, [this]() { this->PlotFormFactors(false); });
	connect(m_plot_ff, &QCustomPlot::mouseMove, this, &FormFactorDlg::FormFactorPlotMouseMove);
	connect(m_plot_ff, &QCustomPlot::mousePress, this, &FormFactorDlg::FormFactorPlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &FormFactorDlg::RescaleFormFactorPlot);
	connect(acSaveFigure, &QAction::triggered, this, &FormFactorDlg::SaveFormFactorPlotFigure);
	connect(acSaveData, &QAction::triggered, this, &FormFactorDlg::SaveFormFactorData);

	// calculation
	connect(m_btnStartStop_ff, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled_ff)
			CalculateFormFactors();
		else
			m_stopRequested_ff = true;
	});

	EnableFormFactorCalculation();
	return panelFFact;
}



/**
 * clears the table of form factor indices
 */
void FormFactorDlg::ClearFormFactorIndices()
{
	m_table_ff->clearContents();
	m_table_ff->setRowCount(0);
}



/**
 * adds a form factor index to the table
 */
void FormFactorDlg::AddFormFactorIndex(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_ff)
		return;

	int row = m_table_ff->rowCount();
	m_table_ff->insertRow(row);

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

	QCheckBox *checkIndex = new QCheckBox(m_table_ff);
	checkIndex->setChecked(enabled);
	connect(checkIndex, &QCheckBox::toggled, [this]() { PlotFormFactors(false); });

	m_table_ff->setItem(row, COL_FF_INDEX, item);
	m_table_ff->setCellWidget(row, COL_FF_ACTIVE, checkIndex);
}



/**
 * verifies if the index' checkbox is checked
 */
bool FormFactorDlg::IsFormFactorIndexEnabled(t_size idx) const
{
	if(!m_table_ff || int(idx) >= m_table_ff->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(m_table_ff->cellWidget(int(idx), COL_FF_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * calculate the data sets and plot the form factors
 */
void FormFactorDlg::PlotFormFactors(bool clear_settings)
{
	if(!m_plot_ff)
		return;

	// keep some settings from previous plot, e.g. the index visibility flags
	std::vector<bool> enabled_indices;
	if(!clear_settings)
	{
		enabled_indices.reserve(m_table_ff->rowCount());
		for(int row = 0; row < m_table_ff->rowCount(); ++row)
			enabled_indices.push_back(IsFormFactorIndexEnabled(t_size(row)));
	}

	ClearFormFactorIndices();
	ClearFormFactorPlot(false);

	if(m_data_ff.size() == 0)
	{
		m_plot_ff->replot();
		return;
	}

	t_size num_Q = m_data_ff.size();
	t_size num_indices = m_data_ff[0].ffacts.size();
	const bool ff_squared = m_ff_squared->isChecked();

	m_plot_ff->yAxis->setLabel(ff_squared
		? "Magnetic Form Factor |F_M(Q)|^2"
		: "Magnetic Form Factor Re(F_M(Q))");

	std::vector<QVector<qreal>> Qs_data_ff{num_indices};
	std::vector<QVector<qreal>> ffact_data_ff{num_indices};

	for(t_size Q_idx = 0; Q_idx < num_Q; ++Q_idx)
	{
		t_real Q = m_data_ff[Q_idx].momentum;

		for(t_size idx = 0; idx < num_indices; ++idx)
		{
			Qs_data_ff[idx].push_back(Q);
			ffact_data_ff[idx].push_back(ff_squared
				? m_data_ff[Q_idx].ffacts2[idx]
				: m_data_ff[Q_idx].ffacts[idx]);
		}
	}

	// sort filtered data by Q
	auto sort_data = [](QVector<qreal>& Qvec, QVector<qreal>& ffvec)
	{
		// sort vectors by Q component
		std::vector<std::size_t> perm = tl2::get_perm(Qvec.size(),
			[&Qvec](std::size_t idx1, std::size_t idx2) -> bool
		{
			return Qvec[idx1] < Qvec[idx2];
		});

		Qvec = tl2::reorder(Qvec, perm);
		ffvec = tl2::reorder(ffvec, perm);
	};

	for(t_size idx = 0; idx < ffact_data_ff.size(); ++idx)
		sort_data(Qs_data_ff[idx], ffact_data_ff[idx]);

	// form factor range
	t_real ffact_min_ff = std::numeric_limits<t_real>::max();
	t_real ffact_max_ff = -ffact_min_ff;

	// how many indices do actually have data?
	t_size num_effective_indices = 0;
	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		if(ffact_data_ff[idx].size() != 0)
			++num_effective_indices;
	}

	// plot form factors per index
	t_size effective_index = 0;
	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		bool enabled = effective_index < enabled_indices.size() ? enabled_indices[effective_index] : true;

		// ignore indices with no data
		if(ffact_data_ff[idx].size() == 0)
			continue;

		QCPCurve *curve = new QCPCurve(m_plot_ff->xAxis, m_plot_ff->yAxis);

		// colour for this form factor
		QPen pen = curve->pen();
		int col[3] = {
			num_effective_indices <= 1 ? 0xff
				: int(std::lerp(1., 0., t_real(effective_index) / t_real(num_effective_indices - 1)) * 255.),
			0x00,
			num_effective_indices <= 1 ? 0x00
				: int(std::lerp(0., 1., t_real(effective_index) / t_real(num_effective_indices - 1)) * 255.),
		};

		//tl2::get_colour<int>(g_colPlot, col);
		const QColor colFull(col[0], col[1], col[2]);
		pen.setColor(colFull);
		pen.setWidthF(2.);

		curve->setPen(pen);
		curve->setLineStyle(QCPCurve::lsLine);
		curve->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 1));
		curve->setAntialiased(true);
		curve->setData(Qs_data_ff[idx], ffact_data_ff[idx]);
		curve->setVisible(enabled);

		if(enabled)
		{
			// form factor range for enabled curves
			auto [min_ffact_iter, max_ffact_iter] = std::minmax_element(ffact_data_ff[idx].begin(), ffact_data_ff[idx].end());
			if(min_ffact_iter != ffact_data_ff[idx].end() && max_ffact_iter != ffact_data_ff[idx].end())
			{
				t_real ffact_range = *max_ffact_iter - *min_ffact_iter;

				ffact_max_ff = std::max<t_real>(ffact_max_ff, *max_ffact_iter + ffact_range*0.05);
				ffact_min_ff = std::min<t_real>(ffact_min_ff, *min_ffact_iter - ffact_range*0.05);
			}
		}

		m_curves_ff.push_back(curve);
		AddFormFactorIndex("#" + tl2::var_to_str(effective_index), colFull, enabled);
		++effective_index;
	}

	// set ranges
	m_plot_ff->xAxis->setRange(m_Q_start_ff->value(), m_Q_end_ff->value());
	m_plot_ff->yAxis->setRange(ffact_min_ff, ffact_max_ff);

	// set font
	m_plot_ff->setFont(font());
	m_plot_ff->xAxis->setLabelFont(font());
	m_plot_ff->yAxis->setLabelFont(font());
	m_plot_ff->xAxis->setTickLabelFont(font());
	m_plot_ff->yAxis->setTickLabelFont(font());

	m_plot_ff->replot();
}



/**
 * calculate the form factors
 */
void FormFactorDlg::CalculateFormFactors()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableFormFactorCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableFormFactorCalculation(false);

	ClearFormFactorPlot(false);

	// get Q coordinates
	t_real Q_start = m_Q_start_ff->value();
	t_real Q_end = m_Q_end_ff->value();

	// keep the scanned Q component in ascending order
	if(Q_start > Q_end)
		std::swap(Q_start, Q_end);

	// get settings
	t_size Q_count = m_num_Q_ff->value();

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested_ff = false;
	m_progress_ff->setMinimum(0);
	m_progress_ff->setMaximum(Q_count);
	m_progress_ff->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count);

	m_data_ff.clear();
	m_data_ff.reserve(Q_count);

	for(t_size Q_idx = 0; Q_idx < Q_count; ++Q_idx)
	{
		auto task = [this, &mtx, &Q_start, &Q_end, Q_idx, Q_count]()
		{
			const t_real Q = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(Q_idx) / t_real(Q_count - 1))
				: Q_start;

			// calculate form factors per index
			FormFactorData data_ff;
			data_ff.momentum = Q;

			std::vector<tl2::ExprParser<t_cplx>> parsers = m_parsers_ff;
			for(tl2::ExprParser<t_cplx>& parser : parsers)
			{
				t_real val_ff = 0.;
				t_real val_ff2 = 0.;

				if(parser)
				{
					parser.register_var("Q", Q);
					parser.register_var("Q2", Q*Q);
					parser.register_var("s", Q / (4.*tl2::pi<t_real>));
					parser.register_var("s2", std::pow(Q / (4.*tl2::pi<t_real>), 2.));

					t_cplx val = parser.eval_noexcept();
					val_ff = val.real();
					val_ff2 = std::norm(val);
				}

				data_ff.ffacts.push_back(val_ff);
				data_ff.ffacts2.push_back(val_ff2);
			}

			std::lock_guard<std::mutex> _lck{mtx};
			m_data_ff.emplace_back(std::move(data_ff));
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

		if(m_stopRequested_ff)
		{
			pool.stop();
			break;
		}

		task->get_future().get();

		if(process_evts || task_idx + 1 == tasks.size())
			m_progress_ff->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested_ff)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort raw unfiltered data by Q
	std::vector<std::size_t> perm_all = tl2::get_perm(m_data_ff.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data_ff[idx1].momentum < m_data_ff[idx2].momentum;
	});

	m_data_ff = tl2::reorder(m_data_ff, perm_all);

	PlotFormFactors(true);
}



/**
 * clears the dispersion graph
 */
void FormFactorDlg::ClearFormFactorPlot(bool replot)
{
	m_curves_ff.clear();

	if(m_plot_ff)
	{
		m_plot_ff->clearPlottables();
		if(replot)
			m_plot_ff->replot();
	}
}



/**
 * show current cursor coordinates
 */
void FormFactorDlg::FormFactorPlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot_ff->xAxis->pixelToCoord(evt->pos().x());
	t_real v = m_plot_ff->yAxis->pixelToCoord(evt->pos().y());

	QString status = m_ff_squared->isChecked()
		? "Q = %1 rlu, |F_M(Q)|^2 = %2."
		: "Q = %1 rlu, F_M(Q) = %2.";
	status = status
		.arg(Q, 0, 'g', g_prec_gui)
		.arg(v, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void FormFactorDlg::FormFactorPlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot_ff)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot_ff->popup(pos);
		evt->accept();
	}
}



/**
 * rescale plot axes to fit the content
 */
void FormFactorDlg::RescaleFormFactorPlot()
{
	if(!m_plot_ff)
		return;

	m_plot_ff->rescaleAxes();
	m_plot_ff->replot();
}



/**
 * save plot as image file
 */
void FormFactorDlg::SaveFormFactorPlotFigure()
{
	if(!m_plot_ff)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("ffact/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("ffact/dir", QFileInfo(filename).path());

	if(!m_plot_ff->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void FormFactorDlg::SaveFormFactorData()
{
	if(m_data_ff.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("ffact/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("ffact/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename).toStdString().c_str());
		return;
	}

	t_size num_indices = m_data_ff[0].ffacts.size();

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
		<< "#\n# Number of form factors: " << num_indices << "\n"
		<< "#\n\n";

	// write column header
	ofstr << std::setw(field_len) << std::left << "# Q" << " ";

	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		std::string ffact = "f_" + tl2::var_to_str(idx);
		ofstr << std::setw(field_len) << std::left << ffact << " ";
	}
	ofstr << "\n";

	// write data
	for(const FormFactorData& data: m_data_ff)
	{
		ofstr << std::setw(field_len) << std::left << data.momentum << " ";

		assert(num_indices == data.ffacts.size());
		for(t_size idx = 0; idx < num_indices; ++idx)
			ofstr << std::setw(field_len) << std::left << data.ffacts[idx] << " ";

		ofstr << "\n";
	}

	ofstr.flush();
}



/**
 * toggle between "calculate" and "stop" button
 */
void FormFactorDlg::EnableFormFactorCalculation(bool enable)
{
	m_calcEnabled_ff = enable;

	if(enable)
	{
		m_btnStartStop_ff->setText("Calculate");
		m_btnStartStop_ff->setToolTip("Start calculation.");
		m_btnStartStop_ff->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop_ff->setText("Stop");
		m_btnStartStop_ff->setToolTip("Stop running calculation.");
		m_btnStartStop_ff->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}

// ============================================================================
