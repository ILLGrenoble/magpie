/**
 * 2d plotter
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

#include "plot2d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/maths.h"



// ============================================================================
// 2d plotter dialog
// ============================================================================

/**
 * sets up the plot dialog
 */
Plot2DDlg::Plot2DDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("2D Plotter");
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
	maingrid->addWidget(CreatePanel(), 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnbox, 1, 3, 1, 1);

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &Plot2DDlg::accept);
}



Plot2DDlg::~Plot2DDlg()
{
}



void Plot2DDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * dialog is closing
 */
void Plot2DDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("plot2d/geo", saveGeometry());
		m_sett->setValue("plot2d/splitter", m_split_plot->saveState());
	}

	QDialog::accept();
}
// ============================================================================



// ============================================================================
// calculate formulas
// ============================================================================

/**
 * parse the formulas and plot them
 */
void Plot2DDlg::Parse()
{
	m_parsers.clear();

	// get formulas
	QString _formulas = m_formulas->toPlainText();
	std::string all_formulas = _formulas.toStdString();

	std::vector<std::string> formulas;
	tl2::get_tokens<std::string, std::string>(all_formulas, ";", formulas);
	if(formulas.size() == 0)
		return;

	for(t_size idx = 0; idx < formulas.size(); ++idx)
	{
		const std::string& formula = formulas[idx];
		
		tl2::ExprParser<t_real> parser;
		parser.SetAutoregisterVariables(false);
		parser.SetInvalid0(false);
		parser.register_var("x", 0.);

		if(!parser.parse_noexcept(formula))
		{
			std::cerr << "Error: Formula #" << idx
				<< ": \"" << formula << "\" could not be parsed."
				<< std::endl;
		}

		m_parsers.emplace_back(std::move(parser));
	}

	Calculate();
}



/**
 * column indices in table for the plot curves
 */
enum : int
{
	COL_FF_INDEX = 0,
	COL_FF_ACTIVE,
	NUM_COLS_FF,
};



/**
 * create the main panel
 */
QWidget* Plot2DDlg::CreatePanel()
{
	QWidget *panelMain = new QWidget(this);

	// plotter
	m_plot = new QCustomPlot(panelMain);
	m_plot->setFont(font());
	m_plot->xAxis->setLabel("x");
	m_plot->yAxis->setLabel("y");
	m_plot->setInteraction(QCP::iRangeDrag, true);
	m_plot->setInteraction(QCP::iRangeZoom, true);
	m_plot->setSelectionRectMode(QCP::srmZoom);
	m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// plot curves
	QWidget *panel_indices = new QWidget(panelMain);
	m_table = new QTableWidget(panel_indices);
	m_table->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table->setShowGrid(true);
	m_table->setSortingEnabled(false);
	m_table->setSelectionBehavior(QTableWidget::SelectRows);
	m_table->setSelectionMode(QTableWidget::SingleSelection);
	m_table->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table->verticalHeader()->setVisible(false);
	m_table->setColumnCount(NUM_COLS_FF);
	m_table->setHorizontalHeaderItem(COL_FF_INDEX, new QTableWidgetItem{"Curve"});
	m_table->setHorizontalHeaderItem(COL_FF_ACTIVE, new QTableWidgetItem{"Act."});
	m_table->setColumnWidth(COL_FF_INDEX, 40);
	m_table->setColumnWidth(COL_FF_ACTIVE, 25);
	m_table->resizeColumnsToContents();

	// splitter for plot and indices list
	m_split_plot = new QSplitter(panelMain);
	m_split_plot->setOrientation(Qt::Horizontal);
	m_split_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot->addWidget(m_plot);
	m_split_plot->addWidget(panel_indices);
	m_split_plot->setCollapsible(0, false);
	m_split_plot->setCollapsible(1, true);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_plot), 24);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(panel_indices), 1);

	// context menu for plotter
	m_menuPlot = new QMenu("Plotter", panelMain);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot);

	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));

	m_menuPlot->addAction(acRescalePlot);
	m_menuPlot->addSeparator();
	m_menuPlot->addAction(acSaveFigure);
	m_menuPlot->addAction(acSaveData);

	// indices panel grid
	int y_indices = 0;
	QGridLayout *grid_indices = new QGridLayout(panel_indices);
	grid_indices->setSpacing(4);
	grid_indices->setContentsMargins(6, 6, 6, 6);
	grid_indices->addWidget(m_table, y_indices++, 0, 1, 1);

	// formulas field
	m_formulas = new QTextEdit(panelMain);
	m_formulas->setPlaceholderText("Enter formulas separated by ';'."
		" The free variable is 'x'.");
	m_formulas->setReadOnly(false);
	m_formulas->setWordWrapMode(QTextOption::WrapAtWordBoundaryOrAnywhere);
	m_formulas->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	// start and stop coordinates
	m_x_start = new QDoubleSpinBox(panelMain);
	m_x_end = new QDoubleSpinBox(panelMain);

	m_x_start->setToolTip("Initial x value.");
	m_x_end->setToolTip("Final x value.");

	m_x_start->setDecimals(4);
	m_x_start->setMinimum(-9999.9999);
	m_x_start->setMaximum(+9999.9999);
	m_x_start->setSingleStep(0.5);
	m_x_start->setValue(0.);
	m_x_start->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	m_x_end->setDecimals(4);
	m_x_end->setMinimum(-9999.9999);
	m_x_end->setMaximum(+9999.9999);
	m_x_end->setSingleStep(0.5);
	m_x_end->setValue(5.);
	m_x_end->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	// number of x points in the plot
	m_num_x = new QSpinBox(panelMain);
	m_num_x->setMinimum(1);
	m_num_x->setMaximum(99999);
	m_num_x->setValue(128);
	m_num_x->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_x->setToolTip("Number of coordinate points to calculate.");

	// progress bar
	m_progress = new QProgressBar(panelMain);
	m_progress->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btnStartStop = new QPushButton("Calculate", panelMain);

	// component grid
	auto grid = new QGridLayout(panelMain);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_split_plot, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Formulas, f_i(x):", panelMain), y++, 0, 1, 4);
	grid->addWidget(m_formulas, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Start x:", panelMain), y, 0, 1, 1);
	grid->addWidget(m_x_start, y, 1, 1, 1);
	grid->addWidget(new QLabel("End x:", panelMain), y, 2, 1, 1);
	grid->addWidget(m_x_end, y++, 3, 1, 1);
	grid->addWidget(new QLabel("x Count:", panelMain), y, 0, 1, 1);
	grid->addWidget(m_num_x, y++, 1, 1, 1);
	grid->addWidget(m_progress, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop, y++, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("plot2d/geo"))
			restoreGeometry(m_sett->value("plot2d/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("plot2d/splitter"))
			m_split_plot->restoreState(m_sett->value("plot2d/splitter").toByteArray());
	}

	// connections
	connect(m_plot, &QCustomPlot::mouseMove, this, &Plot2DDlg::PlotMouseMove);
	connect(m_plot, &QCustomPlot::mousePress, this, &Plot2DDlg::PlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &Plot2DDlg::RescalePlot);
	connect(acSaveFigure, &QAction::triggered, this, &Plot2DDlg::SavePlotFigure);
	connect(acSaveData, &QAction::triggered, this, &Plot2DDlg::SaveData);

	// calculation
	connect(m_btnStartStop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled)
			Parse();
		else
			m_stopRequested = true;
	});

	EnableCalculation();

	return panelMain;
}



/**
 * clears the table of plot curve indices
 */
void Plot2DDlg::ClearIndices()
{
	m_table->clearContents();
	m_table->setRowCount(0);
}



/**
 * adds a plot curve index to the table
 */
void Plot2DDlg::AddIndex(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table)
		return;

	int row = m_table->rowCount();
	m_table->insertRow(row);

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

	QCheckBox *checkIndex = new QCheckBox(m_table);
	checkIndex->setChecked(enabled);
	connect(checkIndex, &QCheckBox::toggled, [this]() { Plot(false); });

	m_table->setItem(row, COL_FF_INDEX, item);
	m_table->setCellWidget(row, COL_FF_ACTIVE, checkIndex);
}



/**
 * verifies if the index' checkbox is checked
 */
bool Plot2DDlg::IsIndexEnabled(t_size idx) const
{
	if(!m_table || int(idx) >= m_table->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(m_table->cellWidget(int(idx), COL_FF_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * calculate the data sets and plot the formulas
 */
void Plot2DDlg::Plot(bool clear_settings)
{
	if(!m_plot)
		return;

	// keep some settings from previous plot, e.g. the index visibility flags
	std::vector<bool> enabled_indices;
	if(!clear_settings)
	{
		enabled_indices.reserve(m_table->rowCount());
		for(int row = 0; row < m_table->rowCount(); ++row)
			enabled_indices.push_back(IsIndexEnabled(t_size(row)));
	}

	ClearIndices();
	ClearPlot(false);

	if(m_data.size() == 0)
	{
		m_plot->replot();
		return;
	}

	t_size num_x = m_data.size();
	t_size num_indices = m_data[0].ys.size();

	std::vector<QVector<qreal>> xs_data{num_indices};
	std::vector<QVector<qreal>> ys_data{num_indices};

	for(t_size x_idx = 0; x_idx < num_x; ++x_idx)
	{
		t_real x = m_data[x_idx].x;

		for(t_size idx = 0; idx < num_indices; ++idx)
		{
			xs_data[idx].push_back(x);
			ys_data[idx].push_back(m_data[x_idx].ys[idx]);
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

	for(t_size idx = 0; idx < ys_data.size(); ++idx)
		sort_data(xs_data[idx], ys_data[idx]);

	// y range
	t_real y_min = std::numeric_limits<t_real>::max();
	t_real y_max = -y_min;

	// how many indices do actually have data?
	t_size num_effective_indices = 0;
	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		if(ys_data[idx].size() != 0)
			++num_effective_indices;
	}

	// plot curves per index
	t_size effective_index = 0;
	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		bool enabled = effective_index < enabled_indices.size() ? enabled_indices[effective_index] : true;

		// ignore indices with no data
		if(ys_data[idx].size() == 0)
			continue;

		QCPCurve *curve = new QCPCurve(m_plot->xAxis, m_plot->yAxis);

		// colour for this curve
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
		curve->setData(xs_data[idx], ys_data[idx]);
		curve->setVisible(enabled);

		if(enabled)
		{
			// y range for enabled curves
			auto [min_y_iter, max_y_iter] = std::minmax_element(ys_data[idx].begin(), ys_data[idx].end());
			if(min_y_iter != ys_data[idx].end() && max_y_iter != ys_data[idx].end())
			{
				t_real y_range = *max_y_iter - *min_y_iter;

				y_max = std::max<t_real>(y_max, *max_y_iter + y_range*0.05);
				y_min = std::min<t_real>(y_min, *min_y_iter - y_range*0.05);
			}
		}

		m_curves.push_back(curve);
		AddIndex("#" + tl2::var_to_str(effective_index + 1), colFull, enabled);
		++effective_index;
	}

	// set ranges
	m_plot->xAxis->setRange(m_x_start->value(), m_x_end->value());
	m_plot->yAxis->setRange(y_min, y_max);

	// set font
	m_plot->setFont(font());
	m_plot->xAxis->setLabelFont(font());
	m_plot->yAxis->setLabelFont(font());
	m_plot->xAxis->setTickLabelFont(font());
	m_plot->yAxis->setTickLabelFont(font());

	m_plot->replot();
}



/**
 * calculate the curves
 */
void Plot2DDlg::Calculate()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	ClearPlot(false);

	// get x coordinates
	t_real x_start = m_x_start->value();
	t_real x_end = m_x_end->value();

	// keep the x component in ascending order
	if(x_start > x_end)
		std::swap(x_start, x_end);

	// get settings
	t_size x_count = m_num_x->value();

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(x_count);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(x_count);

	m_data.clear();
	m_data.reserve(x_count);

	for(t_size x_idx = 0; x_idx < x_count; ++x_idx)
	{
		auto task = [this, &mtx, &x_start, &x_end, x_idx, x_count]()
		{
			const t_real x = x_count > 1
				? tl2::lerp(x_start, x_end, t_real(x_idx) / t_real(x_count - 1))
				: x_start;

			// calculate curves per index
			PlotData data;
			data.x = x;

			std::vector<tl2::ExprParser<t_real>> parsers = m_parsers;
			for(tl2::ExprParser<t_real>& parser : parsers)
			{
				t_real val = 0.;

				if(parser)
				{
					parser.register_var("x", x);
					val = parser.eval_noexcept();
				}

				data.ys.push_back(val);
			}

			std::lock_guard<std::mutex> _lck{mtx};
			m_data.emplace_back(std::move(data));
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
		if(task_idx % std::max<t_size>(tasks.size() / g_stop_check_fraction, 1) == 0)
			qApp->processEvents();

		if(m_stopRequested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
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

	// sort raw unfiltered data by Q
	std::vector<std::size_t> perm_all = tl2::get_perm(m_data.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data[idx1].x < m_data[idx2].x;
	});

	m_data = tl2::reorder(m_data, perm_all);

	Plot(true);
}



/**
 * clears the dispersion graph
 */
void Plot2DDlg::ClearPlot(bool replot)
{
	m_curves.clear();

	if(m_plot)
	{
		m_plot->clearPlottables();
		if(replot)
			m_plot->replot();
	}
}



/**
 * show current cursor coordinates
 */
void Plot2DDlg::PlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot->xAxis->pixelToCoord(evt->pos().x());
	t_real v = m_plot->yAxis->pixelToCoord(evt->pos().y());

	QString status("x = %1, y = %2.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(v, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void Plot2DDlg::PlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot->popup(pos);
		evt->accept();
	}
}



/**
 * rescale plot axes to fit the content
 */
void Plot2DDlg::RescalePlot()
{
	if(!m_plot)
		return;

	m_plot->rescaleAxes();
	m_plot->replot();
}



/**
 * save plot as image file
 */
void Plot2DDlg::SavePlotFigure()
{
	if(!m_plot)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("plot2d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("plot2d/dir", QFileInfo(filename).path());

	if(!m_plot->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void Plot2DDlg::SaveData()
{
	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("plot2d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("plot2d/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename).toStdString().c_str());
		return;
	}

	t_size num_indices = m_data[0].ys.size();

	ofstr.precision(g_prec);
	int field_len = g_prec * 2.5;

	// write meta header
	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr << "#\n"
		<< "# Created by tlibs\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n# Number of curves: " << num_indices << "\n"
		<< "#\n\n";

	// write column header
	ofstr << std::setw(field_len) << std::left << "# x" << " ";

	for(t_size idx = 0; idx < num_indices; ++idx)
	{
		std::string y = "y_" + tl2::var_to_str(idx);
		ofstr << std::setw(field_len) << std::left << y << " ";
	}
	ofstr << "\n";

	// write data
	for(const PlotData& data: m_data)
	{
		ofstr << std::setw(field_len) << std::left << data.x << " ";

		assert(num_indices == data.ys.size());
		for(t_size idx = 0; idx < num_indices; ++idx)
			ofstr << std::setw(field_len) << std::left << data.ys[idx] << " ";

		ofstr << "\n";
	}

	ofstr.flush();
}



/**
 * toggle between "calculate" and "stop" button
 */
void Plot2DDlg::EnableCalculation(bool enable)
{
	m_calcEnabled = enable;

	if(enable)
	{
		m_btnStartStop->setText("Calculate");
		m_btnStartStop->setToolTip("Start calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop->setText("Stop");
		m_btnStartStop->setToolTip("Stop running calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}

// ============================================================================
