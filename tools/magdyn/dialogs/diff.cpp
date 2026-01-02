/**
 * magnetic dynamics -- differentiation / group velocity calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2025
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

#include "diff.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



// ============================================================================
// differentiation dialog
// ============================================================================

/**
 * sets up the differentiation dialog
 */
DiffDlg::DiffDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Differentiation");
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
	maingrid->addWidget(CreateGroupVelocityPanel(), 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnbox, 1, 3, 1, 1);

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &DiffDlg::accept);
}



DiffDlg::~DiffDlg()
{
}



/**
 * set a pointer to the main magdyn kernel
 */
void DiffDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * set the Q start and end points from the main window's dispersion
 */
void DiffDlg::SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend)
{
	m_Qstart = Qstart;
	m_Qend = Qend;
}



void DiffDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * dialog is closing
 */
void DiffDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("diff/geo", saveGeometry());
		m_sett->setValue("diff/splitter", m_split_plot_gv->saveState());
	}

	QDialog::accept();
}
// ============================================================================



// ============================================================================
// calculate group velocity
// ============================================================================

/**
 * column indices in magnon band table for the group velocity
 */
enum : int
{
	COL_GV_BAND = 0,
	COL_GV_ACTIVE,
	NUM_COLS_GV,
};



/**
 * create the panel for the group velocity tab
 */
QWidget* DiffDlg::CreateGroupVelocityPanel()
{
	QWidget *panelGroupVelocity = new QWidget(this);

	// plotter
	m_plot_gv = new QCustomPlot(panelGroupVelocity);
	m_plot_gv->setFont(font());
	m_plot_gv->xAxis->setLabel("Momentum Transfer Q (rlu)");
	m_plot_gv->yAxis->setLabel("Group Velocity v = dE / dq (meV/rlu)");
	m_plot_gv->setInteraction(QCP::iRangeDrag, true);
	m_plot_gv->setInteraction(QCP::iRangeZoom, true);
	m_plot_gv->setSelectionRectMode(QCP::srmZoom);
	m_plot_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// magnon band table
	QWidget *bands_panel = new QWidget(panelGroupVelocity);
	m_table_bands_gv = new QTableWidget(bands_panel);
	m_table_bands_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_bands_gv->setShowGrid(true);
	m_table_bands_gv->setSortingEnabled(false);
	m_table_bands_gv->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_bands_gv->setSelectionMode(QTableWidget::SingleSelection);
	m_table_bands_gv->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_bands_gv->verticalHeader()->setVisible(false);
	m_table_bands_gv->setColumnCount(NUM_COLS_GV);
	m_table_bands_gv->setHorizontalHeaderItem(COL_GV_BAND, new QTableWidgetItem{"Band"});
	m_table_bands_gv->setHorizontalHeaderItem(COL_GV_ACTIVE, new QTableWidgetItem{"Act."});
	m_table_bands_gv->setColumnWidth(COL_GV_BAND, 40);
	m_table_bands_gv->setColumnWidth(COL_GV_ACTIVE, 25);
	m_table_bands_gv->resizeColumnsToContents();

	m_only_pos_E_gv = new QCheckBox("E ≥ 0", bands_panel);
	m_only_pos_E_gv->setChecked(false);
	m_only_pos_E_gv->setToolTip("Ignore magnon annihilation.");

	// splitter for plot and magnon band list
	m_split_plot_gv = new QSplitter(panelGroupVelocity);
	m_split_plot_gv->setOrientation(Qt::Horizontal);
	m_split_plot_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot_gv->addWidget(m_plot_gv);
	m_split_plot_gv->addWidget(bands_panel);
	m_split_plot_gv->setCollapsible(0, false);
	m_split_plot_gv->setCollapsible(1, true);
	m_split_plot_gv->setStretchFactor(m_split_plot_gv->indexOf(m_plot_gv), 24);
	m_split_plot_gv->setStretchFactor(m_split_plot_gv->indexOf(bands_panel), 1);

	// context menu for plotter
	m_menuPlot_gv = new QMenu("Plotter", panelGroupVelocity);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot_gv);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot_gv);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot_gv);

	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));

	m_menuPlot_gv->addAction(acRescalePlot);
	m_menuPlot_gv->addSeparator();
	m_menuPlot_gv->addAction(acSaveFigure);
	m_menuPlot_gv->addAction(acSaveData);

	// bands panel grid
	int y_bands = 0;
	QGridLayout *grid_bands = new QGridLayout(bands_panel);
	grid_bands->setSpacing(4);
	grid_bands->setContentsMargins(6, 6, 6, 6);
	grid_bands->addWidget(m_table_bands_gv, y_bands++, 0, 1, 1);
	grid_bands->addWidget(m_only_pos_E_gv, y_bands++, 0, 1, 1);

	// differentiation order
	m_diff = new QSpinBox(panelGroupVelocity);
	m_diff->setMinimum(1);
	m_diff->setMaximum(9);
	m_diff->setValue(1);
	m_diff->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_diff->setToolTip("Differentiation order, d^n(E)/dq^n.");

	// start and stop coordinates
	m_Q_start_gv[0] = new QDoubleSpinBox(panelGroupVelocity);
	m_Q_start_gv[1] = new QDoubleSpinBox(panelGroupVelocity);
	m_Q_start_gv[2] = new QDoubleSpinBox(panelGroupVelocity);
	m_Q_end_gv[0] = new QDoubleSpinBox(panelGroupVelocity);
	m_Q_end_gv[1] = new QDoubleSpinBox(panelGroupVelocity);
	m_Q_end_gv[2] = new QDoubleSpinBox(panelGroupVelocity);

	m_Q_start_gv[0]->setToolTip("Dispersion initial momentum transfer, h_i (rlu).");
	m_Q_start_gv[1]->setToolTip("Dispersion initial momentum transfer, k_i (rlu).");
	m_Q_start_gv[2]->setToolTip("Dispersion initial momentum transfer, l_i (rlu).");
	m_Q_end_gv[0]->setToolTip("Dispersion final momentum transfer, h_f (rlu).");
	m_Q_end_gv[1]->setToolTip("Dispersion final momentum transfer, k_f (rlu).");
	m_Q_end_gv[2]->setToolTip("Dispersion final momentum transfer, l_f (rlu).");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_start_gv[i]->setDecimals(4);
		m_Q_start_gv[i]->setMinimum(-99.9999);
		m_Q_start_gv[i]->setMaximum(+99.9999);
		m_Q_start_gv[i]->setSingleStep(0.01);
		m_Q_start_gv[i]->setValue(i == 0 ? -1. : 0.);
		m_Q_start_gv[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_start_gv[i]->setPrefix(hklPrefix[i]);

		m_Q_end_gv[i]->setDecimals(4);
		m_Q_end_gv[i]->setMinimum(-99.9999);
		m_Q_end_gv[i]->setMaximum(+99.9999);
		m_Q_end_gv[i]->setSingleStep(0.01);
		m_Q_end_gv[i]->setValue(i == 0 ? 1. : 0.);
		m_Q_end_gv[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_end_gv[i]->setPrefix(hklPrefix[i]);
	}

	// number of Q points in the plot
	m_num_Q_gv = new QSpinBox(panelGroupVelocity);
	m_num_Q_gv->setMinimum(1);
	m_num_Q_gv->setMaximum(99999);
	m_num_Q_gv->setValue(128);
	m_num_Q_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_Q_gv->setToolTip("Number of Q points to calculate.");

	// dispersion Q button
	QPushButton *btnQ = new QPushButton("Set Main Q", panelGroupVelocity);
	btnQ->setToolTip("Set the Q start and end points from the dispersion in the main window.");

	// maximum cutoff for filtering numerical artefacts in group velocity
	m_v_filter_enable_gv = new QCheckBox("Maximum v:", panelGroupVelocity);
	m_v_filter_enable_gv->setChecked(true);
	m_v_filter_enable_gv->setToolTip("Enable maximum cutoff group velocity for filtering numerical artefacts.");

	m_v_filter_gv = new QDoubleSpinBox(panelGroupVelocity);
	m_v_filter_gv->setDecimals(2);
	m_v_filter_gv->setMinimum(0.);
	m_v_filter_gv->setMaximum(999999.99);
	m_v_filter_gv->setSingleStep(1.);
	m_v_filter_gv->setValue(m_v_filter_gv->maximum());
	m_v_filter_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_v_filter_gv->setToolTip("Maximum cutoff group velocity for filtering numerical artefacts.");

	// minimum cutoff for filtering S(Q, E)
	m_S_filter_enable_gv = new QCheckBox("Minimum S(Q, E):", panelGroupVelocity);
	m_S_filter_enable_gv->setChecked(false);
	m_S_filter_enable_gv->setToolTip("Enable minimum S(Q, E).");

	m_S_filter_gv = new QDoubleSpinBox(panelGroupVelocity);
	m_S_filter_gv->setDecimals(5);
	m_S_filter_gv->setMinimum(0.);
	m_S_filter_gv->setMaximum(9999.99999);
	m_S_filter_gv->setSingleStep(0.01);
	m_S_filter_gv->setValue(0.01);
	m_S_filter_gv->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_S_filter_gv->setToolTip("Minimum S(Q, E) to keep.");

	// progress bar
	m_progress_gv = new QProgressBar(panelGroupVelocity);
	m_progress_gv->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btnStartStop_gv = new QPushButton("Calculate", panelGroupVelocity);

	// component grid
	auto grid = new QGridLayout(panelGroupVelocity);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_split_plot_gv, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Differentiation Order:", panelGroupVelocity), y, 0, 1, 1);
	grid->addWidget(m_diff, y++, 1, 1, 1);
	grid->addWidget(new QLabel("Start Q (rlu):", panelGroupVelocity), y, 0, 1, 1);
	grid->addWidget(m_Q_start_gv[0], y, 1, 1, 1);
	grid->addWidget(m_Q_start_gv[1], y, 2, 1, 1);
	grid->addWidget(m_Q_start_gv[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("End Q (rlu):", panelGroupVelocity), y, 0, 1, 1);
	grid->addWidget(m_Q_end_gv[0], y, 1, 1, 1);
	grid->addWidget(m_Q_end_gv[1], y, 2, 1, 1);
	grid->addWidget(m_Q_end_gv[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Count:", panelGroupVelocity), y, 0, 1, 1);
	grid->addWidget(m_num_Q_gv, y, 1, 1, 1);
	grid->addWidget(btnQ, y++, 3, 1, 1);
	grid->addWidget(m_v_filter_enable_gv, y, 0, 1, 1);
	grid->addWidget(m_v_filter_gv, y, 1, 1, 1);
	grid->addWidget(m_S_filter_enable_gv, y, 2, 1, 1);
	grid->addWidget(m_S_filter_gv, y++, 3, 1, 1);
	grid->addWidget(m_progress_gv, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop_gv, y++, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("diff/geo"))
			restoreGeometry(m_sett->value("diff/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("diff/splitter"))
			m_split_plot_gv->restoreState(m_sett->value("diff/splitter").toByteArray());
	}

	// connections
	connect(m_plot_gv, &QCustomPlot::mouseMove, this, &DiffDlg::GroupVelocityPlotMouseMove);
	connect(m_plot_gv, &QCustomPlot::mousePress, this, &DiffDlg::GroupVelocityPlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &DiffDlg::RescaleGroupVelocityPlot);
	connect(acSaveFigure, &QAction::triggered, this, &DiffDlg::SaveGroupVelocityPlotFigure);
	connect(acSaveData, &QAction::triggered, this, &DiffDlg::SaveGroupVelocityData);
	connect(btnQ, &QAbstractButton::clicked, this, &DiffDlg::SetGroupVelocityQ);
	connect(m_v_filter_enable_gv, &QCheckBox::toggled, m_v_filter_gv, &QDoubleSpinBox::setEnabled);
	connect(m_S_filter_enable_gv, &QCheckBox::toggled, m_S_filter_gv, &QDoubleSpinBox::setEnabled);

	// calculation
	connect(m_btnStartStop_gv, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled_gv)
			CalculateGroupVelocity();
		else
			m_stopRequested_gv = true;
	});

	// replotting
	connect(m_only_pos_E_gv, &QCheckBox::toggled, [this]() { PlotGroupVelocity(); });
	connect(m_v_filter_enable_gv, &QCheckBox::toggled, [this]() { PlotGroupVelocity(); });
	connect(m_S_filter_enable_gv, &QCheckBox::toggled, [this]() { PlotGroupVelocity(); });
	connect(m_v_filter_gv, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this]() { PlotGroupVelocity(); });
	connect(m_S_filter_gv, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this]() { PlotGroupVelocity(); });
	connect(m_diff, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this]() { PlotGroupVelocity(); });

	m_v_filter_gv->setEnabled(m_v_filter_enable_gv->isChecked());
	m_S_filter_gv->setEnabled(m_S_filter_enable_gv->isChecked());
	EnableGroupVelocityCalculation();

	return panelGroupVelocity;
}



/**
 * clears the table of magnon bands
 */
void DiffDlg::ClearGroupVelocityBands()
{
	m_table_bands_gv->clearContents();
	m_table_bands_gv->setRowCount(0);
}



/**
 * adds a magnon band to the table
 */
void DiffDlg::AddGroupVelocityBand(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_bands_gv)
		return;

	int row = m_table_bands_gv->rowCount();
	m_table_bands_gv->insertRow(row);

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

	QCheckBox *checkBand = new QCheckBox(m_table_bands_gv);
	checkBand->setChecked(enabled);
	connect(checkBand, &QCheckBox::toggled, [this]() { PlotGroupVelocity(false); });

	m_table_bands_gv->setItem(row, COL_GV_BAND, item);
	m_table_bands_gv->setCellWidget(row, COL_GV_ACTIVE, checkBand);
}



/**
 * verifies if the band's checkbox is checked
 */
bool DiffDlg::IsGroupVelocityBandEnabled(t_size idx) const
{
	if(!m_table_bands_gv || int(idx) >= m_table_bands_gv->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(m_table_bands_gv->cellWidget(int(idx), COL_GV_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * calculate the filtered data sets and plot the group velocities
 */
void DiffDlg::PlotGroupVelocity(bool clear_settings)
{
	if(!m_plot_gv)
		return;

	// keep some settings from previous plot, e.g. the band visibility flags
	std::vector<bool> enabled_bands;
	if(!clear_settings)
	{
		enabled_bands.reserve(m_table_bands_gv->rowCount());
		for(int row = 0; row < m_table_bands_gv->rowCount(); ++row)
			enabled_bands.push_back(IsGroupVelocityBandEnabled(t_size(row)));
	}

	ClearGroupVelocityBands();
	ClearGroupVelocityPlot(false);

	int diff_order = m_diff->value();
	if(diff_order == 1)
		m_plot_gv->yAxis->setLabel("Group Velocity v = dE / dq (meV/rlu)");
	else
		m_plot_gv->yAxis->setLabel(QString("d^%1E / dq^%1").arg(diff_order));

	if(m_data_gv.size() == 0)
	{
		m_plot_gv->replot();
		return;
	}

	// get settings
	t_real max_v = m_v_filter_gv->value();
	if(!m_v_filter_enable_gv->isChecked())
		max_v = -1.;  // disable B filter

	t_real min_S = m_S_filter_gv->value();
	if(!m_S_filter_enable_gv->isChecked())
		min_S = -1.;  // disable S(Q, E) filter

	bool only_creation = m_only_pos_E_gv->isChecked();

	t_size num_Q = m_data_gv.size();
	t_size num_bands = m_data_gv[0].velocities.size();

	// filtered momentum transfer and group velocity per band
	std::vector<QVector<qreal>> Qs_data_gv{num_bands};
	std::vector<QVector<qreal>> vs_data_gv{num_bands};

	for(t_size Q_idx = 0; Q_idx < num_Q; ++Q_idx)
	{
		const t_vec_real& Q = m_data_gv[Q_idx].momentum;

		for(t_size band = 0; band < num_bands; ++band)
		{
			const t_real& v = m_data_gv[Q_idx].velocities[band];

			// filter numerical artefacts in v
			if(max_v >= 0. && std::abs(v) > max_v)
				continue;

			// filter minimum S(Q, E)
			if(min_S >= 0. && std::abs(m_data_gv[Q_idx].weights[band]) <= min_S)
				continue;

			// filter magnon annihilation
			if(only_creation && m_data_gv[Q_idx].energies[band] < 0.)
				continue;

			Qs_data_gv[band].push_back(Q[m_Q_idx_gv]);
			vs_data_gv[band].push_back(v);
		}
	}

	// calculate higher-order differentials if requested
	for(int order = 1; order < diff_order; ++order)
	{
		// differentiate all bands
		for(t_size band = 0; band < num_bands; ++band)
			vs_data_gv[band] = tl2::diff<QVector<qreal>>(Qs_data_gv[band], vs_data_gv[band]);
	}

	// sort filtered data by Q
	auto sort_data = [](QVector<qreal>& Qvec, QVector<qreal>& vvec)
	{
		// sort vectors by Q component
		std::vector<std::size_t> perm = tl2::get_perm(Qvec.size(),
			[&Qvec](std::size_t idx1, std::size_t idx2) -> bool
		{
			return Qvec[idx1] < Qvec[idx2];
		});

		Qvec = tl2::reorder(Qvec, perm);
		vvec = tl2::reorder(vvec, perm);
	};

	for(t_size band = 0; band < vs_data_gv.size(); ++band)
		sort_data(Qs_data_gv[band], vs_data_gv[band]);

	// group velocity range
	t_real v_min_gv = std::numeric_limits<t_real>::max();
	t_real v_max_gv = -v_min_gv;

	// how many bands do actually have data?
	t_size num_effective_bands = 0;
	for(t_size band = 0; band < num_bands; ++band)
	{
		if(vs_data_gv[band].size() != 0)
			++num_effective_bands;
	}

	// plot group velocites per band
	t_size effective_band = 0;
	for(t_size band = 0; band < num_bands; ++band)
	{
		bool enabled = effective_band < enabled_bands.size() ? enabled_bands[effective_band] : true;

		// ignore bands with no data
		if(vs_data_gv[band].size() == 0)
			continue;

		QCPCurve *curve = new QCPCurve(m_plot_gv->xAxis, m_plot_gv->yAxis);

		// colour for this magnon band
		QPen pen = curve->pen();
		int col[3] = {
			num_effective_bands <= 1 ? 0xff
				: int(std::lerp(1., 0., t_real(effective_band) / t_real(num_effective_bands - 1)) * 255.),
			0x00,
			num_effective_bands <= 1 ? 0x00
				: int(std::lerp(0., 1., t_real(effective_band) / t_real(num_effective_bands - 1)) * 255.),
		};

		//tl2::get_colour<int>(g_colPlot, col);
		const QColor colFull(col[0], col[1], col[2]);
		pen.setColor(colFull);
		pen.setWidthF(2.);

		curve->setPen(pen);
		curve->setLineStyle(QCPCurve::lsLine);
		curve->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 1));
		curve->setAntialiased(true);
		curve->setData(Qs_data_gv[band], vs_data_gv[band]);
		curve->setVisible(enabled);

		if(enabled)
		{
			// group velocity range for enabled curves
			auto [min_v_iter, max_v_iter] = std::minmax_element(vs_data_gv[band].begin(), vs_data_gv[band].end());
			if(min_v_iter != vs_data_gv[band].end() && max_v_iter != vs_data_gv[band].end())
			{
				t_real v_range = *max_v_iter - *min_v_iter;

				v_max_gv = std::max(v_max_gv, *max_v_iter + v_range*0.05);
				v_min_gv = std::min(v_min_gv, *min_v_iter - v_range*0.05);
			}
		}

		m_curves_gv.push_back(curve);
		AddGroupVelocityBand("#" + tl2::var_to_str(effective_band + 1), colFull, enabled);
		++effective_band;
	}

	// set labels
	const char* Q_label[]{ "h (rlu)", "k (rlu)", "l (rlu)" };
	m_plot_gv->xAxis->setLabel(QString("Momentum Transfer ") + Q_label[m_Q_idx_gv]);

	// set ranges
	m_plot_gv->xAxis->setRange(m_Q_min_gv, m_Q_max_gv);
	m_plot_gv->yAxis->setRange(v_min_gv, v_max_gv);

	// set font
	m_plot_gv->setFont(font());
	m_plot_gv->xAxis->setLabelFont(font());
	m_plot_gv->yAxis->setLabelFont(font());
	m_plot_gv->xAxis->setTickLabelFont(font());
	m_plot_gv->yAxis->setTickLabelFont(font());

	m_plot_gv->replot();
}



/**
 * calculate the group velocity
 */
void DiffDlg::CalculateGroupVelocity()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableGroupVelocityCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableGroupVelocityCalculation(false);

	ClearGroupVelocityPlot(false);

	// get coordinates
	t_vec_real Q_start = tl2::create<t_vec_real>(
	{
		m_Q_start_gv[0]->value(),
		m_Q_start_gv[1]->value(),
		m_Q_start_gv[2]->value(),
	});

	t_vec_real Q_end = tl2::create<t_vec_real>(
	{
		m_Q_end_gv[0]->value(),
		m_Q_end_gv[1]->value(),
		m_Q_end_gv[2]->value(),
	});

	// get Q component with maximum range
	t_vec_real Q_range = Q_end - Q_start;
	m_Q_idx_gv = 0;
	if(std::abs(Q_range[1]) > std::abs(Q_range[m_Q_idx_gv]))
		m_Q_idx_gv = 1;
	if(std::abs(Q_range[2]) > std::abs(Q_range[m_Q_idx_gv]))
		m_Q_idx_gv = 2;

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx_gv] > Q_end[m_Q_idx_gv])
		std::swap(Q_start, Q_end);

	// Q range
	m_Q_min_gv = Q_start[m_Q_idx_gv];
	m_Q_max_gv = Q_end[m_Q_idx_gv];

	// get settings
	t_size Q_count = m_num_Q_gv->value();
	std::vector<t_size> *perm = nullptr;

	// calculate group velocity
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested_gv = false;
	m_progress_gv->setMinimum(0);
	m_progress_gv->setMaximum(Q_count);
	m_progress_gv->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count);

	m_data_gv.clear();
	m_data_gv.reserve(Q_count);

	for(t_size Q_idx = 0; Q_idx < Q_count; ++Q_idx)
	{
		auto task = [this, &mtx, &dyn, &Q_start, &Q_end, Q_idx, Q_count, perm]()
		{
			const t_vec_real Q1 = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(Q_idx) / t_real(Q_count - 1))
				: Q_start;
			const t_vec_real Q2 = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(Q_idx + 1) / t_real(Q_count - 1))
				: Q_end;

			// calculate group velocities per band
			GroupVelocityData data_gv;
			data_gv.momentum = Q1;
			typename t_magdyn::SofQE S;
			std::tie(data_gv.velocities, S) = dyn.CalcGroupVelocities(Q1, Q2 - Q1, perm);
			t_size num_bands = data_gv.velocities.size();
			data_gv.energies.reserve(num_bands);
			data_gv.weights.reserve(num_bands);

			// calculate energies per band
			//assert(S.E_and_S.size() == num_bands);
			if(S.E_and_S.size() != num_bands)
				return;

			for(t_size band = 0; band < num_bands; ++band)
			{
				data_gv.energies.push_back(S.E_and_S[band].E);
				data_gv.weights.push_back(S.E_and_S[band].weight);
			}

			std::lock_guard<std::mutex> _lck{mtx};
			m_data_gv.emplace_back(std::move(data_gv));
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

		if(m_stopRequested_gv)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress_gv->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested_gv)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort raw unfiltered data by Q
	std::vector<std::size_t> perm_all = tl2::get_perm(m_data_gv.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data_gv[idx1].momentum[m_Q_idx_gv]
			< m_data_gv[idx2].momentum[m_Q_idx_gv];
	});

	m_data_gv = tl2::reorder(m_data_gv, perm_all);

	PlotGroupVelocity(true);
}



/**
 * clears the dispersion graph
 */
void DiffDlg::ClearGroupVelocityPlot(bool replot)
{
	m_curves_gv.clear();

	if(m_plot_gv)
	{
		m_plot_gv->clearPlottables();
		if(replot)
			m_plot_gv->replot();
	}
}



/**
 * show current cursor coordinates
 */
void DiffDlg::GroupVelocityPlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot_gv->xAxis->pixelToCoord(evt->pos().x());
	t_real v = m_plot_gv->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 rlu, B = %2.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(v, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void DiffDlg::GroupVelocityPlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot_gv)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot_gv->popup(pos);
		evt->accept();
	}
}



/**
 * rescale plot axes to fit the content
 */
void DiffDlg::RescaleGroupVelocityPlot()
{
	if(!m_plot_gv)
		return;

	m_plot_gv->rescaleAxes();
	m_plot_gv->replot();
}



/**
 * save plot as image file
 */
void DiffDlg::SaveGroupVelocityPlotFigure()
{
	if(!m_plot_gv)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("diff/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("diff/dir", QFileInfo(filename).path());

	if(!m_plot_gv->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void DiffDlg::SaveGroupVelocityData()
{
	if(m_data_gv.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("diff/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("diff/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename).toStdString().c_str());
		return;
	}

	t_size num_bands = m_data_gv[0].velocities.size();

	ofstr.precision(g_prec);
	int field_len = g_prec * 2.5;

	// write meta header
	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr << "#\n"
		<< "# Created by Magpie\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n# Number of energy bands: " << num_bands << "\n"
		<< "#\n\n";

	// write column header
	ofstr << std::setw(field_len) << std::left << "# h" << " ";
	ofstr << std::setw(field_len) << std::left << "k" << " ";
	ofstr << std::setw(field_len) << std::left << "l" << " ";

	for(t_size band = 0; band < num_bands; ++band)
	{
		std::string E = "E_" + tl2::var_to_str(band);
		std::string S = "Sperp_" + tl2::var_to_str(band);
		std::string v = "v_" + tl2::var_to_str(band);

		ofstr << std::setw(field_len) << std::left << E << " ";
		ofstr << std::setw(field_len) << std::left << S << " ";
		ofstr << std::setw(field_len) << std::left << v << " ";
	}
	ofstr << "\n";

	// write data
	for(const GroupVelocityData& data: m_data_gv)
	{
		ofstr << std::setw(field_len) << std::left << data.momentum[0] << " ";
		ofstr << std::setw(field_len) << std::left << data.momentum[1] << " ";
		ofstr << std::setw(field_len) << std::left << data.momentum[2] << " ";

		assert(num_bands == data.velocities.size());
		for(t_size band = 0; band < num_bands; ++band)
		{
			ofstr << std::setw(field_len) << std::left << data.energies[band] << " ";
			ofstr << std::setw(field_len) << std::left << data.weights[band] << " ";
			ofstr << std::setw(field_len) << std::left << data.velocities[band] << " ";
		}

		ofstr << "\n";
	}

	ofstr.flush();
}



/**
 * toggle between "calculate" and "stop" button
 */
void DiffDlg::EnableGroupVelocityCalculation(bool enable)
{
	m_calcEnabled_gv = enable;

	if(enable)
	{
		m_btnStartStop_gv->setText("Calculate");
		m_btnStartStop_gv->setToolTip("Start calculation.");
		m_btnStartStop_gv->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop_gv->setText("Stop");
		m_btnStartStop_gv->setToolTip("Stop running calculation.");
		m_btnStartStop_gv->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * set the group velocity's Q positions to the main window dispersion Qs
 */
void DiffDlg::SetGroupVelocityQ()
{
	if(m_Qstart.size() < 3 || m_Qend.size() < 3)
		return;

	for(int i = 0; i < 3; ++i)
	{
		m_Q_start_gv[i]->setValue(m_Qstart[i]);
		m_Q_end_gv[i]->setValue(m_Qend[i]);
	}
}
// ============================================================================
