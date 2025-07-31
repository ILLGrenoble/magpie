/**
 * brillouin zone tool -- gui setup
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include "bz.h"

#include <QtCore/QMimeData>
#include <QtWidgets/QApplication>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QStyle>

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;

#include "libs/loadcif.h"


BZDlg::BZDlg(QWidget* pParent) : QDialog{pParent},
	m_sett{new QSettings{"takin", "bz", this}}
{
	setWindowTitle("Brillouin Zones");
	setSizeGripEnabled(true);
	setFont(QFontDatabase::systemFont(QFontDatabase::GeneralFont));

	// restore settings done from takin main settings dialog
	QSettings sett_core("takin", "core");
	if(sett_core.contains("main/font_gen"))
	{
		QString font_str = sett_core.value("main/font_gen").toString();
		QFont font = this->font();
		//font.setPointSize(14);
		if(font.fromString(font_str))
			setFont(font);
	}
	if(sett_core.contains("main/prec"))
	{
		g_prec_bz = sett_core.value("main/prec").toInt();
		g_eps_bz = std::pow(t_real(10), -t_real(g_prec_bz));
	}
	if(sett_core.contains("main/prec_gfx"))
	{
		g_prec_gui_bz = sett_core.value("main/prec_gfx").toInt();
	}

	// create gui
	m_tabs_in = new QTabWidget(this);
	m_tabs_out = new QTabWidget(this);

	// continue creating gui
	CreateSymopsPanel();
	CreateBZPanel();
	CreateFormulasPanel();
	CreateResultsPanel();
	CreateJsonResultsPanel();

	// status bar
	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// splitter for input and output tabs
	m_split_inout = new QSplitter(this);
	m_split_inout->setOrientation(Qt::Horizontal);
	m_split_inout->setChildrenCollapsible(true);
	m_split_inout->addWidget(m_tabs_in);
	m_split_inout->addWidget(m_tabs_out);

	// main grid
	QGridLayout *main_grid = new QGridLayout(this);
	main_grid->setSpacing(4);
	main_grid->setContentsMargins(4,4,4,4);
	main_grid->addWidget(m_split_inout, 0,0, 1,1);
	main_grid->addWidget(m_status, 1,0, 1,1);

	// continue creating gui
	CreateInfoDialog();
	CreateMenuBar(main_grid);

	// restore settings
	if(m_sett)
	{
		// restore window size and position
		if(m_sett->contains("geo"))
			restoreGeometry(m_sett->value("geo").toByteArray());
		else
			resize(800, 600);

		if(m_sett->contains("recent_files"))
			m_recent.SetRecentFiles(m_sett->value("recent_files").toStringList());

		if(m_sett->contains("splitter"))
			m_split_inout->restoreState(m_sett->value("splitter").toByteArray());
	}

	setAcceptDrops(true);

	m_symOpIgnoreChanges = false;
	m_formulaIgnoreChanges = false;
	CalcB(true);

	m_status->setText("Ready.");
}


/**
 * space group and symmetry operations panel
 */
void BZDlg::CreateSymopsPanel()
{
	// already created?
	if(m_symops)
		return;

	QWidget *symopspanel = new QWidget(this);

	m_symops = new QTableWidget(symopspanel);
	m_symops->setShowGrid(true);
	m_symops->setSortingEnabled(true);
	m_symops->setMouseTracking(true);
	m_symops->setSelectionBehavior(QTableWidget::SelectRows);
	m_symops->setSelectionMode(QTableWidget::ContiguousSelection);
	m_symops->setContextMenuPolicy(Qt::CustomContextMenu);
	m_symops->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing()*6 + 4);
	m_symops->verticalHeader()->setVisible(false);
	m_symops->setAlternatingRowColors(true);
	m_symops->setColumnCount(NUM_SYMOP_COLS);
	m_symops->setHorizontalHeaderItem(COL_OP,
		new QTableWidgetItem{"Symmetry Operation"});
	m_symops->setHorizontalHeaderItem(COL_PROP,
		new QTableWidgetItem{"Properties"});
	m_symops->setHorizontalHeaderItem(COL_DET,
		new QTableWidgetItem{"Det."});
	m_symops->setColumnWidth(COL_OP, 150);
	m_symops->setColumnWidth(COL_PROP, 100);
	m_symops->setColumnWidth(COL_DET, 75);

	QToolButton *btnAdd = new QToolButton(symopspanel);
	QToolButton *btnDel = new QToolButton(symopspanel);
	QToolButton *btnUp = new QToolButton(symopspanel);
	QToolButton *btnDown = new QToolButton(symopspanel);
	QPushButton *btnSG = new QPushButton(symopspanel);

	m_symops->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	btnAdd->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnDel->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnUp->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnDown->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnSG->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	btnAdd->setText("Add");
	btnDel->setText("Delete");
	btnUp->setText("Up");
	btnDown->setText("Down");
	btnSG->setText("Get");

	btnAdd->setToolTip("Add a new symmetry operation.");
	btnDel->setToolTip("Delete the selected symmetry operation.");
	btnUp->setToolTip("Move the selected symmetry operation up.");
	btnDown->setToolTip("Move the selected symmetry operation down.");
	btnSG->setToolTip("Get the symmetry operations from the selected space group.");

	btnAdd->setIcon(QIcon::fromTheme("list-add"));
	btnDel->setIcon(QIcon::fromTheme("list-remove"));
	btnUp->setIcon(QIcon::fromTheme("go-up"));
	btnDown->setIcon(QIcon::fromTheme("go-down"));

	m_editA = new QDoubleSpinBox(symopspanel);
	m_editB = new QDoubleSpinBox(symopspanel);
	m_editC = new QDoubleSpinBox(symopspanel);
	m_editAlpha = new QDoubleSpinBox(symopspanel);
	m_editBeta = new QDoubleSpinBox(symopspanel);
	m_editGamma = new QDoubleSpinBox(symopspanel);

	for(QDoubleSpinBox* spin : {m_editA, m_editB, m_editC})
	{
		spin->setMinimum(0.);
		spin->setMaximum(999.);
		spin->setValue(5.);
		spin->setDecimals(3);
		spin->setSingleStep(0.1);
	}

	for(QDoubleSpinBox* spin : {m_editAlpha, m_editBeta, m_editGamma})
	{
		spin->setMinimum(0.);
		spin->setMaximum(360.);
		spin->setValue(90.);
		spin->setDecimals(2);
		spin->setSingleStep(1.);
	}

	m_comboSG = new QComboBox(symopspanel);


	// get space groups and symops
	auto spacegroups = get_sgs<t_mat_bz>();
	m_sg_ops.reserve(spacegroups.size());
	for(auto [sgnum, descr, ops] : spacegroups)
	{
		m_comboSG->addItem(descr.c_str(), m_comboSG->count());
		m_sg_ops.emplace_back(std::move(ops));
	}


	auto tabGrid = new QGridLayout(symopspanel);
	tabGrid->setSpacing(2);
	tabGrid->setContentsMargins(4,4,4,4);
	int y = 0;
	//tabGrid->addWidget(m_plot.get(), y,0,1,4);
	tabGrid->addWidget(m_symops, y,0,1,4);
	tabGrid->addWidget(btnAdd, ++y,0,1,1);
	tabGrid->addWidget(btnDel, y,1,1,1);
	tabGrid->addWidget(btnUp, y,2,1,1);
	tabGrid->addWidget(btnDown, y,3,1,1);

	auto sep1 = new QFrame(symopspanel);
	sep1->setFrameStyle(QFrame::HLine);
	tabGrid->addWidget(sep1, ++y,0, 1,4);

	tabGrid->addWidget(new QLabel("Get Symmetry Operations From Space Group:"), ++y,0,1,4);
	tabGrid->addWidget(m_comboSG, ++y,0,1,3);
	tabGrid->addWidget(btnSG, y,3,1,1);

	auto sep2 = new QFrame(symopspanel);
	sep2->setFrameStyle(QFrame::HLine);
	tabGrid->addWidget(sep2, ++y,0, 1,4);

	tabGrid->addWidget(new QLabel("Lattice (\xe2\x84\xab):"), ++y,0,1,1);
	tabGrid->addWidget(m_editA, y,1,1,1);
	tabGrid->addWidget(m_editB, y,2,1,1);
	tabGrid->addWidget(m_editC, y,3,1,1);
	tabGrid->addWidget(new QLabel("Angles (\xc2\xb0):"), ++y,0,1,1);
	tabGrid->addWidget(m_editAlpha, y,1,1,1);
	tabGrid->addWidget(m_editBeta, y,2,1,1);
	tabGrid->addWidget(m_editGamma, y,3,1,1);


	// context menu
	m_symOpContextMenu = new QMenu(m_symops);
	m_symOpContextMenu->addAction("Add Symmetry Operation Before", this,
		[this]() { this->AddSymOpTabItem(-2); });
	m_symOpContextMenu->addAction("Add Symmetry Operation After", this,
		[this]() { this->AddSymOpTabItem(-3); });
	m_symOpContextMenu->addAction("Clone Symmetry Operation", this,
		[this]() { this->AddSymOpTabItem(-4); });
	m_symOpContextMenu->addAction("Delete Symmetry Operation", this,
		[this]() { BZDlg::DelSymOpTabItem(); });

	// context menu in case nothing is selected
	m_symOpContextMenuNoItem = new QMenu(m_symops);
	m_symOpContextMenuNoItem->addAction("Add Symmetry Operation", this,
		[this]() { this->AddSymOpTabItem(); });
	m_symOpContextMenuNoItem->addAction("Delete Symmetry Operation", this,
		[this]() { BZDlg::DelSymOpTabItem(); });


	// signals
	for(QDoubleSpinBox* spin : { m_editA, m_editB, m_editC,
		m_editAlpha, m_editBeta, m_editGamma })
	{
		connect(spin, static_cast<void (QDoubleSpinBox::*)(double)>(
			&QDoubleSpinBox::valueChanged),
			[this]() { this->CalcB(); });
	}

	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddSymOpTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { BZDlg::DelSymOpTabItem(); });
	connect(btnUp, &QAbstractButton::clicked,
		this, &BZDlg::MoveSymOpTabItemUp);
	connect(btnDown, &QAbstractButton::clicked,
		this, &BZDlg::MoveSymOpTabItemDown);
	connect(btnSG, &QAbstractButton::clicked,
		this, &BZDlg::GetSymOpsFromSG);

	connect(m_symops, &QTableWidget::itemChanged,
		this, &BZDlg::SymOpTableItemChanged);
	connect(m_symops, &QTableWidget::customContextMenuRequested,
		this, &BZDlg::ShowSymOpTableContextMenu);

	m_tabs_in->addTab(symopspanel, "Crystal");
}


/**
 * brillouin zone and cuts panel
 */
void BZDlg::CreateBZPanel()
{
	// already created?
	if(m_bzview)
		return;

	auto bzpanel = new QWidget(this);
	auto grid = new QGridLayout(bzpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(4,4,4,4);

	m_bzscene = new BZCutScene<t_vec_bz, t_real>(bzpanel);
	m_bzview = new BZCutView<t_vec_bz, t_real>(m_bzscene, m_sett);

	for(QDoubleSpinBox** const cut : {
		&m_cutNX, &m_cutNY, &m_cutNZ,
		&m_cutX, &m_cutY, &m_cutZ,
		&m_cutD })
	{
		*cut = new QDoubleSpinBox(bzpanel);
		(*cut)->setMinimum(-99);
		(*cut)->setMaximum(99);
		(*cut)->setDecimals(4);
		(*cut)->setSingleStep(0.1);
		(*cut)->setValue(0);

		// signals
		connect(*cut, static_cast<void (QDoubleSpinBox::*)(double)>(
			&QDoubleSpinBox::valueChanged),
			[this]() { this->CalcBZCut(); });
	}

	m_cutX->setValue(1);
	m_cutNZ->setValue(1);

	m_cutNX->setToolTip("First component of the cutting plane's normal vector, which defines the plane's orientation.");
	m_cutNY->setToolTip("Second component of the cutting plane's normal vector, which defines the plane's orientation..");
	m_cutNZ->setToolTip("Third component of the cutting plane's normal vector, which defines the plane's orientation..");
	m_cutX->setToolTip("First component of the cutting plane's first in-plane vector, which defines the cut's rotation.");
	m_cutY->setToolTip("Second component of the cutting plane's first in-plane vector, which defines the cut's rotation.");
	m_cutZ->setToolTip("Third component of the cutting plane's first in-plane vector, which defines the cut's rotation.");
	m_cutD->setToolTip("Distance of the cutting plane along its normal vector from the origin.");

	int draw_order = 4;
	m_BZDrawOrder = new QSpinBox(bzpanel);
	m_BZDrawOrder->setMinimum(0);
	m_BZDrawOrder->setMaximum(99);
	m_BZDrawOrder->setValue(draw_order);
	m_BZDrawOrder->setToolTip("The maximum order of Brillouin zones to draw.");

	int calc_order = 4;
	m_BZCalcOrder = new QSpinBox(bzpanel);
	m_BZCalcOrder->setMinimum(1);
	m_BZCalcOrder->setMaximum(99);
	m_BZCalcOrder->setValue(calc_order);
	m_BZCalcOrder->setToolTip("The order of Bragg peaks used to calculate the Brillouin zone.");

	QPushButton *btnShowBZ = new QPushButton("3D View...", bzpanel);

	// cuts
	grid->addWidget(m_bzview, 0, 0, 1, 4);
	grid->addWidget(new QLabel("Plane Normal (rlu):"), 1, 0, 1, 1);
	grid->addWidget(m_cutNX, 1, 1, 1, 1);
	grid->addWidget(m_cutNY, 1, 2, 1, 1);
	grid->addWidget(m_cutNZ, 1, 3, 1, 1);
	grid->addWidget(new QLabel("In-Plane Vector (rlu):"), 2, 0, 1, 1);
	grid->addWidget(m_cutX, 2, 1, 1, 1);
	grid->addWidget(m_cutY, 2, 2, 1, 1);
	grid->addWidget(m_cutZ, 2, 3, 1, 1);
	grid->addWidget(new QLabel("Plane Offset (rlu):"), 3, 0, 1, 1);
	grid->addWidget(m_cutD, 3, 1, 1, 1);
	grid->addWidget(new QLabel("Draw Order:"), 3, 2, 1, 1);
	grid->addWidget(m_BZDrawOrder, 3, 3, 1, 1);

	// bz
	auto sep1 = new QFrame(bzpanel); sep1->setFrameStyle(QFrame::HLine);
	grid->addWidget(sep1, 4,0,1,4);
	grid->addWidget(new QLabel("Calc. Order:"), 5,0,1,1);
	grid->addWidget(m_BZCalcOrder, 5,1, 1,1);
	grid->addWidget(btnShowBZ, 5,3, 1,1);

	// signals
	connect(m_BZDrawOrder,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this](int order) { this->SetDrawOrder(order); });

#ifdef BZ_USE_QT_SIGNALS
	connect(m_bzview, &BZCutView<t_vec_bz, t_real>::SignalMouseCoordinates,
		this, &BZDlg::BZCutMouseMoved);
#else
	m_bzview->AddMouseCoordinatesSlot([this](t_real x, t_real y)
	{
		this->BZCutMouseMoved(x, y);
	});
#endif

	connect(m_BZCalcOrder,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this](int order) { this->SetCalcOrder(order); });
	connect(btnShowBZ, &QPushButton::clicked, this, &BZDlg::ShowBZPlot);

	SetCalcOrder(calc_order, false);
	SetDrawOrder(draw_order, false);

	m_tabs_out->addTab(bzpanel, "Brillouin Zone");
}


/**
 * panel for entering formulas to plot
 */
void BZDlg::CreateFormulasPanel()
{
	// already created?
	if(m_formulas)
		return;

	QWidget *formulaspanel = new QWidget(this);

	m_formulas = new QTableWidget(formulaspanel);
	m_formulas->setShowGrid(true);
	m_formulas->setSortingEnabled(true);
	m_formulas->setMouseTracking(true);
	m_formulas->setSelectionBehavior(QTableWidget::SelectRows);
	m_formulas->setSelectionMode(QTableWidget::ContiguousSelection);
	m_formulas->setContextMenuPolicy(Qt::CustomContextMenu);
	m_formulas->verticalHeader()->setVisible(false);
	m_formulas->setAlternatingRowColors(true);
	m_formulas->setColumnCount(NUM_FORMULAS_COLS);
	m_formulas->setHorizontalHeaderItem(COL_FORMULA, new QTableWidgetItem{"Formula to Plot (Variables are Qx and Qy)"});
	m_formulas->setColumnWidth(COL_FORMULA, 325);

	QToolButton *btnAdd = new QToolButton(formulaspanel);
	QToolButton *btnDel = new QToolButton(formulaspanel);
	QToolButton *btnUp = new QToolButton(formulaspanel);
	QToolButton *btnDown = new QToolButton(formulaspanel);

	m_formulas->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	btnAdd->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnDel->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnUp->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	btnDown->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	btnAdd->setText("Add");
	btnDel->setText("Delete");
	btnUp->setText("Up");
	btnDown->setText("Down");

	btnAdd->setToolTip("Add a new formula.");
	btnDel->setToolTip("Delete the selected formula.");
	btnUp->setToolTip("Move the selected formula up.");
	btnDown->setToolTip("Move the selected formula down.");

	btnAdd->setIcon(QIcon::fromTheme("list-add"));
	btnDel->setIcon(QIcon::fromTheme("list-remove"));
	btnUp->setIcon(QIcon::fromTheme("go-up"));
	btnDown->setIcon(QIcon::fromTheme("go-down"));

	auto tabGrid = new QGridLayout(formulaspanel);
	tabGrid->setSpacing(2);
	tabGrid->setContentsMargins(4,4,4,4);
	int y = 0;
	tabGrid->addWidget(m_formulas, y,0,1,4);
	tabGrid->addWidget(btnAdd, ++y,0,1,1);
	tabGrid->addWidget(btnDel, y,1,1,1);
	tabGrid->addWidget(btnUp, y,2,1,1);
	tabGrid->addWidget(btnDown, y,3,1,1);


	// table CustomContextMenu
	m_formulasContextMenu = new QMenu(m_formulas);
	m_formulasContextMenu->addAction("Add Formula Before", this,
		[this]() { this->AddFormulaTabItem(-2); });
	m_formulasContextMenu->addAction("Add Formula After", this,
		[this]() { this->AddFormulaTabItem(-3); });
	m_formulasContextMenu->addAction("Clone Formula", this,
		[this]() { this->AddFormulaTabItem(-4); });
	m_formulasContextMenu->addAction("Delete Formula", this,
		[this]() { BZDlg::DelFormulaTabItem(); });

	// table CustomContextMenu in case nothing is selected
	m_formulasContextMenuNoItem = new QMenu(m_formulas);
	m_formulasContextMenuNoItem->addAction("Add Formula", this,
		[this]() { this->AddFormulaTabItem(); });
	m_formulasContextMenuNoItem->addAction("Delete Formula", this,
		[this]() { BZDlg::DelFormulaTabItem(); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddFormulaTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { BZDlg::DelFormulaTabItem(); });
	connect(btnUp, &QAbstractButton::clicked,
		this, &BZDlg::MoveFormulaTabItemUp);
	connect(btnDown, &QAbstractButton::clicked,
		this, &BZDlg::MoveFormulaTabItemDown);

	connect(m_formulas, &QTableWidget::itemChanged,
		this, &BZDlg::FormulaTableItemChanged);
	connect(m_formulas, &QTableWidget::customContextMenuRequested,
		this, &BZDlg::ShowFormulaTableContextMenu);

	m_tabs_in->addTab(formulaspanel, "Formulas");
}


/**
 * brillouin zone calculation results panel
 */
void BZDlg::CreateResultsPanel()
{
	// already created?
	if(m_bzresults)
		return;

	auto resultspanel = new QWidget(this);
	auto grid = new QGridLayout(resultspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(4,4,4,4);

	m_bzresults = new QPlainTextEdit(resultspanel);
	m_bzresults->setReadOnly(true);
	m_bzresults->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));

	grid->addWidget(m_bzresults, 0,0, 1,4);

	m_tabs_out->addTab(resultspanel, "Results");
}


/**
 * brillouin zone calculation results panel (json)
 */
void BZDlg::CreateJsonResultsPanel()
{
	// already created?
	if(m_bzresultsJSON)
		return;

	auto resultspanel = new QWidget(this);
	auto grid = new QGridLayout(resultspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(4,4,4,4);

	m_bzresultsJSON = new QPlainTextEdit(resultspanel);
	m_bzresultsJSON->setReadOnly(true);
	m_bzresultsJSON->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));

	grid->addWidget(m_bzresultsJSON, 0,0, 1,4);

	m_tabs_out->addTab(resultspanel, "Results (JSON)");
}


void BZDlg::CreateMenuBar(QGridLayout *main_grid)
{
	// already created?
	if(m_menu)
		return;

	m_menu = new QMenuBar(this);
	m_menu->setNativeMenuBar(
		m_sett ? m_sett->value("native_gui", false).toBool() : false);

	auto menuFile = new QMenu("File", m_menu);
	auto menuBZ = new QMenu("Brillouin Zone", m_menu);

	// recent files menu
	m_menuOpenRecent = new QMenu("Open Recent", menuFile);

	m_recent.SetRecentFilesMenu(m_menuOpenRecent);
	m_recent.SetMaxRecentFiles(16);
	m_recent.SetOpenFunc(&m_open_func);

	// file menu
	auto acNew = new QAction("New", menuFile);
	auto acLoad = new QAction("Open...", menuFile);
	auto acSave = new QAction("Save...", menuFile);
	auto acImportCIF = new QAction("Import CIF...", menuFile);
	auto acExit = new QAction("Quit", menuFile);

	// bz menu
	auto ac3DView = new QAction("3D View...", menuFile);
	auto acCutSVG = new QAction("Save Cut to SVG...", menuFile);
	m_acCutHull = new QAction("Calculate Convex Hull for Cut", menuFile);

	// help menu
	auto menuHelp = new QMenu("Help", m_menu);
	auto *acAboutQt = new QAction("About Qt...", menuHelp);
	auto *acAbout = new QAction("About...", menuHelp);

	acNew->setIcon(QIcon::fromTheme("document-new"));
	acLoad->setIcon(QIcon::fromTheme("document-open"));
	acSave->setIcon(QIcon::fromTheme("document-save"));
	acExit->setIcon(QIcon::fromTheme("application-exit"));
	acAboutQt->setIcon(QIcon::fromTheme("help-about"));
	acAbout->setIcon(QIcon::fromTheme("help-about"));
	m_menuOpenRecent->setIcon(QIcon::fromTheme("document-open-recent"));

	acNew->setShortcut(QKeySequence::New);
	acLoad->setShortcut(QKeySequence::Open);
	acSave->setShortcut(QKeySequence::Save);
	acExit->setShortcut(QKeySequence::Quit);

	acExit->setMenuRole(QAction::QuitRole);
	acAboutQt->setMenuRole(QAction::AboutQtRole);
	acAbout->setMenuRole(QAction::AboutRole);

	m_acCutHull->setCheckable(true);
	m_acCutHull->setChecked(true);

	menuFile->addAction(acNew);
	menuFile->addSeparator();
	menuFile->addAction(acLoad);
	menuFile->addMenu(m_menuOpenRecent);
	menuFile->addSeparator();
	menuFile->addAction(acSave);
	menuFile->addSeparator();
	menuFile->addAction(acImportCIF);
	menuFile->addSeparator();
	menuFile->addAction(acExit);

	menuBZ->addAction(m_acCutHull);
	menuBZ->addAction(acCutSVG);
	menuBZ->addSeparator();
	menuBZ->addAction(ac3DView);

	menuHelp->addAction(acAboutQt);
	menuHelp->addAction(acAbout);

	connect(acNew, &QAction::triggered, this, &BZDlg::NewFile);
	connect(acLoad, &QAction::triggered, this,
		static_cast<void (BZDlg::*)()>(&BZDlg::Load));
	connect(acSave, &QAction::triggered, this,
		static_cast<void (BZDlg::*)()>(&BZDlg::Save));
	connect(acImportCIF, &QAction::triggered, this, &BZDlg::ImportCIF);
	connect(acExit, &QAction::triggered, this, &QDialog::close);
	connect(ac3DView, &QAction::triggered, this, &BZDlg::ShowBZPlot);
	connect(acCutSVG, &QAction::triggered, this, &BZDlg::SaveCutSVG);
	connect(m_acCutHull, &QAction::triggered, [this]() { this->CalcBZCut(); });
	connect(acAboutQt, &QAction::triggered, []() { qApp->aboutQt(); });
	connect(acAbout, &QAction::triggered, [this]()
	{
		if(!m_dlgInfo)
			return;

		m_dlgInfo->show();
		m_dlgInfo->raise();
		m_dlgInfo->activateWindow();
	});

	m_menu->addMenu(menuFile);
	m_menu->addMenu(menuBZ);
	m_menu->addMenu(menuHelp);

	main_grid->setMenuBar(m_menu);
}


void BZDlg::CreateInfoDialog()
{
	// already created?
	if(m_dlgInfo)
		return;

	auto infopanel = new QWidget(this);
	auto grid = new QGridLayout(infopanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	auto sep1 = new QFrame(infopanel); sep1->setFrameStyle(QFrame::HLine);
	auto sep2 = new QFrame(infopanel); sep2->setFrameStyle(QFrame::HLine);
	auto sep3 = new QFrame(infopanel); sep3->setFrameStyle(QFrame::HLine);
	auto sep4 = new QFrame(infopanel); sep4->setFrameStyle(QFrame::HLine);

	std::string strBoost = BOOST_LIB_VERSION;
	algo::replace_all(strBoost, "_", ".");

	auto labelTitle = new QLabel("Takin / Brillouin Zone Calculator", infopanel);
	auto fontTitle = labelTitle->font();
	fontTitle.setBold(true);
	labelTitle->setFont(fontTitle);
	labelTitle->setAlignment(Qt::AlignHCenter);

	auto labelAuthor = new QLabel("Written by Tobias Weber <tweber@ill.fr>.", infopanel);
	labelAuthor->setAlignment(Qt::AlignHCenter);

	auto labelDate = new QLabel("May 2022.", infopanel);
	labelDate->setAlignment(Qt::AlignHCenter);

	auto labelPaper = new QLabel("Paper DOI: "
		"<a href=\"https://doi.org/10.1016/j.softx.2023.101471\">10.1016/j.softx.2023.101471</a>.",
		infopanel);
	labelPaper->setWordWrap(true);
	labelPaper->setOpenExternalLinks(true);

	// renderer infos
	for(int i = 0; i < 4; ++i)
	{
		m_labelGlInfos[i] = new QLabel("", infopanel);
		m_labelGlInfos[i]->setSizePolicy(
			QSizePolicy::Ignored,
			m_labelGlInfos[i]->sizePolicy().verticalPolicy());
	}

	int y = 0;
	grid->addWidget(labelTitle, y++,0, 1,1);
	grid->addWidget(labelAuthor, y++,0, 1,1);
	grid->addWidget(labelDate, y++,0, 1,1);

	grid->addItem(new QSpacerItem(16,16,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,1);

	grid->addWidget(labelPaper, y++,0, 1,1);
	grid->addWidget(sep2, y++,0, 1,1);

	grid->addWidget(new QLabel(
		QString("Compiler: ") +
		QString(BOOST_COMPILER) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("C++ Library: ") +
		QString(BOOST_STDLIB) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("Build Date: ") +
		QString(__DATE__) + ", " +
		QString(__TIME__) + ".",
		infopanel), y++,0, 1,1);

	grid->addWidget(sep3, y++,0, 1,1);

	auto labelQt = new QLabel(QString(
		"<a href=\"http://code.qt.io/cgit/\">Qt</a>"
		" Version: %1.").arg(QT_VERSION_STR),
		infopanel);
	labelQt->setOpenExternalLinks(true);
	grid->addWidget(labelQt, y++,0, 1,1);

	auto labelBoost = new QLabel(QString(
		"<a href=\"http://www.boost.org\">Boost</a>"
		" Version: %1.").arg(strBoost.c_str()),
		infopanel);
	labelBoost->setOpenExternalLinks(true);
	grid->addWidget(labelBoost, y++,0, 1,1);

	auto labelGemmi = new QLabel(QString(
		"<a href=\"https://github.com/project-gemmi/gemmi\">Gemmi</a>"
		" Version: %1.").arg(get_gemmi_version<QString>()),
		infopanel);
	labelGemmi->setOpenExternalLinks(true);
	grid->addWidget(labelGemmi, y++,0, 1,1);

	grid->addWidget(sep4, y++,0, 1,1);

	for(int i = 0; i < 4; ++i)
		grid->addWidget(m_labelGlInfos[i], y++,0, 1,1);

	grid->addItem(new QSpacerItem(16,16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0, 1,1);

	// add info panel as a tab
	//m_tabs->addTab(infopanel, "Infos");

	// show info panel as a dialog
	m_dlgInfo = new QDialog(this);
	m_dlgInfo->setWindowTitle("About");
	m_dlgInfo->setSizeGripEnabled(true);
	m_dlgInfo->setFont(this->font());

	QPushButton *infoDlgOk = new QPushButton("OK", m_dlgInfo);
	infoDlgOk->setIcon(style()->standardIcon(QStyle::SP_DialogOkButton));
	connect(infoDlgOk, &QAbstractButton::clicked, m_dlgInfo, &QDialog::accept);

	auto dlgGrid = new QGridLayout(m_dlgInfo);
	dlgGrid->setSpacing(8);
	dlgGrid->setContentsMargins(8, 8, 8, 8);
	dlgGrid->addWidget(infopanel, 0,0, 1,4);
	dlgGrid->addWidget(infoDlgOk, 1,3, 1,1);
}


/**
 * a file is being dragged over the window
 */
void BZDlg::dragEnterEvent(QDragEnterEvent *evt)
{
	if(evt)
		evt->accept();
}


/**
 * a file is being dropped onto the window
 */
void BZDlg::dropEvent(QDropEvent *evt)
{
	const QMimeData *mime = evt->mimeData();
	if(!mime)
		return;

	for(const QUrl& url : mime->urls())
	{
		if(!url.isLocalFile())
			continue;

		Load(url.toLocalFile());
		evt->accept();
		break;
	}
}


void BZDlg::closeEvent(QCloseEvent *)
{
	if(!m_sett)
		return;

	m_recent.TrimEntries();
	m_sett->setValue("recent_files", m_recent.GetRecentFiles());
	m_sett->setValue("geo", saveGeometry());

	if(m_split_inout)
		m_sett->setValue("splitter", m_split_inout->saveState());
}


void BZDlg::UpdateBZDescription()
{
	// brillouin zone description
	std::string descr = m_descrBZ + "\n" + m_descrBZCut;
	m_bzresults->setPlainText(descr.c_str());

	// brillouin zone json description
	descr = "{\n";
	descr = m_descrBZJSON;
	if(m_descrBZCutJSON.size())
		descr += ",\n\n" + m_descrBZCutJSON;
	descr += "\n}";
	m_bzresultsJSON->setPlainText(descr.c_str());
}


/**
 * show the 3d plotter
 */
void BZDlg::ShowBZPlot()
{
	if(!m_dlgPlot)
	{
		m_dlgPlot = new BZPlotDlg(this, m_sett);

		connect(m_dlgPlot, &BZPlotDlg::NeedRecalc,
			[this]() { this->CalcB(); } );
		connect(m_dlgPlot, &BZPlotDlg::GlDeviceInfos,
			[this](const std::string& ver, const std::string& shader_ver,
				const std::string& vendor, const std::string& renderer)
		{
			if(!m_labelGlInfos[0] || !m_labelGlInfos[1] ||
				!m_labelGlInfos[2] || !m_labelGlInfos[3])
				return;

			m_labelGlInfos[0]->setText(QString("GL Version: %1.").arg(ver.c_str()));
			m_labelGlInfos[1]->setText(QString("GL Shader Version: %1.").arg(shader_ver.c_str()));
			m_labelGlInfos[2]->setText(QString("GL Vendor: %1.").arg(vendor.c_str()));
			m_labelGlInfos[3]->setText(QString("GL Device: %1.").arg(renderer.c_str()));
		});
	}

	m_dlgPlot->show();
	m_dlgPlot->raise();
	m_dlgPlot->focusWidget();
}
