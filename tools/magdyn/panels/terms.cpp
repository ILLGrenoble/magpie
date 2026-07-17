/**
 * magnetic dynamics -- gui setup
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#include "magdyn.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>



/**
 * allows the user to specify magnetic couplings between sites
 */
void MagDynDlg::CreateExchangeTermsPanel()
{
	m_termspanel = new QWidget(this);

	m_termstab = new QTableWidget(m_termspanel);
	m_termstab->setShowGrid(true);
	m_termstab->setAlternatingRowColors(true);
	m_termstab->setSortingEnabled(true);
	m_termstab->setMouseTracking(true);
	m_termstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_termstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_termstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_termstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing()*1.25 + 4);
	m_termstab->verticalHeader()->setVisible(true);

	m_termstab->setColumnCount(NUM_XCH_COLS);
	m_termstab->setHorizontalHeaderItem(COL_XCH_NAME, new QTableWidgetItem{"Name"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_ATOM1_IDX, new QTableWidgetItem{"Site 1"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_ATOM2_IDX, new QTableWidgetItem{"Site 2"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DIST_X, new QTableWidgetItem{"Cell Δx"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DIST_Y, new QTableWidgetItem{"Cell Δy"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DIST_Z, new QTableWidgetItem{"Cell Δz"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_SYM_IDX, new QTableWidgetItem{"Sym. Idx."});
	m_termstab->setHorizontalHeaderItem(COL_XCH_INTERACTION, new QTableWidgetItem{"Exch. J"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DMI_X, new QTableWidgetItem{"DMI x"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DMI_Y, new QTableWidgetItem{"DMI y"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_DMI_Z, new QTableWidgetItem{"DMI z"});
	m_termstab->setHorizontalHeaderItem(COL_XCH_RGB, new QTableWidgetItem{"Colour"});

	if(m_allow_general_J)
	{
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_XX, new QTableWidgetItem{"J xx"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_XY, new QTableWidgetItem{"J xy"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_XZ, new QTableWidgetItem{"J xz"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_YX, new QTableWidgetItem{"J yx"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_YY, new QTableWidgetItem{"J yy"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_YZ, new QTableWidgetItem{"J yz"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_ZX, new QTableWidgetItem{"J zx"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_ZY, new QTableWidgetItem{"J zy"});
		m_termstab->setHorizontalHeaderItem(COL_XCH_GEN_ZZ, new QTableWidgetItem{"J zz"});
	}
	else
	{
		m_termstab->setColumnCount(NUM_XCH_COLS - 9);
	}

	m_termstab->setColumnWidth(COL_XCH_NAME, 90);
	m_termstab->setColumnWidth(COL_XCH_ATOM1_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_ATOM2_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_X, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_Y, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_Z, 80);
	m_termstab->setColumnWidth(COL_XCH_SYM_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_INTERACTION, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_X, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_Y, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_Z, 80);
	m_termstab->setColumnWidth(COL_XCH_RGB, 80);

	if(m_allow_general_J)
	{
		m_termstab->setColumnWidth(COL_XCH_GEN_XX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_XY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_XZ, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YZ, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZZ, 80);
	}

	m_termstab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_termspanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_termspanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_termspanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_termspanel);

	btnAdd->setToolTip("Add a coupling between two sites.");
	btnDel->setToolTip("Delete selected coupling(s).");
	btnUp->setToolTip("Move selected coupling(s) up.");
	btnDown->setToolTip("Move selected coupling(s) down.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}


	// couplings from distances
	m_maxdist = new QDoubleSpinBox(m_termspanel);
	m_maxdist->setDecimals(3);
	m_maxdist->setMinimum(0.001);
	m_maxdist->setMaximum(99.999);
	m_maxdist->setSingleStep(0.1);
	m_maxdist->setValue(5);
	m_maxdist->setPrefix("d = ");
	m_maxdist->setToolTip("Maximum distance between sites.");
	m_maxdist->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	m_maxSC = new QSpinBox(m_termspanel);
	m_maxSC->setMinimum(1);
	m_maxSC->setMaximum(99);
	m_maxSC->setValue(4);
	m_maxSC->setPrefix("order = ");
	m_maxSC->setToolTip("Maximum order of supercell to consider.");
	m_maxSC->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	m_maxcouplings = new QSpinBox(m_termspanel);
	m_maxcouplings->setMinimum(-1);
	m_maxcouplings->setMaximum(999);
	m_maxcouplings->setValue(100);
	m_maxcouplings->setPrefix("n = ");
	m_maxcouplings->setToolTip("Maximum number of couplings to generate (-1: no limit).");
	m_maxcouplings->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	QPushButton *btnGenByDist = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_termspanel);
	btnGenByDist->setToolTip("Create possible couplings by distances between sites.");
	btnGenByDist->setFocusPolicy(Qt::StrongFocus);
	btnGenByDist->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});


	// couplings from space group
	QPushButton *btnGenBySG = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_termspanel);
	btnGenBySG->setToolTip("Create couplings from space group"
		" symmetry operators and existing couplings.");
	btnGenBySG->setFocusPolicy(Qt::StrongFocus);
	btnGenBySG->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	// ordering vector
	m_ordering[0] = new QDoubleSpinBox(m_termspanel);
	m_ordering[1] = new QDoubleSpinBox(m_termspanel);
	m_ordering[2] = new QDoubleSpinBox(m_termspanel);

	// normal axis
	m_normaxis[0] = new QDoubleSpinBox(m_termspanel);
	m_normaxis[1] = new QDoubleSpinBox(m_termspanel);
	m_normaxis[2] = new QDoubleSpinBox(m_termspanel);

	for(int i = 0; i < 3; ++i)
	{
		m_ordering[i]->setDecimals(4);
		m_ordering[i]->setMinimum(-9.9999);
		m_ordering[i]->setMaximum(+9.9999);
		m_ordering[i]->setSingleStep(0.01);
		m_ordering[i]->setValue(0.);
		m_ordering[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

		m_normaxis[i]->setDecimals(4);
		m_normaxis[i]->setMinimum(-9.9999);
		m_normaxis[i]->setMaximum(+9.9999);
		m_normaxis[i]->setSingleStep(0.01);
		m_normaxis[i]->setValue(i==0 ? 1. : 0.);
		m_normaxis[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	m_ordering[0]->setPrefix("Oh = ");
	m_ordering[1]->setPrefix("Ok = ");
	m_ordering[2]->setPrefix("Ol = ");

	m_normaxis[0]->setPrefix("Nh = ");
	m_normaxis[1]->setPrefix("Nk = ");
	m_normaxis[2]->setPrefix("Nl = ");


	// grid
	QGridLayout *grid = new QGridLayout(m_termspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_termstab, y++,0,1,4);
	grid->addWidget(btnAdd, y,0,1,1);
	grid->addWidget(btnDel, y,1,1,1);
	grid->addWidget(btnUp, y,2,1,1);
	grid->addWidget(btnDown, y++,3,1,1);

	QFrame *sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	QFrame *sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Generate Possible Coupling Terms"
		" By Distance (\xe2\x84\xab):", m_termspanel), y++,0,1,4);
	grid->addWidget(m_maxdist, y,0,1,1);
	grid->addWidget(m_maxSC, y,1,1,1);
	grid->addWidget(m_maxcouplings, y,2,1,1);
	grid->addWidget(btnGenByDist, y++,3,1,1);
	grid->addWidget(new QLabel("Create Symmetry-Equivalent Couplings:",
		m_termspanel), y,0,1,3);
	grid->addWidget(btnGenBySG, y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Ordering Vector:", m_termspanel), y,0,1,1);
	grid->addWidget(m_ordering[0], y,1,1,1);
	grid->addWidget(m_ordering[1], y,2,1,1);
	grid->addWidget(m_ordering[2], y++,3,1,1);
	grid->addWidget(new QLabel("Rotation Axis:", m_termspanel), y,0,1,1);
	grid->addWidget(m_normaxis[0], y,1,1,1);
	grid->addWidget(m_normaxis[1], y,2,1,1);
	grid->addWidget(m_normaxis[2], y++,3,1,1);

	// table context menu
	QMenu *menuTableContext = new QMenu(m_termstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coupling Before", this,
		[this]() { this->AddTermTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coupling After", this,
		[this]() { this->AddTermTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Coupling", this,
		[this]() { this->AddTermTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coupling", this,
		[this]() { this->DelTabItem(m_termstab); });


	// table context menu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_termstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coupling", this,
		[this]() { this->AddTermTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coupling", this,
		[this]() { this->DelTabItem(m_termstab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddTermTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_termstab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_termstab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_termstab); });
	connect(btnGenByDist, &QAbstractButton::clicked,
		this, &MagDynDlg::GeneratePossibleCouplings);
	connect(btnGenBySG, &QAbstractButton::clicked,
		this, &MagDynDlg::GenerateCouplingsFromSG);

	connect(m_termstab, &QTableWidget::itemSelectionChanged, this, &MagDynDlg::TermsSelectionChanged);
	connect(m_termstab, &QTableWidget::itemChanged, this, &MagDynDlg::TermsTableItemChanged);
	connect(m_termstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_termstab, menuTableContext, menuTableContextNoItem, pt);
	});

	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};

	for(int i = 0; i < 3; ++i)
	{
		connect(m_ordering[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);

		connect(m_normaxis[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}


	m_tabs_setup->addTab(m_termspanel, "Couplings");
}



/**
 * a term has been selected
 */
void MagDynDlg::TermsSelectionChanged()
{
	QList<QTableWidgetItem*> selected = m_termstab->selectedItems();
	if(selected.size() == 0)
		return;

	const QTableWidgetItem* item = *selected.begin();
	m_terms_cursor_row = item->row();
	if(m_terms_cursor_row < 0 ||
		std::size_t(m_terms_cursor_row) >= m_dyn.GetExchangeTermsCount())
	{
		m_status->setText("");
		return;
	}

	const t_term* term = GetTermFromTableIndex(m_terms_cursor_row);
	if(!term)
	{
		m_status->setText("Invalid coupling selected.");
		return;
	}

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr << "Coupling " << term->name
	<< ": length = " << term->length_calc << " \xe2\x84\xab.";
	m_status->setText(ostr.str().c_str());
}
