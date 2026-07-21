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
 * allows the user to specify magnetic sites
 */
void MagDynDlg::CreateSitesPanel()
{
	m_sitespanel = new QWidget(this);

	m_sitestab = new QTableWidget(m_sitespanel);
	m_sitestab->setShowGrid(true);
	m_sitestab->setAlternatingRowColors(true);
	m_sitestab->setSortingEnabled(true);
	m_sitestab->setMouseTracking(true);
	m_sitestab->setSelectionBehavior(QTableWidget::SelectRows);
	m_sitestab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_sitestab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_sitestab->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing()*1.25 + 4);
	m_sitestab->verticalHeader()->setVisible(true);

	m_sitestab->setColumnCount(NUM_SITE_COLS);

	m_sitestab->setHorizontalHeaderItem(COL_SITE_NAME, new QTableWidgetItem{"Name"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_X, new QTableWidgetItem{"x"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_Y, new QTableWidgetItem{"y"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_Z, new QTableWidgetItem{"z"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SYM_IDX, new QTableWidgetItem{"Sym. Idx."});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_X, new QTableWidgetItem{"Spin x"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_Y, new QTableWidgetItem{"Spin y"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_Z, new QTableWidgetItem{"Spin z"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_MAG, new QTableWidgetItem{"Spin |S|"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_FORMFACT_IDX, new QTableWidgetItem{"f_M Idx."});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_RGB, new QTableWidgetItem{"Colour"});

	if(m_allow_ortho_spin)
	{
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_X, new QTableWidgetItem{"Spin ux"});
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_Y, new QTableWidgetItem{"Spin uy"});
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_Z, new QTableWidgetItem{"Spin uz"});
	}
	else
	{
		m_sitestab->setColumnCount(NUM_SITE_COLS - 3);
	}

	m_sitestab->setColumnWidth(COL_SITE_NAME, 90);
	m_sitestab->setColumnWidth(COL_SITE_POS_X, 80);
	m_sitestab->setColumnWidth(COL_SITE_POS_Y, 80);
	m_sitestab->setColumnWidth(COL_SITE_POS_Z, 80);
	m_sitestab->setColumnWidth(COL_SITE_SYM_IDX, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_X, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_Y, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_Z, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_MAG, 80);
	m_sitestab->setColumnWidth(COL_SITE_FORMFACT_IDX, 80);
	m_sitestab->setColumnWidth(COL_SITE_RGB, 80);

	if(m_allow_ortho_spin)
	{
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_X, 80);
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_Y, 80);
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_Z, 80);
	}

	m_sitestab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_sitespanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_sitespanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_sitespanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_sitespanel);

	btnAdd->setToolTip("Add a site.");
	btnDel->setToolTip("Delete selected site(s).");
	btnUp->setToolTip("Move selected site(s) up.");
	btnDown->setToolTip("Move selected site(s) down.");

	QPushButton *btnMirrorAtoms = new QPushButton("Mirror", m_sitespanel);
	QPushButton *btnShowNotes = new QPushButton(
		QIcon::fromTheme("accessories-text-editor"),
		"Notes...", m_sitespanel);
	QPushButton *btnGroundState = new QPushButton("Ground State...", m_sitespanel);
	btnMirrorAtoms->setToolTip("Flip the coordinates of the sites.");
	btnShowNotes->setToolTip("Add notes or comments describing the magnetic structure.");
	btnGroundState->setToolTip("Minimise ground state energy.");
#ifndef __TLIBS2_MAGDYN_USE_MINUIT__
	btnGroundState->setEnabled(false);
#endif

	// extend cell
	const char* idx_names[] = {"x = ", "y = ", "z = "};
	for(int cell_idx = 0; cell_idx < 3; ++cell_idx)
	{
		m_extCell[cell_idx] = new QSpinBox(m_sitespanel);
		m_extCell[cell_idx]->setMinimum(1);
		m_extCell[cell_idx]->setMaximum(99);
		m_extCell[cell_idx]->setValue(cell_idx == 2 ? 2 : 1);
		m_extCell[cell_idx]->setPrefix(idx_names[cell_idx]);
		m_extCell[cell_idx]->setToolTip("Order of supercell.");
		m_extCell[cell_idx]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	QPushButton *btnExtCell = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_sitespanel);
	btnExtCell->setToolTip("Extend the unit cell.");

	QPushButton *btnGenBySG = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_sitespanel);
	btnGenBySG->setToolTip("Create site positions from space group"
		" symmetry operators and existing positions.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown, btnGenBySG, btnExtCell,
		btnGroundState, btnMirrorAtoms, btnShowNotes })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}


	QGridLayout *grid = new QGridLayout(m_sitespanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_sitestab, y, 0, 1, 4);
	grid->addWidget(btnAdd, ++y, 0, 1, 1);
	grid->addWidget(btnDel, y, 1, 1, 1);
	grid->addWidget(btnUp, y, 2, 1, 1);
	grid->addWidget(btnDown, y++, 3, 1, 1);
	grid->addWidget(btnMirrorAtoms, y, 0, 1, 1);
	grid->addWidget(btnShowNotes, y, 2, 1, 1);
	grid->addWidget(btnGroundState, y++, 3, 1, 1);

	QFrame *sep1 = new QFrame(m_sitespanel);
	sep1->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);

	grid->addWidget(new QLabel("Extend Structure, Copying Existing Sites:",
		m_sitespanel), y++, 0, 1, 4);
	grid->addWidget(m_extCell[0], y, 0, 1, 1);
	grid->addWidget(m_extCell[1], y, 1, 1, 1);
	grid->addWidget(m_extCell[2], y, 2, 1, 1);
	grid->addWidget(btnExtCell, y++, 3, 1, 1);

	grid->addWidget(new QLabel("Create Symmetry-Equivalent Sites:",
		m_sitespanel), y, 0, 1, 3);
	grid->addWidget(btnGenBySG, y++, 3, 1, 1);

	// table context menu
	QMenu *menuTableContext = new QMenu(m_sitestab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site Before Current", this,
		[this]() { this->AddSiteTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site After Current", this,
		[this]() { this->AddSiteTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Site", this,
		[this]() { this->AddSiteTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Site",this,
		[this]() { this->DelTabItem(m_sitestab); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Symmetry-Equivalent Sites", this,
		[this]() { this->DelIdentTabItems(m_sitestab, COL_SITE_SYM_IDX); });



	// table context menu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_sitestab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site", this,
		[this]() { this->AddSiteTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Sites", this,
		[this]() { this->DelTabItem(m_sitestab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddSiteTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_sitestab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_sitestab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_sitestab); });
	connect(btnGenBySG, &QAbstractButton::clicked,
		this, &MagDynDlg::GenerateSitesFromSG);
	connect(btnExtCell, &QAbstractButton::clicked,
		this, &MagDynDlg::ExtendStructure);

	connect(btnMirrorAtoms, &QAbstractButton::clicked, this, &MagDynDlg::MirrorAtoms);
	connect(btnShowNotes, &QAbstractButton::clicked, this, &MagDynDlg::ShowNotesDlg);
	connect(btnGroundState, &QAbstractButton::clicked, this, &MagDynDlg::ShowGroundStateDlg);

	connect(m_sitestab, &QTableWidget::itemSelectionChanged, this, &MagDynDlg::SitesSelectionChanged);
	connect(m_sitestab, &QTableWidget::itemChanged, this, &MagDynDlg::SitesTableItemChanged);
	connect(m_sitestab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_sitestab, menuTableContext, menuTableContextNoItem, pt);
	});


	m_tabs_setup->addTab(m_sitespanel, "Sites");
}



/**
 * a site has been selected
 */
void MagDynDlg::SitesSelectionChanged()
{
	QList<QTableWidgetItem*> selected = m_sitestab->selectedItems();
	if(selected.size() == 0)
		return;

	const QTableWidgetItem* item = *selected.begin();
	m_sites_cursor_row = item->row();
	if(m_sites_cursor_row < 0 ||
		std::size_t(m_sites_cursor_row) >= m_dyn.GetMagneticSitesCount())
	{
		m_status->setText("");
		return;
	}

	const t_site* site = GetSiteFromTableIndex(m_sites_cursor_row);
	if(!site)
	{
		m_status->setText("Invalid site selected.");
		return;
	}

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr << "Site " << site->name << ".";
	m_status->setText(ostr.str().c_str());
}
