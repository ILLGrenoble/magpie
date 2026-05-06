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
 * panel for saved favourite Q positions and paths
 */
void MagDynDlg::CreateCoordinatesPanel()
{
	m_coordinatespanel = new QWidget(this);

	// table with saved fields
	m_coordinatestab = new QTableWidget(m_coordinatespanel);
	m_coordinatestab->setShowGrid(true);
	m_coordinatestab->setAlternatingRowColors(true);
	m_coordinatestab->setSortingEnabled(true);
	m_coordinatestab->setMouseTracking(true);
	m_coordinatestab->setSelectionBehavior(QTableWidget::SelectRows);
	m_coordinatestab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_coordinatestab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_coordinatestab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing()*1.25 + 4);
	m_coordinatestab->verticalHeader()->setVisible(true);

	m_coordinatestab->setColumnCount(NUM_COORD_COLS);
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_NAME, new QTableWidgetItem{"Name"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_H, new QTableWidgetItem{"h"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_K, new QTableWidgetItem{"k"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_L, new QTableWidgetItem{"l"});

	m_coordinatestab->setColumnWidth(COL_COORD_NAME, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_H, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_K, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_L, 90);
	m_coordinatestab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAddCoord = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_coordinatespanel);
	QPushButton *btnDelCoord = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_coordinatespanel);
	QPushButton *btnCoordUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_coordinatespanel);
	QPushButton *btnCoordDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_coordinatespanel);

	btnAddCoord->setToolTip("Add a Q coordinate.");
	btnDelCoord->setToolTip("Delete selected Q coordinate.");
	btnCoordUp->setToolTip("Move selected coordinate(s) up.");
	btnCoordDown->setToolTip("Move selected coordinate(s) down.");

	QPushButton *btnSaveMultiDisp = new QPushButton(
		QIcon::fromTheme("text-x-generic"),
		"Save Data...", m_coordinatespanel);
	QPushButton *btnSaveMultiDispScr = new QPushButton(
		QIcon::fromTheme("text-x-script"),
		"Save Script...", m_coordinatespanel);
	btnSaveMultiDisp->setToolTip("Calculate the dispersion paths and save them to a data file.");
	btnSaveMultiDispScr->setToolTip("Calculate the dispersion paths and save them to a script file.");

	QPushButton *btnSetDispersion = new QPushButton("To Dispersion", m_coordinatespanel);
	QPushButton *btnSetHamilton = new QPushButton("To Hamiltonian", m_coordinatespanel);
	btnSetDispersion->setToolTip("Calculate the dispersion relation for the currently selected Q path.");
	btnSetHamilton->setToolTip("Calculate the Hamiltonian for the currently selected Q coordinate.");

	for(QPushButton *btn : { btnAddCoord, btnDelCoord, btnCoordUp, btnCoordDown })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}


	// table context menu
	QMenu *menuTableContext = new QMenu(m_coordinatestab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coordinate Before", this,
		[this]()
	{
		this->AddCoordinateTabItem(-2, "",
			m_Q_start[0]->value(), m_Q_start[1]->value(), m_Q_start[2]->value());
	});
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coordinate After", this,
		[this]()
	{
		this->AddCoordinateTabItem(-3, "",
			m_Q_start[0]->value(), m_Q_start[1]->value(), m_Q_start[2]->value());
	});
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Coordinate", this,
		[this]() { this->AddCoordinateTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coordinate", this,
		[this]() { this->DelTabItem(m_coordinatestab); });
	menuTableContext->addSeparator();
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Calculate Dispersion From This To Next Q", this,
		[this]() { this->SetCurrentCoordinate(0); });
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Calculate Hamiltonian For This Q", this,
		[this]() { this->SetCurrentCoordinate(1); });


	// table context menu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_coordinatestab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Start Q Coordinate", this,
		[this]()
	{
		this->AddCoordinateTabItem(-1, "",
			m_Q_start[0]->value(), m_Q_start[1]->value(), m_Q_start[2]->value());
	});
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add End Q Coordinate", this,
		[this]()
	{
		this->AddCoordinateTabItem(-1, "",
			m_Q_end[0]->value(), m_Q_end[1]->value(), m_Q_end[2]->value());
	});
	menuTableContextNoItem->addSeparator();
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coordinate", this,
		[this]() { this->DelTabItem(m_coordinatestab); });


	QGridLayout *grid = new QGridLayout(m_coordinatespanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Saved Q Coordinates (rlu):", m_coordinatespanel), y++,0,1,4);
	grid->addWidget(m_coordinatestab, y,0,1,4);
	grid->addWidget(btnAddCoord, ++y,0,1,1);
	grid->addWidget(btnDelCoord, y,1,1,1);
	grid->addWidget(btnCoordUp, y,2,1,1);
	grid->addWidget(btnCoordDown, y++,3,1,1);
	grid->addWidget(btnSetDispersion, y,0,1,1);
	grid->addWidget(btnSetHamilton, y,1,1,1);
	grid->addWidget(btnSaveMultiDisp, y,2,1,1);
	grid->addWidget(btnSaveMultiDispScr, y++,3,1,1);

	// signals
	connect(btnAddCoord, &QAbstractButton::clicked,
		[this]()
	{
		this->AddCoordinateTabItem(-1, "",
			m_Q_start[0]->value(), m_Q_start[1]->value(), m_Q_start[2]->value());
	});
	connect(btnDelCoord, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_coordinatestab); });
	connect(btnCoordUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_coordinatestab); });
	connect(btnCoordDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_coordinatestab); });

	connect(btnSetDispersion, &QAbstractButton::clicked,
		[this]() { this->SetCurrentCoordinate(0); });
	connect(btnSetHamilton, &QAbstractButton::clicked,
		[this]() { this->SetCurrentCoordinate(1); });
	connect(btnSaveMultiDisp, &QAbstractButton::clicked,
		[this]() { this->SaveMultiDispersion(false); });
	connect(btnSaveMultiDispScr, &QAbstractButton::clicked,
		[this]() { this->SaveMultiDispersion(true); });

	connect(m_coordinatestab, &QTableWidget::itemSelectionChanged, this,
		&MagDynDlg::CoordinatesSelectionChanged);
	connect(m_coordinatestab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_coordinatestab, menuTableContext,
			menuTableContextNoItem, pt);
	});


	m_tabs_recip->addTab(m_coordinatespanel, "Coordinates");
}



/**
 * a coordinate has been selected
 */
void MagDynDlg::CoordinatesSelectionChanged()
{
	QList<QTableWidgetItem*> selected = m_coordinatestab->selectedItems();
	if(selected.size() == 0)
		return;

	const QTableWidgetItem* item = *selected.begin();
	m_coordinates_cursor_row = item->row();
}
