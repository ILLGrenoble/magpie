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
 * lets the user define variables to be used for the J and DMI parameters
 */
void MagDynDlg::CreateVariablesPanel()
{
	m_varspanel = new QWidget(this);

	m_varstab = new QTableWidget(m_varspanel);
	m_varstab->setShowGrid(true);
	m_varstab->setAlternatingRowColors(true);
	m_varstab->setSortingEnabled(true);
	m_varstab->setMouseTracking(true);
	m_varstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_varstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_varstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_varstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing()*1.25 + 4);
	m_varstab->verticalHeader()->setVisible(true);

	m_varstab->setColumnCount(NUM_VARS_COLS);
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_NAME, new QTableWidgetItem{"Name"});
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_VALUE_REAL, new QTableWidgetItem{"Value (Re)"});
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_VALUE_IMAG, new QTableWidgetItem{"Value (Im)"});

	m_varstab->setColumnWidth(COL_VARS_NAME, 150);
	m_varstab->setColumnWidth(COL_VARS_VALUE_REAL, 150);
	m_varstab->setColumnWidth(COL_VARS_VALUE_IMAG, 150);
	m_varstab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_varspanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_varspanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_varspanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_varspanel);
	QPushButton *btnReplace = new QPushButton(
		QIcon::fromTheme("edit-find-replace"),
		"Replace Values", m_varspanel);

	btnAdd->setToolTip("Add a variable.");
	btnDel->setToolTip("Delete selected variables(s).");
	btnUp->setToolTip("Move selected variable(s) up.");
	btnDown->setToolTip("Move selected variable(s) down.");
	btnReplace->setToolTip("Replace numeric values with variable names.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown, btnReplace })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}


	// grid
	QGridLayout *grid = new QGridLayout(m_varspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_varstab, y++,0,1,4);
	grid->addWidget(btnAdd, y,0,1,1);
	grid->addWidget(btnDel, y,1,1,1);
	grid->addWidget(btnUp, y,2,1,1);
	grid->addWidget(btnDown, y++,3,1,1);
	grid->addWidget(btnReplace, y++,0,1,1);


	// table context menu
	QMenu *menuTableContext = new QMenu(m_varstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable Before Current", this,
		[this]() { this->AddVariableTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable After Current", this,
		[this]() { this->AddVariableTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Variable", this,
		[this]() { this->AddVariableTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Variable(s)", this,
		[this]() { this->DelTabItem(m_varstab); });


	// table context menu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_varstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable", this,
		[this]() { this->AddVariableTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Variable(s)", this,
		[this]() { this->DelTabItem(m_varstab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddVariableTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_varstab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_varstab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_varstab); });
	connect(btnReplace, &QAbstractButton::clicked,
		this, &MagDynDlg::ReplaceValuesWithVariables);

	connect(m_varstab, &QTableWidget::itemSelectionChanged, this, &MagDynDlg::VariablesSelectionChanged);
	connect(m_varstab, &QTableWidget::itemChanged, this, &MagDynDlg::VariablesTableItemChanged);
	connect(m_varstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_varstab, menuTableContext, menuTableContextNoItem, pt);
	});


	m_tabs_setup->addTab(m_varspanel, "Variables");
}



/**
 * a variable has been selected
 */
void MagDynDlg::VariablesSelectionChanged()
{
	QList<QTableWidgetItem*> selected = m_varstab->selectedItems();
	if(selected.size() == 0)
		return;

	const QTableWidgetItem* item = *selected.begin();
	m_variables_cursor_row = item->row();
}
