/**
 * choose magnon bands
 * @author Tobias Weber <tweber@ill.fr>
 * @date 23-august-2026
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

#include "bands.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>



/**
 * column indices in magnon band table for the group velocity
 */
enum : int
{
	COL_BAND = 0,
	COL_ACTIVE,
	NUM_COLS,
};



/**
 * set up the gui
 */
BandsDlg::BandsDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	m_emit = false;
	setWindowTitle("Magnon Bands");
	setSizeGripEnabled(true);


	QPushButton *btnReset = new QPushButton("Reset", this);

	m_table_bands = new QTableWidget(this);
	m_table_bands->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_bands->setShowGrid(true);
	m_table_bands->setSortingEnabled(false);
	m_table_bands->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_bands->setSelectionMode(QTableWidget::SingleSelection);
	m_table_bands->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_bands->verticalHeader()->setVisible(false);
	m_table_bands->setColumnCount(NUM_COLS);
	m_table_bands->setHorizontalHeaderItem(COL_BAND, new QTableWidgetItem{"Band"});
	m_table_bands->setHorizontalHeaderItem(COL_ACTIVE, new QTableWidgetItem{"Active"});
	m_table_bands->setColumnWidth(COL_BAND, 45);
	m_table_bands->setColumnWidth(COL_ACTIVE, 45);
	m_table_bands->resizeColumnsToContents();

	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &BandsDlg::accept);


	// grid
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	grid->addWidget(m_table_bands, y++, 0, 1, 3);
	grid->addWidget(btnReset, y, 0, 1, 1);
	grid->addWidget(btnbox, y++, 2, 1, 1);


	// connections
	connect(btnReset, &QAbstractButton::clicked, this, &BandsDlg::Reset);


	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("bands/geo"))
			restoreGeometry(m_sett->value("bands/geo").toByteArray());
		else
			resize(400, 800);
  }


  Reset();
	m_emit = true;
}



void BandsDlg::EmitStateChanged()
{
	if(m_emit)
		emit StateChanged();
}



/**
 * clears all bands
 */
void BandsDlg::Clear(bool clear_settings)
{
	m_emit = false;

	// keep the previous band visibility flags
	if(!clear_settings)
	{
		m_enabled_bands.resize(m_table_bands->rowCount());
		for(int row = 0; row < m_table_bands->rowCount(); ++row)
			m_enabled_bands[row] = IsChecked(std::size_t(row));
	}
	else
	{
		m_enabled_bands.clear();
	}

	m_table_bands->clearContents();
	m_table_bands->setRowCount(0);

	m_emit = true;
}



/**
 * sets the default configuration
 */
void BandsDlg::Reset()
{
	m_emit = false;

	for(int row = 0; row < m_table_bands->rowCount(); ++row)
		SetChecked(std::size_t(row));

	m_emit = true;
	EmitStateChanged();
}



/**
 * adds a magnon band to the table
 */
void BandsDlg::AddBand(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_bands)
		return;

	int row = m_table_bands->rowCount();
	m_table_bands->insertRow(row);

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

	QCheckBox *checkBand = new QCheckBox(m_table_bands);
	checkBand->setChecked(enabled);
	connect(checkBand, &QCheckBox::toggled, [this]() { EmitStateChanged(); });

	m_table_bands->setItem(row, COL_BAND, item);
	m_table_bands->setCellWidget(row, COL_ACTIVE, checkBand);
}



/**
 * is the given band active?
 */
bool BandsDlg::IsChecked(std::size_t idx) const
{
	if(!m_table_bands || int(idx) >= m_table_bands->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(
		m_table_bands->cellWidget(int(idx), COL_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * sets the given band active
 */
void BandsDlg::SetChecked(std::size_t idx, bool active)
{
	if(!m_table_bands || int(idx) >= m_table_bands->rowCount())
		return;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(
		m_table_bands->cellWidget(int(idx), COL_ACTIVE));
	if(!box)
		return;

	box->setChecked(active);
}



/**
 * was the band previously active?
 */
bool BandsDlg::GetOldChecked(std::size_t idx) const
{
	bool enabled = idx < m_enabled_bands.size()
		? m_enabled_bands[idx]  // previous setting
		: true;                 // default setting

	return enabled;
}



/**
 * close the dialog
 */
void BandsDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("bands/geo", saveGeometry());
	}

	QDialog::accept();
}
