/**
 * choose matrix elements
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-july-2026
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

#include "matrixelems.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>



/**
 * set up the gui
 */
MatrixElemsDlg::MatrixElemsDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	m_emit = false;
	setWindowTitle("Matrix Elements");
	setSizeGripEnabled(true);


	// matrix elements
	const char* comps[] = { "x", "y", "z" };

	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
	{
		m_elems[i*3 + j] = new QCheckBox(this);
		m_elems[i*3 + j]->setText(QString("Re{S%1%2}").arg(comps[i]).arg(comps[j]));

		m_elems[i*3 + j + 3*3] = new QCheckBox(this);
		m_elems[i*3 + j + 3*3]->setText(QString("Im{S%1%2}").arg(comps[i]).arg(comps[j]));
	}


	QPushButton *btnReset = new QPushButton("Reset", this);


	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &MatrixElemsDlg::accept);


	// grid
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
	{
		grid->addWidget(m_elems[i*3 + j],       y + i,     j, 1, 1);
		grid->addWidget(m_elems[i*3 + j + 3*3], y + i + 4, j, 1, 1);
	}

	auto sep1 = new QFrame(this);
	sep1->setFrameStyle(QFrame::HLine);
	grid->addWidget(sep1, y + 3, 0, 1, 3);

	y += 7;

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++, 0, 1, 3);
	grid->addWidget(btnReset, y, 0, 1, 1);
	grid->addWidget(btnbox, y++, 2, 1, 1);


	// connections
	connect(btnReset, &QAbstractButton::clicked, this, &MatrixElemsDlg::Reset);

	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
	{
		connect(m_elems[i*3 + j],       &QCheckBox::toggled, this, &MatrixElemsDlg::EmitStateChanged);
		connect(m_elems[i*3 + j + 3*3], &QCheckBox::toggled, this, &MatrixElemsDlg::EmitStateChanged);
	}


	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("matrixelems/geo"))
			restoreGeometry(m_sett->value("matrixelems/geo").toByteArray());
		else
			resize(400, 400);
  }


  Reset();
	m_emit = true;
}



void MatrixElemsDlg::EmitStateChanged()
{
	if(m_emit)
		emit StateChanged();
}



/**
 * sets the default configuration (diagonal elements)
 */
void MatrixElemsDlg::Reset()
{
	m_emit = false;

	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
	{
		m_elems[i*3 + j]->setChecked(i == j);
		m_elems[i*3 + j + 3*3]->setChecked(false);
	}

	m_emit = true;
	EmitStateChanged();
}



/**
 * is the given channel active?
 */
bool MatrixElemsDlg::IsChecked(std::size_t i, std::size_t j, bool real_elem) const
{
	if(i >= 3)
		i = 3;
	if(j >= 3)
		j = 3;

	std::size_t idx = real_elem ? 0 : 3*3;  // start index
	return m_elems[i*3 + j + idx]->isChecked();
}



/**
 * close the dialog
 */
void MatrixElemsDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("matrixelems/geo", saveGeometry());
	}

	QDialog::accept();
}
