/**
 * magnon dynamics -- assigns multiple couplings
 * @author Tobias Weber <tweber@ill.fr>
 * @date 6-june-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core / magdyn / magpie
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

#include "assign.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>



/**
 * set up the gui
 */
AssignDlg::AssignDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	setWindowTitle("Assign Multiple Couplings");
	setSizeGripEnabled(true);

	auto panel = new QWidget(this);
	auto grid = new QGridLayout(panel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	static const QString comps[] = { "x", "y", "z" };

	m_editJ = new QLineEdit(panel);
	m_editJ->setToolTip("Exchange constant J.");

	for(int i = 0; i < 3; ++i)
	{
		m_editDMI[i] = new QLineEdit(panel);
		m_editDMI[i]->setToolTip(QString("DMI %1 component.").arg(comps[i]));
	}

	for(int i = 0; i < 3*3; ++i)
	{
		m_editJs[i] = new QLineEdit(panel);
		m_editJs[i]->setToolTip(QString("J_%1%2 component of the general exchange matrix.")
		 .arg(comps[i/3]).arg(comps[i%3]));
	}

	m_symmidx = new QSpinBox(panel);
	m_symmidx->setMinimum(0);  // 0 is the invalid index
	m_symmidx->setValue(1);
	m_symmidx->setToolTip("Index of the symmetry group to assign.");

	QPushButton *btnAssignByIdx = new QPushButton("Assign", panel);
	btnAssignByIdx->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	btnAssignByIdx->setToolTip("Assign the given values to all couplings in the same symmetry group.");


	int y = 0;
	grid->addWidget(new QLabel("J:", panel), y, 0, 1, 1);
	grid->addWidget(m_editJ, y++, 1, 1, 1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 4);

	grid->addWidget(new QLabel("DMI Vector:", panel), y, 0, 1, 1);
	for(int i = 0; i < 3; ++i)
		grid->addWidget(m_editDMI[i], y, 1 + i, 1, 1);
	++y;

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 4);

	grid->addWidget(new QLabel("J Matrix:", panel), y, 0, 1, 1);
	for(int i = 0; i < 3*3; ++i)
		grid->addWidget(m_editJs[i], y + i/3, 1 + i%3, 1, 1);
	y += 3;

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 4);
	auto sep1 = new QFrame(panel);
	sep1->setFrameStyle(QFrame::HLine);
	grid->addWidget(sep1, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 4);

	grid->addWidget(new QLabel("Symmetry Index:", panel), y, 0, 1, 1);
	grid->addWidget(m_symmidx, y, 1, 1, 1);
	grid->addWidget(btnAssignByIdx, y++, 3, 1, 1);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++, 0, 1, 4);

	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &AssignDlg::accept);

	auto dlgGrid = new QGridLayout(this);
	dlgGrid->setSpacing(4);
	dlgGrid->setContentsMargins(8, 8, 8, 8);
	dlgGrid->addWidget(panel, 0, 0, 1, 4);
	dlgGrid->addWidget(btnbox, 1, 3, 1, 1);


	// connections
	connect(btnAssignByIdx, &QAbstractButton::clicked, this, &AssignDlg::AssignByIndex);


	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("assign/geo"))
			restoreGeometry(m_sett->value("assign/geo").toByteArray());
		else
			resize(400, 200);
	}
}



/**
 * emits coupling to assign
 */
void AssignDlg::AssignByIndex()
{
	t_size symmidx = m_symmidx->value();
	std::string J = m_editJ->text().toStdString();
	std::string DMI[3];
	std::string Js[3*3];

	for(int i = 0; i < 3*3; ++i)
	{
		if(i < 3)
			DMI[i] = m_editDMI[i]->text().toStdString();
		Js[i] = m_editJs[i]->text().toStdString();
	}

	emit AssignCouplingsBySymmetryIndex(symmidx, &J, DMI, Js);
}



/**
 * close the dialog
 */
void AssignDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("assign/geo", saveGeometry());
	}

	QDialog::accept();
}
