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
 * exports data to different file types
 */
void MagDynDlg::CreateExportPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_exportpanel = new QWidget(this);

	// Q coordinates
	m_exportStartQ[0] = new QDoubleSpinBox(m_exportpanel);
	m_exportStartQ[1] = new QDoubleSpinBox(m_exportpanel);
	m_exportStartQ[2] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[0] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[1] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[2] = new QDoubleSpinBox(m_exportpanel);

	// number of grid points
	for(int i = 0; i < 3; ++i)
	{
		m_exportNumPoints[i] = new QSpinBox(m_exportpanel);
		m_exportNumPoints[i]->setMinimum(1);
		m_exportNumPoints[i]->setMaximum(99999);
		m_exportNumPoints[i]->setValue(128);
		m_exportNumPoints[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	// export
	m_exportFormat = new QComboBox(m_exportpanel);
	m_exportFormat->addItem("Takin Grid File", EXPORT_GRID);
#ifdef USE_HDF5
	m_exportFormat->addItem("HDF5 Data File", EXPORT_HDF5);
#endif
	m_exportFormat->addItem("Text Data File", EXPORT_TEXT);

	QPushButton *btn_export = new QPushButton(
		QIcon::fromTheme("document-save-as"),
		"Export...", m_exportpanel);
	btn_export->setFocusPolicy(Qt::StrongFocus);

	for(int i = 0; i < 3; ++i)
	{
		m_exportStartQ[i]->setDecimals(4);
		m_exportStartQ[i]->setMinimum(-99.9999);
		m_exportStartQ[i]->setMaximum(+99.9999);
		m_exportStartQ[i]->setSingleStep(0.01);
		m_exportStartQ[i]->setValue(-1.);
		m_exportStartQ[i]->setSuffix(" rlu");
		m_exportStartQ[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_exportStartQ[i]->setPrefix(hklPrefix[i]);

		m_exportEndQ[i]->setDecimals(4);
		m_exportEndQ[i]->setMinimum(-99.9999);
		m_exportEndQ[i]->setMaximum(+99.9999);
		m_exportEndQ[i]->setSingleStep(0.01);
		m_exportEndQ[i]->setValue(1.);
		m_exportEndQ[i]->setSuffix(" rlu");
		m_exportEndQ[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_exportEndQ[i]->setPrefix(hklPrefix[i]);
	}


	QGridLayout *grid = new QGridLayout(m_exportpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Export Ranges:", m_exportpanel), y++,0,1,4);
	grid->addWidget(new QLabel("Start Q:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportStartQ[0], y,1,1,1);
	grid->addWidget(m_exportStartQ[1], y,2,1,1);
	grid->addWidget(m_exportStartQ[2], y++,3,1,1);
	grid->addWidget(new QLabel("End Q:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportEndQ[0], y,1,1,1);
	grid->addWidget(m_exportEndQ[1], y,2,1,1);
	grid->addWidget(m_exportEndQ[2], y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	QFrame *sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	grid->addWidget(sep1, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Number of Grid Points per Q Direction:",
		m_exportpanel), y++,0,1,4);
	grid->addWidget(new QLabel("Points:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportNumPoints[0], y,1,1,1);
	grid->addWidget(m_exportNumPoints[1], y,2,1,1);
	grid->addWidget(m_exportNumPoints[2], y++,3,1,1);

	QFrame *sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	QLabel* labelBoseInfo = new QLabel(QString(
		"Info: If this grid file is to be used in Takin's "
		"resolution convolution module (\"Model Source: Uniform Grid\"), "
		"please disable the Bose factor (\"Calculation Options\" -> \"Use Bose Factor\" [off]). "
		"The Bose factor is already managed by the convolution module."),
		m_exportpanel);
	labelBoseInfo->setWordWrap(true);
	grid->addWidget(labelBoseInfo, y++,0,1,4);

	QFrame *sep3 = new QFrame(m_sampleenviropanel);
	sep3->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep3, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0,1,4);

	grid->addWidget(new QLabel("Export Format:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportFormat, y,1,1,1);
	grid->addWidget(btn_export, y++,3,1,1);

	// signals
	connect(btn_export, &QAbstractButton::clicked, this,
		static_cast<void (MagDynDlg::*)()>(&MagDynDlg::ExportSQE));


	m_tabs_out->addTab(m_exportpanel, "Export Data");
}
