/**
 * magnon dynamics -- info dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date june-2024
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

#include "glinfos.h"
#include "defs.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>



/**
 * set up the gui
 */
GlInfoDlg::GlInfoDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	setWindowTitle("About Renderer");
	setSizeGripEnabled(true);

	auto infopanel = new QWidget(this);
	auto grid = new QGridLayout(infopanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	// renderer infos
	for(int i = 0; i < 4; ++i)
	{
		m_labelGlInfos[i] = new QLabel("", infopanel);
		m_labelGlInfos[i]->setSizePolicy(
			QSizePolicy::Ignored,
			m_labelGlInfos[i]->sizePolicy().verticalPolicy());

		if(i == 0)
			m_labelGlInfos[i]->setText("GL renderer not yet initialised.");
	}

	int y = 0;
	for(int i = 0; i < 4; ++i)
		grid->addWidget(m_labelGlInfos[i], y++,0, 1, 1);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0, 1,1);

	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &GlInfoDlg::accept);

	auto dlgGrid = new QGridLayout(this);
	dlgGrid->setSpacing(4);
	dlgGrid->setContentsMargins(8, 8, 8, 8);
	dlgGrid->addWidget(infopanel, 0, 0, 1, 4);
	dlgGrid->addWidget(btnbox, 1, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("glinfos/geo"))
			restoreGeometry(m_sett->value("glinfos/geo").toByteArray());
		else
			resize(400, 200);
	}
}



void GlInfoDlg::SetGlDeviceInfos(const std::string& ver, const std::string& shader_ver,
	const std::string& vendor, const std::string& renderer)
{
	m_labelGlInfos[0]->setText(QString("GL Version: %1.").arg(ver.c_str()));
	m_labelGlInfos[1]->setText(QString("GL Shader Version: %1.").arg(shader_ver.c_str()));
	m_labelGlInfos[2]->setText(QString("GL Vendor: %1.").arg(vendor.c_str()));
	m_labelGlInfos[3]->setText(QString("GL Device: %1.").arg(renderer.c_str()));
}



/**
 * close the dialog
 */
void GlInfoDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("glinfos/geo", saveGeometry());
	}

	QDialog::accept();
}
