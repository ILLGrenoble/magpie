/**
 * magnon dynamics -- transformation calculator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 29-dec-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 7-sep-2023 from my privately developed "gl" project: https://github.com/t-weber/gl .
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "gl" project
 * Copyright (C) 2021-2023  Tobias WEBER (privately developed).
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

#include "trafos.h"
#include "defs.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QLabel>

#include <boost/math/quaternion.hpp>



TrafoCalculator::TrafoCalculator(QWidget* pParent, QSettings *sett)
	: QDialog{pParent}, m_sett(sett)
{
	setWindowTitle("Transformations");
	setSizeGripEnabled(true);

	// tabs
	QTabWidget *tabs = new QTabWidget(this);
	QWidget *rotationPanel = new QWidget(tabs);

	// buttons
	QDialogButtonBox *buttons = new QDialogButtonBox(this);
	buttons->setStandardButtons(QDialogButtonBox::Ok);

	// tab panels
	tabs->addTab(rotationPanel, "Axis Rotation");

	// rotation tab
	QLabel *labelAxis = new QLabel("Axis: ");
	QLabel *labelAngle = new QLabel("Angle (\xc2\xb0): ");
	QLabel *labelVecToRotate = new QLabel("Vector: ");

	m_spinAxis[0] = new QDoubleSpinBox(rotationPanel);
	m_spinAxis[1] = new QDoubleSpinBox(rotationPanel);
	m_spinAxis[2] = new QDoubleSpinBox(rotationPanel);
	m_spinAxis[0]->setValue(0);
	m_spinAxis[1]->setValue(0);
	m_spinAxis[2]->setValue(1);

	m_spinAngle = new QDoubleSpinBox(rotationPanel);
	m_spinAngle->setMinimum(-360.);
	m_spinAngle->setMaximum(360);
	m_spinAngle->setDecimals(3);
	m_spinAngle->setSingleStep(0.1);
	//m_spinAngle->setSuffix("\xc2\xb0");

	m_spinVecToRotate[0] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[1] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[2] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[0]->setValue(1);
	m_spinVecToRotate[1]->setValue(0);
	m_spinVecToRotate[2]->setValue(0);

	for(int i=0; i<3; ++i)
	{
		m_spinAxis[i]->setMinimum(-999.);
		m_spinAxis[i]->setMaximum(999.);
		m_spinAxis[i]->setDecimals(4);
		m_spinAxis[i]->setSingleStep(0.1);

		m_spinVecToRotate[i]->setMinimum(-999.);
		m_spinVecToRotate[i]->setMaximum(999.);
		m_spinVecToRotate[i]->setDecimals(4);
		m_spinVecToRotate[i]->setSingleStep(0.1);
	}

	labelAxis->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	labelAngle->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	labelVecToRotate->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

	m_textRotation = new QTextEdit(rotationPanel);
	m_textRotation->setReadOnly(true);

	// rotation grid
	auto grid_rotation = new QGridLayout(rotationPanel);
	grid_rotation->setSpacing(4);
	grid_rotation->setContentsMargins(6, 6, 6, 6);
	grid_rotation->addWidget(labelAxis, 0, 0, 1, 1);
	grid_rotation->addWidget(m_spinAxis[0], 0, 1, 1, 1);
	grid_rotation->addWidget(m_spinAxis[1], 0, 2, 1, 1);
	grid_rotation->addWidget(m_spinAxis[2], 0, 3, 1, 1);
	grid_rotation->addWidget(labelAngle, 1, 0, 1, 1);
	grid_rotation->addWidget(m_spinAngle, 1, 1, 1, 1);
	grid_rotation->addWidget(labelVecToRotate, 2, 0, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[0], 2, 1, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[1], 2, 2, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[2], 2, 3, 1, 1);
	grid_rotation->addWidget(m_textRotation, 3, 0, 1, 4);

	// main grid
	auto grid_dlg = new QGridLayout(this);
	grid_dlg->setSpacing(4);
	grid_dlg->setContentsMargins(8, 8, 8, 8);
	grid_dlg->addWidget(tabs, 0, 0, 1, 1);
	grid_dlg->addWidget(buttons, 1, 0, 1, 1);

	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("trafocalc/geo"))
			restoreGeometry(m_sett->value("trafocalc/geo").toByteArray());
		else
			resize(500, 500);
	}

	// connections
	connect(buttons, &QDialogButtonBox::accepted,
		this, &TrafoCalculator::accept);
	connect(buttons, &QDialogButtonBox::rejected,
		this, &TrafoCalculator::reject);
	for(QDoubleSpinBox* spin : {
		m_spinAxis[0], m_spinAxis[1], m_spinAxis[2], m_spinAngle,
		m_spinVecToRotate[0], m_spinVecToRotate[1], m_spinVecToRotate[2] })
	{
		connect(spin,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &TrafoCalculator::CalculateRotation);
	}

	CalculateRotation();
}



void TrafoCalculator::CalculateRotation()
{
	using namespace tl2_ops;

	if(!m_spinAngle || !m_textRotation)
		return;

	t_vec_real axis = tl2::create<t_vec_real>({
		(t_real)m_spinAxis[0]->value(),
		(t_real)m_spinAxis[1]->value(),
		(t_real)m_spinAxis[2]->value() });
	t_real angle = tl2::d2r<t_real>(m_spinAngle->value());
	t_vec_real vec = tl2::create<t_vec_real>({
		(t_real)m_spinVecToRotate[0]->value(),
		(t_real)m_spinVecToRotate[1]->value(),
		(t_real)m_spinVecToRotate[2]->value() });

	m_textRotation->clear();

	t_mat_real mat = tl2::rotation<t_mat_real, t_vec_real>(axis, angle, false);
	tl2::set_eps_0(mat, g_eps);


	// print the rotation matrix
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	ostrResult << "<p>Transformation Matrix:\n";
	ostrResult << "<table style=\"border:0px\">\n";
	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		ostrResult << "\t<tr>\n";
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			ostrResult << "\t\t<td style=\"padding-right:8px\">";
			ostrResult << mat(i, j);
			ostrResult << "</td>\n";
		}
		ostrResult << "\t</tr>\n";
	}

	ostrResult << "</table>";
	ostrResult << "</p>\n";

	ostrResult << "<p>As Single-Line String:<br>";
	ostrResult << mat;
	ostrResult << "</p>\n";


	using t_quat = boost::math::quaternion<t_real>;

	t_quat quat = tl2::rot3_to_quat<t_mat_real, t_quat>(mat);
	tl2::set_eps_0(quat, g_eps);
	ostrResult << "<p>As Quaternion:<br>";
	ostrResult << quat;
	ostrResult << "</p>\n";


	// print the rotated test vector
	t_vec_real vec_rot = mat * vec;
	tl2::set_eps_0(vec_rot, g_eps);
	ostrResult << "<p>Rotated Vector:<br>";
	ostrResult << vec_rot;
	ostrResult << "</p>\n";

	m_textRotation->setHtml(ostrResult.str().c_str());
}



/**
 * close the dialog
 */
void TrafoCalculator::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("trafocalc/geo", saveGeometry());
	}

	QDialog::accept();
}
