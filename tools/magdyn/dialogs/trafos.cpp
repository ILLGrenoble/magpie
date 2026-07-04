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



// set kernel from the main window
void TrafoCalculator::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;

	CalculateRotation();
	CalculateProjection();
}



TrafoCalculator::TrafoCalculator(QWidget* pParent, QSettings *sett)
	: QDialog{pParent}, m_sett(sett)
{
	setWindowTitle("Transformations");
	setSizeGripEnabled(true);

	// tabs
	QTabWidget *tabs = new QTabWidget(this);
	QWidget *rotationPanel = CreateRotationPanel();
	QWidget *projectionPanel = CreateProjectionPanel();
	rotationPanel->setParent(tabs);
	projectionPanel->setParent(tabs);

	// buttons
	QDialogButtonBox *buttons = new QDialogButtonBox(this);
	buttons->setStandardButtons(QDialogButtonBox::Ok);

	// tab panels
	tabs->addTab(rotationPanel, "Axis Rotation");
	tabs->addTab(projectionPanel, "Projection");

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
	connect(buttons, &QDialogButtonBox::accepted, this, &TrafoCalculator::accept);
	connect(buttons, &QDialogButtonBox::rejected, this, &TrafoCalculator::reject);

	CalculateRotation();
	CalculateProjection();
}



QWidget* TrafoCalculator::CreateRotationPanel()
{
	// rotation tab (crystal)
	QWidget *rotationPanel = new QWidget(this);

	QLabel *labelAxis = new QLabel("Axis (rlu): ");
	QLabel *labelAngle = new QLabel("Angle (\xc2\xb0): ");
	QLabel *labelVecToRotate = new QLabel("Vector (rlu): ");

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

	m_checkRot = new QCheckBox(rotationPanel);
	m_checkRot->setText("Use Crystal System");
	m_checkRot->setChecked(true);

	QPushButton *btnRecalc = new QPushButton(rotationPanel);
	btnRecalc->setText("Get Crystal");

	m_spinVecToRotate[0] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[1] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[2] = new QDoubleSpinBox(rotationPanel);
	m_spinVecToRotate[0]->setValue(1);
	m_spinVecToRotate[1]->setValue(0);
	m_spinVecToRotate[2]->setValue(0);

	for(int i = 0; i < 3; ++i)
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
	grid_rotation->addWidget(m_checkRot, 1, 2, 1, 1);
	grid_rotation->addWidget(btnRecalc, 1, 3, 1, 1);
	grid_rotation->addWidget(labelVecToRotate, 2, 0, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[0], 2, 1, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[1], 2, 2, 1, 1);
	grid_rotation->addWidget(m_spinVecToRotate[2], 2, 3, 1, 1);
	grid_rotation->addWidget(m_textRotation, 3, 0, 1, 4);

	// connections
	for(QDoubleSpinBox* spin : {
		m_spinAxis[0], m_spinAxis[1], m_spinAxis[2], m_spinAngle,
		m_spinVecToRotate[0], m_spinVecToRotate[1], m_spinVecToRotate[2] })
	{
		connect(spin,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &TrafoCalculator::CalculateRotation);
	}
	connect(m_checkRot, &QCheckBox::toggled, [this, btnRecalc](bool checked)
	{
		btnRecalc->setEnabled(checked);
		CalculateRotation();
	});
	connect(btnRecalc, &QAbstractButton::clicked,
		this, &TrafoCalculator::CalculateRotation);

	return rotationPanel;
}



QWidget* TrafoCalculator::CreateProjectionPanel()
{
	// projection tab (crystal)
	QWidget *projectionPanel = new QWidget(this);

	QLabel *labelAxis = new QLabel("Axis (rlu): ");

	m_spinProjAxis[0] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[1] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[2] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[0]->setValue(0);
	m_spinProjAxis[1]->setValue(0);
	m_spinProjAxis[2]->setValue(1);

	m_checkProj = new QCheckBox(projectionPanel);
	m_checkProj->setText("Use Crystal System");
	m_checkProj->setChecked(true);

	QPushButton *btnRecalc = new QPushButton(projectionPanel);
	btnRecalc->setText("Get Crystal");

	for(int i = 0; i < 3; ++i)
	{
		m_spinProjAxis[i]->setMinimum(-999.);
		m_spinProjAxis[i]->setMaximum(999.);
		m_spinProjAxis[i]->setDecimals(4);
		m_spinProjAxis[i]->setSingleStep(0.1);
	}

	labelAxis->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

	m_textProjection = new QTextEdit(projectionPanel);
	m_textProjection->setReadOnly(true);

	// rotation grid
	auto grid_projection = new QGridLayout(projectionPanel);
	grid_projection->setSpacing(4);
	grid_projection->setContentsMargins(6, 6, 6, 6);
	grid_projection->addWidget(labelAxis, 0, 0, 1, 1);
	grid_projection->addWidget(m_spinProjAxis[0], 0, 1, 1, 1);
	grid_projection->addWidget(m_spinProjAxis[1], 0, 2, 1, 1);
	grid_projection->addWidget(m_spinProjAxis[2], 0, 3, 1, 1);
	grid_projection->addWidget(m_checkProj, 1, 2, 1, 1);
	grid_projection->addWidget(btnRecalc, 1, 3, 1, 1);
	grid_projection->addWidget(m_textProjection, 2, 0, 1, 4);

	// connections
	for(QDoubleSpinBox* spin : { m_spinProjAxis[0], m_spinProjAxis[1], m_spinProjAxis[2] })
	{
		connect(spin,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &TrafoCalculator::CalculateProjection);
	}
	connect(m_checkProj, &QCheckBox::toggled, [this, btnRecalc](bool checked)
	{
		btnRecalc->setEnabled(checked);
		CalculateProjection();
	});
	connect(btnRecalc, &QAbstractButton::clicked,
		this, &TrafoCalculator::CalculateProjection);

	return projectionPanel;
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


	// apply crystal B matrix
	bool use_B = m_checkRot->isChecked() && m_dyn;

	t_mat_real xtalB_inv;
	bool inv_ok = false;
	if(use_B)
	{
		const t_mat_real& xtalB = m_dyn->GetCrystalBTrafo();
		std::tie(xtalB_inv, inv_ok) = tl2::inv(xtalB);
		if(inv_ok)
		{
			axis = xtalB * axis;
			vec = xtalB * vec;
		}
	}

	t_mat_real mat = tl2::rotation<t_mat_real, t_vec_real>(axis, angle, false);
	tl2::set_eps_0(mat, g_eps);


	// print the B and the rotation matrices
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	if(use_B)
	{
		t_mat_real xtalB = m_dyn->GetCrystalBTrafo();
		tl2::set_eps_0(xtalB, g_eps);

		ostrResult << "<p>Crystal B Matrix:\n";
		ostrResult << "<table style=\"border:0px\">\n";
		for(std::size_t i = 0; i < xtalB.size1(); ++i)
		{
			ostrResult << "\t<tr>\n";
			for(std::size_t j = 0; j < xtalB.size2(); ++j)
			{
				ostrResult << "\t\t<td style=\"padding-right:8px\">";
				ostrResult << xtalB(i, j);
				ostrResult << "</td>\n";
			}
			ostrResult << "\t</tr>\n";
		}
		ostrResult << "</table>";
		ostrResult << "</p>\n";

		ostrResult << "<p>B Matrix As Single-Line String:<br>";
		ostrResult << xtalB;
		ostrResult << "</p>\n";
	}

	ostrResult << "<p>Transformation Matrix:\n";
	ostrResult << "<table style=\"border:0px\">\n";
	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		ostrResult << "\t<tr>\n";
		for(std::size_t j = 0; j < mat.size2(); ++j)
		{
			ostrResult << "\t\t<td style=\"padding-right:8px\">";
			ostrResult << mat(i, j);
			ostrResult << "</td>\n";
		}
		ostrResult << "\t</tr>\n";
	}
	ostrResult << "</table>";
	ostrResult << "</p>\n";

	ostrResult << "<p>Trafo As Single-Line String:<br>";
	ostrResult << mat;
	ostrResult << "</p>\n";


	using t_quat = boost::math::quaternion<t_real>;

	t_quat quat = tl2::rot3_to_quat<t_mat_real, t_quat>(mat);
	tl2::set_eps_0(quat, g_eps);
	ostrResult << "<p>Trafo As Quaternion:<br>";
	ostrResult << quat;
	ostrResult << "</p>\n";

	if(use_B && inv_ok)
	{
		tl2::set_eps_0(axis, g_eps);
		ostrResult << "<p>Original Axis (lab): ";
		ostrResult << axis;
		ostrResult << "</p>\n";

		tl2::set_eps_0(vec, g_eps);
		ostrResult << "<p>Original Vector (lab): ";
		ostrResult << vec;
		ostrResult << "</p>\n";
	}

	// print the rotated test vector
	t_vec_real vec_rot = mat * vec;

	tl2::set_eps_0(vec_rot, g_eps);
	ostrResult << "<p>Rotated Vector (rlu): ";
	ostrResult << vec_rot;
	ostrResult << "</p>\n";

	if(use_B && inv_ok)
	{	
		vec_rot = xtalB_inv * vec_rot;
		tl2::set_eps_0(vec_rot, g_eps);
		ostrResult << "<p>Rotated Vector (lab): ";
		ostrResult << vec_rot;
		ostrResult << "</p>\n";
	}


	m_textRotation->setHtml(ostrResult.str().c_str());
}



void TrafoCalculator::CalculateProjection()
{
	using namespace tl2_ops;

	if(!m_textProjection)
		return;

	t_vec_real axis = tl2::create<t_vec_real>({
		(t_real)m_spinProjAxis[0]->value(),
		(t_real)m_spinProjAxis[1]->value(),
		(t_real)m_spinProjAxis[2]->value() });

	m_textProjection->clear();


	// apply crystal B matrix
	bool use_B = m_checkProj->isChecked() && m_dyn;
	if(use_B)
	{
		const t_mat_real& xtalB = m_dyn->GetCrystalBTrafo();
		axis = xtalB * axis;
	}

	t_mat_real matProj = tl2::projector<t_mat_real, t_vec_real>(axis, false);
	t_mat_real matOrthoProj = tl2::ortho_projector<t_mat_real, t_vec_real>(axis, false);
	tl2::set_eps_0(matProj, g_eps);
	tl2::set_eps_0(matOrthoProj, g_eps);


	// print the B and the rotation matrices
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	if(use_B)
	{
		t_mat_real xtalB = m_dyn->GetCrystalBTrafo();
		tl2::set_eps_0(xtalB, g_eps);


		ostrResult << "<p>Crystal B Matrix:\n";
		ostrResult << "<table style=\"border:0px\">\n";
		for(std::size_t i = 0; i < xtalB.size1(); ++i)
		{
			ostrResult << "\t<tr>\n";
			for(std::size_t j = 0; j < xtalB.size2(); ++j)
			{
				ostrResult << "\t\t<td style=\"padding-right:8px\">";
				ostrResult << xtalB(i, j);
				ostrResult << "</td>\n";
			}
			ostrResult << "\t</tr>\n";
		}
		ostrResult << "</table>";
		ostrResult << "</p>\n";

		ostrResult << "<p>B Matrix As Single-Line String:<br>";
		ostrResult << xtalB;
		ostrResult << "</p>\n";
	}


	ostrResult << "<p>Projection Matrix:\n";
	ostrResult << "<table style=\"border:0px\">\n";
	for(std::size_t i = 0; i < matProj.size1(); ++i)
	{
		ostrResult << "\t<tr>\n";
		for(std::size_t j = 0; j < matProj.size2(); ++j)
		{
			ostrResult << "\t\t<td style=\"padding-right:8px\">";
			ostrResult << matProj(i, j);
			ostrResult << "</td>\n";
		}
		ostrResult << "\t</tr>\n";
	}
	ostrResult << "</table>";
	ostrResult << "</p>\n";

	ostrResult << "<p>Projection Matrix Trace: ";
	ostrResult << tl2::trace(matProj);
	ostrResult << "</p>\n";

	ostrResult << "<p>Projection Matrix As Single-Line String:<br>";
	ostrResult << matProj;
	ostrResult << "</p>\n";


	ostrResult << "<p>Orthogonal Projection Matrix:\n";
	ostrResult << "<table style=\"border:0px\">\n";
	for(std::size_t i = 0; i < matOrthoProj.size1(); ++i)
	{
		ostrResult << "\t<tr>\n";
		for(std::size_t j = 0; j < matOrthoProj.size2(); ++j)
		{
			ostrResult << "\t\t<td style=\"padding-right:8px\">";
			ostrResult << matOrthoProj(i, j);
			ostrResult << "</td>\n";
		}
		ostrResult << "\t</tr>\n";
	}
	ostrResult << "</table>";
	ostrResult << "</p>\n";

	ostrResult << "<p>Orthogonal Projection Matrix Trace: ";
	ostrResult << tl2::trace(matOrthoProj);
	ostrResult << "</p>\n";

	ostrResult << "<p>Orthogonal Projection Matrix As Single-Line String:<br>";
	ostrResult << matOrthoProj;
	ostrResult << "</p>\n";


	tl2::set_eps_0(axis, g_eps);
	ostrResult << "<p>Original Axis (lab): ";
	ostrResult << axis;
	ostrResult << "</p>\n";


	m_textProjection->setHtml(ostrResult.str().c_str());
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
