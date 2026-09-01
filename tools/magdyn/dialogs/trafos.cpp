/**
 * magnon dynamics -- transformation calculator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 29-dec-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 7-sep-2023 from my privately developed "gl" project: https://github.com/t-weber/gl .
 *
 * ----------------------------------------------------------------------------
 * magpie & mag-core
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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
	CalculateCrossProduct();
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
	QWidget *crossProdPanel = CreateCrossProductPanel();
	rotationPanel->setParent(tabs);
	projectionPanel->setParent(tabs);
	crossProdPanel->setParent(tabs);

	// buttons
	QDialogButtonBox *buttons = new QDialogButtonBox(this);
	buttons->setStandardButtons(QDialogButtonBox::Ok);

	// tab panels
	tabs->addTab(rotationPanel, "Axis Rotation");
	tabs->addTab(projectionPanel, "Projection");
	tabs->addTab(crossProdPanel, "Normal");

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

	m_checkXtalRot = new QCheckBox(rotationPanel);
	m_checkXtalRot->setText("Use Crystal System");
	m_checkXtalRot->setToolTip("Use the crystallographic B matrix.");
	m_checkXtalRot->setChecked(true);

	QPushButton *btnRecalc = new QPushButton(rotationPanel);
	btnRecalc->setText("Get Crystal");
	btnRecalc->setToolTip("Get the crystallographic B matrix and recalculate.");

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
	grid_rotation->addWidget(m_checkXtalRot, 1, 2, 1, 1);
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
	connect(m_checkXtalRot, &QCheckBox::toggled, [this, btnRecalc](bool checked)
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

	m_spinProjAxis[0] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[1] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[2] = new QDoubleSpinBox(projectionPanel);
	m_spinProjAxis[0]->setValue(0);
	m_spinProjAxis[1]->setValue(0);
	m_spinProjAxis[2]->setValue(1);

	m_spinVecToProj[0] = new QDoubleSpinBox(projectionPanel);
	m_spinVecToProj[1] = new QDoubleSpinBox(projectionPanel);
	m_spinVecToProj[2] = new QDoubleSpinBox(projectionPanel);
	m_spinVecToProj[0]->setValue(1);
	m_spinVecToProj[1]->setValue(0);
	m_spinVecToProj[2]->setValue(0);

	for(int i = 0; i < 3; ++i)
	{
		m_spinProjAxis[i]->setMinimum(-999.);
		m_spinProjAxis[i]->setMaximum(999.);
		m_spinProjAxis[i]->setDecimals(4);
		m_spinProjAxis[i]->setSingleStep(0.1);

		m_spinVecToProj[i]->setMinimum(-999.);
		m_spinVecToProj[i]->setMaximum(999.);
		m_spinVecToProj[i]->setDecimals(4);
		m_spinVecToProj[i]->setSingleStep(0.1);
	}

	m_checkXtalProj = new QCheckBox(projectionPanel);
	m_checkXtalProj->setText("Use Crystal System");
	m_checkXtalProj->setToolTip("Use the crystallographic B matrix.");
	m_checkXtalProj->setChecked(true);

	QPushButton *btnRecalc = new QPushButton(projectionPanel);
	btnRecalc->setText("Get Crystal");
	btnRecalc->setToolTip("Get the crystallographic B matrix and recalculate.");

	QLabel *labelAxis = new QLabel("Axis (rlu): ");
	QLabel *labelVecToProj = new QLabel("Vector (rlu): ");
	labelAxis->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	labelVecToProj->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

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
	grid_projection->addWidget(labelVecToProj, 1, 0, 1, 1);
	grid_projection->addWidget(m_spinVecToProj[0], 1, 1, 1, 1);
	grid_projection->addWidget(m_spinVecToProj[1], 1, 2, 1, 1);
	grid_projection->addWidget(m_spinVecToProj[2], 1, 3, 1, 1);
	grid_projection->addWidget(m_checkXtalProj, 2, 2, 1, 1);
	grid_projection->addWidget(btnRecalc, 2, 3, 1, 1);
	grid_projection->addWidget(m_textProjection, 3, 0, 1, 4);

	// connections
	for(QDoubleSpinBox* spin : {
		m_spinProjAxis[0], m_spinProjAxis[1], m_spinProjAxis[2],
		m_spinVecToProj[0], m_spinVecToProj[1], m_spinVecToProj[2] })
	{
		connect(spin,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &TrafoCalculator::CalculateProjection);
	}
	connect(m_checkXtalProj, &QCheckBox::toggled, [this, btnRecalc](bool checked)
	{
		btnRecalc->setEnabled(checked);
		CalculateProjection();
	});
	connect(btnRecalc, &QAbstractButton::clicked,
		this, &TrafoCalculator::CalculateProjection);

	return projectionPanel;
}



QWidget* TrafoCalculator::CreateCrossProductPanel()
{
	// projection tab (crystal)
	QWidget *crossProdPanel = new QWidget(this);

	m_spinVec1[0] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec1[1] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec1[2] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec1[0]->setValue(1);
	m_spinVec1[1]->setValue(0);
	m_spinVec1[2]->setValue(0);

	m_spinVec2[0] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec2[1] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec2[2] = new QDoubleSpinBox(crossProdPanel);
	m_spinVec2[0]->setValue(0);
	m_spinVec2[1]->setValue(1);
	m_spinVec2[2]->setValue(0);

	for(int i = 0; i < 3; ++i)
	{
		m_spinVec1[i]->setMinimum(-999.);
		m_spinVec1[i]->setMaximum(999.);
		m_spinVec1[i]->setDecimals(4);
		m_spinVec1[i]->setSingleStep(0.1);

		m_spinVec2[i]->setMinimum(-999.);
		m_spinVec2[i]->setMaximum(999.);
		m_spinVec2[i]->setDecimals(4);
		m_spinVec2[i]->setSingleStep(0.1);
	}

	m_checkXtalCrossProd = new QCheckBox(crossProdPanel);
	m_checkXtalCrossProd->setText("Use Crystal System");
	m_checkXtalCrossProd->setToolTip("Use the crystallographic B matrix.");
	m_checkXtalCrossProd->setChecked(true);

	QPushButton *btnRecalc = new QPushButton(crossProdPanel);
	btnRecalc->setText("Get Crystal");
	btnRecalc->setToolTip("Get the crystallographic B matrix and recalculate.");

	QLabel *labelVec1 = new QLabel("Vector 1 (rlu): ");
	QLabel *labelVec2 = new QLabel("Vector 2 (rlu): ");
	labelVec1->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
	labelVec2->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

	m_textCrossProd = new QTextEdit(crossProdPanel);
	m_textCrossProd->setReadOnly(true);

	// rotation grid
	auto grid_projection = new QGridLayout(crossProdPanel);
	grid_projection->setSpacing(4);
	grid_projection->setContentsMargins(6, 6, 6, 6);
	grid_projection->addWidget(labelVec1, 0, 0, 1, 1);
	grid_projection->addWidget(m_spinVec1[0], 0, 1, 1, 1);
	grid_projection->addWidget(m_spinVec1[1], 0, 2, 1, 1);
	grid_projection->addWidget(m_spinVec1[2], 0, 3, 1, 1);
	grid_projection->addWidget(labelVec2, 1, 0, 1, 1);
	grid_projection->addWidget(m_spinVec2[0], 1, 1, 1, 1);
	grid_projection->addWidget(m_spinVec2[1], 1, 2, 1, 1);
	grid_projection->addWidget(m_spinVec2[2], 1, 3, 1, 1);
	grid_projection->addWidget(m_checkXtalCrossProd, 2, 2, 1, 1);
	grid_projection->addWidget(btnRecalc, 2, 3, 1, 1);
	grid_projection->addWidget(m_textCrossProd, 3, 0, 1, 4);

	// connections
	for(QDoubleSpinBox* spin : {
		m_spinVec1[0], m_spinVec1[1], m_spinVec1[2],
		m_spinVec2[0], m_spinVec2[1], m_spinVec2[2] })
	{
		connect(spin,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &TrafoCalculator::CalculateCrossProduct);
	}
	connect(m_checkXtalCrossProd, &QCheckBox::toggled, [this, btnRecalc](bool checked)
	{
		btnRecalc->setEnabled(checked);
		CalculateCrossProduct();
	});
	connect(btnRecalc, &QAbstractButton::clicked, this, &TrafoCalculator::CalculateCrossProduct);

	return crossProdPanel;
}



/**
 * prints a matrix as html table
 */
static void print_matrix(const t_mat33_real& _mat, const std::string& name,
	std::ostream& ostr = std::cout)
{
	using namespace tl2_ops;

	t_mat33_real mat = _mat;
	tl2::set_eps_0(mat, g_eps);

	ostr << "<p>" << name << ":\n";
	ostr << "<table style=\"border:0px\">\n";
	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		ostr << "\t<tr>\n";
		for(std::size_t j = 0; j < mat.size2(); ++j)
		{
			ostr << "\t\t<td style=\"padding-right:8px\">";
			ostr << mat(i, j);
			ostr << "</td>\n";
		}
		ostr << "\t</tr>\n";
	}
	ostr << "</table>";
	ostr << "</p>\n";

	ostr << "<p>" << name << " Trace: ";
	ostr << tl2::trace(mat);
	ostr << "<br>" << name << " As Single-Line String:<br>";
	ostr << mat;
	ostr << "</p>\n";
}



void TrafoCalculator::CalculateRotation()
{
	using namespace tl2_ops;

	if(!m_spinAngle || !m_textRotation)
		return;

	t_vec3_real axis = tl2::create<t_vec3_real>({
		(t_real)m_spinAxis[0]->value(),
		(t_real)m_spinAxis[1]->value(),
		(t_real)m_spinAxis[2]->value() });
	t_real angle = tl2::d2r<t_real>(m_spinAngle->value());
	t_vec3_real vec = tl2::create<t_vec3_real>({
		(t_real)m_spinVecToRotate[0]->value(),
		(t_real)m_spinVecToRotate[1]->value(),
		(t_real)m_spinVecToRotate[2]->value() });

	m_textRotation->clear();


	// apply crystal B matrix
	bool use_B = m_checkXtalRot->isChecked() && m_dyn;
	t_mat33_real xtalB_inv;
	bool inv_ok = false;

	if(use_B)
	{
		const t_mat33_real& xtalB = m_dyn->GetCrystalBTrafo();
		std::tie(xtalB_inv, inv_ok) = tl2::inv(xtalB);

		axis = xtalB * axis;
		vec = xtalB * vec;
	}

	t_mat33_real mat = tl2::rotation<t_mat33_real, t_vec3_real>(axis, angle, false);
	tl2::set_eps_0(mat, g_eps);


	// print the B and the rotation matrices
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	if(use_B)
		print_matrix(m_dyn->GetCrystalBTrafo(), "Crystal B Matrix", ostrResult);

	print_matrix(mat, "Transformation Matrix", ostrResult);

	using t_quat = boost::math::quaternion<t_real>;
	t_quat quat = tl2::rot3_to_quat<t_mat33_real, t_quat>(mat);
	tl2::set_eps_0(quat, g_eps);
	ostrResult << "<p>Trafo As Quaternion:<br>";
	ostrResult << quat;
	ostrResult << "</p>\n";

	if(use_B)
	{
		tl2::set_eps_0(axis, g_eps);
		ostrResult << "<p>Original Axis (lab): ";
		ostrResult << axis;
		ostrResult << "<br>\n";
		tl2::set_eps_0(vec, g_eps);
		ostrResult << "Original Vector (lab): ";
		ostrResult << vec;
		ostrResult << "</p>\n";
	}

	// print the rotated test vector
	t_vec3_real vec_rot = mat * vec;

	tl2::set_eps_0(vec_rot, g_eps);
	ostrResult << "<p>Rotated Vector (lab): ";
	ostrResult << vec_rot;

	if(use_B && inv_ok)
	{
		vec_rot = xtalB_inv * vec_rot;
		tl2::set_eps_0(vec_rot, g_eps);
		ostrResult << "<br>Rotated Vector (rlu): ";
		ostrResult << vec_rot;
	}
	ostrResult << "</p>\n";


	m_textRotation->setHtml(ostrResult.str().c_str());
}



void TrafoCalculator::CalculateProjection()
{
	using namespace tl2_ops;

	if(!m_textProjection)
		return;

	t_vec3_real axis = tl2::create<t_vec3_real>({
		(t_real)m_spinProjAxis[0]->value(),
		(t_real)m_spinProjAxis[1]->value(),
		(t_real)m_spinProjAxis[2]->value() });
	t_vec3_real vec = tl2::create<t_vec3_real>({
		(t_real)m_spinVecToProj[0]->value(),
		(t_real)m_spinVecToProj[1]->value(),
		(t_real)m_spinVecToProj[2]->value() });

	m_textProjection->clear();


	// apply crystal B matrix
	bool use_B = m_checkXtalProj->isChecked() && m_dyn;
	t_mat33_real xtalB_inv;
	bool inv_ok = false;

	if(use_B)
	{
		const t_mat33_real& xtalB = m_dyn->GetCrystalBTrafo();
		std::tie(xtalB_inv, inv_ok) = tl2::inv(xtalB);

		axis = xtalB * axis;
		vec = xtalB * vec;
	}

	t_mat33_real matProj = tl2::projector<t_mat33_real, t_vec3_real>(axis, false);
	t_mat33_real matOrthoProj = tl2::ortho_projector<t_mat33_real, t_vec3_real>(axis, false);
	tl2::set_eps_0(matProj, g_eps);
	tl2::set_eps_0(matOrthoProj, g_eps);


	// print the B and the rotation matrices
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	if(use_B)
		print_matrix(m_dyn->GetCrystalBTrafo(), "Crystal B Matrix", ostrResult);
	
	print_matrix(matProj, "Projection Matrix", ostrResult);
	print_matrix(matOrthoProj, "Orthogonal Projection Matrix", ostrResult);

	// print the original vectors
	tl2::set_eps_0(axis, g_eps);
	ostrResult << "<p>Original Axis (lab): ";
	ostrResult << axis;
	ostrResult << "<br>\n";
	tl2::set_eps_0(vec, g_eps);
	ostrResult << "Original Vector (lab): ";
	ostrResult << vec;
	ostrResult << "</p>\n";

	// print the projected test vector
	t_vec3_real vec_proj = matProj * vec;
	t_vec3_real vec_ortho_proj = matOrthoProj * vec;

	tl2::set_eps_0(vec_proj, g_eps);
	ostrResult << "<p>Projected Vector (lab): ";
	ostrResult << vec_proj;
	ostrResult << "<br>\n";
	tl2::set_eps_0(vec_proj, g_eps);
	ostrResult << "Orthogonally Projected Vector (lab): ";
	ostrResult << vec_ortho_proj;
	ostrResult << "</p>\n";

	if(use_B && inv_ok)
	{
		vec_proj = xtalB_inv * vec_proj;
		vec_ortho_proj = xtalB_inv * vec_ortho_proj;
		tl2::set_eps_0(vec_proj, g_eps);
		tl2::set_eps_0(vec_ortho_proj, g_eps);
		ostrResult << "<p>Projected Vector (rlu): ";
		ostrResult << vec_proj;
		ostrResult << "<br>\n";
		ostrResult << "Orthogonally Projected Vector (rlu): ";
		ostrResult << vec_ortho_proj;
		ostrResult << "</p>\n";
	}


	m_textProjection->setHtml(ostrResult.str().c_str());
}



void TrafoCalculator::CalculateCrossProduct()
{
	using namespace tl2_ops;

	if(!m_textCrossProd)
		return;

	t_vec3_real vec1 = tl2::create<t_vec3_real>({
		(t_real)m_spinVec1[0]->value(),
		(t_real)m_spinVec1[1]->value(),
		(t_real)m_spinVec1[2]->value() });
	t_vec3_real vec2 = tl2::create<t_vec3_real>({
		(t_real)m_spinVec2[0]->value(),
		(t_real)m_spinVec2[1]->value(),
		(t_real)m_spinVec2[2]->value() });

	m_textCrossProd->clear();


	// apply crystal B matrix
	bool use_B = m_checkXtalCrossProd->isChecked() && m_dyn;
	t_mat33_real xtalB_inv;
	bool inv_ok = false;

	if(use_B)
	{
		const t_mat33_real& xtalB = m_dyn->GetCrystalBTrafo();
		std::tie(xtalB_inv, inv_ok) = tl2::inv(xtalB);

		vec1 = xtalB * vec1;
		vec2 = xtalB * vec2;
	}


	// print the B and the rotation matrices
	std::ostringstream ostrResult;
	ostrResult.precision(g_prec);

	if(use_B)
		print_matrix(m_dyn->GetCrystalBTrafo(), "Crystal B Matrix", ostrResult);
	if(use_B && inv_ok)
		print_matrix(xtalB_inv, "Inverse B Matrix", ostrResult);

	tl2::set_eps_0(vec1, g_eps);
	ostrResult << "<p>Vector 1 (lab): ";
	ostrResult << vec1;
	ostrResult << "<br>\n";
	tl2::set_eps_0(vec2, g_eps);
	ostrResult << "Vector 2 (lab): ";
	ostrResult << vec2;
	ostrResult << "</p>\n";

	// print the cross product vector
	t_vec3_real vec_cross = tl2::cross<t_vec3_real>(vec1, vec2);
	t_vec3_real vec_cross_norm = vec_cross / tl2::norm<t_vec3_real>(vec_cross);

	tl2::set_eps_0(vec_cross, g_eps);
	ostrResult << "<p>Normal Vector (lab): ";
	ostrResult << vec_cross;
	ostrResult << "<br>\n";
	tl2::set_eps_0(vec_cross, g_eps);
	ostrResult << "Normalised Normal Vector (lab): ";
	ostrResult << vec_cross_norm;
	ostrResult << "</p>\n";

	if(use_B && inv_ok)
	{
		vec_cross = xtalB_inv * vec_cross;
		vec_cross_norm = xtalB_inv * vec_cross_norm;
		tl2::set_eps_0(vec_cross, g_eps);
		tl2::set_eps_0(vec_cross_norm, g_eps);

		ostrResult << "<p>Normal Vector (rlu): ";
		ostrResult << vec_cross;
		ostrResult << "<br>\n";
		ostrResult << "Normalised Normal Vector (rlu): ";
		ostrResult << vec_cross_norm;
		ostrResult << "</p>\n";
	}


	m_textCrossProd->setHtml(ostrResult.str().c_str());
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
