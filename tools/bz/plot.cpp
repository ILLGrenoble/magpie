/**
 * brillouin zone tool -- 3d plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include <boost/scope_exit.hpp>

#include "plot.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStyle>
#include <QtWidgets/QFileDialog>

#include <sstream>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"

using namespace tl2_ops;


BZPlotDlg::BZPlotDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Brillouin Zone");
	setFont(parent->font());
	setSizeGripEnabled(true);

	m_plot = std::make_shared<tl2::GlPlot>(this);
	m_plot->GetRenderer()->SetRestrictCamTheta(false);
	m_plot->GetRenderer()->SetCull(false);
	m_plot->GetRenderer()->SetBlend(true);
	m_plot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
	m_plot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
	m_plot->GetRenderer()->SetCoordMax(1.);
	m_plot->GetRenderer()->GetCamera().SetDist(2.);
	m_plot->GetRenderer()->GetCamera().SetParalellRange(4.);
	//m_plot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_structplot_fov));
	m_plot->GetRenderer()->GetCamera().UpdateTransformation();
	m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});


	//auto labCoordSys = new QLabel("Coordinate System:", this);
	//auto comboCoordSys = new QComboBox(this);
	//comboCoordSys->addItem("Fractional Units (rlu)");
	//comboCoordSys->addItem("Lab Units (\xe2\x84\xab)");

	m_show_coordcross_lab = new QCheckBox("Lab Basis", this);
	m_show_coordcross_xtal = new QCheckBox("Crystal Basis", this);
	//m_show_labels = new QCheckBox("Labels", this);
	m_show_plane = new QCheckBox("Scattering Plane", this);
	m_show_Qs = new QCheckBox("Q Vertices", this);
	m_show_coordcross_lab->setToolTip("Show or hide the basis vectors of the orthogonal lab system (in units of Å⁻¹).");
	m_show_coordcross_xtal->setToolTip("Show or hide the basis vectors of the generally non-orthogonal crystal system (in units of rlu, expressed in the lab basis).");
	//m_show_labels->setToolTip("Show or hide the object labels.");
	m_show_plane->setToolTip("Show or hide the plane used for the Brillouin zone intersection.");
	m_show_Qs->setToolTip("Show or hide the Brillouin zone vertices and Bragg peaks.");
	m_show_coordcross_lab->setChecked(true);
	m_show_coordcross_xtal->setChecked(false);
	//m_show_labels->setChecked(true);
	m_show_plane->setChecked(true);
	m_show_Qs->setChecked(true);

	m_perspective = new QCheckBox("Perspective Projection", this);
	m_perspective->setToolTip("Switch between perspective and parallel projection.");
	m_perspective->setChecked(true);

	QPushButton *btn_100 = new QPushButton("[100] View", this);
	QPushButton *btn_010 = new QPushButton("[010] View", this);
	QPushButton *btn_001 = new QPushButton("[001] View", this);
	QPushButton *btn_110 = new QPushButton("[110] View", this);
	btn_100->setToolTip("View along [100] axis.");
	btn_010->setToolTip("View along [010] axis.");
	btn_001->setToolTip("View along [001] axis.");
	btn_110->setToolTip("View along [110] axis.");

	m_cam_phi = new QDoubleSpinBox(this);
	m_cam_phi->setRange(0., 360.);
	m_cam_phi->setSingleStep(1.);
	m_cam_phi->setDecimals(std::max(m_prec_gui - 2, 2));
	m_cam_phi->setPrefix("φ = ");
	m_cam_phi->setSuffix("°");
	m_cam_phi->setToolTip("Camera polar rotation angle φ.");

	m_cam_theta = new QDoubleSpinBox(this);
	m_cam_theta->setRange(-180., 180.);
	m_cam_theta->setSingleStep(1.);
	m_cam_theta->setDecimals(std::max(m_prec_gui - 2, 2));
	m_cam_theta->setPrefix("θ = ");
	m_cam_theta->setSuffix("°");
	m_cam_theta->setToolTip("Camera azimuthal rotation angle θ.");

	// status bar
	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);

	// close button
	QPushButton *okbtn = new QPushButton("OK", this);
	okbtn->setIcon(style()->standardIcon(QStyle::SP_DialogOkButton));

	m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	//labCoordSys->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

	// plot context menu
	m_context = new QMenu(this);
	QAction *acSaveImage = new QAction("Save Image...", m_context);
	acSaveImage->setIcon(QIcon::fromTheme("image-x-generic"));
	m_context->addAction(acSaveImage);

	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	grid->addWidget(m_plot.get(), y++, 0, 1, 4);

	//grid->addWidget(labCoordSys, y, 0, 1, 1);
	//grid->addWidget(comboCoordSys, y, 1, 1, 1);
	grid->addWidget(m_show_coordcross_lab, y, 0, 1, 1);
	grid->addWidget(m_show_coordcross_xtal, y, 1, 1, 1);
	//grid->addWidget(m_show_labels, y, 1, 1, 1);
	grid->addWidget(m_show_plane, y, 2, 1, 1);
	grid->addWidget(m_show_Qs, y++, 3, 1, 1);

	grid->addWidget(btn_100, y, 0, 1, 1);
	grid->addWidget(btn_010, y, 1, 1, 1);
	grid->addWidget(btn_001, y, 2, 1, 1);
	grid->addWidget(btn_110, y++, 3, 1, 1);

	grid->addWidget(new QLabel("Camera Angles:", this), y, 0, 1, 1);
	grid->addWidget(m_cam_phi, y, 1, 1, 1);
	grid->addWidget(m_cam_theta, y, 2, 1, 1);
	grid->addWidget(m_perspective, y++, 3, 1, 1);

	grid->addWidget(m_status, y, 0, 1, 3);
	grid->addWidget(okbtn, y++, 3, 1, 1);

	// signals
	connect(m_plot.get(), &tl2::GlPlot::AfterGLInitialisation,
		this, &BZPlotDlg::AfterGLInitialisation);
	connect(m_plot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &BZPlotDlg::CameraHasUpdated);
	connect(m_plot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &BZPlotDlg::PickerIntersection);
	connect(m_plot.get(), &tl2::GlPlot::MouseClick, this, &BZPlotDlg::PlotMouseClick);
	connect(m_plot.get(), &tl2::GlPlot::MouseDown, this, &BZPlotDlg::PlotMouseDown);
	connect(m_plot.get(), &tl2::GlPlot::MouseUp, this, &BZPlotDlg::PlotMouseUp);
	connect(m_show_coordcross_lab, &QCheckBox::toggled, this, &BZPlotDlg::ShowCoordCrossLab);
	connect(m_show_coordcross_xtal, &QCheckBox::toggled, this, &BZPlotDlg::ShowCoordCrossXtal);
	//connect(m_show_labels, &QCheckBox::toggled, this, &BZPlotDlg::ShowLabels);
	connect(m_show_plane, &QCheckBox::toggled, this, &BZPlotDlg::ShowPlane);
	connect(m_show_Qs, &QCheckBox::toggled, this, &BZPlotDlg::ShowQVertices);
	/*connect(comboCoordSys, static_cast<void (QComboBox::*)(int)>(
		&QComboBox::currentIndexChanged),
		this, [this](int val)
	{
		if(this->m_plot)
			this->m_plot->GetRenderer()->SetCoordSys(val);
	});*/
	connect(acSaveImage, &QAction::triggered, this, &BZPlotDlg::SaveImage);
	connect(okbtn, &QAbstractButton::clicked, this, &QDialog::accept);

	connect(btn_100, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(t_real_gl(90.), -t_real_gl(90.));
	});
	connect(btn_010, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(0., -t_real_gl(90.));
	});
	connect(btn_001, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(0., t_real_gl(180.));
	});
	connect(btn_110, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(t_real_gl(45.), -t_real_gl(90.));
	});

	connect(m_perspective, &QCheckBox::toggled, this, &BZPlotDlg::SetPerspectiveProjection);

	connect(m_cam_phi,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl phi) -> void
	{
		this->SetCameraRotation(phi, m_cam_theta->value());
	});

	connect(m_cam_theta,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl theta) -> void
	{
		this->SetCameraRotation(m_cam_phi->value(), theta);
	});

	if(m_sett && m_sett->contains("bz3d/geo"))
		restoreGeometry(m_sett->value("bz3d/geo").toByteArray());
	else
		resize(640, 640);
}


/**
 * dialog is closing
 */
void BZPlotDlg::accept()
{
	if(!m_sett)
		return;

	m_sett->setValue("bz3d/geo", saveGeometry());
	QDialog::accept();
}


/**
 * show or hide the coordinate cross
 */
void BZPlotDlg::ShowCoordCrossLab(bool show)
{
	if(!m_plot)
		return;

	if(auto obj = m_plot->GetRenderer()->GetCoordCross(false); obj)
	{
		m_plot->GetRenderer()->SetObjectVisible(*obj, show);
		m_plot->update();
	}
}


/**
 * show or hide the coordinate cross
 */
void BZPlotDlg::ShowCoordCrossXtal(bool show)
{
	if(!m_plot)
		return;

	if(auto obj = m_plot->GetRenderer()->GetCoordCross(true); obj)
	{
		m_plot->GetRenderer()->SetObjectVisible(*obj, show);
		m_plot->update();
	}
}


/**
 * show or hide the object labels
 */
void BZPlotDlg::ShowLabels(bool show)
{
	if(!m_plot)
		return;

	m_plot->GetRenderer()->SetLabelsVisible(show);
	m_plot->update();
}


/**
 * show or hide the BZ cut plane
 */
void BZPlotDlg::ShowPlane(bool show)
{
	if(!m_plot)
		return;

	m_plot->GetRenderer()->SetObjectVisible(m_plane, show);
	m_plot->update();
}


/**
 * show or hide the Bragg peaks and Voronoi vertices
 */
void BZPlotDlg::ShowQVertices(bool show)
{
	if(!m_plot)
		return;

	for(std::size_t obj : m_objsVoronoi)
		m_plot->GetRenderer()->SetObjectVisible(obj, show);
	for(std::size_t obj : m_objsBragg)
		m_plot->GetRenderer()->SetObjectVisible(obj, show);
	for(std::size_t obj : m_objsLines)
		m_plot->GetRenderer()->SetObjectVisible(obj, show);

	m_plot->update();
}


/**
 * set the crystal matrices
 */
void BZPlotDlg::SetABTrafo(const t_mat_bz& crystA, const t_mat_bz& crystB)
{
	m_crystA = crystA;
	m_crystB = crystB;

	if(!m_plot)
		return;

	t_mat_gl matA{crystA};
	m_plot->GetRenderer()->SetBTrafo(crystB, &matA, false);
}


/**
 * add a vertex to the plot
 */
void BZPlotDlg::AddVertex(const t_vec_bz& pos,
	std::vector<std::size_t>* cont, t_real_gl scale,
	t_real_gl r, t_real_gl g, t_real_gl b)
{
	if(!m_plot)
		return;

	t_real_gl posx = static_cast<t_real_gl>(pos[0]);
	t_real_gl posy = static_cast<t_real_gl>(pos[1]);
	t_real_gl posz = static_cast<t_real_gl>(pos[2]);

	auto obj = m_plot->GetRenderer()->AddLinkedObject(m_sphere, 0,0,0, r,g,b,1);
	//auto obj = m_plot->GetRenderer()->AddSphere(0.05, 0,0,0, r,g,b,1);
	m_plot->GetRenderer()->SetObjectMatrix(obj,
		tl2::hom_translation<t_mat_gl>(posx, posy, posz) *
		tl2::hom_scaling<t_mat_gl>(scale, scale, scale));
	m_plot->GetRenderer()->SetObjectIntersectable(obj, false);

	if(cont)
		cont->push_back(obj);
	m_plot->update();	
}



/**
 * add a voronoi vertex to the plot
 */
void BZPlotDlg::AddVoronoiVertex(const t_vec_bz& pos)
{
	AddVertex(pos, &m_objsVoronoi, 1., 0., 0., 1.);
}


/**
 * add a bragg peak to the plot
 */
void BZPlotDlg::AddBraggPeak(const t_vec_bz& pos)
{
	AddVertex(pos, &m_objsBragg, 1., 1., 0., 0.);
}


/**
 * add a line from a start to an end point
 */
void BZPlotDlg::AddLine(const t_vec_bz& start, const t_vec_bz& end, bool add_vertices)
{
	if(add_vertices)
	{
		AddVertex(start, &m_objsLines, 0.5,  1., 0., 0.);
		AddVertex(end, &m_objsLines, 0.5,  0., 0., 1.);
	}

	t_vec3_gl dir = tl2::convert<t_vec3_gl>(end - start);
	t_vec3_gl mid = tl2::convert<t_vec3_gl>(start + (end - start)*0.5);
	t_real_gl len = tl2::norm(dir);

	auto arrow = m_plot->GetRenderer()->AddArrow(
		0.025, len,  0., 0., 0.5,  1., 1., 1., 1.);
	m_plot->GetRenderer()->SetObjectMatrix(arrow,
		tl2::get_arrow_matrix<t_vec3_gl, t_mat_gl, t_real_gl>(
			dir,                                 // to
			1.,                                  // post-scale
			tl2::create<t_vec3_gl>({ 0, 0, 0 }), // post-translate
			tl2::create<t_vec3_gl>({ 0, 0, 1 }), // from
			1.,                                  // pre-scale
			tl2::convert<t_vec3_gl>(mid)));    // pre-translate

	m_plot->GetRenderer()->SetObjectIntersectable(arrow, false);
	m_objsLines.push_back(arrow);
}


/**
 * add polygons to the plot
 */
void BZPlotDlg::AddTriangles(const std::vector<t_vec_bz>& _vecs,
	const std::vector<std::size_t> *faceindices)
{
	if(!m_plot || _vecs.size() < 3)
		return;
	if(!m_plot->GetRenderer())
		return;

	t_real_gl r = 1, g = 0, b = 0;
	std::vector<t_vec3_gl> vecs, norms;
	vecs.reserve(_vecs.size());
	norms.reserve(_vecs.size());

	for(std::size_t idx = 0; idx < _vecs.size() - 2; idx += 3)
	{
		t_vec3_gl vec1 = tl2::convert<t_vec3_gl>(_vecs[idx]);
		t_vec3_gl vec2 = tl2::convert<t_vec3_gl>(_vecs[idx+1]);
		t_vec3_gl vec3 = tl2::convert<t_vec3_gl>(_vecs[idx+2]);

		t_vec3_gl norm = tl2::cross<t_vec3_gl>(vec2 - vec1, vec3 - vec1);
		norm /= tl2::norm(norm);

		t_vec3_gl mid = (vec1+vec2+vec3)/3.;
		mid /= tl2::norm<t_vec3_gl>(mid);

		// change sign of norm / sense of veritices?
		if(tl2::inner<t_vec3_gl>(norm, mid) < 0.)
		{
			vecs.emplace_back(std::move(vec3));
			vecs.emplace_back(std::move(vec2));
			vecs.emplace_back(std::move(vec1));
			norms.push_back(-norm);
		}
		else
		{
			vecs.emplace_back(std::move(vec1));
			vecs.emplace_back(std::move(vec2));
			vecs.emplace_back(std::move(vec3));
			norms.push_back(norm);
		}
	}

	auto obj = m_plot->GetRenderer()->AddTriangleObject(vecs, norms, r, g, b, 1);
	m_objsBZ.push_back(obj);

	// set (or clear) the face indices of this object's triangles
	if(faceindices)
	{
		m_objFaceIndices.insert(std::make_pair(obj, *faceindices));
	}
	else
	{
		if(auto iter = m_objFaceIndices.find(obj); iter != m_objFaceIndices.end())
			m_objFaceIndices.erase(iter);
	}

	m_plot->update();
}


/**
 * set the brillouin zone cut plane
 */
void BZPlotDlg::SetPlane(const t_vec_bz& _norm, t_real d)
{
	if(!m_plot)
		return;

	t_vec3_gl norm = tl2::convert<t_vec3_gl>(_norm);
	t_vec3_gl norm_old = tl2::create<t_vec3_gl>({ 0, 0, 1 });
	t_vec3_gl rot_vec = tl2::create<t_vec3_gl>({ 1, 0, 0 });

	t_vec3_gl offs = d * norm;
	t_mat_gl rot = tl2::hom_rotation<t_mat_gl, t_vec3_gl>(norm_old, norm, &rot_vec);
	t_mat_gl trans = tl2::hom_translation<t_mat_gl>(offs[0], offs[1], offs[2]);

	m_plot->GetRenderer()->SetObjectMatrix(m_plane, trans*rot);
	m_plot->update();
}


void BZPlotDlg::ClearLines(bool update)
{
	if(!m_plot)
		return;

	for(std::size_t obj : m_objsLines)
		m_plot->GetRenderer()->RemoveObject(obj);

	m_objsLines.clear();

	if(update)
		m_plot->update();
}


void BZPlotDlg::Clear()
{
	if(!m_plot)
		return;

	for(std::size_t obj : m_objsBZ)
		m_plot->GetRenderer()->RemoveObject(obj);
	for(std::size_t obj : m_objsBragg)
		m_plot->GetRenderer()->RemoveObject(obj);
	for(std::size_t obj : m_objsVoronoi)
		m_plot->GetRenderer()->RemoveObject(obj);

	m_objsBZ.clear();
	m_objsBragg.clear();
	m_objsVoronoi.clear();
	m_objFaceIndices.clear();

	ClearLines(false);

	m_plot->update();
}


/**
 * mouse hovers over 3d object
 */
void BZPlotDlg::PickerIntersection(
	const t_vec3_gl *pos, std::size_t objIdx, std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl *posSphere)
{
	if(pos)
		m_curPickedObj = long(objIdx);
	else
		m_curPickedObj = -1;

	if(m_curPickedObj <= 0 || !pos)
	{
		SetStatusMsg("");
		return;
	}

	t_vec_bz QinvA = tl2::convert<t_vec_bz>(*pos);
	t_mat_bz Binv =  tl2::trans(m_crystA) / (t_real(2)*tl2::pi<t_real>);
	m_cur_Qrlu = Binv * QinvA;

	tl2::set_eps_0<t_vec_bz>(QinvA, m_eps);
	tl2::set_eps_0<t_vec_bz>(m_cur_Qrlu, m_eps);

	std::ostringstream ostr;
	ostr.precision(m_prec_gui);

	// see if this object has stored face indices
	if(auto iter = m_objFaceIndices.find(objIdx); iter != m_objFaceIndices.end())
	{
		if(triagIdx < iter->second.size())
			ostr << "Face " << iter->second[triagIdx] << ": ";
	}
	else if(objIdx == m_plane)
	{
		ostr << "Plane: ";
	}

	ostr << "Q = (" << QinvA[0] << ", " << QinvA[1] << ", " << QinvA[2] << ") Å⁻¹";
	ostr << " = (" << m_cur_Qrlu[0] << ", " << m_cur_Qrlu[1] << ", " << m_cur_Qrlu[2] << ") rlu.";

	SetStatusMsg(ostr.str());
}


/**
 * set status label text in 3d dialog
 */
void BZPlotDlg::SetStatusMsg(const std::string& msg)
{
	m_status->setText(msg.c_str());
}


/**
 * mouse button clicked
 */
void BZPlotDlg::PlotMouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(left)
	{
		m_clicked_Q_rlu[0] = m_cur_Qrlu;
	}

	if(right)
	{
		m_clicked_Q_rlu[1] = m_cur_Qrlu;

		const QPointF& _pt = m_plot->GetRenderer()->GetMousePosition();
		QPoint pt = m_plot->mapToGlobal(_pt.toPoint());

		m_context->popup(pt);
	}
}


/**
 * mouse button pressed
 */
void BZPlotDlg::PlotMouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}


/**
 * mouse button released
 */
void BZPlotDlg::PlotMouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}


void BZPlotDlg::AfterGLInitialisation()
{
	if(!m_plot)
		return;

	// reference sphere and plane for linked objects
	m_sphere = m_plot->GetRenderer()->AddSphere(0.05, 0.,0.,0., 1.,1.,1.,1.);
	m_plane = m_plot->GetRenderer()->AddPlane(0.,0.,1., 0.,0.,0., 5., 0.75,0.75,0.75,0.5);
	m_plot->GetRenderer()->SetObjectVisible(m_sphere, false);
	m_plot->GetRenderer()->SetObjectIntersectable(m_sphere, false);
	m_plot->GetRenderer()->SetObjectVisible(m_plane, true);
	m_plot->GetRenderer()->SetObjectPriority(m_plane, 0);

	emit NeedRecalc();

	// GL device info
	auto [ver, shader_ver, vendor, renderer] = m_plot->GetRenderer()->GetGlDescr();
	emit GlDeviceInfos(ver, shader_ver, vendor, renderer);

	ShowCoordCrossLab(m_show_coordcross_lab->isChecked());
	ShowCoordCrossXtal(m_show_coordcross_xtal->isChecked());
	//ShowLabels(m_labels->isChecked());
	SetPerspectiveProjection(m_perspective->isChecked());
	CameraHasUpdated();
}


/**
 * choose between perspective or parallel projection
 */
void BZPlotDlg::SetPerspectiveProjection(bool proj)
{
	m_plot->GetRenderer()->GetCamera().SetPerspectiveProjection(proj);
	m_plot->GetRenderer()->RequestViewportUpdate();
	m_plot->GetRenderer()->GetCamera().UpdateTransformation();
	m_plot->update();
}


/**
 * sets the camera's rotation angles
 */
void BZPlotDlg::SetCameraRotation(t_real_gl phi, t_real_gl theta)
{
	phi = tl2::d2r<t_real_gl>(phi);
	theta = tl2::d2r<t_real_gl>(theta);

	m_plot->GetRenderer()->GetCamera().SetRotation(phi, theta);
	m_plot->GetRenderer()->GetCamera().UpdateTransformation();
	CameraHasUpdated();
	m_plot->update();
}


/**
 * the camera's properties have been updated
 */
void BZPlotDlg::CameraHasUpdated()
{
	auto [phi, theta] = m_plot->GetRenderer()->GetCamera().GetRotation();

	phi = tl2::r2d<t_real_gl>(phi);
	theta = tl2::r2d<t_real_gl>(theta);

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_cam_phi->blockSignals(false);
		this_->m_cam_theta->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_cam_phi->blockSignals(true);
	m_cam_theta->blockSignals(true);

	m_cam_phi->setValue(phi);
	m_cam_theta->setValue(theta);
}


/**
 * save an image of the brillouin zone
 */
void BZPlotDlg::SaveImage()
{
	if(!m_plot)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("bz3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Brillouin Zone Image",
		dirLast, "PNG Files (*.png)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("bz3d/dir", QFileInfo(filename).path());

	m_plot->grabFramebuffer().save(filename, nullptr, 90);
}


void BZPlotDlg::SetEps(t_real eps)
{
	m_eps = eps;
}


void BZPlotDlg::SetPrecGui(int prec)
{
	m_prec_gui = prec;
}


QMenu* BZPlotDlg::GetContextMenu()
{
	return m_context;
}


const t_vec_bz& BZPlotDlg::GetClickedPosition(bool right_button) const
{
	return right_button ? m_clicked_Q_rlu[1] : m_clicked_Q_rlu[0];
}
