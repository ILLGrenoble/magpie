/**
 * 3d plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
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

#include <boost/scope_exit.hpp>

#include <cstdlib>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>

#include "plot3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"


#define SHOW_COORD_CUBE  false   // show coordinate cube by default?



/**
 * sets up the topology dialog
 */
Plot3DDlg::Plot3DDlg(QWidget *parent, std::shared_ptr<QSettings> sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("3D Plotter");
	setSizeGripEnabled(true);

	// create gl plotter
	m_dispplot = new tl2::GlPlot(this);
	m_dispplot->GetRenderer()->SetRestrictCamTheta(false);
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 50, 50, 50 }));
	m_dispplot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -50, -50, -50 }));
	m_dispplot->GetRenderer()->SetCoordMax(50.);
	m_dispplot->GetRenderer()->GetCamera().SetParallelRange(5.);
	m_dispplot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_fov));
	m_dispplot->GetRenderer()->GetCamera().SetDist(2.);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// surfaces table
	QWidget *surfs_panel = new QWidget(this);
	m_table_surfs = new QTableWidget(surfs_panel);
	m_table_surfs->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_surfs->setShowGrid(true);
	m_table_surfs->setSortingEnabled(false);
	m_table_surfs->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_surfs->setSelectionMode(QTableWidget::SingleSelection);
	m_table_surfs->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_surfs->verticalHeader()->setVisible(false);
	m_table_surfs->setColumnCount(NUM_COLS);
	m_table_surfs->setHorizontalHeaderItem(COL_SURF, new QTableWidgetItem{"Surf."});
	m_table_surfs->setHorizontalHeaderItem(COL_ACTIVE, new QTableWidgetItem{"Act."});
	m_table_surfs->setColumnWidth(COL_SURF, 40);
	m_table_surfs->setColumnWidth(COL_ACTIVE, 25);
	m_table_surfs->resizeColumnsToContents();

	// splitter for plot and surface list
	m_split_plot = new QSplitter(this);
	m_split_plot->setOrientation(Qt::Horizontal);
	m_split_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot->addWidget(m_dispplot);
	m_split_plot->addWidget(surfs_panel);
	m_split_plot->setCollapsible(0, false);
	m_split_plot->setCollapsible(1, true);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_dispplot), 24);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(surfs_panel), 1);

	// general plot context menu
	m_context = new QMenu(this);
	QAction *acCentre = new QAction("Centre Camera", m_context);
	QAction *acShowCoords = new QAction("Show Coordinate Planes", m_context);
	QAction *acSaveData = new QAction("Save Data...", m_context);
	QAction *acSaveScript = new QAction("Save Data As Script...", m_context);
	QAction *acSaveImage = new QAction("Save Image...", m_context);
	acShowCoords->setCheckable(true);
	acShowCoords->setChecked(SHOW_COORD_CUBE);
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	acSaveScript->setIcon(QIcon::fromTheme("text-x-script"));
	acSaveImage->setIcon(QIcon::fromTheme("image-x-generic"));
	m_context->addAction(acCentre);
	m_context->addSeparator();
	m_context->addAction(acShowCoords);
	m_context->addSeparator();
	m_context->addAction(acSaveData);
	m_context->addAction(acSaveScript);
	m_context->addSeparator();
	m_context->addAction(acSaveImage);

	// context menu for surfaces
	m_context_surf = new QMenu(this);
	QAction *acCentreOnObject = new QAction("Centre Camera on Surface", m_context_surf);
	m_context_surf->addAction(acCentre);
	m_context_surf->addAction(acCentreOnObject);
	m_context_surf->addSeparator();
	m_context_surf->addAction(acShowCoords);
	m_context_surf->addSeparator();
	m_context_surf->addAction(acSaveData);
	m_context_surf->addAction(acSaveScript);
	m_context_surf->addSeparator();
	m_context_surf->addAction(acSaveImage);

	// formulas
	QGroupBox *groupFormulas = new QGroupBox("Formulas", this);
	m_formulas = new QTextEdit(groupFormulas);
	m_formulas->setPlaceholderText("Enter formulas separated by ';'. Variables are 'x' and 'y'.");
	m_formulas->setReadOnly(false);
	//m_formulas->setFixedHeight(QFontMetrics{m_formulas->font()}.lineSpacing() * 4);
	m_formulas->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Ignored});

	// surfaces
	QGroupBox *groupQ = new QGroupBox("Surface Options", this);
	m_xrange[0] = new QDoubleSpinBox(groupQ);
	m_xrange[1] = new QDoubleSpinBox(groupQ);
	m_yrange[0] = new QDoubleSpinBox(groupQ);
	m_yrange[1] = new QDoubleSpinBox(groupQ);

	m_num_points[0] = new QSpinBox(groupQ);
	m_num_points[1] = new QSpinBox(groupQ);
	m_num_points[0]->setToolTip("Number of grid points along the first axis.");
	m_num_points[1]->setToolTip("Number of grid points along the second axis.");

	static const char* prefix[] = { "x = ", "y = " };
	static const char* prefix_minmax[] = { "min = ", "max = " };
	for(int i = 0; i < 2; ++i)
	{
		m_xrange[i]->setDecimals(4);
		m_xrange[i]->setMinimum(-99.9999);
		m_xrange[i]->setMaximum(+99.9999);
		m_xrange[i]->setSingleStep(0.01);
		m_xrange[i]->setValue(i == 1 ? 1. : 0.);
		m_xrange[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_xrange[i]->setPrefix(prefix_minmax[i]);
		m_xrange[i]->setToolTip("Range of the first axis.");

		m_yrange[i]->setDecimals(4);
		m_yrange[i]->setMinimum(-99.9999);
		m_yrange[i]->setMaximum(+99.9999);
		m_yrange[i]->setSingleStep(0.01);
		m_yrange[i]->setValue(i == 1 ? 1. : 0.);
		m_yrange[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_yrange[i]->setPrefix(prefix_minmax[i]);
		m_yrange[i]->setToolTip("Range of the second axis.");

		m_num_points[i]->setMinimum(1);
		m_num_points[i]->setMaximum(9999);
		m_num_points[i]->setSingleStep(1);
		m_num_points[i]->setPrefix(prefix[i]);
		m_num_points[i]->setValue(64);
	}

	// Q and E scale for plot
	QGroupBox *groupPlotOptions = new QGroupBox("Plot Options", this);
	m_x_scale = new QDoubleSpinBox(groupPlotOptions);
	m_y_scale = new QDoubleSpinBox(groupPlotOptions);
	m_z_scale = new QDoubleSpinBox(groupPlotOptions);

	m_x_scale->setPrefix("x = ");
	m_y_scale->setPrefix("y = ");
	m_z_scale->setPrefix("z = ");

	m_x_scale->setToolTip("Scaling factor along the first axis.");
	m_y_scale->setToolTip("Scaling factor along the second axis.");
	m_z_scale->setToolTip("Scaling factor along the third axis.");

	for(QDoubleSpinBox *box : { m_x_scale, m_y_scale, m_z_scale })
	{
		box->setDecimals(3);
		box->setMinimum(0.001);
		box->setMaximum(999.99);
		box->setSingleStep(0.1);
		box->setValue(1.);
		box->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	// camera view angle
	m_cam_phi = new QDoubleSpinBox(this);
	m_cam_phi->setRange(0., 360.);
	m_cam_phi->setSingleStep(1.);
	m_cam_phi->setDecimals(std::max(g_prec_gui - 2, 2));
	m_cam_phi->setPrefix("φ = ");
	m_cam_phi->setSuffix("°");
	m_cam_phi->setToolTip("Camera polar rotation angle φ.");

	m_cam_theta = new QDoubleSpinBox(this);
	m_cam_theta->setRange(-180., 180.);
	m_cam_theta->setSingleStep(1.);
	m_cam_theta->setDecimals(std::max(g_prec_gui - 2, 2));
	m_cam_theta->setPrefix("θ = ");
	m_cam_theta->setSuffix("°");
	m_cam_theta->setToolTip("Camera azimuthal rotation angle θ.");

	// camera perspective projection
	m_perspective = new QCheckBox("Perspective", this);
	m_perspective->setToolTip("Switch between perspective and parallel projection.");
	m_perspective->setChecked(true);

	// progress bar
	m_progress = new QProgressBar(this);
	m_progress->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btn_start_stop = new QPushButton("Calculate", this);

	// status bar
	m_status = new QLabel(this);
	m_status->setFrameShape(QFrame::Panel);
	m_status->setFrameShadow(QFrame::Sunken);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);

	// close button
	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	QPushButton *btnSaveData = btnbox->addButton("Save Data...", QDialogButtonBox::ActionRole);
	QPushButton *btnSaveScript = btnbox->addButton("Save Script...", QDialogButtonBox::ActionRole);
	btnSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	btnSaveScript->setIcon(QIcon::fromTheme("text-x-script"));
	btnbox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// surfaces panel grid
	int y = 0;
	QGridLayout *grid_surfs = new QGridLayout(surfs_panel);
	grid_surfs->setSpacing(4);
	grid_surfs->setContentsMargins(6, 6, 6, 6);
	grid_surfs->addWidget(m_table_surfs, y++, 0, 1, 1);

	// formulas grid
	y = 0;
	QGridLayout *gridFormulas = new QGridLayout(groupFormulas);
	gridFormulas->setSpacing(4);
	gridFormulas->setContentsMargins(6, 6, 6, 6);
	gridFormulas->addWidget(m_formulas, y++, 0, 1, 1);

	// surface options grid
	y = 0;
	QGridLayout *Qgrid = new QGridLayout(groupQ);
	Qgrid->setSpacing(4);
	Qgrid->setContentsMargins(6, 6, 6, 6);
	Qgrid->addWidget(new QLabel("x Range:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_xrange[0], y, 1, 1, 1);
	Qgrid->addWidget(m_xrange[1], y++, 2, 1, 1);
	Qgrid->addWidget(new QLabel("y Range:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_yrange[0], y, 1, 1, 1);
	Qgrid->addWidget(m_yrange[1], y++, 2, 1, 1);
	Qgrid->addWidget(new QLabel("Grid:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_num_points[0], y, 1, 1, 1);
	Qgrid->addWidget(m_num_points[1], y++, 2, 1, 1);

	// plot options grid
	y = 0;
	QGridLayout *plot_options_grid = new QGridLayout(groupPlotOptions);
	plot_options_grid->setSpacing(4);
	plot_options_grid->setContentsMargins(6, 6, 6, 6);
	plot_options_grid->addWidget(new QLabel("Scale:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_x_scale, y, 1, 1, 1);
	plot_options_grid->addWidget(m_y_scale, y, 2, 1, 1);
	plot_options_grid->addWidget(m_z_scale, y++, 3, 1, 1);
	plot_options_grid->addWidget(new QLabel("Camera Angles:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_cam_phi, y, 1, 1, 1);
	plot_options_grid->addWidget(m_cam_theta, y, 2, 1, 1);
	plot_options_grid->addWidget(m_perspective, y++, 3, 1, 1);

	// status grid
	QWidget *status_panel = new QWidget(this);
	QGridLayout *status_grid = new QGridLayout(status_panel);
	status_grid->setSpacing(0);
	status_grid->setContentsMargins(0, 0, 0, 0);
	status_grid->addWidget(m_status, y, 0, 1, 3);
	status_grid->addWidget(btnbox, y++, 3, 1, 1);

	// main grid
	y = 0;
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(m_split_plot, y++, 0, 1, 4);
	maingrid->addWidget(groupFormulas, y, 0, 1, 2);
	maingrid->addWidget(groupQ, y++, 2, 1, 2);
	maingrid->addWidget(groupPlotOptions, y++, 0, 1, 4);
	maingrid->addWidget(m_progress, y, 0, 1, 3);
	maingrid->addWidget(m_btn_start_stop, y++, 3, 1, 1);
	maingrid->addWidget(status_panel, y++, 0, 1, 4);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("plot3d/geo"))
			restoreGeometry(m_sett->value("plot3d/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("plot3d/splitter"))
			m_split_plot->restoreState(m_sett->value("plot3d/splitter").toByteArray());
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &Plot3DDlg::accept);
	connect(acCentreOnObject, &QAction::triggered, this, &Plot3DDlg::CentrePlotCameraOnObject);
	connect(acCentre, &QAction::triggered, this, &Plot3DDlg::CentrePlotCamera);
	connect(acShowCoords, &QAction::toggled, this, &Plot3DDlg::ShowPlotCoordCube);
	connect(acSaveData, &QAction::triggered, this, &Plot3DDlg::SaveData);
	connect(acSaveScript, &QAction::triggered, this, &Plot3DDlg::SaveScript);
	connect(acSaveImage, &QAction::triggered, this, &Plot3DDlg::SaveImage);
	connect(btnSaveData, &QAbstractButton::clicked, this, &Plot3DDlg::SaveData);
	connect(btnSaveScript, &QAbstractButton::clicked, this, &Plot3DDlg::SaveScript);

	connect(m_dispplot, &tl2::GlPlot::AfterGLInitialisation,
		this, &Plot3DDlg::AfterPlotGLInitialisation);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &Plot3DDlg::PlotPickerIntersection);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &Plot3DDlg::PlotCameraHasUpdated);
	connect(m_dispplot, &tl2::GlPlot::MouseClick, this, &Plot3DDlg::PlotMouseClick);
	connect(m_dispplot, &tl2::GlPlot::MouseDown, this, &Plot3DDlg::PlotMouseDown);
	connect(m_dispplot, &tl2::GlPlot::MouseUp, this, &Plot3DDlg::PlotMouseUp);
	connect(m_perspective, &QCheckBox::toggled, this, &Plot3DDlg::SetPlotPerspectiveProjection);

	for(QDoubleSpinBox *box : { m_x_scale, m_y_scale, m_z_scale })
	{
		connect(box, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]() { Plot(false); });
	}

	connect(m_cam_phi,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl phi) -> void
	{
		this->SetPlotCameraRotation(phi, m_cam_theta->value());
	});

	connect(m_cam_theta,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl theta) -> void
	{
		this->SetPlotCameraRotation(m_cam_phi->value(), theta);
	});

	// calculation
	connect(m_btn_start_stop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calc_enabled)
			Calculate();
		else
			m_stop_requested = true;
	});

	EnableCalculation(true);
}



Plot3DDlg::~Plot3DDlg()
{
}



void Plot3DDlg::ShowError(const QString& msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * toggle between "calculate" and "stop" button
 */
void Plot3DDlg::EnableCalculation(bool enable)
{
	m_calc_enabled = enable;

	if(enable)
	{
		m_btn_start_stop->setText("Calculate");
		m_btn_start_stop->setToolTip("Start surface calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btn_start_stop->setText("Stop");
		m_btn_start_stop->setToolTip("Stop running surface calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * the surface plot's camera properties have been updated
 */
void Plot3DDlg::PlotCameraHasUpdated()
{
	auto [phi, theta] = m_dispplot->GetRenderer()->GetCamera().GetRotation();

	phi = tl2::r2d<t_real>(phi);
	theta = tl2::r2d<t_real>(theta);

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
 * surface plot mouse button clicked
 */
void Plot3DDlg::PlotMouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(right)
	{
		const QPointF& _pt = m_dispplot->GetRenderer()->GetMousePosition();
		QPoint pt = m_dispplot->mapToGlobal(_pt.toPoint());

		if(m_cur_obj && m_surf_objs.find(*m_cur_obj) != m_surf_objs.end())
		{
			// surface selected
			m_context_surf->popup(pt);
		}
		else
		{
			// no surface selected
			m_context->popup(pt);
		}
	}
}



/**
 * surface plot mouse button pressed
 */
void Plot3DDlg::PlotMouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * surface plot mouse button released
 */
void Plot3DDlg::PlotMouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * surface plot has initialised
 */
void Plot3DDlg::AfterPlotGLInitialisation()
{
	if(!m_dispplot)
		return;

	// GL device infos
	auto [ver, shader_ver, vendor, renderer]
		= m_dispplot->GetRenderer()->GetGlDescr();

	emit GlDeviceInfos(ver, shader_ver, vendor, renderer);

	m_dispplot->GetRenderer()->SetCull(false);

	ShowPlotCoordCube(SHOW_COORD_CUBE);
	ShowPlotLabels(false);
	SetPlotPerspectiveProjection(m_perspective->isChecked());
	//SetPlotCoordinateSystem(m_coordsys->currentIndex());

	PlotCameraHasUpdated();
}



/**
 * show or hide the coordinate system
 */
void Plot3DDlg::ShowPlotCoordCube(bool show)
{
	// always hide coordinate cross
	if(auto obj = m_dispplot->GetRenderer()->GetCoordCross(); obj)
		m_dispplot->GetRenderer()->SetObjectVisible(*obj, false);

	// show coordinate cube
	for(auto obj : m_dispplot->GetRenderer()->GetCoordCube())
		m_dispplot->GetRenderer()->SetObjectVisible(obj, show);

	m_dispplot->update();
}



/**
 * show or hide the object labels
 */
void Plot3DDlg::ShowPlotLabels(bool show)
{
	m_dispplot->GetRenderer()->SetLabelsVisible(show);
	m_dispplot->update();
}



/**
 * choose between perspective or parallel projection
 */
void Plot3DDlg::SetPlotPerspectiveProjection(bool proj)
{
	m_dispplot->GetRenderer()->GetCamera().SetPerspectiveProjection(proj);
	m_dispplot->GetRenderer()->RequestViewportUpdate();
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * sets the camera's rotation angles
 */
void Plot3DDlg::SetPlotCameraRotation(t_real_gl phi, t_real_gl theta)
{
	phi = tl2::d2r<t_real>(phi);
	theta = tl2::d2r<t_real>(theta);

	m_dispplot->GetRenderer()->GetCamera().SetRotation(phi, theta);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	PlotCameraHasUpdated();
	m_dispplot->update();
}



/**
 * centre camera on currently selected object
 */
void Plot3DDlg::CentrePlotCameraOnObject()
{
	if(!m_cur_obj)
		return;

	t_mat_gl mat = m_dispplot->GetRenderer()->GetObjectMatrix(*m_cur_obj);

	// selected a surface?
	auto surf_iter = m_surf_objs.find(*m_cur_obj);
	if(surf_iter != m_surf_objs.end())
	{
		// translate camera to the surface's mean z position
		t_real z_mean = GetMeanZ(surf_iter->second) * m_z_scale->value();
		mat(2, 3) = z_mean;
	}

	m_dispplot->GetRenderer()->GetCamera().Centre(mat);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * centre camera on central position
 */
void Plot3DDlg::CentrePlotCamera()
{
	t_mat_gl matCentre = tl2::hom_translation<t_mat_gl>(
		m_cam_centre[0] * m_x_scale->value(),
		m_cam_centre[1] * m_y_scale->value(),
		m_cam_centre[2] * m_z_scale->value());

	m_dispplot->GetRenderer()->GetCamera().Centre(matCentre);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * switch between crystal and lab coordinates
 */
void Plot3DDlg::SetPlotCoordinateSystem(int which)
{
	m_dispplot->GetRenderer()->SetCoordSys(which);
}



/**
 * dialog is closing
 */
void Plot3DDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("plot3d/geo", saveGeometry());
		m_sett->setValue("plot3d/splitter", m_split_plot->saveState());
	}

	QDialog::accept();
}
