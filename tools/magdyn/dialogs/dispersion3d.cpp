/**
 * magnetic dynamics -- 3d dispersion plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "dispersion3d.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"


#define SHOW_COORD_CUBE  false   // show coordinate cube by default?



/**
 * sets up the topology dialog
 */
Dispersion3DDlg::Dispersion3DDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("3D Dispersion");
	setSizeGripEnabled(true);

	// create gl plotter
	m_dispplot = new tl2::GlPlot(this);
	m_dispplot->GetRenderer()->SetRestrictCamTheta(false);
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 50, 50, 50 }));
	m_dispplot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -50, -50, -50 }));
	m_dispplot->GetRenderer()->SetCoordMax(50.);
	m_dispplot->GetRenderer()->GetCamera().SetParalellRange(100.);
	m_dispplot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_structplot_fov));
	m_dispplot->GetRenderer()->GetCamera().SetDist(40.);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// magnon band table
	QWidget *bands_panel = new QWidget(this);
	m_table_bands = new QTableWidget(bands_panel);
	m_table_bands->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_bands->setShowGrid(true);
	m_table_bands->setSortingEnabled(false);
	m_table_bands->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_bands->setSelectionMode(QTableWidget::SingleSelection);
	m_table_bands->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_bands->verticalHeader()->setVisible(false);
	m_table_bands->setColumnCount(NUM_COLS_BC);
	m_table_bands->setHorizontalHeaderItem(COL_BC_BAND, new QTableWidgetItem{"Band"});
	m_table_bands->setHorizontalHeaderItem(COL_BC_ACTIVE, new QTableWidgetItem{"Act."});
	m_table_bands->setColumnWidth(COL_BC_BAND, 40);
	m_table_bands->setColumnWidth(COL_BC_ACTIVE, 25);
	m_table_bands->resizeColumnsToContents();

	m_only_pos_E = new QCheckBox("E ≥ 0", bands_panel);
	m_only_pos_E->setChecked(true);
	m_only_pos_E->setToolTip("Ignore magnon annihilation.");

	// splitter for plot and magnon band list
	m_split_plot = new QSplitter(this);
	m_split_plot->setOrientation(Qt::Horizontal);
	m_split_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot->addWidget(m_dispplot);
	m_split_plot->addWidget(bands_panel);
	m_split_plot->setCollapsible(0, false);
	m_split_plot->setCollapsible(1, true);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_dispplot), 24);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(bands_panel), 1);

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

	// context menu for sites
	m_context_band = new QMenu(this);
	QAction *acCentreOnObject = new QAction("Centre Camera on Band", m_context_band);
	m_context_band->addAction(acCentre);
	m_context_band->addAction(acCentreOnObject);
	m_context_band->addSeparator();
	m_context_band->addAction(acShowCoords);
	m_context_band->addSeparator();
	m_context_band->addAction(acSaveData);
	m_context_band->addAction(acSaveScript);
	m_context_band->addSeparator();
	m_context_band->addAction(acSaveImage);

	// Q coordinates
	QGroupBox *groupQ = new QGroupBox("Dispersion", this);
	m_Q_origin[0] = new QDoubleSpinBox(groupQ);
	m_Q_origin[1] = new QDoubleSpinBox(groupQ);
	m_Q_origin[2] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[0] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[1] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[2] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[0] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[1] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[2] = new QDoubleSpinBox(groupQ);
	m_num_Q_points[0] = new QSpinBox(groupQ);
	m_num_Q_points[1] = new QSpinBox(groupQ);
	m_num_Q_points[0]->setToolTip("Number of grid points along the first momentum axis.");
	m_num_Q_points[1]->setToolTip("Number of grid points along the second momentum axis.");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin[i]->setDecimals(4);
		m_Q_origin[i]->setMinimum(-99.9999);
		m_Q_origin[i]->setMaximum(+99.9999);
		m_Q_origin[i]->setSingleStep(0.01);
		m_Q_origin[i]->setValue(0.);
		//m_Q_origin[i]->setSuffix(" rlu");
		m_Q_origin[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_origin[i]->setPrefix(hklPrefix[i]);
		m_Q_origin[i]->setToolTip("Starting momentum transfer.");

		m_Q_dir1[i]->setDecimals(4);
		m_Q_dir1[i]->setMinimum(-99.9999);
		m_Q_dir1[i]->setMaximum(+99.9999);
		m_Q_dir1[i]->setSingleStep(0.01);
		m_Q_dir1[i]->setValue(i == 0 ? 1. : 0.);
		//m_Q_dir1[i]->setSuffix(" rlu");
		m_Q_dir1[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir1[i]->setPrefix(hklPrefix[i]);
		m_Q_dir1[i]->setToolTip("Direction of momentum transfer along the first axis.");

		m_Q_dir2[i]->setDecimals(4);
		m_Q_dir2[i]->setMinimum(-99.9999);
		m_Q_dir2[i]->setMaximum(+99.9999);
		m_Q_dir2[i]->setSingleStep(0.01);
		m_Q_dir2[i]->setValue(i == 1 ? 1. : 0.);
		//m_Q_dir2[i]->setSuffix(" rlu");
		m_Q_dir2[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir2[i]->setPrefix(hklPrefix[i]);
		m_Q_dir2[i]->setToolTip("Direction of momentum transfer along the second axis.");
	}

	for(int i = 0; i < 2; ++i)
	{
		m_num_Q_points[i]->setMinimum(1);
		m_num_Q_points[i]->setMaximum(9999);
		m_num_Q_points[i]->setSingleStep(1);
		m_num_Q_points[i]->setValue(64);
	}

	// main dispersion button
	QPushButton *btnMainQ = new QPushButton("From Main Q", groupQ);
	btnMainQ->setToolTip("Set the Q origin and directions from the dispersion in the main window.");

	// minimum cutoff for filtering S(Q, E)
	m_S_filter_enable = new QCheckBox("Minimum S(Q, E):", groupQ);
	m_S_filter_enable->setChecked(false);
	m_S_filter_enable->setToolTip("Enable minimum S(Q, E).");

	m_S_filter = new QDoubleSpinBox(groupQ);
	m_S_filter->setDecimals(5);
	m_S_filter->setMinimum(0.);
	m_S_filter->setMaximum(9999.99999);
	m_S_filter->setSingleStep(0.01);
	m_S_filter->setValue(0.01);
	m_S_filter->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_S_filter->setToolTip("Minimum S(Q, E) to keep.");

	// unite degenerate energies
	m_unite_degeneracies = new QCheckBox("Unite degenerate Es.", groupQ);
	m_unite_degeneracies->setChecked(false);
	m_unite_degeneracies->setToolTip("Unite degenerate energies.");

	// Q and E scale for plot
	QGroupBox *groupPlotOptions = new QGroupBox("Plot Options", this);
	m_Q_scale1 = new QDoubleSpinBox(groupPlotOptions);
	m_Q_scale2 = new QDoubleSpinBox(groupPlotOptions);
	m_E_scale = new QDoubleSpinBox(groupPlotOptions);

	m_Q_scale1->setToolTip("Scaling factor along the first momentum axis.");
	m_Q_scale2->setToolTip("Scaling factor along the second momentum axis.");
	m_E_scale->setToolTip("Scaling factor along the energy axis.");

	for(QDoubleSpinBox *box : { m_Q_scale1, m_Q_scale2, m_E_scale })
	{
		box->setDecimals(3);
		box->setMinimum(0.001);
		box->setMaximum(999.99);
		box->setSingleStep(box == m_E_scale ? 0.1 : 0.5);
		box->setValue(box == m_E_scale ? 4. : 32.);
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
	m_perspective = new QCheckBox("Perspective Projection", this);
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

	// bands panel grid
	int y = 0;
	QGridLayout *grid_bands = new QGridLayout(bands_panel);
	grid_bands->setSpacing(4);
	grid_bands->setContentsMargins(6, 6, 6, 6);
	grid_bands->addWidget(m_table_bands, y++, 0, 1, 1);
	grid_bands->addWidget(m_only_pos_E, y++, 0, 1, 1);

	// Q coordinates grid
	y = 0;
	QGridLayout *Qgrid = new QGridLayout(groupQ);
	Qgrid->setSpacing(4);
	Qgrid->setContentsMargins(6, 6, 6, 6);
	Qgrid->addWidget(new QLabel("Q Origin:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_origin[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_origin[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_origin[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Direction 1:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir1[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir1[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir1[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Direction 2:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir2[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir2[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir2[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Grid Points:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_num_Q_points[0], y, 1, 1, 1);
	Qgrid->addWidget(m_num_Q_points[1], y, 2, 1, 1);
	Qgrid->addWidget(btnMainQ, y++, 3, 1, 1);
	Qgrid->addWidget(m_S_filter_enable, y, 0, 1, 1);
	Qgrid->addWidget(m_S_filter, y, 1, 1, 1);
	Qgrid->addWidget(m_unite_degeneracies, y++, 2, 1, 1);

	// plot options grid
	y = 0;
	QGridLayout *plot_options_grid = new QGridLayout(groupPlotOptions);
	plot_options_grid->setSpacing(4);
	plot_options_grid->setContentsMargins(6, 6, 6, 6);
	plot_options_grid->addWidget(new QLabel("Q Scale:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_Q_scale1, y, 1, 1, 1);
	plot_options_grid->addWidget(m_Q_scale2, y, 2, 1, 1);
	plot_options_grid->addWidget(new QLabel("E Scale:", this), y, 3, 1, 1);
	plot_options_grid->addWidget(m_E_scale, y++, 4, 1, 1);
	plot_options_grid->addWidget(new QLabel("Camera Angles:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_cam_phi, y, 1, 1, 1);
	plot_options_grid->addWidget(m_cam_theta, y, 2, 1, 1);
	plot_options_grid->addWidget(m_perspective, y++, 3, 1, 2);

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
	maingrid->addWidget(groupQ, y++, 0, 1, 4);
	maingrid->addWidget(groupPlotOptions, y++, 0, 1, 4);
	maingrid->addWidget(m_progress, y, 0, 1, 3);
	maingrid->addWidget(m_btn_start_stop, y++, 3, 1, 1);
	maingrid->addWidget(status_panel, y++, 0, 1, 4);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("dispersion3d/geo"))
			restoreGeometry(m_sett->value("dispersion3d/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("dispersion3d/splitter"))
			m_split_plot->restoreState(m_sett->value("dispersion3d/splitter").toByteArray());
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &Dispersion3DDlg::accept);
	connect(acCentreOnObject, &QAction::triggered, this, &Dispersion3DDlg::CentrePlotCameraOnObject);
	connect(acCentre, &QAction::triggered, this, &Dispersion3DDlg::CentrePlotCamera);
	connect(acShowCoords, &QAction::toggled, this, &Dispersion3DDlg::ShowPlotCoordCube);
	connect(acSaveData, &QAction::triggered, this, &Dispersion3DDlg::SaveData);
	connect(acSaveScript, &QAction::triggered, this, &Dispersion3DDlg::SaveScript);
	connect(acSaveImage, &QAction::triggered, this, &Dispersion3DDlg::SaveImage);
	connect(btnSaveData, &QAbstractButton::clicked, this, &Dispersion3DDlg::SaveData);
	connect(btnSaveScript, &QAbstractButton::clicked, this, &Dispersion3DDlg::SaveScript);

	connect(m_dispplot, &tl2::GlPlot::AfterGLInitialisation,
		this, &Dispersion3DDlg::AfterPlotGLInitialisation);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &Dispersion3DDlg::PlotPickerIntersection);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &Dispersion3DDlg::PlotCameraHasUpdated);
	connect(m_dispplot, &tl2::GlPlot::MouseClick, this, &Dispersion3DDlg::PlotMouseClick);
	connect(m_dispplot, &tl2::GlPlot::MouseDown, this, &Dispersion3DDlg::PlotMouseDown);
	connect(m_dispplot, &tl2::GlPlot::MouseUp, this, &Dispersion3DDlg::PlotMouseUp);
	connect(m_perspective, &QCheckBox::toggled, this, &Dispersion3DDlg::SetPlotPerspectiveProjection);
	connect(m_only_pos_E, &QCheckBox::toggled, [this]() { Plot(true); });
	connect(btnMainQ, &QAbstractButton::clicked, this, &Dispersion3DDlg::FromMainQ);
	connect(m_S_filter_enable, &QCheckBox::toggled, m_S_filter, &QDoubleSpinBox::setEnabled);

	for(QDoubleSpinBox *box : { m_Q_scale1, m_Q_scale2, m_E_scale })
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

	m_S_filter->setEnabled(m_S_filter_enable->isChecked());
	EnableCalculation(true);
}



Dispersion3DDlg::~Dispersion3DDlg()
{
}



void Dispersion3DDlg::ShowError(const QString& msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * toggle between "calculate" and "stop" button
 */
void Dispersion3DDlg::EnableCalculation(bool enable)
{
	m_calc_enabled = enable;

	if(enable)
	{
		m_btn_start_stop->setText("Calculate");
		m_btn_start_stop->setToolTip("Start dispersion calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btn_start_stop->setText("Stop");
		m_btn_start_stop->setToolTip("Stop running dispersion calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * the dispersion plot's camera properties have been updated
 */
void Dispersion3DDlg::PlotCameraHasUpdated()
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
 * dispersion plot mouse button clicked
 */
void Dispersion3DDlg::PlotMouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(right)
	{
		const QPointF& _pt = m_dispplot->GetRenderer()->GetMousePosition();
		QPoint pt = m_dispplot->mapToGlobal(_pt.toPoint());

		if(m_cur_obj && m_band_objs.find(*m_cur_obj) != m_band_objs.end())
		{
			// band selected
			m_context_band->popup(pt);
		}
		else
		{
			// no band selected
			m_context->popup(pt);
		}
	}
}



/**
 * dispersion plot mouse button pressed
 */
void Dispersion3DDlg::PlotMouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * dispersion plot mouse button released
 */
void Dispersion3DDlg::PlotMouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * dispersion plot has initialised
 */
void Dispersion3DDlg::AfterPlotGLInitialisation()
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
void Dispersion3DDlg::ShowPlotCoordCube(bool show)
{
	// always hide coordinate cross
	if(auto obj = m_dispplot->GetRenderer()->GetCoordCross(); obj)
		m_dispplot->GetRenderer()->SetObjectVisible(*obj, false);

	// show coordinate cube
	if(auto obj = m_dispplot->GetRenderer()->GetCoordCube(); obj)
	{
		m_dispplot->GetRenderer()->SetObjectVisible(*obj, show);
		m_dispplot->update();
	}
}



/**
 * show or hide the object labels
 */
void Dispersion3DDlg::ShowPlotLabels(bool show)
{
	m_dispplot->GetRenderer()->SetLabelsVisible(show);
	m_dispplot->update();
}



/**
 * choose between perspective or parallel projection
 */
void Dispersion3DDlg::SetPlotPerspectiveProjection(bool proj)
{
	m_dispplot->GetRenderer()->GetCamera().SetPerspectiveProjection(proj);
	m_dispplot->GetRenderer()->RequestViewportUpdate();
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * sets the camera's rotation angles
 */
void Dispersion3DDlg::SetPlotCameraRotation(t_real_gl phi, t_real_gl theta)
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
void Dispersion3DDlg::CentrePlotCameraOnObject()
{
	if(!m_cur_obj)
		return;

	t_mat_gl mat = m_dispplot->GetRenderer()->GetObjectMatrix(*m_cur_obj);

	// selected a band?
	auto band_iter = m_band_objs.find(*m_cur_obj);
	if(band_iter != m_band_objs.end())
	{
		// translate camera to the band's mean energy position
		t_real E_mean = GetMeanEnergy(band_iter->second) * m_E_scale->value();
		mat(2, 3) = E_mean;
	}

	m_dispplot->GetRenderer()->GetCamera().Centre(mat);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * centre camera on central position
 */
void Dispersion3DDlg::CentrePlotCamera()
{
	t_mat_gl matCentre = tl2::hom_translation<t_mat_gl>(
		m_cam_centre[0] * m_Q_scale2->value(),
		m_cam_centre[1] * m_Q_scale1->value(),
		m_cam_centre[2] * m_E_scale->value());

	m_dispplot->GetRenderer()->GetCamera().Centre(matCentre);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * switch between crystal and lab coordinates
 */
void Dispersion3DDlg::SetPlotCoordinateSystem(int which)
{
	m_dispplot->GetRenderer()->SetCoordSys(which);
}



/**
 * dialog is closing
 */
void Dispersion3DDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("dispersion3d/geo", saveGeometry());
		m_sett->setValue("dispersion3d/splitter", m_split_plot->saveState());
	}

	QDialog::accept();
}
