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
 * panel for visualising the scattering plane
 */
void MagDynDlg::CreateReciprocalPanel()
{
	m_reciprocalpanel = new QWidget(this);

	// 2d brillouin zone cut
	m_bzscene = new BZCutScene<t_vec_real, t_real>(m_reciprocalpanel);
	m_bzview = new BZCutView<t_vec_real, t_real>(m_bzscene, m_sett);


	// rotation axis
	m_recip_rot_axis[0] = new QDoubleSpinBox(m_reciprocalpanel);
	m_recip_rot_axis[1] = new QDoubleSpinBox(m_reciprocalpanel);
	m_recip_rot_axis[2] = new QDoubleSpinBox(m_reciprocalpanel);

	// rotation angle
	m_recip_rot_angle = new QDoubleSpinBox(m_reciprocalpanel);
	m_recip_rot_angle->setDecimals(3);
	m_recip_rot_angle->setMinimum(-360);
	m_recip_rot_angle->setMaximum(+360);
	m_recip_rot_angle->setSingleStep(0.1);
	m_recip_rot_angle->setValue(1.);
	//m_recip_rot_angle->setSuffix("\xc2\xb0");
	m_recip_rot_angle->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	for(int i = 0; i < 3; ++i)
	{
		m_recip_rot_axis[i]->setDecimals(4);
		m_recip_rot_axis[i]->setMinimum(-99.9999);
		m_recip_rot_axis[i]->setMaximum(+99.9999);
		m_recip_rot_axis[i]->setSingleStep(0.1);
		m_recip_rot_axis[i]->setValue(i == 2 ? 1. : 0.);
		m_recip_rot_axis[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	QPushButton *btn_rotate_ccw = new QPushButton(
		QIcon::fromTheme("object-rotate-left"),
		"Rotate CCW", m_reciprocalpanel);
	QPushButton *btn_rotate_cw = new QPushButton(
		QIcon::fromTheme("object-rotate-right"),
		"Rotate CW", m_reciprocalpanel);
	btn_rotate_ccw->setToolTip("Rotate the magnetic field in the counter-clockwise direction.");
	btn_rotate_cw->setToolTip("Rotate the magnetic field in the clockwise direction.");
	btn_rotate_ccw->setFocusPolicy(Qt::StrongFocus);
	btn_rotate_cw->setFocusPolicy(Qt::StrongFocus);


	// reduce path to first brillouin zone
	QPushButton *btnReduceBZ = new QPushButton("Reduce to First Zone", this);
	btnReduceBZ->setToolTip("Reduce the scan path to the first Brillouin zone.");
	btnReduceBZ->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// show 3d brillouin zone
	QPushButton *btnShowBZ = new QPushButton("3D Brillouin Zone...", this);
	btnShowBZ->setIcon(QIcon::fromTheme("applications-graphics"));
	btnShowBZ->setToolTip("Show a 3D view of the first nuclear Brillouin zone.");
	btnShowBZ->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);


	QFrame *sep1 = new QFrame(m_reciprocalpanel);
	sep1->setFrameStyle(QFrame::HLine);
	QFrame *sep2 = new QFrame(m_reciprocalpanel);
	sep2->setFrameStyle(QFrame::HLine);

	int y = 0;
	QGridLayout *grid = new QGridLayout(m_reciprocalpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);
	grid->addWidget(m_bzview, y++, 0, 1, 4);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);

	grid->addWidget(new QLabel("Rotate Qs:", m_reciprocalpanel), y++, 0, 1, 2);
	grid->addWidget(new QLabel("Axis (rlu):", m_reciprocalpanel), y, 0, 1, 1);
	grid->addWidget(m_recip_rot_axis[0], y, 1, 1, 1);
	grid->addWidget(m_recip_rot_axis[1], y, 2, 1, 1);
	grid->addWidget(m_recip_rot_axis[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Angle (\xc2\xb0):", m_reciprocalpanel), y, 0, 1, 1);
	grid->addWidget(m_recip_rot_angle, y, 1, 1, 1);
	grid->addWidget(btn_rotate_ccw, y, 2, 1, 1);
	grid->addWidget(btn_rotate_cw, y++, 3, 1, 1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep2, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(btnReduceBZ, y, 0, 1, 2);
	grid->addWidget(btnShowBZ, y++, 3, 1, 1);


	// context menu
	QMenu* context = m_bzview->GetContextMenu();
	QAction *acSetQi = new QAction("Set Start Q", context);
	QAction *acSetQf = new QAction("Set End Q", context);

	QAction* acFirst = context->actions().first();
	context->insertAction(acFirst, acSetQi);
	context->insertAction(acFirst, acSetQf);
	context->insertSeparator(acFirst);


	// signals
#ifdef BZ_USE_QT_SIGNALS
	connect(m_bzview, &BZCutView<t_vec, t_real>::SignalMouseCoordinates,
		this, &MagDynDlg::BZCutMouseMoved);
	connect(m_bzview, &BZCutView<t_vec, t_real>::SignalClickCoordinates,
		this, &MagDynDlg::BZCutMouseClicked);
#else
	m_bzview->AddMouseCoordinatesSlot([this](t_real x, t_real y)
	{
		this->BZCutMouseMoved(x, y);
	});
	m_bzview->AddClickCoordinatesSlot([this](int buttons, t_real x, t_real y)
	{
		this->BZCutMouseClicked(buttons, x, y);
	});
#endif

	connect(btnReduceBZ, &QAbstractButton::clicked, this, &MagDynDlg::ReducePathBZ);
	connect(btnShowBZ, &QAbstractButton::clicked, this, &MagDynDlg::ShowBZ3DDlg);

	connect(acSetQi, &QAction::triggered, [this]()
	{
		t_vec_real pos = m_bzview->GetClickedPosition(true);
		auto [QinvA, Qrlu] = m_bz.GetBZCutQ(pos[0], pos[1]);
		if(Qrlu.size() != 3)
			return;

		SetCoordinates(Qrlu, t_vec_real{}, true);
	});

	connect(acSetQf, &QAction::triggered, [this]()
	{
		t_vec_real pos = m_bzview->GetClickedPosition(true);
		auto [QinvA, Qrlu] = m_bz.GetBZCutQ(pos[0], pos[1]);
		if(Qrlu.size() != 3)
			return;

		SetCoordinates(t_vec_real{}, Qrlu, true);
	});

	connect(btn_rotate_ccw, &QAbstractButton::clicked, [this]()
	{
		t_vec_real axis = tl2::create<t_vec_real>(
		{
			(t_real)m_recip_rot_axis[0]->value(),
			(t_real)m_recip_rot_axis[1]->value(),
			(t_real)m_recip_rot_axis[2]->value(),
		});

		t_real angle = tl2::d2r<t_real>(m_recip_rot_angle->value());

		RotateDispersionQs(axis, angle);
	});

	connect(btn_rotate_cw, &QAbstractButton::clicked, [this]()
	{
		t_vec_real axis = tl2::create<t_vec_real>(
		{
			(t_real)m_recip_rot_axis[0]->value(),
			(t_real)m_recip_rot_axis[1]->value(),
			(t_real)m_recip_rot_axis[2]->value(),
		});

		t_real angle = tl2::d2r<t_real>(m_recip_rot_angle->value());

		RotateDispersionQs(axis, -angle);
	});


	m_tabs_recip->addTab(m_reciprocalpanel, "Scattering Plane");
}



/**
 * rotate the Q start and end vectors
 */
void MagDynDlg::RotateDispersionQs(const t_vec_real& axis_rlu, t_real angle)
{
	t_vec_real axis = axis_rlu;
	auto [Q_start, Q_end] = GetDispersionQ();

	const t_mat_real& xtalB = m_dyn.GetCrystalBTrafo();
	auto [xtalB_inv, inv_ok] = tl2::inv(xtalB);
	if(inv_ok)
	{
		axis = xtalB * axis;
		Q_start = xtalB * Q_start;
		Q_end = xtalB * Q_end;
	}

	t_mat_real R = tl2::rotation<t_mat_real, t_vec_real>(axis, angle, false);
	Q_start = R*Q_start;
	Q_end = R*Q_end;

	if(inv_ok)
	{
		Q_start = xtalB_inv * Q_start;
		Q_end = xtalB_inv * Q_end;
	}

	tl2::set_eps_0(Q_start, g_eps);
	tl2::set_eps_0(Q_end, g_eps);

	for(int i = 0; i < 3; ++i)
	{
		m_Q_start[i]->blockSignals(true);
		m_Q_start[i]->setValue(Q_start[i]);
		m_Q_start[i]->blockSignals(false);

		m_Q_end[i]->blockSignals(true);
		m_Q_end[i]->setValue(Q_end[i]);
		m_Q_end[i]->blockSignals(false);
	}

	DispersionQChanged(true);
};
