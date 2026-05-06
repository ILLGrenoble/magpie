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

	// reduce path to first brillouin zone
	QPushButton *btnReduceBZ = new QPushButton("Reduce to First Zone", this);
	btnReduceBZ->setToolTip("Reduce the scan path to the first Brillouin zone.");
	btnReduceBZ->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// show 3d brillouin zone
	QPushButton *btnShowBZ = new QPushButton("3D Brillouin Zone...", this);
	btnShowBZ->setIcon(QIcon::fromTheme("applications-graphics"));
	btnShowBZ->setToolTip("Show a 3D view of the first nuclear Brillouin zone.");
	btnShowBZ->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	QGridLayout *grid = new QGridLayout(m_reciprocalpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	grid->addWidget(m_bzview, 0, 0, 1, 4);
	grid->addWidget(btnReduceBZ, 1, 0, 1, 1);
	grid->addWidget(btnShowBZ, 1, 3, 1, 1);

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

	m_tabs_recip->addTab(m_reciprocalpanel, "Scattering Plane");
}
