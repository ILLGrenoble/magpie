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
 * input for sample environment parameters (field, temperature)
 */
void MagDynDlg::CreateSampleEnvPanel()
{
	m_sampleenviropanel = new QWidget(this);

	// field magnitude
	m_field_mag = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_mag->setDecimals(3);
	m_field_mag->setMinimum(0);
	m_field_mag->setMaximum(+99.999);
	m_field_mag->setSingleStep(0.1);
	m_field_mag->setValue(0.);
	m_field_mag->setPrefix("|B| = ");
	m_field_mag->setSuffix(" T");
	m_field_mag->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	// predefined field directions
	QPushButton *btnDirs = new QPushButton(m_reciprocalpanel);
	btnDirs->setText("Set Direction");
	btnDirs->setToolTip("Set the magnetic field direction to a predefined value.");
	btnDirs->setSizePolicy(QSizePolicy{QSizePolicy::Preferred, QSizePolicy::Preferred});
	QMenu *menuDirs = new QMenu(btnDirs);
	QAction *acDirPlane[] = {
		new QAction("Scattering Plane Vector 1", menuDirs),
		new QAction("Scattering Plane Vector 2", menuDirs),
		new QAction("Scattering Plane Normal", menuDirs) };
	menuDirs->addAction(acDirPlane[0]);
	menuDirs->addAction(acDirPlane[1]);
	menuDirs->addAction(acDirPlane[2]);
	btnDirs->setMenu(menuDirs);

	// field direction
	m_field_dir[0] = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_dir[1] = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_dir[2] = new QDoubleSpinBox(m_sampleenviropanel);

	// align spins along field (field-polarised state)
	m_align_spins = new QCheckBox(
		"Align Spins Along Field", m_sampleenviropanel);
	m_align_spins->setChecked(false);
	m_align_spins->setFocusPolicy(Qt::StrongFocus);

	// align spins along field (field-polarised state)
	m_keep_spin_signs = new QCheckBox(
		"Keep the Spin Senses", m_sampleenviropanel);
	m_keep_spin_signs->setChecked(false);
	m_keep_spin_signs->setFocusPolicy(Qt::StrongFocus);

	// predefined rotation axes
	QPushButton *btnAxes = new QPushButton(m_reciprocalpanel);
	btnAxes->setText("Set Axis");
	btnAxes->setToolTip("Set the magnetic field rotation axis to a predefined value.");
	btnAxes->setSizePolicy(QSizePolicy{QSizePolicy::Preferred, QSizePolicy::Preferred});
	QMenu *menuAxes = new QMenu(btnAxes);
	QAction *acPlane[] = {
		new QAction("Scattering Plane Vector 1", menuAxes),
		new QAction("Scattering Plane Vector 2", menuAxes),
		new QAction("Scattering Plane Normal", menuAxes) };
	menuAxes->addAction(acPlane[0]);
	menuAxes->addAction(acPlane[1]);
	menuAxes->addAction(acPlane[2]);
	btnAxes->setMenu(menuAxes);

	// rotation axis
	m_rot_axis[0] = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_axis[1] = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_axis[2] = new QDoubleSpinBox(m_sampleenviropanel);

	// rotation angle
	m_rot_angle = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_angle->setDecimals(3);
	m_rot_angle->setMinimum(-360);
	m_rot_angle->setMaximum(+360);
	m_rot_angle->setSingleStep(0.1);
	m_rot_angle->setValue(90.);
	//m_rot_angle->setSuffix("\xc2\xb0");
	m_rot_angle->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	QPushButton *btn_rotate_ccw = new QPushButton(
		QIcon::fromTheme("object-rotate-left"),
		"Rotate CCW", m_sampleenviropanel);
	QPushButton *btn_rotate_cw = new QPushButton(
		QIcon::fromTheme("object-rotate-right"),
		"Rotate CW", m_sampleenviropanel);
	btn_rotate_ccw->setToolTip("Rotate the magnetic field in the counter-clockwise direction.");
	btn_rotate_cw->setToolTip("Rotate the magnetic field in the clockwise direction.");
	btn_rotate_ccw->setFocusPolicy(Qt::StrongFocus);
	btn_rotate_cw->setFocusPolicy(Qt::StrongFocus);
	btn_rotate_ccw->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	btn_rotate_cw->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);


	// table with saved fields
	m_fieldstab = new QTableWidget(m_sampleenviropanel);
	m_fieldstab->setShowGrid(true);
	m_fieldstab->setAlternatingRowColors(true);
	m_fieldstab->setSortingEnabled(true);
	m_fieldstab->setMouseTracking(true);
	m_fieldstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_fieldstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_fieldstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_fieldstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing()*1.25 + 4);
	m_fieldstab->verticalHeader()->setVisible(true);

	m_fieldstab->setColumnCount(NUM_FIELD_COLS);
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_H, new QTableWidgetItem{"Bh"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_K, new QTableWidgetItem{"Bk"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_L, new QTableWidgetItem{"Bl"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_MAG, new QTableWidgetItem{"|B|"});

	m_fieldstab->setColumnWidth(COL_FIELD_H, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_K, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_L, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_MAG, 150);
	m_fieldstab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAddField = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_sampleenviropanel);
	QPushButton *btnDelField = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_sampleenviropanel);
	QPushButton *btnFieldUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_sampleenviropanel);
	QPushButton *btnFieldDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_sampleenviropanel);

	btnAddField->setToolTip("Add a magnetic field.");
	btnDelField->setToolTip("Delete selected magnetic field(s).");
	btnFieldUp->setToolTip("Move selected magnetic field(s) up.");
	btnFieldDown->setToolTip("Move selected magnetic field(s) down.");

	QPushButton *btnSetField = new QPushButton("Set Field", m_sampleenviropanel);
	btnSetField->setToolTip("Set the selected field as the currently active one.");

	for(QPushButton *btn : { btnAddField, btnDelField, btnFieldUp, btnFieldDown, btnSetField })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}


	// table context menu
	QMenu *menuTableContext = new QMenu(m_fieldstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field Before Current", this,
		[this]() { this->AddFieldTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field After Current", this,
		[this]() { this->AddFieldTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Field", this,
		[this]() { this->AddFieldTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Field", this,
		[this]() { this->DelTabItem(m_fieldstab); });
	menuTableContext->addSeparator();
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Set As Current Field", this,
		[this]() { this->SetCurrentField(); });


	// table context menu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_fieldstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field", this,
		[this]() { this->AddFieldTabItem(-1,
			m_field_dir[0]->value(),
			m_field_dir[1]->value(),
			m_field_dir[2]->value(),
			m_field_mag->value()); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Field", this,
		[this]() { this->DelTabItem(m_fieldstab); });


	// temperature
	m_temperature = new QDoubleSpinBox(m_sampleenviropanel);
	m_temperature->setDecimals(2);
	m_temperature->setMinimum(0);
	m_temperature->setMaximum(+999.99);
	m_temperature->setSingleStep(0.1);
	m_temperature->setValue(300.);
	m_temperature->setPrefix("T = ");
	m_temperature->setSuffix(" K");
	m_temperature->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	for(int i = 0; i < 3; ++i)
	{
		m_field_dir[i]->setDecimals(4);
		m_field_dir[i]->setMinimum(-99.9999);
		m_field_dir[i]->setMaximum(+99.9999);
		m_field_dir[i]->setSingleStep(0.1);
		m_field_dir[i]->setValue(i == 2 ? 1. : 0.);
		m_field_dir[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

		m_rot_axis[i]->setDecimals(4);
		m_rot_axis[i]->setMinimum(-99.9999);
		m_rot_axis[i]->setMaximum(+99.9999);
		m_rot_axis[i]->setSingleStep(0.1);
		m_rot_axis[i]->setValue(i == 2 ? 1. : 0.);
		m_rot_axis[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	m_field_dir[0]->setPrefix("Bh = ");
	m_field_dir[1]->setPrefix("Bk = ");
	m_field_dir[2]->setPrefix("Bl = ");


	// grid
	QGridLayout *grid = new QGridLayout(m_sampleenviropanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Magnetic Field:", m_sampleenviropanel), y++, 0, 1, 2);
	grid->addWidget(new QLabel("Magnitude:", m_sampleenviropanel), y, 0, 1, 1);
	grid->addWidget(m_field_mag, y, 1, 1, 1);
	grid->addWidget(btnDirs, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Direction (rlu):", m_sampleenviropanel), y, 0, 1, 1);
	grid->addWidget(m_field_dir[0], y, 1, 1, 1);
	grid->addWidget(m_field_dir[1], y, 2, 1, 1);
	grid->addWidget(m_field_dir[2], y++, 3, 1, 1);
	grid->addWidget(m_align_spins, y, 0, 1, 2);
	grid->addWidget(m_keep_spin_signs, y++, 2, 1, 2);

	QFrame *sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	QFrame *sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);
	QFrame *sep3 = new QFrame(m_sampleenviropanel);
	sep3->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep1, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);

	grid->addWidget(new QLabel("Rotate Magnetic Field:", m_sampleenviropanel), y, 0, 1, 2);
	grid->addWidget(btnAxes, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Axis (rlu):", m_sampleenviropanel), y, 0, 1, 1);
	grid->addWidget(m_rot_axis[0], y, 1, 1, 1);
	grid->addWidget(m_rot_axis[1], y, 2, 1, 1);
	grid->addWidget(m_rot_axis[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Angle (\xc2\xb0):", m_sampleenviropanel), y, 0, 1, 1);
	grid->addWidget(m_rot_angle, y, 1, 1, 1);
	grid->addWidget(btn_rotate_ccw, y, 2, 1, 1);
	grid->addWidget(btn_rotate_cw, y++, 3, 1, 1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep2, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);

	grid->addWidget(new QLabel("Saved Fields:", m_sampleenviropanel), y++, 0, 1, 4);
	grid->addWidget(m_fieldstab, y, 0, 1, 4);
	grid->addWidget(btnAddField, ++y, 0, 1, 1);
	grid->addWidget(btnDelField, y, 1, 1, 1);
	grid->addWidget(btnFieldUp, y, 2, 1, 1);
	grid->addWidget(btnFieldDown, y++, 3, 1, 1);
	grid->addWidget(btnSetField, y++, 3, 1, 1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep3, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);

	grid->addWidget(new QLabel("Temperature:", m_sampleenviropanel), y, 0, 1, 1);
	grid->addWidget(m_temperature, y++, 1, 1, 1);

	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};


	// signals
	connect(m_field_mag,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		calc_all);

	for(int i = 0; i < 3; ++i)
	{
		connect(m_field_dir[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}

	connect(m_temperature,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		calc_all);

	connect(m_align_spins, &QCheckBox::toggled, calc_all);
	connect(m_keep_spin_signs, &QCheckBox::toggled, calc_all);

	connect(btn_rotate_ccw, &QAbstractButton::clicked, [this]()
	{
		t_vec_real axis = tl2::create<t_vec_real>(
		{
			(t_real)m_rot_axis[0]->value(),
			(t_real)m_rot_axis[1]->value(),
			(t_real)m_rot_axis[2]->value(),
		});

		t_real angle = tl2::d2r<t_real>(m_rot_angle->value());

		RotateField(axis, angle);
	});

	connect(btn_rotate_cw, &QAbstractButton::clicked, [this]()
	{
		t_vec_real axis = tl2::create<t_vec_real>(
		{
			(t_real)m_rot_axis[0]->value(),
			(t_real)m_rot_axis[1]->value(),
			(t_real)m_rot_axis[2]->value(),
		});

		t_real angle = tl2::d2r<t_real>(m_rot_angle->value());

		RotateField(axis, -angle);
	});

	connect(btnAddField, &QAbstractButton::clicked,
		[this]() { this->AddFieldTabItem(-1,
			m_field_dir[0]->value(),
			m_field_dir[1]->value(),
			m_field_dir[2]->value(),
			m_field_mag->value()); });
	connect(btnDelField, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_fieldstab); });
	connect(btnFieldUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_fieldstab); });
	connect(btnFieldDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_fieldstab); });

	connect(btnSetField, &QAbstractButton::clicked,
		[this]() { this->SetCurrentField(); });

	connect(m_fieldstab, &QTableWidget::itemSelectionChanged, this, &MagDynDlg::FieldsSelectionChanged);
	connect(m_fieldstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_fieldstab, menuTableContext, menuTableContextNoItem, pt);
	});


	for(int i = 0; i < 3; ++i)
	{
		connect(acDirPlane[i], &QAction::triggered, [this, i, calc_all]()
		{
			const t_vec_real& vec = m_dyn.GetScatteringPlane()[i];

			for(int j = 0; j < 3; ++j)
			{
				m_field_dir[j]->blockSignals(true);
				m_field_dir[j]->setValue(vec[j]);
				m_field_dir[j]->blockSignals(false);
			}

			calc_all();
		});

		connect(acPlane[i], &QAction::triggered, [this, i]()
		{
			const t_vec_real& vec = m_dyn.GetScatteringPlane()[i];

			m_rot_axis[0]->setValue(vec[0]);
			m_rot_axis[1]->setValue(vec[1]);
			m_rot_axis[2]->setValue(vec[2]);
		});
	}


	m_tabs_in->insertTab(1, m_sampleenviropanel, "Environment");
}



/**
 * a field value has been selected
 */
void MagDynDlg::FieldsSelectionChanged()
{
	QList<QTableWidgetItem*> selected = m_fieldstab->selectedItems();
	if(selected.size() == 0)
		return;

	const QTableWidgetItem* item = *selected.begin();
	m_fields_cursor_row = item->row();
}



/**
 * rotate the direction of the magnetic field
 */
void MagDynDlg::RotateField(const t_vec_real& axis_rlu, t_real angle)
{
	t_vec_real axis = axis_rlu;

	t_vec_real B = tl2::create<t_vec_real>(
	{
		(t_real)m_field_dir[0]->value(),
		(t_real)m_field_dir[1]->value(),
		(t_real)m_field_dir[2]->value(),
	});

	const t_mat_real& xtalB = m_dyn.GetCrystalBTrafo();
	auto [xtalB_inv, inv_ok] = tl2::inv(xtalB);
	if(inv_ok)
	{
		axis = xtalB * axis;
		B = xtalB * B;
	}

	t_mat_real R = tl2::rotation<t_mat_real, t_vec_real>(axis, angle, false);
	B = R*B;

	if(inv_ok)
		B = xtalB_inv * B;

	tl2::set_eps_0(B, g_eps);

	for(int i = 0; i < 3; ++i)
	{
		m_field_dir[i]->blockSignals(true);
		m_field_dir[i]->setValue(B[i]);
		m_field_dir[i]->blockSignals(false);
	}

	if(m_autocalc->isChecked())
		CalcAll();
};
