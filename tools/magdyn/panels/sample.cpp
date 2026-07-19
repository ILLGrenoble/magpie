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
 * allows the user to specify the sample properties
 */
void MagDynDlg::CreateSamplePanel()
{
	m_samplepanel = new QWidget(this);

	// crystal lattice and angles
	static const char* latticestr[] = { "a = ", "b = ", "c = " };
	for(int i = 0; i < 3; ++i)
	{
		m_xtallattice[i] = new QDoubleSpinBox(m_samplepanel);
		m_xtallattice[i]->setDecimals(3);
		m_xtallattice[i]->setMinimum(0.001);
		m_xtallattice[i]->setMaximum(99.999);
		m_xtallattice[i]->setSingleStep(0.1);
		m_xtallattice[i]->setValue(5);
		m_xtallattice[i]->setPrefix(latticestr[i]);
		//m_xtallattice[i]->setSuffix(" \xe2\x84\xab");
		m_xtallattice[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	static const char* anlesstr[] = { "α = ", "β = ", "γ = " };
	for(int i = 0; i < 3; ++i)
	{
		m_xtalangles[i] = new QDoubleSpinBox(m_samplepanel);
		m_xtalangles[i]->setDecimals(2);
		m_xtalangles[i]->setMinimum(0.01);
		m_xtalangles[i]->setMaximum(180.);
		m_xtalangles[i]->setSingleStep(0.1);
		m_xtalangles[i]->setValue(90);
		m_xtalangles[i]->setPrefix(anlesstr[i]);
		//m_xtalangles[i]->setSuffix("\xc2\xb0");
		m_xtalangles[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	// space groups
	m_comboSG = new QComboBox(m_samplepanel);
	m_comboSG->setFocusPolicy(Qt::StrongFocus);
	m_comboSG->setToolTip("List of space groups.");

	m_checkFilterSG = new QCheckBox("Filter Space Groups:", m_samplepanel);
	m_checkFilterSG->setChecked(true);

	m_editFilterSG = new QLineEdit(m_samplepanel);
	m_editFilterSG->setPlaceholderText("Space group filter.");

	// scattering plane
	static const char* recipstr[] = { "h = ", "k = ", "l = " };
	for(int i = 0; i < 6; ++i)
	{
		m_scatteringplane[i] = new QDoubleSpinBox(m_samplepanel);
		m_scatteringplane[i]->setDecimals(3);
		m_scatteringplane[i]->setMinimum(-99.999);
		m_scatteringplane[i]->setMaximum(99.999);
		m_scatteringplane[i]->setSingleStep(1.);
		m_scatteringplane[i]->setValue(0);
		m_scatteringplane[i]->setPrefix(recipstr[i % 3]);
		//m_scatteringplane[i]->setSuffix(" rlu");
		m_scatteringplane[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}
	m_scatteringplane[0]->setValue(1);
	m_scatteringplane[4]->setValue(1);

	// magnetic form factor
	m_ffact = new QPlainTextEdit(m_samplepanel);
	m_ffact->setPlaceholderText("Enter magnetic form factor formula."
		" The free variable is 'Q' or 's'.");
	m_ffact->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// maximum number of form factors
	m_num_ffacts = new QSpinBox(m_samplepanel);
	m_num_ffacts->setMinimum(1);
	m_num_ffacts->setMaximum(10000);
	m_num_ffacts->setValue(1);
	m_num_ffacts->setPrefix("num = ");
	m_num_ffacts->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_ffacts->setToolTip("Maximum number of magnetic form factors.");

	// index of currently shown form factor
	m_cur_ffact = new QSpinBox(m_samplepanel);
	m_cur_ffact->setMinimum(0);
	m_cur_ffact->setMaximum(0);
	m_cur_ffact->setValue(0);
	m_cur_ffact->setPrefix("cur = ");
	m_cur_ffact->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_cur_ffact->setToolTip("Index of currently shown magnetic form factor.");

	// form factor table
	QPushButton *btn_set_ffact = nullptr;
	if(m_ff.GetFormfactorCount())
	{
		m_combo_ffacts = new QComboBox(m_samplepanel);
		m_combo_ffacts->setFocusPolicy(Qt::StrongFocus);
		m_combo_ffacts->setToolTip("List of magnetic form factors.");

		m_editFilterFFacts = new QLineEdit(m_samplepanel);
		m_editFilterFFacts->setPlaceholderText("Form factor filter.");

		PopulateFormFactors();

		btn_set_ffact = new QPushButton("Set", m_samplepanel);
		btn_set_ffact->setToolTip("Set the form factor term from the currently selected magnetic ion.");
	}

	QPushButton *btn_ffact_j0 = new QPushButton("<j0> Templ.", m_samplepanel);
	QPushButton *btn_ffact_j2 = new QPushButton("<j0-2> Templ.", m_samplepanel);
	QPushButton *btn_ffact_j4 = new QPushButton("<j0-4> Templ.", m_samplepanel);
	QPushButton *btn_ffact_plot = new QPushButton("Plot...", m_samplepanel);
	btn_ffact_j0->setToolTip("Add a template <j0> term.");
	btn_ffact_j2->setToolTip("Add template <j0> and <j2> terms.");
	btn_ffact_j4->setToolTip("Add template <j0>, <j2>, and <j4> terms.");
	btn_ffact_plot->setToolTip("Plot the form factors.");

	for(QPushButton *btn : { btn_set_ffact, btn_ffact_j0, btn_ffact_j2, btn_ffact_j4, btn_ffact_plot })
	{
		if(!btn)
			continue;
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	QGridLayout *grid = new QGridLayout(m_samplepanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;

	// crystal
	grid->addWidget(new QLabel("Crystal Definition", m_samplepanel), y++, 0, 1, 4);
	grid->addWidget(new QLabel("Lattice (\xe2\x84\xab):", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_xtallattice[0], y, 1, 1, 1);
	grid->addWidget(m_xtallattice[1], y, 2, 1, 1);
	grid->addWidget(m_xtallattice[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Angles (\xc2\xb0):", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_xtalangles[0], y, 1, 1, 1);
	grid->addWidget(m_xtalangles[1], y, 2, 1, 1);
	grid->addWidget(m_xtalangles[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Space Group:", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_comboSG, y++, 1, 1, 3);
	grid->addWidget(m_checkFilterSG, y, 0, 1, 1);
	grid->addWidget(m_editFilterSG, y++, 1, 1, 3);

	// separator
	QFrame *sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);
	grid->addWidget(sep1, y++,0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);

	// scattering plane
	grid->addWidget(new QLabel("Scattering Plane", m_samplepanel), y++, 0, 1, 4);
	grid->addWidget(new QLabel("Vector 1 (rlu):", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_scatteringplane[0], y, 1, 1, 1);
	grid->addWidget(m_scatteringplane[1], y, 2, 1, 1);
	grid->addWidget(m_scatteringplane[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Vector 2 (rlu):", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_scatteringplane[3], y, 1, 1, 1);
	grid->addWidget(m_scatteringplane[4], y, 2, 1, 1);
	grid->addWidget(m_scatteringplane[5], y++, 3, 1, 1);

	// separator
	QFrame *sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);
	grid->addWidget(sep2, y++, 0, 1, 4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++, 0, 1, 1);

	// magnetic form factor formula
	grid->addWidget(new QLabel("Magnetic Form Factors", m_samplepanel), y, 0, 1, 1);
	grid->addWidget(m_num_ffacts, y, 2, 1, 1);
	grid->addWidget(m_cur_ffact, y++, 3, 1, 1);
	grid->addWidget(new QLabel("Formula, F_M(Q or s) = ", m_samplepanel), y++, 0, 1, 4);
	grid->addWidget(m_ffact, y++, 0, 1, 4);
	if(m_combo_ffacts && btn_set_ffact)
	{
		grid->addWidget(new QLabel("Select Ion:", m_samplepanel), y, 0, 1, 1);
		grid->addWidget(m_combo_ffacts, y, 1, 1, 1);
		grid->addWidget(m_editFilterFFacts, y, 2, 1, 1);
		grid->addWidget(btn_set_ffact, y++, 3, 1, 1);
	}
	grid->addWidget(btn_ffact_j0, y, 0, 1, 1);
	grid->addWidget(btn_ffact_j2, y, 1, 1, 1);
	grid->addWidget(btn_ffact_j4, y, 2, 1, 1);
	grid->addWidget(btn_ffact_plot, y++, 3, 1, 1);
	//grid->addItem(new QSpacerItem(8, 8,
	//	QSizePolicy::Minimum, QSizePolicy::Expanding),
	//	y++, 0, 1, 1);

	auto calc_all = [this]()
	{
		m_needsBZCalc = true;

		if(this->m_autocalc->isChecked())
		{
			this->CalcAll();
			this->CalcBZ();
			this->DispersionQChanged(false);
		}
	};

	connect(m_comboSG,
		static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		this, &MagDynDlg::CalcBZ);

#if QT_VERSION < QT_VERSION_CHECK(6, 7, 0)
	connect(m_checkFilterSG, &QCheckBox::stateChanged, [this](int checked)
#else
	connect(m_checkFilterSG, &QCheckBox::checkStateChanged, [this](Qt::CheckState checked)
#endif
	{
		m_editFilterSG->setEnabled(int(checked) != 0);
		PopulateSpaceGroups();
	});

	connect(m_editFilterSG, &QLineEdit::textChanged, [this]() { PopulateSpaceGroups(); });

	for(int i = 0; i < 2*3; ++i)
	{
		if(i < 3)
		{
			connect(m_xtallattice[i],
				static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
				calc_all);
			connect(m_xtalangles[i],
				static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
				calc_all);
		}

		connect(m_scatteringplane[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}

	// magnetic form factors
	m_ffacts.resize(m_num_ffacts->value());

	connect(m_editFilterFFacts, &QLineEdit::textChanged, [this]() { PopulateFormFactors(); });
	connect(btn_ffact_plot, &QPushButton::clicked, this, &MagDynDlg::ShowFormFactorDlg);
	connect(btn_ffact_j0, &QPushButton::clicked, [this]()
	{
		// for formulas (and coefficients to use), see: https://mcphase.github.io/webpage/manual/node164.html .
		m_ffact->setPlainText(
			"A = 1; a = 0;\n"
			"B = 0; b = 0;\n"
			"C = 0; c = 0;\n"
			"D = 0;\n"
			"A*exp(-a*s2) + B*exp(-b*s2) + C*exp(-c*s2) + D");
	});
	connect(btn_ffact_j2, &QPushButton::clicked, [this]()
	{
		// for formulas (and coefficients to use), see: https://mcphase.github.io/webpage/manual/node164.html .
		m_ffact->setPlainText(
			"A0 = 1; a0 = 0;\n"
			"B0 = 0; b0 = 0;\n"
			"C0 = 0; c0 = 0;\n"
			"D0 = 0;\n"
			"A2 = 0; a2 = 0;\n"
			"B2 = 0; b2 = 0;\n"
			"C2 = 0; c2 = 0;\n"
			"D2 = 0; g = 1;\n"
			" A0*exp(-a0*s2) + B0*exp(-b0*s2) + C0*exp(-c0*s2) + D0 +\n"
			"(A2*exp(-a2*s2) + B2*exp(-b2*s2) + C2*exp(-c2*s2) + D2) * (2/g - 1) * s2");
	});
	connect(btn_ffact_j4, &QPushButton::clicked, [this]()
	{
		// for formulas (and coefficients to use), see: https://mcphase.github.io/webpage/manual/node164.html .
		m_ffact->setPlainText(
			"A0 = 1; a0 = 0;\n"
			"B0 = 0; b0 = 0;\n"
			"C0 = 0; c0 = 0;\n"
			"D0 = 0;\n"
			"A2 = 0; a2 = 0;\n"
			"B2 = 0; b2 = 0;\n"
			"C2 = 0; c2 = 0;\n"
			"D2 = 0; g = 1;\n"
			"A4 = 0; a4 = 0;\n"
			"B4 = 0; b4 = 0;\n"
			"C4 = 0; c4 = 0;\n"
			"D4 = 0;\n"
			" A0*exp(-a0*s2) + B0*exp(-b0*s2) + C0*exp(-c0*s2) + D0 +\n"
			"(A2*exp(-a2*s2) + B2*exp(-b2*s2) + C2*exp(-c2*s2) + D2) * (2/g - 1) * s2 +\n"
			"(A4*exp(-a4*s2) + B4*exp(-b4*s2) + C4*exp(-c4*s2) + D4) * (2/g - 1) * s2");
	});

	// maximum number of form factors changed
	connect(m_num_ffacts,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this](int num)
	{
		m_ffacts.resize(num);
		m_cur_ffact->setMaximum(num - 1);

		PlotFormFactors();
	});

	// current form factor changed
	connect(m_cur_ffact,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this](int idx)
	{
		// display selected form factor
		if(idx < (int)m_ffacts.size())
			m_ffact->setPlainText(m_ffacts[idx].c_str());
		else
			m_ffact->clear();
	});

	connect(m_ffact, &QPlainTextEdit::textChanged, [this, calc_all]()
	{
		int cur_idx = m_cur_ffact->value();
		if(cur_idx >= (int)m_ffacts.size())
			m_ffacts.resize(cur_idx + 1);

		// save current form factor text
		m_ffacts[cur_idx] = m_ffact->toPlainText().toStdString();

		PlotFormFactors();
		calc_all();
	});

	// set the selected ion's form factor
	if(m_combo_ffacts && btn_set_ffact)
	{
		connect(btn_set_ffact, &QAbstractButton::clicked, [this]()
		{
			std::string ion_name = m_combo_ffacts->currentText().toStdString();
			const auto* ion = m_ff.GetFormfactor(ion_name);
			if(!ion)
			{
				m_ffact->setPlainText("");
				ShowError(("Could not find " + ion_name + " ion.").c_str(), true);
				return;
			}

			m_ffact->setPlainText(ion->to_string().c_str());
		});
	}

	m_tabs_setup->addTab(m_samplepanel, "Crystal");
}

