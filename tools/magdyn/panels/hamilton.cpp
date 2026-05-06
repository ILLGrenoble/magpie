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
 * shows the hamilton operator for a given Q position
 */
void MagDynDlg::CreateHamiltonPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_hamiltonianpanel = new QWidget(this);

	// hamiltonian
	m_hamiltonian = new QTextEdit(m_hamiltonianpanel);
	m_hamiltonian->setReadOnly(true);
	m_hamiltonian->setWordWrapMode(QTextOption::NoWrap);
	m_hamiltonian->setLineWrapMode(QTextEdit::NoWrap);
	m_hamiltonian->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// Q coordinates
	m_Q[0] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[1] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[2] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[0]->setToolTip("Momentum transfer component h (rlu).");
	m_Q[1]->setToolTip("Momentum transfer component k (rlu).");
	m_Q[2]->setToolTip("Momentum transfer component l (rlu).");

	for(int i = 0; i < 3; ++i)
	{
		m_Q[i]->setDecimals(4);
		m_Q[i]->setMinimum(-99.9999);
		m_Q[i]->setMaximum(+99.9999);
		m_Q[i]->setSingleStep(0.01);
		m_Q[i]->setValue(0.);
		m_Q[i]->setSuffix(" rlu");
		m_Q[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q[i]->setPrefix(hklPrefix[i]);
	}

	QGridLayout *grid = new QGridLayout(m_hamiltonianpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_hamiltonian, y++,0,1,4);
	grid->addWidget(new QLabel("Q:", m_hamiltonianpanel), y,0,1,1);
	grid->addWidget(m_Q[0], y,1,1,1);
	grid->addWidget(m_Q[1], y,2,1,1);
	grid->addWidget(m_Q[2], y++,3,1,1);

	// signals
	for(int i = 0; i < 3; ++i)
	{
		connect(m_Q[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			if(this->m_autocalc->isChecked())
				this->CalcHamiltonian();
		});
	}


	m_tabs_out->addTab(m_hamiltonianpanel, "Hamiltonian");
}

