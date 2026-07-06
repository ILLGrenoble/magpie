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
 * plots the dispersion relation for a given Q path
 */
void MagDynDlg::CreateDispersionPanel()
{
	m_disppanel = new QWidget(this);

	// plotter
	m_plot = new QCustomPlot(m_disppanel);
	m_plot->setFont(this->font());
	m_plot->xAxis->setLabel("Q (rlu)");
	m_plot->yAxis->setLabel("E (meV)");
	m_plot->setInteraction(QCP::iRangeDrag, true);
	m_plot->setInteraction(QCP::iRangeZoom, true);
	m_plot->setSelectionRectMode(QCP::srmZoom);
	m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// start and stop coordinates
	m_Q_start[0] = new QDoubleSpinBox(m_disppanel);
	m_Q_start[1] = new QDoubleSpinBox(m_disppanel);
	m_Q_start[2] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[0] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[1] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[2] = new QDoubleSpinBox(m_disppanel);

	m_Q_start[0]->setToolTip("Dispersion initial momentum transfer, h_i (rlu).");
	m_Q_start[1]->setToolTip("Dispersion initial momentum transfer, k_i (rlu).");
	m_Q_start[2]->setToolTip("Dispersion initial momentum transfer, l_i (rlu).");
	m_Q_end[0]->setToolTip("Dispersion final momentum transfer, h_f (rlu).");
	m_Q_end[1]->setToolTip("Dispersion final momentum transfer, k_f (rlu).");
	m_Q_end[2]->setToolTip("Dispersion final momentum transfer, l_f (rlu).");

	// number of Q points in the plot
	m_num_points = new QSpinBox(m_disppanel);
	m_num_points->setMinimum(1);
	m_num_points->setMaximum(99999);
	m_num_points->setValue(512);
	m_num_points->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_num_points->setToolTip("Number of Q points in the plot.");

	// scaling factor for weights
	for(auto** comp : { &m_weight_scale, &m_weight_min, &m_weight_max })
	{
		*comp = new QDoubleSpinBox(m_disppanel);
		(*comp)->setDecimals(4);
		(*comp)->setMinimum(0.);
		(*comp)->setMaximum(+9999.9999);
		(*comp)->setSingleStep(0.1);
		(*comp)->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	m_weight_scale->setValue(1.);
	m_weight_min->setValue(0.);
	m_weight_max->setValue(99.);
	m_weight_min->setMinimum(-1.);	// -1: disable clamping
	m_weight_max->setMinimum(-1.);	// -1: disable clamping
	m_weight_min->setToolTip("Minimum spectral weight for clamping.");
	m_weight_max->setToolTip("Maximum spectral weight for clamping.");
	m_weight_scale->setToolTip("Spectral weight scaling factor.");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_start[i]->setDecimals(4);
		m_Q_start[i]->setMinimum(-99.9999);
		m_Q_start[i]->setMaximum(+99.9999);
		m_Q_start[i]->setSingleStep(0.01);
		m_Q_start[i]->setValue(i == 0 ? -1. : 0.);
		//m_Q_start[i]->setSuffix(" rlu");
		m_Q_start[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_start[i]->setPrefix(hklPrefix[i]);

		m_Q_end[i]->setDecimals(4);
		m_Q_end[i]->setMinimum(-99.9999);
		m_Q_end[i]->setMaximum(+99.9999);
		m_Q_end[i]->setSingleStep(0.01);
		m_Q_end[i]->setValue(i == 0 ? 1. : 0.);
		//m_Q_end[i]->setSuffix(" rlu");
		m_Q_end[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_end[i]->setPrefix(hklPrefix[i]);
	}

	QGridLayout *grid = new QGridLayout(m_disppanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	if(m_plot)
		grid->addWidget(m_plot, y++,0,1,4);
	grid->addWidget(new QLabel("Start Q (rlu):", m_disppanel), y,0,1,1);
	grid->addWidget(m_Q_start[0], y,1,1,1);
	grid->addWidget(m_Q_start[1], y,2,1,1);
	grid->addWidget(m_Q_start[2], y++,3,1,1);
	grid->addWidget(new QLabel("End Q (rlu):", m_disppanel), y,0,1,1);
	grid->addWidget(m_Q_end[0], y,1,1,1);
	grid->addWidget(m_Q_end[1], y,2,1,1);
	grid->addWidget(m_Q_end[2], y++,3,1,1);
	grid->addWidget(new QLabel("Q Count:", m_disppanel), y,0,1,1);
	grid->addWidget(m_num_points, y,1,1,1);
	grid->addWidget(new QLabel("Weight Scale:", m_disppanel), y,2,1,1);
	grid->addWidget(m_weight_scale, y++,3,1,1);
	grid->addWidget(new QLabel("Min. Weight:", m_disppanel), y,0,1,1);
	grid->addWidget(m_weight_min, y,1,1,1);
	grid->addWidget(new QLabel("Max. Weight:", m_disppanel), y,2,1,1);
	grid->addWidget(m_weight_max, y++,3,1,1);

	// signals
	for(int i = 0; i < 3; ++i)
	{
		connect(m_Q_start[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
			{
				DispersionQChanged(true);
			});

		connect(m_Q_end[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
			{
				DispersionQChanged(true);
			});
	}

	connect(m_num_points,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this]()
		{
			DispersionQChanged(true);
		});
	
	for(auto* comp : {m_weight_scale, m_weight_min, m_weight_max})
	{
		connect(comp,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			// update graph weights
			for(GraphWithWeights* graph : m_graphs)
				graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
			if(m_plot)
				m_plot->replot();
		});
	}

	if(m_plot)
	{
		connect(m_plot, &QCustomPlot::mouseMove, this, &MagDynDlg::PlotMouseMove);
		connect(m_plot, &QCustomPlot::mousePress, this, &MagDynDlg::PlotMousePress);
	}


	m_tabs_out->addTab(m_disppanel, "Dispersion");
}

