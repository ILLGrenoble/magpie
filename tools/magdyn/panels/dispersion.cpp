/**
 * magnetic dynamics -- plotting of the dispersion
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

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



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
		[this](int val)
		{
			m_Qidx->setMaximum(val - 1);
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



/**
 * draw the calculated dispersion curve
 */
void MagDynDlg::PlotDispersion()
{
	if(!m_plot)
		return;

	m_plot->clearPlottables();
	m_graphs.clear();

	const bool plot_channels = m_plot_channels->isChecked();
	const bool plot_degeneracies = m_plot_degeneracies->isChecked();

	// create plot curves for each matrix element
	if(plot_channels)
	{
		// TODO: handle case if channels AND degeneracies are plotted

		const QColor colChannel[2*3*3]
		{
			QColor(std::lerp(0x00, 0xff, 1.),   0x00, 0x00),  // xx, real
			QColor(std::lerp(0x00, 0xff, 0.75), 0x00, 0x00),  // xy, real
			QColor(std::lerp(0x00, 0xff, 0.5),  0x00, 0x00),  // xz, real

			QColor(0x00, std::lerp(0x00, 0xff, 0.5),  0x00),  // yx, real
			QColor(0x00, std::lerp(0x00, 0xff, 1.),   0x00),  // yy, real
			QColor(0x00, std::lerp(0x00, 0xff, 0.57), 0x00),  // yz, real

			QColor(0x00, 0x00, std::lerp(0x00, 0xff, 0.5)),   // zx, real
			QColor(0x00, 0x00, std::lerp(0x00, 0xff, 0.75)),  // zy, real
			QColor(0x00, 0x00, std::lerp(0x00, 0xff, 1.)),    // zz, real

			QColor(std::lerp(0x00, 0xff, 1.),   0x22, 0x22),  // xx, imag
			QColor(std::lerp(0x00, 0xff, 0.75), 0x22, 0x22),  // xy, imag
			QColor(std::lerp(0x00, 0xff, 0.5),  0x22, 0x22),  // xz, imag

			QColor(0x22, std::lerp(0x00, 0xff, 0.5),  0x22),  // yx, imag
			QColor(0x22, std::lerp(0x00, 0xff, 1.),   0x22),  // yy, imag
			QColor(0x22, std::lerp(0x00, 0xff, 0.57), 0x22),  // yz, imag

			QColor(0x22, 0x22, std::lerp(0x00, 0xff, 0.5)),   // zx, imag
			QColor(0x22, 0x22, std::lerp(0x00, 0xff, 0.75)),  // zy, imag
			QColor(0x22, 0x22, std::lerp(0x00, 0xff, 1.)),    // zz, imag
		};

		// iterate channels
		for(t_size imag_elem = 0; imag_elem < 2; ++imag_elem)
		for(t_size i = 0; i < 3; ++i)
		for(t_size j = 0; j < 3; ++j)
		{
			int idx = imag_elem ? 3*3 : 0;  // channel index
			idx += i*3 + j;

			m_matrixelems_dlg->SetActive(i, j, imag_elem == 0,
				!tl2::equals_0(m_ws_total_channel[idx], g_eps));

			if(!m_matrixelems_dlg->IsChecked(i, j, imag_elem == 0))
				continue;

			GraphWithWeights *graph = new GraphWithWeights(m_plot->xAxis, m_plot->yAxis);
			QPen pen = graph->pen();
			pen.setColor(colChannel[idx]);
			pen.setWidthF(1.);
			graph->setPen(pen);
			graph->setBrush(QBrush(pen.color(), Qt::SolidPattern));
			graph->setLineStyle(QCPGraph::lsNone);
			graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, m_weight_scale->value()));
			graph->setAntialiased(true);
			graph->setData(m_qs_data_channel[idx], m_Es_data_channel[idx], true /*already sorted*/);
			graph->SetWeights(m_ws_data_channel[idx]);
			graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
			graph->SetWeightAsPointSize(m_plot_weights_pointsize->isChecked());
			graph->SetWeightAsAlpha(m_plot_weights_alpha->isChecked());
			m_graphs.push_back(graph);
		}
	}
	else
	{
		GraphWithWeights *graph = new GraphWithWeights(m_plot->xAxis, m_plot->yAxis);

		// dispersion colour
		int col_comp[3] = { 0, 0, 0xff };           // default colour
		tl2::get_colour<int>(g_colPlot, col_comp);  // get actual colour
		const QColor colFull(col_comp[0], col_comp[1], col_comp[2]);

		QPen pen = graph->pen();
		pen.setColor(colFull);

		// colour degeneracies
		if(plot_degeneracies)
		{
			// colour for degenerate dispersion points
			int col_comp_degen[3] = { 0xff, 0, 0 };      // default colour
			tl2::get_colour<int>(g_colPlotDegen, col_comp_degen);  // get actual colour
			const QColor colDegen(col_comp_degen[0], col_comp_degen[1], col_comp_degen[2]);

			graph->AddColour(colFull);
			graph->AddColour(colDegen);

			graph->SetColourIndices(m_degen_data);
		}

		pen.setWidthF(1.);
		graph->setPen(pen);
		graph->setBrush(QBrush(pen.color(), Qt::SolidPattern));
		graph->setLineStyle(QCPGraph::lsNone);
		graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, m_weight_scale->value()));
		graph->setAntialiased(true);
		graph->setData(m_qs_data, m_Es_data, true /*already sorted*/);
		graph->SetWeights(m_ws_data);
		graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
		graph->SetWeightAsPointSize(m_plot_weights_pointsize->isChecked());
		graph->SetWeightAsAlpha(m_plot_weights_alpha->isChecked());
		m_graphs.push_back(graph);
	}

	// set labels
	const char* Q_label[]{ "h (rlu)", "k (rlu)", "l (rlu)" };
	m_plot->xAxis->setLabel(Q_label[m_Q_idx]);

	// set Q plot range
	m_plot->xAxis->setRange(m_Q_min, m_Q_max);

	// set E plot range
	auto [min_E_iter, max_E_iter] = std::minmax_element(m_Es_data.begin(), m_Es_data.end());
	if(min_E_iter != m_Es_data.end() && max_E_iter != m_Es_data.end())
	{
		t_real E_range = *max_E_iter - *min_E_iter;

		m_E_min = *min_E_iter - E_range*0.05;
		m_E_max = *max_E_iter + E_range*0.05;
		//if(ignore_annihilation)
		//	m_E_min = t_real(0);

	}
	else
	{
		m_E_min = 0.;
		m_E_max = 1.;
	}

	m_plot->yAxis->setRange(m_E_min, m_E_max);

	// set font
	m_plot->setFont(font());
	m_plot->xAxis->setLabelFont(font());
	m_plot->yAxis->setLabelFont(font());
	m_plot->xAxis->setTickLabelFont(font());
	m_plot->yAxis->setTickLabelFont(font());

	m_plot->replot();

	// - change of E range (due to the "ignore annihilation" option) is not handled in DispersionQChanged()
	// - new E range is only available after CalcDispersion()
	if(m_powder_dlg)
	{
		auto [ E_start, E_end ] = GetDispersionE();
		m_powder_dlg->SetDispersionE(E_start, E_end);
	}
}



/**
 * mouse move event of the plot
 */
void MagDynDlg::PlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	auto [Q1, Q2] = GetDispersionQ();

	// plot (Q, E)
	const t_real Q = m_plot->xAxis->pixelToCoord(evt->pos().x());
	const t_real E = m_plot->yAxis->pixelToCoord(evt->pos().y());

	// lerp parameter
	t_size Q_idx = m_Q_idx;
	if(Q_idx > 2)
		Q_idx = 2;
	const t_real t = (Q - Q1[Q_idx]) / (Q2[Q_idx] - Q1[Q_idx]);
	const t_vec_real Qvec = tl2::create<t_vec_real>({
		std::lerp(Q1[0], Q2[0], t),
		std::lerp(Q1[1], Q2[1], t),
		std::lerp(Q1[2], Q2[2], t) });

	const t_mat_real& B = m_dyn.GetCrystalBTrafo();
	const t_vec_real Qvec_invA = B * Qvec;
	const t_real Q_invA = tl2::norm(Qvec_invA);

	// write status
	QString status("Q = (%1, %2, %3) rlu, |Q| = %4 Å⁻¹, E = %5 meV.");
	status = status
		.arg(Qvec[0], 0, 'g', g_prec_gui)
		.arg(Qvec[1], 0, 'g', g_prec_gui)
		.arg(Qvec[2], 0, 'g', g_prec_gui)
		.arg(Q_invA, 0, 'g', g_prec_gui)
		.arg(E, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * mouse button has been pressed in the plot
 */
void MagDynDlg::PlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuDisp)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuDisp->popup(pos);
		evt->accept();
	}
}
