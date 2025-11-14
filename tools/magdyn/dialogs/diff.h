/**
 * magnetic dynamics -- differentiation / group velocity calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 202t
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

#ifndef __MAG_DYN_DIFF_DLG_H__
#define __MAG_DYN_DIFF_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QProgressBar>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	#include <QtWidgets/QAction>
#else
	#include <QtGui/QAction>
#endif

#include <qcustomplot.h>
#include <vector>

#include "gui_defs.h"



struct GroupVelocityData
{
	t_vec_real momentum{};

	std::vector<t_real> energies{};   // E
	std::vector<t_real> weights{};    // S_perp
	std::vector<t_real> velocities{}; // group velocities, dE/dq
};



/**
 * topology dialog
 */
class DiffDlg : public QDialog
{ Q_OBJECT
public:
	DiffDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~DiffDlg();

	DiffDlg(const DiffDlg&) = delete;
	DiffDlg& operator=(const DiffDlg&) = delete;

	// set kernel and Q path from main window
	void SetKernel(const t_magdyn* dyn);
	void SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);


	// ------------------------------------------------------------------------
	// group velocity panel
	QWidget* CreateGroupVelocityPanel();

	// plot functions
	void RescaleGroupVelocityPlot();
	void ClearGroupVelocityPlot(bool replot = true);
	void PlotGroupVelocity(bool clear_settings = true);
	void SaveGroupVelocityPlotFigure();
	void GroupVelocityPlotMouseMove(QMouseEvent *evt);
	void GroupVelocityPlotMousePress(QMouseEvent *evt);

	// band table functions
	void ClearGroupVelocityBands();
	void AddGroupVelocityBand(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsGroupVelocityBandEnabled(t_size idx) const;

	// calculation functions
	void EnableGroupVelocityCalculation(bool enable = true);
	void CalculateGroupVelocity();
	void SaveGroupVelocityData();

	void SetGroupVelocityQ();
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	// from main dialog
	const t_magdyn *m_dyn{};            // main calculation kernel
	t_vec_real m_Qstart{}, m_Qend{};    // Qs from main window

	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// group velocity panel
	std::vector<GroupVelocityData> m_data_gv{};  // all (non-filtered) groupvelocity data

	QCustomPlot *m_plot_gv{};           // group velocity plotter
	std::vector<QCPCurve*> m_curves_gv{};  // group velocity plot curves
	t_size m_Q_idx_gv{};                // index of dominant Q component
	t_real m_Q_min_gv{}, m_Q_max_gv{};  // range of dominant Q component

	QSplitter *m_split_plot_gv{};
	QTableWidget *m_table_bands_gv{};   // table listing the magnon bands

	QCheckBox *m_only_pos_E_gv{};       // ignore magnon annihilation?
	QSpinBox *m_diff{};                 // differentiation order

	QDoubleSpinBox *m_Q_start_gv[3]{};  // Q start coordinate
	QDoubleSpinBox *m_Q_end_gv[3]{};    // Q end coordinate
	QSpinBox *m_num_Q_gv{};             // number of Q coordinates

	QCheckBox *m_v_filter_enable_gv{};  // switch to enable maximum B value
	QDoubleSpinBox *m_v_filter_gv{};    // maximum B value

	QCheckBox *m_S_filter_enable_gv{};  // switch to enable minimum S(Q,E) value
	QDoubleSpinBox *m_S_filter_gv{};    // minimum S(Q,E) value

	QPushButton *m_btnStartStop_gv{};   // start/stop calculation
	bool m_calcEnabled_gv{};            // enable calculations
	bool m_stopRequested_gv{};          // stop running calculations

	QProgressBar *m_progress_gv{};      // progress bar
	QMenu *m_menuPlot_gv{};             // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
