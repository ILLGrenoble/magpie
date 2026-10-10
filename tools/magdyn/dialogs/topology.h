/**
 * magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * magpie & mag-core
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

#ifndef __MAG_DYN_TOPO_DLG_H__
#define __MAG_DYN_TOPO_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QTabWidget>
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
#include <tuple>

#include "gui_defs.h"



struct BerryCurvatureData
{
	t_vec3_real momentum{};

	std::vector<t_real> energies{};    // E
	std::vector<t_real> weights{};     // S_perp
	std::vector<t_cplx> curvatures{};  // B
};



/**
 * topology dialog
 */
class TopologyDlg : public QDialog
{ Q_OBJECT
public:
	TopologyDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~TopologyDlg();

	TopologyDlg(const TopologyDlg&) = delete;
	TopologyDlg& operator=(const TopologyDlg&) = delete;

	// set kernel and Q path from main window
	void SetKernel(const t_magdyn* dyn);
	void SetDispersionQ(const t_vec3_real& Qstart, const t_vec3_real& Qend);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);

	// ------------------------------------------------------------------------
	// berry curvature tab
	QWidget* CreateBerryCurvaturePanel();

	// plot functions
	void RescaleBerryCurvaturePlot();
	void ClearBerryCurvaturePlot(bool replot = true);
	void PlotBerryCurvature(bool clear_settings = true);
	void SaveBerryCurvaturePlotFigure();
	void BerryCurvaturePlotMouseMove(QMouseEvent *evt);
	void BerryCurvaturePlotMousePress(QMouseEvent *evt);

	// band table functions
	void ClearBerryCurvatureBands();
	void AddBerryCurvatureBand(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsBerryCurvatureBandEnabled(t_size idx) const;

	// calculation functions
	void EnableBerryCurvatureCalculation(bool enable = true);
	void CalculateBerryCurvature();
	void SaveBerryCurvatureData();

	void SetBerryCurvatureQ();
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// berry curvature map tab
	QWidget* CreateBerryCurvatureMapPanel();

	// helper functions
	std::tuple<t_vec3_real, t_vec3_real, t_vec3_real> GetMapQVectors() const;
	t_vec3_real GetMapQFromIndices(std::size_t idx1, std::size_t idx2) const;
	void MapQFromMainQ();

	// plot functions
	void ClearBerryCurvatureMapPlot(bool replot = true);
	void RescaleBerryCurvatureMapPlot();
	void PlotBerryCurvatureMap();
	void BerryCurvatureMapPlotMouseMove(QMouseEvent *evt);
	void BerryCurvatureMapPlotMousePress(QMouseEvent *evt);
	void SaveBerryCurvatureMapData();
	void SaveBerryCurvatureMapPlotFigure();

	// calculation functions
	void EnableBerryCurvatureMapCalculation(bool enable = true);
	void CalculateBerryCurvatureMap();
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	// from main dialog
	const t_magdyn *m_dyn{};            // main calculation kernel
	t_vec3_real m_Qstart{}, m_Qend{};   // Qs from main window

	QTabWidget *m_tabs{};               // tabs
	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// berry curvature tab
	std::vector<BerryCurvatureData> m_data_bc{};  // all (non-filtered) berry curvature data

	QCustomPlot *m_plot_bc{};           // berry curvature plotter
	std::vector<QCPCurve*> m_curves_bc{};  // berry curvature plot curves
	t_size m_Q_idx_bc{};                // index of dominant Q component
	t_real m_Q_min_bc{}, m_Q_max_bc{};  // range of dominant Q component

	QSplitter *m_split_plot_bc{};
	QTableWidget *m_table_bands_bc{};   // table listing the magnon bands

	QDoubleSpinBox *m_Q_start_bc[3]{};  // Q start coordinate
	QDoubleSpinBox *m_Q_end_bc[3]{};    // Q end coordinate
	QSpinBox *m_num_Q_bc{};             // number of Q coordinates

	QCheckBox *m_B_filter_enable_bc{};  // switch to enable maximum B value
	QDoubleSpinBox *m_B_filter_bc{};    // maximum B value

	QCheckBox *m_S_filter_enable_bc{};  // switch to enable minimum S(Q,E) value
	QDoubleSpinBox *m_S_filter_bc{};    // minimum S(Q,E) value

	QSpinBox *m_coords_bc[2]{};         // berry curvature component indices
	QAction *m_imag_bc{};               // imaginary or real components?
	QCheckBox *m_only_pos_E_bc{};       // ignore magnon annihilation?
	//QCheckBox *m_enforce_commutator_bc{};

	QPushButton *m_btnStartStop_bc{};   // start/stop calculation
	bool m_calcEnabled_bc{};            // enable calculations
	std::atomic<bool> m_stopRequested_bc{};  // stop running calculations

	QProgressBar *m_progress_bc{};      // progress bar
	QMenu *m_menuPlot_bc{};             // context menu for plot
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// berry curvature map tab
	QCustomPlot *m_plot_bcm{};            // berry curvature plotter
	QCPColorScale *m_plot_colour_bcm{};   // plot colour scale
	QCPColorMap *m_plot_map_bcm{};        // plot colour map

	QDoubleSpinBox *m_Q_origin_bcm[3]{};  // Q start coordinate
	QDoubleSpinBox *m_Q_dir1_bcm[3]{};    // Q direction vector 1
	QDoubleSpinBox *m_Q_dir2_bcm[3]{};    // Q direction vector 2
	QSpinBox *m_num_Q_bcm[2]{};           // number of Q points

	QCheckBox *m_B_filter_enable_bcm{};   // switch to enable maximum B value
	QDoubleSpinBox *m_B_filter_bcm{};     // maximum B value

	QCheckBox *m_S_filter_enable_bcm{};   // switch to enable minimum S(Q,E) value
	QDoubleSpinBox *m_S_filter_bcm{};     // minimum S(Q,E) value

	QSpinBox *m_coords_bcm[2]{};          // berry curvature component indices
	QSpinBox *m_band_bcm{};               // band index
	QAction *m_imag_bcm{};                // imaginary or real components?
	QCheckBox *m_only_pos_E_bcm{};        // ignore magnon annihilation?

	QPushButton *m_btnStartStop_bcm{};    // start/stop calculation
	bool m_calcEnabled_bcm{};             // enable calculations
	std::atomic<bool> m_stopRequested_bcm{};  // stop running calculations

	QProgressBar *m_progress_bcm{};       // progress bar
	QMenu *m_menuPlot_bcm{};              // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
