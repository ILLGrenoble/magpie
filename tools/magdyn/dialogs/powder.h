/**
 * magnetic dynamics -- powder excitations
 * @author Tobias Weber <tweber@ill.fr>
 * @date July 2026
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

#ifndef __MAG_DYN_POWDER_DLG_H__
#define __MAG_DYN_POWDER_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>
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



struct PowderData
{
	t_real momentum{};
	typename t_magdyn::t_histo histogram{};
};



/**
 * differentiation dialog
 */
class PowderDlg : public QDialog
{ Q_OBJECT
public:
	PowderDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~PowderDlg();

	PowderDlg(const PowderDlg&) = delete;
	PowderDlg& operator=(const PowderDlg&) = delete;

	// set kernel and Q path from main window
	void SetKernel(const t_magdyn* dyn);
	void SetDispersionQE(const t_vec_real& Qstart, const t_vec_real& Qend, t_real Estart, t_real Eend);
	void SetDispersionE(t_real Estart, t_real Eend);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);


	// ------------------------------------------------------------------------
	// powder panel
	QWidget* CreatePowderPanel();

	// plot functions
	void RescalePowderPlot();
	void ClearPowderPlot(bool replot = true);
	void PlotPowder();
	void SavePowderPlotFigure();
	void PowderPlotMouseMove(QMouseEvent *evt);
	void PowderPlotMousePress(QMouseEvent *evt);

	// calculation functions
	void EnablePowderCalculation(bool enable = true);
	void CalculatePowder();
	void SavePowderData();

	void SetPowderQE();
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	// from main dialog
	const t_magdyn *m_dyn{};            // main calculation kernel
	t_vec_real m_Qstart{}, m_Qend{};    // Qs from main window
	t_real m_Estart{0}, m_Eend{1};      // Es from main window

	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// powder panel
	std::vector<PowderData> m_data_powder{}; // all powder data

	QCustomPlot *m_plot_powder{};           // powder plotter
	QCPColorScale *m_plot_colour{};         // plot colour scale
	QCPColorMap *m_plot_map{};              // plot colour map

	t_real m_Q_min_powder{}, m_Q_max_powder{};  // range of Q component
	t_real m_E_min_powder{}, m_E_max_powder{};  // range of E component

	QDoubleSpinBox *m_Q_start_powder{};     // Q start
	QDoubleSpinBox *m_Q_end_powder{};       // Q end
	QDoubleSpinBox *m_E_start_powder{};     // E start
	QDoubleSpinBox *m_E_end_powder{};       // E end
	QSpinBox *m_num_Q_powder{};             // number of Q coordinates
	QSpinBox *m_num_Qvecs_powder{};         // number of Q vectors per coordinate
	QSpinBox *m_num_E_powder{};             // number of E coordinates
	QCheckBox *m_use_proj{};                // use orthogonal projector

	QPushButton *m_btnStartStop_powder{};   // start/stop calculation
	bool m_calcEnabled_powder{};            // enable calculations
	bool m_stopRequested_powder{};          // stop running calculations

	QProgressBar *m_progress_powder{};      // progress bar
	QMenu *m_menuPlot_powder{};             // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
