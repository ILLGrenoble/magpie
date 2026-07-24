/**
 * magnetic dynamics -- form factor plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date May 2026
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

#ifndef __MAG_DYN_FFACT_DLG_H__
#define __MAG_DYN_FFACT_DLG_H__

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



struct FormFactorData
{
	t_real momentum{};              // Q
	std::vector<t_real> ffacts{};   // Re(f_M(Q))
	std::vector<t_real> ffacts2{};  // |f_M(Q)|^2
};



/**
 * form factor dialog
 */
class FormFactorDlg : public QDialog
{ Q_OBJECT
public:
	FormFactorDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~FormFactorDlg();

	FormFactorDlg(const FormFactorDlg&) = delete;
	FormFactorDlg& operator=(const FormFactorDlg&) = delete;

	void PlotFormFactors(const std::vector<std::string>& ffacts);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);


	// ------------------------------------------------------------------------
	// form factor panel
	QWidget* CreateFormFactorPanel();

	// plot functions
	void RescaleFormFactorPlot();
	void ClearFormFactorPlot(bool replot = true);
	void PlotFormFactors(bool clear_settings = true);
	void SaveFormFactorPlotFigure();
	void FormFactorPlotMouseMove(QMouseEvent *evt);
	void FormFactorPlotMousePress(QMouseEvent *evt);

	// index table functions
	void ClearFormFactorIndices();
	void AddFormFactorIndex(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsFormFactorIndexEnabled(t_size idx) const;

	// calculation functions
	void EnableFormFactorCalculation(bool enable = true);
	void CalculateFormFactors();
	void SaveFormFactorData();
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	// from main dialog
	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// formulas for the magnetic form factors
	std::vector<tl2::ExprParser<t_cplx>> m_parsers_ff{};
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// form factors panel
	std::vector<FormFactorData> m_data_ff{};  // all form factor data

	QCustomPlot *m_plot_ff{};           // form factor plotter
	std::vector<QCPCurve*> m_curves_ff{};  // form factor plot curves

	QSplitter *m_split_plot_ff{};
	QTableWidget *m_table_ff{};         // table listing individual form factors

	QDoubleSpinBox *m_Q_start_ff{};     // Q start coordinate
	QDoubleSpinBox *m_Q_end_ff{};       // Q end coordinate
	QSpinBox *m_num_Q_ff{};             // number of Q coordinates
	QCheckBox *m_ff_squared{};          // plot |F|^2 instead of F

	QPushButton *m_btnStartStop_ff{};   // start/stop calculation
	bool m_calcEnabled_ff{};            // enable calculations
	bool m_stopRequested_ff{};          // stop running calculations

	QProgressBar *m_progress_ff{};      // progress bar
	QMenu *m_menuPlot_ff{};             // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
