/**
 * 2d plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date May 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * Magpie
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

#ifndef __TL2_PLOT2D_DLG_H__
#define __TL2_PLOT2D_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QSlider>
#include <QtWidgets/QProgressBar>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	#include <QtWidgets/QAction>
#else
	#include <QtGui/QAction>
#endif

#include <qcustomplot.h>
#include <vector>

#include "defs.h"

#include "tlibs2/libs/expr.h"



struct PlotData
{
	t_real x{};
	std::vector<t_real> ys{};
};



/**
 * 2d plot dialog
 */
class Plot2DDlg : public QDialog
{ Q_OBJECT
public:
	Plot2DDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~Plot2DDlg();

	Plot2DDlg(const Plot2DDlg&) = delete;
	Plot2DDlg& operator=(const Plot2DDlg&) = delete;


protected:
	virtual void accept() override;

	void ShowError(const char* msg);


	// ------------------------------------------------------------------------
	// main panel
	QWidget* CreatePanel();

	// plot functions
	void RescalePlot();
	void ClearPlot(bool replot = true);
	void Plot(bool clear_settings = true);
	void SavePlotFigure();
	void PlotMouseMove(QMouseEvent *evt);
	void PlotMousePress(QMouseEvent *evt);

	// index table functions
	void ClearIndices();
	void AddIndex(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsIndexEnabled(t_size idx) const;

	// calculation functions
	void Parse();
	void EnableCalculation(bool enable = true);
	void Calculate();
	void SaveData();
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// formulas
	QTextEdit *m_formulas{};            // formulas for the curves
	std::vector<tl2::ExprParser<t_real>> m_parsers{};
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// plot panel
	std::vector<PlotData> m_data{};  // all plot curve data

	QCustomPlot *m_plot{};           // plotter
	std::vector<QCPCurve*> m_curves{};  // plot curves

	QSplitter *m_split_plot{};
	QTableWidget *m_table{};         // table listing individual curves

	QDoubleSpinBox *m_x_start{};     // x start coordinate
	QDoubleSpinBox *m_x_end{};       // x end coordinate
	QSpinBox *m_num_x{};             // number of x coordinates

	QDoubleSpinBox *m_t_start{};     // t start parameter
	QDoubleSpinBox *m_t_end{};       // t end parameter
	QSpinBox *m_num_t{};             // number of t subdivisions
	QSlider *m_slider_t{};           // t parameter position

	QPushButton *m_btnStartStop{};   // start/stop calculation
	bool m_calcEnabled{};            // enable calculations
	bool m_stopRequested{};          // stop running calculations

	QProgressBar *m_progress{};      // progress bar
	QMenu *m_menuPlot{};             // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
