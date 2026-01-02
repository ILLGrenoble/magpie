/**
 * 3d plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
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

#ifndef __TL2_PLOT3D_DLG_H__
#define __TL2_PLOT3D_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>

#include <vector>
#include <array>
#include <tuple>
#include <unordered_map>

#include "tlibs2/libs/qt/glplot.h"
#include "defs.h"



/**
 * 3d plotter dialog
 */
class Plot3DDlg : public QDialog
{ Q_OBJECT
private:
	// column indices in surface table
	enum : int
	{
		COL_BC_BAND = 0,
		COL_BC_ACTIVE,
		NUM_COLS_BC,
	};


protected:
	using t_data_Q = std::tuple<t_vec /*0: Q*/, t_real /*1: E*/, t_real /*2: S*/,
		t_size /*3: Q_idx_1*/, t_size /*4: Q_idx_2*/,
		t_size /*5: degeneracy*/, bool /*6: valid*/>;
	using t_data_Qs = std::vector<t_data_Q>;
	using t_data_bands = std::vector<t_data_Qs>;


public:
	Plot3DDlg(QWidget *parent, std::shared_ptr<QSettings> sett);
	virtual ~Plot3DDlg();

	Plot3DDlg(const Plot3DDlg&) = delete;
	Plot3DDlg& operator=(const Plot3DDlg&) = delete;


protected:
	virtual void accept() override;

	// calculation functions
	void EnableCalculation(bool enable = true);
	void Calculate();
	void Plot(bool clear_settings = true);

	// calculation helper functions
	std::pair<t_size, t_size> NumValid(const t_data_Qs& data) const;
	bool IsValid(const t_data_Qs& data) const;
	t_real GetMeanZ(const t_data_Qs& data) const;
	t_real GetMeanZ(t_size band_idx) const;
	std::array<int, 3> GetBranchColour(t_size branch_idx, t_size num_branches) const;
	void ShowError(const QString& msg);

	// surface table functions
	void ClearSurfaces();
	void AddSurface(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsSurfaceEnabled(t_size idx) const;

	// export data
	void WriteHeader(std::ostream& ostr) const;
	void SaveData();
	void SaveScript();
	void SaveImage();

	// ------------------------------------------------------------------------
	// plotter interface
	void AfterPlotGLInitialisation();
	void PlotPickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl* posSphere);

	void PlotCameraHasUpdated();
	void CentrePlotCamera();
	void CentrePlotCameraOnObject();

	void PlotMouseClick(bool left, bool mid, bool right);
	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);

	void SetPlotCoordinateSystem(int which);
	void ShowPlotCoordCube(bool show);
	void ShowPlotLabels(bool show);
	void SetPlotPerspectiveProjection(bool proj);
	void SetPlotCameraRotation(t_real_gl phi, t_real_gl theta);
	// ------------------------------------------------------------------------


private:
	t_size m_x_count{}, m_y_count{};     // number of points along the two directions
	t_data_bands m_data{};               // data for all surfaces
	std::array<t_vec, 2> m_minmax_x{};   // minimum and maximum x values
	std::array<t_vec, 2> m_minmax_y{};   // minimum and maximum y values
	std::array<t_real, 2> m_minmax_z{};  // minimum and maximum z values

	std::shared_ptr<QSettings> m_sett{}; // program settings

	// dispersion
	tl2::GlPlot *m_dispplot{};           // 3d plotter
	std::optional<std::size_t> m_cur_obj{};
	std::unordered_map<std::size_t /*plot object*/, t_size /*surface index*/> m_band_objs{};
	t_vec_gl m_cam_centre{tl2::zero<t_vec_gl>(3)};

	QSplitter *m_split_plot{};
	QTableWidget *m_table_bands{};       // table listing the surfaces

	// dispersion options
	QDoubleSpinBox *m_xrange[2]{}, *m_yrange[2]{};
	QSpinBox *m_num_points[2]{};         // number of points on the (x, y) grid

	// context menus
	QMenu *m_context{};                  // general plot context menu
	QMenu *m_context_band{};             // context menu for the surfaces

	// plot options
	QDoubleSpinBox *m_x_scale{}, *m_y_scale{}, *m_z_scale{};
	QDoubleSpinBox *m_cam_phi{}, *m_cam_theta{};
	QCheckBox *m_perspective{};

	QPushButton *m_btn_start_stop{};     // start/stop calculation
	QProgressBar *m_progress{};          // progress bar
	QLabel *m_status{};                  // status bar

	bool m_calc_enabled{};               // enable calculations
	bool m_stop_requested{};             // stop running calculations


signals:
	void GlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);
};


#endif
