/**
 * brillouin zone tool -- 3d plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __BZTOOL_PLOTTER_H__
#define __BZTOOL_PLOTTER_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QMenu>

#include <memory>
#include <vector>
#include <unordered_map>
#include <string>

#include "globals.h"

#include "tlibs2/libs/qt/glplot.h"


class BZPlotDlg : public QDialog
{ Q_OBJECT
public:
	BZPlotDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	~BZPlotDlg() = default;

	BZPlotDlg(const BZPlotDlg&) = delete;
	BZPlotDlg& operator=(const BZPlotDlg&) = delete;

	void Clear();

	void SetABTrafo(const t_mat_bz& crystA, const t_mat_bz& crystB);
	void AddVoronoiVertex(const t_vec_bz& pos);
	void AddBraggPeak(const t_vec_bz& pos);
	void AddTriangles(const std::vector<t_vec_bz>& vecs,
		const std::vector<std::size_t> *faceindices = nullptr);
	void SetPlane(const t_vec_bz& norm, t_real d);

	void SetEps(t_real eps);
	void SetPrecGui(int prec);

	QMenu* GetContextMenu();
	const t_vec_bz& GetClickedPosition(bool right_button = false) const;


protected:
	void SetStatusMsg(const std::string& msg);
	void ShowCoordCrossLab(bool show);
	void ShowCoordCrossXtal(bool show);
	void ShowLabels(bool show);
	void ShowPlane(bool show);
	void ShowQVertices(bool show);
	void SetCameraRotation(t_real_gl phi, t_real_gl theta);
	void SetPerspectiveProjection(bool proj);

	void CameraHasUpdated();

	void PlotMouseClick(bool left, bool mid, bool right);
	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);
	void PickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl* posSphere);
	void AfterGLInitialisation();

	void SaveImage();

	virtual void accept() override;


private:
	t_real m_eps{ 1e-7 };
	int m_prec_gui{ 4 };

	t_mat_bz m_crystA{ tl2::unit<t_mat_bz>(3) };  // crystal A matrix
	t_mat_bz m_crystB{ tl2::unit<t_mat_bz>(3) };  // crystal B matrix

	QSettings *m_sett{};

	std::shared_ptr<tl2::GlPlot> m_plot{};
	std::size_t m_sphere{ 0 }, m_plane{ 0 };

	QLabel *m_status{};
	QCheckBox *m_show_coordcross_lab{}, *m_show_coordcross_xtal{};
	QCheckBox *m_show_labels{}, *m_show_plane{}, *m_show_Qs{};
	QCheckBox *m_perspective{};
	QDoubleSpinBox *m_cam_phi{}, *m_cam_theta{};

	QMenu *m_context{};                       // plot context menu

	long m_curPickedObj{ -1 };                // current 3d bz object
	std::vector<std::size_t> m_objsBragg{};   // Bragg peak plot objects
	std::vector<std::size_t> m_objsVoronoi{}; // Voronoi vertex plot objects
	std::vector<std::size_t> m_objsBZ{};      // BZ triangle plot objects

	t_vec_bz m_cur_Qrlu{};                    // current cursor position
	t_vec_bz m_clicked_Q_rlu[2]{};            // clicked position

	// maps an object to its triangles' BZ face indices
	std::unordered_map<std::size_t, std::vector<std::size_t>> m_objFaceIndices{};


signals:
	void NeedRecalc();

	void GlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);
};


#endif
