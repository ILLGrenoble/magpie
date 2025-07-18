/**
 * brillouin zone tool
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

#include "bz.h"

#include <QtWidgets/QMessageBox>

#include <iostream>
#include <sstream>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/expr.h"
#include "tlibs2/libs/qt/helper.h"

using namespace tl2_ops;


/**
 * precalculates Q vectors for BZ cut calculation
 */
void BZDlg::SetDrawOrder(int order, bool recalc)
{
	m_bzcalc.CalcPeaks(order, false /*invA*/, true /*cut*/);

	if(recalc)
		CalcBZCut();
}


/**
 * precalculates Q vectors for BZ calculation
 */
void BZDlg::SetCalcOrder(int order, bool recalc)
{
	m_bzcalc.CalcPeaks(order, false /*invA*/, false /*cut*/);

	if(recalc)
		CalcBZ();
}


/**
 * calculate crystal B matrix
 */
bool BZDlg::CalcB(bool full_recalc)
{
	if(m_ignoreCalc)
		return false;

	t_real a = tl2::stoval<t_real>(m_editA->text().toStdString());
	t_real b = tl2::stoval<t_real>(m_editB->text().toStdString());
	t_real c = tl2::stoval<t_real>(m_editC->text().toStdString());
	t_real alpha = tl2::stoval<t_real>(m_editAlpha->text().toStdString());
	t_real beta = tl2::stoval<t_real>(m_editBeta->text().toStdString());
	t_real gamma = tl2::stoval<t_real>(m_editGamma->text().toStdString());

	if(tl2::equals<t_real>(a, 0., g_eps_bz) || a <= 0. ||
		tl2::equals<t_real>(b, 0., g_eps_bz) || b <= 0. ||
		tl2::equals<t_real>(c, 0., g_eps_bz) || c <= 0. ||
		tl2::equals<t_real>(alpha, 0., g_eps_bz) || alpha <= 0. ||
		tl2::equals<t_real>(beta, 0., g_eps_bz) || beta <= 0. ||
		tl2::equals<t_real>(gamma, 0., g_eps_bz) || gamma <= 0.)
	{
		//QMessageBox::critical(this, "Brillouin Zones", "Error: Invalid lattice.");
		m_status->setText("<font color=\"red\">Error: Invalid lattice.</font>");
		return false;
	}

	t_mat_bz crystB = tl2::B_matrix<t_mat_bz>(a, b, c,
		tl2::d2r<t_real>(alpha), tl2::d2r<t_real>(beta), tl2::d2r<t_real>(gamma));

	bool ok = true;
	t_mat_bz crystA = tl2::unit<t_mat_bz>(3);
	std::tie(crystA, ok) = tl2::inv(crystB);
	if(!ok)
	{
		QMessageBox::critical(this, "Brillouin Zones", "Error: Cannot invert B matrix.");
		return false;
	}

	m_crystA = crystA * t_real(2)*tl2::pi<t_real>;
	m_crystB = std::move(crystB);

	if(m_dlgPlot)
		m_dlgPlot->SetABTrafo(m_crystA, m_crystB);

	m_status->setText((std::string("<font color=\"green\">")
		+ "Crystal B matrix calculated successfully."
		+ std::string("</font>")).c_str());

	if(full_recalc)
		CalcBZ();

	return true;
}


/**
 * calculate brillouin zone
 */
bool BZDlg::CalcBZ(bool full_recalc)
{
	if(m_ignoreCalc)
		return false;

	// set up bz calculator
	m_bzcalc.SetEps(g_eps_bz);
	m_bzcalc.SetSymOps(GetSymOps(true), true);
	m_bzcalc.SetCrystalA(m_crystA);
	m_bzcalc.SetCrystalB(m_crystB);
	m_bzcalc.CalcPeaksInvA();

	// calculate bz
	bool ok = m_bzcalc.CalcBZ();

	if(m_dlgPlot)
	{
		m_dlgPlot->SetEps(g_eps_bz);
		m_dlgPlot->SetPrecGui(g_prec_gui_bz);

		// clear old plot
		m_dlgPlot->Clear();

		// add gamma point
		std::size_t idx000 = m_bzcalc.Get000Peak();
		const std::vector<t_vec_bz>& Qs_invA = m_bzcalc.GetPeaks(true);
		if(idx000 < Qs_invA.size())
			m_dlgPlot->AddBraggPeak(Qs_invA[idx000]);

		// add voronoi vertices forming the vertices of the BZ
		for(const t_vec_bz& voro : m_bzcalc.GetVertices())
			m_dlgPlot->AddVoronoiVertex(voro);

		// add voronoi bisectors
		m_dlgPlot->AddTriangles(m_bzcalc.GetAllTriangles(),
			&m_bzcalc.GetAllTrianglesFaceIndices());
	}

	// set bz description string
	m_descrBZ = m_bzcalc.Print(g_prec_bz);
	m_descrBZJSON = m_bzcalc.PrintJSON(g_prec_bz);

	m_status->setText((std::string("<font color=\"green\">")
		+ "Brillouin zone calculated successfully."
		+ std::string("</font>")).c_str());

	bool cut_ok = true;
	if(full_recalc)
		cut_ok = CalcBZCut();
	else
		UpdateBZDescription();

	return ok && cut_ok;
}


/**
 * calculate brillouin zone cut
 */
bool BZDlg::CalcBZCut()
{
	if(m_ignoreCalc || !m_bzcalc.GetTriangles().size())
		return false;

	// get plane coordinate system
	t_vec_bz vec1_rlu = tl2::create<t_vec_bz>({
		m_cutX->value(),
		m_cutY->value(),
		m_cutZ->value()
	});
	t_vec_bz norm_rlu = tl2::create<t_vec_bz>({
		m_cutNX->value(),
		m_cutNY->value(),
		m_cutNZ->value()
	});

	t_real d_rlu = m_cutD->value();

	// calculate cut
	bool ok = m_bzcalc.CalcBZCut(vec1_rlu, norm_rlu, d_rlu, m_acCutHull->isChecked());

	// draw cut
	m_bzscene->SetEps(g_eps_bz);
	m_bzscene->SetPrecGui(g_prec_gui_bz);
	m_bzscene->ClearAll();
	m_bzscene->AddCut(m_bzcalc.GetCutLines(false));
	m_bzscene->AddPeaks(m_bzcalc.GetPeaksOnPlane(true), &m_bzcalc.GetPeaksOnPlane(false));
	m_bzview->Centre();

	// get description of the cutting plane
	m_descrBZCut = m_bzcalc.PrintCut(g_prec_bz);

	// update calculation results
	if(m_dlgPlot)
	{
		m_dlgPlot->SetPlane(
			tl2::col<t_mat_bz, t_vec_bz>(m_bzcalc.GetCutPlane(), 2),  // normal
			m_bzcalc.GetCutPlaneD());                           // distance
	}

	UpdateBZDescription();
	bool formulas_ok = CalcFormulas();

	if(formulas_ok)
	{
		m_status->setText((std::string("<font color=\"green\">")
			+ "Brillouin zone cut calculated successfully."
			+ std::string("</font>")).c_str());
	}

	return ok;
}


/**
 * evaluate the formulas in the table and plot them
 */
bool BZDlg::CalcFormulas()
{
	auto [ min_x, min_y, max_x, max_y ] = m_bzcalc.GetCutMinMax();

	m_bzscene->ClearCurves();
	if(max_x < min_x)
		return false;

	t_real plane_d = m_cutD->value() * m_bzcalc.GetCutNormScale();

	bool all_ok = true;
	std::vector<std::string> formulas = GetFormulas();
	for(std::size_t formula_idx = 0; formula_idx < formulas.size(); ++formula_idx)
	{
		const std::string& formula = formulas[formula_idx];

		try
		{
			tl2::ExprParser<t_real> parser;
			//parser.SetDebug(true);
			parser.SetAutoregisterVariables(false);
			parser.register_var("x", 0.);
			parser.register_var("Qx", 0.);
			parser.register_var("Qy", 0.);

			if(bool ok = parser.parse(formula); !ok)
				continue;

			int num_pts = 512;
			t_real x_delta = (max_x - min_x) / t_real(num_pts);

			std::vector<t_vec_bz> curve;
			curve.reserve(num_pts);

			for(t_real x = min_x; x <= max_x; x += x_delta)
			{
				t_vec_bz QinvA = m_bzcalc.GetCutPlane() * tl2::create<t_vec_bz>({ x, 0., plane_d });
				//std::cout << x << ": " << QinvA << std::endl;

				parser.register_var("x", x);
				parser.register_var("Qx", QinvA[0]);
				parser.register_var("Qy", QinvA[1]);

				t_real y = parser.eval();
				if(std::isnan(y) || std::isinf(y))
					continue;
				if(y < min_y || y > max_y)
					continue;

				curve.emplace_back(tl2::create<t_vec_bz>({ x, y }));
			}

			m_bzscene->AddCurve(curve);
		}
		catch(const std::exception& ex)
		{
			m_status->setText((std::string("<font color=\"red\">")
				+ "Formula " + tl2::var_to_str(formula_idx + 1) + ": "
				+ ex.what() + std::string("</font>")).c_str());
			all_ok = false;
		}
	}

	if(all_ok)
	{
		m_status->setText((std::string("<font color=\"green\">")
			+ "Formulas parsed successfully."
			+ std::string("</font>")).c_str());
	}

	return all_ok;
}


/**
 * calculate reciprocal coordinates of the cursor position
 */
void BZDlg::BZCutMouseMoved(t_real x, t_real y)
{
	t_real plane_d = m_cutD->value() * m_bzcalc.GetCutNormScale();

	t_vec_bz QinvA = m_bzcalc.GetCutPlane() * tl2::create<t_vec_bz>({ x, y, plane_d });
	t_mat_bz B_inv = m_crystA / (t_real(2)*tl2::pi<t_real>);
	t_vec_bz Qrlu = B_inv * QinvA;

	tl2::set_eps_0(QinvA, g_eps_bz);
	tl2::set_eps_0(Qrlu, g_eps_bz);

	std::ostringstream ostr;
	ostr.precision(g_prec_gui_bz);

	ostr << "Q = (" << QinvA[0] << ", " << QinvA[1] << ", " << QinvA[2] << ") Å⁻¹";
	ostr << " = (" << Qrlu[0] << ", " << Qrlu[1] << ", " << Qrlu[2] << ") rlu.";
	m_status->setText(ostr.str().c_str());
}
