/**
 * magnetic dynamics -- calculations for sites and coupling terms
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "magdyn.h"

#include <iostream>
#include <boost/scope_exit.hpp>

using t_numitem = tl2::NumericTableWidgetItem<t_real>;



/**
 * flip the coordinates of the magnetic site positions
 * (e.g. to get the negative phase factor for weights)
 */
void MagDynDlg::MirrorAtoms()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// iterate the magnetic sites
	for(int row=0; row<m_sitestab->rowCount(); ++row)
	{
		auto *pos_x = static_cast<t_numitem*>(m_sitestab->item(row, COL_SITE_POS_X));
		auto *pos_y = static_cast<t_numitem*>(m_sitestab->item(row, COL_SITE_POS_Y));
		auto *pos_z = static_cast<t_numitem*>(m_sitestab->item(row, COL_SITE_POS_Z));

		if(!pos_x || !pos_y || !pos_z)
		{
			std::cerr << "Invalid entry in sites table row "
				<< row << "." << std::endl;
			continue;
		}

		pos_x->SetValue(-pos_x->GetValue());
		pos_y->SetValue(-pos_y->GetValue());
		pos_z->SetValue(-pos_z->GetValue());
	}
}



/**
 * rotate the direction of the magnetic field
 */
void MagDynDlg::RotateField(bool ccw)
{
	t_vec_real axis = tl2::create<t_vec_real>(
	{
		(t_real)m_rot_axis[0]->value(),
		(t_real)m_rot_axis[1]->value(),
		(t_real)m_rot_axis[2]->value(),
	});

	t_vec_real B = tl2::create<t_vec_real>(
	{
		(t_real)m_field_dir[0]->value(),
		(t_real)m_field_dir[1]->value(),
		(t_real)m_field_dir[2]->value(),
	});

	t_real angle = tl2::d2r<t_real>(m_rot_angle->value());
	if(!ccw)
		angle = -angle;

	t_mat_real R = tl2::rotation<t_mat_real, t_vec_real>(
		axis, angle, false);
	B = R*B;
	tl2::set_eps_0(B, g_eps);

	for(int i = 0; i < 3; ++i)
	{
		m_field_dir[i]->blockSignals(true);
		m_field_dir[i]->setValue(B[i]);
		m_field_dir[i]->blockSignals(false);
	}

	if(m_autocalc->isChecked())
		CalcAll();
};



/**
 * set selected field as current
 */
void MagDynDlg::SetCurrentField()
{
	if(m_fields_cursor_row < 0 || m_fields_cursor_row >= m_fieldstab->rowCount())
		return;

	const auto* Bh = static_cast<t_numitem*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_H));
	const auto* Bk = static_cast<t_numitem*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_K));
	const auto* Bl = static_cast<t_numitem*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_L));
	const auto* Bmag = static_cast<t_numitem*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_MAG));

	if(!Bh || !Bk || !Bl || !Bmag)
		return;

	m_ignoreCalc = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	m_field_dir[0]->setValue(Bh->GetValue());
	m_field_dir[1]->setValue(Bk->GetValue());
	m_field_dir[2]->setValue(Bl->GetValue());
	m_field_mag->setValue(Bmag->GetValue());
}



/**
 * generate magnetic sites form the space group symmetries
 */
void MagDynDlg::GenerateSitesFromSG()
{
	try
	{
		const auto& symops = GetSymOpsForCurrentSG();

		SyncToKernel();
		m_dyn.SymmetriseMagneticSites(symops);
		SyncSitesFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		ShowError(ex.what());
	}
}



/**
 * generate exchange terms from space group symmetries
 */
void MagDynDlg::GenerateCouplingsFromSG()
{
	try
	{
		const auto& symops = GetSymOpsForCurrentSG();

		SyncToKernel();
		m_dyn.SymmetriseExchangeTerms(symops);
		SyncTermsFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		ShowError(ex.what());
	}
}



/**
 * extend the unit cell by copying the existing elements
 */
void MagDynDlg::ExtendStructure()
{
	try
	{
		t_size x_size = m_extCell[0]->value();
		t_size y_size = m_extCell[1]->value();
		t_size z_size = m_extCell[2]->value();

		bool remove_duplicates = true;
		bool flip_spin = false;

		SyncToKernel();
		m_dyn.ExtendStructure(x_size, y_size, z_size,
			remove_duplicates, flip_spin);
		SyncSitesFromKernel();
		SyncTermsFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		ShowError(ex.what());
	}
}



/**
 * generate possible couplings up to a certain distance
 */
void MagDynDlg::GeneratePossibleCouplings()
{
	try
	{
		const auto& symops = GetSymOpsForCurrentSG();

		t_real dist_max = m_maxdist->value();
		t_size sc_max = m_maxSC->value();
		t_size couplings_max = m_maxcouplings->value();

		SyncToKernel();
		m_dyn.GeneratePossibleExchangeTerms(dist_max, sc_max, couplings_max);
		m_dyn.CalcSymmetryIndices(symops);
		SyncTermsFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		ShowError(ex.what());
	}
}



/**
 * get the symmetry operators for the currently selected space group
 */
const std::vector<t_mat_real>& MagDynDlg::GetSymOpsForCurrentSG(bool show_err) const
{
	// current space group index
	int sgidx = m_comboSG->itemData(m_comboSG->currentIndex()).toInt();

	if(sgidx < 0 || t_size(sgidx) >= m_SGops.size())
	{
		if(show_err)
			ShowError("Invalid space group selected.");

		// return empty symop list
		static const std::vector<t_mat_real> nullvec{};
		return nullvec;
	}

	return m_SGops[sgidx];
}



/**
 * open the table import dialog
 */
void MagDynDlg::ShowTableImporter()
{
	if(!m_table_import_dlg)
	{
		m_table_import_dlg = new TableImportDlg(this, m_sett);

		connect(m_table_import_dlg, &TableImportDlg::SetAtomsSignal,
			this, &MagDynDlg::ImportAtoms);
		connect(m_table_import_dlg, &TableImportDlg::SetCouplingsSignal,
			this, &MagDynDlg::ImportCouplings);
	}

	m_table_import_dlg->show();
	m_table_import_dlg->raise();
	m_table_import_dlg->activateWindow();
}



/**
 * import magnetic site positions from table dialog
 */
void MagDynDlg::ImportAtoms(const std::vector<TableImportAtom>& atompos_vec,
	bool clear_existing)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// remove existing sites
	if(clear_existing)
		DelTabItem(m_sitestab, -1);

	for(const TableImportAtom& atompos : atompos_vec)
	{
		std::string S[3] = { "0", "0", "1" };
		std::string Smag = "1";
		t_size sym_idx = 0;  // TODO

		for(std::size_t i = 0; i < 3; ++i)
		{
			if(atompos.S[i] != "")
				S[i] = atompos.S[i];
			if(atompos.Smag != "")
				Smag = atompos.Smag;
		}

		AddSiteTabItem(-1, atompos.name, sym_idx,
			atompos.pos[0], atompos.pos[1], atompos.pos[2],
			S[0], S[1], S[2], Smag);
	}
}



/**
 * import magnetic couplings from table dialog
 */
void MagDynDlg::ImportCouplings(const std::vector<TableImportCoupling>& couplings,
	bool clear_existing)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// remove existing couplings
	if(clear_existing)
		DelTabItem(m_termstab, -1);

	for(const TableImportCoupling& coupling : couplings)
	{
		t_size atom_1 = 0, atom_2 = 0;
		t_size sym_idx = 0;  // TODO

		if(coupling.atomidx1) atom_1 = *coupling.atomidx1;
		if(coupling.atomidx2) atom_2 = *coupling.atomidx2;

		std::string atom_1_name = tl2::var_to_str(atom_1);
		std::string atom_2_name = tl2::var_to_str(atom_2);

		// get the site names from the table
		if((int)atom_1 < m_sitestab->rowCount())
		{
			if(auto *name = m_sitestab->item(atom_1, COL_SITE_NAME); name)
				atom_1_name = name->text().toStdString();
		}
		if((int)atom_2 < m_sitestab->rowCount())
		{
			if(auto *name = m_sitestab->item(atom_2, COL_SITE_NAME); name)
				atom_2_name = name->text().toStdString();
		}

		AddTermTabItem(-1, coupling.name, sym_idx,
			atom_1_name, atom_2_name,
			coupling.d[0], coupling.d[1], coupling.d[2], coupling.J,
			coupling.dmi[0], coupling.dmi[1], coupling.dmi[2],
			coupling.Jgen[0], coupling.Jgen[1], coupling.Jgen[2],
			coupling.Jgen[3], coupling.Jgen[4], coupling.Jgen[5],
			coupling.Jgen[6], coupling.Jgen[7], coupling.Jgen[8]);
	}
}



/**
 * assign symmetry groups to sites and couplings
 */
void MagDynDlg::CalcSymmetryIndices()
{
	m_dyn.CalcSymmetryIndices(GetSymOpsForCurrentSG());
	SyncSymmetryIndicesFromKernel();
}



/**
 * sort couplings by their lengths
 */
void MagDynDlg::SortTerms()
{
	m_dyn.SortExchangeTerms();
	SyncTermsFromKernel();
}



/**
 * calculate brillouin zone and cut
 */
void MagDynDlg::CalcBZ()
{
	if(m_ignoreCalc)
		return;

	const t_mat_real& xtalA = m_dyn.GetCrystalATrafo();
	const t_mat_real& xtalB = m_dyn.GetCrystalBTrafo();

	m_bz.SetEps(g_eps);
	m_bz.SetSymOps(GetSymOpsForCurrentSG(), false);
	m_bz.SetCrystalA(xtalA);
	m_bz.SetCrystalB(xtalB);

	m_bz.CalcPeaks(g_bz_calc_order, false /*invA*/, false /*cut*/);
	m_bz.CalcPeaks(g_bz_draw_order, false /*invA*/, true /*cut*/);
	m_bz.CalcPeaksInvA();

	// get plane coordinate system
	const t_vec_real* plane = m_dyn.GetScatteringPlane();
	t_real plane_d = 0.;

	// calculate brillouin zone
	if(!m_bz.CalcBZ())
	{
		std::cerr << "Error calculating Brillouin zone." << std::endl;
		return;
	}

	// calculate brillouin zone cut
	if(!m_bz.CalcBZCut(plane[0], plane[2], plane_d, true))
	{
		std::cerr << "Error calculating Brillouin zone cut." << std::endl;
		return;
	}

	// draw brillouin zone
	if(m_bz_dlg)
	{
		m_bz_dlg->SetABTrafo(xtalA, xtalB);
		m_bz_dlg->SetEps(g_eps);
		m_bz_dlg->SetPrecGui(g_prec_gui);

		// clear old plot
		m_bz_dlg->Clear();

		// add gamma point
		std::size_t idx000 = m_bz.Get000Peak();
		const std::vector<t_vec_bz>& Qs_invA = m_bz.GetPeaks(true);
		if(idx000 < Qs_invA.size())
			m_bz_dlg->AddBraggPeak(Qs_invA[idx000]);

		// add voronoi vertices forming the vertices of the BZ
		for(const t_vec_bz& voro : m_bz.GetVertices())
			m_bz_dlg->AddVoronoiVertex(voro);

		// add voronoi bisectors
		m_bz_dlg->AddTriangles(m_bz.GetAllTriangles(), &m_bz.GetAllTrianglesFaceIndices());

		// scattering plane
		m_bz_dlg->SetPlane(
			tl2::col<t_mat_real, t_vec_real>(m_bz.GetCutPlane(), 2),  // normal
			m_bz.GetCutPlaneD());                                     // distance, here: 0
	}

	// draw brillouin zone cut
	if(m_bzscene)
	{
		m_bzscene->SetEps(g_eps);
		m_bzscene->SetPrecGui(g_prec_gui);
		m_bzscene->ClearAll();
		m_bzscene->AddCut(m_bz.GetCutLines(false));
		m_bzscene->AddPeaks(m_bz.GetPeaksOnPlane(true), &m_bz.GetPeaksOnPlane(false));
		m_bzview->Centre();
	}
}



/**
 * calculate reciprocal coordinates of the cursor position
 */
std::pair<t_vec_real, t_vec_real> MagDynDlg::GetBZCutQ(t_real x, t_real y) const
{
	t_real plane_d = 0. * m_bz.GetCutNormScale();  // plane always goes through Gamma point

	t_vec_real QinvA = m_bz.GetCutPlane() * tl2::create<t_vec_real>({ x, y, plane_d });
	t_mat_real B_inv = tl2::trans(m_dyn.GetCrystalATrafo()) / (t_real(2)*tl2::pi<t_real>);
	//auto [B_inv, ok] = tl2::inv(m_dyn.GetCrystalBTrafo());
	//t_mat_real UB_inv = m_dyn.GetCrystalUBTrafo(true);
	t_vec_real Qrlu = B_inv * QinvA;

	tl2::set_eps_0(QinvA, g_eps);
	tl2::set_eps_0(Qrlu, g_eps);

	return std::make_pair(QinvA, Qrlu);
}



/**
 * mouse cursor moved
 */
void MagDynDlg::BZCutMouseMoved(t_real x, t_real y)
{
	auto [QinvA, Qrlu] = GetBZCutQ(x, y);

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);

	ostr << "Q = (" << QinvA[0] << ", " << QinvA[1] << ", " << QinvA[2] << ") Å⁻¹";
	ostr << " = (" << Qrlu[0] << ", " << Qrlu[1] << ", " << Qrlu[2] << ") rlu.";
	m_status->setText(ostr.str().c_str());
}



/**
 * mouse clicked
 */
void MagDynDlg::BZCutMouseClicked(int buttons, t_real x, t_real y)
{
	if(buttons & 4)  // right button
	{
		auto [QinvA, Qrlu] = GetBZCutQ(x, y);
		m_bz_cur_pos = std::move(Qrlu);

		QPointF _pt{x*m_bzscene->GetScale(), y*m_bzscene->GetScale()};
		QPoint pt = m_bzview->mapToGlobal(m_bzview->mapFromScene(_pt));

		m_bzcontext->popup(pt);
	}
}
