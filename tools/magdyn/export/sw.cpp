/**
 * magnetic dynamics -- exporting the magnetic structure to other magnon tools
 * @author Tobias Weber <tweber@ill.fr>
 * @date 26-june-2024
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

#include <boost/scope_exit.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "../magdyn.h"

#include <QtCore/QString>

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <cstdlib>

#include "tlibs2/libs/str.h"
#include "tlibs2/libs/file.h"
#include "tlibs2/libs/units.h"
#include "libs/symops.h"


// precision
extern int g_prec;
extern t_real g_eps;


extern std::string get_str_var(const std::string& var, bool add_brackets = false);



/**
 * export the magnetic structure to spinw
 *   (https://github.com/SpinW/spinw)
 */
void MagDynDlg::ExportToSpinW()
{
	QString dirLast = m_sett->value("dir_export_sw", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save As m File", dirLast, "m files (*.m)");
	if(filename == "")
		return;

	if(ExportToSpinW(filename))
		m_sett->setValue("dir_export_sw", QFileInfo(filename).path());
}



/**
 * export the magnetic structure to spinw
 *   (https://github.com/SpinW/spinw)
 */
bool MagDynDlg::ExportToSpinW(const QString& _filename)
{
	// make sure the symmetry indices are up-to-date
	CalcSymmetryIndices();
	std::string filename = _filename.toStdString();

	std::ofstream ofstr(filename);
	if(!ofstr)
	{
		ShowError("Cannot open file for writing.");
		return false;
	}

	ofstr.precision(g_prec);

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr	<< "%\n"
		<< "% Created by Magpie\n"
		<< "% URL: https://github.com/ILLGrenoble/magpie\n"
		<< "% DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "% User: " << user << "\n"
		<< "% Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "%\n\n";

	ofstr << "tic();\n";
	ofstr << "sw_obj = spinw();\n\n";


	// --------------------------------------------------------------------
	t_real h1 = (t_real)m_Q_start[0]->value();
	t_real k1 = (t_real)m_Q_start[1]->value();
	t_real l1 = (t_real)m_Q_start[2]->value();
	t_real h2 = (t_real)m_Q_end[0]->value();
	t_real k2 = (t_real)m_Q_end[1]->value();
	t_real l2 = (t_real)m_Q_end[2]->value();

	ofstr << "% variables\n";

	// internal constants and variables
	ofstr << "g_e     = " << tl2::g_e<t_real> << ";\n";
	ofstr << "Qstart  = [ " << h1 << " " << k1 << " " << l1 << " ];\n";
	ofstr << "Qend    = [ " << h2 << " " << k2 << " " << l2 << " ];\n";
	ofstr << "Qpts    = " << m_num_points->value() << ";\n";

	// user variables
	for(const auto &var : m_dyn.GetVariables())
	{
		ofstr << var.name << " = " << var.value.real();
		if(!tl2::equals_0<t_real>(var.value.imag(), g_eps))
			ofstr << " + var.value.imag()" << "j";
		ofstr << ";\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n% xtal lattice\n";

	const auto& symops = GetSymOpsForCurrentSG();

	ofstr << "symops = \'";
	for(t_size opidx = 0; opidx < symops.size(); ++opidx)
	{
		std::string symops_xyz =
			symop_to_xyz<t_mat_real, t_real>(symops[opidx], g_prec, g_eps);

		ofstr << symops_xyz;

		if(opidx < symops.size() - 1)
			ofstr << "; ";
	}
	ofstr << "\';\n";

	const auto& xtal = m_dyn.GetCrystalLattice();

	ofstr << "sw_obj.genlattice("
		<< "\"lat_const\", [ " << xtal[0] << " " << xtal[1] << " " << xtal[2] << " ], "
		<< "\"angled\", [ "
		<< tl2::r2d<t_real>(xtal[3]) << " "
		<< tl2::r2d<t_real>(xtal[4]) << " "
		<< tl2::r2d<t_real>(xtal[5]) << " ], "
		<< "\"sym\", symops);\n";

	ofstr << "\n% magnetic sites\n";

	std::unordered_set<t_size> seen_site_sym_indices;
	for(const t_site& site : m_dyn.GetMagneticSites())
	{
		// only emit one position of a symmetry group
		if(seen_site_sym_indices.find(site.sym_idx) != seen_site_sym_indices.end())
			continue;
		seen_site_sym_indices.insert(site.sym_idx);

		ofstr << "sw_obj.addatom(\"r\", "
			//<< "\"label\", \"" << site.name << "\", ";
			<< "[ "
			<< get_str_var(site.pos[0]) << " "
			<< get_str_var(site.pos[1]) << " "
			<< get_str_var(site.pos[2]) << " ]"
			<< ", \"S\", " << get_str_var(site.spin_mag)
			<< ");";
		ofstr << " % " << site.name << "\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n% spin directions\n";

	const auto& field = m_dyn.GetExternalField();
	const t_vec_real& prop = m_dyn.GetOrderingWavevector();
	const t_vec_real& axis = m_dyn.GetRotationAxis();

	ofstr << "spins = [ ";
	for(int i = 0; i < 3; ++i)
	{
		for(const t_site& site : m_dyn.GetMagneticSites())
		{
			if(field.align_spins)
				ofstr << field.dir[i] << " ";
			else
				ofstr << get_str_var(site.spin_dir[i]) << " ";
		}

		if(i < 2)
			ofstr << "; ";
	}
	ofstr << "];\n";

	std::string mode = m_dyn.IsIncommensurate() ? "\"helical\"" : "\"direct\"";
	ofstr << "sw_obj.genmagstr(\"mode\", " << mode << ", \"S\", spins, "
		<< "\"k\", [ " << prop[0] << " " << prop[1] << " " << prop[2] << " ], "
		<< "\"n\", [ " << axis[0] << " " << axis[1] << " " << axis[2] << " ]"
		<< ");\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n% magnetic couplings\n";
	ofstr << "sw_obj.gencoupling();\n";

	std::unordered_set<t_size> seen_term_sym_indices;
	for(const t_term& term : m_dyn.GetExchangeTerms())
	{
		t_size idx1 = m_dyn.GetMagneticSiteIndex(term.site1);
		t_size idx2 = m_dyn.GetMagneticSiteIndex(term.site2);
		bool is_aniso = (idx1 == idx2 && tl2::equals_0(term.dist_calc, g_eps));

		if(!is_aniso)
		{
			// only emit one coupling of a symmetry group
			if(seen_term_sym_indices.find(term.sym_idx) != seen_term_sym_indices.end())
				continue;
			seen_term_sym_indices.insert(term.sym_idx);
		}

		ofstr << "% " << term.name << "\n";

		if(!tl2::equals_0(term.J_calc, g_eps))
		{
			std::string J = "\'J_" + term.name + "\'";

			ofstr << "sw_obj.addmatrix(\"label\", "
				<< J << ", \"value\", "
				<< get_str_var(term.J)
				<< ");\n";

			if(is_aniso)
			{
				ofstr << "sw_obj.addaniso(" << J << ", " << idx1 + 1 << ");\n";
			}
			else
			{
				ofstr << "sw_obj.addcoupling(\"mat\", " << J << ", \"bond\", "
					<< term.sym_idx << ");\n";
			}
		}

		if(!tl2::equals_0(term.dmi_calc, g_eps))
		{
			std::string DMI = "\'DMI_" + term.name + "\'";

			ofstr << "sw_obj.addmatrix(\"label\", "
				<< DMI << ", \"value\", [ "
				<< get_str_var(term.dmi[0]) << " "
				<< get_str_var(term.dmi[1]) << " "
				<< get_str_var(term.dmi[2])
				<< " ]);\n";

			if(is_aniso)
			{
				ofstr << "sw_obj.addaniso(" << DMI << ", " << idx1 + 1 << ");\n";
			}
			else
			{
				ofstr << "sw_obj.addcoupling(\"mat\", " << DMI << ", \"bond\", "
					<< term.sym_idx << ");\n";
			}
		}

		if(!tl2::equals_0(term.Jgen_calc, g_eps))
		{
			std::string GEN = "\'GEN_" + term.name + "\'";

			ofstr << "sw_obj.addmatrix(\"label\", "
				<< GEN << ", \"value\", [ "
				<< get_str_var(term.Jgen[0][0]) << " "
				<< get_str_var(term.Jgen[0][1]) << " "
				<< get_str_var(term.Jgen[0][2]) << "; "
				<< get_str_var(term.Jgen[1][0]) << " "
				<< get_str_var(term.Jgen[1][1]) << " "
				<< get_str_var(term.Jgen[1][2]) << "; "
				<< get_str_var(term.Jgen[2][0]) << " "
				<< get_str_var(term.Jgen[2][1]) << " "
				<< get_str_var(term.Jgen[2][2]) << " "
				<< " ]);\n";

			if(is_aniso)
			{
				//std::string site_name = "\'" + m_dyn.GetMagneticSite(idx1).name + "\'";
				ofstr << "sw_obj.addaniso(" << GEN << ", " << /*site_name*/ idx1 + 1 << ");\n";
			}
			else
			{
				ofstr << "sw_obj.addcoupling(\"mat\", " << GEN << ", \"bond\", "
					<< term.sym_idx << ");\n";
			}
		}
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	if(m_dyn.GetTemperature() >= 0.)
	{
		ofstr << "\nsw_obj.temperature(" << m_dyn.GetTemperature() << ");\n";
	}

	if(!tl2::equals_0<t_real>(field.mag, g_eps))
	{
		ofstr << "\nsw_obj.field([ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ] * " << field.mag
			<< ");\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n% spin-wave calculation\n";

	std::string Scomp = m_use_projector->isChecked() ? "\'Sperp\'" : "\'Sxx+Syy+Szz\'";
	ofstr << "calc = sw_neutron(sw_obj.spinwave({ Qstart, Qend, Qpts }, \"hermit\", false));\n";
	ofstr << "bins = sw_egrid(calc, \"component\", " << Scomp << ");\n";
	ofstr << "toc();\n";

	ofstr << "\n% plotting\n";
	ofstr << "figure();\n";
	ofstr << "sw_plotspec(bins, \"mode\", 3, \"dE\", 0.1);\n";
	// --------------------------------------------------------------------

	return true;
}

