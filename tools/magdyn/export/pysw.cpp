/**
 * magnetic dynamics -- exporting the magnetic structure to other magnon tools
 * @author Tobias Weber <tweber@ill.fr>
 * @date 25-oct-2026
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

#include <boost/scope_exit.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "../magdyn.h"

#include <QtCore/QString>

#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <cstdlib>

#include "tlibs2/libs/str.h"
#include "tlibs2/libs/file.h"
#include "tlibs2/libs/units.h"

#include "libs/sym/symops.h"


// precision
extern int g_prec;
extern t_real g_eps;


extern std::string get_str_var(const std::string& var, bool add_brackets = false);



/**
 * export the magnetic structure to pyspinw
 *   (https://github.com/SpinW/pySpinW)
 */
void MagDynDlg::ExportToPySpinW()
{
	QString dirLast = m_sett->value("dir_export_pysw", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save As py File", dirLast, "py files (*.py)");
	if(filename == "")
		return;

	if(ExportToPySpinW(filename))
		m_sett->setValue("dir_export_pysw", QFileInfo(filename).path());
}



/**
 * export the magnetic structure to pyspinw
 *   (https://github.com/SpinW/pySpinW)
 */
bool MagDynDlg::ExportToPySpinW(const QString& _filename)
{
	// make sure the symmetry indices are up-to-date
	CalcSymmetryIndices();
	std::string filename = _filename.toStdString();

	std::ofstream ofstr(filename);
	if(!ofstr)
	{
		ShowError(QString("Cannot open file \"%1\" for writing.").arg(filename.c_str()));
		return false;
	}

	ofstr.precision(g_prec);

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr	<< "#\n"
		<< "# Created by Magpie " << MAGPIE_VER << "\n"
		<< "# Author: Tobias Weber\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n\n";

	ofstr << "import pyspinw as sw\n";
	ofstr << "import time\n\n";

	ofstr << "time_begin = time.perf_counter()\n\n";


	// --------------------------------------------------------------------
	t_real h1 = (t_real)m_Q_start[0]->value();
	t_real k1 = (t_real)m_Q_start[1]->value();
	t_real l1 = (t_real)m_Q_start[2]->value();
	t_real h2 = (t_real)m_Q_end[0]->value();
	t_real k2 = (t_real)m_Q_end[1]->value();
	t_real l2 = (t_real)m_Q_end[2]->value();
	t_real peak1x = (t_real)m_scatteringplane[0]->value();
	t_real peak1y = (t_real)m_scatteringplane[1]->value();
	t_real peak1z = (t_real)m_scatteringplane[2]->value();
	t_real peak2x = (t_real)m_scatteringplane[3]->value();
	t_real peak2y = (t_real)m_scatteringplane[4]->value();
	t_real peak2z = (t_real)m_scatteringplane[5]->value();

	ofstr << "# variables\n";

	// internal constants and variables
	ofstr << "g_e     = " << tl2::g_e<t_real> << "\n";
	ofstr << "Qstart  = [ " << h1 << ", " << k1 << ", " << l1 << " ]\n";
	ofstr << "Qend    = [ " << h2 << ", " << k2 << ", " << l2 << " ]\n";
	ofstr << "Qpts    = " << m_num_points->value() << "\n";
	ofstr << "plane1  = [ " << peak1x << ", " << peak1y << ", " << peak1z << " ]\n";
	ofstr << "plane2  = [ " << peak2x << ", " << peak2y << ", " << peak2z << " ]\n";

	// user variables
	for(const auto &var : m_dyn.GetVariables())
	{
		ofstr << var.name << " = " << var.value.real();
		if(!tl2::equals_0<t_real>(var.value.imag(), g_eps))
			ofstr << " + var.value.imag()" << "j";
		ofstr << "\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	/*const auto& symops = GetSymOpsForCurrentSG();

	ofstr << "symops = \'";
	for(t_size opidx = 0; opidx < symops.size(); ++opidx)
	{
		std::string symops_xyz =
			sym::symop_to_xyz<t_mat44_real, t_real>(symops[opidx], g_prec, g_eps);

		ofstr << symops_xyz;

		if(opidx < symops.size() - 1)
			ofstr << "; ";
	}
	ofstr << "\'\n";*/

	bool only_unique = false;  // only emit symmetry-unique sites and couplings


	const auto& xtal = m_dyn.GetCrystalLattice();
	ofstr << "\n# xtal lattice and magnetic sites\n";
	ofstr << "sites = sw.Structure([\n";

	std::unordered_set<t_size> seen_site_sym_indices;
	std::unordered_map<std::string, t_size> local_site_indices;
	t_size local_site_idx = 0;
	for(const t_site& site : m_dyn.GetMagneticSites())
	{
		if(only_unique)
		{
			// only emit one position of a symmetry group
			if(seen_site_sym_indices.find(site.sym_idx) != seen_site_sym_indices.end())
				continue;
			seen_site_sym_indices.insert(site.sym_idx);
		}

		ofstr << "\tsw.Site("
			<< get_str_var(site.pos[0]) << ", "
			<< get_str_var(site.pos[1]) << ", "  
			<< get_str_var(site.pos[2]) << ", "
			<< get_str_var(site.spin_dir[0]) << " * " << get_str_var(site.spin_mag) << ", "
			<< get_str_var(site.spin_dir[1]) << " * " << get_str_var(site.spin_mag) << ", "  
			<< get_str_var(site.spin_dir[2]) << " * " << get_str_var(site.spin_mag) << ", "
			<< "name = \"" << site.name << "\"),\n";

		local_site_indices.insert(std::make_pair(site.name, local_site_idx));
		++local_site_idx;
	}

	ofstr << "\t], sw.UnitCell("
		<< xtal[0] << ", " << xtal[1] << ", " << xtal[2] << ", "
		<< tl2::r2d<t_real>(xtal[3]) << ", "
		<< tl2::r2d<t_real>(xtal[4]) << ", "
		<< tl2::r2d<t_real>(xtal[5]) << ")\n";

	ofstr << ")\n";
	// --------------------------------------------------------------------


	/*// --------------------------------------------------------------------
	// TODO
	ofstr << "\n# spin directions\n";

	const auto& field = m_dyn.GetExternalField();
	const t_vec3_real& field_dir = m_dyn.GetExternalFieldLab(false);
	const t_vec3_real& prop = m_dyn.GetOrderingWavevector();
	const t_vec3_real& axis = m_dyn.GetRotationAxis();

	ofstr << "spins = [ ";
	for(int i = 0; i < 3; ++i)
	{
		for(const t_site& site : m_dyn.GetMagneticSites())
		{
			if(field.align_spins && field.dir)
				ofstr << field_dir[i] << " ";
			else
				ofstr << get_str_var(site.spin_dir[i]) << " ";
		}

		if(i < 2)
			ofstr << "; ";
	}
	ofstr << "]\n";

	std::string mode = m_dyn.IsIncommensurate() ? "\"helical\"" : "\"direct\"";
	ofstr << "sw_obj.genmagstr(\"mode\", " << mode << ", \"S\", spins, "
		<< "\"k\", [ " << prop[0] << " " << prop[1] << " " << prop[2] << " ], "
		<< "\"n\", [ " << axis[0] << " " << axis[1] << " " << axis[2] << " ]"
		<< ")\n";*/
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";
	ofstr << "couplings = [\n";

	std::unordered_set<t_size> seen_term_sym_indices;
	for(const t_term& term : m_dyn.GetExchangeTerms())
	{
		t_size idx1 = m_dyn.GetMagneticSiteIndex(term.site1);
		t_size idx2 = m_dyn.GetMagneticSiteIndex(term.site2);
		bool is_aniso = (idx1 == idx2 && tl2::equals_0(term.dist_calc, g_eps));

		if(is_aniso)
			continue;

		if(only_unique)
		{
			// only emit one coupling of a symmetry group
			if(seen_term_sym_indices.find(term.sym_idx) != seen_term_sym_indices.end())
				continue;
			seen_term_sym_indices.insert(term.sym_idx);
		}

		if(!tl2::equals_0(term.J_calc, g_eps))
		{
			ofstr << "\tsw.HeisenbergExchange(sites.sites[" << local_site_indices[term.site1] << "], "
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = ("
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2]
				<< "), j = " << get_str_var(term.J) << "),"
				<< "\t\t# " << term.name << "\n";
		}


		if(term.dmi_calc && !tl2::equals_0(*term.dmi_calc, g_eps))
		{
			ofstr << "\tsw.DMExchange(sites.sites[" << local_site_indices[term.site1] << "], "
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = ("
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2] << "), "
				<< "d_x = " << get_str_var(term.dmi[0]) << ", "
				<< "d_y = " << get_str_var(term.dmi[1]) << ", "
				<< "d_z = " << get_str_var(term.dmi[2]) << "),"
				<< "\t\t# " << term.name << "\n";
		}

		if(term.Jgen_calc && !tl2::equals_0(*term.Jgen_calc, g_eps))
		{
			ofstr << "\tsw.Exchange(sites.sites[" << local_site_indices[term.site1] << "], "
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = ("
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2]
				<< "), exchange_matrix = [[ "
				<< get_str_var(term.Jgen[0][0]) << ", "
				<< get_str_var(term.Jgen[0][1]) << ", "
				<< get_str_var(term.Jgen[0][2]) << " ], [ "
				<< get_str_var(term.Jgen[1][0]) << ", "
				<< get_str_var(term.Jgen[1][1]) << ", "
				<< get_str_var(term.Jgen[1][2]) << " ], [ "
				<< get_str_var(term.Jgen[2][0]) << ", "
				<< get_str_var(term.Jgen[2][1]) << ", "
				<< get_str_var(term.Jgen[2][2]) << " ]]),"
				<< "\t\t# " << term.name << "\n";
		}
	}

	ofstr << "]\n\n";


	// TODO
	ofstr << "sia_couplings = []\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// TODO
	/*if(m_dyn.GetTemperature() >= 0.)
	{
	}

	if(!tl2::equals_0<t_real>(field.mag, g_eps) && field.dir)
	{
	}*/
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# spin-wave calculation\n";
	ofstr << "H = sw.Hamiltonian(sites, couplings, sia_couplings)\n";

	std::string Scomp = m_use_projector->isChecked() ? "\"Sperp\"" : "\"Sxx+Syy+Szz\"";
	ofstr << "H.spaghetti_plot("
		<< "sw.Path([ Qstart, Qend ], Qpts), dE = 0.1, "
		<< "components = " << Scomp << ")\n";
	// --------------------------------------------------------------------


	ofstr << "\ntime_total = time.perf_counter() - time_begin\n";
	ofstr << "print(\"Calculation took %.6f seconds.\" % time_total)\n";

	return true;
}
