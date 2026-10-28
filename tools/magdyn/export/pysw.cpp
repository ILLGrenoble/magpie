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
#include <boost/algorithm/string/replace.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;
namespace algo = boost::algorithm;

#include "../magdyn.h"

#include <QtCore/QString>

#include <fstream>
#include <sstream>
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


static std::string get_str_var_repl(const std::string& _var, const std::string& S_mag, bool add_brackets = false)
{
	std::string var = _var;
	if(S_mag != "")
		algo::replace_all(var, "_S_mag", S_mag);

	return get_str_var(var, add_brackets);
}



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
	ofstr << "import numpy as np\n";
	ofstr << "import time\n\n";


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


	ofstr << "# options\n";
	ofstr << "only_pos_E    = True     # ignore magnon annihilation?\n";
	ofstr << "print_results = False    # output final results\n";
	ofstr << "weight_scale  = 16.      # S(Q, E) scaling factor for plotting\n";
	ofstr << "\n";


	ofstr << "# variables\n";

	// internal constants and variables
	ofstr << "g_e     = " << std::abs(tl2::g_e<t_real>) << "\n";
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


	ofstr << "\ntime_begin = time.perf_counter()\n";


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
	const auto& field = m_dyn.GetExternalField();
	const t_vec3_real& field_dir = m_dyn.GetExternalFieldLab(false);
	const t_vec3_real& prop = m_dyn.GetOrderingWavevector();
	const t_vec3_real& axis = m_dyn.GetRotationAxis();

	ofstr << "\n# xtal lattice and magnetic sites\n";
	ofstr << "sites = sw.Structure(sites = [\n";

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
			<< "i = " << get_str_var(site.pos[0]) << ", "
			<< "j = " << get_str_var(site.pos[1]) << ", "  
			<< "k = " << get_str_var(site.pos[2]) << ", ";
		if(field.align_spins && field.dir)
		{
			ofstr
				<< "sx = " << -field_dir[0] << " * " << get_str_var(site.spin_mag) << ", "
				<< "sy = " << -field_dir[1] << " * " << get_str_var(site.spin_mag) << ", "
				<< "sz = " << -field_dir[2] << " * " << get_str_var(site.spin_mag) << ", ";
		}
		else
		{
			ofstr
				<< "sx = " << get_str_var(site.spin_dir[0]) << " * " << get_str_var(site.spin_mag) << ", "
				<< "sy = " << get_str_var(site.spin_dir[1]) << " * " << get_str_var(site.spin_mag) << ", "  
				<< "sz = " << get_str_var(site.spin_dir[2]) << " * " << get_str_var(site.spin_mag) << ", ";
		}
		ofstr << "g = g_e, ";
		ofstr << "name = \"" << site.name << "\"),\n";

		local_site_indices.insert(std::make_pair(site.name, local_site_idx));
		++local_site_idx;
	}

	ofstr << "\t]";

	if(m_dyn.IsIncommensurate())
	{
		ofstr << ",\n\tsupercell = sw.RotationSupercell("
			<< "[ " << axis[0] << ", " << axis[1] << ", " << axis[2] << " ], "
			<< "[ " << prop[0] << ", " << prop[1] << ", " << prop[2] << " ])";
	}

	ofstr << ",\n\tunit_cell = sw.UnitCell("
		<< xtal[0] << ", " << xtal[1] << ", " << xtal[2] << ", "
		<< tl2::r2d<t_real>(xtal[3]) << ", "
		<< tl2::r2d<t_real>(xtal[4]) << ", "
		<< tl2::r2d<t_real>(xtal[5]) << ")";

	ofstr << ", spacegroup = None";

	ofstr << ")\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";
	ofstr << "couplings = [\n";

	std::unordered_set<t_size> seen_term_sym_indices;
	for(const t_term& term : m_dyn.GetExchangeTerms())
	{
		const t_size idx1 = m_dyn.GetMagneticSiteIndex(term.site1);
		const t_size idx2 = m_dyn.GetMagneticSiteIndex(term.site2);
		const bool is_aniso = (idx1 == idx2 && tl2::equals_0(term.dist_calc, g_eps));
	
		std::string S1_mag;
		if(is_aniso)
			S1_mag = get_str_var(m_dyn.GetMagneticSites()[idx1].spin_mag);

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
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = [ "
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2]
				<< " ], j = " << get_str_var_repl(term.J, S1_mag) << "),"
				<< "\t\t# " << term.name << "\n";
		}


		if(term.dmi_calc && !tl2::equals_0(*term.dmi_calc, g_eps))
		{
			ofstr << "\tsw.DMExchange(sites.sites[" << local_site_indices[term.site1] << "], "
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = [ "
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2] << " ], "
				<< "d_x = " << get_str_var_repl(term.dmi[0], S1_mag) << ", "
				<< "d_y = " << get_str_var_repl(term.dmi[1], S1_mag) << ", "
				<< "d_z = " << get_str_var_repl(term.dmi[2], S1_mag) << "),"
				<< "\t\t# " << term.name << "\n";
		}

		if(term.Jgen_calc && !tl2::equals_0(*term.Jgen_calc, g_eps))
		{
			ofstr << "\tsw.exchange.Exchange(sites.sites[" << local_site_indices[term.site1] << "], "
				<< "sites.sites[" << local_site_indices[term.site2] << "], " << "cell_offset = [ "
				<< term.dist_calc[0]
				<< ", " << term.dist_calc[1]
				<< ", " << term.dist_calc[2]
				<< " ], exchange_matrix = np.array([[ "
				<< get_str_var_repl(term.Jgen[0][0], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[0][1], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[0][2], S1_mag) << " ], [ "
				<< get_str_var_repl(term.Jgen[1][0], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[1][1], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[1][2], S1_mag) << " ], [ "
				<< get_str_var_repl(term.Jgen[2][0], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[2][1], S1_mag) << ", "
				<< get_str_var_repl(term.Jgen[2][2], S1_mag) << " ]])),"
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
	}*/
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# spin-wave calculation\n";
	ofstr << "H = sw.Hamiltonian(sites, couplings, sia_couplings)\n";
	ofstr << "# H.print_summary()\n";

	ofstr << "Qs = sw.Path([ Qstart, Qend ], Qpts)\n";
	std::string Scomp = m_use_projector->isChecked() ? "\"Sperp\"" : "\"Sxx+Syy+Szz\"";

	std::ostringstream ostrField;
	ostrField.precision(g_prec);
	
	if(!tl2::equals_0<t_real>(field.mag, g_eps) && field.dir)
	{
		ostrField << ", field = np.array([ "
			<< field_dir[0] << ", "
			<< field_dir[1] << ", "
			<< field_dir[2] << " ]) * " << field.mag;
	}

	ofstr << "# H.spaghetti_plot(Qs, dE = 0.1, components = " << Scomp;
	ofstr << ostrField.str() << ")\n";
	ofstr << "[ allEs, allws ] = H.energies_and_intensities(Qs, components = " << Scomp;
	ofstr << ostrField.str() << ")\n";

	ofstr << R"BLOCK(
if print_results:
	print("%-15s %-15s %-15s %-15s %-15s" % ("# h (rlu)", "k (rlu)", "l (rlu)", "E (meV)", "S (a.u.)"))

hs, ks, ls, Es, ws = [], [], [], [], []
for Q, _Es, _ws in zip(Qs.q_points(), allEs, allws):
	for E, w in zip(_Es, _ws):
		E = E.real
		if only_pos_E and E < 0.:
			continue

		hs.append(Q[0])
		ks.append(Q[1])
		ls.append(Q[2])
		Es.append(E)
		ws.append(w * weight_scale)

		if print_results:
			print("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f" % (Q[0], Q[1], Q[2], E, w))

time_total = time.perf_counter() - time_begin
print("Calculation took %.6f seconds." % time_total)

# choose the Q component with the most change
qs = hs
if np.std(ks) > np.std(qs):
	qs = ks
if np.std(ls) > np.std(qs):
	qs = ls

# plotting
print("Plotting dispersion...")
import matplotlib.pyplot as plt

plt.plot()
plt.xlabel("q (rlu)")
plt.ylabel("E (meV)")
plt.scatter(qs, Es, marker = '.', s = ws)
plt.show()
)BLOCK";
	// --------------------------------------------------------------------


	return true;
}
