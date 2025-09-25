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
 * export the magnetic structure to the sunny tool
 *   (https://github.com/SunnySuite/Sunny.jl)
 */
void MagDynDlg::ExportToSunny()
{
	QString dirLast = m_sett->value("dir_export_sun", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save As Jl File", dirLast, "jl files (*.jl)");
	if(filename == "")
		return;

	if(ExportToSunny(filename))
		m_sett->setValue("dir_export_sun", QFileInfo(filename).path());
}



/**
 * export the magnetic structure to the sunny tool
 *   (https://github.com/SunnySuite/Sunny.jl)
 */
bool MagDynDlg::ExportToSunny(const QString& _filename)
{
	std::string filename = _filename.toStdString();
	std::string dispname_abs = tl2::get_file_noext(filename) + ".dat";
	std::string dispname_rel = tl2::get_file_nodir(dispname_abs);

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

	ofstr	<< "#\n"
		<< "# Created by Magpie\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n\n";

	ofstr << "using Sunny\nusing Printf\n";


	// --------------------------------------------------------------------
	ofstr << "\n# options\n";
	ofstr << "calc_groundstate = false\n";
	ofstr << "plot_structure   = true\n";
	ofstr << "plot_dynamics    = true\n";
	ofstr << "save_dynamics    = true\n";
	ofstr << "datfile          = \"" << dispname_rel << "\"\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	t_real h1 = (t_real)m_Q_start[0]->value();
	t_real k1 = (t_real)m_Q_start[1]->value();
	t_real l1 = (t_real)m_Q_start[2]->value();
	t_real h2 = (t_real)m_Q_end[0]->value();
	t_real k2 = (t_real)m_Q_end[1]->value();
	t_real l2 = (t_real)m_Q_end[2]->value();

	ofstr << "\n# variables\n";

	// internal constants and variables
	ofstr << "g_e     = " << tl2::g_e<t_real> << "\n";
	ofstr << "Qstart  = [ " << h1 << ", " << k1 << ", " << l1 << " ]\n";
	ofstr << "Qend    = [ " << h2 << ", " << k2 << ", " << l2 << " ]\n";
	ofstr << "Qpts    = " << m_num_points->value() << "\n";

	// user variables
	for(const auto &var : m_dyn.GetVariables())
	{
		ofstr << var.name << " = " << var.value.real();
		if(!tl2::equals_0<t_real>(var.value.imag(), g_eps))
			ofstr << " + var.value.imag()" << "im";
		ofstr << "\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic sites and xtal lattice\n";
	ofstr << "@printf(\"Setting up magnetic sites...\\n\")\n";

	const auto& xtal = m_dyn.GetCrystalLattice();
	ofstr << "magsites = Crystal(\n"
		<< "\tlattice_vectors("
		<< xtal[0] << ", "
		<< xtal[1] << ", "
		<< xtal[2] << ", "
		<< tl2::r2d<t_real>(xtal[3]) << ", "
		<< tl2::r2d<t_real>(xtal[4]) << ", "
		<< tl2::r2d<t_real>(xtal[5]) << "),\n\t[\n";

	ofstr << "\t\t# site list\n";
	for(const t_site &site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t\t[ "
			<< get_str_var(site.pos[0]) << ", "
			<< get_str_var(site.pos[1]) << ", "
			<< get_str_var(site.pos[2]) << " ],"
			<< " # " << site.name << "\n";
	}

	// save as the P1 space group, as we have already performed the symmetry operations
	// (you can also manually set the crystal's space group and delete all
	//  symmetry-equivalent positions and couplings in the generated file)
	ofstr << "\t], 1)\n";

	ofstr << "num_sites = length(magsites.positions)\n\n";


	ofstr << "# spin magnitudes and magnetic system\n";
	ofstr << "magsys = System(magsites, #( 1, 1, 1 ),\n\t[\n";

	t_size site_idx = 1;
	for(const t_site& site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t\t" << site_idx << " => Moment("
			<< "s = " << get_str_var(site.spin_mag) << ", "
			<< "g = -[ g_e 0 0; 0 g_e 0; 0 0 g_e ]),"
			<< " # " << site.name << "\n";
		++site_idx;
	}
	ofstr << "\t], :dipole)\n\n";


	ofstr << "# spin directions\n";
	const auto& field = m_dyn.GetExternalField();
	if(field.align_spins)
	{
		// set all spins to field direction
		ofstr << "polarize_spins!(magsys, [ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ])\n";
	}
	else
	{
		// set individual spins
		site_idx = 1;
		for(const t_site& site : m_dyn.GetMagneticSites())
		{
			ofstr << "set_dipole!(magsys, [ "
				<< get_str_var(site.spin_dir[0]) << ", "
				<< get_str_var(site.spin_dir[1]) << ", "
				<< get_str_var(site.spin_dir[2]) << " ], "
				<< "( 1, 1, 1, " << site_idx << " ))"
				<< " # " << site.name << "\n";
			++site_idx;
		}
	}

	ofstr << "\n@printf(\"%s\", magsites)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";
	ofstr << "@printf(\"Setting up magnetic couplings...\\n\")\n";

	for(const t_term& term : m_dyn.GetExchangeTerms())
	{
		t_size idx1 = m_dyn.GetMagneticSiteIndex(term.site1) + 1;
		t_size idx2 = m_dyn.GetMagneticSiteIndex(term.site2) + 1;
		bool is_aniso = (idx1 == idx2 && tl2::equals_0(term.dist_calc, g_eps));

		if(is_aniso)
		{
			ofstr << "set_onsite_coupling!(magsys, S -> "
			    << get_str_var(term.J, true) << "*(S[1]^2 + S[2]^2 + S[3]^2), "
				<< idx1 << ");\n";

			// TODO: also treat dmi vector and general matrix
		}
		else
		{
			ofstr << "set_exchange!(magsys," << " # " << term.name
				<< "\n\t[\n"
				<< "\t\t" << get_str_var(term.J, true)             // 0,0
				<< "   " << get_str_var(term.dmi[2], true)         // 0,1
				<< "  -" << get_str_var(term.dmi[1], true) << ";"  // 0,2
				<< "\n\t\t-" << get_str_var(term.dmi[2], true)     // 1,0
				<< "   " << get_str_var(term.J, true)              // 1,1
				<< "   " << get_str_var(term.dmi[0], true) << ";"  // 1,2
				<< "\n\t\t" << get_str_var(term.dmi[1], true)      // 2,0
				<< "  -" << get_str_var(term.dmi[0], true)         // 2,1
				<< "   " << get_str_var(term.J, true)              // 2,2
				<< "\n\t]";

			if(!tl2::equals_0(term.Jgen_calc, g_eps))
			{
				ofstr << " +\n\t[\n"
					<< "\t\t" << get_str_var(term.Jgen[0][0], true)
					<< "  " << get_str_var(term.Jgen[0][1], true)
					<< "  " << get_str_var(term.Jgen[0][2], true) << ";"
					<< "\n\t\t" << get_str_var(term.Jgen[1][0], true)
					<< "  " << get_str_var(term.Jgen[1][1], true)
					<< "  " << get_str_var(term.Jgen[1][2], true) << ";"
					<< "\n\t\t" << get_str_var(term.Jgen[2][0], true)
					<< "  " << get_str_var(term.Jgen[2][1], true)
					<< "  " << get_str_var(term.Jgen[2][2], true)
					<< "\n\t]";
			}

			ofstr << ", Bond(" << idx1 << ", " << idx2 << ", [ "
				<< get_str_var(term.dist[0]) << ", "
				<< get_str_var(term.dist[1]) << ", "
				<< get_str_var(term.dist[2])
				<< " ]))\n";
		}
	}


	if(!tl2::equals_0<t_real>(field.mag, g_eps))
	{
		ofstr << "\n# external field\n";
		ofstr << "phys_units = Units(:meV, :angstrom)\n";
		ofstr << "set_field!(magsys, -[ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ] * " << field.mag
			<< " * phys_units.T"
			<< ")\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# optionally calculate the ground state\n";
	ofstr << "if calc_groundstate\n";
	ofstr << "\t@printf(\"Calculating ground state...\\n\")\n";
	ofstr << "\trandomize_spins!(magsys)\n";
	ofstr << "\tminimize_energy!(magsys; maxiters = 1024)\n";
	ofstr << "end\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# optionally plot nuclear and magnetic structure\n";
	ofstr << "if plot_structure\n";
	ofstr << "\t@printf(\"Plotting structure...\\n\")\n";
	ofstr << "\tusing GLMakie\n";
	ofstr << "\tview_crystal(magsys, refbonds = 15, compass = true)\n";
	ofstr << "\tplot_spins(magsys, compass = true)\n";
	ofstr << "end\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	if(m_dyn.IsIncommensurate())
	{
		ofstr << "\n# supercell for incommensurate structure\n";

		const t_vec_real& prop = m_dyn.GetOrderingWavevector();
		const t_vec_real& axis = m_dyn.GetRotationAxis();
		t_vec_real s0 = tl2::cross(prop, axis);
		s0 /= tl2::norm(s0);

		int sc_x = tl2::equals_0(prop[0], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[0]));
		int sc_y = tl2::equals_0(prop[1], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[1]));
		int sc_z = tl2::equals_0(prop[2], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[2]));

		ofstr << "magsys = reshape_supercell(magsys, [ "
			<< sc_x << " 0 0; 0 " << sc_y << " 0; 0 0 " << sc_z
			<< " ])\n";

		ofstr << "set_spiral_order!(magsys; "
			<< "k = [ " << prop[0] << ", " << prop[1] << ", " << prop[2] << " ], "
			<< "axis = [ " << axis[0] << ", " << axis[1] << ", " << axis[2] << " ], "
			<< "S0 = [ " << s0[0] << ", " << s0[1] << ", " << s0[2] << " ])\n";
	}

	ofstr << "\n@printf(\"%s\\n\", magsys)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# spin-wave calculation\n";
	ofstr << "@printf(\"Calculating S(Q, E)...\\n\")\n";

	ofstr << "cholesky_eps = 1e-8\n";
	std::string proj = m_use_projector->isChecked() ? "ssf_perp" : "ssf_trace";
	ofstr << "calc = SpinWaveTheory(magsys; measure = " << proj << "(magsys), regularization = cholesky_eps)\n";

	//ofstr << "momenta = collect(range(Qstart, Qend, Qpts))\n";
	ofstr << "momenta = q_space_path(magsys.crystal, [ Qstart, Qend ], Qpts)\n";
	ofstr << "bands = intensities_bands(calc, momenta)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# plot the dispersion\n";
	ofstr << "if plot_dynamics\n";
	ofstr << "\t@printf(\"Plotting dispersion...\\n\")\n";
	ofstr << "\tusing GLMakie\n";
	ofstr << "\tplot_intensities(bands; fwhm = 0.1, units = phys_units)\n";
	ofstr << "end\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# output the dispersion and spin-spin correlation\n";
	ofstr << "if save_dynamics\n";
	ofstr << "\t@printf(\"Outputting dispersion data to \\\"%s\\\", plot with (adapting x index):\\n"
		<< "\\tgnuplot -p -e \\\"plot \\\\\\\"%s\\\\\\\" u 1:4:(\\\\\\$5) w p pt 7 ps var\\\"\\n\", "
		<< "datfile, datfile)\n";

	ofstr << "\tenergies = bands.disp\n";
	ofstr << "\tcorrelations = bands.data\n";

	ofstr << "\topen(datfile, \"w\") do ostr\n";
	ofstr <<
		R"BLOCK(		@printf(ostr, "# %8s %10s %10s %10s %10s\n",
			"h (rlu)", "k (rlu)", "l (rlu)", "E (meV)", "S(Q, E)")
		for q_idx in 1:length(momenta.qs)
			for e_idx in 1:length(energies[:, q_idx])
				@printf(ostr, "%10.4f %10.4f %10.4f %10.4f %10.4f\n",
					momenta.qs[q_idx][1], momenta.qs[q_idx][2], momenta.qs[q_idx][3],
					energies[e_idx, q_idx],
					correlations[e_idx, q_idx] / num_sites)
			end
		end
	end
end
)BLOCK";
	// --------------------------------------------------------------------

	return true;
}

