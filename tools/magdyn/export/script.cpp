/**
 * magnetic dynamics -- exporting the magnetic structure to other magnon tools
 * @author Tobias Weber <tweber@ill.fr>
 * @date 26-june-2024, 3-april-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (McClarty 2022) https://doi.org/10.1146/annurev-conmatphys-031620-104715
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
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
#include <cstdlib>

#include "tlibs2/libs/str.h"
#include "tlibs2/libs/file.h"
#include "tlibs2/libs/units.h"
#include "libs/symops.h"


// precision
extern int g_prec;
extern t_real g_eps;


std::string get_str_var(const std::string& var, bool add_brackets = false)
{
	if(var == "")
	{
		return "0";
	}
	else
	{
		if(add_brackets)
			return "(" + var + ")";
		return var;
	}
}



/*
 * export the magnetic structure to a self-contained script
 */
void MagDynDlg::ExportToScript()
{
	if(m_dyn.IsIncommensurate())
	{
		ShowError("Incommensurate structure export is not yet supported.", true);
		return;
	}

	QString dirLast = m_sett->value("dir_export_script", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save As Py File", dirLast, "py files (*.py)");
	if(filename == "")
		return;

	if(ExportToScript(filename))
		m_sett->setValue("dir_export_script", QFileInfo(filename).path());
}



/**
 * export the magnetic structure to a self-contained script
 */
bool MagDynDlg::ExportToScript(const QString& _filename)
{
	std::string scr = R"BLOCK(
only_pos_E    = True     # hide magnon annihilation?
verbose_print = False    # print intermediate results
weight_scale  = 16.      # S(q, E) scaling factor for plotting


# debug output
def print_infos(str):
	if verbose_print:
		print(str)

# skew-symmetric (cross-product) matrix
def skew(vec):
	return np.array([
		[      0.,   vec[2],  -vec[1] ],
		[ -vec[2],       0.,   vec[0] ],
		[  vec[1],  -vec[0],       0. ] ])

# calculate structure properties
def init(sites, couplings):
	# calculate spin rotations towards ferromagnetic order along [001]
	for site in sites:
		zdir = np.array([ 0., 0., 1. ])
		Sdir = np.array(site["Sdir"]) / la.norm(site["Sdir"])
		rotaxis = np.array([ 0., 1., 0. ])
		s = 0.

		if np.allclose(Sdir, zdir):
			# spin and z axis parallel
			c = 1.
		elif np.allclose(Sdir, -zdir):
			# spin and z axis anti-parallel
			c = -1.
		else:
			# sine and cosine of the angle between spin and z axis
			rotaxis = np.cross(Sdir, zdir)
			s = la.norm(rotaxis)
			c = np.dot(Sdir, zdir)
			rotaxis /= s

		# rotation via rodrigues' formula, see (Arens 2015), p. 718 and p. 816
		rot = (1. - c) * np.outer(rotaxis, rotaxis) + np.diag([ c, c, c ]) - skew(rotaxis)*s
		site["u"] = rot[0, :] + 1j * rot[1, :]
		site["v"] = rot[2, :]

		print_infos(np.dot(rot, Sdir))
		print_infos("\nrot = \n%s\nu = %s\nv = %s" % (rot, site["u"], site["v"]))

	# calculate real interaction matrices
	for coupling in couplings:
		J = coupling["J"]
		coupling["J_real"] = np.diag([ J, J, J ]) + skew(coupling["DMI"])
		print_infos("\nJ_real =\n%s" % coupling["J_real"])

# get the energies of the dispersion at the momentum transfer Qvec
def get_energies(Qvec, sites, couplings, field = None):
	print_infos("\n\nQ = %s" % Qvec)

	# fourier transform interaction matrices
	num_sites = len(sites)
	J_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)
	J0_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)

	for coupling in couplings:
		dist = np.array(coupling["dist"])
		J_real = coupling["J_real"]
		site1 = coupling["sites"][0]
		site2 = coupling["sites"][1]

		J_ft = J_real * np.exp(-2j*np.pi * np.dot(dist, Qvec))
		J_fourier[site1, site2] += J_ft
		J_fourier[site2, site1] += J_ft.transpose().conj()
		J0_fourier[site1, site2] += J_real
		J0_fourier[site2, site1] += J_real.transpose().conj()

	print_infos("\nJ_fourier =\n%s\n\nJ0_fourier =\n%s" % (J_fourier, J0_fourier))

	# hamiltonian
	H = np.zeros((2*num_sites, 2*num_sites), dtype = complex)

	for i in range(num_sites):
		S_i = sites[i]["S"]
		u_i = sites[i]["u"]
		v_i = sites[i]["v"]

		for j in range(num_sites):
			S_j = sites[j]["S"]
			u_j = sites[j]["u"]
			v_j = sites[j]["v"]
			S = 0.5 * np.sqrt(S_i * S_j)

			H[            i,             j] += S   * np.dot(u_i,        np.dot(J_fourier[i, j],  u_j.conj()))
			H[            i,             i] -= S_j * np.dot(v_i,        np.dot(J0_fourier[i, j], v_j))
			H[num_sites + i, num_sites + j] += S   * np.dot(u_i.conj(), np.dot(J_fourier[i, j],  u_j))
			H[num_sites + i, num_sites + i] -= S_j * np.dot(v_i,        np.dot(J0_fourier[i, j], v_j))
			H[            i, num_sites + j] += S   * np.dot(u_i,        np.dot(J_fourier[i, j],  u_j))
			H[num_sites + i,             j] += (S  * np.dot(u_j,        np.dot(J_fourier[j, i],  u_i))).conj()

		if field != None:
			zee = mu_B * g_e * field["mag"]*field["dir"] @ v_i

			H[            i,             i] -= zee
			H[num_sites + i, num_sites + i] -= zee.conjugate()
	print_infos("\nH =\n%s" % H)

	# trafo
	C = la.cholesky(H, upper = True)
	signs = np.diag(np.concatenate((np.repeat(1, num_sites), np.repeat(-1, num_sites))))
	H_trafo = np.dot(C, np.dot(signs, C.transpose().conj()))
	print_infos("\nC =\n%s\n\nH_trafo =\n%s" % (C, H_trafo))

	# the eigenvalues of H give the energies
	#Es = la.eigvalsh(H_trafo)
	Es, states = la.eigh(H_trafo)
	print_infos("\nEs = %s" % Es)

	# sort by the energies in descending order
	Es = np.flip(Es)
	states = np.flip(states, axis = 1)
	return Es, states, H_trafo, C, signs

# get the spin-spin correlation at the momentum transfer Qvec
def get_correlations(Qvec, states, H, C, signs, sites):
	num_sites = len(sites)

	energy_mat = np.dot(states.transpose().conj(), np.dot(H, states))
	#Es = np.diag(energy_mat)
	E_sqrt = np.sqrt(np.dot(signs, energy_mat))
	boson_ops = np.dot(la.inv(C), np.dot(states, E_sqrt))

	S_mats = np.zeros((2*num_sites, 3, 3), dtype = complex)
	for x in range(3):
		for y in range(3):
			M = np.zeros((2*num_sites, 2*num_sites), dtype = complex)
			for i in range(num_sites):
				S_i = sites[i]["S"]
				u_i = g_e * sites[i]["u"]
				for j in range(num_sites):
					S_j = sites[j]["S"]
					u_j = g_e * sites[j]["u"]

					S = np.sqrt(S_i * S_j)
					e = np.exp(2j*np.pi * np.dot(Qvec, (np.array(sites[j]["pos"]) - np.array(sites[i]["pos"]))))

					M[            i,             j] = e * S * u_i[x]        * u_j[y].conj()
					M[            i, num_sites + j] = e * S * u_i[x]        * u_j[y]
					M[num_sites + i,             j] = e * S * u_i[x].conj() * u_j[y].conj()
					M[num_sites + i, num_sites + j] = e * S * u_i[x].conj() * u_j[y]

			M = np.dot(boson_ops.transpose().conj(), np.dot(M, boson_ops))
			for E_idx in range(num_sites * 2):
				S_mats[E_idx, x, y] += M[E_idx, E_idx] / (2. * num_sites)

	proj = np.eye(3) - np.outer(Qvec, Qvec) / la.norm(Qvec)**2.
	weights = []
	for E_idx in range(num_sites * 2):
		S_mats[E_idx, :, :] = np.dot(proj, np.dot(S_mats[E_idx, :, :], proj))
		weights.append(np.abs(S_mats[E_idx, :, :].trace().real))

	return weights
)BLOCK";

	std::string filename = _filename.toStdString();

	std::ofstream ofstr(filename);
	if(!ofstr)
	{
		ShowError(QString("Cannot open file \"%1\" for writing.").arg(filename.c_str()));
		return false;
	}

	ofstr.precision(g_prec);

	// header
	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr << "#\n"
		<< "# Generated by Magpie\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n\n";

	// imports
	ofstr << "import numpy as np\n";
	ofstr << "import numpy.linalg as la\n\n";


	// --------------------------------------------------------------------
	t_real h1 = (t_real)m_Q_start[0]->value();
	t_real k1 = (t_real)m_Q_start[1]->value();
	t_real l1 = (t_real)m_Q_start[2]->value();
	t_real h2 = (t_real)m_Q_end[0]->value();
	t_real k2 = (t_real)m_Q_end[1]->value();
	t_real l2 = (t_real)m_Q_end[2]->value();

	ofstr << "\n# constants and variables\n";
	ofstr << "g_e     = " << tl2::g_e<t_real> << "\n";
	ofstr << "mu_B    = " << t_real(tl2::mu_B<t_real> / tl2::meV<t_real> * tl2::tesla<t_real>) << "\n\n";

	ofstr << "Qstart  = np.array([ " << h1 << ", " << k1 << ", " << l1 << " ])\n";
	ofstr << "Qend    = np.array([ " << h2 << ", " << k2 << ", " << l2 << " ])\n";
	ofstr << "Qpts    = " << m_num_points->value() << "\n\n";

	// user variables
	for(const auto &var : m_dyn.GetVariables())
	{
		ofstr << var.name << " = " << var.value.real();
		if(!tl2::equals_0<t_real>(var.value.imag(), g_eps))
			ofstr << " + var.value.imag()" << "j";
		ofstr << "\n";
	}
	// --------------------------------------------------------------------


	// script functions
	ofstr << scr << "\n";


	// --------------------------------------------------------------------
	const auto& field = m_dyn.GetExternalField();
	bool use_field = !tl2::equals_0<t_real>(field.mag, g_eps);

	if(use_field)
	{
		ofstr << "\n# external field\n";
		ofstr << "field = { ";
		ofstr << "\"dir\" : -np.array([ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ]), "
			<< "\"mag\" : " << field.mag
			<< " }\n";
	}
	else
	{
		ofstr << "\n# no external field\n";
		ofstr << "field = None\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# xtal lattice and magnetic sites\n";

	const auto& xtal = m_dyn.GetCrystalLattice();
	ofstr << "xtal_lattice = [ "
		<< xtal[0] << ", " << xtal[1] << ", " << xtal[2] << " ]\n";
	ofstr << "xtal_angles = [ "
		<< tl2::r2d<t_real>(xtal[0]) << ", "
		<< tl2::r2d<t_real>(xtal[1]) << ", "
		<< tl2::r2d<t_real>(xtal[2]) << " ]\n";

	ofstr << "sites = [\n";
	for(const t_site &site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t{ \"name\" : \"" << site.name << "\", "
			<< "\"S\" : " << get_str_var(site.spin_mag);

		if(field.align_spins)
		{
			ofstr << ", \"Sdir\" : [ "
				<< -field.dir[0] << ", "
				<< -field.dir[1] << ", "
				<< -field.dir[2] << " ]";
		}
		else
		{
			ofstr << ", \"Sdir\" : [ "
				<< get_str_var(site.spin_dir[0]) << ", "
				<< get_str_var(site.spin_dir[1]) << ", "
				<< get_str_var(site.spin_dir[2]) << " ]";
		}
		
		ofstr << ", \"pos\" : [ "
			<< get_str_var(site.pos[0]) << ", "
			<< get_str_var(site.pos[1]) << ", "
			<< get_str_var(site.pos[2]) << " ]"
			<< " },\n";
	}
	ofstr << "]\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";

	ofstr << "couplings = [\n";
	for(const t_term& term : m_dyn.GetExchangeTerms())
	{
		t_size idx1 = m_dyn.GetMagneticSiteIndex(term.site1);
		t_size idx2 = m_dyn.GetMagneticSiteIndex(term.site2);

		ofstr << "\t{ \"name\" : \"" << term.name << "\", "
			<< "\"sites\" : [ " << idx1 << ", " << idx2 << " ], "
			<< "\"dist\" : [ "
			<< get_str_var(term.dist[0]) << ", "
			<< get_str_var(term.dist[1]) << ", "
			<< get_str_var(term.dist[2]) << " ], "
			<< "\"J\" : " << get_str_var(term.J, true) << ", "
			<< "\"DMI\" : [ "
			<< get_str_var(term.dmi[0], true) << ", "
			<< get_str_var(term.dmi[1], true) << ", "
			<< get_str_var(term.dmi[2], true) << " ]"
			<< " },\n";
	}
	ofstr << "]\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	if(m_dyn.IsIncommensurate())
	{
		ofstr << "\n# incommensurability\n";

		const t_vec_real& prop = m_dyn.GetOrderingWavevector();
		const t_vec_real& axis = m_dyn.GetRotationAxis();
		t_vec_real s0 = tl2::cross(prop, axis);
		s0 /= tl2::norm(s0);

		ofstr << "k = [ " << prop[0] << ", " << prop[1] << ", " << prop[2] << " ]\n"
			<< "axis = [ " << axis[0] << ", " << axis[1] << ", " << axis[2] << " ]\n"
			<< "S0 = [ " << s0[0] << ", " << s0[1] << ", " << s0[2] << " ]\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	
	ofstr <<  R"BLOCK(
# create magnetic structure
init(sites, couplings)

# plot the dispersion branch
import matplotlib.pyplot as plt

hs, ks, ls, Es, ws = [], [], [], [], []
for Qidx in range(Qpts):
	try:
		Qvec = Qstart + (Qend - Qstart) * Qidx / Qpts
		allEs, allstates, H, C, signs = get_energies(Qvec, sites, couplings, field)
		allws = get_correlations(Qvec, allstates, H, C, signs, sites)
		for E, state, w in zip(allEs, allstates, allws):
			if only_pos_E and E < 0.:
				continue
			hs.append(Qvec[0])
			ks.append(Qvec[1])
			ls.append(Qvec[2])
			Es.append(E)
			ws.append(w * weight_scale)
	except la.LinAlgError:
		pass

# choose the Q component with the most change
qs = hs
if np.std(ks) > np.std(qs):
	qs = ks
if np.std(ls) > np.std(qs):
	qs = ls

plt.plot()
plt.xlabel("q (rlu)")
plt.ylabel("E (meV)")
plt.scatter(qs, Es, marker = '.', s = ws)
plt.show()
)BLOCK";
	// --------------------------------------------------------------------

	return true;
}
