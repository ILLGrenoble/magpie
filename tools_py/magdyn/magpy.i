/**
 * swig interface for magdyn
 * @author Tobias Weber <tweber@ill.fr>
 * @date 12-oct-2023
 * @license see 'LICENSE' file
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

%module magpy
%{
	#ifndef MAGDYN_COMPAT
		#define MAGDYN_COMPAT
	#endif
	#include "../../libs/magdyn.h"
%}

#define MAGDYN_COMPAT

%include "std_vector.i"
%include "std_list.i"
%include "std_array.i"
%include "std_string.i"
%include "std_complex.i"
%include "std_shared_ptr.i"


//%template(CplxD) std::complex<double>;
%template(VecD) std::vector<double>;
%template(VecCplxD) std::vector<std::complex<double>>;

%template(LstStr) std::list<std::string>;

%template(ArrD3) std::array<double, 3>;
%template(ArrD4) std::array<double, 4>;
%template(ArrCplx3) std::array<std::complex<double>, 3>;
%template(ArrCplx4) std::array<std::complex<double>, 4>;
%template(ArrStr3) std::array<std::string, 3>;
%template(ArrStr33) std::array<std::array<std::string, 3>, 3>;


// ----------------------------------------------------------------------------
// matrix and vector containers
// ----------------------------------------------------------------------------
//%include "maths.h"

// TODO: add these as soon as swig supports concepts

//%template(VectorD) tl2::vec<double>;
//%template(VectorCplx) tl2::vec<std::complex<double>>;
//%template(MatrixD) tl2::mat<double>;
//%template(MatrixCplx) tl2::mat<std::complex<double>>;

//%template(Vector3D) tl2::vec_static<double, 3>;
//%template(Vector3Cplx) tl2::vec_static<std::complex<double>, 3>;
//%template(Matrix33D) tl2::mat_static<double, 3, 3>;
//%template(Matrix33Cplx) tl2::mat_static<std::complex<double>, 3, 3>;
//%template(Vector4D) tl2::vec_static<double, 4>;
//%template(Vector4Cplx) tl2::vec_static<std::complex<double>, 4>;
//%template(Matrix44D) tl2::mat_static<double, 4, 4>;
//%template(Matrix44Cplx) tl2::mat_static<std::complex<double>, 4, 4>;
// ----------------------------------------------------------------------------


%include "../../libs/magdyn.h"
%include "../../libs/magdyn/magdyn.h"
%include "../../libs/magdyn/structs.h"
%include "../../libs/magdyn/helpers.h"
%include "../../libs/magdyn/checks.h"
%include "../../libs/magdyn/getters.h"
%include "../../libs/magdyn/generators.h"
%include "../../libs/magdyn/file.h"
%include "../../libs/magdyn/configuration.h"
%include "../../libs/magdyn/groundstate.h"
%include "../../libs/magdyn/precalc.h"
%include "../../libs/magdyn/hamilton.h"
%include "../../libs/magdyn/correlation.h"
%include "../../libs/magdyn/polarisation.h"
%include "../../libs/magdyn/dispersion.h"
%include "../../libs/magdyn/powder.h"
%include "../../libs/magdyn/topology.h"
%include "../../libs/magdyn/diff.h"


// ----------------------------------------------------------------------------
// input- and output structs (and vectors of them)
// ----------------------------------------------------------------------------
%template(MagneticSite) magdyn::t_MagneticSite<
	tl2::mat_static<std::complex<double>, 3, 3>,
	tl2::vec_static<std::complex<double>, 3>,
	tl2::vec_static<double, 3>,
	std::size_t, double>;
%template(VecMagneticSite) std::vector<
	magdyn::t_MagneticSite<
		tl2::mat_static<std::complex<double>, 3, 3>,
		tl2::vec_static<std::complex<double>, 3>,
		tl2::vec_static<double, 3>,
		std::size_t, double>
	>;
%template(VecMagneticSitePtr) std::vector<
	const magdyn::t_MagneticSite<
		tl2::mat_static<std::complex<double>, 3, 3>,
		tl2::vec_static<std::complex<double>, 3>,
		tl2::vec_static<double, 3>,
		std::size_t, double>*
	>;

%template(ExchangeTerm) magdyn::t_ExchangeTerm<
	tl2::mat_static<std::complex<double>, 3, 3>,
	tl2::vec_static<std::complex<double>, 3>,
	tl2::vec_static<double, 3>,
	std::size_t, std::complex<double>, double>;
%template(VecExchangeTerm) std::vector<
	magdyn::t_ExchangeTerm<
		tl2::mat_static<std::complex<double>, 3, 3>,
		tl2::vec_static<std::complex<double>, 3>,
		tl2::vec_static<double, 3>,
		std::size_t, std::complex<double>, double>
	>;

%template(ExternalField) magdyn::t_ExternalField<
	tl2::vec_static<double, 3>, double>;

%template(EnergyAndWeight) magdyn::t_EnergyAndWeight<
	tl2::mat_static<std::complex<double>, 3, 3>,
	tl2::vec<std::complex<double>>,
	double, std::size_t, std::complex<double>>;
%template(VecEnergyAndWeight) std::vector<
	magdyn::t_EnergyAndWeight<
		tl2::mat_static<std::complex<double>, 3, 3>,
		tl2::vec<std::complex<double>>,
		double, std::size_t, std::complex<double>>
	>;

%template(SofQE) magdyn::t_SofQE<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	double, std::size_t, std::complex<double>,
	tl2::vec_static<double, 3>,
	tl2::mat_static<std::complex<double>, 3, 3>>;
%template(VecSofQE) std::vector<
	magdyn::t_SofQE<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		double, std::size_t, std::complex<double>,
		tl2::vec_static<double, 3>,
		tl2::mat_static<std::complex<double>, 3, 3>>
	>;

%template(Variable) magdyn::t_Variable<
	std::complex<double>>;
%template(VecVariable) std::vector<
	magdyn::t_Variable<
		std::complex<double>>
	>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
/**
 * main magdyn class
 */
%template(MagDyn) magdyn::MagDyn<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	std::complex<double>, double, std::size_t>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// additional functions
// ----------------------------------------------------------------------------
%inline
%{
	#include "../../libs/sym/loadcif.h"
	#include "../../libs/vers.h"


	/**
	 * types
	 */
	using t_str = std::string;
	using t_real = double;
	using t_cplx = std::complex<double>;
	using t_mat = tl2::mat<std::complex<double>>;
	using t_vec = tl2::vec<std::complex<double>>;
	//using t_mat_real = tl2::mat<double>;
	//using t_vec_real = tl2::vec<double>;

	using t_MagDyn = magdyn::MagDyn<t_mat, t_vec, t_cplx, t_real, std::size_t>;

	using t_vec3_real = typename t_MagDyn::t_vec3_real;
	using t_mat44_real = typename t_MagDyn::t_mat44_real;
	using t_SofQE = typename t_MagDyn::SofQE;


	/**
	 * magpie version number
	 */
	t_str get_version()
	{
		return MAGPIE_VER;
	}


	/**
	 * adds a variable
	 */
	void add_variable(t_MagDyn& magdyn,
		const t_str& name,
		const t_cplx& val)
	{
		typename t_MagDyn::Variable var{};

		var.name = name;
		var.value = val;

		// set existing or add new variable
		magdyn.SetVariable(std::move(var));
	}


	/**
	 * sets an external magnetic field
	 */
	void set_field(t_MagDyn& magdyn,
		t_real Bx, t_real By, t_real Bz, t_real Bmag = 1.,
		bool align_spins = true)
	{
		typename t_MagDyn::ExternalField field{};

		field.dir = tl2::create<t_vec3_real>({ Bx, By, Bz });
		field.mag = Bmag;
		field.align_spins = align_spins;

		magdyn.SetExternalField(field);
		magdyn.CalcExternalField();
	}


	/**
	 * sets the temperature
	 */
	void set_temperature(t_MagDyn& magdyn, t_real T)
	{
		magdyn.SetTemperature(T);
	}


	/**
	 * sets up an ordering wave vector
	 */
	void set_ordering(t_MagDyn& magdyn, t_real x, t_real y, t_real z)
	{
		t_vec3_real vec = tl2::create<t_vec3_real>({ x, y, z });
		magdyn.SetOrderingWavevector(vec);
	}


	/**
	 * sets up the rotation axis for the ordering wave vector
	 */
	void set_rotation(t_MagDyn& magdyn, t_real x, t_real y, t_real z)
	{
		t_vec3_real vec = tl2::create<t_vec3_real>({ x, y, z });
		magdyn.SetRotationAxis(vec);
	}


	/**
	 * adds a magnetic site
	 * (using strings with expressions)
	 */
	void add_site(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& x, const t_str& y, const t_str& z,
		const t_str& sx = "0", const t_str& sy = "0", const t_str& sz = "1",
		const t_str& S = "1")
	{
		typename t_MagDyn::MagneticSite site{};

		site.name = name;

		site.pos[0] = x;
		site.pos[1] = y;
		site.pos[2] = z;

		site.spin_dir[0] = sx;
		site.spin_dir[1] = sy;
		site.spin_dir[2] = sz;
		site.spin_mag = S;

		magdyn.CalcMagneticSite(site);
		magdyn.AddMagneticSite(std::move(site));
	}


	/**
	 * adds a magnetic site
	 * (using scalar types)
	 * TODO: remove this function and use std::variant<std::string, t_real> as soon as swig supports this
	 */
	void add_site(t_MagDyn& magdyn,
		const t_str& name,
		t_real x, t_real y, t_real z,
		t_real sx = 0., t_real sy = 0., t_real sz = 1.,
		t_real S = 1.)
	{
		typename t_MagDyn::MagneticSite site{};

		site.name = name;

		site.pos[0] = tl2::var_to_str(x);
		site.pos[1] = tl2::var_to_str(y);
		site.pos[2] = tl2::var_to_str(z);

		site.spin_dir[0] = tl2::var_to_str(sx);
		site.spin_dir[1] = tl2::var_to_str(sy);
		site.spin_dir[2] = tl2::var_to_str(sz);
		site.spin_mag = tl2::var_to_str(S);

		magdyn.CalcMagneticSite(site);
		magdyn.AddMagneticSite(std::move(site));
	}


	/**
	 * sets a site's spin
	 */
	void set_spin(t_MagDyn& magdyn,
		const t_str& name,
		t_real sx = 0., t_real sy = 0., t_real sz = 1.)
	{
		std::size_t idx = magdyn.GetMagneticSiteIndex(name);
		if(idx >= magdyn.GetMagneticSitesCount())
			return;
		typename t_MagDyn::MagneticSite& site = magdyn.GetMagneticSite(idx);

		site.spin_dir[0] = tl2::var_to_str(sx);
		site.spin_dir[1] = tl2::var_to_str(sy);
		site.spin_dir[2] = tl2::var_to_str(sz);
		//site.spin_mag = tl2::var_to_str(S);

		magdyn.CalcMagneticSite(site);
	}


	/**
	 * adds a coupling term between two magnetic sites
	 * (using strings with expressions)
	 */
	void add_coupling(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& site1, const t_str& site2,
		const t_str& dx, const t_str& dy, const t_str& dz,
		const t_str& J,
		const t_str& dmix = "", const t_str& dmiy = "", const t_str& dmiz = "",
		const t_str& Jxx = "", const t_str& Jxy = "", const t_str& Jxz = "",
		const t_str& Jyx = "", const t_str& Jyy = "", const t_str& Jyz = "",
		const t_str& Jzx = "", const t_str& Jzy = "", const t_str& Jzz = "")
	{
		typename t_MagDyn::ExchangeTerm coupling{};

		coupling.name = name;

		coupling.site1 = site1;
		coupling.site2 = site2;

		coupling.dist[0] = dx;
		coupling.dist[1] = dy;
		coupling.dist[2] = dz;

		coupling.J = J;

		coupling.dmi[0] = dmix;
		coupling.dmi[1] = dmiy;
		coupling.dmi[2] = dmiz;

		coupling.Jgen[0][0] = Jxx;
		coupling.Jgen[0][1] = Jxy;
		coupling.Jgen[0][2] = Jxz;
		coupling.Jgen[1][0] = Jyx;
		coupling.Jgen[1][1] = Jyy;
		coupling.Jgen[1][2] = Jyz;
		coupling.Jgen[2][0] = Jzx;
		coupling.Jgen[2][1] = Jzy;
		coupling.Jgen[2][2] = Jzz;

		magdyn.CalcExchangeTerm(coupling);
		magdyn.AddExchangeTerm(std::move(coupling));
	}


	/**
	 * adds a coupling term between two magnetic sites
	 * (using scalar types)
	 * TODO: remove this function and use std::variant<std::string, t_real> as soon as swig supports this
	 */
	void add_coupling(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& site1, const t_str& site2,
		t_real dx, t_real dy, t_real dz,
		t_real J,
		t_real dmix = 0., t_real dmiy = 0., t_real dmiz = 0.,
		t_real Jxx = 0., t_real Jxy = 0., t_real Jxz = 0.,
		t_real Jyx = 0., t_real Jyy = 0., t_real Jyz = 0.,
		t_real Jzx = 0., t_real Jzy = 0., t_real Jzz = 0.)
	{
		typename t_MagDyn::ExchangeTerm coupling{};

		coupling.name = name;

		coupling.site1 = site1;
		coupling.site2 = site2;

		coupling.dist[0] = tl2::var_to_str(dx);
		coupling.dist[1] = tl2::var_to_str(dy);
		coupling.dist[2] = tl2::var_to_str(dz);

		coupling.J = tl2::var_to_str(J);

		coupling.dmi[0] = tl2::var_to_str(dmix);
		coupling.dmi[1] = tl2::var_to_str(dmiy);
		coupling.dmi[2] = tl2::var_to_str(dmiz);

		coupling.Jgen[0][0] = tl2::var_to_str(Jxx);
		coupling.Jgen[0][1] = tl2::var_to_str(Jxy);
		coupling.Jgen[0][2] = tl2::var_to_str(Jxz);
		coupling.Jgen[1][0] = tl2::var_to_str(Jyx);
		coupling.Jgen[1][1] = tl2::var_to_str(Jyy);
		coupling.Jgen[1][2] = tl2::var_to_str(Jyz);
		coupling.Jgen[2][0] = tl2::var_to_str(Jzx);
		coupling.Jgen[2][1] = tl2::var_to_str(Jzy);
		coupling.Jgen[2][2] = tl2::var_to_str(Jzz);

		magdyn.CalcExchangeTerm(coupling);
		magdyn.AddExchangeTerm(std::move(coupling));
	}


	/**
	 * calculates sites and couplings
	 */
	void calc(t_MagDyn& magdyn)
	{
		magdyn.CalcExternalField();
		magdyn.CalcMagneticSites();
		magdyn.CalcExchangeTerms();
	}


	/**
	 * get a list of all space groups
	 */
	std::list<std::string> get_spacegroups()
	{
		std::list<std::string> sg_names;

		for(const auto& sg : sym::get_sgs<t_mat44_real>(true, false))
			sg_names.push_back(std::get<1>(sg));

		return sg_names;
	}


	/**
	 * create symmetry equivalent site positions using a space group
	 */
	 void symmetrise_sites(t_MagDyn& magdyn, const std::string& spacegroup)
	 {
		std::vector<t_mat44_real> ops = sym::get_sg_ops<t_mat44_real>(spacegroup);
		if(ops.size() == 0)
		{
			std::cerr << "Error: Space group \"" << spacegroup << "\" was not found. "
				<< "Please use get_spacegroups() to get a list."
				<< std::endl;
			return;
		}

		magdyn.SymmetriseMagneticSites(ops);
	 }


	 /**
	  * output the list of magnetic sites
	  */
	 void print_sites(t_MagDyn& magdyn)
	 {
		int w = magdyn.GetPrecision()*1.5;
		std::cout.precision(magdyn.GetPrecision());

		std::cout
			<< std::setw(24) << std::left << "site" << " "
			<< std::setw(w) << std::left << "x" << " "
			<< std::setw(w) << std::left << "y" << " "
			<< std::setw(w) << std::left << "z" << " "
			<< std::setw(w) << std::left << "sx" << " "
			<< std::setw(w) << std::left << "sy" << " "
			<< std::setw(w) << std::left << "sz" << " "
			<< std::setw(w) << std::left << "|S|"
			<< std::endl;

		for(const auto& site : magdyn.GetMagneticSites())
		{
			std::cout
				<< std::setw(24) << std::left << site.name << " "
				<< std::setw(w) << std::left << site.pos_calc[0] << " "
				<< std::setw(w) << std::left << site.pos_calc[1] << " "
				<< std::setw(w) << std::left << site.pos_calc[2] << " "
				<< std::setw(w) << std::left << site.spin_dir_calc[0] << " "
				<< std::setw(w) << std::left << site.spin_dir_calc[1] << " "
				<< std::setw(w) << std::left << site.spin_dir_calc[2] << " "
				<< std::setw(w) << std::left << site.spin_mag_calc
				<< std::endl;
		}
	}


	/**
	 * create symmetry equivalent exchange couplings using a space group
	 */
	 void symmetrise_couplings(t_MagDyn& magdyn, const std::string& spacegroup)
	 {
		std::vector<t_mat44_real> ops = sym::get_sg_ops<t_mat44_real>(spacegroup);
		if(ops.size() == 0)
		{
			std::cerr << "Error: Space group \"" << spacegroup << "\" was not found. "
				<< "Please use get_spacegroups() to get a list."
				<< std::endl;
			return;
		}

		magdyn.SymmetriseExchangeTerms(ops);
	 }


 	/**
	 * create couplings up to a certain distance
	 */
	 void generate_couplings(t_MagDyn& magdyn,
		t_real dist_max, std::size_t sc_max, std::size_t couplings_max,
		const std::string& spacegroup)
	 {
		std::vector<t_mat44_real> ops = sym::get_sg_ops<t_mat44_real>(spacegroup);
		if(ops.size() == 0)
		{
			std::cerr << "Error: Space group \"" << spacegroup << "\" was not found. "
				<< "Please use get_spacegroups() to get a list."
				<< std::endl;
			return;
		}

		magdyn.GeneratePossibleExchangeTerms(dist_max, sc_max, couplings_max);
		magdyn.CalcSymmetryIndices(ops);
	 }


	 /**
	  * output the list of magnetic couplings
	  */
	 void print_couplings(t_MagDyn& magdyn)
	 {
		int w = magdyn.GetPrecision()*1.5;
		std::cout.precision(magdyn.GetPrecision());

		std::cout
			<< std::setw(24) << std::left << "coupling" << " "
			<< std::setw(w) << std::left << "site1" << " "
			<< std::setw(w) << std::left << "site2" << " "
			<< std::setw(w) << std::left << "dx" << " "
			<< std::setw(w) << std::left << "dy" << " "
			<< std::setw(w) << std::left << "dz" << " "
			<< std::setw(w) << std::left << "length" << " "
			<< std::setw(w) << std::left << "J" << " "
			<< std::endl;

		for(const auto& coupling : magdyn.GetExchangeTerms())
		{
			std::cout
				<< std::setw(24) << std::left << coupling.name << " "
				<< std::setw(w) << std::left << coupling.site1_calc << " "
				<< std::setw(w) << std::left << coupling.site2_calc << " "
				<< std::setw(w) << std::left << coupling.dist_calc[0] << " "
				<< std::setw(w) << std::left << coupling.dist_calc[1] << " "
				<< std::setw(w) << std::left << coupling.dist_calc[2] << " "
				<< std::setw(w) << std::left << coupling.length_calc << " "
				<< std::setw(w) << std::left << coupling.J_calc << " "
				<< std::endl;
		}
	}


 	/**
	 * assign the couplings correponsing to a given symmetry index
	 */
	 void assign_couplings_by_symmetry(t_MagDyn& magdyn, std::size_t sym_idx, const std::string& J)
	 {
		magdyn.AssignCouplingsBySymmetryIndex(sym_idx, &J);
	 }


	 /**
	  * helper function to access vector components which are not (yet) seen by swig
	  * TODO: remove this once swig understands concepts
	  */
	 t_real get_h(const t_SofQE& S)
	 {
		return S.Q_rlu[0];
	 }


	 /**
	  * helper function to access vector components which are not (yet) seen by swig
	  * TODO: remove this once swig understands concepts
	  */
	 t_real get_k(const t_SofQE& S)
	 {
		return S.Q_rlu[1];
	 }


	 /**
	  * helper function to access vector components which are not (yet) seen by swig
	  * TODO: remove this once swig understands concepts
	  */
	 t_real get_l(const t_SofQE& S)
	 {
		return S.Q_rlu[2];
	 }
%}
// ----------------------------------------------------------------------------
