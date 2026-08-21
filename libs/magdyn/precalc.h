/**
 * magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAGDYN_PRECALC_H__
#define __MAGDYN_PRECALC_H__


#include "tlibs2/libs/units.h"

#include "magdyn.h"


// --------------------------------------------------------------------
// pre-calculation functions
// --------------------------------------------------------------------
/**
 * calculate the rotation matrix for the external field
 */
MAGDYN_TEMPL void MAGDYN_INST::CalcExternalField()
{
	const bool use_field_rot =
		(!tl2::equals_0<t_real>(m_field.mag, m_eps) || m_field.align_spins)
		&& !!m_field.dir;
	if(!use_field_rot)
		return;

	// rotate field to [001] direction
	m_rot_field = tl2::convert<t_mat33>(
		tl2::trans<t_mat33_real>(
			tl2::rotation<t_mat33_real, t_vec3_real>(
				-*m_field.dir, m_zdir, &m_rotaxis, m_eps)));

	// rotate -field to [001] direction
	if(m_field.keep_signs)
	{
		m_rot_negfield = tl2::convert<t_mat33>(
			tl2::trans<t_mat33_real>(
				tl2::rotation<t_mat33_real, t_vec3_real>(
					*m_field.dir, m_zdir, &m_rotaxis, m_eps)));
	}

#ifdef __MAGDYN_DEBUG_OUTPUT__
	std::cout << "Field rotation from:\n";
	tl2::niceprint(std::cout, -(*m_field.dir), 1e-4, 4);
	std::cout << "\nto:\n";
	tl2::niceprint(std::cout, m_zdir, 1e-4, 4);
	std::cout << "\nmatrix:\n";
	tl2::niceprint(std::cout, m_rot_field, 1e-4, 4);
	std::cout << std::endl;
#endif
}



/**
 * calculate the spin rotation trafo for the magnetic sites
 * and parse any given expressions
 */
MAGDYN_TEMPL void MAGDYN_INST::CalcMagneticSite(MagneticSite& site)
{
	try
	{
		tl2::ExprParser<t_cplx> parser = GetExprParser();

		// spin direction and orthogonal plane
		bool has_explicit_trafo = true;

		// defaults
		site.pos_calc = tl2::zero<t_vec3_real>(3);
		site.spin_dir_calc = tl2::zero<t_vec3_real>(3);
		site.trafo_z_calc = tl2::zero<t_vec3>(3);
		site.trafo_plane_calc = tl2::zero<t_vec3>(3);
		site.trafo_plane_conj_calc = tl2::zero<t_vec3>(3);
		if(!site.g_e)
			site.g_e = tl2::g_e<t_real> * tl2::unit<t_mat33>(3);

		// spin magnitude
		if(parser.parse_noexcept(site.spin_mag))
		{
			site.spin_mag_calc = parser.eval_noexcept().real();
		}
		else
		{
			MAGDYN_CERR_OPT << "Magdyn error: Parsing spin magnitude \""
				<< site.spin_mag << "\""
				<< " for site \"" << site.name << "\""
				<< "." << std::endl;
		}

		for(std::uint8_t idx = 0; idx < 3; ++idx)
		{
			// position
			if(site.pos[idx] != "")
			{
				if(parser.parse_noexcept(site.pos[idx]))
				{
					site.pos_calc[idx] = parser.eval_noexcept().real();
				}
				else
				{
					MAGDYN_CERR_OPT << "Magdyn error: Parsing position \""
						<< site.pos[idx] << "\""
						<< " for site \"" << site.name << "\""
						<< " and component " << int(idx)
						<< "." << std::endl;
				}
			}

			// spin direction
			if(site.spin_dir[idx] != "")
			{
				if(parser.parse_noexcept(site.spin_dir[idx]))
				{
					site.spin_dir_calc[idx] = parser.eval_noexcept().real();
				}
				else
				{
					MAGDYN_CERR_OPT << "Magdyn error: Parsing spin direction \""
						<< site.spin_dir[idx] << "\""
						<< " for site \"" << site.name << "\""
						<< " and component " << int(idx)
						<< "." << std::endl;
				}
			}

			// orthogonal spin direction
			if(site.spin_ortho[idx] != "")
			{
				if(parser.parse_noexcept(site.spin_ortho[idx]))
				{
					site.trafo_plane_calc[idx] = parser.eval_noexcept();
					site.trafo_plane_conj_calc[idx] =
						std::conj(site.trafo_plane_calc[idx]);
				}
				else
				{
					has_explicit_trafo = false;

					MAGDYN_CERR_OPT << "Magdyn error: Parsing spin orthogonal plane \""
						<< site.spin_ortho[idx] << "\""
						<< " for site \"" << site.name << "\""
						<< " and component " << int(idx)
						<< "." << std::endl;
				}
			}
			else
			{
				has_explicit_trafo = false;
			}
		}  // idx

		// spin rotation of equation (9) from (Toth 2015)
		if(m_field.align_spins && m_field.dir)
		{
			if(m_field.keep_signs && tl2::inner(site.spin_dir_calc, *m_field.dir) > 0.)
			{
				std::tie(site.trafo_plane_calc, site.trafo_z_calc) =
					rot_to_trafo(m_rot_negfield);
			}
			else
			{
				std::tie(site.trafo_plane_calc, site.trafo_z_calc) =
					rot_to_trafo(m_rot_field);
			}
		}
		else
		{
			if(!has_explicit_trafo)
			{
				// calculate u and v from the spin rotation
				std::tie(site.trafo_plane_calc, site.trafo_z_calc) =
					spin_to_trafo(site.spin_dir_calc);
			}

			// TODO: normalise the v vector as well as the real and imaginary u vectors
			// in case they are explicitly given

#ifdef __MAGDYN_DEBUG_OUTPUT__
			std::cout << "Site " << site.name << ": u = "
				<< site.trafo_plane_calc[0] << " "
				<< site.trafo_plane_calc[1] << " "
				<< site.trafo_plane_calc[2]
				<< std::endl;
			std::cout << "Site " << site.name << ": v = "
				<< site.trafo_z_calc[0] << " "
				<< site.trafo_z_calc[1] << " "
				<< site.trafo_z_calc[2]
				<< std::endl;
#endif
		}

		if(!has_explicit_trafo)
			site.trafo_plane_conj_calc = tl2::conj(site.trafo_plane_calc);

		// multiply g factor
		site.ge_trafo_z_calc = (*site.g_e) * site.trafo_z_calc;
		site.ge_trafo_plane_calc = (*site.g_e) * site.trafo_plane_calc;
		site.ge_trafo_plane_conj_calc = (*site.g_e) * site.trafo_plane_conj_calc;

	}
	catch(const std::exception& ex)
	{
		MAGDYN_CERR_OPT << "Magdyn error: Calculating site \"" << site.name << "\"."
			<< " Reason: " << ex.what()
			<< std::endl;
	}
}



/**
 * calculate the spin rotation trafo for the magnetic sites
 * and parse any given expressions
 */
MAGDYN_TEMPL void MAGDYN_INST::CalcMagneticSites()
{
	for(MagneticSite& site : GetMagneticSites())
		CalcMagneticSite(site);
}



/**
 * parse the exchange term expressions and calculate all properties
 */
MAGDYN_TEMPL void MAGDYN_INST::CalcExchangeTerm(MAGDYN_TYPE::ExchangeTerm& term)
{
	try
	{
		tl2::ExprParser<t_cplx> parser = GetExprParser();

		// defaults
		term.dist_calc = tl2::zero<t_vec3_real>(3);  // distance
		term.dmi_calc  = std::nullopt;
		term.Jgen_calc = std::nullopt;

		// get site indices
		term.site1_calc = GetMagneticSiteIndex(term.site1);
		term.site2_calc = GetMagneticSiteIndex(term.site2);

#ifdef __MAGDYN_DEBUG_OUTPUT__
		std::cout << "Coupling: "
			<< term.name << ": " << term.site1 << " (" << term.site1_calc << ") -> "
			<< term.site2 << " (" << term.site2_calc << ")." << std::endl;
#endif

		if(term.site1_calc >= GetMagneticSitesCount())
		{
			MAGDYN_CERR_OPT << "Magdyn error: Unknown site 1 name \"" << term.site1 << "\"."
				<< " in coupling \"" << term.name << "\"."
				<< std::endl;
			return;
		}
		if(term.site2_calc >= GetMagneticSitesCount())
		{
			MAGDYN_CERR_OPT << "Magdyn error: Unknown site 2 name \"" << term.site2 << "\"."
				<< " in coupling \"" << term.name << "\"."
				<< std::endl;
			return;
		}

		// register the spin magnitude as parser variable
		const MagneticSite& site1 = GetMagneticSite(term.site1_calc);
		const MagneticSite& site2 = GetMagneticSite(term.site2_calc);
		parser.register_var("_S1_mag", site1.spin_mag_calc);
		parser.register_var("_S2_mag", site2.spin_mag_calc);
		if(term.site1_calc == term.site2_calc)
			parser.register_var("_S_mag", site1.spin_mag_calc);

		// symmetric interaction
		if(term.J == "")
		{
			term.J_calc = t_real(0);
		}
		else if(parser.parse_noexcept(term.J))
		{
			term.J_calc = parser.eval_noexcept();
		}
		else
		{
			MAGDYN_CERR_OPT << "Magdyn error: Parsing J term \""
				<< term.J << "\"." << std::endl;
		}

		// check if dmi is enabled
		auto _has_dmi = [&term]() -> bool
		{
			for(std::uint8_t i = 0; i < 3; ++i)
				if(term.dist[i] != "")
					return true;
			return false;
		};

		// check if the general J matrix is enabled
		auto _has_Jgen = [&term]() -> bool
		{
			for(std::uint8_t i = 0; i < 3; ++i)
				for(std::uint8_t j = 0; j < 3; ++j)
					if(term.Jgen[i][j] != "")
						return true;
			return false;
		};

		const bool has_dmi = _has_dmi();
		const bool has_Jgen = _has_Jgen();

		// enable dmi
		if(has_dmi)
			term.dmi_calc  = tl2::zero<t_vec3>(3);
		// enable general J
		if(has_Jgen)
			term.Jgen_calc = tl2::zero<t_mat33>(3, 3);

		for(std::uint8_t i = 0; i < 3; ++i)
		{
			// distance
			if(term.dist[i] != "")
			{
				if(parser.parse_noexcept(term.dist[i]))
				{
					term.dist_calc[i] = parser.eval_noexcept().real();
				}
				else
				{
					MAGDYN_CERR_OPT << "Magdyn error: Parsing distance term \""
						<< term.dist[i]
						<< "\" (index " << int(i) << ")"
						<< "." << std::endl;
				}
			}

			// dmi
			if(has_dmi && term.dmi[i] != "")
			{
				if(parser.parse_noexcept(term.dmi[i]))
				{
					(*term.dmi_calc)[i] = parser.eval_noexcept();
				}
				else
				{
					MAGDYN_CERR_OPT << "Magdyn error: Parsing DMI term \""
						<< term.dmi[i]
						<< "\" (index " << int(i) << ")"
						<< "." << std::endl;
				}
			}

				// general exchange interaction
			if(has_Jgen)
			{
				for(std::uint8_t j = 0; j < 3; ++j)
				{
					if(term.Jgen[i][j] == "")
						continue;

					if(parser.parse_noexcept(term.Jgen[i][j]))
					{
						(*term.Jgen_calc)(i, j) = parser.eval_noexcept();
					}
					else
					{
						MAGDYN_CERR_OPT << "Magdyn error: Parsing general term \""
							<< term.Jgen[i][j]
							<< "\" (indices " << int(i) << ", "
							<< int(j) << ")" << "." << std::endl;
					}
				}  // j iteration
			}  // has_Jgen
		}  // i iteration

		const t_vec3_real& pos1_uc = GetMagneticSite(term.site1_calc).pos_calc;
		const t_vec3_real& pos2_uc = GetMagneticSite(term.site2_calc).pos_calc;
		t_vec3_real pos2_sc = pos2_uc + term.dist_calc;

		// transform to lab units for correct distance
		t_vec3_real pos1_uc_lab = m_xtalA * pos1_uc;
		t_vec3_real pos2_sc_lab = m_xtalA * pos2_sc;

		term.length_calc = tl2::norm<t_vec3_real>(pos2_sc_lab - pos1_uc_lab);
	}
	catch(const std::exception& ex)
	{
		MAGDYN_CERR_OPT << "Magdyn error: Calculating coupling \""
			<< term.name << "\"."
			<< " Reason: " << ex.what() << "."
			<< std::endl;
	}
}



/**
 * parse all exchange term expressions and calculate all properties
 */
MAGDYN_TEMPL void MAGDYN_INST::CalcExchangeTerms()
{
	for(ExchangeTerm& term : GetExchangeTerms())
		CalcExchangeTerm(term);
}



/**
 * converts the rotation matrix rotating the local spins to ferromagnetic
 * [001] directions into the vectors comprised of the matrix columns
 * @see equation (9) and (51) from (Toth 2015)
 */
MAGDYN_TEMPL
std::tuple<typename MAGDYN_INST::t_vec3, typename MAGDYN_INST::t_vec3>
MAGDYN_INST::rot_to_trafo(const MAGDYN_INST::t_mat33& R)
{
	const t_vec3 xy_plane = tl2::col<t_mat33, t_vec3>(R, 0)
		+ s_imag * tl2::col<t_mat33, t_vec3>(R, 1);
	const t_vec3 z = tl2::col<t_mat33, t_vec3>(R, 2);

	return std::make_tuple(xy_plane, z);
}



/**
 * rotate local spin to ferromagnetic [001] direction
 * @see equations (7) and (9) from (Toth 2015)
 */
MAGDYN_TEMPL
std::tuple<typename MAGDYN_INST::t_vec3, typename MAGDYN_INST::t_vec3>
MAGDYN_INST::spin_to_trafo(const MAGDYN_INST::t_vec3_real& spin_dir)
{
	const t_mat33_real _rot = tl2::trans(tl2::rotation<t_mat33_real, t_vec3_real>(
		spin_dir, m_zdir, &m_rotaxis, m_eps));

	const t_mat33 rot = tl2::convert<t_mat33, t_mat33_real>(_rot);
	return rot_to_trafo(rot);
}
// --------------------------------------------------------------------



// --------------------------------------------------------------------
// pre-calculation functions
// --------------------------------------------------------------------
/**
 * rotate a Q vector given in rlu
 */
MAGDYN_TEMPL
typename MAGDYN_INST::t_vec3_real
MAGDYN_INST::RotateQ(const typename MAGDYN_INST::t_vec3_real& Q_rlu,
	const typename MAGDYN_INST::t_vec3_real& axis_rlu, t_real angle) const
{
	t_vec3_real axis = axis_rlu;
	t_vec3_real Q = Q_rlu;

	const t_mat33_real& xtalB = GetCrystalBTrafo();
	auto [xtalB_inv, inv_ok] = tl2::inv(xtalB);

	// get Q and the axis in A^(-1)
	if(inv_ok)
	{
		axis = xtalB * axis;
		Q = xtalB * Q_rlu;
	}

	// rotate Q in the orthonormal A^(-1) system
	t_mat33_real R = tl2::rotation<t_mat33_real, t_vec3_real>(axis, angle, false);
	Q = R * Q;

	// convert Q back to rlu
	if(inv_ok)
		Q = xtalB_inv * Q;

	tl2::set_eps_0(Q, m_eps);
	return Q;
};



/**
 * rotate a Q vector given in rlu
 */
MAGDYN_TEMPL
std::array<t_real, 3>
MAGDYN_INST::RotateQ(t_real h, t_real k, t_real l,
	t_real axis_h, t_real axis_k, t_real axis_l, t_real angle) const
{
	const t_vec3_real Q = tl2::create<t_vec3_real>({ h, k, l });
	const t_vec3_real axis = tl2::create<t_vec3_real>({ axis_h, axis_k, axis_l });

	const t_vec3_real Qrot = RotateQ(Q, axis, angle);

	return std::array<t_real, 3>({ Qrot[0], Qrot[1], Qrot[2] });
}
// --------------------------------------------------------------------

#endif
