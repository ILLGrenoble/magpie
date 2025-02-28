/**
 * tlibs2
 * physical units library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015-2021
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __TLIBS2_UNITS__
#define __TLIBS2_UNITS__

#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/io.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#include <boost/math/constants/constants.hpp>

namespace units = boost::units;
namespace co = boost::units::si::constants::codata;


namespace tl2 {


template<typename T = double> constexpr T __pi = boost::math::constants::pi<T>();


// general quantities
template<class Sys, class T = double> using t_length = units::quantity<units::unit<units::length_dimension, Sys>, T>;
template<class Sys, class T = double> using t_momentum = units::quantity<units::unit<units::momentum_dimension, Sys>, T>;
template<class Sys, class T = double> using t_wavenumber = units::quantity<units::unit<units::wavenumber_dimension, Sys>, T>;
template<class Sys, class T = double> using t_velocity = units::quantity<units::unit<units::velocity_dimension, Sys>, T>;
template<class Sys, class T = double> using t_frequency = units::quantity<units::unit<units::frequency_dimension, Sys>, T>;
template<class Sys, class T = double> using t_energy = units::quantity<units::unit<units::energy_dimension, Sys>, T>;
template<class Sys, class T = double> using t_angle = units::quantity<units::unit<units::plane_angle_dimension, Sys>, T>;
template<class Sys, class T = double> using t_temperature = units::quantity<units::unit<units::temperature_dimension, Sys>, T>;
template<class Sys, class T = double> using t_mass = units::quantity<units::unit<units::mass_dimension, Sys>, T>;
template<class Sys, class T = double> using t_time = units::quantity<units::unit<units::time_dimension, Sys>, T>;
template<class Sys, class T = double> using t_flux = units::quantity<units::unit<units::magnetic_flux_density_dimension, Sys>, T>;
template<class Sys, class T = double> using t_area = units::quantity<units::unit<units::area_dimension, Sys>, T>;
template<class Sys, class T = double> using t_volume = units::quantity<units::unit<units::volume_dimension, Sys>, T>;

template<class Sys, class T = double> using t_length_inverse =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, T>;
template<class Sys, class T = double> using t_length_square =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, 2>::type, Sys>, T>;
template<class Sys, class T = double> using t_momentum_square =
	units::quantity<units::unit<units::derived_dimension<units::momentum_dimension, 2>::type, Sys>, T>;
template<class Sys, class T = double> using t_action =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2, units::time_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_energy_per_temperature =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2,
	units::time_base_dimension,-2, units::temperature_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_energy_per_field =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,1, units::length_base_dimension,2>::type, Sys>, T>;
template<class Sys, class T = double> using t_inv_flux_time =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,1, units::time_base_dimension,1, units::mass_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_dimensionless =
	units::quantity<units::unit<units::dimensionless_type, Sys>, T>;


// synonyms
template<class Sys, class T = double> using t_freq = t_frequency<Sys, T>;
template<class Sys, class T = double> using t_temp = t_temperature<Sys, T>;


// si quantities -- partial specialisations
template<class T = double> using t_length_si = t_length<units::si::system, T>;
template<class T = double> using t_length_inverse_si = t_length_inverse<units::si::system, T>;
template<class T = double> using t_momentum_si = t_momentum<units::si::system, T>;
template<class T = double> using t_wavenumber_si = t_wavenumber<units::si::system, T>;
template<class T = double> using t_velocity_si = t_velocity<units::si::system, T>;
template<class T = double> using t_frequency_si = t_frequency<units::si::system, T>;
template<class T = double> using t_energy_si = t_energy<units::si::system, T>;
template<class T = double> using t_angle_si = t_angle<units::si::system, T>;
template<class T = double> using t_temperature_si = t_temperature<units::si::system, T>;
template<class T = double> using t_mass_si = t_mass<units::si::system, T>;
template<class T = double> using t_time_si = t_time<units::si::system, T>;
template<class T = double> using t_flux_si = t_flux<units::si::system, T>;
template<class T = double> using t_area_si = t_area<units::si::system, T>;
template<class T = double> using t_action_si = t_action<units::si::system, T>;
template<class T = double> using t_energy_per_temperature_si = t_energy_per_temperature<units::si::system, T>;
template<class T = double> using t_inv_flux_time_si = t_inv_flux_time<units::si::system, T>;


// constants
template<class T = double> constexpr t_mass_si<T> kg = T(1) * units::si::kilograms;
template<class T = double> constexpr t_temperature_si<T> kelvin = T(1) * units::si::kelvins;
template<class T = double> constexpr t_temperature_si<T> kelvins = kelvin<T>;
template<class T = double> constexpr t_angle_si<T> deg = T(__pi<T>/T(180.)) * units::si::radians;
template<class T = double> constexpr t_angle_si<T> radian = T(1) * units::si::radians;
template<class T = double> constexpr t_angle_si<T> radians = radian<T>;
template<class T = double> constexpr t_flux_si<T> tesla = T(1) * units::si::teslas;
template<class T = double> constexpr t_flux_si<T> teslas = T(1)*units::si::teslas;
template<class T = double> constexpr t_flux_si<T> kilogauss = T(0.1) * units::si::teslas;
template<class T = double> constexpr t_length_si<T> meter = T(1) * units::si::meters;
template<class T = double> constexpr t_length_si<T> meters = meter<T>;
template<class T = double> constexpr t_length_si<T> cm = meters<T>/T(100);
template<class T = double> constexpr t_time_si<T> seconds = T(1) * units::si::seconds;
template<class T = double> constexpr t_time_si<T> sec = seconds<T>;
template<class T = double> constexpr t_time_si<T> ps = T(1e-12)*seconds<T>;
template<class T = double> constexpr t_area_si<T> barns = T(1e-28) * units::si::meters*units::si::meters;;
template<class T = double> constexpr t_area_si<T> barn = barns<T>;
template<class T = double> constexpr t_energy_si<T> meV = T(1e-3) * T(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts;;
template<class T = double> constexpr t_energy_si<T> eV = T(1) * T(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts;
template<class T = double> constexpr t_length_si<T> angstrom = T(1e-10) * units::si::meters;
template<class T = double> constexpr t_length_si<T> angstroms = angstrom<T>;
template<class T = double> constexpr t_mass_si<T> m_n = T(co::m_n/units::si::kilograms)*units::si::kilograms;
template<class T = double> constexpr t_mass_si<T> m_e = T(co::m_e/units::si::kilograms)*units::si::kilograms;
template<class T = double> constexpr t_mass_si<T> amu = T(co::m_u/units::si::kilograms)*units::si::kilograms;
template<class T = double> constexpr t_action_si<T> hbar = T(co::hbar/units::si::joules/units::si::seconds)*units::si::joules*units::si::seconds;
template<class T = double> constexpr t_action_si<T> h = hbar<T> * T(2)*__pi<T>;
template<class T = double> constexpr t_velocity_si<T> c = T(co::c/units::si::meters*units::si::seconds)*units::si::meters/units::si::seconds;

template<class T = double> constexpr t_length_si<T> r_e = T(co::r_e/units::si::meters)*units::si::meters;
template<class T = double> constexpr T g_n = T(co::g_n.value());
template<class T = double> constexpr T g_e = T(co::g_e.value());
template<class T = double> constexpr t_inv_flux_time<units::si::system, T> gamma_n = T(co::gamma_n*units::si::tesla*units::si::seconds)/units::si::tesla/units::si::seconds;
template<class T = double> constexpr t_inv_flux_time<units::si::system, T> gamma_e = T(co::gamma_e*units::si::tesla*units::si::seconds)/units::si::tesla/units::si::seconds;
template<class T = double> constexpr t_energy_per_temperature<units::si::system, T> kB = T(co::k_B*units::si::kelvin/units::si::joules)/units::si::kelvin*units::si::joules;
template<class T = double> constexpr t_energy_per_field<units::si::system, T> mu_B = T(co::mu_B/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla;
template<class T = double> constexpr t_energy_per_field<units::si::system, T> mu_n = T(co::mu_n/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla;
template<class T = double> constexpr t_energy_per_field<units::si::system, T> mu_N = T(co::mu_N/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla;
template<class T = double> constexpr t_energy_per_field<units::si::system, T> mu_e = T(co::mu_e/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla;


// helper functions
template<class t_quant>
t_quant my_units_sqrt(const decltype(t_quant() * t_quant())& val)
{
	using t_quant_sq = decltype(t_quant() * t_quant());
	using T = typename t_quant::value_type;

	t_quant one_quant = t_quant::from_value(T(1));
	t_quant_sq one_quant_sq = t_quant_sq::from_value(T(1));

	T valsq = T(val / one_quant_sq);
	return std::sqrt(valsq) * one_quant;
}


template<class t_quant>
decltype(t_quant()*t_quant()) my_units_pow2(const t_quant& val)
{
	return val*val;
}


template<class t_elem, template<class...> class t_vec>
t_elem my_units_norm2(const t_vec<t_elem>& vec)
{
	using t_elem_sq = decltype(t_elem()*t_elem());
	t_elem_sq tRet = t_elem_sq();

	for(std::size_t i=0; i<vec.size(); ++i)
		tRet += vec[i]*vec[i];

	return my_units_sqrt<t_elem>(tRet);
}

}
#endif
