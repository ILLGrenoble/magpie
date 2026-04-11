/**
 * system benchmark
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-april-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "tlibs2/libs/magdyn.h"
#include "tlibs2/libs/algos.h"
#include "defs.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;



/**
 * runs an individual benchmark for the given number of Q points and threads
 */
static bool benchmark(t_magdyn& magdyn,
  t_real h_start, t_real k_start, t_real l_start,
  t_real h_end, t_real k_end, t_real l_end,
	t_size num_Qs, t_size num_threads)
{
	// ignore results
	std::function<void(const typename t_magdyn::SofQE*)> results_fkt =
		[](const typename t_magdyn::SofQE*) -> void {};

	//magdyn.CalcEnergies(h, k, l, false);
	magdyn.CalcDispersion(h_start, k_start, l_start,
		h_end, k_end, l_end, num_Qs, num_threads,
		true, nullptr, &results_fkt);

	return true;
}



/**
 * runs all benchmarks
 */
bool benchmark(const std::string& model_file,
	t_size max_threads = 4, t_size max_Qs = 65536)
{
	if(max_threads == 0)
		max_threads = 4;
	if(max_Qs == 0)
		max_Qs = 65536;

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// magnon calculator
	t_magdyn magdyn{};

	// load model from input file
	if(!magdyn.Load(model_file))
		return false;

	// settings
	magdyn.SetSilent(true);
	magdyn.SetPerformChecks(false);

	stopwatch.stop();
	const t_real startup_time = stopwatch.GetDur();


	// print header
	std::cout << "# magpie version " << MAGPIE_VER << std::endl;
	std::cout << "# setup time (ms): " << startup_time * 1000. << std::endl;


	// get the configuration options
	pt::ptree root_node;
	std::ifstream ifstr{model_file};
	pt::read_xml(ifstr, root_node);
	const auto &magdyn_node = root_node.get_child("magdyn");

	// calculate the dispersion along the given Q_i and Q_f positions
	t_real h_start = magdyn_node.get<t_real>("config.h_start", 0.);
	t_real k_start = magdyn_node.get<t_real>("config.k_start", 0.);
	t_real l_start = magdyn_node.get<t_real>("config.l_start", 0.);
	t_real h_end = magdyn_node.get<t_real>("config.h_end", 1.);
	t_real k_end = magdyn_node.get<t_real>("config.k_end", 0.);
	t_real l_end = magdyn_node.get<t_real>("config.l_end", 0.);


	// calculate benchmarks for different Q and thread counts
	int field_w = g_prec * 2.5;
	for(t_size num_threads = 1; num_threads <= max_threads; num_threads *= 2)
	{
		std::cout << "\n# benchmark for " << num_threads << " thread(s)" << std::endl;
		std::cout
			<< std::left << std::setw(field_w) << "# Qs" << " "
			<< std::left << std::setw(field_w) << "calc t (ms)"
			<< std::left << std::setw(field_w) << "full t (ms)"
			<< std::endl;

		for(t_size num_Qs = 1; num_Qs <= max_Qs; num_Qs *= 2)
		{
			stopwatch.start();

			if(!benchmark(magdyn,
				h_start, k_start, l_start,
				h_end, k_end, l_end,
				num_Qs, num_threads))
				return false;

			stopwatch.stop();

			std::cout
				<< std::left << std::setw(field_w) << num_Qs << " "
				<< std::left << std::setw(field_w) << stopwatch.GetDur() * 1000.
				<< std::left << std::setw(field_w) << (stopwatch.GetDur() + startup_time) * 1000.
				<< std::endl;
		}
	}

	return true;
}
