/**
 * magnetic dynamics -- entry point
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef DONT_USE_QT
	#include "magdyn.h"
	#include "tlibs2/libs/qt/helper.h"
#endif

#include "tlibs2/libs/algos.h"
#include "defs.h"

#include <iostream>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace args = boost::program_options;
namespace pt = boost::property_tree;


#ifdef DONT_USE_QT
static int gui_main(int argc, char** argv, const std::string& model_file,
	const t_vec_real& Qi, const t_vec_real& Qf, t_size num_Q_pts)
{
	std::cerr << "Error: The GUI is not available in this version." << std::endl;
	return -1;
}
#else
extern int gui_main(int argc, char** argv, const std::string& model_file,
	const t_vec_real& Qi, const t_vec_real& Qf, t_size num_Q_pts);
#endif



/**
 * starts the cli program
 */
static int cli_main(const std::string& model_file, const std::string& results_file,
	const t_vec_real& Qi, const t_vec_real& Qf, t_size num_Q_pts)
{
	using namespace tl2_ops;

	if(model_file == "")
	{
		std::cerr << "Error: No magnetic model given." << std::endl;
		return -1;
	}


	// load model from input file
	t_magdyn magdyn;
	if(magdyn.Load(model_file))
	{
		std::cout << "Loaded magnetic model from file \"" <<
			model_file << "\"." << std::endl;
	}
	else
	{
		std::cerr << "Error: Failed loading magnetic model \"" <<
			model_file << "\"." << std::endl;
		return -1;
	}


	// print some infos about the model
	std::cout << "Model infos:" << std::endl;

	std::cout << "\t" << magdyn.GetMagneticSitesCount() << " magnetic site(s) and "
		<< magdyn.GetExchangeTermsCount() << " coupling(s) defined." << std::endl;

	if(magdyn.IsIncommensurate())
	{
		std::cout << "\tSystem is incommensurate with ordering vector: "
			<< magdyn.GetOrderingWavevector() << "." << std::endl;
	}
	else
	{
		std::cout << "\tSystem is commensurate." << std::endl;
	}

	t_real T = magdyn.GetTemperature();
	if(T < 0.)
		std::cout << "\tTemperature disabled." << std::endl;
	else
		std::cout << "\tTemperature: " << T << "." << std::endl;

	const auto& field = magdyn.GetExternalField();
	std::cout << "\tMagnetic field magnitude: " << field.mag << "." << std::endl;
	std::cout << "\tMagnetic field direction: " << field.dir << "." << std::endl;
	if(field.align_spins)
		std::cout << "\tAligning spins to field." << std::endl;
	else
		std::cout << "\tNot aligning spins to field." << std::endl;

	std::cout << "\tGround state energy: " << magdyn.CalcGroundStateEnergy()
		<< " meV." << std::endl;

	// get output stream for results
	std::ostream* postr = &std::cout;
	std::unique_ptr<std::ofstream> ofstr;

	if(results_file != "")
	{
		ofstr = std::make_unique<std::ofstream>(results_file);
		postr = ofstr.get();
	}
	else
	{
		std::cerr << "Warning: No output file given, using standard output."
			<< std::endl;
	}


	// get the configuration options
	pt::ptree root_node;
	std::ifstream ifstr{model_file};
	pt::read_xml(ifstr, root_node);
	const auto &magdyn_node = root_node.get_child("magdyn");

	t_real h_start =  magdyn_node.get<t_real>("config.h_start", 0.);
	t_real k_start = magdyn_node.get<t_real>("config.k_start", 0.);
	t_real l_start = magdyn_node.get<t_real>("config.l_start", 0.);
	t_real h_end = magdyn_node.get<t_real>("config.h_end", 1.);
	t_real k_end = magdyn_node.get<t_real>("config.k_end", 0.);
	t_real l_end = magdyn_node.get<t_real>("config.l_end", 0.);

	if(num_Q_pts == 0)
		num_Q_pts = magdyn_node.get<t_size>("config.num_Q_points", 128);

	// get the override options
	if(Qi.size() == 3)
	{
		h_start = Qi[0];
		k_start = Qi[1];
		l_start = Qi[2];
	}
	// get the override options
	if(Qf.size() == 3)
	{
		h_end = Qf[0];
		k_end = Qf[1];
		l_end = Qf[2];
	}


	// calculate the dispersion
	std::cout << "\nCalculating dispersion from"
		<< " Q_i = (" << h_start << ", " << k_start << ", " << l_start << ") to"
		<< " Q_f = (" << h_end << ", " << k_end << ", " << l_end << ")"
		<< " in " << num_Q_pts << " steps and " << g_num_threads << " threads..."
		<< std::endl;

	magdyn.SaveDispersion(*postr,
		h_start, k_start, l_start,
		h_end, k_end, l_end,
		num_Q_pts, g_num_threads);

	if(results_file != "")
		std::cout << "Wrote results to \"" << results_file << "\"." << std::endl;

	return 0;
}



static inline void set_locales()
{
#ifndef DONT_USE_QT
	tl2::set_locales();
#else
	std::ios_base::sync_with_stdio(false);
	::setlocale(LC_ALL, "C");
	std::locale::global(std::locale("C"));
#endif
}



int main(int argc, char** argv)
{
	try
	{
		set_locales();

		bool show_help = false;
		bool use_cli = false;
		bool show_timing = false;

#ifdef DONT_USE_QT
		use_cli = true;
#endif

		t_size num_Q_pts = 0;
		std::string model_file, results_file;

		args::options_description arg_descr("Takin/Magdyn arguments");
		arg_descr.add_options()
			("help,h", args::bool_switch(&show_help), "show help")
#ifndef DONT_USE_QT
			("cli,c", args::bool_switch(&use_cli), "use command-line interface")
#endif
			("input,i", args::value(&model_file), "input magnetic model file (.magdyn)")
			("output,o", args::value(&results_file), "output results file (in cli mode)")
			("timing", args::bool_switch(&show_timing), "show time needed for calculation")
			("threads,t", args::value(&g_num_threads), "number of threads for calculation")
			("points,p", args::value(&num_Q_pts), "number of Q points")
			("hi", args::value<t_real>(), "initial h coordinate")
			("ki", args::value<t_real>(), "initial k coordinate")
			("li", args::value<t_real>(), "initial l coordinate")
			("hf", args::value<t_real>(), "final h coordinate")
			("kf", args::value<t_real>(), "final k coordinate")
			("lf", args::value<t_real>(), "final l coordinate");

		args::positional_options_description posarg_descr;
		posarg_descr.add("input", 1);

		auto argparser = args::command_line_parser{argc, argv};
		argparser.options(arg_descr);
		argparser.positional(posarg_descr);
		argparser.allow_unregistered();
		auto parsedArgs = argparser.run();

		args::variables_map mapArgs;
		args::store(parsedArgs, mapArgs);
		args::notify(mapArgs);

		if(show_help)
		{
			std::cout << "This is Takin/Magdyn by Tobias Weber <tweber@ill.fr>.\n\n"
				<< arg_descr
				<< R"BLOCK(
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
)BLOCK";
			std::cout << std::endl;
			return 0;
		}


		// get Qi and Qf override values
		auto iterHi = mapArgs.find("hi");
		auto iterKi = mapArgs.find("ki");
		auto iterLi = mapArgs.find("li");
		auto iterHf = mapArgs.find("hf");
		auto iterKf = mapArgs.find("kf");
		auto iterLf = mapArgs.find("lf");

		t_vec_real Qi, Qf;
		if(iterHi != mapArgs.end() &&iterKi != mapArgs.end() && iterLi != mapArgs.end())
		{
			Qi = tl2::create<t_vec_real>({
				iterHi->second.as<t_real>(),
				iterKi->second.as<t_real>(),
				iterLi->second.as<t_real>()
			});
		}
		if(iterHf != mapArgs.end() &&iterKf != mapArgs.end() && iterLf != mapArgs.end())
		{
			Qf = tl2::create<t_vec_real>({
				iterHf->second.as<t_real>(),
				iterKf->second.as<t_real>(),
				iterLf->second.as<t_real>()
			});
		}


		// either start the cli or the gui program
		std::unique_ptr<tl2::Stopwatch<t_real>> stopwatch;

		if(show_timing)
		{
			stopwatch = std::make_unique<tl2::Stopwatch<t_real>>();
			stopwatch->start();
		}

		int ret = 0;
		if(use_cli)
			ret = cli_main(model_file, results_file, Qi, Qf, num_Q_pts);
		else
			ret = gui_main(argc, argv, model_file, Qi, Qf, num_Q_pts);

		if(show_timing)
		{
			stopwatch->stop();
			std::cout << "\n================================================================================\n"
				<< "Magdyn start time: " << stopwatch->GetStartTimeStr() << ".\n"
				<< "Magdyn stop time:  " << stopwatch->GetStopTimeStr() << ".\n"
				<< "Elapsed time:      " << stopwatch->GetDur() << " s.\n"
				<< "================================================================================"
				<< std::endl;
		}

		return ret;
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
