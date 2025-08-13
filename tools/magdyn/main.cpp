/**
 * magnetic dynamics -- entry point
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2025
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
static int gui_main(int, char**, const std::string&, const t_vec_real&, const t_vec_real&, t_size)
{
	std::cerr << "Error: The GUI is not available in this version." << std::endl;
	return -1;
}
#else
extern int gui_main(int argc, char** argv, const std::string& model_file,
	const t_vec_real& Qi, const t_vec_real& Qf, t_size num_Q_pts);
#endif


#define CERR_OPT if(!silent) std::cerr
#define COUT_OPT if(!silent) std::cout


/**
 * starts the cli program
 */
static int cli_main(const std::string& model_file, const std::string& results_file,
	const t_vec_real& Qi, const t_vec_real& Qf, const t_vec_real& Qf2,
	t_size num_Q_pts, t_size num_Q_pts2, const std::string& Qlist_file,
	bool as_py = false, bool as_bin = false, bool calc_weights = true, bool silent = false/*,
	t_real Emin = 1., t_real Emax = -1*/)
{
	using namespace tl2_ops;

	// load model from input file
	t_magdyn magdyn;
	if(magdyn.Load(model_file))
	{
		COUT_OPT << "Loaded magnetic model from file \"" <<
			model_file << "\"." << std::endl;
	}
	else
	{
		CERR_OPT << "Error: Failed loading magnetic model \"" <<
			model_file << "\"." << std::endl;
		return -1;
	}

	// calculation options
	//if(Emax >= Emin)
	//	magdyn.SetFilterEnergies(Emin, Emax);
	magdyn.SetSilent(silent);

	if(!silent)
	{
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
	}

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
		CERR_OPT << "Warning: No output file given, using standard output."
			<< std::endl;
	}


	// get the configuration options
	pt::ptree root_node;
	std::ifstream ifstr{model_file};
	pt::read_xml(ifstr, root_node);
	const auto &magdyn_node = root_node.get_child("magdyn");


	if(Qlist_file != "")
	{
		std::ifstream ifstr_Qs(Qlist_file);
		if(!ifstr_Qs)
		{
			CERR_OPT << "Error: Cannot open Q list file \""
				<< Qlist_file << "\"." << std::endl;
			return -1;
		}

		// calculate the dispersion for the given Q points
		std::vector<t_vec_real> Qs;

		// read one h k l triplet per line
		std::string line;
		std::size_t lineno = 0;
		while(std::getline(ifstr_Qs, line))
		{
			++lineno;
			tl2::trim(line);
			if(line == "")
				continue;  // ignore empty lines
			if(line.size() > 0 && line[0] == '#')
				continue;  // ignore comments

			t_vec_real Q;
			tl2::get_tokens<t_real, std::string, t_vec_real>(line, std::string(" \t;,"), Q);
			if(Q.size() != 3)
			{
				CERR_OPT << "Error: Invalid Q vector length in line " << lineno
					<< " of \"" << Qlist_file
					<< "\", three components, h, k, and l, are required."
					<< std::endl;
				continue;
			}
			Qs.emplace_back(std::move(Q));
		}

		COUT_OPT << "\nCalculating dispersion for custom Q points in "
			<< g_num_threads << " threads..."
			<< std::endl;

		magdyn.SaveDispersion(*postr, Qs, g_num_threads,
			as_py, as_bin, calc_weights);
	}
	else
	{
		// calculate the dispersion along the given Q_i and Q_f positions
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

		bool use_2d_map = false;
		if(Qf2.size() == 3 && num_Q_pts2 > 0)
			use_2d_map = true;

		if(!use_2d_map)
		{
			// 1d dispersion
			COUT_OPT << "\nCalculating 1d dispersion from"
				<< " Q_i = (" << h_start << ", " << k_start << ", " << l_start << ") to"
				<< " Q_f = (" << h_end << ", " << k_end << ", " << l_end << ")"
				<< " in " << num_Q_pts << " steps and " << g_num_threads << " threads..."
				<< std::endl;

			magdyn.SaveDispersion(*postr,
				h_start, k_start, l_start,
				h_end, k_end, l_end,
				num_Q_pts, g_num_threads,
				as_py, as_bin, calc_weights);
		}
		else
		{
			// 2d dispersion
			// generate Q points
			std::vector<t_vec_real> Qs;
			Qs.reserve(num_Q_pts * num_Q_pts2);

			t_vec_real Q0 = tl2::create<t_vec_real>({ h_start, k_start, l_start });
			t_vec_real Qdir1 = tl2::create<t_vec_real>({ h_end, k_end, l_end }) - Q0;
			t_vec_real Qdir2 = Qf2 - Q0;
			Qdir1 /= t_real(num_Q_pts);
			Qdir2 /= t_real(num_Q_pts2);

			//tl2::Stopwatch<t_real> stopwatch;
			//stopwatch.start();
			for(t_size i = 0; i < num_Q_pts; ++i)
			for(t_size j = 0; j < num_Q_pts2; ++j)
			{
				t_vec_real Q = Q0 + Qdir1 * t_real(i) + Qdir2 * t_real(j);
				Qs.emplace_back(std::move(Q));
			}
			//stopwatch.stop();
			//COUT_OPT << "Generated " << num_Q_pts * num_Q_pts2 << " Q points in "
			//	<< stopwatch.GetDur() << " s." << std::endl;

			COUT_OPT << "\nCalculating 2d dispersion from"
				<< " Q_i = (" << h_start << ", " << k_start << ", " << l_start << ") to"
				<< " Q_f = (" << h_end << ", " << k_end << ", " << l_end << ") and"
				<< " Q_f2 = (" << Qf2[0] << ", " << Qf2[1] << ", " << Qf2[2] << ")"
				<< " in " << num_Q_pts << "*" << num_Q_pts2 << " steps and "
				<< g_num_threads << " threads..."
				<< std::endl;

			magdyn.SaveDispersion(*postr, Qs, g_num_threads,
				as_py, as_bin, calc_weights);
		}
	}

	if(results_file != "")
		COUT_OPT << "Wrote results to \"" << results_file << "\"." << std::endl;

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
		bool silent = false;
		bool use_cli = false;
		bool show_timing = false;
		bool as_py = false;
		bool as_bin = false;
		bool no_weights = false;
		//t_real Emin = 1., Emax = -1.;

#ifdef DONT_USE_QT
		use_cli = true;
#endif

		t_size num_Q_pts = 0, num_Q_pts2 = 0;
		std::string model_file, results_file;
		std::string Qlist_file;

		args::options_description arg_descr("Magpie arguments");
		arg_descr.add_options()
			("help,h", args::bool_switch(&show_help), "show help")
#ifndef DONT_USE_QT
			("cli,c", args::bool_switch(&use_cli), "use command-line interface")
#endif
			("silent", args::bool_switch(&silent), "disable console output")
			("input,i", args::value(&model_file), "input magnetic model file (.magpie)")
			("output,o", args::value(&results_file), "output results file (in cli mode)")
			("qlist", args::value(&Qlist_file), "input file containing Q points")
			("py", args::bool_switch(&as_py), "create py script")
			("bin", args::bool_switch(&as_bin), "create binary output")
			("no_weights", args::bool_switch(&no_weights), "don't calculate spectral weights")
			("timing", args::bool_switch(&show_timing), "show time needed for calculation")
			("threads,t", args::value(&g_num_threads), "number of threads for calculation")
			("points,p", args::value(&num_Q_pts), "number of Q points")
			("points2", args::value(&num_Q_pts2), "number of Q points in second direction")
			("hi", args::value<t_real>(), "initial h coordinate")
			("ki", args::value<t_real>(), "initial k coordinate")
			("li", args::value<t_real>(), "initial l coordinate")
			("hf", args::value<t_real>(), "final h coordinate")
			("kf", args::value<t_real>(), "final k coordinate")
			("lf", args::value<t_real>(), "final l coordinate")
			("hf2", args::value<t_real>(), "second final h coordinate")
			("kf2", args::value<t_real>(), "second final k coordinate")
			("lf2", args::value<t_real>(), "second final l coordinate")
			/*("Emin", args::value<t_real>(&Emin), "minimum energy to calculate")
			("Emax", args::value<t_real>(&Emax), "maximum energy to calculate")*/;

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
			std::cout << "This is Magpie " << MAGPIE_VER
				<< ", written 2022-2025 by Tobias Weber <tweber@ill.fr>.\n\n"
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
		auto iterHf2 = mapArgs.find("hf2");
		auto iterKf2 = mapArgs.find("kf2");
		auto iterLf2 = mapArgs.find("lf2");

		t_vec_real Qi, Qf, Qf2;
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
		if(iterHf2 != mapArgs.end() &&iterKf2 != mapArgs.end() && iterLf2 != mapArgs.end())
		{
			Qf2 = tl2::create<t_vec_real>({
				iterHf2->second.as<t_real>(),
				iterKf2->second.as<t_real>(),
				iterLf2->second.as<t_real>()
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
		{
			if(model_file == "")
			{
				CERR_OPT << "Error: No magnetic model given.\n\n"
					<< arg_descr << std::endl;
				return -1;
			}

			ret = cli_main(model_file, results_file, Qi, Qf, Qf2,
				num_Q_pts, num_Q_pts2, Qlist_file,
				as_py, as_bin, !no_weights, silent/*,
				Emin, Emax*/);
		}
		else
		{
			ret = gui_main(argc, argv, model_file, Qi, Qf, num_Q_pts);
		}

		if(show_timing)
		{
			stopwatch->stop();
			std::cout << "\n================================================================================\n"
				<< "Magpie start time: " << stopwatch->GetStartTimeStr() << ".\n"
				<< "Magpie stop time:  " << stopwatch->GetStopTimeStr() << ".\n"
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
