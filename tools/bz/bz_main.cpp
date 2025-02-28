/**
 * brillouin zone tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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
	#include "bz.h"
	#include "tlibs2/libs/qt/helper.h"

	#include <QtCore/QDir>
	#include <QtWidgets/QApplication>
#endif

#include "bz_lib.h"
#include "bz_conf.h"

#include <iostream>
#include <fstream>


#ifndef DONT_USE_BOOTS_PROGOPTS
	#include <boost/program_options.hpp>
	namespace args = boost::program_options;
#endif



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



/**
 * starts the cli program
 */
static int cli_main(const std::string& cfg_file, const std::string& results_file, bool use_stdin)
{
	try
	{
		BZConfig cfg = load_bz_config(cfg_file, use_stdin);

		BZCalc<t_mat, t_vec, t_real> bzcalc;
		bzcalc.SetEps(g_eps);
		bzcalc.SetSymOps(cfg.symops, false);
		if(cfg.xtal_a && cfg.xtal_b && cfg.xtal_c &&
			cfg.xtal_alpha && cfg.xtal_beta && cfg.xtal_gamma)
		{
			bzcalc.SetCrystal(*cfg.xtal_a, *cfg.xtal_b, *cfg.xtal_c,
				*cfg.xtal_alpha, *cfg.xtal_beta, *cfg.xtal_gamma);
		}
		bzcalc.CalcPeaks(cfg.order ? *cfg.order : 5, true);

		if(!bzcalc.CalcBZ())
		{
			std::cerr << "Error calculating Brillouin zone." << std::endl;
			return -1;
		}

		// get calculated bz
		std::string results = bzcalc.PrintJSON(g_prec);

		if(results_file == "")
		{
			// output results to console
			std::cout << results << std::endl;
		}
		else
		{
			// output results to file
			std::ofstream ofstrResults{results_file};
			ofstrResults << results << std::endl;
		}

		return 0;
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << std::endl;
		return -1;
	}
}



#ifndef DONT_USE_QT
/**
 * starts the gui program
 */
static int gui_main(int argc, char** argv, const std::string& cfg_file, bool use_stdin)
{
	tl2::set_gl_format(1, _GL_MAJ_VER, _GL_MIN_VER, 8);

	// application
	QApplication::addLibraryPath(QString(".") + QDir::separator() + "qtplugins");
	auto app = std::make_unique<QApplication>(argc, argv);

	// main window
	auto dlg = std::make_unique<BZDlg>(nullptr);
	dlg->show();

	// if a configuration file is given, load it
	if(cfg_file != "" || use_stdin)
		dlg->Load(cfg_file.c_str(), use_stdin);

	return app->exec();
}
#endif



#ifndef DONT_USE_BOOTS_PROGOPTS
/**
 * starts the cli or the gui program
 */
int main(int argc, char** argv)
{
	set_locales();

#ifndef DONT_USE_QT
	bool use_cli = false;
#endif
	bool show_help = false;
	bool use_stdin = false;
	t_real eps = -1.;
	std::string cfg_file, results_file;

	args::options_description arg_descr("Takin/BZ arguments");
	arg_descr.add_options()
		("help,h", args::bool_switch(&show_help), "show help")
#ifndef DONT_USE_QT
		("cli,c", args::bool_switch(&use_cli), "use command-line interface")
#endif
		("stdin,s", args::bool_switch(&use_stdin), "load configuration file from standard input")
		("eps,e", args::value(&eps), "set epsilon value")
		("input,i", args::value(&cfg_file), "input configuration file")
		("output,o", args::value(&results_file), "output results file");

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
		std::cout << arg_descr << std::endl;
		return 0;
	}

	if(eps >= 0.)
		set_eps(eps);

#ifndef DONT_USE_QT
	// either start the cli or the gui program
	if(use_cli)
		return cli_main(cfg_file, results_file, use_stdin);
	return gui_main(argc, argv, cfg_file, use_stdin);
#else
	// only start the cli program
	if(argc <= 1)
	{
		std::cout << arg_descr << std::endl;
		return -1;
	}
	return cli_main(cfg_file, results_file, use_stdin);
#endif
}



#else   // DONT_USE_BOOTS_PROGOPTS
/**
 * starts the cli program using a minimal interface
 */
int main(int argc, char** argv)
{
	set_locales();
	std::string cfg_file, results_file;
	bool has_cfg_file = false;

	// find input and output file names
	for(int arg=1; arg<argc; ++arg)
	{
		// ignore switches
		if(argv[arg][0] == '-')
			continue;

		if(!has_cfg_file)
		{
			cfg_file = argv[arg];
			has_cfg_file = true;
		}
		else
		{
			results_file = argv[arg];
			break;
		}
	}

	if(argc < 2 || !has_cfg_file)
	{
		std::cout << "Usage: " << argv[0] << "<input.xml> [output.json]" << std::endl;
		return -1;
	}

	return cli_main(cfg_file, results_file, false);
}
#endif  // DONT_USE_BOOTS_PROGOPTS
