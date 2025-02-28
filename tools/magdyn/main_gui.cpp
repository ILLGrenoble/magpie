/**
 * magnetic dynamics -- gui entry point
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

#include "magdyn.h"
#include "tlibs2/libs/qt/gl.h"
#include "tlibs2/libs/qt/helper.h"
#include "tlibs2/libs/algos.h"

#include <QtWidgets/QApplication>

#include <iostream>
#include <memory>



/**
 * starts the gui program
 */
int gui_main(int argc, char** argv, const std::string& model_file,
	const t_vec_real& Qi, const t_vec_real& Qf, t_size num_Q_pts)
{
	tl2::set_gl_format(1, _GL_MAJ_VER, _GL_MIN_VER, 8);

	// application
	QApplication::addLibraryPath(QString(".") + QDir::separator() + "qtplugins");
	auto app = std::make_unique<QApplication>(argc, argv);

	// main window
	auto magdyn = std::make_unique<MagDynDlg>(nullptr);
	magdyn->show();

	// if a configuration file is given, load it
	if(model_file != "")
	{
		QString qfile = model_file.c_str();
		if(magdyn->Load(qfile, false))
			magdyn->SetCurrentFileAndDir(qfile);
	}

	// override the dispersion branch to plot
	if(Qi.size() == 3 && Qf.size() == 3)
		magdyn->SetCoordinates(Qi, Qf, false);
	if(num_Q_pts > 0)
		magdyn->SetNumQPoints(num_Q_pts);

	magdyn->CalcDispersion();
	magdyn->CalcHamiltonian();

	return app->exec();
}
