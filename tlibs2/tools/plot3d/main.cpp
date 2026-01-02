/**
 * 3d plotter -- entry point
 * @author Tobias Weber <tweber@ill.fr>
 * @date january 2026
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include <QtCore/QDir>
#include <QtCore/QSettings>
#include <QtWidgets/QApplication>

#include <iostream>
#include <memory>

#include "plot3d.h"
#include "defs.h"

#include "tlibs2/libs/qt/gl.h"
#include "tlibs2/libs/qt/helper.h"
#include "tlibs2/libs/algos.h"



int main(int argc, char** argv)
{
	tl2::set_gl_format(1, _GL_MAJ_VER, _GL_MIN_VER, 8);

	// application
	QApplication::addLibraryPath(QString(".") + QDir::separator() + "qtplugins");
	auto app = std::make_unique<QApplication>(argc, argv);

	// re-set locales
	tl2::set_locales();

	// main window
	auto sett = std::make_shared<QSettings>("tlibs", "plot3d");
	auto plt = std::make_unique<Plot3DDlg>(nullptr, sett);
	plt->show();

	return app->exec();
}
