/**
 * space group browser
 * @author Tobias Weber <tweber@ill.fr>
 * @date Apr-2018
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools" project
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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

#include <QtCore/QSettings>
#include <QtWidgets/QApplication>

#include <locale>
#include <memory>

#include "browser.h"
#include "tlibs2/libs/qt/helper.h"


// ----------------------------------------------------------------------------


int main(int argc, char** argv)
{
	QSettings sett("takin", "sgbrowser", nullptr);

	auto app = std::make_unique<QApplication>(argc, argv);
	tl2::set_locales();

	auto dlg = std::make_unique<SgBrowserDlg>(nullptr, &sett);
	dlg->show();

	return app->exec();
}

// ----------------------------------------------------------------------------
