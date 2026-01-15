/**
 * magnetic dynamics c library interface test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15-january-2026
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

#include <stdio.h>
#include "magpie_c.h"


int main(int argc, char** argv)
{
	if(argc < 2)
	{
		fprintf(stderr, "Please give a magpie model file.\n");
		return -1;
	}

	void* mag = magpie_create();

	if(!magpie_load(mag, argv[1]))
	{
		fprintf(stderr, "Cannot load magpie model file \"%s\".\n", argv[1]);
		return -2;
	}

	magpie_free(mag);
	return 0;
}
