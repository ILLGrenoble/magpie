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
#include <stdlib.h>

#include "magpie_c.h"


int main(int argc, char** argv)
{
	if(argc < 2)
	{
		fprintf(stderr, "Please give a magpie model file.\n");
		return -1;
	}

	const char *model_file = argv[1];
	const char *disp_file = "disp.dat";


	// initialise magpie
	void* mag = magpie_create();
	if(!mag)
	{
		fprintf(stderr, "Cannot initialise magpie.\n");
		return -2;
	}


	// load a magnetic model
	fprintf(stderr, "Loading model file \"%s\"...\n", model_file);

	if(!magpie_load(mag, model_file))
	{
		fprintf(stderr, "Cannot load magpie model file \"%s\".\n", model_file);
		return -3;
	}


	// directly save a dispersion branch to a text data file
	t_magpie_real h_step = 0.05;
	fprintf(stderr, "Saving dispersion branch to file \"%s\"...\n", disp_file);

	if(!magpie_save_dispersion(mag,
		disp_file,
		h_step, 0., 0.,
		1. - h_step, 0., 0.,
		128))
	{
		fprintf(stderr, "Cannot save dispersion branch to file  \"%s\".\n", disp_file);
		return -4;
	}


	// clean up magpie
	magpie_free(mag);

	return 0;
}
