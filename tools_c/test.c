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

#define CALC_DISP 1
#define SAVE_DISP 0


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


	t_magpie_real h_step = 0.05;

#if CALC_DISP
	// calculate the energies and weights for the given momentum transfers
	unsigned int max_branches = 2*magpie_site_count(mag);
	t_magpie_real *Es = malloc(max_branches * sizeof(t_magpie_real));
	t_magpie_real *ws = malloc(max_branches * sizeof(t_magpie_real));

	t_magpie_real k = 0., l = 0.;

	printf("%10s %10s %10s %10s %10s\n", "h (rlu)", "k (rlu)", "l (rlu)", "E (meV)", "S (a.u.)");
	for(t_magpie_real h = h_step; h < 1.; h += h_step)
	{
		unsigned int num_branches = magpie_calc_energies(mag, h, k, l, Es, ws);
		for(unsigned int branch_idx = 0; branch_idx < num_branches; ++branch_idx)
			printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n", h, k, l, Es[branch_idx], ws[branch_idx]);
	}

	free(Es);
	free(ws);
#endif


#if SAVE_DISP
	// directly save a dispersion branch to a text data file
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
#endif


	// clean up magpie
	magpie_free(mag);

	return 0;
}
