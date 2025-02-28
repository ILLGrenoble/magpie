/**
 * tlibs2 -- (C-only) math library test
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2020
 * @license GPLv2 or GPLv3, see 'LICENSE' file
 *
 * gcc -Wall -Wextra -std=c99 -o mathlib mathlib.c ../libs/mathlib.c
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "../libs/mathlib.h"
#include <stdio.h>
#include <stdlib.h>


void tst1()
{
	struct tl2_list* vecs = tl2_lst_create(0);
	struct tl2_list* probs = tl2_lst_create(0);
	vecs->elem = calloc(3, sizeof(double));
	probs->elem = calloc(1, sizeof(double));

	((double*)vecs->elem)[0] = 1.;
	((double*)vecs->elem)[1] = 2.;
	((double*)vecs->elem)[2] = 3.;
	*((double*)probs->elem) = 1.;


	struct tl2_list* vecs1 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs1 = tl2_lst_append(probs, 0);
	vecs1->elem = calloc(3, sizeof(double));
	probs1->elem = calloc(1, sizeof(double));

	((double*)vecs1->elem)[0] = 3.;
	((double*)vecs1->elem)[1] = 2.;
	((double*)vecs1->elem)[2] = 1.;
	*((double*)probs1->elem) = 1.;


	struct tl2_list* vecs2 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs2 = tl2_lst_append(probs, 0);
	vecs2->elem = calloc(3, sizeof(double));
	probs2->elem = calloc(1, sizeof(double));

	((double*)vecs2->elem)[0] = 5.;
	((double*)vecs2->elem)[1] = -2.;
	((double*)vecs2->elem)[2] = 0.5;
	*((double*)probs2->elem) = 0.5;


	double *mean = calloc(3, sizeof(double));
	//tl2_vec_mean(vecs, probs, mean, 3);

	double *cov = calloc(3*3, sizeof(double));
	tl2_covariance(vecs, probs, cov, mean, 3);
	printf("mean: %g %g %g\n", mean[0], mean[1], mean[2]);
	printf("cov:\n\t%g %g %g\n\t%g %g %g\n\t%g %g %g\n",
		cov[0], cov[1], cov[2],
		cov[3], cov[4], cov[5],
		cov[6], cov[7], cov[8]);
	free(mean);
	free(cov);


	tl2_lst_free(vecs);
	tl2_lst_free(probs);
}


void tst2()
{
	struct tl2_list* vecs = tl2_lst_create(0);
	struct tl2_list* probs = tl2_lst_create(0);
	vecs->elem = calloc(4, sizeof(double));
	probs->elem = calloc(1, sizeof(double));

	((double*)vecs->elem)[0] = 1.;
	((double*)vecs->elem)[1] = 2.;
	((double*)vecs->elem)[2] = 3.;
	((double*)vecs->elem)[3] = 4.;
	*((double*)probs->elem) = 1.;


	struct tl2_list* vecs1 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs1 = tl2_lst_append(probs, 0);
	vecs1->elem = calloc(4, sizeof(double));
	probs1->elem = calloc(1, sizeof(double));

	((double*)vecs1->elem)[0] = 0.8;
	((double*)vecs1->elem)[1] = 2.2;
	((double*)vecs1->elem)[2] = 3.1;
	((double*)vecs1->elem)[3] = 3.5;
	*((double*)probs1->elem) = 1.;


	struct tl2_list* vecs2 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs2 = tl2_lst_append(probs, 0);
	vecs2->elem = calloc(4, sizeof(double));
	probs2->elem = calloc(1, sizeof(double));

	((double*)vecs2->elem)[0] = 1.2;
	((double*)vecs2->elem)[1] = 2.5;
	((double*)vecs2->elem)[2] = 2.9;
	((double*)vecs2->elem)[3] = 4.1;
	*((double*)probs2->elem) = 0.5;


	struct tl2_list* vecs3 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs3 = tl2_lst_append(probs, 0);
	vecs3->elem = calloc(4, sizeof(double));
	probs3->elem = calloc(1, sizeof(double));

	((double*)vecs3->elem)[0] = 0.9;
	((double*)vecs3->elem)[1] = 1.9;
	((double*)vecs3->elem)[2] = 3.7;
	((double*)vecs3->elem)[3] = 4.1;
	*((double*)probs3->elem) = 0.5;


	struct tl2_list* vecs4 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs4 = tl2_lst_append(probs, 0);
	vecs4->elem = calloc(4, sizeof(double));
	probs4->elem = calloc(1, sizeof(double));

	((double*)vecs4->elem)[0] = 1.2;
	((double*)vecs4->elem)[1] = 1.9;
	((double*)vecs4->elem)[2] = 4.2;
	((double*)vecs4->elem)[3] = 3.9;
	*((double*)probs4->elem) = 0.5;


	double *cov = calloc(4*4, sizeof(double));
	double *reso = calloc(4*4, sizeof(double));
	tl2_reso(vecs, probs, cov, reso);
	printf("cov:\n\t%g %g %g %g\n\t%g %g %g %g\n\t%g %g %g %g\n\t%g %g %g %g\n",
		cov[0], cov[1], cov[2], cov[3],
		cov[4], cov[5], cov[6], cov[7],
		cov[8], cov[9], cov[10], cov[11],
		cov[12], cov[13], cov[14], cov[15]);
	printf("reso:\n\t%g %g %g %g\n\t%g %g %g %g\n\t%g %g %g %g\n\t%g %g %g %g\n",
		reso[0], reso[1], reso[2], reso[3],
		reso[4], reso[5], reso[6], reso[7],
		reso[8], reso[9], reso[10], reso[11],
		reso[12], reso[13], reso[14], reso[15]);
	free(cov);
	free(reso);


	tl2_lst_free(vecs);
	tl2_lst_free(probs);
}


int main()
{
	tst1();
	tst2();
	return 0;
}
