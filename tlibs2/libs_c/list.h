/**
 * tlibs2 -- (C-only) linked list
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2020
 * @note forked on 20-nov-2020 from the runtime library of my private "matrix_calc" project (https://github.com/t-weber/matrix_calc/blob/master/src/runtime.c).
 * @license GPLv2 or GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * matrix_calc
 * Copyright (C) 2020       Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 or version 3 of the License.
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

#ifndef __TLIBS2_C_LIST_H__
#define __TLIBS2_C_LIST_H__


// ----------------------------------------------------------------------------
// linked list
// ----------------------------------------------------------------------------
struct tl2_list
{
	struct tl2_list *next;
	void *elem;
};

extern struct tl2_list* tl2_lst_create(void *elem);
extern struct tl2_list* tl2_lst_append(struct tl2_list *lst, void *elem);
extern void tl2_lst_remove(struct tl2_list *lst, void *elem);
extern void tl2_lst_free(struct tl2_list *lst);
// ----------------------------------------------------------------------------


#endif
