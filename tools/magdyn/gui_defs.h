/**
 * magnetic dynamics -- definitions for gui
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#ifndef __MAGDYN_GUI_DEFS_H__
#define __MAGDYN_GUI_DEFS_H__

#include "defs.h"



/**
 * user data for table items
 */
enum : int
{
	DATA_ROLE_NAME = Qt::UserRole + 0,
	//DATA_ROLE_IDX  = Qt::UserRole + 1,  // TODO: index used in kernel
};



/**
 * export file types
 */
enum : int
{
	EXPORT_HDF5 = 0,
	EXPORT_GRID = 1,
	EXPORT_TEXT = 2,
};



/**
 * columns of the sites table
 */
enum : int
{
	COL_SITE_NAME = 0,                                 // label
	COL_SITE_POS_X, COL_SITE_POS_Y, COL_SITE_POS_Z,    // site position
	COL_SITE_SYM_IDX,                                  // symmetry index
	COL_SITE_SPIN_X, COL_SITE_SPIN_Y, COL_SITE_SPIN_Z, // spin direction
	COL_SITE_SPIN_MAG,                                 // spin magnitude
	COL_SITE_RGB,                                      // colour
	COL_SITE_SPIN_ORTHO_X, COL_SITE_SPIN_ORTHO_Y, COL_SITE_SPIN_ORTHO_Z,

	NUM_SITE_COLS
};



/**
 * columns of the exchange terms table
 */
enum : int
{
	COL_XCH_NAME = 0,                                // label
	COL_XCH_ATOM1_IDX, COL_XCH_ATOM2_IDX,            // site names or indices
	COL_XCH_DIST_X, COL_XCH_DIST_Y, COL_XCH_DIST_Z,  // unit cell distance
	COL_XCH_SYM_IDX,                                 // symmetry index
	COL_XCH_INTERACTION,                             // isotropic exchange interaction
	COL_XCH_DMI_X, COL_XCH_DMI_Y, COL_XCH_DMI_Z,     // antisymmetric DMI
	COL_XCH_RGB,                                     // colour
	COL_XCH_GEN_XX, COL_XCH_GEN_XY, COL_XCH_GEN_XZ,  //
	COL_XCH_GEN_YX, COL_XCH_GEN_YY, COL_XCH_GEN_YZ,  // general interaction matrix
	COL_XCH_GEN_ZX, COL_XCH_GEN_ZY, COL_XCH_GEN_ZZ,  //

	NUM_XCH_COLS
};



/**
 * columns of the variables table
 */
enum : int
{
	COL_VARS_NAME = 0,
	COL_VARS_VALUE_REAL, COL_VARS_VALUE_IMAG,

	NUM_VARS_COLS
};



/**
 * columns of the table with saved fields
 */
enum : int
{
	COL_FIELD_H = 0, COL_FIELD_K, COL_FIELD_L,
	COL_FIELD_MAG,

	NUM_FIELD_COLS
};



/**
 * columns of the table with saved Q coordinates
 */
enum : int
{
	COL_COORD_NAME = 0,
	COL_COORD_H, COL_COORD_K, COL_COORD_L,

	NUM_COORD_COLS
};


#endif
