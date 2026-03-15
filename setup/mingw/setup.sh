#!/bin/bash
#
# installs the mingw packages
# @author Tobias Weber <tweber@ill.fr>
# @date 15-mar-2026
# @license GPLv3, see 'LICENSE' file
#
# -----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                     Grenoble, France).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
#

# install mingw libraries
dnf install mingw64-cmake mingw64-gcc-c++ mingw64-gcc-gfortran \
	mingw64-boost mingw64-boost-static \
	mingw64-qt6* 

# compile external libraries
pushd ../externals
./build_gemmi.sh --mingw
./build_lapacke.sh --mingw
./build_minuit.sh --mingw
./build_qcp.sh --mingw
./build_qhull.sh --mingw
popd
