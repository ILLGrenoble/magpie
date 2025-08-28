#!/bin/bash
#
# builds external libraries
# @author Tobias Weber <tweber@ill.fr>
# @date 28-august-2025
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magpie
# Copyright (C) 2022-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
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
# ----------------------------------------------------------------------------
#

TMPDIR=tmp

rm -rf "${TMPDIR}"
mkdir -p "${TMPDIR}"
pushd "${TMPDIR}"


if ! ../setup/externals/build_qhull.sh; then
	echo -e "Error: Could not build QHull."
	exit -1
fi

if ! ../setup/externals/build_gemmi.sh; then
echo -e "Error: Could not build Gemmi."
	exit -1
fi

cp -v ../setup/externals/CMakeLists_qcp.txt .
if ! ../setup/externals/build_qcp.sh; then
	echo -e "Error: Could not build QCustomPlot."
	exit -1
fi


popd
