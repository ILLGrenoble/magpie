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

# create temporary directory for building
TMPDIR=tmp
rm -rf "${TMPDIR}"
mkdir -p "${TMPDIR}"
pushd "${TMPDIR}"


#MAGPIE_ROOT=$(dirname $0)
MAGPIE_ROOT=..
echo -e "Magpie root directory: ${MAGPIE_ROOT}"


# parse command-line arguments
skip_takin_libs=0
for((arg_idx=1; arg_idx<=$#; ++arg_idx)); do
	if [[ "${!arg_idx}" == "--skip-takin-libs" ]]; then
		skip_takin_libs=1
		echo -e "Skipping libraries already built for Takin."
	fi
done


# build libraries
if ! ${MAGPIE_ROOT}/setup/externals/build_qhull.sh; then
	echo -e "Error: Could not build QHull."
	exit -1
fi

if ! ${MAGPIE_ROOT}/setup/externals/build_gemmi.sh; then
echo -e "Error: Could not build Gemmi."
	exit -1
fi

cp -v ${MAGPIE_ROOT}/setup/externals/CMakeLists_qcp.txt .
if ! ${MAGPIE_ROOT}/setup/externals/build_qcp.sh; then
	echo -e "Error: Could not build QCustomPlot."
	exit -1
fi


# build libraries that are also needed by takin
if [ $skip_takin_libs == 0 ]; then
	if ! ${MAGPIE_ROOT}/setup/externals/build_lapacke.sh; then
	echo -e "Error: Could not build Lapack(e)."
		exit -1
	fi

	if ! ${MAGPIE_ROOT}/setup/externals/build_minuit.sh; then
	echo -e "Error: Could not build Minuit."
		exit -1
	fi
fi


echo -e "Building external libraries for Magpie finished successfully."
popd
