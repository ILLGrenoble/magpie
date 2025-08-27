#!/bin/bash
#
# builds qcp
# @author Tobias Weber <tweber@ill.fr>
# @date oct-2024
# @note forked from TAS-Paths: https://github.com/ILLGrenoble/taspaths
# @note thanks to J. Krüger for cleaning up this script
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

NUM_CORES=$(nproc)


BUILD_FOR_MINGW=0
if [ "$1" == "--mingw" ]; then
	BUILD_FOR_MINGW=1
fi


QCP_REMOTE=https://www.qcustomplot.com/release/2.1.1/QCustomPlot-source.tar.gz
QCP_LOCAL=${QCP_REMOTE##*[/\\]}


rm -f "${QCP_LOCAL}"


if ! wget ${QCP_REMOTE}; then
	echo -e "Could not download ${QCP_REMOTE}."
	exit -1
fi


if ! tar xzvf ${QCP_LOCAL}
then
	echo -e "Error extracting QCustomPlot.";
	exit -1;
fi


cp -v CMakeLists_qcp.txt qcustomplot-source/CMakeLists.txt
cd qcustomplot-source


if [ $BUILD_FOR_MINGW -ne 0 ]; then
	mkdir build_lib && cd build_lib
	mingw64-cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=True ..
	mingw64-make -j${NUM_CORES} && sudo mingw64-make install/strip
else
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=True -B build_lib .
	cmake --build build_lib --parallel ${NUM_CORES} && sudo cmake --install build_lib --strip
fi
