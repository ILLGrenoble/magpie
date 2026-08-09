#!/bin/bash
#
# builds openblas
# @author Tobias Weber <tweber@ill.fr>
# @date 9-aug-2026
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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
TMP_DIR=tmp


BUILD_FOR_MINGW=0
if [ "$1" == "--mingw" ]; then
	BUILD_FOR_MINGW=1
fi


OPENBLAS_VER=0.3.34
OPENBLAS_REMOTE="https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${OPENBLAS_VER}.tar.gz"
OPENBLAS_LOCAL="openblas-${OPENBLAS_VER}"


mkdir -v "${TMP_DIR}"
rm -f "${TMP_DIR}/${OPENBLAS_LOCAL}"


if ! wget ${OPENBLAS_REMOTE}; then
	echo -e "Could not download ${OPENBLAS_REMOTE}."
	exit -1
fi


rm -rf build_openblas
OPENBLAS_LOCAL_TAR=${OPENBLAS_REMOTE##*[/\\]}
mv -v "${OPENBLAS_LOCAL_TAR}" "${TMP_DIR}"
cd "${TMP_DIR}"
tar xzvf "${OPENBLAS_LOCAL_TAR}"


if [ $BUILD_FOR_MINGW -ne 0 ]; then
	mkdir build_openblas
	cd build_openblas
	mingw64-cmake -DCMAKE_BUILD_TYPE=Release -DDYNAMIC_ARCH=TRUE \
		-DBUILD_WITHOUT_LAPACK=FALSE -DBUILD_WITHOUT_LAPACKE=FALSE -DBUILD_LAPACK_DEPRECATED=FALSE \
		-DBUILD_TESTING=FALSE -DBUILD_BENCHMARKS=FALSE \
		-DBUILD_STATIC_LIBS=TRUE -DBUILD_SHARED_LIBS=FALSE \
		"../${OPENBLAS_LOCAL}"
	mingw64-make -j${NUM_CORES} && sudo mingw64-make install/strip
else
	cmake -DCMAKE_BUILD_TYPE=Release -DDYNAMIC_ARCH=TRUE \
		-DBUILD_WITHOUT_LAPACK=FALSE -DBUILD_WITHOUT_LAPACKE=FALSE -DBUILD_LAPACK_DEPRECATED=FALSE \
		-DBUILD_TESTING=FALSE -DBUILD_BENCHMARKS=FALSE \
		-DBUILD_STATIC_LIBS=TRUE -DBUILD_SHARED_LIBS=FALSE \
		-B build_openblas "${OPENBLAS_LOCAL}"
	cmake --build build_openblas --parallel ${NUM_CORES} && sudo cmake --install build_openblas --strip
fi
