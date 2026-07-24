#!/bin/bash
#
# @author Tobias Weber <tweber@ill.fr>
# @date 18-dec-19
# @license GPLv3
#
# ----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

PROGDIR=~/.wine/drive_c/magpie
MINGWDIR=/usr/x86_64-w64-mingw32/sys-root/mingw/bin
MINGWQTDIR=/usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt5

mkdir ${PROGDIR}
mkdir -p ${PROGDIR}/qtplugins/platforms

cp -v build/magpie.exe ${PROGDIR}/

cp -v ${MINGWDIR}/Qt5Core.dll		${PROGDIR}/
cp -v ${MINGWDIR}/Qt5Gui.dll		${PROGDIR}/
cp -v ${MINGWDIR}/Qt5Widgets.dll		${PROGDIR}/
cp -v ${MINGWQTDIR}/plugins/platforms/*.dll	${PROGDIR}/qtplugins/platforms/

cp -v ${MINGWDIR}/libgcc*.dll		${PROGDIR}/
cp -v ${MINGWDIR}/libstdc++-*.dll	${PROGDIR}/
cp -v ${MINGWDIR}/libglib*.dll		${PROGDIR}/
cp -v ${MINGWDIR}/libwinpthread*.dll	${PROGDIR}/

cp -v ${MINGWDIR}/zlib*.dll		${PROGDIR}/
cp -v ${MINGWDIR}/libbz2*.dll		${PROGDIR}/

cp -v ${MINGWDIR}/libfreetype*.dll	${PROGDIR}/
cp -v ${MINGWDIR}/libpng*.dll		${PROGDIR}/

cp -v ${MINGWDIR}/iconv.dll		${PROGDIR}/
cp -v ${MINGWDIR}/libpcre*.dll		${PROGDIR}/
cp -v ${MINGWDIR}/libharfbuzz*.dll	${PROGDIR}/
cp -v ${MINGWDIR}/libintl*.dll		${PROGDIR}/

#echo -e "[paths]\nplugins=qtplugins\n" > ${PROGDIR}/qt.conf

find ${PROGDIR} -name "*.exe" -or -name "*.dll" -exec strip {} \; -print
