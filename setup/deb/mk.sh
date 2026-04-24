#!/bin/bash
#
# creates a debian package
# @author Tobias Weber <tweber@ill.fr>
# @date 2016, 23-may-2021, 11-mar-2026
# @license GPLv3
#
# -----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# TAS-Paths (part of the Takin software suite)
# Copyright (C) 2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

# options
create_appdir=1
create_deb=1

# defines
APPNAME="magpie"
APPDIRNAME="tmp/${APPNAME}"
APPDEBNAME="${APPNAME}.deb"

# install directories
BINDIR=/usr/local/bin
SHAREDIR=/usr/local/share
PY_DISTDIR=/usr/local/lib/python3.12/dist-packages


#
# create application directories
#
if [ $create_appdir -ne 0 ]; then
	# clean up old dir
	rm -rfv "${APPDIRNAME}"

	# directories
	mkdir -pv ${APPDIRNAME}${BINDIR}
	mkdir -pv ${APPDIRNAME}${SHAREDIR}/${APPNAME}/res
	mkdir -pv ${APPDIRNAME}${PY_DISTDIR}
	mkdir -pv ${APPDIRNAME}/usr/share/applications
	mkdir -pv ${APPDIRNAME}/DEBIAN

	# package control file
	echo -e "Package: ${APPNAME}\nVersion: 0.9.6" >         ${APPDIRNAME}/DEBIAN/control
	echo -e "Architecture: $(dpkg --print-architecture)" >> ${APPDIRNAME}/DEBIAN/control
	echo -e "Section: base\nPriority: optional" >>          ${APPDIRNAME}/DEBIAN/control
	echo -e "Description: A Graphical Magnon Software" >>   ${APPDIRNAME}/DEBIAN/control
	echo -e "Maintainer: tweber@ill.fr" >>                  ${APPDIRNAME}/DEBIAN/control

	if [ "$1" == "resolute" ]; then
		echo -e "Choosing debendencies for Resolute..."

		echo -e "Depends:" \
			"libstdc++6 (>=10.0.0)," \
			"libboost-program-options1.90.0 (>=1.90.0)," "libboost-container1.90.0 (>=1.90.0)," \
			"libqhull-r8.0 (>=2020.2)," "libqhullcpp8.0 (>=2020.2)," \
			"libqt6core6 (>=6.0.0)," "libqt6gui6 (>=6.0.0)," \
			"libqt6widgets6 (>=6.0.0)," "libqt6svg6 (>=6.0.0)," \
			"libqt6opengl6 (>=6.0.0)," "libqt6openglwidgets6 (>=6.0.0)," \
			"libqt6printsupport6 (>=6.0.0),"\
			"libqcustomplot2.1 (>=2.0.0)," \
			"libhdf5-310," "libhdf5-cpp-310," \
			"python3 (>=3.10.0)," "python3-numpy," "python3-matplotlib," \
			"liblapacke (>=3.11)," \
			"libgmp10 (>=2:6.2)," \
			"libopengl0 (>=1.3.0)\n" \
				>> ${APPDIRNAME}/DEBIAN/control
	elif [ "$1" == "noble" ] || [  "$1" == "" ]; then
		echo -e "Choosing debendencies for Noble..."

		echo -e "Depends:" \
			"libstdc++6 (>=10.0.0)," \
			"libboost-program-options1.83.0 (>=1.83.0)," "libboost-container1.83.0 (>=1.83.0)," \
			"libqhull-r8.0 (>=2020.2)," "libqhullcpp8.0 (>=2020.2)," \
			"libqt6core6t64 (>=6.0.0)," "libqt6gui6t64 (>=6.0.0)," \
			"libqt6widgets6t64 (>=6.0.0)," "libqt6svg6 (>=6.0.0)," \
			"libqt6opengl6t64 (>=6.0.0)," "libqt6openglwidgets6t64 (>=6.0.0)," \
			"libqt6printsupport6t64 (>=6.0.0),"\
			"libqcustomplot2.1 (>=2.0.0)," \
			"libhdf5-103-1t64," "libhdf5-cpp-103-1t64," \
			"python3 (>=3.10.0)," "python3-numpy," "python3-matplotlib," \
			"liblapacke (>=3.11)," \
			"libgmp10 (>=2:6.2)," \
			"libopengl0 (>=1.3.0)\n" \
				>> ${APPDIRNAME}/DEBIAN/control
	elif [ "$1" == "jammy" ]; then
		echo -e "Choosing debendencies for Jammy..."

		echo -e "Depends:" \
			"libstdc++6 (>=10.0.0)," \
			"libboost-program-options1.74.0 (>=1.74.0)," "libboost-container1.74.0 (>=1.74.0)," \
			"libqhull-r8.0 (>=2020.2)," "libqhullcpp8.0 (>=2020.2)," \
			"libqt6core6 (>=6.0.0)," "libqt6gui6 (>=6.0.0)," \
			"libqt6widgets6 (>=6.0.0)," "libqt6svg6 (>=6.0.0)," \
			"libqt6opengl6 (>=6.0.0)," "libqt6openglwidgets6 (>=6.0.0)," \
			"libqt6printsupport6 (>=6.0.0)," \
			"libqcustomplot2.0 (>=2.0.0)," \
			"libhdf5-103-1," "libhdf5-cpp-103-1," \
			"python3 (>=3.10.0)," "python3-numpy," "python3-matplotlib," \
			"libgmp10 (>=2:6.2)," \
			"libopengl0 (>=1.3.0)\n" \
				>> ${APPDIRNAME}/DEBIAN/control
	else
		echo -e "Invalid target system: ${1}."
		exit -1
	fi

	# program files (changing bin name to not interfere with takin's bundled version)
	cp -v  build/magpie/magpie          ${APPDIRNAME}${BINDIR}/themagpie
	cp -v  AUTHORS                      ${APPDIRNAME}${SHAREDIR}/${APPNAME}/
	cp -v  LICENSE                      ${APPDIRNAME}${SHAREDIR}/${APPNAME}/
	cp -v  LICENSES                     ${APPDIRNAME}${SHAREDIR}/${APPNAME}/
	cp -rv res/*                        ${APPDIRNAME}${SHAREDIR}/${APPNAME}/res/
	cp -v  setup/deb/themagpie.desktop  ${APPDIRNAME}/usr/share/applications
	#cp -v  setup/deb/magpie.desktop     ${APPDIRNAME}/usr/share/applications

	# py interface
	cp -v  build/tools_py/magdyn/_magpy.so  ${APPDIRNAME}${PY_DISTDIR}
	cp -v  build/tools_py/magdyn/magpy.py   ${APPDIRNAME}${PY_DISTDIR}

	# cleanups
	chmod a+x  ${APPDIRNAME}${BINDIR}/themagpie
	strip -v   ${APPDIRNAME}${BINDIR}/themagpie
	strip -v   ${APPDIRNAME}${PY_DISTDIR}/_magpy.so
fi


#
# create deb package
#
if [ $create_deb -ne 0 ]; then
	cd tmp
	chmod -R 775 ${APPNAME}
	dpkg-deb -v --root-owner-group --build ${APPNAME}
fi
