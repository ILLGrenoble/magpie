#!/bin/bash
#
# creates a mingw distribution
# @author Tobias Weber <tweber@ill.fr>
# @date 2016, may-2021, nov-2025
# @license GPLv3, see 'LICENSE' file
#
# -----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2025  Tobias WEBER (Institut Laue-Langevin (ILL),
#                     Grenoble, France).
# TAS-Paths
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

APPDIRNAME="magpie_mingw"


# icon
APPICON="magpie.svg"
APPICON_ICO="${APPICON%\.svg}.ico"


# mingw directories
MINGW_ROOT=/usr/x86_64-w64-mingw32/sys-root/mingw
MINGW_QT_ROOT=/usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt6


# third-party libraries
EXT_LIBS=( \
	Qt6Core.dll Qt6Gui.dll Qt6Widgets.dll \
	Qt6OpenGL.dll Qt6OpenGLWidgets.dll \
	Qt6DBus.dll Qt6PrintSupport.dll Qt6Svg.dll \
	libboost_program_options-mt-x64.dll \
	libMinuit2.dll libMinuit2Math.dll \
	libqhull_r.dll libqcustomplot.dll \
	liblapack.dll liblapacke.dll libblas.dll libgfortran-5.dll libquadmath-0.dll \
	libstdc++-6.dll libwinpthread-1.dll libglib-2.0-0.dll libgcc_s_seh-1.dll \
	libbz2-1.dll zlib1.dll libpng16-16.dll \
	libfreetype-6.dll libfontconfig-1.dll \
	libpcre2-16-0.dll libpcre2-8-0.dll libexpat-1.dll \
	libharfbuzz-0.dll iconv.dll libintl-8.dll \
	icui18n76.dll icuuc76.dll icudata76.dll \
)


# qt plugins
QT_PLUGINS=( \
	platforms/*.dll \
	styles/*.dll \
)


# create program directory
cd "$(dirname $0)/../.."
mkdir -p "${APPDIRNAME}"


#
# create a png icon with the specified size out of an svg
#
svg_to_png() {
	local ICON_SVG="$1"
	local ICON_PNG="${APPICON%\.svg}.png"
	local ICON_SIZE="$2"

	echo -e "${ICON_SVG} -> ${ICON_PNG} (size: ${ICON_SIZE}x${ICON_SIZE})..."
	magick -resize "${ICON_SIZE}x${ICON_SIZE}" \
		-antialias -channel rgba \
		-background "#ffffff00" -alpha background \
		"${ICON_SVG}" "${ICON_PNG}"

	svg_to_png_result="${ICON_PNG}"
}


# create the application icon
svg_to_png "${APPICON}" 128
APPICON_PNG="${svg_to_png_result}"

echo -e "${APPICON_PNG} -> ${APPICON_ICO}..."
magick "${APPICON_PNG}" "${APPICON_ICO}"


# copy program files
cp -v  tools/magdyn/build/*.exe  "${APPDIRNAME}/"
cp -rv tools/magdyn/test         "${APPDIRNAME}/"
cp -rv tools/magdyn/scripts      "${APPDIRNAME}/"
cp -v  AUTHORS                   "${APPDIRNAME}/"
cp -v  LICENSE                   "${APPDIRNAME}/"
cp -rv 3rdparty_licenses         "${APPDIRNAME}/"


# copy third-party libraries
for THELIB in ${EXT_LIBS[@]}; do
	cp -v ${MINGW_ROOT}/bin/${THELIB}  "${APPDIRNAME}/"
done


# copy qt plugins
for THELIB in ${QT_PLUGINS[@]}; do
	LIBDIRNAME=$(dirname ${THELIB})
	mkdir -p "${APPDIRNAME}/qt_plugins/${LIBDIRNAME}"
	cp -v "${MINGW_QT_ROOT}/plugins/${THELIB}" "${APPDIRNAME}/qt_plugins/${LIBDIRNAME}"
done


# create qt config file
echo -e "[Paths]\n\tPlugins = qt_plugins\n" > "${APPDIRNAME}/qt.conf"


# strip binaries
find ${APPDIRNAME} -type f \( -name "*.exe" -o -name "*.dll" \) -exec strip -v {} \;
