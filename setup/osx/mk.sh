#!/bin/bash
#
# creates an app bundle and a dmg file
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2019, apr-2021, mar-2026
# @license GPLv3, see 'LICENSE' file
#
# -----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
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

# settings
create_icon=1
create_appdir=1
create_dmg=1
strip_binaries=1
clean_frameworks=1


# tools
STRIP=/opt/homebrew/Cellar/llvm/22.1.0/bin/llvm-strip


# application name
APPNAME="Magpie"
APPDIRNAME="${APPNAME}.app"
APPDMGNAME="${APPNAME}.dmg"
TMPFILE="${APPNAME}_tmp.dmg"
APPICON="res/magpie.svg"
APPICON_ICNS="${APPICON%\.svg}.icns"


# directories
LOCAL_DIR=/opt/homebrew
QT_PLUGIN_DIR=${LOCAL_DIR}/share/qt/plugins
GCC_LIB_DIR=${LOCAL_DIR}/Cellar/gcc/15.2.0_1/lib/gcc/current/
LOCAL_FRAMEWORKS_DIR=${LOCAL_DIR}/Frameworks
ICU_LIB_DIR=${LOCAL_DIR}/Cellar/icu4c@78/78.2/lib


# qt libraries
declare -a QT_LIBS=(
	QtCore QtGui QtWidgets
	QtOpenGL QtOpenGLWidgets
	QtDBus QtPrintSupport QtSvg
)


# libraries
declare -a LOCAL_LIBS=(
	libboost_program_options.dylib libboost_container.dylib
	libqhull_r.8.0.dylib
	libhdf5.320.dylib libhdf5_cpp.320.dylib
	libdouble-conversion.3.dylib
	libfreetype.6.dylib libharfbuzz.0.dylib
	libpng16.16.dylib libzstd.1.dylib
	libpcre2-8.0.dylib libpcre2-16.0.dylib
	libglib-2.0.0.dylib libgthread-2.0.0.dylib libintl.8.dylib
	libdbus-1.3.dylib libsz.2.dylib libaec.0.dylib
	libb2.1.dylib libmd4c.0.dylib libgraphite2.3.dylib
)


# libraries with full path
declare -a MISC_LIBS=(
	# gcc libraries
	${GCC_LIB_DIR}/libgfortran.5.dylib
	${GCC_LIB_DIR}/libquadmath.0.dylib
	${GCC_LIB_DIR}/libgcc_s.1.1.dylib

	# icu libraries
	${ICU_LIB_DIR}/libicui18n.78.dylib
	${ICU_LIB_DIR}/libicuuc.78.dylib
	${ICU_LIB_DIR}/libicudata.78.dylib
)


# ansi colours
COL_ERR="\033[1;31m"
COL_WARN="\033[1;31m"
COL_NORM="\033[0m"


#
# tests if the given file is a binary image
#
function is_binary()
{
	local binary="$1"

	# does the file exist?
	if [ ! -f $binary ]; then
		return 0
	fi

	# is it a binary?
	if [[ "$(file ${binary})" != *"Mach-O"* ]]; then
		return 0
	fi

	return 1
}


#
# tests if the given binary still has remaining bindings to LOCAL_DIR
#
function check_local_bindings()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "${COL_WARN}Warning: ${cfile} does not exist.${COL_NORM}"
		return
	fi

	local local_binding=$(otool -L ${cfile} | grep --color ${LOCAL_DIR})
	local num_local_bindings=$(echo -e "${local_binding}" | wc -l)

	if [ "${num_local_bindings}" -gt "1" ]; then
		echo -e "${COL_WARN}Warning: Possible local binding remaining, please check \"tmp/local_bindings.txt\"!${COL_NORM}"
	fi

	echo -e "--------------------------------------------------------------------------------" >> tmp/local_bindings.txt
	echo -e "${num_local_bindings} local binding(s) in binary \"${cfile}\":" >> tmp/local_bindings.txt
	echo -e "${local_binding}" >> tmp/local_bindings.txt
	echo -e "--------------------------------------------------------------------------------\n" >> tmp/local_bindings.txt
}


#
# changes a LOCAL_DIR linker path to an @rpath
#
function change_to_rpath()
{
	local binary="$1"

	# local lib paths
	local old_paths1=$(otool -L ${binary} | \
		grep "${LOCAL_DIR}" | \
		sed -e "s/(.*)//p" -n | sed -e "s/\t//p" -n)

	# locally compiled files beginning with "lib/"
	local old_paths2=$(otool -L ${binary} | \
		grep "lib/" | \
		sed -e "s/(.*)//p" -n | sed -e "s/\t//p" -n | \
		grep "^lib/")

	local old_paths=("$old_paths1[@] $old_paths2[@]")

	is_binary ${binary}
	if [[ $? == 0 ]]; then
		return
	fi

	for old_path in $old_paths; do
		case "${old_path}" in
			*".framework"*)
				local new_base=$(basename ${old_path})
				local cur_path=$(dirname ${old_path})
				local new_path=${new_base}

				while true; do
					local next_base=$(basename ${cur_path})
					local cur_path=$(dirname ${cur_path})
					local new_path="${next_base}/${new_path}"

					if [[ "${next_base}" == *".framework"* ]]; then
						break
					fi
				done

				local new_path=@rpath/${new_path}
				;;

			*)
				local new_path=@rpath/$(basename $old_path)
				;;
		esac

		echo -e "Changing linker path: $old_path -> $new_path"
		install_name_tool -change ${old_path} ${new_path} ${binary}
	done
}


#
# creates a png icon with the specified size out of an svg
#
function svg_to_png()
{
	local ICON_SVG="$1"
	local ICON_PNG="${APPICON%\.svg}.png"
	local ICON_SIZE="$2"

	echo -e "${ICON_SVG} -> ${ICON_PNG} (size: ${ICON_SIZE}x${ICON_SIZE})..."
	magick convert -resize "${ICON_SIZE}x${ICON_SIZE}" \
		-antialias -channel rgba \
		-background "#ffffff00" -alpha background \
		"${ICON_SVG}" "${ICON_PNG}"

	svg_to_png_result="${ICON_PNG}"
}


#
# creates the application icon
#
if [ $create_icon -ne 0 ] && [ -e ${APPICON} ]; then
	echo -e "\nCreating icons from ${APPICON}..."

	svg_to_png "${APPICON}" 512
	APPICON_PNG="${svg_to_png_result}"

	echo -e "${APPICON_PNG} -> ${APPICON_ICNS}..."
	makeicns -512 -in "${APPICON_PNG}" -out "${APPICON_ICNS}"

	echo -e "--------------------------------------------------------------------------------"
fi


#
# creates the application directory
#
if [ $create_appdir -ne 0 ]; then
	echo -e "\nCleaning and (re)creating directories..."
	rm -rfv "${APPDIRNAME}"

	mkdir -pv "${APPDIRNAME}/Contents/MacOS"
	mkdir -pv "${APPDIRNAME}/Contents/Resources"
	mkdir -pv "${APPDIRNAME}/Contents/Libraries"
	mkdir -pv "${APPDIRNAME}/Contents/Libraries/Qt_Plugins"
	mkdir -pv "${APPDIRNAME}/Contents/Frameworks"

	ln -sf "Libraries/Qt_Plugins" "${APPDIRNAME}/Contents/PlugIns"
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nCopying files to ${APPDIRNAME}..."

	# program files
	cp -v setup/osx/Info.plist  "${APPDIRNAME}/Contents/"
	cp -v build/magpie/magpie   "${APPDIRNAME}/Contents/MacOS/"

	# resources
	cp -v "${APPICON_ICNS}"     "${APPDIRNAME}/Contents/Resources/"
	cp -v "${APPICON}"          "${APPDIRNAME}/Contents/Resources/"
	cp -v AUTHORS               "${APPDIRNAME}/Contents/Resources/AUTHORS.txt"
	cp -v LICENSE               "${APPDIRNAME}/Contents/Resources/LICENSE.txt"
	cp -v LICENSES              "${APPDIRNAME}/Contents/Resources/LICENSES.txt"
	cp -v LITERATURE            "${APPDIRNAME}/Contents/Resources/LITERATURE.txt"

	# local libraries
	for (( libidx=0; libidx<${#LOCAL_LIBS[@]}; ++libidx )); do
		LOCAL_LIB=${LOCAL_LIBS[$libidx]}

		cp -v ${LOCAL_DIR}/lib/${LOCAL_LIB} "${APPDIRNAME}/Contents/Libraries/"
	done

	# libraries with full path
	for (( libidx=0; libidx<${#MISC_LIBS[@]}; ++libidx )); do
		MISC_LIB=${MISC_LIBS[$libidx]}

		cp -v ${MISC_LIB} "${APPDIRNAME}/Contents/Libraries/"
	done

	# frameworks
	for (( libidx=0; libidx<${#QT_LIBS[@]}; ++libidx )); do
		QT_LIB=${QT_LIBS[$libidx]}

		cp -rv ${LOCAL_FRAMEWORKS_DIR}/${QT_LIB}.framework "${APPDIRNAME}/Contents/Frameworks/"
	done

	# remove unnecessary files from frameworks
	if [ $clean_frameworks -ne 0 ]; then
		echo -e "\nCleaning frameworks..."
		find ${APPDIRNAME}/Contents/Frameworks/ -type d -name "Headers" -exec rm -rfv {} \;
		find ${APPDIRNAME}/Contents/Frameworks/ -type l -name "Headers" -exec rm -rv {} \;
		find ${APPDIRNAME} -type d -name "_CodeSignature" -exec rm -rfv {} \;
		echo -e "--------------------------------------------------------------------------------"
	fi

	# qt plugins
	cp -rv ${QT_PLUGIN_DIR}/platforms     "${APPDIRNAME}/Contents/Libraries/Qt_Plugins/"
	cp -rv ${QT_PLUGIN_DIR}/styles        "${APPDIRNAME}/Contents/Libraries/Qt_Plugins/"
	cp -rv ${QT_PLUGIN_DIR}/imageformats  "${APPDIRNAME}/Contents/Libraries/Qt_Plugins/"
	cp -rv ${QT_PLUGIN_DIR}/iconengines   "${APPDIRNAME}/Contents/Libraries/Qt_Plugins/"

	rm -fv ${APPDIRNAME}/Contents/Libraries/Qt_Plugins/platforms/libqoffscreen.dylib
	rm -fv ${APPDIRNAME}/Contents/Libraries/Qt_Plugins/platforms/libqwebgl.dylib

	rm -fv ${APPDIRNAME}/Contents/Libraries/Qt_Plugins/imageformats/libq[^sj][^vp][^g]*.dylib
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nChanging linked names..."
	rm -fv tmp/local_bindings.txt

	# binaries
	for binary in $(ls "${APPDIRNAME}/Contents/MacOS/"); do
		echo -e "\nProcessing binary ${binary}..."

		install_name_tool \
			-add_rpath @executable_path/../Libraries \
			-add_rpath @executable_path/../Frameworks \
			"${APPDIRNAME}/Contents/MacOS/${binary}"

		change_to_rpath      "${APPDIRNAME}/Contents/MacOS/${binary}"
		check_local_bindings "${APPDIRNAME}/Contents/MacOS/${binary}"
		chmod a+xr-w         "${APPDIRNAME}/Contents/MacOS/${binary}"

		if [ $strip_binaries -ne 0 ]; then
			${STRIP}           "${APPDIRNAME}/Contents/MacOS/${binary}"
		fi
	done

	# libraries and frameworks
	for library in $(find "${APPDIRNAME}/Contents/Libraries/" -type f && \
		find "${APPDIRNAME}/Contents/Frameworks/" -type f)
	do
		is_binary ${library}
		if [[ $? == 0 ]]; then
			continue
		fi

		echo -e "\nProcessing library ${library}..."

		change_to_rpath      "${library}"
		check_local_bindings "${library}"
		chmod a-wx+r         "${library}"

		if [ $strip_binaries -ne 0 ]; then
			${STRIP}           "${library}"
		fi
	done

	echo -e "--------------------------------------------------------------------------------"
fi


#
# creates a dmg image
#
if [ $create_dmg -ne 0 ]; then
	echo -e "\nCreating ${APPDMGNAME} from ${APPDIRNAME}..."
	rm -fv "${APPDMGNAME}"
	rm -fv "${TMPFILE}"
	if ! hdiutil create "${APPDMGNAME}" -srcfolder "${APPDIRNAME}" \
		-fs UDF -format "UDRW" -volname "${APPNAME}"
	then
		echo -e "${COL_ERR}Error: Cannot create ${APPDMGNAME}.${COL_NORM}"
		exit -1
	fi
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nMounting ${APPDMGNAME}..."
	if ! hdiutil attach "${APPDMGNAME}" -readwrite; then
		echo -e "${COL_ERR}Error: Cannot mount ${APPDMGNAME}.${COL_NORM}"
		exit -1
	fi

	echo -e "\nAdding files to ${APPDMGNAME}..."
	ln -sf /Applications "/Volumes/${APPNAME}/Install by dragging ${APPDIRNAME} here."

	echo -e "\nUnmounting ${APPDMGNAME}..."
	if ! hdiutil detach "/Volumes/${APPNAME}"; then
		echo -e "${COL_ERR}Error: Cannot detach ${APPDMGNAME}.${COL_NORM}"
		exit -1
	fi
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nCompressing ${APPDMGNAME} into ${TMPFILE}..."
	if ! hdiutil convert "${APPDMGNAME}" -o "${TMPFILE}" -format "UDBZ"
	then
		echo -e "${COL_ERR}Error: Cannot compress ${APPDMGNAME}.${COL_NORM}"
		exit -1
	fi
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nCopying ${APPDMGNAME}..."
	mv -v "${TMPFILE}" "${APPDMGNAME}"

	echo -e "\nSuccessfully created "${APPDMGNAME}"."
	echo -e "--------------------------------------------------------------------------------"
fi
