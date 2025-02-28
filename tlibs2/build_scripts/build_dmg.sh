#!/bin/bash
#
# creates an app bundle and a dmg file
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2019, apr-2021
# @license GPLv3, see 'LICENSE' file
#
# @note Forked on 16-july-2022 from my TAS-Paths project (https://code.ill.fr/scientific-software/takin/paths).
#
# -----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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

create_appdir=1
create_dmg=1
strip_binaries=1
clean_frameworks=1


# gets application name from the command line
APPNAME="takin_app"
if [[ $# > 0 ]]; then
	APPNAME="$1"
fi

APPDIRNAME="${APPNAME}.app"
APPDMGNAME="${APPNAME}.dmg"
TMPFILE="${APPNAME}_tmp.dmg"


# libraries
declare -a QT_LIBS=(
	QtCore
	QtGui
	QtWidgets
	QtPrintSupport
	QtDBus
	QtSvg
)

declare -a QT_PLUGIN_LIBS=(
	platforms
	styles
	#renderers
	#imageformats
	#iconengines
	#printsupport
	#platformthemes
)

declare -a OTHER_LIBS=(
	/usr/local/lib/libqhull_r.8.0.dylib
	#/usr/local/opt/lapack/lib/liblapacke.3.dylib
	#/usr/local/opt/lapack/lib/liblapack.3.dylib
)


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
# tests if the given binary still has remaining bindings to /usr/local/
#
function check_local_bindings()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "${COL_WARN}Warning: ${cfile} does not exist.${COL_NORM}"
		return
	fi

	local local_binding=$(otool -L ${cfile} | grep --color /usr/local)
	local num_local_bindings=$(echo -e "${local_binding}" | wc -l)

	if [ "${num_local_bindings}" -gt "1" ]; then
		echo -e "${COL_WARN}Warning: Possible local binding remaining, please check \"local_bindings.txt\"!${COL_NORM}"
	fi

	echo -e "--------------------------------------------------------------------------------" >> local_bindings.txt
	echo -e "${num_local_bindings} local binding(s) in binary \"${cfile}\":" >> local_bindings.txt
	echo -e "${local_binding}" >> local_bindings.txt
	echo -e "--------------------------------------------------------------------------------\n" >> local_bindings.txt
}



#
# changes a /usr/local linker path to an @rpath
#
function change_to_rpath()
{
	local binary="$1"
	local old_paths=$(otool -L ${binary} | grep /usr/local | sed -e "s/(.*)//p" -n | sed -e "s/\t//p" -n)

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
	cp -v Info.plist "${APPDIRNAME}/Contents/"
	cp -v ${APPNAME} "${APPDIRNAME}/Contents/MacOS/"

	# resources
	cp -rv 3rdparty_licenses "${APPDIRNAME}/Contents/Resources/"

	# libraries
	for library in ${OTHER_LIBS[*]}; do
		cp -v ${library} "${APPDIRNAME}/Contents/Libraries/"
	done

	# frameworks
	for (( libidx=0; libidx<${#QT_LIBS[@]}; ++libidx )); do
		QT_LIB=${QT_LIBS[$libidx]}

		cp -rv /usr/local/Frameworks/${QT_LIB}.framework/ \
			"${APPDIRNAME}/Contents/Frameworks/"
	done

	# remove unnecessary files from frameworks
	if [ $clean_frameworks -ne 0 ]; then
		echo -e "\nCleaning frameworks..."
		find ${APPDIRNAME}/Contents/Frameworks/ -type d -name "Headers" -exec rm -rfv {} \;
		find ${APPDIRNAME}/Contents/Frameworks/ -type l -name "Headers" -exec rm -rv {} \;
		echo -e "--------------------------------------------------------------------------------"
	fi

	# qt plugins
	for library in ${QT_PLUGIN_LIBS[*]}; do
		cp -rv /usr/local/opt/qt@5/plugins/${library} "${APPDIRNAME}/Contents/Libraries/Qt_Plugins/"
	done

	rm -fv ${APPDIRNAME}/Contents/Libraries/Qt_Plugins/platforms/libqoffscreen.dylib
	rm -fv ${APPDIRNAME}/Contents/Libraries/Qt_Plugins/platforms/libqwebgl.dylib
	echo -e "--------------------------------------------------------------------------------"


	echo -e "\nChanging linked names..."
	rm -fv local_bindings.txt

	# binaries
	for binary in $(ls "${APPDIRNAME}/Contents/MacOS/"); do
		echo -e "\nProcessing ${binary}..."

		install_name_tool \
			-add_rpath @executable_path/../Libraries \
			-add_rpath @executable_path/../Frameworks \
			"${APPDIRNAME}/Contents/MacOS/${binary}"

		change_to_rpath "${APPDIRNAME}/Contents/MacOS/${binary}"
		check_local_bindings "${APPDIRNAME}/Contents/MacOS/${binary}"

		if [ $strip_binaries -ne 0 ]; then
			/usr/local/Cellar/llvm/13.0.1_1/bin/llvm-strip "${APPDIRNAME}/Contents/MacOS/${binary}"
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

		echo -e "\nProcessing ${library}..."

		change_to_rpath "${library}"
		check_local_bindings "${library}"

		if [ $strip_binaries -ne 0 ]; then
			/usr/local/Cellar/llvm/13.0.1_1/bin/llvm-strip "${library}"
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
	ln -sf /Applications \
		"/Volumes/${APPNAME}/Drag ${APPDIRNAME} here to install."

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
