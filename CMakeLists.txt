#
# build all magpie tools
# @author Tobias Weber <tweber@ill.fr>
# @date 17-jun-2019
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magpie
# Copyright (C) 2022-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# takin / magcore
# Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

cmake_minimum_required(VERSION 3.10)
project(magpie)


option(BUILD_MAGPIE        "build main program"         TRUE)
option(BUILD_MODULES       "build tool programs"        TRUE)
option(BUILD_PY_MODULES    "build py scripting modules" TRUE)
option(ONLY_BUILD_FINISHED "only build finished tools"  TRUE)
option(USE_QT6             "use Qt version 6"           TRUE)


# mingw settings
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	set(Boost_USE_STATIC_LIBS FALSE)
endif()


# main program
if(BUILD_MAGPIE)
	# magpie -- magnetic dynamics tool
	message("\n================================================================================")
	message("Configuring magpie...")
	message("================================================================================")
	add_subdirectory(tools/magdyn magpie)
endif()


# additional modules (used by takin)
if(BUILD_MODULES)
	# Brillouin zone tool
	message("\n================================================================================")
	message("Configuring bz...")
	message("================================================================================")
	add_subdirectory(tools/bz)

	# cif converter
	message("\n================================================================================")
	message("Configuring cif2xml...")
	message("================================================================================")
	add_subdirectory(tools/cif2xml)

	# nuclear and magnetic structure factor tools
	message("\n================================================================================")
	message("Configuring structfact...")
	message("================================================================================")
	add_subdirectory(tools/structfact)


	if(NOT ONLY_BUILD_FINISHED)
		# scan browser tool
		message("\n================================================================================")
		message("Configuring scanbrowser...")
		message("================================================================================")
		add_subdirectory(tools/scanbrowser)

		# setup tools for magnetic space groups
		message("\n================================================================================")
		message("Configuring setup...")
		message("================================================================================")
		add_subdirectory(tools/magsgbrowser/setup)

		# magnetic space group browser
		message("\n================================================================================")
		message("Configuring magsgbrowser...")
		message("================================================================================")
		add_subdirectory(tools/magsgbrowser)

		# molecular dynamics tool
		message("\n================================================================================")
		message("Configuring moldyn...")
		message("================================================================================")
		add_subdirectory(tools/moldyn)

		# tlibs2 script compiler
		message("\n================================================================================")
		message("Configuring mcalc...")
		message("================================================================================")
		add_subdirectory(tlibs2/script)
	endif()
endif()


# scripting modules
if(BUILD_PY_MODULES)
	# magnetic dynamics
	message("\n================================================================================")
	message("Configuring py-magdyn...")
	message("================================================================================")
	add_subdirectory(tools_py/magdyn)

	# brillouin zones
	message("\n================================================================================")
	message("Configuring py-bz...")
	message("================================================================================")
	add_subdirectory(tools_py/bz)

	# instrument data file loader
	message("\n================================================================================")
	message("Configuring py-instr...")
	message("================================================================================")
	add_subdirectory(tools_py/instr)
endif()
