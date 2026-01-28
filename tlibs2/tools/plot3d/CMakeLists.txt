#
# 3d plotter
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2026
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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
project(plot3d)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra -Weffc++)


# -----------------------------------------------------------------------------
#options
# -----------------------------------------------------------------------------
option(USE_QT6 "use Qt version 6"  TRUE)

if(USE_QT6)
	set(QT_VER 6)
else()
	set(QT_VER 5)
endif()
# -----------------------------------------------------------------------------


list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../../../setup/cmake"
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_compile_options(-g -ggdb)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)

	message("Debug building for ${CMAKE_SYSTEM_NAME} systems.")

elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
	add_compile_options(-g -ggdb3 -pg -O0)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)

	message("Profile building for ${CMAKE_SYSTEM_NAME} systems.")

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_compile_options("-DNDEBUG")
	add_compile_options("-Wno-#pragma-messages")

	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		# see: https://gcc.gnu.org/onlinedocs/gcc/Developer-Options.html
		add_compile_options(-fcompare-debug-second)
	endif()

	message("Release building for ${CMAKE_SYSTEM_NAME} systems.")
endif()


set(MINGW_EXTRALIBS)
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	# pass linker --subsystem option
	add_compile_options(-Wl,--subsystem,windows)
	set(MINGW_EXTRALIBS "ws2_32")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_options(-mmacosx-version-min=10.15)
endif()


# -----------------------------------------------------------------------------
# packages
# -----------------------------------------------------------------------------
find_package(Threads REQUIRED)

find_package(Boost REQUIRED COMPONENTS program_options)
set(Boost_USE_MULTITHREADED TRUE)
add_definitions(${Boost_CXX_FLAGS})

find_package(Qhull REQUIRED)
if(Qhull_FOUND)
	add_definitions(-DUSE_QHULL)
	include_directories("${Qhull_INCLUDE_DIRS}")
endif()
# -----------------------------------------------------------------------------


include_directories("${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../.." "${PROJECT_SOURCE_DIR}/../../..")

include_directories(SYSTEM "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/..")


message("Trying to select Qt version ${QT_VER}.")

if(QT_VER EQUAL 6)
	find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGLWidgets PrintSupport Svg)
	include_directories("${Qt6Core_INCLUDE_DIRS}/..")
	list(APPEND QtLibraries Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGLWidgets Qt6::PrintSupport Qt6::Svg)
elseif(QT_VER EQUAL 5)
	find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets PrintSupport Svg)
	include_directories("${Qt5Core_INCLUDE_DIRS}/..")
	list(APPEND QtLibraries Qt5::Core Qt5::Gui Qt5::Widgets Qt5::PrintSupport Qt5::Svg)
else()
	message(FATAL_ERROR "Tried to select a non-supported Qt version ${QT_VER}.")
endif()

set(CMAKE_AUTOUIC TRUE)
set(CMAKE_AUTOMOC TRUE)

add_executable(plot3d
	main.cpp
	plot3d.cpp plot3d_file.cpp
	plot3d_calc.cpp plot3d.h
	defs.cpp defs.h

	# library
	../../libs/qt/gl.cpp ../../libs/qt/gl.h
	../../libs/qt/glplot.cpp ../../libs/qt/glplot.h
	../../libs/qt/recent.cpp ../../libs/qt/recent.h
)

target_link_libraries(plot3d
	Threads::Threads
	${QtLibraries}
	${Qhull_LIBRARIES}
	${Boost_LIBRARIES}
	${MINGW_EXTRALIBS}
)


# -----------------------------------------------------------------------------
# rpath
# -----------------------------------------------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_link_options(plot3d
		# create an __info_plist section in the binary
		PRIVATE LINKER:-sectcreate,__TEXT,__info_plist,${PROJECT_SOURCE_DIR}/../../../setup/plists/InfoPlot3D.plist
	)

	add_custom_command(TARGET plot3d POST_BUILD
		COMMAND install_name_tool -add_rpath /usr/local/lib $<TARGET_FILE:plot3d> || true
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)
endif()
# -----------------------------------------------------------------------------
