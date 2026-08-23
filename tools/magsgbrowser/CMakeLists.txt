#
# @author Tobias Weber
# @date Apr-2018
# @license GPLv3, see 'LICENSE' file
# @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
#
# ----------------------------------------------------------------------------
# mag-core / magpie
# Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# "magtools" project
# Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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
project(magsgbrowser)

set(CMAKE_CXX_STANDARD 20)
add_definitions(-std=c++20)


find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)


find_package(Boost REQUIRED)
set(Boost_USE_MULTITHREADED TRUE)
add_definitions(${Boost_CXX_FLAGS})


find_package(Gemmi REQUIRED)


include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/../.."
	"${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	gemmi::headers)


# -----------------------------------------------------------------------------
# main program
add_executable(magsgbrowser
	main.cpp
	browser.cpp browser.h
)

target_link_libraries(magsgbrowser
	${Boost_LIBRARIES}
	gemmi::gemmi_cpp
	Qt6::Core Qt6::Gui Qt6::Widgets
)
# -----------------------------------------------------------------------------

