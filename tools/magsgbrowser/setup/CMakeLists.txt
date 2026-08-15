#
# @author Tobias Weber
# @date Apr-2018
# @license GPLv3, see 'LICENSE' file
# @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

project(magsg_setup)
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
add_definitions(-std=c++20)
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/../../..")


message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}.")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}.")
message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}.")


find_package(Boost REQUIRED)
set(Boost_USE_MULTITHREADED TRUE)
add_definitions(${Boost_CXX_FLAGS})
include_directories("${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/..")


# -----------------------------------------------------------------------------
# setup tool
add_executable(magsg_setup convmag.cpp)
target_link_libraries(magsg_setup ${Boost_LIBRARIES})
# -----------------------------------------------------------------------------
