#
# magdyn -- unit tests
# @author Tobias Weber <tweber@ill.fr>
# @date 8-jun-2020
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core & Magpie
# Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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
project(tlibs2_tests)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
enable_testing()

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../setup/cmake")


message("Build type: ${CMAKE_BUILD_TYPE}")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()


set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra)

add_definitions(-D__MAGDYN_USE_MINUIT__)
add_definitions(-D__TLIBS2_USE_MINUIT__)


# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_STACKTRACE_USE_NOOP)
find_package(Boost REQUIRED)
add_definitions(${Boost_CXX_FLAGS})
include_directories(SYSTEM "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/..")

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Lapacke
find_package(GFortran)
find_package(Lapacke REQUIRED)
add_definitions(-DUSE_LAPACK)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# further packages
find_package(Minuit2 REQUIRED)
find_package(Threads REQUIRED)
#find_package(OpenMP REQUIRED)
# -----------------------------------------------------------------------------


include_directories(
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/.."
	"${Minuit2_INCLUDE_DIRS}"
	"${Lapacke_INCLUDE_DIRS}"
)


# -----------------------------------------------------------------------------
# build individual test cases
add_executable(magdyn magdyn.cpp)
add_executable(magdyn2 magdyn2.cpp)

target_link_libraries(magdyn ${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES}
	${GFortran_LIBRARIES} ${Minuit2_LIBRARIES} Threads::Threads)
target_link_libraries(magdyn2 ${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES}
	${GFortran_LIBRARIES} ${Minuit2_LIBRARIES} Threads::Threads)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# run individual test cases
add_test(magdyn magdyn)
add_test(magdyn2 magdyn2)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# add rpaths
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	foreach(tool IN ITEMS magdyn magdyn2)
		message(STATUS "Setting rpath for ${tool}...")

		add_custom_command(TARGET ${tool} POST_BUILD
			COMMAND install_name_tool -add_rpath /usr/local/lib $<TARGET_FILE:${tool}>
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	endforeach()
endif()
# -----------------------------------------------------------------------------
