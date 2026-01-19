#
# magpie py bindings
# @author Tobias Weber <tweber@ill.fr>
# @date 2023 - 2026
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
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
project(magpie_c)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../setup/cmake")


message("Project: ${PROJECT_NAME}, build type: ${CMAKE_BUILD_TYPE}.")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()



# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS filesystem iostreams OPTIONAL_COMPONENTS system)

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Lapack(e)
find_package(GFortran)
find_package(Lapacke REQUIRED)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Gemmi
find_package(Gemmi REQUIRED)
# -----------------------------------------------------------------------------



include_directories(
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/libs"
	"${PROJECT_SOURCE_DIR}/.." "${PROJECT_SOURCE_DIR}/../tlibs2"
	"${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${Lapacke_INCLUDE_DIRS}"
	gemmi::headers
)



# -----------------------------------------------------------------------------
# magpie c library
add_library(magpie_c SHARED magpie_c.cpp magpie_c.h)

set_target_properties(magpie_c PROPERTIES CMAKE_CXX_STANDARD 20)
target_compile_options(magpie_c PRIVATE -std=c++20)
target_compile_options(magpie_c PUBLIC -Wall -Wextra)

target_compile_definitions(magpie_c PRIVATE -DUSE_LAPACK=1)
target_compile_definitions(magpie_c PRIVATE ${Boost_CXX_FLAGS})

target_link_libraries(magpie_c
	${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES}
	#/usr/local/lib/libgemmi_cpp.a
	gemmi::gemmi_cpp
)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# test programs
add_executable(test_calc_disp tests/calc_disp.c)

set_target_properties(test_calc_disp PROPERTIES CMAKE_C_STANDARD 99)
target_compile_options(test_calc_disp PRIVATE -std=c99)

target_link_libraries(test_calc_disp magpie_c)


add_executable(test_save_disp tests/save_disp.c)

set_target_properties(test_save_disp PROPERTIES CMAKE_C_STANDARD 99)
target_compile_options(test_save_disp PRIVATE -std=c99)

target_link_libraries(test_save_disp magpie_c)
# -----------------------------------------------------------------------------
