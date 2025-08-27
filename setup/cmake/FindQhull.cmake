#
# finds the qhull libs
# @author Tobias Weber <tweber@ill.fr>
# @date 13-aug-2020
# @license GPLv3
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
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

find_path(Qhull_INCLUDE_DIRS
	NAMES Qhull.h
	PATH_SUFFIXES libqhullcpp libqhull_r qhull
	HINTS ${PROJECT_SOURCE_DIR}/externals/qhull-inst/usr/local/include
		/usr/local/include/ /usr/include/ /opt/local/include
	DOC "Qhull include directories"
)


find_library(Qhull_c_LIBRARY
	NAMES libqhullstatic_r.a libqhull_r.a qhull_r
	HINTS ${PROJECT_SOURCE_DIR}/externals/qhull-inst/usr/local/lib
		/usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib
	DOC "Qhull C library"
)


find_library(Qhull_cpp_LIBRARY
	NAMES libqhullcpp.a qhullcpp
	HINTS ${PROJECT_SOURCE_DIR}/externals/qhull-inst/usr/local/lib
		/usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib
	DOC "Qhull C++ library"
)


if(Qhull_INCLUDE_DIRS AND Qhull_c_LIBRARY AND Qhull_cpp_LIBRARY)
	set(Qhull_FOUND TRUE)

	list(APPEND Qhull_LIBRARIES "${Qhull_c_LIBRARY}")
	list(APPEND Qhull_LIBRARIES "${Qhull_cpp_LIBRARY}")

	message("Qhull include directories: ${Qhull_INCLUDE_DIRS}")
	message("Qhull libraries: ${Qhull_LIBRARIES}")
else()
	set(Qhull_FOUND FALSE)

	message("Error: Qhull could not be found!")
endif()
