#
# finds the minuit libs
# @author Tobias Weber <tweber@ill.fr>
# @date 2016
# @license GPLv3
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

find_path(Minuit2_INCLUDE_DIRS
	NAMES MnMigrad.h
	PATH_SUFFIXES root Minuit2 Minuit root/Minuit2 root/Minuit Minuit2/Minuit2
	HINTS /usr/local/include/ /usr/include/ /opt/local/include
	DOC "Root/Minuit2 include directories"
)

# also include root base dir
list(APPEND Minuit2_INCLUDE_DIRS "${Minuit2_INCLUDE_DIRS}/..")


find_library(Minuit2_LIBRARIES
	NAMES Minuit2
	HINTS /usr/local/lib64 /usr/local/lib64/root /usr/local/lib /usr/local/lib/root /usr/lib64 /usr/lib64/root /usr/lib /usr/lib/root /opt/local/lib /opt/local/lib/root /usr/lib32 /usr/lib32/root
	DOC "Root/Minuit2 library"
)

if(Minuit2_INCLUDE_DIRS AND Minuit2_LIBRARIES)
	set(Minuit2_FOUND TRUE)

	message("Minuit include directories: ${Minuit2_INCLUDE_DIRS}")
	message("Minuit library: ${Minuit2_LIBRARIES}")
else()
	set(Minuit2_FOUND FALSE)

	message("Error: Minuit2 could not be found!")
endif()
