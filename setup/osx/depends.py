#!/usr/bin/env python3
#
# gets all dependent libraries of a binary
# @author Tobias Weber <tweber@ill.fr>
# @date 7-mar-2026
# @license GPLv3, see 'LICENSE' file
#
# -----------------------------------------------------------------------------
# Magpie
# Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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
# -----------------------------------------------------------------------------
#

import re
import os
import sys
import subprocess


#
# get the registered @rpaths of a binary
#
def get_rpaths(bin):
	res = subprocess.run([ "otool", "-l", bin ], capture_output = True, text = True)
	if res.returncode != 0:
		return []

	bin_dir = os.path.dirname(os.path.abspath(bin))

	rpaths = set()
	has_rpath = False
	skip_lines = 0
	for line in res.stdout.splitlines():
		if skip_lines > 0:
			skip_lines -= 1
			continue

		if has_rpath and skip_lines == 0:
			# get @rpath
			path = re.split("\\(offset", line)[0].strip()[5:]
			# resolve @loader_path
			path = path.replace("@loader_path", bin_dir)
			rpaths.add(path)
			has_rpath = False
			continue

		if line.find("LC_RPATH") >= 0:
			has_rpath = True
			skip_lines = 1

	return rpaths


# set of found library dependencies
libs = set()


#
# add a new library to the set
#
def add_lib(lib):
	lib = os.path.abspath(lib)
	if lib in libs:
		return False  # already seen
	libs.add(lib)

	# find dependencies for this library
	get_depends(lib)
	return True


#
# get the library dependencies of a binary
#
def get_depends(bin):
	print("Processing \"%s\"..." % bin)
	rpaths = get_rpaths(bin)
	if len(rpaths) > 0:
		print("@rpaths: %s." % rpaths)

	res = subprocess.run([ "otool", "-L", bin ], capture_output = True, text = True)
	if res.returncode != 0:
		print("Could not run otool for %s." % bin, file = sys.stderr)
		#print(res.stderr, file = sys.stderr)
		return

	for line in res.stdout.splitlines()[1:]:
		lib = re.split("\\(compatibility", line)[0].strip()

		# is there an unresolved @rpath?
		if lib.find("@rpath") >= 0:
			for rpath in rpaths:
				new_lib = lib.replace("@rpath", rpath)

				if not add_lib(new_lib):
					continue
		else:
			if not add_lib(lib):
				continue


#
# main
#
if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Please give a program binary.", file = sys.stderr)
		exit(-1)

	for arg in sys.argv[1:]:
		get_depends(arg)

	# restrict to existing files
	libs = { lib for lib in libs if os.path.isfile(lib) }
	libs = sorted(libs)

	print("\nFound %d dependencies:" % len(libs))
	for lib in libs:
		print("\t%s" % lib)
