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


# dependent libraries
libs = set()


def get_depends(bin):
	print("Processing \"%s\"..." % bin)

	res = subprocess.run([ "otool", "-L", bin ], capture_output = True, text = True)
	if res.returncode != 0:
		print("Could not run otool for %s." % bin, file = sys.stderr)
		#print(res.stderr, file = sys.stderr)
		return

	lines = res.stdout.splitlines()

	for line in lines[1:]:
		lib = re.split("\\(compatibility", line)[0].strip()
		if lib in libs:
			continue  # already seen

		libs.add(lib)
		get_depends(lib)


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
