#
# convert all data files in a directory to a simplified format
# @author Tobias Weber <tweber@ill.fr>
# @date 1-apr-2024
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

import os
import numpy
import instrhelper


#
# convert all data files in a given directory
#
def convert_data_files(in_dir, out_dir):
	for datfile in os.listdir(in_dir):
		instrhelper.save_data(out_dir + "/" + datfile,
			instrhelper.load_data(in_dir + "/" + datfile))
		print()


def main(argv):
	if len(argv) < 3:
		print("Please give an input and an output directory.")
		exit(-1)

	in_dir = argv[1]
	out_dir = argv[2]

	convert_data_files(in_dir, out_dir)


if __name__ == "__main__":
	import sys
	main(sys.argv)
