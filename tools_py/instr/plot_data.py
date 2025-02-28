#
# plots data files
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

import numpy
import instrhelper


#
# load scan files
#
def main(argv):
	if len(argv) < 2:
		print("Please give a scan file.")
		exit(-1)

	scanfile = argv[1]
	mergefiles = argv[2:]
	[ hs, ks, ls, Es, Is, Is_err, instrdat ] = instrhelper.load_data(scanfile, mergefiles)
	h = numpy.mean(hs)
	k = numpy.mean(ks)
	l = numpy.mean(ls)

	# save a simplified version of the data file
	#instrhelper.save_data("tst.dat", [ hs, ks, ls, Es, Is, Is_err, instrdat ])


	# plot the scan data
	import matplotlib.pyplot as plt

	plt.figure()
	plt.title("Q = (%.2f %.2f %.2f)" % (h, k, l))
	plt.xlabel("E (meV)")
	plt.ylabel("S (a.u.)")
	plt.errorbar(Es, Is, Is_err, marker="o", markersize=6, capsize=4, ls="none")
	plt.tight_layout()
	plt.show()
	plt.close()


if __name__ == "__main__":
	import sys
	main(sys.argv)
