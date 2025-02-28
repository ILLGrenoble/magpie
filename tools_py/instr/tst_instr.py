#
# tlibs2 py interface test
# @author Tobias Weber <tweber@ill.fr>
# @date 9-jun-2020
# @license GPLv3, see 'LICENSE' file
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

import math
import instr

import instrhelper


#
# load an instrument data file
#
def load_data(datfile):
	#print("Loading \"%s\"." % (datfile))
	dat = instr.FileInstrBaseD.LoadInstr(datfile)
	if dat == None:
		return

	cnt = dat.GetCountVar()
	mon = dat.GetMonVar()
	cntcol = dat.GetCol(cnt)
	moncol = dat.GetCol(mon)

	for point_idx in range(cntcol.size()):
		(h, k, l, ki, kf) = dat.GetScanHKLKiKf(point_idx)
		E = instrhelper.get_E(ki, kf)

		counts = cntcol[point_idx]
		counts_err = math.sqrt(counts)
		mon_counts = moncol[point_idx]
		intensity = counts/mon_counts
		intensity_err = counts_err / mon_counts

		print("{0:12.4g} {1:12.4g} {2:12.4g} {3:12.4g} {4:12.4g} {5:12.4g}".format(h, k, l, E, intensity, intensity_err))

		#print("Q = (%.4f %.4f %.4f), E = %.4f: Monitor: %d, Counts: %d +- %d, Counts/Monitor: %.5g +- %.5g" \
		#	% (h, k, l, E, mon_counts, counts, counts_err, intensity, intensity_err))
	#print()


#
# load all instrument data files in a directory
#
def load_all(dir):
	import os

	for datfile in os.listdir(dir):
		load_data(dir + "/" + datfile)


def main(argv):
	if len(argv) < 2:
		print("No scan directory given.")
		exit(-1)

	print("#          h            k            l            E            S        S_err")
	load_all(argv[1])


if __name__ == "__main__":
	import sys
	main(sys.argv)
