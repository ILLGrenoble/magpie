#
# brillouin zone scripting test
# @author Tobias Weber <tweber@ill.fr>
# @date April-2023
# @license GPLv3, see 'LICENSE' file
#

import os
import sys
import magpy

sys.path.append(os.getcwd())


# create brillouin zone calculator
bz = magpy.BZCalc()
bz.SetEps(1e-6)


# set-up crystal
sg_name = "F d -3 m"
bz.SetCrystal(5, 5, 5, 90, 90, 90)
num_ops = bz.SetSymOpsFromSpaceGroup(sg_name)
print("Using %d centring symops." % num_ops)


# calculate Bragg peaks
num_peaks = bz.CalcPeaks(4, True)
print("Using %d reflections." % num_peaks)


# calculate brillouin zone
calc_ok = bz.CalcBZ()
if calc_ok:
	print("Brillouin zone calculation successful.")
else:
	print("Brillouin zone calculation failed.")


# output results
if calc_ok:
	print("\nJSON Output:")
	json = bz.PrintJSON(6)
	print(json)
