#
# helper functions for instrument data file loader
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

import instr
import numpy
import scipy.constants as co


#
# normalise detector (y) by monitor (m) counts
# y_new = y / m
# dy_new = 1/m dy - y/m^2 dm
#
def norm_counts_to_mon(y, dy, m, dm):
	val = y / m;
	err = numpy.sqrt((dy/m)**2. + (dm*y/(m*m))**2.)

	return [val, err]


#
# get energy transfer from ki and kf
#
def get_E(ki, kf):
	E_to_k2 = 2.*co.neutron_mass/co.hbar**2. * co.elementary_charge/1000. * 1e-20
	return (ki**2. - kf**2.) / E_to_k2


#
# load an instrument data file
#
def load_data(datfile, mergefiles = [], intensity_scale = 1.):
	print("Loading \"%s\"." % (datfile))
	instrdat = instr.FileInstrBaseD.LoadInstr(datfile)
	if instrdat == None:
		return

	for mergefile in mergefiles:
		tomerge = instr.FileInstrBaseD.LoadInstr(mergefile)
		if tomerge != None:
			print("Merging with \"%s\"." % (mergefile))
			instrdat.MergeWith(tomerge)


	hs = []; ks = []; ls = []; Es = []
	Is = []; Is_err = []

	cntcol = instrdat.GetCol(instrdat.GetCountVar())
	moncol = instrdat.GetCol(instrdat.GetMonVar())

	# iterate scan points
	for point_idx in range(cntcol.size()):
		( h, k, l, ki, kf ) = instrdat.GetScanHKLKiKf(point_idx)
		E = get_E(ki, kf)

		hs.append(h)
		ks.append(k)
		ls.append(l)
		Es.append(E)

		counts = cntcol[point_idx]
		mon = moncol[point_idx]
		counts_err = numpy.sqrt(counts)
		mon_err = numpy.sqrt(mon)

		if counts == 0:
			counts_err = 1.
		if mon == 0:
			mon_err = 1.

		( I, I_err ) = norm_counts_to_mon(counts, counts_err, mon, mon_err)
		Is.append(I * intensity_scale)
		Is_err.append(I_err * intensity_scale)

	return (hs, ks, ls, Es, Is, Is_err, instrdat)


#
# save the normalised scans to simplified data files
#
def save_data(datfile, dat):
	if dat == None:
		return
	print("Saving normalised scan to \"%s\"." % (datfile))

	hs = dat[0]
	ks = dat[1]
	ls = dat[2]
	Es = dat[3]
	Is = dat[4]
	Is_err = dat[5]
	instrdat = dat[6]

	header = """
sample_a     = %f
sample_b     = %f
sample_c     = %f

sample_alpha = %f
sample_beta  = %f
sample_gamma = %f

mono_d       = %f
ana_d        = %f

sense_m      = %d
sense_s      = %d
sense_a      = %d

orient1_x    = %f
orient1_y    = %f
orient1_z    = %f

orient2_x    = %f
orient2_y    = %f
orient2_z    = %f

is_ki_fixed  = %d
k_fix        = %f

col_h        = 1
col_k        = 2
col_l        = 3
col_E        = 4
cols_scanned = 4
col_ctr      = 5
col_ctr_err  = 6
	""" % (
		instrdat.GetSampleLattice()[0],
		instrdat.GetSampleLattice()[1],
		instrdat.GetSampleLattice()[2],
		instrdat.GetSampleAngles()[0] / numpy.pi * 180.,
		instrdat.GetSampleAngles()[1] / numpy.pi * 180.,
		instrdat.GetSampleAngles()[2] / numpy.pi * 180.,
		instrdat.GetMonoAnaD()[0], instrdat.GetMonoAnaD()[1],
		1 if instrdat.GetScatterSenses()[0] else -1,
		1 if instrdat.GetScatterSenses()[1] else -1,
		1 if instrdat.GetScatterSenses()[2] else -1,
		instrdat.GetScatterPlane0()[0],
		instrdat.GetScatterPlane0()[1],
		instrdat.GetScatterPlane0()[2],
		instrdat.GetScatterPlane1()[0],
		instrdat.GetScatterPlane1()[1],
		instrdat.GetScatterPlane1()[2],
		1 if instrdat.IsKiFixed() else 0,
		instrdat.GetKFix())

	numpy.savetxt(datfile, numpy.transpose((hs, ks, ls, Es, Is, Is_err)), header = header)
