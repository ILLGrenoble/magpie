#
# magpy interface demo -- creating and plotting a model
# @author Tobias Weber <tweber@ill.fr>
# @date june-2025
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magpie
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
# ----------------------------------------------------------------------------
#

import numpy
import magpy


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
save_config_file       = False  # save magpie file
save_dispersion        = False  # write dispersion to file
print_dispersion       = False  # write dispersion to console
plot_dispersion        = True   # show dispersion plot
only_positive_energies = True   # ignore magnon annihilation?

use_pointwise          = False  # use point-wise dispersion calculation
use_threadpool         = True   # parallelise calculation
max_threads            = 4      # number of worker threads, 0: automatic determination
num_Q_points           = 1024   # number of Qs to calculate on a dispersion direction

S_scale                = 4.     # weight scaling and clamp factors
S_clamp_min            = 0.     #
S_clamp_max            = 9990.  #
S_filter_min           = -1.    # don't filter

dispfile               = "disp.dat"
modelfile              = "MnTe.magpie"

# dispersion plotting range
hkl_start              = numpy.array([ 1./3., 0., 0. ])
hkl_end                = numpy.array([ 1./3., 0., 2. ])
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Create the magnetic model
# -----------------------------------------------------------------------------
# The given magnetic model and its parameters are from this paper:
#	https://doi.org/10.1103/PhysRevLett.133.156702
#	(which is also available here: https://arxiv.org/abs/2408.16490).
#
# Thanks to N. Heinsdorf for providing his spinw model as a check.
#

# create the magpy object
mag = magpy.MagDyn()

spacegroup = "P 63/m m c"
mag.SetCrystalLattice(4.1, 4.1, 6.62,
	90./180.*numpy.pi, 90./180.*numpy.pi, 120./180.*numpy.pi)


#
# add coupling variables (values from https://arxiv.org/abs/2408.16490)
#
magpy.add_variable(mag, "J_afm1", 3.99)
magpy.add_variable(mag, "J_afm2", 0.47)

magpy.add_variable(mag, "J_fm1", -0.12)

magpy.add_variable(mag, "J_am1",  0.068)
magpy.add_variable(mag, "J_am2", -0.022)


#
# set up the sample environment
#
magpy.set_temperature(mag, 50.)
#magpy.set_field(mag,  1., -2., 0.,  2.5,  False)


#
# create magnetic sites
#
print("Creating sites...")

#                   name    x    y    z     sx   sy   sz    S
magpy.add_site(mag, "Mn",   0.,  0.,  0.,   0., +1.,  0.,   5./2.)
magpy.symmetrise_sites(mag, spacegroup)
magpy.set_spin(mag, "1", 0., -1., 0.)
magpy.print_sites(mag)
print()


#
# create couplings between sites
#
print("Creating couplings...")

#                              d  sc  cnt
magpy.generate_couplings(mag, 10., 5, 200, spacegroup)
#magpy.symmetrise_couplings(mag, spacegroup)
#mag.SortExchangeTerms()

magpy.assign_couplings_by_symmetry(mag,  1, "J_afm1")
magpy.assign_couplings_by_symmetry(mag,  2, "J_fm1")
magpy.assign_couplings_by_symmetry(mag,  3, "J_afm2")
magpy.assign_couplings_by_symmetry(mag, 10, "J_am1")
magpy.assign_couplings_by_symmetry(mag, 11, "J_am2")

mag.RemoveUnusedExchangeTerms()
magpy.print_couplings(mag)
print()
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate sites and couplings
# -----------------------------------------------------------------------------
print("Calculating structure...")
magpy.calc(mag)


# minimum energy
print("\nEnergy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))


if save_config_file:
	# save the magpy configuration file for use in the gui program or the convolution fitter
	print("\nSaving configuration to {}...".format(modelfile))
	if not mag.Save(modelfile):
		print("Error saving configuration to {}".format(modelfile))


if save_dispersion:
	# directly calculate a dispersion and write it to a file
	print("\nSaving dispersion to {}...".format(dispfile))
	mag.SaveDispersion(dispfile,
		hkl_start[0], hkl_start[1], hkl_start[2],
		hkl_end[0], hkl_end[1], hkl_end[2],
		num_Q_points)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# manually calculate the dispersion
# -----------------------------------------------------------------------------
print("\nCalculating dispersion...")
if print_dispersion:
	print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))


#
# add a data point
#
data_Q = [ [], [], [] ]
data_E, data_S = [], []

def append_data(h, k, l, E, S):
	weight = S * S_scale

	if weight < S_filter_min:
		return

	if weight < S_clamp_min:
		weight = S_clamp_min
	elif weight > S_clamp_max:
		weight = S_clamp_max

	data_Q[0].append(h)
	data_Q[1].append(k)
	data_Q[2].append(l)
	data_E.append(E)
	data_S.append(weight)


# classical way to calculate the dispersion along straight lines
if not use_pointwise:
	data_disp = mag.CalcDispersion(
		hkl_start[0], hkl_start[1], hkl_start[2],
		hkl_end[0], hkl_end[1], hkl_end[2],
		num_Q_points, max_threads)
	for S in data_disp:
		for data_EandS in S.E_and_S:
			if only_positive_energies and data_EandS.E < 0.:
				continue

			append_data(magpy.get_h(S), magpy.get_k(S), magpy.get_l(S),
				data_EandS.E, data_EandS.weight)

			if print_dispersion:
				print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
					magpy.get_h(S), magpy.get_k(S), magpy.get_l(S),
					data_EandS.E, data_EandS.weight))

# alternate way to calculate the dispersion branches in a point-wise fashion
else:
	if use_threadpool:
		import concurrent.futures as fut
		executor = fut.ThreadPoolExecutor
		#executor = fut.ProcessPoolExecutor
		print(f"Using {max_threads} threads.")

		#
		# calculate the energies and weights for a Q point
		#
		def calc_Es(h, k, l):
			Es_dicts = mag.CalcEnergies(h, k, l, False)

			Es = []
			for Es_dict in Es_dicts:
				Es.append(( h, k, l, Es_dict.E, Es_dict.weight ))

			return Es


		with executor(max_workers = max_threads) as exe:
			Es_futures = []

			# submit tasks
			for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
				Es_futures.append(
					exe.submit(calc_Es, hkl[0], hkl[1], hkl[2]))

			# get results from tasks
			for Es_future in Es_futures:
				if Es_future == None:
					continue
				result = Es_future.result()
				if result == None:
					continue

				for (h, k, l, E, weight) in result:
					if only_positive_energies and E < 0.:
						continue

					append_data(h, k, l, E, weight)

					if print_dispersion:
						print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
							h, k, l, E, weight))

	else:  # no thread pool
		for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
			for S in mag.CalcEnergies(hkl[0], hkl[1], hkl[2], False):
				if only_positive_energies and S.E < 0.:
					continue

				append_data(hkl[0], hkl[1], hkl[2], S.E, S.weight)

				if print_dispersion:
					print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
						hkl[0], hkl[1], hkl[2], S.E, S.weight))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the results
# -----------------------------------------------------------------------------
if plot_dispersion:
	# find dominant q index
	q_idx = 0
	q_label = "h (rlu)"
	q_dst = numpy.abs(hkl_end[q_idx] - hkl_start[q_idx])
	if numpy.abs(hkl_end[1] - hkl_start[1]) > q_dst:
		q_idx = 1
		q_dst = numpy.abs(hkl_end[1] - hkl_start[1])
		q_label = "k (rlu)"
	if numpy.abs(hkl_end[2] - hkl_start[2]) > q_dst:
		q_idx = 2
		q_label = "l (rlu)"

	import matplotlib.pyplot as plot
	print("Plotting dispersion...")

	fig = plot.figure()

	plt = fig.add_subplot(1, 1, 1)
	plt.set_xlabel(q_label)
	plt.set_ylabel("E (meV)")
	plt.scatter(data_Q[q_idx], data_E, marker = '.', s = data_S)

	plot.tight_layout()
	plot.show()
# -----------------------------------------------------------------------------
