#
# magpy interface demo -- plotting multiple dispersion branches
# @author Tobias Weber <tweber@ill.fr>
# @date 18-sep-2024
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

import os
import time
import numpy
import numpy.linalg
import magpy


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
print_dispersion       = False  # write dispersion to console
only_positive_energies = True   # ignore magnon annihilation?

use_pointwise          = False  # use point-wise dispersion calculation
use_threadpool         = True   # parallelise calculation
max_threads            = 0      # number of worker threads or processes, 0: automatic determination

show_dividers          = False  # show vertical bars between dispersion branches
show_markers           = False  # mark scan positions

use_tick_labels        = True   # label q points (False: use raw q labels)
use_custom_labels      = True   # show custom dispersion labels (branch_labels)
use_colours            = True   # use dispersion colours (branch_colours)

num_Q_points           = 256    # number of Qs on a dispersion branch
S_scale                = 64.    # weight scaling and clamp factors
S_clamp_min            = 1.     #
S_clamp_max            = 500.   #
S_filter_min           = -1.    # don't filter

modelfile              = "model.magpie"
plotfile               = ""     # file to save plot to

print_progress         = False  # show progress of calculation

# cubic high-symmetry points
pt_G  = numpy.array([ 0.0, 0.0, 0.0 ])
pt_X1 = numpy.array([ 0.5, 0.0, 0.0 ])
pt_X2 = numpy.array([ 0.0, 0.5, 0.0 ])
pt_X3 = numpy.array([ 0.0, 0.0, 0.5 ])
pt_M1 = numpy.array([ 0.5, 0.5, 0.0 ])
pt_M2 = numpy.array([ 0.5, 0.0, 0.5 ])
pt_M3 = numpy.array([ 0.0, 0.5, 0.5 ])
pt_R  = numpy.array([ 0.5, 0.5, 0.5 ])

# dispersion branches to plot
dispersion     = [ pt_G, pt_X1, pt_M1, pt_R ]
branch_labels  = [ "Γ", "X", "M", "R" ]
branch_colours = [ "#ffffff", "#eeeeee", "#ffffff" ]

# scan markers
marker_colour = "#00000055"
marker_Qs     = [ [ ], [ 0.25 ], [ ] ]
marker_dQs    = [ 0.025, 0.025, 0.025 ]
marker_Es     = [ [ ], [ 9 ], [ ] ]
marker_dEs    = [ [ ], [ 7 ], [ ] ]
marker_labels = [ [ ], [ "(a)" ], [ ] ]

width_ratios = []                   # lengths from one dispersion point to the next
num_branches = len(dispersion) - 1  # number of dispersion braches

# determine number of threads to use
if max_threads == 0:
	max_threads = int(os.cpu_count() / 2)
if max_threads < 1:
	max_threads = 1
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# load the magnetic model
# -----------------------------------------------------------------------------
#
# calculate the energies and weights for a Q point
#
def calc_Es(mag, h, k, l):
	Es_dicts = mag.CalcEnergies(h, k, l, False)

	Es = []
	for Es_dict in Es_dicts:
		Es.append(( h, k, l, Es_dict.E, Es_dict.weight ))

	return Es


def setup_struct():
	# create the magpy object
	mag = magpy.MagDyn()

	# load the model file
	print("Loading {}...".format(modelfile))
	if not mag.Load(modelfile):
		print("Failed loading {}.".format(modelfile))
		exit(-1)

	# minimum energy
	print("\nEnergy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
	print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))

	return mag
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate the dispersion branches along straight lines
# -----------------------------------------------------------------------------
def calc_disp(mag):
	# calculate the dispersion branches
	print("\nCalculating %d dispersion branches..." % num_branches)
	if print_dispersion:
		print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

	dispersion_plot_indices = []  # Q components to use for plotting branches
	data = []                     # data for all branches

	for branch_idx in range(num_branches):
		hkl_start = dispersion[branch_idx]
		hkl_end = dispersion[branch_idx + 1]
		width_ratios.append(numpy.linalg.norm(hkl_end - hkl_start))

		print("[%d/%d] Calculating %s  ->  %s branch..." %
		   (branch_idx + 1, num_branches, hkl_start, hkl_end))

		# find scan axis
		Q_diff = [
			numpy.abs(hkl_start[0] - hkl_end[0]),
			numpy.abs(hkl_start[1] - hkl_end[1]),
			numpy.abs(hkl_start[2] - hkl_end[2]) ]

		axis_idx = 0
		if Q_diff[1] > Q_diff[axis_idx]:
			axis_idx = 1
		elif Q_diff[2] > Q_diff[axis_idx]:
			axis_idx = 2
		dispersion_plot_indices.append(axis_idx)


		#
		# add a data point
		#
		data_h, data_k, data_l = [], [], []
		data_E, data_S = [], []

		def append_data(h, k, l, E, S):
			#print(f"S(Q = ({h} {k} {l}) rlu, E = {E} meV) = {S}")
			weight = S * S_scale

			if weight < S_filter_min:
				return

			if weight < S_clamp_min:
				weight = S_clamp_min
			elif weight > S_clamp_max:
				weight = S_clamp_max

			data_h.append(h)
			data_k.append(k)
			data_l.append(l)
			data_E.append(E)
			data_S.append(weight)


		# calculate the dispersion
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

		data.append([ branch_idx, data_h, data_k, data_l, data_E, data_S ])

	return (data, dispersion_plot_indices)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# alternate way to calculate the dispersion branches in a point-wise fashion
# -----------------------------------------------------------------------------
def calc_disp_pointwise(mag):
	# calculate the dispersion branches
	print("\nCalculating %d dispersion branches..." % num_branches)
	if print_dispersion:
		print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

	dispersion_plot_indices = []  # Q components to use for plotting branches
	data = []                     # data for all branches

	for branch_idx in range(num_branches):
		hkl_start = dispersion[branch_idx]
		hkl_end = dispersion[branch_idx + 1]
		width_ratios.append(numpy.linalg.norm(hkl_end - hkl_start))

		print("[%d/%d] Calculating %s  ->  %s branch..." %
		   (branch_idx + 1, num_branches, hkl_start, hkl_end))

		# find scan axis
		Q_diff = [
			numpy.abs(hkl_start[0] - hkl_end[0]),
			numpy.abs(hkl_start[1] - hkl_end[1]),
			numpy.abs(hkl_start[2] - hkl_end[2]) ]

		axis_idx = 0
		if Q_diff[1] > Q_diff[axis_idx]:
			axis_idx = 1
		elif Q_diff[2] > Q_diff[axis_idx]:
			axis_idx = 2
		dispersion_plot_indices.append(axis_idx)


		#
		# add a data point
		#
		data_h, data_k, data_l = [], [], []
		data_E, data_S = [], []

		def append_data(h, k, l, E, S):
			#print(f"S(Q = ({h} {k} {l}) rlu, E = {E} meV) = {S}")
			weight = S * S_scale

			if weight < S_filter_min:
				return

			if weight < S_clamp_min:
				weight = S_clamp_min
			elif weight > S_clamp_max:
				weight = S_clamp_max

			data_h.append(h)
			data_k.append(k)
			data_l.append(l)
			data_E.append(E)
			data_S.append(weight)


		if use_threadpool:
			import concurrent.futures as fut

			if print_progress:
				print(f"Using {max_threads} threads.")

			executor = fut.ThreadPoolExecutor
			#executor = fut.ProcessPoolExecutor

			# calculate a branch
			with executor(max_workers = max_threads) as exe:
				Es_futures = []

				# submit tasks
				for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
					Es_futures.append(
						exe.submit(calc_Es, mag, hkl[0], hkl[1], hkl[2]))

				# print progress
				last_num_finished = -1
				while print_progress:
					num_futures = len(Es_futures)
					num_finished = len([future for future in Es_futures if future.done()])

					if last_num_finished != num_finished:
						print(f"{num_finished}/{num_futures} finished.")
						last_num_finished = num_finished

					if num_futures == num_finished:
						break
					time.sleep(0.01)

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

		data.append([ branch_idx, data_h, data_k, data_l, data_E, data_S ])

	return (data, dispersion_plot_indices)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the results
# -----------------------------------------------------------------------------
def plot_disp(data, dispersion_plot_indices):
	import matplotlib.pyplot as plot

	plot.rcParams.update({
		"font.sans-serif" : "DejaVu Sans",
		"font.family"     : "sans-serif",
		"font.size"       : 16,
		"axes.linewidth"  : 2.,
		"figure.figsize"  : (6.4, 4.8)
	})

	print("\nPlotting dispersion branches...")

	(plt, axes) = plot.subplots(nrows = 1, ncols = num_branches,
		width_ratios = width_ratios, sharey = True)

	# in case there's only one sub-plot
	if type(axes) != numpy.ndarray:
		axes = [ axes ]

	for ( branch_idx, data_h, data_k, data_l, data_E, data_S ) in data:
		# dispersion branch start and end points
		b1 = dispersion[branch_idx]
		b2 = dispersion[branch_idx + 1]

		plot_idx = dispersion_plot_indices[branch_idx]
		if plot_idx == 0:
			data_x = data_h
		elif plot_idx == 1:
			data_x = data_k
		elif plot_idx == 2:
			data_x = data_l

		# ticks and labels
		axes[branch_idx].set_xlim(data_x[0], data_x[-1])

		if use_colours:
			axes[branch_idx].set_facecolor(branch_colours[branch_idx])

		if use_tick_labels:
			if use_custom_labels:
				tick_labels = [
					branch_labels[branch_idx],
					branch_labels[branch_idx + 1] ]
			else:
				tick_labels = [
					"(%.4g %.4g %.4g)" % (b1[0], b1[1], b1[2]),
					"(%.4g %.4g %.4g)" % (b2[0], b2[1], b2[2]) ]

		if branch_idx == 0:
			axes[branch_idx].set_ylabel("E (meV)")
		else:
			axes[branch_idx].get_yaxis().set_visible(False)
			if not show_dividers:
				axes[branch_idx].spines["left"].set_visible(False)

			if use_tick_labels:
				tick_labels[0] = ""

		if not show_dividers and branch_idx != num_branches - 1:
			axes[branch_idx].spines["right"].set_visible(False)

		if use_tick_labels:
			axes[branch_idx].set_xticks([data_x[0], data_x[-1]], labels = tick_labels)

		if branch_idx == num_branches / 2 - 1:
			axes[branch_idx].set_xlabel("Q (rlu)")

		# plot the dispersion branch
		axes[branch_idx].scatter(data_x, data_E, marker = '.', s = data_S)

		if show_markers:
			for marker_Q, marker_E, marker_dE, marker_label in zip(marker_Qs[branch_idx], \
				marker_Es[branch_idx], marker_dEs[branch_idx], marker_labels[branch_idx]):
				# scan marker
				axes[branch_idx].add_patch(plot.Rectangle(
					xy = (marker_Q - marker_dQs[branch_idx]*0.5,
						marker_E - marker_dE*0.5),
					width = marker_dQs[branch_idx], height = marker_dE,
					color = marker_colour))

				# marker label
				axes[branch_idx].annotate(marker_label, xy = (
					marker_Q - marker_dQs[branch_idx]*0.5,
					marker_E - marker_dE*0.6))


	#plt.legend(handletextpad = 0.1)
	plt.tight_layout()
	plt.subplots_adjust(wspace = 0)

	if plotfile != "":
		plot.savefig(plotfile)
	plot.show()
# -----------------------------------------------------------------------------


if __name__ == "__main__":
	mag = setup_struct()

	tick = time.time()
	if use_pointwise:
		(data, dispersion_plot_indices) = calc_disp_pointwise(mag)
	else:
		(data, dispersion_plot_indices) = calc_disp(mag)
	time_needed = time.time() - tick
	print("Calculation took %.4g s." % time_needed)

	plot_disp(data, dispersion_plot_indices)
