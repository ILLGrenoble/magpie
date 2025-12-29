/**
 * magnetic dynamics -- 3d dispersion plot -- file related functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <QtWidgets/QFileDialog>

#include "dispersion3d.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * write the meta data header for text file exports
 */
void Dispersion3DDlg::WriteHeader(std::ostream& ostr) const
{
	using namespace tl2_ops;

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	const t_size num_bands = m_data.size();
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	ostr << "#\n"
		<< "# Created by Magpie " << MAGPIE_VER << "\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n# branch_count: " << num_bands << "\n"
		<< "# Q_origin: " << Q_origin << "\n"
		<< "# Q_direction_1: " << Q_dir_1 << "\n"
		<< "# Q_direction_2: " << Q_dir_2 << "\n"
		<< "#\n\n";
}



/**
 * save the dispersion as a text data file
 */
void Dispersion3DDlg::SaveData()
{
	bool skip_invalid_points = true;
	bool use_weights = m_S_filter_enable->isChecked();

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("dispersion3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Dispersion Data",
		dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("dispersion3d/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename));
		return;
	}

	ofstr.precision(g_prec);
	WriteHeader(ofstr);

	// write column header
	int field_len = g_prec * 2.5;
	ofstr << std::setw(field_len) << std::left << "# h" << " ";
	ofstr << std::setw(field_len) << std::left << "k" << " ";
	ofstr << std::setw(field_len) << std::left << "l" << " ";
	ofstr << std::setw(field_len) << std::left << "E" << " ";
	if(!skip_invalid_points)
		ofstr << std::setw(field_len) << std::left << "valid" << " ";
	if(use_weights)
		ofstr << std::setw(field_len) << std::left << "S" << " ";
	ofstr << std::setw(field_len) << std::left << "band" << " ";
	ofstr << std::setw(field_len) << std::left << "Qidx1" << " ";
	ofstr << std::setw(field_len) << std::left << "Qidx2" << " ";
	ofstr << std::setw(field_len) << std::left << "degen" << "\n";

	const t_size num_bands = m_data.size();
	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		for(t_data_Q& data : m_data[band_idx])
		{
			const t_vec_real& Q = std::get<0>(data);
			t_real E = std::get<1>(data);
			t_real S = std::get<2>(data);
			t_size Qidx1 = std::get<3>(data);
			t_size Qidx2 = std::get<4>(data);
			t_size degen = std::get<5>(data);
			bool valid = std::get<6>(data);

			if(skip_invalid_points && !valid)
				continue;

			ofstr << std::setw(field_len) << std::left << Q[0] << " ";
			ofstr << std::setw(field_len) << std::left << Q[1] << " ";
			ofstr << std::setw(field_len) << std::left << Q[2] << " ";
			ofstr << std::setw(field_len) << std::left << E << " ";
			if(!skip_invalid_points)
				ofstr << std::setw(field_len) << std::left << valid << " ";
			if(use_weights)
				ofstr << std::setw(field_len) << std::left << S << " ";
			ofstr << std::setw(field_len) << std::left << band_idx << " ";
			ofstr << std::setw(field_len) << std::left << Qidx1 << " ";
			ofstr << std::setw(field_len) << std::left << Qidx2 << " ";
			ofstr << std::setw(field_len) << std::left << degen << "\n";
		}
	}

	ofstr.flush();
}



/**
 * save the dispersion as a script file
 */
void Dispersion3DDlg::SaveScript()
{
	using namespace tl2_ops;

	bool skip_invalid_points = true;
	bool use_weights = m_S_filter_enable->isChecked();

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("dispersion3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Dispersion Data As Script",
		dirLast, "Py Files (*.py)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("dispersion3d/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename));
		return;
	}

	ofstr.precision(g_prec);
	WriteHeader(ofstr);

	// ------------------------------------------------------------------------
	// plot script
	std::string pyscr = R"RAW(import sys
import numpy

# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
plot_file      = ""     # file to save plot to
only_pos_E     = %%ONLY_POS_E%%  # ignore magnon annihilation?
S_filter_min   = 1e-5   # cutoff minimum spectral weight

h_column       =  0     # h column index in data files
k_column       =  1     # k column index
l_column       =  2     # l column index
E_column       =  3     # E column index
S_column       = %%S_INDEX%%     # S column index, -1: not available
# -----------------------------------------------------------------------------

# get the magnons with a certain branch index
def get_branch(data, branch_data, E_branch_idx = 0, Q_idx1 = 0, Q_idx2 = 1):
	# iterate energy branches
	Q1_minmax = [ +999999999., -999999999. ]
	Q2_minmax = [ +999999999., -999999999. ]
	E_minmax  = [ +999999999., -999999999. ]

	# filter data for given branch
	data_S = []
	data_Q = [
		[ row[h_column + Q_idx1] for (row, branch_idx) in zip(data, branch_data) \
			if branch_idx == E_branch_idx ],
		[ row[h_column + Q_idx2] for (row, branch_idx) in zip(data, branch_data) \
			if branch_idx == E_branch_idx ]
	]
	data_E = [ row[E_column] for (row, branch_idx) in zip(data, branch_data) \
		if branch_idx == E_branch_idx ]
	if S_column >= 0:
		data_S = [ row[S_column] for (row, branch_idx) in zip(data, branch_data) \
			if branch_idx == E_branch_idx ]

	if only_pos_E:
		# ignore magnon annihilation
		data_Q[0] = [ Q for (Q, E) in zip(data_Q[0], data_E) if E >= 0. ]
		data_Q[1] = [ Q for (Q, E) in zip(data_Q[1], data_E) if E >= 0. ]
		if S_column >= 0:
			data_S = [ S for (S, E) in zip(data_S, data_E) if E >= 0. ]
		data_E = [ E for E in data_E if E >= 0. ]

	if S_column >= 0 and S_filter_min >= 0.:
		# filter weights below cutoff
		data_Q[0] = [ Q for (Q, S) in zip(data_Q[0], data_S) if S >= S_filter_min ]
		data_Q[1] = [ Q for (Q, S) in zip(data_Q[1], data_S) if S >= S_filter_min ]
		data_E = [ E for (E, S) in zip(data_E, data_S) if S >= S_filter_min ]
		data_S = [ S for S in data_S if S >= S_filter_min ]

	if len(data_E) < 1:
		return None

	# data ranges
	Q1_minmax[0] = numpy.min([ numpy.min(data_Q[0]), Q1_minmax[0] ])
	Q1_minmax[1] = numpy.max([ numpy.max(data_Q[0]), Q1_minmax[1] ])
	Q2_minmax[0] = numpy.min([ numpy.min(data_Q[1]), Q2_minmax[0] ])
	Q2_minmax[1] = numpy.max([ numpy.max(data_Q[1]), Q2_minmax[1] ])
	E_minmax[0] = numpy.min([ numpy.min(data_E), E_minmax[0] ])
	E_minmax[1] = numpy.max([ numpy.max(data_E), E_minmax[1] ])

	return [data_Q, data_E, data_S, Q1_minmax, Q2_minmax, E_minmax]

# plot using mayavi
def plot_disp_mvi(data, branch_data, degen_data, branch_colours, Q_idx1 = 0, Q_idx2 = 1):
	# convert a colour from a string like "#ff0000" to a tuple like (1, 0., 0)
	def conv_col(col_str):
		r = float(int(col_str[1:3], 16)) / 255.
		g = float(int(col_str[3:5], 16)) / 255.
		b = float(int(col_str[5:7], 16)) / 255.
		return (r, g, b)

	from mayavi import mlab
	fig = mlab.figure(size = (800, 600), fgcolor = (0., 0., 0.), bgcolor = (1., 1., 1.))
	mlab.view(azimuth = -(%%AZIMUTH%%) + 270., elevation = %%ELEVATION%% + 90., figure = fig)
	fig.scene.parallel_projection = %%ORTHO_PROJ%%
	axis_extents = [0., 1., 0., 1., 0., 1.]

	# iterate energy branches
	Q1_minmax = [ +999999999., -999999999. ]
	Q2_minmax = [ +999999999., -999999999. ]
	E_minmax  = [ +999999999., -999999999. ]
	E_branch_eff_idx = 0             # effective index of actually plotted bands
	E_branch_max = max(branch_data)  # maximum branch index
	for E_branch_idx in range(0, E_branch_max + 1):
		branch = get_branch(data, branch_data, E_branch_idx, Q_idx1, Q_idx2)
		if branch == None:
			continue

		[data_Q, data_E, data_S, _Q1_minmax, _Q2_minmax, _E_minmax] = branch
		Q1_minmax[0] = min(Q1_minmax[0], _Q1_minmax[0])
		Q1_minmax[1] = max(Q1_minmax[1], _Q1_minmax[1])
		Q2_minmax[0] = min(Q2_minmax[0], _Q2_minmax[0])
		Q2_minmax[1] = max(Q2_minmax[1], _Q2_minmax[1])
		E_minmax[0] = min(E_minmax[0], _E_minmax[0])
		E_minmax[1] = max(E_minmax[1], _E_minmax[1])

		colour_idx = E_branch_eff_idx
		if only_pos_E:
			colour_idx *= 2  # skip every other colour if E >= 0

		surf_repr = "surface" # "wireframe"
		zval = float(colour_idx) / float(E_branch_max)
		points = mlab.points3d(data_Q[0], data_Q[1], data_E, [zval]*len(data_E),
			mode = "point", opacity = 0.5, figure = fig, extent = axis_extents)
		triags = mlab.pipeline.delaunay2d(points, figure = fig)
		surface = mlab.pipeline.surface(triags, representation = surf_repr,
			figure = fig, extent = axis_extents, opacity = 1.,
			name = ("branch_%d" % E_branch_idx), line_width = 1.,
			color = conv_col(branch_colours[colour_idx]))
		E_branch_eff_idx += 1

	labels = [ "h (rlu)", "k (rlu)", "l (rlu)" ]
	mlab.xlabel(labels[Q_idx1])
	mlab.ylabel(labels[Q_idx2])
	mlab.zlabel("E (meV)")
	axes = mlab.axes(figure = fig, ranges = [*Q1_minmax, *Q2_minmax, *E_minmax],
		extent = axis_extents, line_width = 2.)
	axes.axes.font_factor = 1.75
	mlab.outline(figure = fig, extent = axis_extents, line_width = 2.)

	if plot_file != "":
		mlab.savefig(plot_file, figure = fig)
	mlab.show()
	mlab.close(all = True)

# plot using matplotlib
def plot_disp_mpl(data, branch_data, degen_data, branch_colours, Q_idx1 = 0, Q_idx2 = 1):
	from matplotlib import colors
	from matplotlib import pyplot
	pyplot.rcParams.update({
		"font.sans-serif" : "DejaVu Sans",
		"font.family" : "sans-serif",
		"font.size" : 16,
	})

	(plt, axis) = pyplot.subplots(nrows = 1, ncols = 1,
		width_ratios = None, sharey = True,
		subplot_kw = { "projection" : "3d", "proj_type" : "%%PROJ_TYPE%%",
			%%FOCAL_LEN%% "computed_zorder" : False,
			"azim" : %%AZIMUTH%%, "elev" : %%ELEVATION%% })
	light = colors.LightSource(azdeg = %%AZIMUTH%%, altdeg = 45.)

	# iterate energy branches
	Q1_minmax = [ +999999999., -999999999. ]
	Q2_minmax = [ +999999999., -999999999. ]
	E_minmax  = [ +999999999., -999999999. ]
	E_branch_eff_idx = 0             # effective index of actually plotted bands
	E_branch_max = max(branch_data)  # maximum branch index
	for E_branch_idx in range(0, E_branch_max + 1):
		branch = get_branch(data, branch_data, E_branch_idx, Q_idx1, Q_idx2)
		if branch == None:
			continue

		[data_Q, data_E, data_S, _Q1_minmax, _Q2_minmax, _E_minmax] = branch
		Q1_minmax[0] = min(Q1_minmax[0], _Q1_minmax[0])
		Q1_minmax[1] = max(Q1_minmax[1], _Q1_minmax[1])
		Q2_minmax[0] = min(Q2_minmax[0], _Q2_minmax[0])
		Q2_minmax[1] = max(Q2_minmax[1], _Q2_minmax[1])
		E_minmax[0] = min(E_minmax[0], _E_minmax[0])
		E_minmax[1] = max(E_minmax[1], _E_minmax[1])

		colour_idx = E_branch_eff_idx
		if only_pos_E:
			colour_idx *= 2  # skip every other colour if E >= 0

		axis.plot_trisurf(data_Q[0], data_Q[1], data_E,
			color = branch_colours[colour_idx], alpha = 1., shade = True,
			lightsource = light, antialiased = False,
			zorder = E_branch_max - E_branch_idx)
		E_branch_eff_idx += 1

	labels = [ "h (rlu)", "k (rlu)", "l (rlu)" ]
	axis.set_xlabel(labels[Q_idx1], labelpad = 12)
	axis.set_ylabel(labels[Q_idx2], labelpad = 12)
	axis.set_zlabel("E (meV)", labelpad = 12)

	axis.set_xlim([Q1_minmax[0], Q1_minmax[1]])
	axis.set_ylim([Q2_minmax[0], Q2_minmax[1]])
	axis.set_zlim([E_minmax[0], E_minmax[1]])
	axis.set_box_aspect([ 1., 1., 1. ])

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0, bottom = 0.15)

	if plot_file != "":
		pyplot.savefig(plot_file)
	pyplot.show()
	pyplot.close()

if __name__ == "__main__":
	h_data = %%H_DATA%%
	k_data = %%K_DATA%%
	l_data = %%L_DATA%%
	E_data = %%E_DATA%%
	S_data = %%S_DATA%%
	branch_data = %%BRANCH_DATA%%
	degen_data = %%DEGEN_DATA%%
	colours = %%BRANCH_COLOURS%%

	if len(S_data) == len(E_data):
		data = numpy.array([ h_data, k_data, l_data, E_data, S_data ]).T
	else:
		data = numpy.array([ h_data, k_data, l_data, E_data ]).T

	try:
		# try to plot using mayavi...
		plot_disp_mvi(data, branch_data, degen_data, colours, %%Q_IDX_1%%, %%Q_IDX_2%%)
	except ModuleNotFoundError:
		# ... otherwise resort to matplotlib
		plot_disp_mpl(data, branch_data, degen_data, colours, %%Q_IDX_1%%, %%Q_IDX_2%%)
)RAW";
	// ------------------------------------------------------------------------

	const t_size num_bands = m_data.size();

	// create data arrays
	std::ostringstream h_data, k_data, l_data, E_data, S_data, bandidx_data, degen_data, colours;
	for(std::ostringstream* ostr : { &h_data, &k_data, &l_data, &E_data, &S_data,
		&bandidx_data, &degen_data, &colours })
	{
		ostr->precision(g_prec);
	}

	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		// band data
		for(t_data_Q& data : m_data[band_idx])
		{
			const t_vec_real& Q = std::get<0>(data);
			t_real E = std::get<1>(data);
			t_real S = std::get<2>(data);
			//t_size Qidx1 = std::get<3>(data);
			//t_size Qidx2 = std::get<4>(data);
			t_size degen = std::get<5>(data);
			bool valid = std::get<6>(data);

			if(skip_invalid_points && !valid)
				continue;

			h_data << Q[0] << ", ";
			k_data << Q[1] << ", ";
			l_data << Q[2] << ", ";
			E_data << E << ", ";
			if(use_weights)
				S_data << S << ", ";

			bandidx_data << band_idx << ", ";
			degen_data << degen << ", ";
		}

		// band colour
		std::array<int, 3> col = GetBranchColour(band_idx, num_bands);
		colours << "\"#" << std::hex
			<< std::setw(2) << std::setfill('0') << col[0]
			<< std::setw(2) << std::setfill('0') << col[1]
			<< std::setw(2) << std::setfill('0') << col[2]
			<< "\", ";
	}

	std::ostringstream focal_len;
	focal_len.precision(g_prec);
	if(m_perspective->isChecked())
	{
		// see: https://en.wikipedia.org/wiki/Focal_length
		focal_len << "\"focal_length\" : 1. / numpy.tan(0.5 * "
			<< g_structplot_fov << "/180.*numpy.pi),\n";
	}

	auto [Q_idx_1, Q_idx_2] = GetQIndices();

	// TODO: for the moment the azimuth only matches for Q directions in a right-handed system -> add a trafo
	algo::replace_all(pyscr, "%%H_DATA%%", "[ " + h_data.str() + "]");
	algo::replace_all(pyscr, "%%K_DATA%%", "[ " + k_data.str() + "]");
	algo::replace_all(pyscr, "%%L_DATA%%", "[ " + l_data.str() + "]");
	algo::replace_all(pyscr, "%%E_DATA%%", "[ " + E_data.str() + "]");
	algo::replace_all(pyscr, "%%S_DATA%%", "[ " + S_data.str() + "]");
	algo::replace_all(pyscr, "%%S_INDEX%%", use_weights ? " 4" : "-1");
	algo::replace_all(pyscr, "%%BRANCH_DATA%%", "[ " + bandidx_data.str() + "]");
	algo::replace_all(pyscr, "%%DEGEN_DATA%%", "[ " + degen_data.str() + "]");
	algo::replace_all(pyscr, "%%BRANCH_COLOURS%%", "[ " + colours.str() + "]");
	algo::replace_all(pyscr, "%%ONLY_POS_E%%", m_only_pos_E->isChecked() ? "True " : "False");
	algo::replace_all(pyscr, "%%PROJ_TYPE%%", m_perspective->isChecked() ? "persp" : "ortho");
	algo::replace_all(pyscr, "%%ORTHO_PROJ%%", m_perspective->isChecked() ? "False" : "True");
	algo::replace_all(pyscr, "%%FOCAL_LEN%%", focal_len.str());
	algo::replace_all(pyscr, "%%AZIMUTH%%", tl2::var_to_str(-90. - m_cam_phi->value(), g_prec));
	algo::replace_all(pyscr, "%%ELEVATION%%", tl2::var_to_str(90. + m_cam_theta->value(), g_prec));
	algo::replace_all(pyscr, "%%Q_IDX_1%%", tl2::var_to_str(Q_idx_1, g_prec));
	algo::replace_all(pyscr, "%%Q_IDX_2%%", tl2::var_to_str(Q_idx_2, g_prec));

	ofstr << pyscr << std::endl;
}



/**
 * save an image of the plot
 */
void Dispersion3DDlg::SaveImage()
{
	if(!m_dispplot)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("dispersion3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Plot Image",
		dirLast, "PNG Files (*.png)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("dispersion3d/dir", QFileInfo(filename).path());

	if(!m_dispplot->grabFramebuffer().save(filename, nullptr, 90))
		ShowError(QString("Could not save plot image to file \"%1\".").arg(filename));
}
