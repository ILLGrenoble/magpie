/**
 * 3d plotter -- file related functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "plot3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/ver.h"



/**
 * write the meta data header for text file exports
 */
void Plot3DDlg::WriteHeader(std::ostream& ostr) const
{
	using namespace tl2_ops;

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ostr << "#\n"
		<< "# Created by tlibs " << TLIBS2_VER << "\n"
		<< "# URL: https://github.com/ILLGrenoble/magpie\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.16180814\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n\n";
}



/**
 * save the surfaces as a text data file
 */
void Plot3DDlg::SaveData()
{
	bool skip_invalid_points = true;

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("plot3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Surface Data",
		dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("plot3d/dir", QFileInfo(filename).path());

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
	ofstr << std::setw(field_len) << std::left << "# x" << " ";
	ofstr << std::setw(field_len) << std::left << "y" << " ";
	ofstr << std::setw(field_len) << std::left << "z" << " ";
	if(!skip_invalid_points)
		ofstr << std::setw(field_len) << std::left << "valid" << " ";
	ofstr << std::setw(field_len) << std::left << "surf." << "\n";

	const t_size num_surfs = m_data.size();
	for(t_size surf_idx = 0; surf_idx < num_surfs; ++surf_idx)
	{
		for(t_data_Q& data : m_data[surf_idx])
		{
			const t_vec& xy = std::get<0>(data);
			t_real z = std::get<1>(data);
			bool valid = std::get<4>(data);

			if(skip_invalid_points && !valid)
				continue;

			ofstr << std::setw(field_len) << std::left << xy[0] << " ";
			ofstr << std::setw(field_len) << std::left << xy[1] << " ";
			ofstr << std::setw(field_len) << std::left << z << " ";
			if(!skip_invalid_points)
				ofstr << std::setw(field_len) << std::left << valid << " ";
			ofstr << std::setw(field_len) << std::left << surf_idx << "\n";
		}
	}

	ofstr.flush();
}



/**
 * save the surfaces as a script file
 */
void Plot3DDlg::SaveScript()
{
	using namespace tl2_ops;

	bool skip_invalid_points = true;

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("plot3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Surface Data As Script",
		dirLast, "Py Files (*.py)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("plot3d/dir", QFileInfo(filename).path());

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
plot_file = ""  # file to save plot to
# -----------------------------------------------------------------------------

# get the surface with a certain index
def get_surf(data, surf_data, z_surf_idx = 0):
	# iterate surfaces
	x_minmax = [ +999999999., -999999999. ]
	y_minmax = [ +999999999., -999999999. ]
	z_minmax = [ +999999999., -999999999. ]

	# filter data for given surface
	data_Q = [
		[ row[0] for (row, surf_idx) in zip(data, surf_data) \
			if surf_idx == z_surf_idx ],
		[ row[1] for (row, surf_idx) in zip(data, surf_data) \
			if surf_idx == z_surf_idx ]
	]
	data_E = [ row[2] for (row, surf_idx) in zip(data, surf_data) \
		if surf_idx == z_surf_idx ]

	if len(data_E) < 1:
		return None

	# data ranges
	x_minmax[0] = numpy.min([ numpy.min(data_Q[0]), x_minmax[0] ])
	x_minmax[1] = numpy.max([ numpy.max(data_Q[0]), x_minmax[1] ])
	y_minmax[0] = numpy.min([ numpy.min(data_Q[1]), y_minmax[0] ])
	y_minmax[1] = numpy.max([ numpy.max(data_Q[1]), y_minmax[1] ])
	z_minmax[0] = numpy.min([ numpy.min(data_E), z_minmax[0] ])
	z_minmax[1] = numpy.max([ numpy.max(data_E), z_minmax[1] ])

	return [data_Q, data_E, x_minmax, y_minmax, z_minmax]

# plot using mayavi
def plot_disp_mvi(data, surf_data, surf_colours):
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

	# iterate surfaces
	x_minmax = [ +999999999., -999999999. ]
	y_minmax = [ +999999999., -999999999. ]
	z_minmax = [ +999999999., -999999999. ]
	z_surf_eff_idx = 0           # effective index of actually plotted surfaces
	z_surf_max = max(surf_data)  # maximum surface index
	for z_surf_idx in range(0, z_surf_max + 1):
		surf = get_surf(data, surf_data, z_surf_idx)
		if surf == None:
			continue

		[data_Q, data_E, _x_minmax, _y_minmax, _z_minmax] = surf
		x_minmax[0] = min(x_minmax[0], _x_minmax[0])
		x_minmax[1] = max(x_minmax[1], _x_minmax[1])
		y_minmax[0] = min(y_minmax[0], _y_minmax[0])
		y_minmax[1] = max(y_minmax[1], _y_minmax[1])
		z_minmax[0] = min(z_minmax[0], _z_minmax[0])
		z_minmax[1] = max(z_minmax[1], _z_minmax[1])

		colour_idx = z_surf_eff_idx

		surf_repr = "surface" # "wireframe"
		zval = float(colour_idx) / float(z_surf_max)
		points = mlab.points3d(data_Q[0], data_Q[1], data_E, [zval]*len(data_E),
			mode = "point", opacity = 0.5, figure = fig, extent = axis_extents)
		triags = mlab.pipeline.delaunay2d(points, figure = fig)
		surface = mlab.pipeline.surface(triags, representation = surf_repr,
			figure = fig, extent = axis_extents, opacity = 1.,
			name = ("surf_%d" % z_surf_idx), line_width = 1.,
			color = conv_col(surf_colours[colour_idx]))
		z_surf_eff_idx += 1

	mlab.xlabel("x")
	mlab.ylabel("y")
	mlab.zlabel("z")
	axes = mlab.axes(figure = fig, ranges = [*x_minmax, *y_minmax, *z_minmax],
		extent = axis_extents, line_width = 2.)
	axes.axes.font_factor = 1.75
	mlab.outline(figure = fig, extent = axis_extents, line_width = 2.)

	if plot_file != "":
		mlab.savefig(plot_file, figure = fig)
	mlab.show()
	mlab.close(all = True)

# plot using matplotlib
def plot_disp_mpl(data, surf_data, surf_colours):
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

	# iterate surfaces
	x_minmax = [ +999999999., -999999999. ]
	y_minmax = [ +999999999., -999999999. ]
	z_minmax = [ +999999999., -999999999. ]
	z_surf_eff_idx = 0           # effective index of actually plotted surfaces
	z_surf_max = max(surf_data)  # maximum surface index
	for z_surf_idx in range(0, z_surf_max + 1):
		surf = get_surf(data, surf_data, z_surf_idx)
		if surf == None:
			continue

		[data_Q, data_E, _x_minmax, _y_minmax, _z_minmax] = surf
		x_minmax[0] = min(x_minmax[0], _x_minmax[0])
		x_minmax[1] = max(x_minmax[1], _x_minmax[1])
		y_minmax[0] = min(y_minmax[0], _y_minmax[0])
		y_minmax[1] = max(y_minmax[1], _y_minmax[1])
		z_minmax[0] = min(z_minmax[0], _z_minmax[0])
		z_minmax[1] = max(z_minmax[1], _z_minmax[1])

		colour_idx = z_surf_eff_idx

		axis.plot_trisurf(data_Q[0], data_Q[1], data_E,
			color = surf_colours[colour_idx], alpha = 1., shade = True,
			lightsource = light, antialiased = False,
			zorder = z_surf_max - z_surf_idx)
		z_surf_eff_idx += 1

	axis.set_xlabel("x", labelpad = 12)
	axis.set_ylabel("y", labelpad = 12)
	axis.set_zlabel("z", labelpad = 12)

	axis.set_xlim([x_minmax[0], x_minmax[1]])
	axis.set_ylim([y_minmax[0], y_minmax[1]])
	axis.set_zlim([z_minmax[0], z_minmax[1]])
	axis.set_box_aspect([ 1., 1., 1. ])

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0, bottom = 0.15)

	if plot_file != "":
		pyplot.savefig(plot_file)
	pyplot.show()
	pyplot.close()

if __name__ == "__main__":
	x_data = %%X_DATA%%
	y_data = %%Y_DATA%%
	z_data = %%Z_DATA%%
	surf_data = %%SURF_DATA%%
	colours = %%SURF_COLOURS%%

	data = numpy.array([ x_data, y_data, z_data ]).T

	try:
		# try to plot using mayavi...
		plot_disp_mvi(data, surf_data, colours)
	except ModuleNotFoundError:
		# ... otherwise resort to matplotlib
		plot_disp_mpl(data, surf_data, colours)
)RAW";
	// ------------------------------------------------------------------------

	const t_size num_surfs = m_data.size();

	// create data arrays
	std::ostringstream x_data, y_data, z_data, surfidx_data, colours;
	for(std::ostringstream* ostr : { &x_data, &y_data, &z_data, &surfidx_data, &colours })
	{
		ostr->precision(g_prec);
	}

	for(t_size surf_idx = 0; surf_idx < num_surfs; ++surf_idx)
	{
		// surface data
		for(t_data_Q& data : m_data[surf_idx])
		{
			const t_vec& xy = std::get<0>(data);
			t_real z = std::get<1>(data);
			bool valid = std::get<4>(data);

			if(skip_invalid_points && !valid)
				continue;

			x_data << xy[0] << ", ";
			y_data << xy[1] << ", ";
			z_data << z << ", ";

			surfidx_data << surf_idx << ", ";
		}

		// surface colour
		std::array<int, 3> col = GetBranchColour(surf_idx, num_surfs);
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
			<< g_fov << "/180.*numpy.pi),\n";
	}

	algo::replace_all(pyscr, "%%X_DATA%%", "[ " + x_data.str() + "]");
	algo::replace_all(pyscr, "%%Y_DATA%%", "[ " + y_data.str() + "]");
	algo::replace_all(pyscr, "%%Z_DATA%%", "[ " + z_data.str() + "]");
	algo::replace_all(pyscr, "%%SURF_DATA%%", "[ " + surfidx_data.str() + "]");
	algo::replace_all(pyscr, "%%SURF_COLOURS%%", "[ " + colours.str() + "]");
	algo::replace_all(pyscr, "%%PROJ_TYPE%%", m_perspective->isChecked() ? "persp" : "ortho");
	algo::replace_all(pyscr, "%%ORTHO_PROJ%%", m_perspective->isChecked() ? "False" : "True");
	algo::replace_all(pyscr, "%%FOCAL_LEN%%", focal_len.str());
	algo::replace_all(pyscr, "%%AZIMUTH%%", tl2::var_to_str(-90. - m_cam_phi->value(), g_prec));
	algo::replace_all(pyscr, "%%ELEVATION%%", tl2::var_to_str(90. + m_cam_theta->value(), g_prec));

	ofstr << pyscr << std::endl;
}



/**
 * save an image of the plot
 */
void Plot3DDlg::SaveImage()
{
	if(!m_dispplot)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("plot3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Plot Image",
		dirLast, "PNG Files (*.png)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("plot3d/dir", QFileInfo(filename).path());

	if(!m_dispplot->grabFramebuffer().save(filename, nullptr, 90))
		ShowError(QString("Could not save plot image to file \"%1\".").arg(filename));
}
