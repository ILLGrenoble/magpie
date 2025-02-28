/**
 * peak finder test
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-2021
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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
 *
 * g++-10 -std=c++20 -I.. -o peaks peaks.cpp
 */

#include "libs/maths.h"
#include "libs/fit.h"


using t_real = double;


int main()
{
	const int spline_order = 3;
	const int spline_pts = 512;

	std::vector<t_real> x = tl2::linspace<t_real>(0, 360., 256);

	std::vector<t_real> y, y_neg;
	y.reserve(x.size());
	y_neg.reserve(x.size());

	for(t_real val : x)
	{
		t_real val_fkt = std::sin(val / 180. * tl2::pi<t_real>);

		y.push_back(val_fkt);
		y_neg.push_back(-val_fkt);
	}


	std::vector<t_real> peaks_x, peaks_sizes, peaks_widths;
	std::vector<bool> peaks_minima;

	tl2::find_peaks(x.size(), x.data(), y.data(), spline_order,
		peaks_x, peaks_sizes, peaks_widths, peaks_minima,
		spline_pts, 1e-8);

	std::cout << "local extrema:" << std::endl;
	for(std::size_t i=0; i<peaks_x.size(); ++i)
	{
		if(peaks_minima[i])
			std::cout << "\tminimum: ";
		else
			std::cout << "\tmaximum: ";

		std::cout
			<< "x = " << peaks_x[i]
			<< ", h = " << peaks_sizes[i]
			<< ", w = " << peaks_widths[i]
			<< std::endl;
	}

	return 0;
}
