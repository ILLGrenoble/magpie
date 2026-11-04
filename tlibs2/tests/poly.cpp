/**
 * 2d polygon test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 1-nov-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
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
 */

// g++ -std=c++20 -O2 -march=native -I.. -I/opt/homebrew/include -o poly poly.cpp

#include <iostream>

#include "libs/maths.h"
#include "libs/cam.h"


using namespace tl2_ops;


void poly1()
{
	using t_int = int;
	using t_vec = tl2::vec_static<t_int, 2>;

	std::vector<t_vec> poly
	{
		tl2::create<t_vec>({ 10, 1 }),
		tl2::create<t_vec>({ 30, 1 }),
		tl2::create<t_vec>({ 20, 8 }),
	};

	for(int y = 0; y < 10; ++y)
	{
		for(int x = 0; x < 50; ++x)
		{
			t_vec pt = tl2::create<t_vec>({ x, y });

			bool inside = tl2::in_poly_2d(pt, poly);
			std::cout << (inside ? "*" : " ");
		}

		std::cout << std::endl;
	}
}


void poly2()
{
	using t_real = double;
	using t_vec3 = tl2::vec_static<t_real, 3>;
	using t_vec4 = tl2::vec_static<t_real, 4>;
	using t_mat44 = tl2::mat_static<t_real, 4, 4>;

	int w = 60, h = 20;

	tl2::Camera<t_mat44, t_vec4, t_vec3, t_real> cam;
	cam.SetScreenDimensions(w, h);
	cam.SetRotation(0., 0.);
	cam.SetPosition(tl2::create<t_vec3>({ 0., 0., 0. }));
	cam.SetDist(0.);
	cam.UpdateViewport();
	cam.UpdatePerspective();
	cam.UpdateTransformation();

	std::cout << "viewp: " << cam.GetViewport() << std::endl;
	std::cout << "perps: " << cam.GetPerspective() << std::endl;
	std::cout << "trafo: " << cam.GetTransformation() << std::endl;

	char* screen = new char[w*h];
	t_real* zbuf = new t_real[w*h];
	for(int y = 0.; y < h; ++y)
	{
		for(int x = 0.; x < w; ++x)
		{
			screen[y*w + x] = ' ';
			zbuf[y*w + x] = 999.;
		}
	}


	auto draw = [w, h, screen, zbuf, &cam](const std::vector<t_vec4>& poly, char sym)
	{
		std::vector<t_vec4> poly_proj = cam.ToScreenCoords(poly);
		std::cout << poly_proj[0] << std::endl;
		std::cout << poly_proj[1] << std::endl;
		std::cout << poly_proj[2] << std::endl;

		t_real z = poly_proj[0][2];

		for(int y = 0.; y < h; ++y)
		{
			for(int x = 0.; x < w; ++x)
			{
				if(x == 0 || y == 0 || x == w - 1 || y == h - 1)
					screen[y*w + x] = '#';

				t_vec4 pt = tl2::create<t_vec4>({ t_real(x), t_real(y), 0., 0. });

				bool inside = tl2::in_poly_2d(pt, poly_proj, false, &z);
				//std::cout << z << " " << poly_proj[0][2] << std::endl;
				if(inside && z < zbuf[y*w + x])
				{
					screen[y*w + x] = sym;
					zbuf[y*w + x] = z;
				}
			}
		}
	};

	std::vector<t_vec4> poly
	{
		tl2::create<t_vec4>({ -10., 5., -15., 1. }),
		tl2::create<t_vec4>({ 10., 5., -15., 1. }),
		tl2::create<t_vec4>({ 0., -5., -15., 1. }),
	};
	std::vector<t_vec4> poly2
	{
		tl2::create<t_vec4>({ -16., -9., -16., 1. }),
		tl2::create<t_vec4>({ 16., -9., -16., 1. }),
		tl2::create<t_vec4>({ 0., 10., -14., 1. }),
	};
	std::vector<t_vec4> poly3
	{
		tl2::create<t_vec4>({ -20., -13., -17., 1. }),
		tl2::create<t_vec4>({ -20., 13., -17., 1. }),
		tl2::create<t_vec4>({ 20., 13., -17., 1. }),
		tl2::create<t_vec4>({ 20., -13., -17., 1. }),
	};

	draw(poly, '*');
	draw(poly2, '.');
	draw(poly3, '-');

	for(int y = 0.; y < h; ++y)
	{
		for(int x = 0.; x < w; ++x)
			std::cout << screen[y*w + x];
		std::cout << std::endl;
	}

	delete[] screen;
	delete[] zbuf;
}


void poly3()
{
	using t_real = double;
	using t_vec3 = tl2::vec_static<t_real, 3>;
	using t_vec4 = tl2::vec_static<t_real, 4>;
	using t_mat33 = tl2::mat_static<t_real, 3, 3>;
	using t_mat44 = tl2::mat_static<t_real, 4, 4>;

	int w = 80, h = 40;

	tl2::Camera<t_mat44, t_vec4, t_vec3, t_real> cam;
	cam.SetScreenDimensions(w, h);
	cam.SetRotation(0., tl2::pi<t_real>/4.9);
	cam.SetPosition(tl2::create<t_vec3>({ 0., 0., 0. }));
	cam.SetDist(12.);
	cam.UpdateViewport();
	cam.UpdatePerspective();
	cam.UpdateTransformation();

	std::cout << "viewp: " << cam.GetViewport() << std::endl;
	std::cout << "perps: " << cam.GetPerspective() << std::endl;
	std::cout << "trafo: " << cam.GetTransformation() << std::endl;

	char* screen = new char[w*h];
	t_real* zbuf = new t_real[w*h];
	for(int y = 0.; y < h; ++y)
	{
		for(int x = 0.; x < w; ++x)
		{
			screen[y*w + x] = ' ';
			zbuf[y*w + x] = 999.;
		}
	}


	auto draw = [w, h, screen, zbuf, &cam](const std::vector<t_vec4>& poly, char sym)
	{
		std::vector<t_vec4> poly_proj = cam.ToScreenCoords(poly);
		std::cout << poly_proj[0] << std::endl;
		std::cout << poly_proj[1] << std::endl;
		std::cout << poly_proj[2] << std::endl;

		t_real z = poly_proj[0][2];

		for(int y = 0.; y < h; ++y)
		{
			for(int x = 0.; x < w; ++x)
			{
				if(x == 0 || y == 0 || x == w - 1 || y == h - 1)
					screen[y*w + x] = '#';

				t_vec4 pt = tl2::create<t_vec4>({ t_real(x), t_real(y), 0., 0. });

				bool inside = tl2::in_poly_2d(pt, poly_proj, false, &z);
				//std::cout << z << " " << poly_proj[0][2] << std::endl;
				if(inside && z < zbuf[y*w + x])
				{
					screen[y*w + x] = sym;
					zbuf[y*w + x] = z;
				}
			}
		}
	};

	auto solid = tl2::create_cuboid<t_vec3>(5., 5., 5.);

	char fill[] = { '*', '-', '~', '+', '.', '|' };
	std::size_t face_idx = 0;
	for(const auto& face : std::get<1>(solid))
	{
		std::vector<t_vec4> poly;
		for(std::size_t vert_idx : face)
		{
			t_vec3 vert = std::get<0>(solid)[vert_idx];
			vert = tl2::givens<t_mat33>(3, 2, 0, tl2::pi<t_real>/4.) * vert;
			poly.push_back(tl2::create<t_vec4>({ vert[0], vert[1], vert[2], 1. }));
		}

		draw(poly, fill[face_idx % sizeof(fill)]);
		++face_idx;
	}

	for(int y = 0.; y < h; ++y)
	{
		for(int x = 0.; x < w; ++x)
			std::cout << screen[y*w + x];
		std::cout << std::endl;
	}

	delete[] screen;
	delete[] zbuf;
}


int main()
{
	poly1();
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	poly2();
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	poly3();

	return 0;
}
