/**
 * 3d plotter calculation
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

#include <QtWidgets/QApplication>

#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <cstdlib>
#include <mutex>
#include <sstream>

#include "plot3d.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * calculate the surfaces via the given formulas
 */
void Plot3DDlg::Calculate()
{
	// get formulas
	QString _formulas = m_formulas->toPlainText();

	// if no variables are given, treat it as data grid
	if(_formulas.indexOf("x") == -1 && _formulas.indexOf("y") == -1)
	{
		CalculateGrid();
		return;
	}

	std::string all_formulas = _formulas.toStdString();
	std::vector<std::string> formulas;
	tl2::get_tokens<std::string, std::string>(all_formulas, ";", formulas);
	if(formulas.size() == 0)
		return;

	// parse formulas
	std::vector<tl2::ExprParser<t_real>> parsers;
	for(const std::string& formula : formulas)
	{
		tl2::ExprParser<t_real> parser;

		parser.SetAutoregisterVariables(false);
		parser.register_var("x", 0.);
		parser.register_var("y", 0.);
		parser.register_var("t", 0.);

		if(!parser.parse_noexcept(formula))
		{
			std::cerr << "Error parsing formula \"" << formula << "\"." << std::endl;
			continue;
		}

		parsers.emplace_back(std::move(parser));
	}

	if(parsers.size() == 0)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	m_minmax_z[0] = +std::numeric_limits<t_real>::max();
	m_minmax_z[1] = -std::numeric_limits<t_real>::max();
	m_minmax_x[0] = m_minmax_x[1] = 0.;
	m_minmax_y[0] = m_minmax_y[1] = 0.;
	m_data.clear();

	m_x_count = m_num_points[0]->value();
	m_y_count = m_num_points[1]->value();

	t_vec_real start = tl2::create<t_vec_real>({ (t_real)m_xrange[0]->value(), (t_real)m_yrange[0]->value() });
	t_vec_real end = tl2::create<t_vec_real>({ (t_real)m_xrange[1]->value(), (t_real)m_yrange[1]->value() });

	m_minmax_x[0] = start[0];
	m_minmax_x[1] = end[0];
	m_minmax_y[0] = start[1];
	m_minmax_y[1] = end[1];

	// get t parameter
	int t_idx = m_slider_t->value();
	t_real t = tl2::lerp(m_trange[0]->value(), m_trange[1]->value(),
		t_real(t_idx) / t_real(m_slider_t->maximum() - m_slider_t->minimum()));

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stop_requested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(m_x_count * m_y_count);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(m_x_count * m_y_count);

	m_data.resize(parsers.size());

	for(t_size x_idx = 0; x_idx < m_x_count; ++x_idx)
	for(t_size y_idx = 0; y_idx < m_y_count; ++y_idx)
	{
		auto task = [this, &parsers, &mtx, &start, &end, x_idx, y_idx, t]()
		{
			// calculate the surface at the given coordinate
			t_vec_real pos = start;
			pos[0] += (end[0] - start[0]) / t_real(m_x_count) * t_real(x_idx);
			pos[1] += (end[1] - start[1]) / t_real(m_y_count) * t_real(y_idx);

			// iterate the z components for this position
			for(t_size surf_idx = 0; surf_idx < parsers.size(); ++surf_idx)
			{
				tl2::ExprParser<t_real> localparser = parsers[surf_idx];
				localparser.register_var("x", pos[0]);
				localparser.register_var("y", pos[1]);
				localparser.register_var("t", t);

				bool valid = true;
				t_real z = localparser.eval();
				if(std::isnan(z) || std::isinf(z))
					valid = false;

				if(valid)
				{
					m_minmax_z[0] = std::min(m_minmax_z[0], z);
					m_minmax_z[1] = std::max(m_minmax_z[1], z);
				}

				// generate and add data point
				t_data dat{std::make_tuple(pos, z, x_idx, y_idx, valid)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[surf_idx].emplace_back(std::move(dat));
			}
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		// process events to see if the stop button was clicked
		// only do this for a fraction of the points to avoid gui overhead
		if(task_idx % std::max<t_size>(tasks.size() / std::sqrt(g_stop_check_fraction), 1) == 0)
			qApp->processEvents();

		if(m_stop_requested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	// finish parallel calculations
	pool.join();

	// get sorting of data
	for(t_size surf_idx = 0; surf_idx < m_data.size(); ++surf_idx)
	{
		std::vector<std::size_t> perm = tl2::get_perm(m_data[surf_idx].size(),
			[this, surf_idx](std::size_t idx1, std::size_t idx2) -> bool
		{
			// sorting by indices
			t_size x_idx_1 = std::get<2>(m_data[surf_idx][idx1]);
			t_size x_idx_2 = std::get<3>(m_data[surf_idx][idx1]);
			t_size y_idx_1 = std::get<2>(m_data[surf_idx][idx2]);
			t_size y_idx_2 = std::get<3>(m_data[surf_idx][idx2]);

			if(x_idx_1 != y_idx_1)
				return x_idx_1 < y_idx_1;
			return x_idx_2 < y_idx_2;
		});

		m_data[surf_idx] = tl2::reorder(m_data[surf_idx], perm);
	}

	// remove fully invalid surfaces
	for(auto iter = m_data.begin(); iter != m_data.end();)
	{
		if(!IsValid(*iter))
			iter = m_data.erase(iter);
		else
			++iter;
	}

	// sort surfaces in descending order
	std::stable_sort(m_data.begin(), m_data.end(), [this](const t_data_vec& dat1, const t_data_vec& dat2)
	{
		return GetMeanZ(dat1) >= GetMeanZ(dat2);
	});

	// show elapsed time
	stopwatch.stop();
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stop_requested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	Plot(true);
}



/**
 * determine if a surface is valid or only contains invalid points
 */
bool Plot3DDlg::IsValid(const t_data_vec& data) const
{
	for(const t_data& data : data)
	{
		// data point valid?
		if(std::get<4>(data))
			return true;
	}

	// all points invalid
	return false;
}



/**
 * count the number of valid and invalid points in a surface
 */
std::pair<t_size, t_size> Plot3DDlg::NumValid(const t_data_vec& data) const
{
	t_size valid = 0, invalid = 0;

	for(const t_data& data : data)
	{
		// data point valid?
		if(std::get<4>(data))
			++valid;
		else
			++invalid;
	}

	return std::make_pair(valid, invalid);
}



/**
 * calculate the mean surface z value
 */
t_real Plot3DDlg::GetMeanZ(const Plot3DDlg::t_data_vec& data) const
{
	t_real E_mean = 0.;
	t_size num_pts = 0;

	for(const t_data& data : data)
	{
		// data point invalid?
		if(!std::get<4>(data))
			continue;

		E_mean += std::get<1>(data);
		++num_pts;
	}

	if(num_pts)
		E_mean /= static_cast<t_real>(num_pts);
	return E_mean;
}



/**
 * calculate the mean surface z value
 */
t_real Plot3DDlg::GetMeanZ(t_size surf_idx) const
{
	if(surf_idx >= m_data.size())
		return 0.;

	return GetMeanZ(m_data[surf_idx]);
}



/**
 * clears the table of surfaces
 */
void Plot3DDlg::ClearSurfaces()
{
	m_surf_objs.clear();
	m_cur_obj = std::nullopt;

	m_table_surfs->clearContents();
	m_table_surfs->setRowCount(0);
}



/**
 * adds a surface to the table
 */
void Plot3DDlg::AddSurface(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_surfs)
		return;

	int row = m_table_surfs->rowCount();
	m_table_surfs->insertRow(row);

	QTableWidgetItem *item = new QTableWidgetItem{name.c_str()};
	item->setFlags(item->flags() & ~Qt::ItemIsEditable);

	QBrush bg = item->background();
	bg.setColor(colour);
	bg.setStyle(Qt::SolidPattern);
	item->setBackground(bg);

	QBrush fg = item->foreground();
	fg.setColor(QColor{0xff, 0xff, 0xff});
	fg.setStyle(Qt::SolidPattern);
	item->setForeground(fg);

	QCheckBox *checkSurf = new QCheckBox(m_table_surfs);
	checkSurf->setChecked(enabled);
	connect(checkSurf, &QCheckBox::toggled, [this]() { Plot(false); });

	m_table_surfs->setItem(row, COL_SURF, item);
	m_table_surfs->setCellWidget(row, COL_ACTIVE, checkSurf);
}



/**
 * verifies if the surface's checkbox is checked
 */
bool Plot3DDlg::IsSurfaceEnabled(t_size idx) const
{
	if(!m_table_surfs || int(idx) >= m_table_surfs->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(
		m_table_surfs->cellWidget(int(idx), COL_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}
