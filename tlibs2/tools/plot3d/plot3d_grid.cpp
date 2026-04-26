/**
 * 3d grid data plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date 26 April 2026
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
#include "tlibs2/libs/instr/raw.h"



/**
 * calculate the surfaces via the given grid points
 */
void Plot3DDlg::CalculateGrid()
{
	// get data
	std::string all_data = m_formulas->toPlainText().toStdString();
	std::vector<std::string> datavec;
	tl2::get_tokens<std::string, std::string>(all_data, ";", datavec);

	std::vector<tl2::DatFile<t_real>> datvec;
	datvec.reserve(datavec.size());
	for(const std::string& dat : datavec)
	{
		std::istringstream istr(dat);
		tl2::DatFile<t_real> raw;
		if(raw.Load(istr))
			datvec.emplace_back(std::move(raw));
	}

	if(datvec.size() == 0)
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

	// get maximum grid size
	m_x_count = datvec[0].GetColumnCount();
	m_y_count = datvec[0].GetRowCount();
	for(t_size surf_idx = 1; surf_idx < datvec.size(); ++surf_idx)
	{
		m_x_count = std::max(datvec[surf_idx].GetColumnCount(), m_x_count);
		m_y_count = std::max(datvec[surf_idx].GetRowCount(), m_y_count);
	}
	//std::cout << "Grid: " << m_x_count << " x " << m_y_count << std::endl;

	t_vec_real start = tl2::create<t_vec_real>({ m_xrange[0]->value(), m_yrange[0]->value() });
	t_vec_real end = tl2::create<t_vec_real>({ m_xrange[1]->value(), m_yrange[1]->value() });

	m_minmax_x[0] = start[0];
	m_minmax_x[1] = end[0];
	m_minmax_y[0] = start[1];
	m_minmax_y[1] = end[1];

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

	m_data.resize(datvec.size());

	for(t_size x_idx = 0; x_idx < m_x_count; ++x_idx)
	for(t_size y_idx = 0; y_idx < m_y_count; ++y_idx)
	{
		auto task = [this, &datvec, &mtx, &start, &end, x_idx, y_idx]()
		{
			// calculate the surface at the given coordinate
			t_vec_real Q = start;
			Q[0] += (end[0] - start[0]) / t_real(m_x_count) * t_real(x_idx);
			Q[1] += (end[1] - start[1]) / t_real(m_y_count) * t_real(y_idx);

			// iterate the energies for this Q point
			for(t_size surf_idx = 0; surf_idx < datvec.size(); ++surf_idx)
			{
				bool valid = true;
				t_real z = datvec[surf_idx].GetData(y_idx, x_idx);
				if(std::isnan(z) || std::isinf(z))
					valid = false;

				if(valid)
				{
					m_minmax_z[0] = std::min(m_minmax_z[0], z);
					m_minmax_z[1] = std::max(m_minmax_z[1], z);
				}

				// generate and add data point
				t_data dat{std::make_tuple(Q, z, x_idx, y_idx, valid)};

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
