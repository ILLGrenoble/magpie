/**
 * brillouin zone tool -- 2d plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __BZCUT_H__
#define __BZCUT_H__

#ifndef BZ_USE_QT_SIGNALS
	// qt signals can't be emitted from a template class
	// TODO: remove this as soon as this is supported
	#include <boost/signals2/signal.hpp>
#endif

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGraphicsItem>
#include <QtGui/QMouseEvent>
#include <QtGui/QWheelEvent>
#include <QtWidgets/QApplication>

#include <sstream>
#include <tuple>
#include <array>



// --------------------------------------------------------------------------------
template<class t_vec, class t_real = typename t_vec::value_type>
class BZCutScene : public QGraphicsScene
{
public:
	BZCutScene(QWidget *parent = nullptr) : QGraphicsScene(parent)
	{
	}


	virtual ~BZCutScene()
	{
		ClearAll();
	}


	/**
	 * adds the line segments of a brillouin zone cut
	 * @arg [start, end, Q]
	 */
	void AddCut(const std::vector<
		// [x, y, Q]
		std::tuple<t_vec, t_vec, std::array<t_real, 3>>>& lines)
	{
		// (000) brillouin zone
		std::vector<const std::tuple<t_vec, t_vec, std::array<t_real, 3>>*> lines000;

		QPen pen;
		pen.setCosmetic(true);
		pen.setColor(qApp->palette().color(QPalette::WindowText));
		pen.setWidthF(1.);

		m_bzcut.reserve(lines.size() * 2);

		// draw brillouin zones
		for(const auto& line : lines)
		{
			const auto& Q = std::get<2>(line);

			// (000) BZ?
			if(tl2::equals_0(Q[0], m_eps) &&
				tl2::equals_0(Q[1], m_eps) &&
				tl2::equals_0(Q[2], m_eps))
			{
				lines000.push_back(&line);
				continue;
			}

			QGraphicsLineItem *plot_line = addLine(QLineF(
				std::get<0>(line)[0]*m_scale, std::get<0>(line)[1]*m_scale,
				std::get<1>(line)[0]*m_scale, std::get<1>(line)[1]*m_scale),
				pen);

			m_bzcut.push_back(plot_line);
		}


		// draw (000) brillouin zone
		pen.setColor(QColor(0xff, 0x00, 0x00));
		pen.setWidthF(2.);

		for(const auto* line : lines000)
		{
			QGraphicsLineItem *plot_line = addLine(QLineF(
				std::get<0>(*line)[0]*m_scale, std::get<0>(*line)[1]*m_scale,
				std::get<1>(*line)[0]*m_scale, std::get<1>(*line)[1]*m_scale),
				pen);

			m_bzcut.push_back(plot_line);
		}
	}


	/**
	 * adds bragg peaks
	 */
	void AddPeaks(const std::vector<t_vec>& peaks, const std::vector<t_vec>* peaks_rlu = nullptr)
	{
		QColor col(0x00, 0x99, 0x00);

		QPen pen;
		pen.setCosmetic(true);
		pen.setColor(col);

		QBrush brush;
		brush.setStyle(Qt::SolidPattern);
		brush.setColor(col);

		m_peaks.reserve(peaks.size());
		const t_real w = 6.;

		for(std::size_t Q_idx = 0; Q_idx < peaks.size(); ++Q_idx)
		{
			const t_vec& Q = peaks[Q_idx];

			QGraphicsEllipseItem *ell = addEllipse(
				Q[0]*m_scale - w/2., Q[1]*m_scale - w/2., w, w,
				pen, brush
			);

			if(peaks_rlu)
			{
				const t_vec& Q_rlu = (*peaks_rlu)[Q_idx];

				std::ostringstream ostr;
				ostr.precision(m_prec_gui);

				ostr << "(" << Q_rlu[0] << " " << Q_rlu[1] << " " << Q_rlu[2] << ")";
				ell->setToolTip(ostr.str().c_str());
			}

			m_peaks.push_back(ell);
		}
	}


	/**
	 * adds a plot curve from a set of points
	 */
	void AddCurve(const std::vector<t_vec>& points)
	{
		if(points.size() < 2)
			return;

		QPen pen;
		pen.setCosmetic(true);
		pen.setColor(QColor(0x00, 0x00, 0xff));
		pen.setWidthF(2.);

		for(std::size_t i = 0; i < points.size() - 1; ++i)
		{
			std::size_t j = i + 1;

			QGraphicsLineItem *plot_line = addLine(QLineF(
				points[i][0]*m_scale, points[i][1]*m_scale,
				points[j][0]*m_scale, points[j][1]*m_scale),
				pen);

			m_curves.push_back(plot_line);
		}
	}


	void ClearAll()
	{
		ClearCut();
		ClearPeaks();
		ClearCurves();
		clear();
	}


	void ClearCut()
	{
		for(QGraphicsItem *item : m_bzcut)
			delete item;
		m_bzcut.clear();
	}


	void ClearPeaks()
	{
		for(QGraphicsItem *item : m_peaks)
			delete item;
		m_peaks.clear();
	}


	void ClearCurves()
	{
		for(QGraphicsItem *item : m_curves)
			delete item;
		m_curves.clear();
	}


	/**
	 * get the centre of the brillouin zone
	 */
	QPointF GetCentre() const
	{
		QPointF centre{0, 0};

		for(const QGraphicsItem* item : m_bzcut)
			centre += item->pos();
		centre /= qreal(m_bzcut.size());

		return centre;
	}


	void SetEps(t_real eps)
	{
		m_eps = eps;
	}


	void SetPrecGui(int prec)
	{
		m_prec_gui = prec;
	}


	void SetScale(t_real scale)
	{
		m_scale = scale;
	}


	t_real GetScale() const
	{
		return m_scale;
	}


protected:
	t_real m_eps{ 1e-7 };
	int m_prec_gui{ 4 };

	t_real m_scale{ 100. };

	std::vector<QGraphicsItem*> m_bzcut{};
	std::vector<QGraphicsItem*> m_peaks{};
	std::vector<QGraphicsItem*> m_curves{};
};
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
template<class t_vec, class t_real = typename t_vec::value_type>
class BZCutView : public QGraphicsView
{
#ifdef BZ_USE_QT_SIGNALS
	Q_OBJECT
#endif
public:
	BZCutView(BZCutScene<t_vec, t_real>* scene)
		: QGraphicsView(scene, static_cast<QWidget*>(scene->parent())),
		  m_scene(scene)
	{
		setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
		setDragMode(QGraphicsView::ScrollHandDrag);
		setInteractive(true);
		setMouseTracking(true);
		scale(1., -1.);
	}


	BZCutView(const BZCutView<t_vec, t_real>&) = delete;
	BZCutView<t_vec, t_real>& operator=(const BZCutView<t_vec, t_real>&) = delete;


	virtual ~BZCutView()
	{
	}


	/**
	 * centre the view
	 */
	void Centre()
	{
		QPointF centre{0, 0};
		if(m_scene)
			centre = m_scene->GetCentre();

		centerOn(centre);
	}


protected:
	virtual void mouseMoveEvent(QMouseEvent *evt) override
	{
		QPointF pos = mapToScene(evt->pos());
		t_real scale = m_scene->GetScale();
#ifdef BZ_USE_QT_SIGNALS
		emit SignalMouseCoordinates(pos.x()/scale, pos.y()/scale);
#else
		m_sigMouseCoordinates(pos.x()/scale, pos.y()/scale);
#endif

		QGraphicsView::mouseMoveEvent(evt);
	}


	virtual void wheelEvent(QWheelEvent *evt) override
	{
		t_real sc = std::pow(2., evt->angleDelta().y()/8.*0.01);
		QGraphicsView::scale(sc, sc);
	}


#ifdef BZ_USE_QT_SIGNALS
signals:
	void SignalMouseCoordinates(t_real x, t_real y);
#else
	boost::signals2::signal<void(t_real, t_real)> m_sigMouseCoordinates{ };

public:
	template<class t_slot>
	boost::signals2::connection AddMouseCoordinatesSlot(const t_slot& slot)
	{
		return m_sigMouseCoordinates.connect(slot);
	}
#endif


protected:
	BZCutScene<t_vec, t_real>* m_scene{ };
};
// --------------------------------------------------------------------------------


#endif
