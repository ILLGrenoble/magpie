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

#include <QtCore/QSettings>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGraphicsItem>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMenu>
#include <QtGui/QMouseEvent>
#include <QtGui/QWheelEvent>
#include <QtSvg/QSvgGenerator>

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

		m_bzcut.reserve(m_bzcut.size() + lines.size() * 2);

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

		m_peaks.reserve(m_peaks.size() + peaks.size());
		const t_real w = 6.;

		for(std::size_t Q_idx = 0; Q_idx < peaks.size(); ++Q_idx)
		{
			const t_vec& Q = peaks[Q_idx];

			QGraphicsEllipseItem *ell = addEllipse(
				Q[0]*m_scale - w/2., Q[1]*m_scale - w/2., w, w,
				pen, brush);

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

		m_curves.reserve(m_curves.size() + points.size() - 1);

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


	/**
	 * adds a line from a start to an end point
	 */
	void AddLine(const t_vec& start, const t_vec& end,
		bool add_verts = false,
		const std::string& line_tooltip = "",
		const std::string& start_tooltip = "",
		const std::string& end_tooltip = "")
	{
		// draw line
		QPen pen;
		pen.setCosmetic(true);
		pen.setColor(QColor(0x00, 0x00, 0xff));
		pen.setWidthF(2.);

		QGraphicsLineItem *plot_line = addLine(QLineF(
			start[0]*m_scale, start[1]*m_scale,
			end[0]*m_scale, end[1]*m_scale),
			pen);
		plot_line->setToolTip(line_tooltip.c_str());

		m_lines.reserve(m_lines.size() + 1);
		m_lines.push_back(plot_line);

		// draw vertices
		const t_real w = 6.;
		if(add_verts)
		{
			QColor colStart(0xff, 0x00, 0x00);
			QColor colEnd(0x00, 0x00, 0xff);

			QPen penStart, penEnd;
			penStart.setCosmetic(true);
			penStart.setColor(colStart);
			penEnd.setCosmetic(true);
			penEnd.setColor(colEnd);

			QBrush brushStart, brushEnd;
			brushStart.setStyle(Qt::SolidPattern);
			brushStart.setColor(colStart);
			brushEnd.setStyle(Qt::SolidPattern);
			brushEnd.setColor(colEnd);

			QGraphicsEllipseItem *ellStart = addEllipse(
				start[0]*m_scale - w/2., start[1]*m_scale - w/2., w, w,
				penStart, brushStart);
			QGraphicsEllipseItem *ellEnd = addEllipse(
				end[0]*m_scale - w/2., end[1]*m_scale - w/2., w, w,
				penEnd, brushEnd);

			ellStart->setToolTip(start_tooltip.c_str());
			ellEnd->setToolTip(end_tooltip.c_str());

			m_line_verts.reserve(m_line_verts.size() + 2);
			m_line_verts.push_back(ellStart);
			m_line_verts.push_back(ellEnd);
		}
	}


	void ClearAll()
	{
		ClearCut();
		ClearPeaks();
		ClearCurves();
		ClearLines();
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


	void ClearLines()
	{
		for(QGraphicsItem *item : m_lines)
			delete item;
		for(QGraphicsItem *item : m_line_verts)
			delete item;
		m_lines.clear();
		m_line_verts.clear();
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


	/**
	 * save the current plot as image
	 */
	bool SaveImage(const QString& filename)
	{
		if(filename == "")
			return false;

		QSvgGenerator svg;
		//svg.setViewBox(sceneRect());
		svg.setSize(QSize{ int(sceneRect().width()), int(sceneRect().height()) });
		svg.setFileName(filename);
		svg.setTitle("Brillouin zone cut.");
		svg.setDescription("Created with Magpie (https://doi.org/10.5281/zenodo.16180814).");

		QPainter painter;
		painter.begin(&svg);
		render(&painter);
		painter.end();

		return true;
	}


private:
	t_real m_eps{ 1e-7 };
	int m_prec_gui{ 4 };

	t_real m_scale{ 100. };

	std::vector<QGraphicsItem*> m_bzcut{}, m_peaks{};
	std::vector<QGraphicsItem*> m_curves{};
	std::vector<QGraphicsItem*> m_lines{}, m_line_verts{};
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
	BZCutView(BZCutScene<t_vec, t_real>* scene, QSettings *sett = nullptr)
		: QGraphicsView(scene, static_cast<QWidget*>(scene->parent())),
		  m_scene{scene}, m_sett{sett}
	{
		// context menu
		m_context = new QMenu(this);
		QAction *acSaveImage = new QAction("Save Image...", m_context);
		acSaveImage->setIcon(QIcon::fromTheme("image-x-generic"));
		m_context->addAction(acSaveImage);

		// connections
		connect(acSaveImage, &QAction::triggered, this, &BZCutView::SaveImage);

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


	QMenu* GetContextMenu()
	{
		return m_context;
	}


	const t_vec& GetClickedPosition(bool right_button = false) const
	{
		return right_button ? m_cur_pos[1] : m_cur_pos[0];
	}


	/**
	 * save the current plot as image
	 */
	bool SaveImage()
	{
		if(!m_scene)
			return false;

		QString dirLast;
		if(m_sett)
			dirLast = m_sett->value("bz2d/dir", "").toString();
		QString filename = QFileDialog::getSaveFileName(
			this, "Save Reciprocal Space Image",
			dirLast, "SVG Files (*.svg)");
		if(filename == "")
			return false;

		if(!m_scene->SaveImage(filename))
			return false;

		if(m_sett)
			m_sett->setValue("bz2d/dir", QFileInfo(filename).path());
		return true;
	}


protected:
	virtual void mouseMoveEvent(QMouseEvent *evt) override
	{
		QPointF pos = mapToScene(evt->pos());
		t_real scale = m_scene ? m_scene->GetScale() : 1.;

#ifdef BZ_USE_QT_SIGNALS
		emit SignalMouseCoordinates(pos.x()/scale, pos.y()/scale);
#else
		m_sigMouseCoordinates(pos.x()/scale, pos.y()/scale);
#endif

		QGraphicsView::mouseMoveEvent(evt);
	}


	virtual void mousePressEvent(QMouseEvent *evt) override
	{
		QPointF pos = mapToScene(evt->pos());
		t_real scale = m_scene ? m_scene->GetScale() : 1.;

		int buttons = 0;
		if(evt->buttons() & Qt::LeftButton)
			buttons |= 1;
		if(evt->buttons() & Qt::MiddleButton)
			buttons |= 2;
		if(evt->buttons() & Qt::RightButton)
			buttons |= 4;

		// show context menu
		if(buttons & 1)  // left button
		{
			m_cur_pos[0] = tl2::create<t_vec>({ pos.x() / scale, pos.y() / scale });
		}

		if(buttons & 4)  // right button
		{
			m_cur_pos[1] = tl2::create<t_vec>({ pos.x() / scale, pos.y() / scale });

			QPointF _pt{ pos.x(), pos.y() };
			QPoint pt = mapToGlobal(mapFromScene(_pt));

			m_context->popup(pt);
		}

#ifdef BZ_USE_QT_SIGNALS
		emit SignalClickCoordinates(buttons, pos.x()/scale, pos.y()/scale);
#else
		m_sigClickCoordinates(buttons, pos.x()/scale, pos.y()/scale);
#endif

		QGraphicsView::mousePressEvent(evt);
	}


	virtual void wheelEvent(QWheelEvent *evt) override
	{
		t_real sc = std::pow(2., evt->angleDelta().y()/8.*0.01);
		QGraphicsView::scale(sc, sc);
	}


public:
#ifdef BZ_USE_QT_SIGNALS
signals:
	void SignalMouseCoordinates(t_real x, t_real y);
	void SignalClickCoordinates(int buttons, t_real x, t_real y);

#else

	boost::signals2::signal<void(t_real, t_real)> m_sigMouseCoordinates{ };
	boost::signals2::signal<void(int, t_real, t_real)> m_sigClickCoordinates{ };

	template<class t_slot>
	boost::signals2::connection AddMouseCoordinatesSlot(const t_slot& slot)
	{
		return m_sigMouseCoordinates.connect(slot);
	}


	template<class t_slot>
	boost::signals2::connection AddClickCoordinatesSlot(const t_slot& slot)
	{
		return m_sigClickCoordinates.connect(slot);
	}
#endif


private:
	BZCutScene<t_vec, t_real>* m_scene{ };
	QSettings *m_sett{};

	QMenu *m_context{};     // right-click context menu
	t_vec m_cur_pos[2]{};   // current position when clicking on plot
};
// --------------------------------------------------------------------------------


#endif
