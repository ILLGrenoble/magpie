/**
 * magnon dynamics -- info dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date june-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * magpie & mag-core
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

#ifndef __MAGDYN_INFOS_H__
#define __MAGDYN_INFOS_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>



class InfoDlg : public QDialog
{ Q_OBJECT
public:
	InfoDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~InfoDlg() = default;

	InfoDlg(const InfoDlg&) = delete;
	const InfoDlg& operator=(const InfoDlg&) = delete;


private:
	QSettings *m_sett{};


protected slots:
	virtual void accept() override;
};


#endif
