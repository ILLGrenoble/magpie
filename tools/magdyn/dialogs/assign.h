/**
 * magnon dynamics -- assigns multiple couplings
 * @author Tobias Weber <tweber@ill.fr>
 * @date 6-june-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core / magdyn / magpie
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

#ifndef __MAGDYN_ASSIGN_H__
#define __MAGDYN_ASSIGN_H__


#include <QtCore/QSettings>

#include <QtWidgets/QDialog>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QSpinBox>

#include "defs.h"



class AssignDlg : public QDialog
{ Q_OBJECT
public:
	AssignDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~AssignDlg() = default;

	AssignDlg(const AssignDlg&) = delete;
	const AssignDlg& operator=(const AssignDlg&) = delete;


private:
	QSettings *m_sett{};

	QLineEdit *m_editJ{};
	QLineEdit *m_editDMI[3]{};
	QLineEdit *m_editJs[3*3]{};

	QSpinBox *m_symmidx{};


protected:
	void AssignByIndex();


protected slots:
	virtual void accept() override;


signals:
	void AssignCouplingsBySymmetryIndex(t_size symmidx,
		const std::string* J, const std::string* DMI, const std::string* Js);
};


#endif
