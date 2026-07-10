/**
 * choose matrix elements
 * @author Tobias Weber <tweber@ill.fr>
 * @date 10-july-2026
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
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

#ifndef __MAGDYN_MATRIXELEMS_H__
#define __MAGDYN_MATRIXELEMS_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QCheckBox>

#include <string>



class MatrixElemsDlg : public QDialog
{ Q_OBJECT
public:
	MatrixElemsDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~MatrixElemsDlg() = default;

	MatrixElemsDlg(const MatrixElemsDlg&) = delete;
	const MatrixElemsDlg& operator=(const MatrixElemsDlg&) = delete;


private:
	QSettings *m_sett{};
	QCheckBox *m_elems[2*3*3]{};
	bool m_emit{true};


public:
	void Reset();
	bool IsChecked(std::size_t i, std::size_t j, bool real_elem = true) const;


protected slots:
	void EmitStateChanged();
	virtual void accept() override;


signals:
	void StateChanged();
};


#endif
