/**
 * magnon dynamics -- transformation calculator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 29-dec-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 7-sep-2023 from my privately developed "gl" project: https://github.com/t-weber/gl .
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "gl" project
 * Copyright (C) 2021-2023  Tobias WEBER (privately developed).
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

#ifndef __MAGDYN_TRAFOCALC_H__
#define __MAGDYN_TRAFOCALC_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QDoubleSpinBox>



class TrafoCalculator : public QDialog
{ Q_OBJECT
public:
	TrafoCalculator(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~TrafoCalculator() = default;

	TrafoCalculator(const TrafoCalculator&) = delete;
	const TrafoCalculator& operator=(const TrafoCalculator&) = delete;


private:
	QSettings *m_sett{};

	QTextEdit *m_textRotation{};
	QDoubleSpinBox *m_spinAxis[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_spinAngle{};
	QDoubleSpinBox *m_spinVecToRotate[3]{nullptr, nullptr, nullptr};


protected slots:
	virtual void accept() override;

	void CalculateRotation();
};


#endif
