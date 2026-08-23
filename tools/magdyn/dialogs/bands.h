/**
 * choose magnon bands
 * @author Tobias Weber <tweber@ill.fr>
 * @date 23-august-2026
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

#ifndef __MAGDYN_BANDS_H__
#define __MAGDYN_BANDS_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QTableWidget>

#include <string>
#include <vector>



class BandsDlg : public QDialog
{ Q_OBJECT
public:
	BandsDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~BandsDlg() = default;

	BandsDlg(const BandsDlg&) = delete;
	const BandsDlg& operator=(const BandsDlg&) = delete;


private:
	QSettings *m_sett{};
	bool m_emit{true};

	QTableWidget *m_table_bands{};        // table listing the magnon bands
	std::vector<bool> m_enabled_bands{};  // previous band settings


public:
	void Clear(bool clear_settings = true);
	void Reset();
	
	bool GetOldChecked(std::size_t i) const;
	bool IsChecked(std::size_t i) const;
	void SetChecked(std::size_t i, bool active = true);
	void AddBand(const std::string& name, const QColor& colour, bool enabled = true);


protected slots:
	void EmitStateChanged();

	virtual void accept() override;


signals:
	void StateChanged();
};


#endif
