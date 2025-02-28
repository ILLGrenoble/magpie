/**
 * magnetic dynamics -- minimise ground state energy
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAG_DYN_GROUNDSTATE_DLG_H__
#define __MAG_DYN_GROUNDSTATE_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>

#include <thread>
#include <atomic>
#include <memory>

#include "gui_defs.h"




/**
 * ground state dialog
 */
class GroundStateDlg : public QDialog
{ Q_OBJECT
public:
	GroundStateDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~GroundStateDlg();

	GroundStateDlg(const GroundStateDlg&) = delete;
	GroundStateDlg& operator=(const GroundStateDlg&) = delete;

	void SetKernel(const t_magdyn* dyn);
	void SyncFromKernel(const t_magdyn *dyn = nullptr,
		const std::unordered_set<std::string> *fixed_spins = nullptr);


protected:
	void SyncToKernel();

	void EnableMinimisation(bool enable = true);
	void Minimise();
	void CalcGroundStateEnergy();

	void ShowError(const char* msg);

	void SpinsTableItemChanged(QTableWidgetItem *item);
	void UpdateSpinFromTable(int row = -1);
	void UpdateSpinsFromTable();

	virtual void accept() override;


private:
	const t_magdyn *m_dyn_kern{};     // main kernel
	std::optional<t_magdyn> m_dyn{};  // local copy to work on

	QSettings *m_sett{};
	QTableWidget *m_spinstab{};
	QPushButton *m_btnFromKernel{}, *m_btnToKernel{};
	QPushButton *m_btnMinimise{};
	QLabel *m_status{};

	std::unique_ptr<std::thread> m_thread{};  // minimiser thread
	bool m_stop_request{false};               // stop minimisation
	std::atomic<bool> m_running{false};       // is minimisation running?


signals:
	void SpinsUpdated(const t_magdyn* dyn);
};


#endif
