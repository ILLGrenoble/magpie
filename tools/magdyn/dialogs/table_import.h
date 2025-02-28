/**
 * import magnetic structure from a table
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __MAGDYN_TABLE_IMPORT_H__
#define __MAGDYN_TABLE_IMPORT_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QSpinBox>

#include <optional>

#include "defs.h"



struct TableImportAtom
{
	std::string name{};
	std::string pos[3]{};
	std::string S[3]{}, Smag{};
};



struct TableImportCoupling
{
	std::string name{};
	std::optional<t_size> atomidx1{}, atomidx2{};
	std::string d[3]{};
	std::string J{}, dmi[3]{};
	std::string Jgen[3*3]{};
};



class TableImportDlg : public QDialog
{ Q_OBJECT
public:
	TableImportDlg(QWidget* parent = nullptr, QSettings* sett = nullptr);
	virtual ~TableImportDlg();

	TableImportDlg(const TableImportDlg&) = delete;
	const TableImportDlg& operator=(const TableImportDlg&) = delete;


protected:
	virtual void closeEvent(QCloseEvent *evt) override;

	void ImportAtoms();
	void ImportCouplings();

	void ShowHelp();

	bool HasSymmetricCoupling(const std::vector<TableImportCoupling>& couplings,
		const TableImportCoupling& coupling);


private:
	QSettings *m_sett{};

	// edit boxes
	QTextEdit *m_editAtoms{};
	QTextEdit *m_editCouplings{};

	// atom table column indices
	QSpinBox *m_spinAtomName{};
	QSpinBox *m_spinAtomX{}, *m_spinAtomY{}, *m_spinAtomZ{};
	QSpinBox *m_spinAtomSX{}, *m_spinAtomSY{}, *m_spinAtomSZ{};
	QSpinBox *m_spinAtomSMag{};

	// coupling table column
	QSpinBox *m_spinCouplingName{};
	QSpinBox *m_spinCouplingAtom1{}, *m_spinCouplingAtom2{};
	QSpinBox *m_spinCouplingDX{}, *m_spinCouplingDY{}, *m_spinCouplingDZ{};
	QSpinBox *m_spinCouplingJ{}, *m_spinCouplingJGen{};
	QSpinBox *m_spinCouplingDMIX{}, *m_spinCouplingDMIY{}, *m_spinCouplingDMIZ{};

	// options
	QCheckBox *m_checkIndices1Based{};
	QCheckBox *m_checkUniteIncompleteTokens{};
	QCheckBox *m_checkIgnoreSymmetricCoupling{};
	QCheckBox *m_checkClearExisting{};


signals:
	void SetAtomsSignal(const std::vector<TableImportAtom>&, bool clear_existing);
	void SetCouplingsSignal(const std::vector<TableImportCoupling>&, bool clear_existing);
};


#endif
