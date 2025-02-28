/**
 * magnon dynamics -- notes
 * @author Tobias Weber <tweber@ill.fr>
 * @date june-2024
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

#ifndef __MAGDYN_NOTES_H__
#define __MAGDYN_NOTES_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QPlainTextEdit>

#include <string>



class NotesDlg : public QDialog
{ Q_OBJECT
public:
	NotesDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~NotesDlg() = default;

	NotesDlg(const NotesDlg&) = delete;
	const NotesDlg& operator=(const NotesDlg&) = delete;

	void ClearNotes();
	void SetNotes(const std::string& notes);
	std::string GetNotes() const;


private:
	QSettings *m_sett{};
	QPlainTextEdit *m_notes{};


protected slots:
	virtual void accept() override;
};


#endif
