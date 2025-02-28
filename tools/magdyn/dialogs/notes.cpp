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

#include "notes.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>



// prefix to base64-encoded strings
static const std::string g_b64_prefix = "__base64__";



/**
 * set up the gui
 */
NotesDlg::NotesDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	setWindowTitle("Notes");
	setSizeGripEnabled(true);

	// notes/comments
	m_notes = new QPlainTextEdit(this);
	m_notes->setReadOnly(false);
	m_notes->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &NotesDlg::accept);

	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	grid->addWidget(new QLabel(QString("Comments / Notes:"), this), y++, 0, 1, 1);
	grid->addWidget(m_notes, y++, 0, 1, 4);
	grid->addWidget(btnbox, y++, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("notes/geo"))
			restoreGeometry(m_sett->value("notes/geo").toByteArray());
		else
			resize(500, 500);
        }
}



void NotesDlg::ClearNotes()
{
	m_notes->clear();
}



/**
 * set the notes string
 */
void NotesDlg::SetNotes(const std::string& notes)
{
	if(notes.starts_with(g_b64_prefix))
	{
		m_notes->setPlainText(QByteArray::fromBase64(notes.substr(g_b64_prefix.length()).data(),
			QByteArray::Base64Encoding | QByteArray::AbortOnBase64DecodingErrors));
	}
	else
	{
		m_notes->setPlainText(notes.data());
	}
}



/**
 * get the notes as string
 */
std::string NotesDlg::GetNotes() const
{
	return g_b64_prefix + m_notes->toPlainText().toUtf8().toBase64(
		QByteArray::Base64Encoding).constData();
}



/**
 * close the dialog
 */
void NotesDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("notes/geo", saveGeometry());
	}

	QDialog::accept();
}
