/**
 * recent files
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2021
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * TAS-Paths (part of the Takin software suite)
 * Copyright (C) 2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                     Grenoble, France).
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

#include "recent.h"

#include <QtCore/QObject>
#include <QtCore/QFile>
#include <QtCore/QDir>


namespace tl2 {


RecentFiles::~RecentFiles()
{
	Clear();
}


void RecentFiles::Clear()
{
	if(!m_menuOpenRecent)
		return;

	for(QAction *action : m_menuOpenRecent->actions())
		delete action;

	m_menuOpenRecent->clear();
}


/**
 * adds a file to the recent files menu
 */
void RecentFiles::AddRecentFile(const QString &file)
{
	// remove the file if it's already in the list
	if(int idx = m_recentFiles.indexOf(file); idx >= 0 && idx < m_recentFiles.size())
		m_recentFiles.removeAt(idx);

	m_recentFiles.push_back(file);
	RebuildRecentFiles();
}


/**
 * sets the recent file menu
 */
void RecentFiles::SetRecentFiles(const QStringList &files)
{
	m_recentFiles = files;
	RebuildRecentFiles();
}


/**
 * creates the "recent files" sub-menu
 */
void RecentFiles::RebuildRecentFiles()
{
	Clear();

	std::size_t num_recent_files = 0;
	for(auto iter = m_recentFiles.rbegin(); iter != m_recentFiles.rend();)
	{
		QString filename = *iter;

		// remove recent file entries which do not exist anymore
		if(!QFile::exists(filename) || IsFileInForbiddenDir(filename))
		{
			// current iterator position
			std::ptrdiff_t positer = iter - m_recentFiles.rbegin();

			// get corresponding forward iterator
			auto iter_fwd = (iter+1).base();
			iter_fwd = m_recentFiles.erase(iter_fwd);

			// restore iterator position
			iter = m_recentFiles.rbegin() + positer;

			// more files to check?
			if(iter != m_recentFiles.rend())
				continue;
			break;
		}

		auto *acFile = new QAction(
			QIcon::fromTheme("document"), filename, m_menuOpenRecent);

		// action when clicked on file name
		QObject::connect(acFile, &QAction::triggered, [this, filename]()
		{
			if(!m_open_func)
				return;

			if((*m_open_func)(filename))
			{
				AddRecentFile(filename);  // move file to the front of the sub-menu
				m_curFile = filename;
			}
		});

		m_menuOpenRecent->addAction(acFile);

		if(++num_recent_files >= m_maxRecentFiles)
			break;

		++iter;
	}
}


/**
 * remove superfluous entries
 */
void RecentFiles::TrimEntries()
{
	m_recentFiles.removeDuplicates();

	// remove superfluous entries and save the recent files list
	while((std::size_t)m_recentFiles.size() > m_maxRecentFiles)
		m_recentFiles.pop_front();
}


/**
 * tests if the file is in the list of forbidden directories
 */
bool RecentFiles::IsFileInForbiddenDir(const QString& file) const
{
	QString file_abs = QDir(file).absolutePath();

	for(const QString& dir : m_forbiddenDirs)
	{
		QString dir_abs = QDir(dir).absolutePath();

		// is file path shorter than directory?
		if(file_abs.length() < dir_abs.length())
			continue;

		if(std::equal(dir_abs.begin(), dir_abs.end(), file_abs.begin()))
			return true;
	}

	return false;
}


}
