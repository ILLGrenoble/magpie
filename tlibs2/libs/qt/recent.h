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

#ifndef __TL2_RECENT_FILES_H__
#define __TL2_RECENT_FILES_H__

#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>

#include <functional>


namespace tl2 {


class RecentFiles
{
public:
	RecentFiles() = default;
	~RecentFiles();

	RecentFiles(const RecentFiles& other) = default;
	RecentFiles& operator=(const RecentFiles& other) = default;

	void Clear();

	// adds a file to the recent files menu
	void AddRecentFile(const QString &file);

	// sets the recent file menu
	void SetRecentFiles(const QStringList &files);
	const QStringList& GetRecentFiles() const { return m_recentFiles; }

	// creates the "recent files" sub-menu
	void RebuildRecentFiles();

	// remove superfluous entries
	void TrimEntries();

	// set the function to be called when the menu element is clicked
	void SetOpenFunc(const std::function<bool(const QString& filename)>* func)
	{ m_open_func = func; }

	void SetRecentFilesMenu(QMenu *menu) { m_menuOpenRecent = menu; };
	QMenu* GetRecentFilesMenu() { return m_menuOpenRecent; }

	void SetCurFile(const QString& file) { m_curFile = file; }
	const QString& GetCurFile() const { return m_curFile; }

	void AddForbiddenDir(const QString& dir) { m_forbiddenDirs << dir; }
	bool IsFileInForbiddenDir(const QString& file) const;

	void SetMaxRecentFiles(std::size_t num) { m_maxRecentFiles = num; }


private:
	// maximum number of recent files
	std::size_t m_maxRecentFiles{ 16 };

	// recent file menu
	QMenu* m_menuOpenRecent{};

	// recent file list
	QStringList m_recentFiles{};

	// currently active file
	QString m_curFile{};

	// list of directories for which recent files shouldn't be saved
	QStringList m_forbiddenDirs{};

	// function to be called when the menu element is clicked
	const std::function<bool(const QString& filename)>* m_open_func{};
};


}
#endif
