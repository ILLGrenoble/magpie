/**
 * magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2026
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

 #ifndef __MAGDYN_SITE_ITEM_H__
 #define __MAGDYN_SITE_ITEM_H__


#include <QtWidgets/QComboBox>
#include <QtWidgets/QTableWidgetItem>

#include <iostream>



/**
 * a table widget item with an associated combo box
 */
struct SitesWidgetItem : public QTableWidgetItem
{
	SitesWidgetItem(QComboBox *combo) : m_combo{combo}
	{
	}


	virtual ~SitesWidgetItem() = default;

	SitesWidgetItem(const SitesWidgetItem&) = delete;
	SitesWidgetItem& operator=(const SitesWidgetItem&) = delete;


	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		const QComboBox* combo = dynamic_cast<const SitesWidgetItem*>(&item)->GetComboBox();
		if(!combo || !GetComboBox())
		{
			std::cerr << "Error: Sites table widget item does not have an associated combo box." << std::endl;
			return true;
		}

		return GetComboBox()->currentIndex() < combo->currentIndex();
	}


	QComboBox* GetComboBox() const
	{
		return m_combo;
	}


private:
	QComboBox *m_combo{};
};



/**
 * combo box showing the magnetic sites and sorting according to their index
 */
struct SitesComboBox : public QComboBox
{
	SitesComboBox(QWidget* parent = nullptr)
		: QComboBox(parent), m_item{new SitesWidgetItem{this}}
	{
	}


	virtual ~SitesComboBox()
	{
		clear();
	}


	SitesComboBox(const SitesComboBox&) = delete;
	SitesComboBox& operator=(const SitesComboBox&) = delete;


	void clear()
	{
		if(m_item)
		{
			delete m_item;
			m_item = nullptr;
		}
	}


	SitesWidgetItem* GetItem() const
	{
		//if(!m_item)
		//	m_item = new SitesWidgetItem{this};
		return m_item;
	}


private:
	SitesWidgetItem *m_item{};
};


#endif
