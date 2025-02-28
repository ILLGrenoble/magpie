/**
 * settings dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date apr-2021, nov-2023
 * @note The present version was forked on 28-Nov-2023 from TAS-Paths (https://github.com/ILLGrenoble/taspaths).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * TAS-Paths
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

#ifndef __MAGDYN_SETTINGS__
#define __MAGDYN_SETTINGS__

#include <QtCore/QSettings>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QStyleFactory>
#include <QtWidgets/QFontDialog>
#include <QtWidgets/QHeaderView>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/qt/numerictablewidgetitem.h"

#include <type_traits>
#include <string>
#include <variant>
#include <unordered_map>

#ifndef MAGDYN_SETTINGS_USE_QT_SIGNALS
	// qt signals can't be emitted from a template class
	// TODO: remove this as soon as this is supported
	#include <boost/signals2/signal.hpp>
#endif

#include "defs.h"



// ----------------------------------------------------------------------------
// settings variable struct
// ----------------------------------------------------------------------------
enum class SettingsVariableEditor
{
	NONE,
	YESNO,
	COMBOBOX,
};



struct SettingsVariable
{
	using t_variant = std::variant<t_real, int, unsigned int, std::string>;
	using t_variant_ptr = std::variant<t_real*, int*, unsigned int*, std::string*>;

	const char* description{};
	const char* key{};

	t_variant_ptr value{};
	bool is_angle{false};

	SettingsVariableEditor editor{SettingsVariableEditor::NONE};
	const char* editor_config{};
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
enum class SettingsColumn : int
{
	SETTING = 0,
	TYPE    = 1,
	VALUE   = 2,
};



/**
 * settings dialog
 */
template<std::size_t num_settingsvariables,
	const std::array<SettingsVariable,
		num_settingsvariables> *settingsvariables>
class SettingsDlg : public QDialog
{
#ifdef MAGDYN_SETTINGS_USE_QT_SIGNALS
	Q_OBJECT
#endif

public:
	/**
	 * constructor
	 */
	SettingsDlg(QWidget* parent = nullptr, QSettings *sett = nullptr)
		: QDialog{parent}, m_sett{sett}
	{
		InitGui();
	}



	/**
	 * set-up the settings dialog gui
	 */
	void InitGui()
	{
		setWindowTitle("Preferences");
		setSizeGripEnabled(true);

		// table column widths
		int col_sett_w = 200;
		int col_type_w = 100;
		int col_value_w = 150;

		if(m_sett)
		{
			// restore dialog geometry
			if(m_sett->contains("settings/geo"))
				restoreGeometry(m_sett->value("settings/geo").toByteArray());
			else
				resize(512, 425);

			// restore table column widths
			if(m_sett->contains("settings/col_setting_width"))
				col_sett_w = m_sett->value("settings/col_setting_width").toInt();
			if(m_sett->contains("settings/col_type_width"))
				col_type_w = m_sett->value("settings/col_type_width").toInt();
			if(m_sett->contains("settings/col_value_width"))
				col_value_w = m_sett->value("settings/col_value_width").toInt();
		}

		// general settings
		QWidget *panelGeneral = new QWidget(this);
		QGridLayout *gridGeneral = new QGridLayout(panelGeneral);
		gridGeneral->setSpacing(4);
		gridGeneral->setContentsMargins(6, 6, 6, 6);

		// create the settings table
		m_table = new QTableWidget(panelGeneral);
		m_table->setShowGrid(true);
		m_table->setAlternatingRowColors(true);
		m_table->setSortingEnabled(false);
		m_table->setMouseTracking(false);
		m_table->setSelectionBehavior(QTableWidget::SelectRows);
		m_table->setSelectionMode(QTableWidget::SingleSelection);

		// table headers
		m_table->horizontalHeader()->setDefaultSectionSize(125);
		m_table->verticalHeader()->setDefaultSectionSize(32);
		m_table->verticalHeader()->setVisible(false);
		m_table->setColumnCount(3);
		m_table->setColumnWidth((int)SettingsColumn::SETTING, col_sett_w);
		m_table->setColumnWidth((int)SettingsColumn::TYPE, col_type_w);
		m_table->setColumnWidth((int)SettingsColumn::VALUE, col_value_w);
		m_table->setHorizontalHeaderItem((int)SettingsColumn::SETTING, new QTableWidgetItem{"Setting"});
		m_table->setHorizontalHeaderItem((int)SettingsColumn::TYPE, new QTableWidgetItem{"Type"});
		m_table->setHorizontalHeaderItem((int)SettingsColumn::VALUE, new QTableWidgetItem{"Value"});

		// table contents
		PopulateSettingsTable();


		// search field
		QLabel *labelSearch = new QLabel("Search:", panelGeneral);
		QLineEdit *editSearch = new QLineEdit(panelGeneral);

		labelSearch->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
		editSearch->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

		gridGeneral->addWidget(m_table, 0, 0, 1, 2);
		gridGeneral->addWidget(labelSearch, 1, 0, 1, 1);
		gridGeneral->addWidget(editSearch, 1, 1, 1, 1);


		// gui settings
		QWidget *panelGui = new QWidget(this);
		QGridLayout *gridGui = new QGridLayout(panelGui);
		gridGui->setSpacing(4);
		gridGui->setContentsMargins(6, 6, 6, 6);
		int yGui = 0;

		// theme
		if(s_theme)
		{
			QLabel *labelTheme = new QLabel("Style:", panelGui);
			labelTheme->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
			m_comboTheme = new QComboBox(panelGui);
			m_comboTheme->addItems(QStyleFactory::keys());

			get_setting<QString>(m_sett, "theme", s_theme);
			if(*s_theme != "")
			{
				int idxTheme = m_comboTheme->findText(*s_theme);
				if(idxTheme >= 0 && idxTheme < m_comboTheme->count())
					m_comboTheme->setCurrentIndex(idxTheme);
			}

			gridGui->addWidget(labelTheme, yGui,0,1,1);
			gridGui->addWidget(m_comboTheme, yGui++,1,1,2);
		}

		// font
		if(s_font)
		{
			QLabel *labelFont = new QLabel("Font:", panelGui);
			labelFont->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

			m_editFont = new QLineEdit(panelGui);
			m_editFont->setReadOnly(true);

			QPushButton *btnFont = new QPushButton("Select...", panelGui);
			connect(btnFont, &QPushButton::clicked, [this]()
			{
				// current font
				QFont font = QApplication::font();

				// select a new font
				bool okClicked = false;
				font = QFontDialog::getFont(&okClicked, font, this);
				if(okClicked)
				{
					*s_font = font.toString();
					if(*s_font == "")
						*s_font = QApplication::font().toString();
					m_editFont->setText(*s_font);
				}

				// hack for the QFontDialog hiding the settings dialog
				this->show();
				this->raise();
				this->activateWindow();
			});

			get_setting<QString>(m_sett, "font", s_font);
			if(*s_font == "")
				*s_font = QApplication::font().toString();
			m_editFont->setText(*s_font);

			gridGui->addWidget(labelFont, yGui,0,1,1);
			gridGui->addWidget(m_editFont, yGui,1,1,1);
			gridGui->addWidget(btnFont, yGui++,2,1,1);
		}

		// native menubar
		if(s_use_native_menubar)
		{
			m_checkMenubar = new QCheckBox("Use native menubar.", panelGui);
			get_setting<int>(m_sett, "native_menubar", s_use_native_menubar);
			m_checkMenubar->setChecked(*s_use_native_menubar != 0);

			gridGui->addWidget(m_checkMenubar, yGui++,0,1,3);
		}

		// native dialogs
		if(s_use_native_dialogs)
		{
			m_checkDialogs = new QCheckBox("Use native dialogs.", panelGui);
			get_setting<int>(m_sett, "native_dialogs", s_use_native_dialogs);
			m_checkDialogs->setChecked(*s_use_native_dialogs != 0);

			gridGui->addWidget(m_checkDialogs, yGui++,0,1,3);
		}

		QSpacerItem *spacer_end = new QSpacerItem(1, 1,
			QSizePolicy::Minimum, QSizePolicy::Expanding);
		gridGui->addItem(spacer_end, yGui++,0,1,3);


		// main grid
		QGridLayout *grid = new QGridLayout(this);
		grid->setSpacing(4);
		grid->setContentsMargins(8, 8, 8, 8);
		int y = 0;

		QTabWidget *tab = new QTabWidget(this);
		tab->addTab(panelGeneral, "General");
		tab->addTab(panelGui, "GUI");
		grid->addWidget(tab, y++,0,1,1);

		QLabel *labelRestart = new QLabel(
			"Important: Applying all settings requires a program restart.",
			this);
		labelRestart->setWordWrap(true);
		QFont fontRestart = labelRestart->font();
		fontRestart.setBold(true);
		labelRestart->setFont(fontRestart);
		grid->addWidget(labelRestart, y++,0,1,1);

		QDialogButtonBox *buttons = new QDialogButtonBox(this);
		buttons->setStandardButtons(
			QDialogButtonBox::Ok |
			QDialogButtonBox::Apply |
			QDialogButtonBox::RestoreDefaults |
			QDialogButtonBox::Cancel);
		grid->addWidget(buttons, y++,0,1,1);


		// connections
		connect(buttons, &QDialogButtonBox::accepted, this, &SettingsDlg::accept);
		connect(buttons, &QDialogButtonBox::rejected, this, &SettingsDlg::reject);
		connect(buttons, &QDialogButtonBox::clicked, [this, buttons](QAbstractButton* btn)
		{
			// apply button was pressed
			if(btn == static_cast<QAbstractButton*>(buttons->button(QDialogButtonBox::Apply)))
				ApplySettings();
			else if(btn == static_cast<QAbstractButton*>(buttons->button(QDialogButtonBox::RestoreDefaults)))
				RestoreDefaultSettings();
		});

		// search items
		connect(editSearch, &QLineEdit::textChanged, [editSearch, this]()
		{
			QList<QTableWidgetItem*> items =
				m_table->findItems(editSearch->text(), Qt::MatchContains);

			// scroll to first found item
			if(items.size())
				m_table->setCurrentItem(items[0]);
		});
	}



	/**
	 * destructor
	 */
	virtual ~SettingsDlg()
	{
	}



	/**
	 * copy constructor
	 */
	SettingsDlg(const SettingsDlg&) = delete;
	const SettingsDlg& operator=(const SettingsDlg&) = delete;



	/**
	 * read the settings and set the global variables
	 */
	static void ReadSettings(QSettings* sett)
	{
		// save the initial values as default settings
		static bool first_run = true;
		if(first_run)
		{
			SaveDefaultSettings();
			first_run = false;
		}

		if(!sett)
			return;

		auto seq = std::make_index_sequence<num_settingsvariables>();
		get_settings_loop(seq, sett);

		get_setting<QString>(sett, "theme", s_theme);
		get_setting<QString>(sett, "font", s_font);
		get_setting<int>(sett, "native_menubar", s_use_native_menubar);
		get_setting<int>(sett, "native_dialogs", s_use_native_dialogs);

		ApplyGuiSettings();
	}



	// common gui settings
	static void SetGuiTheme(QString* str) { s_theme = str; }
	static void SetGuiFont(QString* str) { s_font = str; }
	static void SetGuiUseNativeMenubar(int *i) { s_use_native_menubar = i; }
	static void SetGuiUseNativeDialogs(int *i) { s_use_native_dialogs = i; }



	/**
	 * save the current setting values as default values
	 */
	static void SaveDefaultSettings()
	{
		// general program settings
		auto seq = std::make_index_sequence<num_settingsvariables>();
		save_default_values_loop(seq, s_defaults);

		// gui settings
		if(s_theme)
			s_defaults.insert_or_assign("<theme>", s_theme->toStdString());
		if(s_font)
			s_defaults.insert_or_assign("<font>", s_font->toStdString());
		if(s_use_native_menubar)
			s_defaults.insert_or_assign("<native_menubar>", *s_use_native_menubar);
		if(s_use_native_dialogs)
			s_defaults.insert_or_assign("<native_dialogs>", *s_use_native_dialogs);
	}



protected:
	/**
	 * 'OK' was clicked
	 */
	virtual void accept() override
	{
		ApplySettings();

		if(m_sett)
		{
			// save dialog geometry
			m_sett->setValue("settings/geo", saveGeometry());

			// save table column widths
			if(m_table)
			{
				m_sett->setValue("settings/col_setting_width",
					m_table->columnWidth((int)SettingsColumn::SETTING));
				m_sett->setValue("settings/col_type_width",
					m_table->columnWidth((int)SettingsColumn::TYPE));
				m_sett->setValue("settings/col_value_width",
					m_table->columnWidth((int)SettingsColumn::VALUE));
			}
		}

		QDialog::accept();
	}



	/**
	 * populate the settings table using the global settings items
	 */
	void PopulateSettingsTable()
	{
		m_table->clearContents();
		m_table->setRowCount((int)num_settingsvariables);

		auto seq = std::make_index_sequence<num_settingsvariables>();
		add_table_item_loop(seq, m_table);

		// set value field editable
		for(int row = 0; row < m_table->rowCount(); ++row)
		{
			// setting key
			m_table->item(row, (int)SettingsColumn::SETTING)->setFlags(
				m_table->item(row, (int)SettingsColumn::SETTING)->flags() & ~Qt::ItemIsEditable);

			// type
			m_table->item(row, (int)SettingsColumn::TYPE)->setFlags(
				m_table->item(row, (int)SettingsColumn::TYPE)->flags() & ~Qt::ItemIsEditable);

			// value
			QTableWidgetItem *itemVal = m_table->item(row, (int)SettingsColumn::VALUE);
			if(itemVal)
			{
				m_table->item(row, (int)SettingsColumn::VALUE)->setFlags(
					itemVal->flags() | Qt::ItemIsEditable);
			}
		}
	}



	/**
	 * 'Restore Defaults' was clicked, restore original settings
	 */
	void RestoreDefaultSettings()
	{
		// general program settings
		auto seq = std::make_index_sequence<num_settingsvariables>();
		restore_default_values_loop(seq, s_defaults);

		// re-populate the settings table
		PopulateSettingsTable();

		// gui settings
		if(s_theme)
		{
			if(auto iter = s_defaults.find(std::string{"<theme>"}); iter != s_defaults.end())
			{
				if(std::holds_alternative<std::string>(iter->second))
					*s_theme = std::get<std::string>(iter->second).c_str();
				if(m_comboTheme && *s_theme != "")
				{
					int idxTheme = m_comboTheme->findText(*s_theme);
					if(idxTheme >= 0 && idxTheme < m_comboTheme->count())
						m_comboTheme->setCurrentIndex(idxTheme);
				}
			}
		}
		if(s_font)
		{
			if(auto iter = s_defaults.find(std::string{"<font>"}); iter != s_defaults.end())
			{
				if(std::holds_alternative<std::string>(iter->second))
					*s_font = std::get<std::string>(iter->second).c_str();
				if(m_editFont && *s_font != "")
					m_editFont->setText(*s_font);
			}
		}

		int* vars[] =
		{
			s_use_native_menubar,
			s_use_native_dialogs,
		};

		std::string idents[] =
		{
			"<native_menubar>",
			"<native_dialogs>",
		};

		QCheckBox* checks[] =
		{
			m_checkMenubar,
			m_checkDialogs,
		};

		for(std::size_t idx=0; idx<sizeof(vars)/sizeof(*vars); ++idx)
		{
			if(!vars[idx])
				continue;

			if(auto iter = s_defaults.find(idents[idx]);
			   iter != s_defaults.end())
			{
				if(std::holds_alternative<int>(iter->second))
					*vars[idx] = std::get<int>(iter->second);
				if(checks[idx])
					checks[idx]->setChecked(*vars[idx] != 0);
			}
		}
	}



	/**
	 * 'Apply' was clicked, write the settings from the global variables
	 */
	void ApplySettings()
	{
		auto seq = std::make_index_sequence<num_settingsvariables>();
		apply_settings_loop(seq, m_table, m_sett);

		// set the global variables
		if(s_theme)
			*s_theme = m_comboTheme->currentText();
		if(s_font)
			*s_font = m_editFont->text();
		if(s_use_native_menubar)
			*s_use_native_menubar = m_checkMenubar->isChecked();
		if(s_use_native_dialogs)
			*s_use_native_dialogs = m_checkDialogs->isChecked();

		// write out the settings
		if(m_sett)
		{
			if(s_theme)
				m_sett->setValue("theme", *s_theme);
			if(s_font)
				m_sett->setValue("font", *s_font);
			if(s_use_native_menubar)
				m_sett->setValue("native_menubar",
					*s_use_native_menubar);
			if(s_use_native_dialogs)
				m_sett->setValue("native_dialogs",
					*s_use_native_dialogs);
		}

		ApplyGuiSettings();

#ifdef MAGDYN_SETTINGS_USE_QT_SIGNALS
		emit SettingsHaveChanged();
#else
		m_sigSettingsHaveChanged();
#endif
	}



	static void ApplyGuiSettings()
	{
		// set gui theme
		if(s_theme && *s_theme != "")
		{
			if(QStyle* theme = QStyleFactory::create(*s_theme); theme)
				QApplication::setStyle(theme);
		}

		// set gui font
		if(s_font && *s_font != "")
		{
			QFont font;
			if(font.fromString(*s_font))
				QApplication::setFont(font);
		}

		// set native menubar
		if(s_use_native_menubar)
		{
			QApplication::setAttribute(
				Qt::AA_DontUseNativeMenuBar, !*s_use_native_menubar);
		}

		// set native dialogs
		if(s_use_native_dialogs)
		{
			QApplication::setAttribute(
				Qt::AA_DontUseNativeDialogs, !*s_use_native_dialogs);
		}
	}



	// ------------------------------------------------------------------------
	// helpers
	// ------------------------------------------------------------------------
	/**
	 * get type names
	 */
	template<class T>
	static constexpr const char* get_type_str()
	{
		if constexpr (std::is_same_v<T, t_real>)
			return "Real";
		else if constexpr (std::is_same_v<T, int>)
			return "Integer";
		else if constexpr (std::is_same_v<T, unsigned int>)
			return "Integer, unsigned";
		else if constexpr (std::is_same_v<T, std::string>)
			return "String";
		else
			return "Unknown";
	}



	/**
	 * adds a settings item from a global variable to the table
	 */
	template<std::size_t idx>
	static void add_table_item(QTableWidget *table)
	{
		constexpr const SettingsVariable& var = std::get<idx>(*settingsvariables);
		constexpr const auto* value = std::get<var.value.index()>(var.value);
		using t_value = std::decay_t<decltype(*value)>;

		table->setItem((int)idx, (int)SettingsColumn::SETTING,
			new QTableWidgetItem{var.description});
		table->setItem((int)idx, (int)SettingsColumn::TYPE,
			new QTableWidgetItem{get_type_str<t_value>()});

		if constexpr (!std::is_same_v<t_value, std::string>)
		{  // non-string items
			t_value finalval = *value;

			if(var.is_angle)
				finalval = tl2::r2d<t_real>(finalval);

			if(var.editor == SettingsVariableEditor::YESNO)
			{
				QComboBox *combo = new QComboBox(table);
				combo->addItem("No");
				combo->addItem("Yes");

				combo->setCurrentIndex(finalval == 0 ? 0 : 1);
				table->setCellWidget((int)idx, (int)SettingsColumn::VALUE, combo);
			}
			else if(var.editor == SettingsVariableEditor::COMBOBOX)
			{
				std::vector<std::string> config_tokens;
				tl2::get_tokens_seq<std::string, std::string>(
					var.editor_config, ";;", config_tokens, true);

				QComboBox *combo = new QComboBox(table);
				for(const std::string& config_token : config_tokens)
					combo->addItem(config_token.c_str());

				combo->setCurrentIndex((int)finalval);
				table->setCellWidget((int)idx, (int)SettingsColumn::VALUE, combo);
			}
			else
			{
				QTableWidgetItem *item = new tl2::NumericTableWidgetItem<t_value>(finalval, 10);
				table->setItem((int)idx, (int)SettingsColumn::VALUE, item);
			}
		}
		else  // string items
		{
			QTableWidgetItem *item = new QTableWidgetItem{value->c_str()};
			table->setItem((int)idx, (int)SettingsColumn::VALUE, item);
		}
	}



	/**
	 * adds all settings items from the global variables to the table
	 */
	template<std::size_t ...seq>
	static constexpr void add_table_item_loop(
		const std::index_sequence<seq...>&, QTableWidget *table)
	{
		// a sequence of function calls
		( (add_table_item<seq>(table)), ... );
	}



	/**
	 * get a value from a QSettings object
	 */
	template<class T>
	static void get_setting(QSettings* sett, const char* key, T* val)
	{
		if(sett->contains(key))
		{
			if constexpr (!std::is_same_v<T, std::string>)
				*val = sett->value(key, *val).template value<T>();
			else
				*val = sett->value(key, val->c_str()).template value<QString>().toStdString();
		}
	}



	/**
	 * gets a settings item from the QSettings object and saves it
	 * to the corresponding global variable
	 */
	template<std::size_t idx>
	static void get_settings_item(QSettings *sett)
	{
		constexpr const SettingsVariable& var = std::get<idx>(*settingsvariables);
		constexpr auto* value = std::get<var.value.index()>(var.value);
		using t_value = std::decay_t<decltype(*value)>;

		get_setting<t_value>(sett, var.key, value);
	}



	/**
	 * gets all settings items from the QSettings object and saves them
	 * to the global variables
	 */
	template<std::size_t ...seq>
	static constexpr void get_settings_loop(
		const std::index_sequence<seq...>&, QSettings *sett)
	{
		// a sequence of function calls
		( (get_settings_item<seq>(sett)), ... );
	}



	/**
	 * save the current settings value as its default value
	 */
	template<std::size_t idx>
	static void save_default_value(
		std::unordered_map<std::string, SettingsVariable::t_variant>& map)
	{
		constexpr const SettingsVariable& var = std::get<idx>(*settingsvariables);
		constexpr auto* value = std::get<var.value.index()>(var.value);

		map.insert_or_assign(var.key, *value);
	}



	/**
	 * save the current settings values as their default values
	 */
	template<std::size_t ...seq>
	static constexpr void save_default_values_loop(
		const std::index_sequence<seq...>&,
		std::unordered_map<std::string, SettingsVariable::t_variant>& map)
	{
		// a sequence of function calls
		( (save_default_value<seq>(map)), ... );
	}



	/**
	 * restore the current settings value from its default value
	 */
	template<std::size_t idx>
	static void restore_default_value(
		const std::unordered_map<std::string, SettingsVariable::t_variant>& map)
	{
		constexpr const SettingsVariable& var = std::get<idx>(*settingsvariables);
		constexpr auto* value = std::get<var.value.index()>(var.value);
		using t_value = std::decay_t<decltype(*value)>;

		// does the defaults map have this key?
		if(auto iter = map.find(std::string{var.key}); iter!=map.end())
		{
			// does the settings variable variant hold this type?
			if(std::holds_alternative<t_value>(iter->second))
				*value = std::get<t_value>(iter->second);
		}
	}



	/**
	 * restore the current settings values from their default values
	 */
	template<std::size_t ...seq>
	static constexpr void restore_default_values_loop(
		const std::index_sequence<seq...>&,
		const std::unordered_map<std::string, SettingsVariable::t_variant>& map)
	{
		// a sequence of function calls
		( (restore_default_value<seq>(map)), ... );
	}



	/**
	 * reads a settings item from the table and saves it to the
	 * corresponding global variable and to the QSettings object
	 */
	template<std::size_t idx>
	static void apply_settings_item(QTableWidget *table, QSettings *sett)
	{
		constexpr const SettingsVariable& var = std::get<idx>(*settingsvariables);
		constexpr auto* value = std::get<var.value.index()>(var.value);
		using t_value = std::decay_t<decltype(*value)>;

		t_value finalval{};

		if constexpr (!std::is_same_v<t_value, std::string>)
		{  // non-string values
			// use the value from the editor widget if available
			if(var.editor == SettingsVariableEditor::YESNO ||
				var.editor == SettingsVariableEditor::COMBOBOX)
			{
				QComboBox *combo = static_cast<QComboBox*>(
					table->cellWidget((int)idx, (int)SettingsColumn::VALUE));

				finalval = (t_value)combo->currentIndex();
			}
			else
			{
				finalval = dynamic_cast<tl2::NumericTableWidgetItem<t_value>*>(
					table->item((int)idx,
					(int)SettingsColumn::VALUE))->GetValue();

				if(var.is_angle)
					finalval = tl2::d2r<t_real>(finalval);
			}
		}
		else  // string values
		{
			finalval = dynamic_cast<QTableWidgetItem*>(
				table->item((int)idx,
				(int)SettingsColumn::VALUE))->text().toStdString();
		}

		// set the global variable
		*value = finalval;

		// write out the settings
		if(sett)
		{
			if constexpr (!std::is_same_v<t_value, std::string>)
				sett->setValue(var.key, *value);
			else
				sett->setValue(var.key, value->c_str());				
		}
	}



	/**
	 * reads all settings items from the table and saves them to the
	 * global variables and to the QSettings object
	 */
	template<std::size_t ...seq>
	static constexpr void apply_settings_loop(
		const std::index_sequence<seq...>&, QTableWidget* table, QSettings *sett)
	{
		// a sequence of function calls
		( (apply_settings_item<seq>(table, sett)), ... );
	}
	// ------------------------------------------------------------------------



private:
	QSettings *m_sett{nullptr};
	QTableWidget *m_table{nullptr};

	QComboBox *m_comboTheme{nullptr};
	QLineEdit *m_editFont{nullptr};
	QCheckBox *m_checkMenubar{nullptr};
	QCheckBox *m_checkDialogs{nullptr};

	// common gui settings
	static QString *s_theme;           // gui theme
	static QString *s_font;            // gui font
	static int *s_use_native_menubar;  // use native menubar?
	static int *s_use_native_dialogs;  // use native dialogs?

	// default setting values
	static std::unordered_map<std::string, SettingsVariable::t_variant> s_defaults;



#ifdef MAGDYN_SETTINGS_USE_QT_SIGNALS
signals:
	// signal emitted when settings are applied
	void SettingsHaveChanged();

#else
	boost::signals2::signal<void()> m_sigSettingsHaveChanged{};

public:
	template<class t_slot>
	boost::signals2::connection AddChangedSettingsSlot(const t_slot& slot)
	{
		return m_sigSettingsHaveChanged.connect(slot);
	}
#endif

};



// initialise static variable holding defaults
template<
	std::size_t num_settingsvariables,
	const std::array<SettingsVariable, num_settingsvariables> *settingsvariables>
std::unordered_map<std::string, SettingsVariable::t_variant>
SettingsDlg<num_settingsvariables, settingsvariables>::s_defaults{};


// initialise static variables for common gui settings
template<
	std::size_t num_settingsvariables,
	const std::array<SettingsVariable, num_settingsvariables> *settingsvariables>
QString *SettingsDlg<num_settingsvariables, settingsvariables>::s_theme{nullptr};

template<
	std::size_t num_settingsvariables,
	const std::array<SettingsVariable, num_settingsvariables> *settingsvariables>
QString *SettingsDlg<num_settingsvariables, settingsvariables>::s_font{nullptr};

template<
	std::size_t num_settingsvariables,
	const std::array<SettingsVariable, num_settingsvariables> *settingsvariables>
int *SettingsDlg<num_settingsvariables, settingsvariables>::s_use_native_menubar{nullptr};

template<
	std::size_t num_settingsvariables,
	const std::array<SettingsVariable, num_settingsvariables> *settingsvariables>
int *SettingsDlg<num_settingsvariables, settingsvariables>::s_use_native_dialogs{nullptr};

// ----------------------------------------------------------------------------

#endif
