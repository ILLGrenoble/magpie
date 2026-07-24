/**
 * magnetic dynamics -- outputting of the hamiltonian
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2026
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * magpie & mag-core
 * Copyright (C) 2018-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// these need to be included before all other things on mingw
#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include "magdyn.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"

using namespace tl2_ops;



/**
 * shows the hamilton operator for a given Q position
 */
void MagDynDlg::CreateHamiltonPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_hamiltonianpanel = new QWidget(this);

	// hamiltonian
	m_hamiltonian = new QTextEdit(m_hamiltonianpanel);
	m_hamiltonian->setReadOnly(true);
	m_hamiltonian->setWordWrapMode(QTextOption::NoWrap);
	m_hamiltonian->setLineWrapMode(QTextEdit::NoWrap);
	m_hamiltonian->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// Q coordinates
	m_Q[0] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[1] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[2] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[0]->setToolTip("Momentum transfer component h (rlu).");
	m_Q[1]->setToolTip("Momentum transfer component k (rlu).");
	m_Q[2]->setToolTip("Momentum transfer component l (rlu).");

	for(int i = 0; i < 3; ++i)
	{
		m_Q[i]->setDecimals(4);
		m_Q[i]->setMinimum(-99.9999);
		m_Q[i]->setMaximum(+99.9999);
		m_Q[i]->setSingleStep(0.01);
		m_Q[i]->setValue(0.);
		//m_Q[i]->setSuffix(" rlu");
		m_Q[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q[i]->setPrefix(hklPrefix[i]);
	}

	// main Q index
	m_Qidx = new QSpinBox(m_hamiltonianpanel);
	m_Qidx->setMinimum(0);
	m_Qidx->setMaximum(m_num_points->value() - 1);
	m_Qidx->setValue(0);
	m_Qidx->setPrefix("Q_idx = ");
	m_Qidx->setToolTip("Q index from the dispersion panel.");
	m_Qidx->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	QPushButton *btnFromDispersion = new QPushButton(m_hamiltonianpanel);
	btnFromDispersion->setText("From Dispersion");
	btnFromDispersion->setToolTip("Set the Q position with the given index from the dispersion panel.");
	btnFromDispersion->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});

	QGridLayout *grid = new QGridLayout(m_hamiltonianpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_hamiltonian, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Q (rlu):", m_hamiltonianpanel), y, 0, 1, 1);
	grid->addWidget(m_Q[0], y, 1, 1, 1);
	grid->addWidget(m_Q[1], y, 2, 1, 1);
	grid->addWidget(m_Q[2], y++, 3, 1, 1);
	grid->addWidget(m_Qidx, y, 2, 1, 1);
	grid->addWidget(btnFromDispersion, y++, 3, 1, 1);

	// signals
	for(int i = 0; i < 3; ++i)
	{
		connect(m_Q[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			if(this->m_autocalc->isChecked())
				this->CalcHamiltonian();
		});
	}

	connect(btnFromDispersion, &QAbstractButton::clicked, [this]()
	{
		// Q range from main dispersion
		auto [Q_start, Q_end] = GetDispersionQ();

		int idx = m_Qidx->value();
		int num_pts = m_num_points->value();

		// interpolate Q
		const t_vec_real Q = num_pts > 1
			? tl2::create<t_vec_real>(
			{
				std::lerp(Q_start[0], Q_end[0], t_real(idx) / t_real(num_pts - 1)),
				std::lerp(Q_start[1], Q_end[1], t_real(idx) / t_real(num_pts - 1)),
				std::lerp(Q_start[2], Q_end[2], t_real(idx) / t_real(num_pts - 1)),
			})
			: tl2::create<t_vec_real>({ Q_start[0], Q_start[1], Q_start[2] });

		// set Q
		for(int i = 0; i < 3; ++i)
			m_Q[i]->setValue(Q[i]);
	});


	m_tabs_out->addTab(m_hamiltonianpanel, "Hamiltonian");
}



/**
 * calculate the hamiltonian for a single Q value
 */
void MagDynDlg::CalcHamiltonian()
{
	if(m_ignoreCalc)
		return;

	// options
	const bool only_energies = !m_use_weights->isChecked();
	const bool use_projector = m_use_projector->isChecked();
	const bool use_polcoords = m_use_polcoords->isChecked();;
	const bool ignore_annihilation = m_ignore_annihilation->isChecked();
	const bool unite_degeneracies = m_unite_degeneracies->isChecked();
	const bool force_incommensurate = m_force_incommensurate->isChecked();

	m_dyn.SetUniteDegenerateEnergies(unite_degeneracies);
	m_dyn.SetForceIncommensurate(force_incommensurate);

	m_hamiltonian->clear();

	const t_vec_real Q = tl2::create<t_vec_real>(
	{
		(t_real)m_Q[0]->value(),
		(t_real)m_Q[1]->value(),
		(t_real)m_Q[2]->value(),
	});

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);


	// print hamiltonian
	auto print_H = [&ostr](const t_mat& H, const t_vec_real& Qvec,
		const std::string& Qstr = "Q", const std::string& mQstr = "-Q",
		const std::string& title = "")
	{
		ostr << "<p><h3>Hamiltonian at " << Qstr <<  " = ("
			<< Qvec[0] << ", " << Qvec[1] << ", " << Qvec[2] << ")"
			<< title << "</h3>";
		ostr << "<table style=\"border:0px\">";

		// horizontal header
		ostr << "<tr><th/>";
		for(std::size_t col = 0; col < H.size2()/2; ++col)
		{
			ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (col + 1)
				<< "</sub>(" << Qstr << ")" << "</th>";
		}
		for(std::size_t col = H.size2()/2; col < H.size2(); ++col)
		{
			ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (col - H.size2()/2 + 1)
				<< "</sub><sup>&#x2020;</sup>(" << mQstr << ")" << "</th>";
		}
		ostr << "</tr>";

		for(std::size_t row = 0; row < H.size1(); ++row)
		{
			ostr << "<tr>";

			// vertical header
			if(row < H.size1() / 2)
			{
				ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (row + 1)
					<< "</sub><sup>&#x2020;</sup>(" << Qstr << ")" << "</th>";
			}
			else
			{
				ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (row - H.size1()/2 + 1)
					<< "</sub>(" << mQstr << ")" << "</th>";
			}

			// components
			for(std::size_t col = 0; col < H.size2(); ++col)
			{
				t_cplx elem = H(row, col);
				tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
				ostr << "<td style=\"padding-right:8px\">"
					<< elem << "</td>";
			}

			ostr << "</tr>";
		}
		ostr << "</table></p>";
	};


	// get hamiltonian at Q
	t_mat H = m_dyn.CalcHamiltonian(Q);
	const bool is_comm = !m_dyn.IsIncommensurate();
	if(m_hamiltonian_comp[0]->isChecked() || is_comm)  // always calculate commensurate case
		print_H(H, Q, "Q", "-Q");


	// print shifted hamiltonians for incommensurate case
	bool print_incomm_p = false;
	bool print_incomm_m = false;
	t_vec_real O;

	if(!is_comm)
	{
		// ordering wave vector
		O = tl2::create<t_vec_real>(
		{
			(t_real)m_ordering[0]->value(),
			(t_real)m_ordering[1]->value(),
			(t_real)m_ordering[2]->value(),
		});

		if(!tl2::equals_0<t_vec_real>(O, g_eps))
		{
			if(m_hamiltonian_comp[1]->isChecked())
			{
				// get hamiltonian at Q + ordering vector
				t_mat H_p = m_dyn.CalcHamiltonian(Q + O);
				print_H(H_p, Q + O, "Q + O", "Q - O");

				print_incomm_p = true;
			}

			if(m_hamiltonian_comp[2]->isChecked())
			{
				// get hamiltonian at Q - ordering vector
				t_mat H_m = m_dyn.CalcHamiltonian(Q - O);
				print_H(H_m, Q - O, "Q - O", "Q + O");

				print_incomm_m = true;
			}
		}
	}


	// get energies and correlation functions
	using t_E_and_S = typename decltype(m_dyn)::EnergyAndWeight;
	typename t_magdyn::SofQE S;

	if(is_comm)
	{
		// commensurate case
		S = m_dyn.CalcEnergiesFromHamiltonian(H, Q, only_energies);
		if(!only_energies)
			m_dyn.CalcIntensities(S);
		if(unite_degeneracies)
			S = m_dyn.UniteEnergies(S);
	}
	else
	{
		// incommensurate case
		S = m_dyn.CalcEnergies(Q, only_energies);
	}


	ostr << "<hr>";
	print_H(S.H_comm, Q, "Q", "-Q", ", Correct Commutators");
	if(print_incomm_p)
		print_H(S.H_comm_p, Q + O, "Q + O", "Q - O", ", Correct Commutators");
	if(print_incomm_m)
		print_H(S.H_comm_m, Q - O, "Q - O", "Q + O", ", Correct Commutators");
	ostr << "<hr>";


	if(only_energies)  // print energies
	{
		// split into positive and negative energies
		std::vector<t_magdyn::EnergyAndWeight> Es_neg, Es_pos;
		for(const t_E_and_S& E_and_S : S.E_and_S)
		{
			t_real E = E_and_S.E;

			if(E < t_real(0))
				Es_neg.push_back(E_and_S);
			else
				Es_pos.push_back(E_and_S);
		}

		std::stable_sort(Es_neg.begin(), Es_neg.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return std::abs(E1) < std::abs(E2);
		});

		std::stable_sort(Es_pos.begin(), Es_pos.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return std::abs(E1) < std::abs(E2);
		});

		ostr << "<p><h3>Energies</h3>";
		ostr << "<table style=\"border:0px\">";
		ostr << "<tr>";
		ostr << "<th style=\"padding-right:8px\">Creation</th>";
		for(const t_E_and_S& E_and_S : Es_pos)
		{
			t_real E = E_and_S.E;
			tl2::set_eps_0(E);

			ostr << "<td style=\"padding-right:8px\">"
				<< E << " meV" << "</td>";
		}
		ostr << "</tr>";

		if(!ignore_annihilation)
		{
			ostr << "<tr>";
			ostr << "<th style=\"padding-right:8px\">Annihilation</th>";
			for(const t_E_and_S& E_and_S : Es_neg)
			{
				t_real E = E_and_S.E;
				tl2::set_eps_0(E);

				ostr << "<td style=\"padding-right:8px\">"
					<< E << " meV" << "</td>";
			}
			ostr << "</tr>";
		}

		ostr << "</table></p>";
	}
	else  // print energies and weights
	{
		std::stable_sort(S.E_and_S.begin(), S.E_and_S.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return E1 < E2;
		});

		ostr << "<p><h3>Spectrum</h3>";
		ostr << "<table style=\"border:0px\">";
		ostr << "<tr>";
		ostr << "<th style=\"padding-right:16px\">Energy E</td>";
		ostr << "<th style=\"padding-right:16px\">Correlation S(Q, E)</td>";
		ostr << "<th style=\"padding-right:16px\">Projection S<sub>&#x27C2;</sub>(Q, E)</td>";
		ostr << "<th style=\"padding-right:16px\">Weight</td>";
		ostr << "</tr>";

		for(const t_E_and_S& E_and_S : S.E_and_S)
		{
			t_real E = E_and_S.E;
			if(ignore_annihilation && E < t_real(0))
				continue;

			const t_mat& S = use_polcoords ? E_and_S.S_pol : E_and_S.S;
			const t_mat& S_perp = use_polcoords ? E_and_S.S_pol_perp : E_and_S.S_perp;
			t_real weight = use_polcoords ? E_and_S.weight_pol_perp : E_and_S.weight_perp;

			if(!use_projector)
				weight = tl2::trace<t_mat>(S).real();

			tl2::set_eps_0(E);
			tl2::set_eps_0(weight);

			// E
			ostr << "<tr>";
			ostr << "<td style=\"padding-right:16px\">"
				<< E << " meV" << "</td>";

			// S(Q, E)
			ostr << "<td style=\"padding-right:16px\">";
			ostr << "<table style=\"border:0px\">";
			for(std::size_t i = 0; i < S.size1(); ++i)
			{
				ostr << "<tr>";
				for(std::size_t j = 0; j < S.size2(); ++j)
				{
					t_cplx elem = S(i, j);
					tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
					ostr << "<td style=\"padding-right:8px\">"
						<< elem << "</td>";
				}
				ostr << "</tr>";
			}
			ostr << "</table>";
			ostr << "</td>";

			// S_perp(Q, E)
			ostr << "<td style=\"padding-right:16px\">";
			ostr << "<table style=\"border:0px\">";
			for(std::size_t i = 0; i < S_perp.size1(); ++i)
			{
				ostr << "<tr>";
				for(std::size_t j = 0; j < S_perp.size2(); ++j)
				{
					t_cplx elem = S_perp(i, j);
					tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
					ostr << "<td style=\"padding-right:8px\">"
						<< elem << "</td>";
				}
				ostr << "</tr>";
			}
			ostr << "</table>";
			ostr << "</td>";

			// tr(S_perp(Q, E))
			ostr << "<td style=\"padding-right:16px\">" << weight << "</td>";

			ostr << "</tr>";
		}
		ostr << "</table></p>";
	}


	// print eigenstates
	if(S.E_and_S.size() && S.E_and_S[0].state.size())
	{
		ostr << "<hr>";

		ostr << "<p><h3>Eigenstates</h3>";
		ostr << "<table style=\"border:0px\">";
		ostr << "<tr>";
		ostr << "<th style=\"padding-right:16px\">Energy E</td>";
		ostr << "<th style=\"padding-right:16px\">Degeneracy</td>";
		ostr << "<th style=\"padding-right:16px\">State |s></td>";
		ostr << "</tr>";

		for(const t_E_and_S& E_and_S : S.E_and_S)
		{
			t_size degen = E_and_S.degeneracy;
			t_real E = E_and_S.E;
			if(ignore_annihilation && E < t_real(0))
				continue;

			t_vec state = E_and_S.state;

			tl2::set_eps_0(E);
			tl2::set_eps_0(state);

			// energy
			ostr << "<tr>";
			ostr << "<td style=\"padding-right:16px\">"
				<< E << " meV" << "</td>";

			// degeneracy
			ostr << "<td style=\"padding-right:16px\">"
				<< degen << "</td>";

			// state
			ostr << "<td style=\"padding-right:16px\">";
			for(t_size idx = 0; idx < state.size(); ++idx)
			{
				ostr << state[idx];
				if(idx < state.size() - 1)
					ostr << ", ";
			}
			ostr << "</td>";
			ostr << "</tr>";
		}

		ostr << "</table></p>";
	}


	m_hamiltonian->setHtml(ostr.str().c_str());
}
