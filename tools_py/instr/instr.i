/**
 * tlibs2 -- swig interface for instrument file loader
 * @author Tobias Weber <tweber@ill.fr>
 * @date 4-jun-2020
 * @license see 'LICENSE' file
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

%module instr
%{
	#include "../../tlibs2/libs/instr.h"
%}

%include "std_vector.i"
%include "std_array.i"
%include "std_map.i"
%include "std_unordered_map.i"
%include "std_pair.i"
%include "std_string.i"
%include "std_shared_ptr.i"

%shared_ptr(tl2::FileInstrBase<double>);

%template(VecStr) std::vector<std::string>;
%template(VecD) std::vector<double>;
%template(ArrD2) std::array<double, 2>;
%template(ArrD3) std::array<double, 3>;
%template(ArrD4) std::array<double, 4>;
%template(ArrD5) std::array<double, 5>;
%template(ArrB3) std::array<bool, 3>;
%template(MapStrStr) std::unordered_map<std::string, std::string>;

%include "../../tlibs2/libs/instr.h"
%include "../../tlibs2/libs/instr/base.h"
%include "../../tlibs2/libs/instr/base_loader.h"
%include "../../tlibs2/libs/instr/frm.h"
%include "../../tlibs2/libs/instr/psi.h"
%include "../../tlibs2/libs/instr/macs.h"
%include "../../tlibs2/libs/instr/tax.h"
%include "../../tlibs2/libs/instr/trisp.h"
%include "../../tlibs2/libs/instr/raw.h"

%template(FileInstrBaseD) tl2::FileInstrBase<double>;
