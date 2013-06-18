/* ============================================================================

   recurrences_ntl.h:  header for recurrences_ntl.cpp

   This file is part of hypellfrob (version 2.1.1).

   Copyright (C) 2007, 2008, David Harvey

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

============================================================================ */


#include <NTL/ZZ.h>
#include <vector>


namespace hypellfrob {


template <typename SCALAR, typename POLY, typename POLYMODULUS,
          typename VECTOR, typename MATRIX, typename FFTREP>
void ntl_interval_products(std::vector<MATRIX>& output,
                           const MATRIX& M0, const MATRIX& M1,
                           const std::vector<NTL::ZZ>& target);


};

// ----------------------- end of file
