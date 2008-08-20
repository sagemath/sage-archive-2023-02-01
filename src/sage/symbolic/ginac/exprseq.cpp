/** @file exprseq.cpp
 *
 *  Implementation of GiNaC's exprseq. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "exprseq.h"

namespace GiNaC {

template <> GINAC_IMPLEMENT_REGISTERED_CLASS_OPT_T(exprseq, basic,
  print_func<print_context>(&exprseq::do_print).
  print_func<print_tree>(&exprseq::do_print_tree))

/** Specialization of container::info() for exprseq. */
template <> bool exprseq::info(unsigned inf) const
{
	if (inf == info_flags::exprseq)
		return true;
	else
		return inherited::info(inf);
}

} // namespace GiNaC
