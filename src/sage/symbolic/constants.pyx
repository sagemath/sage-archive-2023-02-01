###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "../libs/ginac/decl.pxi"

from sage.symbolic.expression cimport new_Expression_from_GEx


pi = new_Expression_from_GEx(g_Pi)
catalan = new_Expression_from_GEx(g_Catalan)
euler = new_Expression_from_GEx(g_Euler)

