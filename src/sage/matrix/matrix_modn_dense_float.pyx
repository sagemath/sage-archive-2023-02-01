"""
Dense matrices over `\ZZ/n\ZZ` for `n < 2^{11}` using LinBox's ``Modular<float>``

AUTHORS:
- Burcin Erocal
- Martin Albrecht
"""
###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "../ext/stdsage.pxi"
include "../ext/interrupt.pxi"

# randstate in template needs this
include '../ext/random.pxi'

from sage.libs.linbox.modular cimport ModFloatField as ModField, ModFloatFieldElement as ModFieldElement

from sage.libs.linbox.fflas cimport ModFloat_fgemm as Mod_fgemm, ModFloat_fgemv as Mod_fgemv, \
        ModFloatDet as ModDet, \
        ModFloatRank as ModRank, ModFloat_echelon as Mod_echelon, \
        ModFloat_applyp as Mod_applyp, \
        ModFloat_MinPoly as Mod_MinPoly, \
        ModFloat_CharPoly as Mod_CharPoly

# LinBox supports up to 2^11 using float but that's double dog slow,
# so we pick a smaller value for crossover
MAX_MODULUS = 2**8

include "matrix_modn_dense_template.pxi"

cdef class Matrix_modn_dense_float(Matrix_modn_dense_template):
    """
    Dense matrices over `\ZZ/n\ZZ` for `n < 2^{11}` using LinBox's ``Modular<float>``
    """
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        """
        Set the (i,j) entry of self to the int value.
        """
        self._matrix[i][j] = <float>value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        # doesn't do bounds or any other checks; assumes x is in self._base_ring
        self._matrix[i][j] = <float>(<IntegerMod_int>x).ivalue

    cdef IntegerMod_int get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        # doesn't do checks
        return IntegerMod_int(self._base_ring, <mod_int>(<Matrix_modn_dense_template>self)._matrix[i][j])
