r"""
Victor Shoup's NTL C++ Library

SAGE provides an interface to Victor Shoup's C++ library NTL.
Features of this library include {\em incredibly fast} arithmetic with
polynomials and assymptotically fast factorization of polynomials.
"""

__doc_exclude = []  # to include everything

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.ntl.ntl_ZZ import (
                 ntl_setSeed, \
                 ntl_ZZ as ZZ,
                 randomBnd as ZZ_random,
                 randomBits as ZZ_random_bits )

from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext as ZZ_pContext

from sage.libs.ntl.ntl_ZZ_p import (
                 ntl_ZZ_p as ZZ_p,
                 set_ZZ_p_modulus as set_modulus,
                 ntl_ZZ_p_random as ZZ_p_random )

from sage.libs.ntl.ntl_ZZX import (
                 ntl_ZZX as ZZX,
                 zero_ZZX, one_ZZX )

from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX as ZZ_pX

from sage.libs.ntl.ntl_mat_ZZ import ntl_mat_ZZ as mat_ZZ

from sage.libs.ntl.ntl_GF2X import (
                 ntl_GF2X as GF2X
                  )


from sage.libs.ntl.ntl_GF2E import (
                 ntl_GF2E as GF2E, \
                 ntl_GF2E_random as GF2E_random, \
                 ntl_GF2E_modulus as GF2E_modulus, \
                 ntl_GF2E_modulus_degree as GF2E_degree, \
                 ntl_GF2E_sage as GF2E_sage,
                 GF2X_hex_repr as hex_output )

from sage.libs.ntl.ntl_GF2EX import ntl_GF2EX as GF2EX

from sage.libs.ntl.ntl_mat_GF2E import ntl_mat_GF2E as mat_GF2E
