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

from sage.libs.ntl.ntl import (make_new_ZZ as ZZ,
                 ntl_ZZ as ZZ_class,
                               randomBnd as ZZ_random,
                               randomBits as ZZ_random_bits,
                 make_new_ZZ_p as ZZ_p,
                 ntl_ZZ_p as ZZ_p_class, \
                 set_ZZ_p_modulus as set_modulus, \
                 ntl_ZZ_p_random as ZZ_p_random,

                 make_new_ZZX as ZZX, \
                 ntl_ZZX as ZZX_class, \
                 zero_ZZX, one_ZZX, \

                 ntl_setSeed, \

                 make_new_ZZ_pX as ZZ_pX, \
                 ntl_ZZ_pX as ZZ_pX_class, \

                 ntl_mat_ZZ as mat_ZZ, \

                 make_new_GF2X as GF2X, \
                 ntl_GF2X as GF2X_class, \

                 make_new_GF2E as GF2E, \
                 ntl_GF2E as GF2E_class, \
                 ntl_GF2E_random as GF2E_random, \
                 ntl_GF2E_modulus as GF2E_modulus, \
                 ntl_GF2E_modulus_degree as GF2E_degree, \
                 ntl_GF2E_sage as GF2E_sage, \
                 GF2X_hex_repr as hex_output, \

                 make_new_GF2EX as GF2EX, \
                 ntl_GF2EX as GF2EX_class, \


                 ntl_mat_GF2E as mat_GF2E, \

                 )

mat_ZZ_class = mat_ZZ
mat_GF2E_class = mat_GF2E
