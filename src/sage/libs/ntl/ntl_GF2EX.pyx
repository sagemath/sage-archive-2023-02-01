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

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

##############################################################################
#
# ntl_GF2EX: Polynomials over GF(2) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de> 2006-01: initial version
#
##############################################################################

cdef class ntl_GF2EX:
    r"""
    """
    def __init__(self, x=[]):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(3))
            sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: ntl.GF2EX('[[1 0] [2 1]]')
            [[1] [0 1]]
        """
        s = str(x)
        _sig_on
        GF2EX_from_str(&self.x, s)
        _sig_off

    def __new__(self, x=[]):
        GF2EX_construct(&self.x)

    def __dealloc__(self):
        GF2EX_destruct(&self.x)

    def __reduce__(self):
        raise NotImplementedError

    def __repr__(self):
        return GF2EX_to_PyString(&self.x)

    def __mul__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        cdef ntl_GF2EX r = ntl_GF2EX()
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        GF2EX_mul(r.x, self.x, y.x)
        _sig_off
        return r

    def __sub__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        cdef ntl_GF2EX r = ntl_GF2EX()
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        GF2EX_sub(r.x, self.x, y.x)
        _sig_off
        return r

    def __add__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        cdef ntl_GF2EX r = ntl_GF2EX()
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        GF2EX_add(r.x, self.x, y.x)
        _sig_off
        return r

    def __neg__(ntl_GF2EX self):
        cdef ntl_GF2EX r = ntl_GF2EX()
        _sig_on
        GF2EX_negate(r.x, self.x)
        _sig_off
        return r

    def __pow__(ntl_GF2EX self, long e, ignored):
        cdef ntl_GF2EX r = ntl_GF2EX()
        _sig_on
        GF2EX_power(r.x, self.x, e)
        _sig_off
        return r
