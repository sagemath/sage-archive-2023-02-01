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

from ntl_ZZ import unpickle_class_value

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
    Minimal wrapper of NTL's GF2EX class.
    """
    def __init__(self, x=[]):
        """
        Minimal wrapper of NTL's GF2EX class.

        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
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
        """
        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0 1] [1 0 0 1] [1]]')
            sage: f == loads(dumps(f))
            True
        """
        return unpickle_class_value, (ntl_GF2EX, self.__repr__())

    def __cmp__(self, other):
        """
        Compare self to other.

        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0 1] [1 0 0 1] [1]]')
            sage: g = ntl.GF2EX('[[1 0 1] [1 1] [1] [0 0 1]]')
            sage: f == f
            True
            sage: f == g
            False
        """
        ## TODO: this is shady. fix that.
        if (type(self) != type(other)):
            return cmp(type(self), type(other))
        return cmp(self.__repr__(), other.__repr__())

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: m = ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: ntl.GF2EX('[[1 0] [2 1]]').__repr__()
            '[[1] [0 1]]'
        """
        return GF2EX_to_PyString(&self.x)

    def __mul__(ntl_GF2EX self, other):
        """
        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0] [2 1]]')
            sage: g = ntl.GF2EX('[[1 0 1 1] [0 1 1 0 1] [1 0 1]]')
            sage: f*g ## indirect doctest
            [[1 0 1 1] [0 0 1 1] [1 0 0 1 0 1] [0 1 0 1]]
        """
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
        """
        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0] [2 1]]')
            sage: g = ntl.GF2EX('[[1 0 1 1] [0 1 1 0 1] [1 0 1]]')
            sage: f-g ## indirect doctest
            [[0 0 1 1] [0 0 1 0 1] [1 0 1]]
        """
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
        """
        EXAMPLES:
            sage: ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0] [2 1]]')
            sage: g = ntl.GF2EX('[[1 0 1 1] [0 1 1 0 1] [1 0 1]]')
            sage: f+g ## indirect doctest
            [[0 0 1 1] [0 0 1 0 1] [1 0 1]]
        """
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
        """
        EXAMPLES:
            sage: m = ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0] [2 1]]')
            sage: -f ## indirect doctest
            [[1] [0 1]]
        """
        cdef ntl_GF2EX r = ntl_GF2EX()
        _sig_on
        GF2EX_negate(r.x, self.x)
        _sig_off
        return r

    def __pow__(ntl_GF2EX self, long e, ignored):
        """
        EXAMPLES:
            sage: m = ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1]))
            sage: f = ntl.GF2EX('[[1 0] [2 1]]')
            sage: f**2 ## indirect doctest
            [[1] [] [0 0 1]]
        """
        cdef ntl_GF2EX r = ntl_GF2EX()
        _sig_on
        GF2EX_power(r.x, self.x, e)
        _sig_off
        return r
