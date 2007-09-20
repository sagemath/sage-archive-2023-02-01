"""
FLINT wrapper

AUTHORS:
   - Robert Bradshaw (2007-09-15) Initial version.
"""


#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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


include "../../ext/python_sequence.pxi"

from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer

cdef class Fmpz_poly(SageObject):

    def __new__(self, v=None):
        fmpz_poly_init(self.poly)

    def __init__(self, v):
        """
        Construct a new fmpz_poly from a sequence, constant coefficient,
        or string (in the same format as it prints).

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: Fmpz_poly([1,2,3])
            3  1 2 3
            sage: Fmpz_poly(5)
            1  5
            sage: Fmpz_poly(str(Fmpz_poly([3,5,7])))
            3  3 5 7
        """
        cdef Py_ssize_t i
        cdef long c
        if PY_TYPE_CHECK(v, str):
            if fmpz_poly_from_string(self.poly, v):
                return
            else:
                raise ValueError, "Unable to create Fmpz_poly from that string."
        if not PySequence_Check(v):
            v = [v]
        try:
            for i from 0 <= i < len(v):
                fmpz_poly_set_coeff_si(self.poly, i, v[i])
        except OverflowError:
            raise ValueError, "No fmpz_poly_set_coeff_mpz() method."

    def __dealloc__(self):
        fmpz_poly_clear(self.poly)

    def __setitem__(self, i, value):
        """
        Set the $i$-th item of self, which is the coefficient of the $x^i$ term.

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly(range(10))
            sage: f[7] = 100; f
            10  0 1 2 3 4 5 6 100 8 9
        """
        fmpz_poly_set_coeff_si(self.poly, i, value)

    def __getitem__(self, i):
        """
        Return the $i$-th item of self, which is the coefficient of the $x^i$ term.

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly(range(100))
            sage: f[13]
            13
            sage: f[200]
            0
        """
        cdef Integer res = <Integer>PY_NEW(Integer)
        fmpz_poly_get_coeff_mpz(res.value, self.poly, i)
        return res

    def __repr__(self):
        """
        Print self according to the native FLINT format.

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([0,1]); f^7
            8  0 0 0 0 0 0 0 1
        """
        cdef char* ss = fmpz_poly_to_string(self.poly)
        s = ss
        free(ss)
        return s

    def degree(self):
        """
        The degree of self.

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([1,2,3]); f
            3  1 2 3
            sage: f.degree()
            2
            sage: Fmpz_poly(range(1000)).degree()
            999
            sage: Fmpz_poly([2,0]).degree()
            0
        """
        return fmpz_poly_degree(self.poly)

    def list(self):
        """
        Return self as a list of coefficients, lowest terms first.

        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([2,1,0,-1])
            sage: f.list()
            [2, 1, 0, -1]
        """
        return [self[i] for i in xrange(self.degree()+1)]

    def __add__(left, right):
        if not PY_TYPE_CHECK(left, Fmpz_poly) or not PY_TYPE_CHECK(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>PY_NEW(Fmpz_poly)
        raise NotImplementedError, "no fmpz_poly_add" # fmpz_poly_add(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __sub__(left, right):
        if not PY_TYPE_CHECK(left, Fmpz_poly) or not PY_TYPE_CHECK(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>PY_NEW(Fmpz_poly)
        raise NotImplementedError, "no fmpz_poly_sub" # fmpz_poly_sub(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __neg__(self):
        cdef Fmpz_poly res = <Fmpz_poly>PY_NEW(Fmpz_poly)
        raise NotImplementedError, "no _fmpz_poly_neg" #_fmpz_poly_neg(res.poly, self.poly)
        return res

    def __mul__(left, right):
        """
        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([0,1]); g = Fmpz_poly([2,3,4])
            sage: f*g
            4  0 2 3 4
            sage: f = Fmpz_poly([1,0,-1]); g = Fmpz_poly([2,3,4])
            sage: f*g
            5  2 3 2 -3 -4
        """
        if not PY_TYPE_CHECK(left, Fmpz_poly) or not PY_TYPE_CHECK(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>PY_NEW(Fmpz_poly)
        fmpz_poly_mul(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __pow__(self, n, dummy):
        """
        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: f**6
            7  1 6 15 20 15 6 1
            sage: f = Fmpz_poly([2])
            sage: f^150
            1  1427247692705959881058285969449495136382746624
            sage: 2^150
            1427247692705959881058285969449495136382746624
        """
        cdef long nn = n
        if not PY_TYPE_CHECK(self, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>PY_NEW(Fmpz_poly)
        fmpz_poly_power(res.poly, (<Fmpz_poly>self).poly, nn)
        return res

    def truncate(self, n):
        """
        EXAMPLES:
            sage: from sage.libs.flint.fmpz_poly import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: g = f**10; g
            11  1 10 45 120 210 252 210 120 45 10 1
            sage: g.truncate(5); g
            5  1 10 45 120 210
        """
        cdef long nn = n
        fmpz_poly_truncate(self.poly, nn) # mutating!

