"""
Free algebra elements

AUTHOR: David Kohel, 2005-09
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU Genral Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
from sage.misc.misc import repr_lincomb
from sage.rings.integer import Integer
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.algebra_element import AlgebraElement

class FreeAlgebraElement(AlgebraElement):
    def __init__(self, A, x):
        """
        Create the element x of the FreeAlgebra A.
        """
        if isinstance(x, FreeAlgebraElement) and x.parent() == A:
            return x
        AlgebraElement.__init__(self, A)
        R = A.base_ring()
        if isinstance(x, AlgebraElement): #and x.parent() == A.base_ring():
            self.__monomial_coefficients = { A.monoid()(1):R(x) }
        elif isinstance(x, FreeMonoidElement):
            self.__monomial_coefficients = { x:R(1) }
        elif True:
            self.__monomial_coefficients = x
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def __repr__(self):
        v = list(self.__monomial_coefficients.iteritems())
        v.sort()
        mons = [ m for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def __add__(self, y):
        A = self.parent()
        if isinstance(y, (int, long, Integer)):
            z_elt = dict(self.__monomial_coefficients)
            e = A.monoid()(1)
            if z_elt.has_key(e):
                z_elt[e] += A.base_ring()(y)
            else:
                z_elt[e] = A.base_ring()(y)
            z = A(0)
            z.__monomial_coefficients = z_elt
            return z
        if not isinstance(y, FreeAlgebraElement) or not A == y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z_elt = dict(self.__monomial_coefficients)
        for m, c in y.__monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] + c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = c
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def __neg__(self):
        y = self.parent()(0)
        y_elt = {}
        for m, c in self.__monomial_coefficients.iteritems():
            y_elt[m] = -c
        y.__monomial_coefficients = y_elt
        return y

    def __sub__(self, y):
        A = self.parent()
        if isinstance(y, (int, long, Integer)):
            z_elt = dict(self.__monomial_coefficients)
            e = A.monoid()(1)
            if z_elt.has_key(e):
                z_elt[e] += A.base_ring()(-y)
            else:
                z_elt[e] = A.base_ring()(-y)
            z = A(0)
            z.__monomial_coefficients = z_elt
            return z
        if not isinstance(y, FreeAlgebraElement) or not A == y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z_elt = dict(self.__monomial_coefficients)
        for m, c in y.__monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] - c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = -c
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def __mul__(self, y):
        A = self.parent()
        if not isinstance(y, FreeAlgebraElement) or not A == y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z_elt = {}
        for mx, cx in self.__monomial_coefficients.iteritems():
            for my, cy in y.__monomial_coefficients.iteritems():
                z_elt[mx*my] = cx*cy
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def __pow__(self, n):
        if not isinstance(n, (int, long, Integer)):
            raise TypeError, "Argument n (= %s) must be an integer."%n
        if n < 0:
            raise IndexError, "Argument n (= %s) must be positive."%n
        elif n == 0:
            return self.parent()(1)
        elif n == 1:
            return self
        elif n == 2:
            return self * self
        k = n//2
        return self**k * self**(n-k)

