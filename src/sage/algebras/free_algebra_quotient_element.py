"""
Free algebra quotient elements

AUTHOR: David Kohel, 2005-09
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
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
from sage.rings.coerce import bin_op
from sage.misc.misc import repr_lincomb
from sage.rings.ring import Ring
from sage.rings.ring_element import RingElement
from sage.ext.integer import Integer
from sage.algebras.algebra_element import AlgebraElement
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import FreeModuleElement
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.free_algebra_element import FreeAlgebraElement
from sage.structure.element import Element_cmp_

def is_FreeAlgebraQuotientElement(x):
    return isinstance(x, FreeAlgebraQuotientElement)

class FreeAlgebraQuotientElement(AlgebraElement, Element_cmp_):
    def __init__(self, A, x):
        """
        Create the element x of the FreeAlgebraQuotient A.
        """
        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() is A:
            return x
        AlgebraElement.__init__(self, A)
        Q = self.parent()
        if isinstance(x, (int, long, Integer)):
            self.__vector = Q.module().gen(0) * x
            return
        elif isinstance(x, FreeModuleElement) and x.parent() is Q.module():
            self.__vector = x
            return
        elif isinstance(x, FreeModuleElement) and x.parent() == A.module():
            self.__vector = x
            return
        R = A.base_ring()
        M = A.module()
        F = A.monoid()
        B = A.monomial_basis()
        if isinstance(x, (int, long, Integer)):
            self.__vector = x*M.gen(0)
        elif isinstance(x, RingElement) and not isinstance(x, AlgebraElement) and x in R:
            self.__vector = x * M.gen(0)
        elif isinstance(x, FreeMonoidElement) and x.parent() is F:
            if x in B:
                self.__vector = M.gen(B.index(x))
            else:
                raise AttributeError, \
                      "Argument x (= %s) is not in monomial basis"%x
        elif isinstance(x, list) and len(x) == A.dimension():
            try:
                self.__vector = M(x)
            except TypeError:
                raise TypeError, "Argument x (= %s) is of the wrong type."%x
        elif isinstance(x, FreeAlgebraElement) and x.parent() is A.free_algebra():
            # Need to do more work here to include monomials not
            # represented in the monomial basis.
            self.__vector = M(0)
            for m, c in x._FreeAlgebraElement__monomial_coefficients.iteritems():
                self.__vector += c*M.gen(B.index(m))
        elif isinstance(x, dict):
            self.__vector = M(0)
            for m, c in x.iteritems():
                self.__vector += c*M.gen(B.index(m))
        elif isinstance(x, AlgebraElement) and x.parent().ambient_algebra() is A:
            self.__vector = x.ambient_algebra_element().vector()
        else:
            print "x =", x
            print "type(x) =", type(x)
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def __repr__(self):
        Q = self.parent()
        cffs = list(self.__vector)
        mons = Q.monomial_basis()
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def vector(self):
        return self.__vector

    def __add__(self, y):
        if not isinstance(y, FreeAlgebraQuotientElement):
            return bin_op(self, y, operator.add)
        A = self.parent()
        if not A is y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z = A(0)
        z.__vector = self.__vector + y.__vector
        return z

    def __radd__(self, y):
        """
        Can be deleted once this moves to the base ring class.
        """
        return bin_op(y, self, operator.add)

    def __neg__(self):
        y = self.parent()(0)
        y.__vector = -self.__vector
        return y

    def __sub__(self, y):
        if not isinstance(y, FreeAlgebraQuotientElement):
            return bin_op(self, y, operator.sub)
        A = self.parent()
        if not A is y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z = A(0)
        z.__vector = self.__vector - y.__vector
        return z

    def __mul__(self, y):
        if not isinstance(y, FreeAlgebraQuotientElement):
            return bin_op(self, y, operator.mul)
        A = self.parent()
        if not A is y.parent():
            raise TypeError, "Argument y (= %s) is of the wrong type."%y
        def monomial_product(X,w,m):
            mats = X._FreeAlgebraQuotient__matrix_action
            for (j,k) in m._element_list:
                M = mats[int(j)]
                for l in range(k): w *= M
            return w
        u = self.__vector.copy()
        v = y.__vector
        z = A(0)
        B = A.monomial_basis()
        for i in range(A.dimension()):
            c = v[i]
            if c != 0: z.__vector += monomial_product(A,c*u,B[i])
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


