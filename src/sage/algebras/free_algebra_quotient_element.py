"""
Free algebra quotient elements

AUTHORS:

- David Kohel (2005-09)
"""

from __future__ import with_statement

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
from sage.misc.misc import repr_lincomb
from sage.rings.ring import Ring
from sage.rings.ring_element import RingElement
from sage.rings.integer import Integer
from sage.algebras.algebra_element import AlgebraElement
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import FreeModuleElement
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.free_algebra_element import FreeAlgebraElement
from sage.structure.parent_gens import localvars

def is_FreeAlgebraQuotientElement(x):
    return isinstance(x, FreeAlgebraQuotientElement)

class FreeAlgebraQuotientElement(AlgebraElement):
    def __init__(self, A, x):
        """
        Create the element x of the FreeAlgebraQuotient A.
        """
        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() is A:
            return x

        AlgebraElement.__init__(self, A)
        Q = self.parent()

        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() == Q:
            self.__vector = Q.module()(x.vector())
            return
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
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def _repr_(self):
        Q = self.parent()
        M = Q.monoid()
        with localvars(M, Q.variable_names()):
            cffs = list(self.__vector)
            mons = Q.monomial_basis()
            x = repr_lincomb(mons, cffs).replace("*1 "," ")
            if x[len(x)-2:] == "*1":
                return x[:len(x)-2]
            else:
                return x

    def _latex_(self):
        Q = self.parent()
        M = Q.monoid()
        with localvars(M, Q.variable_names()):
            cffs = list(self.__vector)
            mons = Q.monomial_basis()
            x = repr_lincomb(mons, cffs, True).replace("*1 "," ")
            if x[len(x)-2:] == "*1":
                return x[:len(x)-2]
            else:
                return x

    def vector(self):
        return self.__vector

    def __cmp__(self, right):
        return cmp(self.vector(),right.vector())

    def __neg__(self):
        y = self.parent()(0)
        y.__vector = -self.__vector
        return y

    def _add_(self, y):
        A = self.parent()
        z = A(0)
        z.__vector = self.__vector + y.__vector
        return z

    def _sub_(self, y):
        A = self.parent()
        z = A(0)
        z.__vector = self.__vector - y.__vector
        return z

    def _mul_(self, y):
        A = self.parent()
        def monomial_product(X,w,m):
            mats = X._FreeAlgebraQuotient__matrix_action
            for (j,k) in m._element_list:
                M = mats[int(j)]
                for l in range(k): w *= M
            return w
        u = self.__vector.__copy__()
        v = y.__vector
        z = A(0)
        B = A.monomial_basis()
        for i in range(A.dimension()):
            c = v[i]
            if c != 0: z.__vector += monomial_product(A,c*u,B[i])
        return z

    def _rmul_(self, c):
        return self.parent([c*a for a in self.__vector])

    def _lmul_(self, c):
        return self.parent([a*c for a in self.__vector])




