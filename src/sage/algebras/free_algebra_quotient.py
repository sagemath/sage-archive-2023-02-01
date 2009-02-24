"""
Free algebra quotients

TESTS::

    sage: n = 2
    sage: A = FreeAlgebra(QQ,n,'x')
    sage: F = A.monoid()
    sage: i, j = F.gens()
    sage: mons = [ F(1), i, j, i*j ]
    sage: r = len(mons)
    sage: M = MatrixSpace(QQ,r)
    sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]), M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]) ]
    sage: H2.<i,j> = A.quotient(mons,mats)
    sage: H2 == loads(dumps(H2))
    True
    sage: i == loads(dumps(i))
    True
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

from sage.rings.integer import Integer
from sage.modules.free_module import FreeModule
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra import is_FreeAlgebra
from sage.algebras.free_algebra_quotient_element import FreeAlgebraQuotientElement
from sage.structure.parent_gens import ParentWithGens

class FreeAlgebraQuotient(Algebra, object):
    def __init__(self, A, mons, mats, names):
        """
        Returns a quotient algebra defined via the action of a free algebra
        A on a (finitely generated) free module. The input for the quotient
        algebra is a list of monomials (in the underlying monoid for A)
        which form a free basis for the module of A, and a list of
        matrices, which give the action of the free generators of A on this
        monomial basis.

        EXAMPLES:

        Quaternion algebra defined in terms of three generators::

            sage: n = 3
            sage: A = FreeAlgebra(QQ,n,'i')
            sage: F = A.monoid()
            sage: i, j, k = F.gens()
            sage: mons = [ F(1), i, j, k ]
            sage: M = MatrixSpace(QQ,4)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]),  M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]),  M([0,0,0,1, 0,0,-1,0, 0,1,0,0, -1,0,0,0]) ]
            sage: H3.<i,j,k> = FreeAlgebraQuotient(A,mons,mats)
            sage: x = 1 + i + j + k
            sage: x
            1 + i + j + k
            sage: x**128
            -170141183460469231731687303715884105728 + 170141183460469231731687303715884105728*i + 170141183460469231731687303715884105728*j + 170141183460469231731687303715884105728*k

        Same algebra defined in terms of two generators, with some penalty
        on already slow arithmetic.

        ::

            sage: n = 2
            sage: A = FreeAlgebra(QQ,n,'x')
            sage: F = A.monoid()
            sage: i, j = F.gens()
            sage: mons = [ F(1), i, j, i*j ]
            sage: r = len(mons)
            sage: M = MatrixSpace(QQ,r)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]), M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]) ]
            sage: H2.<i,j> = A.quotient(mons,mats)
            sage: k = i*j
            sage: x = 1 + i + j + k
            sage: x
            1 + i + j + i*j
            sage: x**128
            -170141183460469231731687303715884105728 + 170141183460469231731687303715884105728*i + 170141183460469231731687303715884105728*j + 170141183460469231731687303715884105728*i*j
        """
        if not is_FreeAlgebra(A):
            raise TypeError, "Argument A must be an algebra."
        R = A.base_ring()
#        if not R.is_field():  # TODO: why?
#            raise TypeError, "Base ring of argument A must be a field."
        n = A.ngens()
        assert n == len(mats)
        self.__free_algebra = A
        self.__ngens = n
        self.__dim = len(mons)
        self.__module = FreeModule(R,self.__dim)
        self.__matrix_action = mats
        self.__monomial_basis = mons # elements of free monoid
        ParentWithGens.__init__(self, R, names, normalize=True)

    def __eq__(self,right):
        return type(self) == type(right) and \
               self.ngens() == right.ngens() and \
               self.rank() == right.rank() and \
               self.module() == right.module() and \
               self.matrix_action() == right.matrix_action() and \
               self.monomial_basis() == right.monomial_basis()


    def __call__(self, x):
        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() is self:
            return x
        return FreeAlgebraQuotientElement(self,x)

    def _coerce_impl(self, x):
        """
        Return the coercion of x into this free algebra quotient.

        The algebras that coerce into this quotient ring canonically, are:

        - this quotient algebra

        - anything that coerces into the algebra of which this is the
          quotient
        """
        return self._coerce_try(x, [self.__free_algebra, self.base_ring()])

    def _repr_(self):
        R = self.base_ring()
        n = self.__ngens
        r = self.__module.dimension()
        x = self.variable_names()
        return "Free algebra quotient on %s generators %s and dimension %s over %s"%(n,x,r,R)

    def __contains__(self, x):
        return isinstance(x, FreeAlgebraQuotientElement) and x.parent() == self

    def gen(self,i):
        """
        The i-th generator of the algebra.
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        R = self.base_ring()
        F = self.__free_algebra.monoid()
        n = self.__ngens
        return FreeAlgebraQuotientElement(self,{F.gen(i):R(1)})

    def ngens(self):
        """
        The number of generators of the algebra.
        """
        return self.__ngens

    def dimension(self):
        """
        The rank of the algebra (as a free module).
        """
        return self.__dim

    def matrix_action(self):
        return self.__matrix_action

    def monomial_basis(self):
        return self.__monomial_basis

    def rank(self):
        """
        The rank of the algebra (as a free module).
        """
        return self.__dim

    def module(self):
        """
        The free module of the algebra.
        """
        return self.__module

    def monoid(self):
        """
        The free monoid of generators of the algebra.
        """
        return self.__free_algebra.monoid()

    def monomial_basis(self):
        """
        The free monoid of generators of the algebra as elements of a free
        monoid.
        """
        return self.__monomial_basis

    def free_algebra(self):
        """
        The free algebra generating the algebra.
        """
        return self.__free_algebra


