"""
Free algebra quotients
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

class FreeAlgebraQuotient(Algebra, object):
    def __init__(self, A, mons, mats, names = None): #, gens_function = None):
        """
        Returns a quotient algebra defined via the action of a free algebra A
        on a (finitely generated) free module.  The input for the quotient algebra
        is a list monomials (in the underlying monoid for A) which form a free
        basis for the module of A, and a list of matrices, which give the action
        of the free generators of A on this monomial basis.

        EXAMPLES:

        Quaternion algebra defined in terms of three generators:

            sage: n = 3
            sage: A = FreeAlgebra(QQ,n); B = A.gens()
            sage: F = A.monoid()
            sage: i, j, k = F.gens()
            sage: mons = [ F(1), i, j, k ]
            sage: M = MatrixSpace(QQ,4)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]),  M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]),  M([0,0,0,1, 0,0,-1,0, 0,1,0,0, -1,0,0,0]) ]
            sage: H3 = FreeAlgebraQuotient(A,mons,mats)
            sage: H3.assign_names(["i","j","k"])
            sage: i, j, k = H3.gens()
            sage: x = 1 + i + j + k
            sage: x
            1 + i + j + k
            sage: x**128
            -170141183460469231731687303715884105728 + 170141183460469231731687303715884105728*i + 170141183460469231731687303715884105728*j + 170141183460469231731687303715884105728*k

        Same algebra defined in terms of two generators, with some
        penalty on already slow arithmetic.

            sage: n = 2
            sage: A = FreeAlgebra(QQ,n); B = A.gens()
            sage: F = A.monoid()
            sage: i, j = F.gens()
            sage: mons = [ F(1), i, j, i*j ]
            sage: r = len(mons)
            sage: M = MatrixSpace(QQ,r)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]), M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]) ]
            sage: H2 = FreeAlgebraQuotient(A,mons,mats)
            sage: i, j = H2.gens(); k = i*j
            sage: H2.assign_names(["i","j"])
            sage: x = 1 + i + j + k
            sage: x
            1 + i + j + i*j
            sage: x**128
            -170141183460469231731687303715884105728 + 170141183460469231731687303715884105728*i + 170141183460469231731687303715884105728*j + 170141183460469231731687303715884105728*i*j

        """
        if not is_FreeAlgebra(A):
            raise TypeError, "Argument A must be an algebra."
        R = A.base_ring()
        if not R.is_field():
            raise TypeError, "Base ring of argument A must be a field."
        n = A.ngens()
        assert n == len(mats)
        self.__free_algebra = A
        self.__base_ring = R
        self.__ngens = n
        self.__dim = len(mons)
        self.__module = FreeModule(R,self.__dim)
        self.__matrix_action = mats
        self.__monomial_basis = mons # elements of free monoid
        self.assign_names(names)

    def __call__(self, x):
        if isinstance(x, FreeAlgebraQuotientElement) and x.parent() is self:
            return x
        return FreeAlgebraQuotientElement(self,x)

    def __repr__(self):
        R = self.__base_ring
        n = self.__ngens
        r = self.__module.dimension()
        x = self.__free_algebra.variable_names()
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
        R = self.__base_ring
        F = self.__free_algebra.monoid()
        n = self.__ngens
        return FreeAlgebraQuotientElement(self,{F.gen(i):R(1)})

    def ngens(self):
        """
        The number of generators of the algebra.
        """
        return self.__ngens

    def base_ring(self):
        return self.__base_ring

    def dimension(self):
        """
        The rank of the algebra (as a free module).
        """
        return self.__dim

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
        The free monoid of generators of the algebra as elements
        of a free monoid.
        """
        return self.__monomial_basis

    def free_algebra(self):
        """
        The free algebra generating the algebra.
        """
        return self.__free_algebra

    def assign_names(self,names):
        """
        Assign the printing names for the generators; this will have the unfortunate
        effect of overwriting the names for the covering algebra; this also does not
        overwrite the return value of names() for the Algebra.
        """
        self.monoid().assign_names(names)

    def variable_names(self):
        """
        I override this in order to return to names of the underlying monoid.
        """
        return self.monoid().variable_names()

