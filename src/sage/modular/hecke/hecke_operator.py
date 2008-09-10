"""
Hecke operators
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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


import math
import operator

import sage.algebras.algebra_element
import sage.categories.all as cat
import sage.misc.latex as latex
import sage.misc.misc as misc
import sage.modules.module
import sage.modules.free_module_morphism as free_module_morphism
import sage.modular.dims as dims
import sage.rings.arith as arith
from   sage.rings.integer import Integer

import algebra
import morphism


def is_HeckeOperator(x):
    return isinstance(x, HeckeOperator)

def is_HeckeAlgebraElement(x):
    return isinstance(x, HeckeAlgebraElement)

class HeckeAlgebraElement(sage.algebras.algebra_element.AlgebraElement):
    def __init__(self, parent):
        if not algebra.is_HeckeAlgebra(parent):
            raise TypeError, "parent (=%s) must be a Hecke algebra"%parent
        sage.algebras.algebra_element.AlgebraElement.__init__(self, parent)

    def domain(self):
        return self.parent().module()

    def codomain(self):
        return self.parent().module()

    def hecke_module_morphism(self):
        """
        Return the endomorphism of Hecke modules defined
        by the matrix attached to this Hecke operator.

        EXAMPLES:
            sage: M = ModularSymbols(Gamma1(13))
            sage: t = M.hecke_operator(2)
            sage: t
            Hecke operator T_2 on Modular Symbols space of dimension 15 for Gamma_1(13) of weight 2 with sign 0 and over Rational Field
            sage: t.hecke_module_morphism()
            Hecke module morphism T_2 defined by the matrix
            (not printing 15 x 15 matrix)
            Domain: Modular Symbols space of dimension 15 for Gamma_1(13) of weight ...
            Codomain: Modular Symbols space of dimension 15 for Gamma_1(13) of weight ...
        """
        try:
            return self.__hecke_module_morphism
        except AttributeError:
            T = self.matrix()
            M = self.domain()
            H = cat.End(M)
            if isinstance(self, HeckeOperator):
                name = "T_%s"%self.index()
            else:
                name = ""
            self.__hecke_module_morphism = morphism.HeckeModuleMorphism_matrix(H, T, name)
            return self.__hecke_module_morphism

    def __is_compatible(self, other):
        return isinstance(other, HeckeAlgebraElement) and self.parent() == other.parent()

    def _add_(self, other):
        """
        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: t = M.hecke_operator(2)
            sage: t
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: t + t
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [ 6  0 -2]
            [ 0 -4  0]
            [ 0  0 -4]

        We can also add Hecke operators with different indexes:

            sage: M = ModularSymbols(Gamma1(6),4)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2 - t3
            Hecke operator on Modular Symbols space of dimension 6 for Gamma_1(6) of weight 4 with sign 0 and over Rational Field defined by:
            (not printing 6 x 6 matrix)
            sage: (t2 - t3).charpoly('x')
            x^6 + 36*x^5 + 104*x^4 - 3778*x^3 + 7095*x^2 - 3458*x
        """
        return self.parent()(self.matrix() + other.matrix())

    def __call__(self, x):
        """
        Apply this Hecke operator to $x$.

        EXAMPLES:
            sage: M = ModularSymbols(11); t2 = M.hecke_operator(2)
            sage: t2(M.gen(0))
            3*(1,0) - (1,9)

            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t3(t2(M.gen(0)))
            12*(1,0) - 2*(1,9)
            sage: (t3*t2)(M.gen(0))
            12*(1,0) - 2*(1,9)
        """
        T = self.hecke_module_morphism()
        return T(x)

    def __rmul__(self, left):
        """
        EXAMPLES:
            sage: M = ModularSymbols(11); t2 = M.hecke_operator(2)
            sage: 2*t2
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [ 6  0 -2]
            [ 0 -4  0]
            [ 0  0 -4]
        """
        return self.parent()(left * self.matrix())

    def _sub_(self, other):
        """
        Compute the difference of self and other.

        EXAMPLES:
            sage: M = ModularSymbols(Gamma1(6),4)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2 - t3
            Hecke operator on Modular Symbols space of dimension 6 for Gamma_1(6) of weight 4 with sign 0 and over Rational Field defined by:
            (not printing 6 x 6 matrix)
        """
        return self.parent()(self.matrix() - other.matrix())

    def apply_sparse(self, x):
        """
        Apply this Hecke operator to x, where we avoid computing the matrix
        of x if possible.

        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: T = M.hecke_operator(23)
            sage: T.apply_sparse(M.gen(0))
            24*(1,0) - 5*(1,9)
        """
        if not x in self.domain():
            raise TypeError, "x (=%s) must be in %s"%(x, self.domain())
        T = self.hecke_module_morphism()
        return T(x)

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this Hecke operator.

        INPUT:
            var -- string (default: 'x')

        OUTPUT:
            a monic polynomial in the given variable.

        EXAMPLES:
            sage: M = ModularSymbols(Gamma1(6),4)
            sage: M.hecke_operator(2).charpoly('x')
            x^6 - 14*x^5 + 29*x^4 + 172*x^3 - 124*x^2 - 320*x + 256
        """
        return self.matrix().charpoly(var)

    def decomposition(self):
        """
        Decompose the Hecke module under the action of this Hecke
        operator.

        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: t2 = M.hecke_operator(2)
            sage: t2.decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            ]

            sage: M = ModularSymbols(33, sign=1).new_submodule()
            sage: T = M.hecke_operator(2)
            sage: T.decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field
            ]
        """
        try:
            return self.__decomposition
        except AttributeError:
            pass
        if isinstance(self, HeckeOperator) and \
               arith.gcd(self.index(), self.domain().level()) == 1:
            D = self.hecke_module_morphism().decomposition(is_diagonalizable=True)
        else:
            # TODO: There are other weaker hypotheses that imply diagonalizability.
            D = self.hecke_module_morphism().decomposition()
        D.sort()
        D.set_immutable()
        self.__decomposition = D
        return D

    def det(self):
        """
        Return the determinant of this Hecke operator.

        EXAMPLES:
            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.det()
            100
        """
        return self.hecke_module_morphism().det()

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of
        this Hecke operator.

        EXAMPLES:
            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
        """
        return self.hecke_module_morphism().fcp(var)

    def image(self):
        """
        Return the image of this Hecke operator.

        EXAMPLES:
            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
            sage: T.image()
            Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T-4).image()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T**2-5).image()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_module_morphism().image()

    def kernel(self):
        """
        Return the kernel of this Hecke operator.

        EXAMPLES:
            sage: M = ModularSymbols(23)
            sage: T = M.hecke_operator(3)
            sage: T.fcp('x')
            (x - 4) * (x^2 - 5)^2
            sage: T.kernel()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T-4).kernel()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: (T**2-5).kernel()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_module_morphism().kernel()

    def trace(self):
        """
        Return the trace of this Hecke operator.
            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2)
            sage: T.trace()
            2001
        """
        return self.hecke_module_morphism().trace()

    def __getitem__(self, ij):
        """
        EXAMPLE:
            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2).matrix_form()
            sage: T[0,0]
            -24
        """
        return self.matrix()[ij]


class HeckeAlgebraElement_matrix(HeckeAlgebraElement):
    def __init__(self, parent, A):
        HeckeAlgebraElement.__init__(self, parent)
        if not sage.matrix.all.is_Matrix(A):
            raise TypeError, "A must be a matrix"
        self.__matrix = A

    def __cmp__(self, other):
        if not isinstance(other, HeckeAlgebraElement_matrix):
            if isinstance(other, HeckeOperator):
                return cmp(self, other.matrix_form())
            return sage.rings.coerce.cmp(self, other)
        c = cmp(self.parent(), other.parent())
        if c: return c
        return cmp(self.__matrix, other.__matrix)

    def _repr_(self):
        if max(self.__matrix.nrows(),self.__matrix.ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.__matrix.nrows(), self.__matrix.ncols())
        else:
            mat = str(self.__matrix)
        return "Hecke operator on %s defined by:\n%s"%(self.parent().module(), mat)

    def _latex_(self):
        return latex(self.__matrix)

    def matrix(self):
        """
        Return the matrix that defines this Hecke algebra element.

        EXAMPLES:
            sage: M = ModularSymbols(1,12)
            sage: T = M.hecke_operator(2).matrix_form()
            sage: T.matrix()
            [ -24    0    0]
            [   0  -24    0]
            [4860    0 2049]
        """
        return self.__matrix

    def _mul_(self, other):
        return self.parent()(other.matrix() * self.matrix())



class HeckeOperator(HeckeAlgebraElement):
    def __init__(self, parent, n):
        """
        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: M.hecke_operator(2005)
            Hecke operator T_2005 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field

        We create a Hecke operator of large index (greater than 32 bits):
            sage: M1 =  ModularSymbols(21,2)
            sage: M1.hecke_operator(13^9)
            Hecke operator T_10604499373 on Modular Symbols space of dimension 5 for Gamma_0(21) of weight 2 with sign 0 over Rational Field
        """
        HeckeAlgebraElement.__init__(self, parent)
        if not isinstance(n, (int,long,Integer)):
            raise TypeError, "n must be an int"
        self.__n = int(n)

    def __cmp__(self, other):
        if not isinstance(other, HeckeOperator):
            if isinstance(other, HeckeAlgebraElement_matrix):
                return cmp(self.matrix_form(), other)
            return sage.rings.coerce.cmp(self, other)
        c = cmp(self.parent(), other.parent())
        if c: return c
        if self.__n == other.__n:
            return 0
        return cmp(self.matrix(), other.matrix())

    def _repr_(self):
        return "Hecke operator T_%s on %s"%(self.__n, self.domain())

    def _latex_(self):
        return "T_{%s}"%self.__n

    def __mul__(self, other):
        """
        EXAMPLES:
        We create the space of modular symbols of level $11$ and weight $2$, then compute
        $T_2$ and $T_3$ on it, along with their composition.

            sage: M = ModularSymbols(11)
            sage: t2 = M.hecke_operator(2); t3 = M.hecke_operator(3)
            sage: t2*t3
            Hecke operator T_6 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: t3.matrix() * t2.matrix()
            [12  0 -2]
            [ 0  2  0]
            [ 0  0  2]
            sage: (t2*t3).matrix()
            [12  0 -2]
            [ 0  2  0]
            [ 0  0  2]

        When we compute $T_2^2$ the result is not (easily seen to be)
        a Hecke operator of the form $T_n$, so it is returned as a
        Hecke module homomorphism defined as a matrix:

            sage: t2**5
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [243   0 -55]
            [  0 -32   0]
            [  0   0 -32]
        """
        if isinstance(other, HeckeOperator) and other.parent() == self.parent():
            n = None
            if arith.gcd(self.__n, other.__n) == 1:
                n = self.__n * other.__n
            else:
                P = set(arith.prime_divisors(self.domain().level()))
                if P.issubset(set(arith.prime_divisors(self.__n))) and \
                   P.issubset(set(arith.prime_divisors(other.__n))):
                    n = self.__n * other.__n
            if n:
                return HeckeOperator(self.parent(), n)
        # otherwise
        return self.matrix_form() * other

    def index(self):
        """
        Return the index of this Hecke operator, i.e., if
        this Hecke operator is $T_n$, return the int $n$.

        EXAMPLES:
            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.index()
            17
        """
        return self.__n

    def matrix(self):
        """
        Return the matrix underlying this Hecke operator.

        EXAMPLES:
            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.matrix()
            [18  0 -4]
            [ 0 -2  0]
            [ 0  0 -2]
        """
        try:
            return self.__matrix
        except AttributeError:
            self.__matrix = self.parent().hecke_matrix(self.__n)
            return self.__matrix

    def matrix_form(self):
        """
        Return the matrix form of this element of a Hecke algebra.
            sage: T = ModularSymbols(11).hecke_operator(17)
            sage: T.matrix_form()
            Hecke operator on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field defined by:
            [18  0 -4]
            [ 0 -2  0]
            [ 0  0 -2]
        """
        try:
            return self.__matrix_form
        except AttributeError:
            self.__matrix_form = self.parent()(self.matrix())
            return self.__matrix_form


