"""
Hecke algebras and modules
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
import weakref

import sage.rings.arith as arith
import sage.rings.infinity
import sage.misc.latex as latex
import sage.categories.all
import module
import hecke_operator
import sage.rings.commutative_algebra
from sage.misc.misc import verbose
from sage.matrix.constructor import matrix

def is_HeckeAlgebra(x):
    return isinstance(x, HeckeAlgebra_base)

# The basis_matrix stuff here is a workaround for a subtle bug discovered by
# me (David Loeffler) 2009-04-13. The problem is that if one creates two
# subspaces of a Hecke module which are equal as subspaces but have different
# bases, then the caching machinery needs to distinguish between them.

# The AttributeError occurs in two distinct ways: if M is not a free module
# over its base ring, it might not have a basis_matrix method; and for
# SupersingularModule objects, the basis_matrix method exists but raises an
# error -- this is a known bug (#4306).

# See the doctest for the __call__ method below, which tests that this caching
# is working as it should.

_anemic_cache = {}
def AnemicHeckeAlgebra(M):
    try:
        k = (M, M.basis_matrix())
    except AttributeError:
        k = M
    if _anemic_cache.has_key(k):
        T = _anemic_cache[k]()
        if not (T is None):
            return T
    T = HeckeAlgebra_anemic(M)
    _anemic_cache[k] = weakref.ref(T)
    return T

_cache = {}
def HeckeAlgebra(M):
    try:
        k = (M, M.basis_matrix())
    except AttributeError:
        k = M
    if _cache.has_key(k):
        T = _cache[k]()
        if not (T is None):
            return T
    T = HeckeAlgebra_full(M)
    _cache[k] = weakref.ref(T)
    return T


class HeckeAlgebra_base(sage.rings.commutative_algebra.CommutativeAlgebra):
    """
    An algebra of Hecke operators on a fixed Hecke module
    """
    def __init__(self, M):
        """
        INPUT:


        -  ``M`` - a Hecke module
        """
        if not module.is_HeckeModule(M):
            raise TypeError, "M (=%s) must be a HeckeModule"%M
        self.__M = M
        sage.rings.commutative_algebra.CommutativeAlgebra.__init__(self, M.base_ring())

    def _repr_(self):
        return "Hecke algebra acting on %s"%self.__M

    def __call__(self, x):
        r"""
        Convert x into an element of this Hecke algebra.

        TESTS:

        We test that coercion is OK between the Hecke algebras associated to two submodules which are equal but have different bases::

            sage: M = CuspForms(Gamma0(57))
            sage: f1,f2,f3 = M.newforms()
            sage: N1 = M.submodule(M.free_module().submodule_with_basis([f1.element().element(), f2.element().element()]))
            sage: N2 = M.submodule(M.free_module().submodule_with_basis([f1.element().element(), (f1.element() + f2.element()).element()]))
            sage: N1.hecke_operator(5).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
            sage: N2.hecke_operator(5).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [-4  1]
            sage: N1.hecke_algebra()(N2.hecke_operator(5)).matrix_form()
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
            sage: N1.hecke_algebra()(N2.hecke_operator(5).matrix_form())
            Hecke operator on Modular Forms subspace of dimension 2 of ... defined by:
            [-3  0]
            [ 0  1]
        """
        try:
            if x.parent() is self:
                return x
            elif hecke_operator.is_HeckeOperator(x):
                if x.parent() == self:
                    return hecke_operator.HeckeOperator(self, x.index())
                else:
                    raise TypeError
            elif hecke_operator.is_HeckeAlgebraElement(x):
                if x.parent() == self:
                    if x.parent().module().basis_matrix() == self.module().basis_matrix():
                        return hecke_operator.HeckeAlgebraElement_matrix(self, x.matrix())
                    else:
                        A = matrix([self.module().coordinate_vector(x.parent().module().gen(i)) \
                            for i in xrange(x.parent().module().rank())])
                        return hecke_operator.HeckeAlgebraElement_matrix(self, ~A * x.matrix() * A)
                else:
                    raise TypeError
            else:
                return hecke_operator.HeckeAlgebraElement_matrix(self, self.__matrix_space()(x))

        except TypeError:
            raise TypeError, "coercion of %s into Hecke algebra not defined."%x

    def _coerce_impl(self, x):
        return self._coerce_try(x, self.__matrix_space())

    def gen(self, n):
        """
        Return the `n`-th Hecke operator.

        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_algebra()
            sage: T.gen(2)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
        """
        return self.hecke_operator(n)

    def ngens(self):
        return sage.rings.infinity.infinity

    def is_noetherian(self):
        """
        Return True if this Hecke algebra is Noetherian as a ring.
        """
        return self.base_ring().is_noetherian()

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: T = ModularSymbols(11).hecke_algebra()
            sage: T.gen(2) in T
            True
            sage: 5 in T
            False
        """
        try:
            return x.parent() == self
        except AttributeError:
            return False

    def __matrix_space(self):
        try:
            return self.__matrix_space_cache
        except AttributeError:
            self.__matrix_space_cache = sage.matrix.matrix_space.MatrixSpace(
                self.base_ring(), self.module().rank())
            return self.__matrix_space_cache

    def _latex_(self):
        return "\\mathbf{T}_{%s}"%latex(self.__M)

    def level(self):
        return self.module().level()

    def module(self):
        """
        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.module()
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        return self.__M

    def rank(self):
        raise NotImplementedError

    def basis(self):
        raise NotImplementedError

    def discriminant(self):
        raise NotImplementedError

    def gens(self):
        r"""
        Return a generator over all Hecke operator `T_n` for
        `n = 1, 2, 3, \ldots`. This is infinite.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: g = T.gens()
            sage: g.next()
            Hecke operator T_1 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
            sage: g.next()
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        n = 1
        while True:
            yield self.hecke_operator(n)
            n += 1

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator `T_n`.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        try:
            return self.__hecke_operator[n]
        except AttributeError:
            self.__hecke_operator = {}
        except KeyError:
            pass
        n = int(n)
        T = self.__M._hecke_operator_class()(self, n)
        self.__hecke_operator[n] = T
        return T

    def hecke_matrix(self, n):
        """
        Return the matrix of the n-th Hecke operator `T_n`.

        EXAMPLES::

            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: T.hecke_matrix(2)
            [ -24    0    0]
            [   0  -24    0]
            [4860    0 2049]
        """
        return self.__M.hecke_matrix(n)


class HeckeAlgebra_full(HeckeAlgebra_base):
    def _repr_(self):
        return "Full Hecke algebra acting on %s"%self.module()

    def __cmp__(self, other):
        if not isinstance(other, HeckeAlgebra_full):
            return -1
        return cmp(self.module(), other.module())

    def is_anemic(self):
        """
        Return True if this is an anemic Hecke algebra.

        EXAMPLES:
        """
        return False


class HeckeAlgebra_anemic(HeckeAlgebra_base):
    def _repr_(self):
        return "Anemic Hecke algebra acting on %s"%self.module()

    def _latex_(self):
        return "\\mathbf{T}'_{%s}"%latex(self.__M)

    def __cmp__(self, other):
        if not isinstance(HeckeAlgebra_anemic):
            return -1
        return cmp(self.module(), other.module())

    def hecke_operator(self, n):
        """
        Return the `n`-th Hecke operator, for `n` any
        positive integer coprime to the level.

        EXAMPLES::

            sage: T = ModularSymbols(Gamma1(5),3).anemic_hecke_algebra()
            sage: T.hecke_operator(2)
            Hecke operator T_2 on Modular Symbols space of dimension 4 for Gamma_1(5) of weight 3 with sign 0 and over Rational Field
            sage: T.hecke_operator(5)
            Traceback (most recent call last):
            ...
            IndexError: Hecke operator T_5 not defined in the anemic Hecke algebra
        """
        n = int(n)
        if arith.gcd(self.module().level(), n) != 1:
            raise IndexError, "Hecke operator T_%s not defined in the anemic Hecke algebra"%n
        return self.module()._hecke_operator_class()(self, n)

    def is_anemic(self):
        return True

    def gens(self):
        """
        Return a generator over all Hecke operator `T_n` for
        `n = 1, 2, 3, \ldots`, with `n` coprime to the
        level. This is an infinite sequence.

        EXAMPLES::

            sage: T = ModularSymbols(12,2).anemic_hecke_algebra()
            sage: g = T.gens()
            sage: g.next()
            Hecke operator T_1 on Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Rational Field
            sage: g.next()
            Hecke operator T_5 on Modular Symbols space of dimension 5 for Gamma_0(12) of weight 2 with sign 0 over Rational Field
        """
        level = self.level()
        n = 1
        while True:
            if arith.gcd(n, level) == 1:
                yield self.hecke_operator(n)
            n += 1




