r"""
Incidence Algebras
"""
#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
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

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.algebras import Algebras
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.matrix_space import MatrixSpace

from copy import copy

class IncidenceAlgebra(CombinatorialFreeModule):
    r"""
    The incidence algebra of a poset.
    """
    def __init__(self, R, P, prefix='I'):
        """
        Initialize ``self``.

        TESTS::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: TestSuite(I).run()
        """
        cat = Algebras(R).WithBasis()
        if P in FiniteEnumeratedSets():
            cat = cat.FiniteDimensional()
        self._poset = P
        CombinatorialFreeModule.__init__(self, R, map(tuple, P.intervals()),
                                         prefix=prefix, category=cat)

    def _repr_term(self, A):
        """
        Return a string representation of the term labeled by ``A``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I._repr_term((4, 12))
            'I[4, 12]'
        """
        return self.prefix() + str(list(A))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: P.incidence_algebra(QQ)
            Incidence algebra of Finite lattice containing 16 elements
             over Rational Field
        """
        return "Incidence algebra of {} over {}".format(self._poset, self.base_ring())

    def poset(self):
        """
        Return the defining poset of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.poset()
            Finite lattice containing 16 elements
            sage: I.poset() == P
            True
        """
        return self._poset

    def product_on_basis(self, A, B):
        r"""
        Return the product of basis elements indexed by ``A`` and ``B``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.product_on_basis((1, 3), (3, 11))
            I[1, 11]
            sage: I.product_on_basis((1, 3), (2, 2))
            0
        """
        if A[1] == B[0]:
            return self.monomial((A[0], B[1]))
        return self.zero()

    @cached_method
    def one(self):
        """
        Return the element `1` in ``self`` (which is the Kronecker
        delta `\delta(x, y)`).

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.one()
            I[0, 0] + I[1, 1] + I[2, 2] + I[3, 3] + I[4, 4] + I[5, 5]
             + I[6, 6] + I[7, 7] + I[8, 8] + I[9, 9] + I[10, 10]
             + I[11, 11] + I[12, 12] + I[13, 13] + I[14, 14] + I[15, 15]
        """
        return self.sum_of_monomials((x, x) for x in self._poset)

    delta = one

    @cached_method
    def zeta(self):
        r"""
        Return the `\zeta` function in ``self``.

        The `\zeta` function on a poset `P` is given by

        .. MATH::

            \zeta(x, y) = \begin{cases}
            1 & x \leq y, \\
            0 & x \not\leq y.
            \end{cases}

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I.zeta() * I.mobius() == I.one()
            True
        """
        return self.sum(self.basis())

    @cached_method
    def mobius(self):
        """
        Return the Mobius function of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(2)
            sage: I = P.incidence_algebra(QQ)
            sage: I.mobius()
            I[0, 0] - I[0, 1] - I[0, 2] + I[0, 3] + I[1, 1]
             - I[1, 3] + I[2, 2] - I[2, 3] + I[3, 3]
        """
        mu = self._poset.mobius_function
        R = self.base_ring()
        return self.sum_of_terms((A, R(mu(*A))) for A in self.basis().keys())

    def __getitem__(self, A):
        """
        Return the basis element indexed by ``A``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I[2, 6]
            I[2, 6]
            sage: I[(2, 6)]
            I[2, 6]
            sage: I[[2, 6]]
            I[2, 6]
            sage: I[2]
            I[2, 2]

            sage: I[2, 5]
            Traceback (most recent call last):
            ...
            ValueError: not an interval

            sage: I[-1]
            Traceback (most recent call last):
            ...
            ValueError: not an element of the poset
        """
        if not isinstance(A, (list, tuple)):
            if A not in self._poset.list():
                raise ValueError("not an element of the poset")
            A = (A, A)
        else:
            A = tuple(A)
            if len(A) != 2 or A not in self.basis().keys():
                raise ValueError("not an interval")
        return self.monomial(A)

    @lazy_attribute
    def _linear_extension(self):
        """
        Return a fixed linear extension of the defining poset of ``self``.

        EXAMPLES::

            sage: P = posets.BooleanLattice(4)
            sage: I = P.incidence_algebra(QQ)
            sage: I._linear_extension
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
        """
        return tuple(self._poset.linear_extension())

    class Element(CombinatorialFreeModule.Element):
        """
        An element of an incidence algebra.
        """
        def __call__(self, x, y):
            """
            Return ``self(x, y)``.

            EXAMPLES::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.mobius()
                sage: mu(0, 12)
                1
                sage: mu(0, 7)
                -1
                sage: mu(15, 0)
                0
            """
            P = self.parent()._poset
            x = P(x)
            y = P(y)
            return self[x,y]

        @cached_method
        def to_matrix(self):
            r"""
            Return ``self`` as a matrix.

            We define a matrix `M_{xy} = \alpha(x, y)` for some element
            `\alpha \in I_P` in the incidence algebra `I_P` and we order
            the elements `x,y \in P` by some linear extension of `P`. This
            defines an algebra (iso)morphism; in particular, multiplication
            in the incidence algebra goes to matrix multiplication.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: I.mobius().to_matrix()
                [ 1 -1 -1  1]
                [ 0  1  0 -1]
                [ 0  0  1 -1]
                [ 0  0  0  1]
                sage: I.zeta().to_matrix()
                [1 1 1 1]
                [0 1 0 1]
                [0 0 1 1]
                [0 0 0 1]

            TESTS:

            We check that this is an algebra (iso)morphism::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.mobius()
                sage: (mu*mu).to_matrix() == mu.to_matrix() * mu.to_matrix()
                True
            """
            P = self.parent()
            MS = MatrixSpace(P.base_ring(), P._poset.cardinality(), sparse=True)
            L = P._linear_extension
            M = copy(MS.zero())
            for i,c in self:
                M[L.index(i[0]), L.index(i[1])] = c
            M.set_immutable()
            return M

        def is_unit(self):
            """
            Return if ``self`` is a unit.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.mobius()
                sage: mu.is_unit()
                True
                sage: zeta = I.zeta()
                sage: zeta.is_unit()
                True
                sage: x = mu - I.zeta() + I[2,2]
                sage: x.is_unit()
                False
                sage: y = I.mobius() + I.zeta()
                sage: y.is_unit()
                True

            This depends on the base ring::

                sage: I = P.incidence_algebra(ZZ)
                sage: y = I.mobius() + I.zeta()
                sage: y.is_unit()
                False
            """
            return all(self[x,x].is_unit() for x in self.parent()._poset)

        def __invert__(self):
            """
            Return the inverse of ``self`` if it exists.

            EXAMPLES::

                sage: P = posets.BooleanLattice(2)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.mobius()
                sage: ~mu
                I[0, 0] + I[0, 1] + I[0, 2] + I[0, 3] + I[1, 1]
                 + I[1, 3] + I[2, 2] + I[2, 3] + I[3, 3]
                sage: x = mu - I.zeta() + I[2,2]
                sage: ~x
                Traceback (most recent call last):
                ...
                ValueError: element is not invertible

            TESTS::

                sage: P = posets.BooleanLattice(4)
                sage: I = P.incidence_algebra(QQ)
                sage: mu = I.mobius()
                sage: ~mu == I.zeta()
                True
                sage: ~I.one() == I.one()
                True
            """
            M = self.to_matrix()
            if not M.is_invertible():
                raise ValueError("element is not invertible")
            inv = ~M
            L = self.parent()._linear_extension
            return self.parent().sum_of_terms( ((L[i], L[j]), inv[i,j])
                       for i,j in inv.nonzero_positions(copy=False) )

