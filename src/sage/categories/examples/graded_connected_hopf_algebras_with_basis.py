r"""
Examples of graded connected Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.functions.other import binomial
from sage.misc.cachefunc import cached_method
from sage.sets.non_negative_integers import NonNegativeIntegers


class GradedConnectedHopfAlgebraOfInteger(CombinatorialFreeModule):
    r"""
    This class illustrates an implementation of a graded hopf algebra
    with basis: the polynomial Hopf algebra of one variable.

    """
    def __init__(self, base_ring):
        """
        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: TestSuite(H).run()

        """
        CombinatorialFreeModule.__init__(self, base_ring, NonNegativeIntegers(),
                                         category=GradedHopfAlgebrasWithBasis(base_ring).Connected())


    @cached_method
    def one_basis(self):
        """
        Returns 0, which index the unit of the hopf algebra.

        EXAMPLES::r

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: H.one_basis()
            0
            sage: H.one()
            P0

        """
        return self.basis().keys()(0)

    def degree_on_basis(self, t):
        """
        The degree of an integer is the integer

        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: H.degree_on_basis(45)
            45

        """
        return t

    def _repr_(self):
        """
        Print representation

        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: GradedConnectedHopfAlgebraOfInteger(QQ)
            An example of a graded connected hopf algebra with basis over Rational Field

        """
        return "An example of a graded connected hopf algebra with basis over %s" % self.base_ring()

    def _repr_term(self, t):
        """
        Print representation for the basis element represented by the
        integer ``t``.

        EXAMPLES::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: H._repr_term(45)
            'P45'
        """
        return 'P' + repr(t)

    def product_on_basis(self, i, j):
        """
        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: H.monomial(4) * H.monomial(5)
            P9

        """
        return self.monomial(i+j)

    def coproduct_on_basis(self, i):
        """
        TESTS::

            sage: from sage.categories.examples.graded_connected_hopf_algebras_with_basis import GradedConnectedHopfAlgebraOfInteger
            sage: H = GradedConnectedHopfAlgebraOfInteger(QQ)
            sage: H.monomial(3).coproduct()
            P0 # P3 + 3*P1 # P2 + 3*P2 # P1 + P3 # P0

        """
        return self.sum_of_terms(
            ((i-j, j), binomial(i, j))
            for j in range(i+1)
        )

Example = GradedConnectedHopfAlgebraOfInteger
