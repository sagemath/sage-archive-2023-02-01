r"""
Examples of a Lie algebra with basis
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import LieAlgebras
from sage.combinat.free_module import CombinatorialFreeModule

class AbelianLieAlgebra(CombinatorialFreeModule):
    r"""
    An example of a Lie algebra: the abelian Lie algebra.

    This class illustrates a minimal implementation of a Lie algebra with
    a distinguished basis.
    """
    def __init__(self, R, gens):
        """
        EXAMPLES::

            sage: A = LieAlgebras(QQ).example(); A
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
            sage: TestSuite(A).run()
        """
        cat = LieAlgebras(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, gens, category=cat)

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        # TODO: Implement using a combinatorial free module with a free abelian monoid
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self.base_ring(), self.variable_names())

    def _repr_(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).example()
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
        """
        return "An example of a Lie algebra: the abelian Lie algebra on the" \
               " generators indexed by {} over {}".format(
                        self.basis().keys(), self.base_ring())

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L.lie_algebra_generators()
            (B['a'], B['b'], B['c'])
        """
        return self.basis()

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket on basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L.bracket_on_basis('a', 'c')
            0
        """
        return self.zero()

    def is_solvable(self):
        """
        Return if ``self`` is a solvable Lie algebra.
        """
        return True

    def is_nilpotent(self):
        """
        Return if ``self`` is a nilpotent Lie algebra.
        """
        return True

    class Element(CombinatorialFreeModule.Element):
        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).example()
                sage: elt = L.an_element()
                sage: elt.lift()
                2*a + 2*b + 3*c
            """
            UEA = self.parent().universal_enveloping_algebra()
            gens_dict = UEA.gens_dict()
            return UEA.sum(c * gens_dict[t] for t, c in self)

Example = AbelianLieAlgebra

