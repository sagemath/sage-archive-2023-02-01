r"""
Examples of a finite dimensional Lie algebra with basis
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
from sage.modules.free_module import FreeModule
from sage.structure.parent import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

class AbelianLieAlgebra(Parent, UniqueRepresentation):
    r"""
    An example of a finite dimensional Lie algebra with basis:
    the abelian Lie algebra.

    This class illustrates a minimal implementation of a finite dimensional
    Lie algebra with basis.
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.
        """
        return super(AbelianLieAlgebra, cls).__classcall__(cls, R, tuple(names))

    def __init__(self, R, names):
        """
        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example(); L
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
            sage: TestSuite(L).run()
        """
        self._M = FreeModule(R, len(names))
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        Parent.__init__(self, names=names, category=cat)

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self.base_ring(), self.variable_names())

    def _repr_(self):
        """
        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
        """
        return "An example of a finite dimensional Lie algebra with basis:"\
               " the abelian Lie algebra on the generators "\
               " {} over {}".format(self.lie_algebra_generators(), self.base_ring())

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.basis()
        """
        names = self.variable_names()
        d = {names[i]: self.element_class(self, b)
             for i,b in enumerate(self._M.basis())}
        return Family(d)

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.lie_algebra_generators()
        """
        return self.basis()

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket on basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
            sage: L.bracket_on_basis('a', 'c')
            0
        """
        return self.zero()

    class Element(ElementWrapper):
        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: elt = L.an_element()
                sage: elt.lift()
                2*a + 2*b + 3*c
            """
            UEA = self.parent().universal_enveloping_algebra()
            gens_dict = UEA.gens_dict()
            return UEA.sum(c * gens_dict[t] for t, c in self)

Example = AbelianLieAlgebra

