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
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.examples.lie_algebras import LieAlgebraFromAssociative as BaseExample

class AbelianLieAlgebra(Parent, UniqueRepresentation):
    r"""
    An example of a finite dimensional Lie algebra with basis:
    the abelian Lie algebra.

    This class illustrates a minimal implementation of a finite dimensional
    Lie algebra with basis.
    """
    @staticmethod
    def __classcall_private__(cls, R, names, M=None):
        """
        Normalize input to ensure a unique representation.
        """
        if isinstance(names, str):
            names = names.split(',')
        if M is None:
            M = FreeModule(R, len(names))
        elif len(names) != M.dimension():
            raise ValueError("number of generators is not correct")
        else:
            M = M.change_ring(R)
        return super(AbelianLieAlgebra, cls).__classcall__(cls, R, tuple(names), M)

    def __init__(self, R, names, M):
        """
        EXAMPLES::

            sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example(); L
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
            sage: TestSuite(L).run()
        """
        self._M = M
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        Parent.__init__(self, base=R, names=names, category=cat)

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
        ret = "An example of a finite dimensional Lie algebra with basis:" \
              " the abelian Lie algebra with generators {!r} over {}".format(
                         self.variable_names(), self.base_ring())
        B = self._M.basis_matrix()
        if not B.is_one():
            ret += " with basis matrix:\n{!r}".format(B)
        return ret

    def subalgebra(self, gens, names='x'):
        """
        Return a subalgebra of ``self``.
        """
        if isinstance(names, str):
            names = names.split(',')
        if len(names) == 1 and len(gens) != 1:
            names = tuple( names[0] + str(i) for i in range(len(gens)) )
        N = self._M.subspace([g.value for g in gens])
        return AbelianLieAlgebra(self.base_ring(), names, N)

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

    def gens(self):
        """
        Return the generators of ``self``.
        """
        G = self.lie_algebra_generators()
        return tuple(G[i] for i in self.variable_names())

    def free_module(self):
        """
        Return ``self`` as a free module.
        """
        return self._M

    class Element(BaseExample.Element):
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.
            """
            return self.parent().zero()

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
            gens = UEA.gens()
            return UEA.sum(c * gens[i] for i, c in self.value.iteritems())

        def to_vector(self):
            """
            Return ``self`` as a vector.
            """
            return self.value

Example = AbelianLieAlgebra

