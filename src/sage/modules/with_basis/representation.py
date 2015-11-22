"""
Representations Of A Semigroup

AUTHORS:

- Travis Scrimshaw (2015-11-21): Initial version
"""

####################################################################################
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
####################################################################################

from sage.structure.element import Element
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.groups import Groups
from sage.categories.semigroups import Semigroups
from sage.categories.modules import Modules
from sage.algebras.group_algebra import GroupAlgebra

class Representation(CombinatorialFreeModule):
    """
    Representation of a semigroup.

    INPUT:

    - ``semigroup`` -- a semigroup
    - ``module`` -- a module with a basis
    - ``on_basis`` -- function which takes as input ``g``, ``m``, where
      ``g`` is an element of the semigroup and ``m`` is an element of the
      indexing set for the basis
    - ``left_repr`` -- boolean; whether this is a left or a right
      representation

    EXAMPLES:

    We construct the sign representation of the permutation group::

        sage: G = SymmetricGroup(4)
        sage: M = CombinatorialFreeModule(QQ, ['v'])
        sage: from sage.modules.with_basis.representation import Representation
        sage: on_basis = lambda g,m: M.term(m, g.sign())
        sage: R = Representation(G, M, on_basis)
        sage: x = R.an_element(); x
        2*B['v']
        sage: c,s = G.gens()
        sage: c,s
        ((1,2,3,4), (1,2))
        sage: c * x
        -2*B['v']
        sage: s * x
        -2*B['v']
        sage: c * s * x
        2*B['v']
        sage: (c * s) * x
        2*B['v']

    REFERENCES:

    - :wikipedia:`Group_representation`
    """
    def __init__(self, semigroup, module, on_basis, left_repr=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: R = Representation(G, M, on_basis)
            sage: R._test_representation()
        """
        self._semigroup = semigroup
        self._semigroup_algebra = semigroup.algebra(module.base_ring())
        self._on_basis = on_basis
        self._left_repr = left_repr
        self._module = module
        indices = module.basis().keys()
        cat = Modules(module.base_ring()).WithBasis()
        if 'FiniteDimensional' in module.category().axioms():
            cat = cat.FiniteDimensional()
        CombinatorialFreeModule.__init__(self, module.base_ring(), indices,
                                         category=cat, **module.print_options())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P = Permutations(4)
            sage: M = CombinatorialFreeModule(QQ, ['v'])
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.term(m, g.sign())
            sage: Representation(P, M, on_basis)
            Representation of Standard permutations of 4 indexed by {'v'}
             over Rational Field
        """
        return "Representation of {} indexed by {} over {}".format(
            self._semigroup, self.basis().keys(), self.base_ring())

    def _repr_term(self, b):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._repr_term(b) == SGA._repr_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._repr_term(b)

    def _latex_term(self, b):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: R = SGA.regular_representation()
            sage: all(R._latex_term(b) == SGA._latex_term(b) for b in SGA.basis().keys())
            True
        """
        return self._module._latex_term(b)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: A = G.algebra(ZZ)
            sage: R = A.regular_representation()
            sage: x = A.an_element(); x
            B[()] + 4*B[(1,2,3,4)] + 2*B[(1,4)(2,3)]
            sage: R(x)
            B[()] + 4*B[(1,2,3,4)] + 2*B[(1,4)(2,3)]
        """
        if isinstance(x, Element) and x.parent() is self._module:
            return self._from_dict(x.monomial_coefficients(copy=False), remove_zeros=False)
        return super(Representation, self)._element_constructor_(x)

    def _test_representation(self, **options):
        """
        Check that ``self`` is a representation.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()
        """
        if self._semigroup in Semigroups().Finite():
            sgelts = list(self._semigroup)
        else:
            sgelts = self._semigroup.some_elements()
        tester = self._tester(**options)
        S = tester.some_elements()
        for x in sgelts:
            for y in sgelts:
                for elt in S:
                    if self._left_repr:
                        tester.assertEqual(x*(y*elt), (x*y)*elt)
                    else:
                        tester.assertEqual((elt*y)*x, elt*(y*x))

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: G = groups.misc.WeylGroup(['B',2], prefix='s')
                sage: R = G.regular_representation()
                sage: s1,s2 = G.gens()
                sage: x = R.an_element(); x
                2*B[s2*s1*s2] + B[s1*s2] + 3*B[s2] + B[1]
                sage: s1 * x
                2*B[s2*s1*s2*s1] + 3*B[s1*s2] + B[s1] + B[s2]
                sage: s2 * x
                B[s2*s1*s2] + 2*B[s1*s2] + B[s2] + 3*B[1]
            """
            if isinstance(scalar, Element):
                P = self.parent()
                if scalar.parent() is P._semigroup:
                    if self_on_left == P._left_repr:
                        scalar = ~scalar
                    return P.linear_combination(((P._on_basis(scalar, m), c)
                                                 for m,c in self), not self_on_left)
                R = self.base_ring()
                if scalar.parent() is P._semigroup_algebra:
                    ret = self
                    for ms,cs in scalar:
                        if self_on_left == P._left_repr:
                            ms = ~ms
                        return P.linear_combination(((P._on_basis(ms, m), cs)
                                                     for m,c in self), not self_on_left)
            return CombinatorialFreeModule.Element._acted_upon_(scalar, self_on_left)

        _rmul_ = _lmul_ = _acted_upon_

class RegularRepresentation(Representation):
    """
    The regular representation of a semigroup.

    REFERENCES:

    - :wikipedia:`Regular_representation`
    """
    def __init__(self, semigroup, base_ring, left_repr=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: TestSuite(R).run()
        """
        if left_repr:
            on_basis = self._left_on_basis
        else:
            on_basis = self._right_on_basis
        module = semigroup.algebra(base_ring)
        Representation.__init__(self, semigroup, module, on_basis, left_repr)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.regular_representation()
            Left Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
            sage: G.regular_representation(left_repr=False)
            Right Regular Representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        if self._left_repr:
            base = "Left Regular Representation"
        else:
            base = "Right Regular Representation"
        return base + " of {} over {}".format(self._semigroup, self.base_ring())

    def _left_on_basis(self, g, m):
        """
        Return the left action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation()
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(g*m)

    def _right_on_basis(self, g, m):
        """
        Return the right action of ``g`` on ``m``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: R = G.regular_representation(left_repr=False)
            sage: R._test_representation()  # indirect doctest
        """
        return self.monomial(m*g)

class TrivialRepresentation(CombinatorialFreeModule):
    """
    The trivial representation of a semigroup.

    REFERENCES:

    - :wikipedia:`Trivial_representation`
    """
    def __init__(self, semigroup, base_ring):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = groups.permutation.PGL(2, 3)
            sage: V = G.trivial_representation()
            sage: TestSuite(V).run()
        """
        self._semigroup = semigroup
        self._semigroup_algebra = semigroup.algebra(base_ring)
        cat = Modules(base_ring).WithBasis().FiniteDimensional()
        CombinatorialFreeModule.__init__(self, base_ring, ['v'], category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = groups.permutation.Dihedral(4)
            sage: G.trivial_representation()
            Trivial representation of Dihedral group of order 8
             as a permutation group over Integer Ring
        """
        return "Trivial representation of {} over {}".format(self._semigroup,
                                                             self.base_ring())

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 3)
                sage: V = SGA.trivial_representation()
                sage: x = V.an_element()
                sage: all(x * b == x for b in SGA.basis())
                True
                sage: all(b * x == x for b in SGA.basis())
                True
            """
            if isinstance(scalar, Element):
                if scalar.parent() is self.parent()._semigroup:
                    return self
                if scalar.parent() is self.parent()._semigroup_algebra:
                    d = self.monomial_coefficients(copy=True)
                    d['v'] *= sum(scalar.coefficients())
                    return self.parent()._from_dict(d)
            return CombinatorialFreeModule.Element._acted_upon_(scalar, self_on_left)

        _rmul_ = _lmul_ = _acted_upon_

