r"""
Associated Graded Algebras To Filtered Algebras

AUTHORS:

- Travis Scrimshaw (2014-10-08): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from copy import copy

from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.free_module import CombinatorialFreeModule

class AssociatedGradedAlgebra(CombinatorialFreeModule):
    r"""
    The associated graded algebra `\operatorname{gr} A`
    of a filtered algebra with basis `A`.

    Let `A` be a filtered algebra with basis over a commutative
    ring `R`. Let `(F_i)_{i \in I}` be the filtration of `A`, and
    define

    .. MATH::

        G_i = F_i / \sum_{j < i} F_j

    and then

    .. MATH::

        \operatorname{gr} A = \bigoplus_{i \in I} G_i.

    There are canonical projections `p_i : F_i \to G_i` for
    every `i \in I`. Moreover `\operatorname{gr} A` is naturally a
    graded module with `G_i` being the `i`-th graded component.

    Let `u \in G_i` and `v \in G_j` with lifts `u' \in F_i`
    and `v' \in F_j` respectively. Therefore we define
    multiplication `*` on `\operatorname{gr} A` (not to be mistaken
    for the multiplication of the original algebra `A`) by

    .. MATH::

        u * v = p_{i+j}(u' v').

    The *associated graded algebra* (or, for short, just *graded algebra*)
    of `A` is the graded algebra `\operatorname{gr} A`.

    In particular, let `(b_x)_{x \in X}` be the basis of `A`,
    and consider the partition of the set `X = \bigsqcup_{i \in I} X_i`,
    which is part of the data of a filtered algebra with basis.
    We know (see
    :class:`~sage.categories.filtered_modules_with_basis.FilteredModulesWithBasis`)
    that `A` (being a filtered `R`-module with basis) is canonically
    isomorphic to `\operatorname{gr} A` as an `R`-module. Therefore
    the `k`-th graded component `G_k` can be identified with
    the span of `(b_x)_{x \in X_k}`, or equivalently the
    `k`-th homogeneous component of `A`. Suppose
    that `u' v' = \sum_{k \leq i+j} m_k` where `m_k \in G_k` (which
    has been identified with the `k`-th homogeneous component of `A`).
    Then `u * v = m_{i+j}`. We also note that the choice of
    identification of `G_k` with the `k`-th homogeneous component
    depends on the given basis.

    INPUT:

    - ``A`` -- a filtered algebra

    EXAMPLES::

        sage: A = Algebras(QQ).WithBasis().Filtered().example()
        sage: grA = A.graded_algebra()
        sage: grA.category()
        Category of graded algebras with basis over Rational Field
        sage: x,y,z = map(lambda s: grA.algebra_generators()[s], ['x','y','z'])
        sage: x
        bar(U['x'])
        sage: y * x + z
        bar(U['x']*U['y']) + bar(U['z'])
        sage: A(y) * A(x) + A(z)
        U['x']*U['y']

    We note that the conversion between ``A`` and ``grA`` is
    the canonical ``QQ``-module isomorphism stemming from the
    fact that the underlying ``QQ``-modules of ``A`` and
    ``grA`` are the same::

        sage: grA(A.an_element())
        bar(U['x']^2*U['y']^2*U['z']^3)
        sage: elt = A.an_element() + A.algebra_generators()['x'] + 2
        sage: grelt = grA(elt); grelt
        bar(U['x']^2*U['y']^2*U['z']^3) + bar(U['x']) + 2*bar(1)
        sage: A(grelt) == elt
        True

    .. TODO::

        The algebra ``A`` must currently be an instance of (a subclass of)
        :class:`CombinatorialFreeModule`. This should work with any algebra
        with a basis.

    .. TODO::

        Implement a version of for filtered algebras without a
        distinguished basis.

    REFERENCES:

    - :wikipedia:`Filtered_algebra#Associated_graded_algebra`
    """
    def __init__(self, A, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: TestSuite(grA).run(elements=[prod(grA.algebra_generators())])
        """
        if A not in AlgebrasWithBasis(A.base_ring()).Filtered():
            raise ValueError("the base algebra must be filtered")
        self._A = A

        if category is None:
            category = A.category().Graded()
        opts = copy(A.print_options())
        if not opts['prefix'] and not opts['bracket']:
            opts['bracket'] = '('
        opts['prefix'] = opts['prefix'] + 'bar'

        CombinatorialFreeModule.__init__(self, A.base_ring(), A.indices(),
                                         category=category, **opts)

        # Setup the conversion back
        phi = self.module_morphism(lambda x: A.monomial(x), codomain=A)
        self._A.register_conversion(phi)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: A.graded_algebra()
            Graded Algebra of An example of a filtered algebra with basis:
             the universal enveloping algebra of Lie algebra of RR^3
             with cross product over Rational Field
        """
        return "Graded Algebra of {}".format(self._A)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: latex(A.graded_algebra())
            \operatorname{gr} ...
        """
        from sage.misc.latex import latex
        return "\\operatorname{gr} " + latex(self._A)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        This constructs an element from the filtered algebra `A`
        by the canonical module isomorphism (stemming from the
        fact that `A` and the associated graded algebra `A`
        have the same underlying `R`-module).

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: grA(A.an_element())
            bar(U['x']^2*U['y']^2*U['z']^3)
            sage: grA(A.an_element() + A.algebra_generators()['x'] + 2)
            bar(U['x']^2*U['y']^2*U['z']^3) + bar(U['x']) + 2*bar(1)
        """
        if isinstance(x, CombinatorialFreeModule.Element):
            if x.parent() is self._A:
                return self._from_dict(dict(x))
        return super(AssociatedGradedAlgebra, self)._element_constructor_(x)

    def gen(self, *args, **kwds):
        """
        Return a generator of ``self``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: grA.gen('x')
            bar(U['x'])
        """
        try:
            x = self._A.gen(*args, **kwds)
        except AttributeError:
            x = self._A.algebra_generators()[args[0]]
        return self(x)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: grA.algebra_generators()
            Finite family {'y': bar(U['y']), 'x': bar(U['x']), 'z': bar(U['z'])}
        """
        G = self._A.algebra_generators()
        return Family(G.keys(), lambda x: self(G[x]), name="generator")

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: grA.one_basis()
            1
        """
        return self._A.one_basis()

    def product_on_basis(self, x, y):
        """
        Return the product on basis elements given by the
        indices ``x`` and ``y``.

        EXAMPLES::

            sage: A = Algebras(QQ).WithBasis().Filtered().example()
            sage: grA = A.graded_algebra()
            sage: G = grA.algebra_generators()
            sage: x,y,z = G['x'], G['y'], G['z']
            sage: x * y # indirect doctest
            bar(U['x']*U['y'])
            sage: y * x
            bar(U['x']*U['y'])
            sage: z * y * x
            bar(U['x']*U['y']*U['z'])
        """
        ret = self._A.product_on_basis(x, y)
        deg = self._A.degree_on_basis(x) + self._A.degree_on_basis(y)
        return self.sum_of_terms([(i,c) for i,c in ret
                                     if self._A.degree_on_basis(i) == deg],
                                 distinct=True)

