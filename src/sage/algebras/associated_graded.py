r"""
Associated Graded Algebras

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
    of a filtered-algebra-with-basis `A`.

    Let `A` be a filtered algebra with basis over a
    commutative ring `R`. Let `(F_0, F_1, F_2, \ldots)` be
    its filtration, let `(b_i)_{i \in I}` be its basis,
    and consider the partition of the set `I` into subsets
    `I_0, I_1, I_2, \ldots` which is part of the data of
    a filtered algebra with basis. The *associated graded
    algebra* (or, for short, just *graded algebra*) of
    `A` is a graded algebra with basis defined as follows:

    We know (see :class:`FilteredModulesWithBasis`) that
    `A` (being a filtered `R`-module with basis) canonically
    becomes a graded `R`-module with basis. (Its `n`-th
    graded component, for every `n \in \NN`, is the
    `R`-submodule of `A` spanned by `(b_i)_{i \in I_n}`.)
    We define a multiplication `*` on this graded `R`-module
    `A` (not to be mistaken for the multiplication of the
    original algebra `A`) by requiring that

    .. MATH::

        b_i * b_j = \left( \text{the }
        (n+m)\text{-th homogeneous component of }
        b_i b_j \right)
        \qquad \text{for all } n, m \in \NN \text{ and }
        i \in I_n \text{ and } j \in I_m

    (or, equivalently,

    .. MATH::

        u * v = \left( \text{the }
        (n+m)\text{-th homogeneous component of }
        uv \right)
        \qquad \text{for all } n, m \in \NN
        \text{ and any homogeneous elements } u
        \text{ and } v \text{ of respective degrees }
        n \text{ and } m

    ).

    Thus, `(A, *)` is a graded `R`-algebra with basis.
    This is called the associated graded algebra of `A`,
    and denoted by `\operatorname{gr} A`.

    Notice that the multiplication `*` of this associated
    graded algebra depends not only on the filtered algebra
    `A`, but also on its basis.

    INPUT:

    - ``A`` -- a filtered algebra

    EXAMPLES:

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

