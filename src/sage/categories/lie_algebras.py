r"""
Lie Algebras

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.distributive_magmas_and_additive_magmas import DistributiveMagmasAndAdditiveMagmas
from sage.categories.finite_enuemrated_sets import FiniteEnumeratedSets
from sage.categories.modules import Modules
from sage.categories.sets_cat import Sets
from sage.categories.homset import Hom
from sage.categories.morphism import Morphism
from sage.structure.sage_object import have_same_parent
from sage.structure.element import get_coercion_model, coerce_binop

class LieAlgebras(Category_over_base_ring):
    """
    The category of Lie algebras.

    EXAMPLES::

        sage: C = LieAlgebras(QQ); C
        Category of Lie algebras over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of vector spaces over Rational Field]

    We construct a typical parent in this category, and do some
    computations with it::

        sage: A = C.example(); A
        An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field

        sage: A.category()
        Category of Lie algebras with basis over Rational Field

        sage: A.base_ring()
        Rational Field
        sage: A.basis().keys()
        {'a', 'b', 'c'}

        sage: (a,b,c) = A.gens()
        sage: a^3, b^2
        (a^3, b^2)
        sage: a*c*b
        a*b*c

        sage: A.product
        <bound method FreeAlgebra_with_category._product_from_product_on_basis_multiply of An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field>
        sage: A.product(a*b,b)
        a*b^2

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_characteristic() . . . pass
        running ._test_distributivity() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.lie_algebras.AbelianLieAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.algeras.lie_algebras.AbelianLieAlgebra_with_category.element_class'>

    Please see the source code of `A` (with ``A??``) for how to
    implement other Lie algebras.

    TESTS::

        sage: TestSuite(LieAlgebras(QQ)).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).super_categories()
            [Category of vector spaces over Rational Field]
        """
        # We do not also derive from (Magmatic) algebras since we don't want *
        #   to be our Lie bracket
        # Also this doesn't inherit the ability to add axioms like Associative
        #   and Unital, both of which do not make sense for Lie algebras
        return [Modules(self.base_ring())]

    def example(self, gens=None):
        """
        Return an example of a Lie algebra as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: LieAlgebras(QQ).example()
            An example of a Lie algebra: the abelian Lie algebra on the generators ('a', 'b', 'c') over Rational Field

        Another set of generators can be specified as an optional argument::

            sage: LieAlgebras(QQ).example(('x','y','z'))
            An example of a Lie algebra: the abelian Lie algebra on the generators ('a', 'b', 'c') over Rational Field
        """
        if gens is None:
            from sage.combinat.partition import Partitions
            gens = Partitions()
        from sage.categories.examples.lie_algebras import Example
        return Example(self.base_ring(), gens)

    WithBasis = LazyImport('sage.categories.lie_algebras_with_basis',
                           'LieAlgebrasWithBasis', as_name='WithBasis')

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        WithBasis = LazyImport('sage.categories.finite_dimensional_lie_algebras_with_basis',
                               'FiniteDimensionalLieAlgebrasWithBasis', as_name='WithBasis')

        def extra_super_categories(self):
            """
            Implements the fact that a finite dimensional Lie algebra over
            a finite ring is finite.

            EXAMPLES::

                sage: Modules(IntegerModRing(4)).FiniteDimensional().extra_super_categories()
                [Category of finite sets]
                sage: Modules(ZZ).FiniteDimensional().extra_super_categories()
                []
                sage: Modules(GF(5)).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(ZZ).FiniteDimensional().is_subcategory(Sets().Finite())
                False
            """
            if self.base_ring() in Sets().Finite():
                return [Sets().Finite()]
            return []

    class ParentMethods:
        def bracket(self, lhs, rhs):
            """
            Return the Lie bracket ``[lhs, rhs]`` after coercing ``lhs`` and
            ``rhs`` into elements of ``self``.
            """
            return self(lhs)._bracket_(self(rhs))

        # Do not override this, instead implement _construct_UEA() in order
        #   to automatically setup the coercion
        def universal_enveloping_algebra(self):
            """
            Return the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
                sage: L.universal_enveloping_algebra()
                Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            """
            return self.lift.codomain()

        @abstract_method(optional=True)
        def _construct_UEA(self):
            """
            Construct the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
                sage: L.universal_enveloping_algebra() # indirect doctest
                Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            """

        @lazy_attribute
        def lift(self):
            """
            Construct the lift morphism from ``self`` to the universal
            enveloping algebra of ``self``.
            """
            M = LiftMorphism(self, self._construct_UEA())
            M.register_as_coercion()
            return M

        @abstract_method(optional=True)
        def killing_form(self, x, y):
            """
            Return the Killing form of ``x`` and ``y``.
            """

        def is_abelian(self):
            r"""
            Return ``True`` if this Lie algebra is abelian.

            A Lie algebra `\mathfrak{g}` is abelian if `[x, y] = 0` for all
            `x, y \in \mathfrak{g}`.

            EXAMPLES::

                sage: L.<x> = LieAlgebra(QQ,1)
                sage: L.is_abelian()
                True
                sage: L.<x,y> = LieAlgebra(QQ,2)
                sage: L.is_abelian()
                False
            """
            G = self.lie_algebra_generators()
            if G not in FiniteEnumeratedSets():
                raise NotImplementedError("infinite number of generators")
            zero = self.zero()
            return all(x._bracket_(y) == zero for x in G for y in G)

        def is_commutative(self):
            """
            Return if ``self`` is commutative. This is equivalent to ``self``
            being abelian.

            EXAMPLES::

                sage: L.<x> = LieAlgebra(QQ, 1)
                sage: L.is_commutative()
                True
            """
            return self.is_abelian()

        def is_field(self, proof=True):
            """
            Return ``False`` since Lie algebras are never a field since
            they are not associative and antisymmetric.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 1, 'x', abelian=True)
                sage: L.is_field()
                False
            """
            return False

        @abstract_method(optional=True)
        def is_solvable(self):
            """
            Return if ``self`` is a solvable Lie algebra.
            """

        @abstract_method(optional=True)
        def is_nilpotent(self):
            """
            Return if ``self`` is a nilpotent Lie algebra.
            """

        def _test_jacobi_identity(self, **options):
            """
            Test that the Jacobi identity is satisfied.
            """
            tester = self._tester(**options)
            elts = self.some_elements()
            jacobi = lambda x, y, z: self.bracket(x, self.bracket(y, z)) + \
                self.bracket(y, self.bracket(z, x)) + \
                self.bracket(z, self.bracket(x, y))
            zero = self.zero()
            for x in elts:
                for y in elts:
                    if x == y:
                        continue
                    for z in elts:
                        if x == z or y == z: # Trivial
                            continue
                        tester.assert_(jacobi(x, y, z) == zero)

        def _test_distributivity(self, **options):
            r"""
            Test the distributivity of the Lie bracket `[,]` on `+` on (not
            necessarily all) elements of this set.

            INPUT::

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES:

            By default, this method runs the tests only on the
            elements returned by ``self.some_elements()``::

                sage: LieAlgebra(QQ, 3, 'x,y,z').some_elements()
                [x]
                sage: LieAlgebra(QQ, 3, 'x,y,z')._test_distributivity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
                sage: h1 = L.gen(0)
                sage: h2 = L.gen(1)
                sage: e2 = L.gen(3)
                sage: L._test_distributivity(elements=[h1, h2, e2])

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            P = S[0].parent()
            from sage.combinat.cartesian_product import CartesianProduct
            for x,y,z in tester.some_elements(CartesianProduct(S,S,S)):
                # left distributivity
                tester.assert_(P.bracket(x, (y + z)) == P.bracket(x, y) + P.bracket(x, z))
                # right distributivity
                tester.assert_(P.bracket((x + y), z) == P.bracket(x, z) + P.bracket(y, z))

    class ElementMethods:
        @coerce_binop
        def bracket(self, rhs):
            """
            Return the Lie bracket ``[self, rhs]``.
            """
            return self._bracket_(rhs)

        # Implement this in order to avoid having to deal with the coercions
        @abstract_method(optional=True)
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::
            """

        def lift(self):
            """
            Lift ``self`` into an element of the universal enveloping algebra.

            EXAMPLES::

                sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
                sage: x.lift()
                x
            """
            return self.parent().lift(self)

        def killing_form(self, x):
            """
            Return the Killing form of ``self`` and ``x``.
            """
            return self.parent().killing_form(self, x)

class LiftMorphism(Morphism):
    """
    The natural lifting morphism from a Lie algebra to its universal
    enveloping algebra.
    """
    def __init__(self, domain, codomain):
        """
        Initialize ``self``.
        """
        Morphism.__init__(self, Hom(domain, codomain))

    def _call_(self, x):
        """
        Lift ``x`` to the universal enveloping algebra.
        """
        return x.lift()

    def section(self):
        """
        Return the section map of ``self``.
        """
        raise NotImplementedError

