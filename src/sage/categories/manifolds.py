r"""
Manifolds
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
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
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields

class Manifolds(Category_over_base_ring):
    r"""
    The category of manifolds over any field.

    Let `k` be a topological field. A `d`-dimensional `k`-*manifold* `M`
    is a second countable Hausdorff space such that the neighborhood of
    any point `x \in M` is homeomorphic to `k^d`.

    EXAMPLES::

        sage: from sage.categories.manifolds import Manifolds
        sage: C = Manifolds(RR); C
        Category of manifolds over Real Field with 53 bits of precision
        sage: C.super_categories()
        [Category of topological spaces]

    TESTS::

        sage: TestSuite(C).run()
    """
    def __init__(self, base, name=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: C = Manifolds(RR)
            sage: TestSuite(C).run()
        """
        if base not in Fields().Topological():
            raise ValueError("base must be a topological field")
        Category_over_base_ring.__init__(self, base, name)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(RR).super_categories()
            [Category of topological spaces]
        """
        return [Sets().Topological()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of manifolds defines no new
        structure: a morphism of metric spaces between manifolds
        is a manifold morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(RR).additional_structure()
        """
        return None

    class ParentMethods:
        @abstract_method
        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: M = Manifolds(RR).example()
                sage: M.dimension()
                3
            """

    class SubcategoryMethods:
        @cached_method
        def Connected(self):
            """
            Return the full subcategory of the connected objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Connected()
                Category of connected manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Connected()).run()
                sage: Manifolds(RR).Connected.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Connected')

        @cached_method
        def FiniteDimensional(self):
            """
            Return the full subcategory of the finite dimensional
            objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: C = Manifolds(RR).Connected().FiniteDimensional(); C
                Category of finite dimensional connected manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: from sage.categories.manifolds import Manifolds
                sage: C = Manifolds(RR).Connected().FiniteDimensional()
                sage: TestSuite(C).run()
                sage: Manifolds(RR).Connected().FiniteDimensional.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('FiniteDimensional')

        @cached_method
        def Real(self):
            """
            Return the subcategory of manifolds over `\RR` of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Real()
                Category of real manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Real()).run()
                sage: Manifolds(RR).Real.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Real')

        @cached_method
        def Complex(self):
            """
            Return the subcategory of manifolds over `\CC` of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Complex()
                Category of complex manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Complex()).run()
                sage: Manifolds(RR).Complex.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Complex')

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        """
        Category of finite dimensional manifolds.
        """

    class Connected(CategoryWithAxiom_over_base_ring):
        """
        The category of connected manifolds.
        """

    class Real(CategoryWithAxiom_over_base_ring):
        """
        The category of manifolds over `\RR`.
        """
        class SubcategoryMethods:
            @cached_method
            def Complex(self):
                r"""
                Raise an error as a manifold over `\RR` is not a manifold
                over `\CC`.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Complex()
                    Traceback (most recent call last):
                    ...
                    TypeError: a real manifold is not a complex manifold
                """
                raise TypeError("a real manifold is not a complex manifold")

            @cached_method
            def Differentiable(self):
                """
                Return the subcategory of the differentiable objects
                of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Differentiable()
                    Category of differentiable real manifolds
                     over Real Field with 53 bits of precision

                TESTS::

                    sage: TestSuite(Manifolds(RR).Real().Differentiable()).run()
                    sage: Manifolds(RR).Real().Differentiable.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('Differentiable')

            @cached_method
            def Smooth(self):
                """
                Return the subcategory of the smooth objects of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Smooth()
                    Category of smooth real manifolds
                     over Real Field with 53 bits of precision

                TESTS::

                    sage: TestSuite(Manifolds(RR).Real().Smooth()).run()
                    sage: Manifolds(RR).Real().Smooth.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('Smooth')

            @cached_method
            def Analytic(self):
                """
                Return the subcategory of the analytic objects of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Analytic()
                    Category of analytic real manifolds
                     over Real Field with 53 bits of precision

                TESTS::

                    sage: TestSuite(Manifolds(RR).Real().Analytic()).run()
                    sage: Manifolds(RR).Real().Analytic.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('Analytic')

            @cached_method
            def AlmostComplex(self):
                """
                Return the subcategory of the almost complex objects
                of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().AlmostComplex()
                    Category of almost complex real manifolds
                     over Real Field with 53 bits of precision

                TESTS::

                    sage: TestSuite(Manifolds(RR).Real().AlmostComplex()).run()
                    sage: Manifolds(RR).Real().AlmostComplex.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('AlmostComplex')

        class Differentiable(CategoryWithAxiom_over_base_ring):
            """
            The category of differentiable manifolds over `\RR`.

            A `d`-dimensional differentiable manifold is a manifold whose
            underlying vector space is `\RR^d` and differentiable atlas.
            """

        class Smooth(CategoryWithAxiom_over_base_ring):
            """
            The category of smooth manifolds over `\RR`.

            A `d`-dimensional differentiable manifold is a manifold whose
            underlying vector space is `\RR^d` and smooth atlas.
            """
            def extra_super_categories(self):
                """
                Return the extra super categories of ``self``.

                A smooth manifold is differentiable.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Smooth().super_categories() # indirect doctest
                    [Category of differentiable real manifolds
                     over Real Field with 53 bits of precision]
                """
                return [Manifolds(self.base()).Real().Differentiable()]

        class Analytic(CategoryWithAxiom_over_base_ring):
            r"""
            The category of complex manifolds.

            A `d`-dimensional analytic manifold is a manifold whose underlying
            vector space is `\RR^d` and an analytic atlas.
            """
            def extra_super_categories(self):
                """
                Return the extra super categories of ``self``.

                An analytic manifold is smooth.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Analytic().super_categories() # indirect doctest
                    [Category of smooth real manifolds
                     over Real Field with 53 bits of precision]
                """
                return [Manifolds(self.base()).Real().Smooth()]

        class AlmostComplex(CategoryWithAxiom_over_base_ring):
            r"""
            The category of almost complex manifolds.

            A `d`-dimensional almost complex manifold `M` is a manifold
            whose underlying vector space is `\RR^d` with a smooth tensor
            field `J` of rank `(1, 1)` such that `J^2 = -1` when regarded as a
            vector bundle isomorphism `J : TM \to TM` on the tangent bundle.
            The tensor field `J` is called the almost complex structure of `M`.
            """
            def extra_super_categories(self):
                """
                Return the extra super categories of ``self``.

                An analytic manifold is smooth.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Real().Analytic().super_categories() # indirect doctest
                    [Category of smooth real manifolds
                     over Real Field with 53 bits of precision]
                """
                return [Manifolds(self.base()).Real().Smooth()]

    class Complex(CategoryWithAxiom_over_base_ring):
        r"""
        The category of complex manifolds.

        A `d`-dimensional complex manifold is a manifold whose underlying
        vector space is `\CC^d` and a holomorphic atlas.
        """
        class SubcategoryMethods:
            @cached_method
            def Real(self):
                r"""
                Raise an error as a manifold over `\RR` is not a manifold
                over `\CC`.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds(RR).Complex().Real()
                    Traceback (most recent call last):
                    ...
                    TypeError: a complex manifold is not a real manifold
                """
                raise TypeError("a complex manifold is not a real manifold")

