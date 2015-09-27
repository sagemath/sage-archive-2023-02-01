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
from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields

class Manifolds(Category_singleton):
    r"""
    The category of manifolds over any field.

    Let `k` be a topological field. A `d`-dimensional `k`-*manifold* `M`
    is a second countable Hausdorff space such that the neighborhood of
    any point `x \in M` is homeomorphic to `k^d`.

    EXAMPLES::

        sage: from sage.categories.manifolds import Manifolds
        sage: C = Manifolds(); C
        Category of manifolds
        sage: C.super_categories()
        [Category of topological spaces]

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds().super_categories()
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
            sage: Manifolds().additional_structure()
        """
        return None

    class ParentMethods:
        @abstract_method
        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: M = Manifolds().example()
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
                sage: Manifolds().Connected()
                Category of connected manifolds

            TESTS::

                sage: TestSuite(Manifolds().Connected()).run()
                sage: Manifolds().Connected.__module__
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
                sage: C = Manifolds().Connected().FiniteDimensional(); C
                Category of finite dimensional connected manifolds

            TESTS::

                sage: from sage.categories.manifolds import Manifolds
                sage: C = Manifolds().Connected().FiniteDimensional()
                sage: TestSuite(C).run()
                sage: Manifolds().Connected().FiniteDimensional.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('FiniteDimensional')

        @cached_method
        def Real(self):
            """
            Return the subcategory of manifolds over `\RR` of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds().Real()
                Category of real manifolds

            TESTS::

                sage: TestSuite(Manifolds().Real()).run()
                sage: Manifolds().Real.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Real')

        @cached_method
        def Complex(self):
            """
            Return the subcategory of manifolds over `\CC` of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds().Complex()
                Category of complex manifolds

            TESTS::

                sage: TestSuite(Manifolds().Complex()).run()
                sage: Manifolds().Complex.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Complex')

    class FiniteDimensional(CategoryWithAxiom):
        """
        Category of finite dimensional manifolds.
        """

    class Connected(CategoryWithAxiom):
        """
        The category of connected manifolds.
        """

    class Real(CategoryWithAxiom):
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
                    sage: Manifolds().Real().Complex()
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
                    sage: Manifolds().Real().Differentiable()
                    Category of differentiable real manifolds

                TESTS::

                    sage: TestSuite(Manifolds().Real().Differentiable()).run()
                    sage: Manifolds().Real().Differentiable.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('Differentiable')

            @cached_method
            def Smooth(self):
                """
                Return the subcategory of the smooth objects of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds().Real().Smooth()
                    Category of smooth real manifolds

                TESTS::

                    sage: TestSuite(Manifolds().Real().Smooth()).run()
                    sage: Manifolds().Real().Smooth.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('Smooth')

            @cached_method
            def Analytic(self):
                """
                Return the subcategory of the analytic objects of ``self``.

                EXAMPLES::

                    sage: from sage.categories.manifolds import Manifolds
                    sage: Manifolds().Real().Analytic()
                    Category of analytic real manifolds

                TESTS::

                    sage: TestSuite(Manifolds().Real().Analytic()).run()
                    sage: Manifolds().Real().Analytic.__module__
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
                    sage: Manifolds().Real().AlmostComplex()
                    Category of almost complex real manifolds

                TESTS::

                    sage: TestSuite(Manifolds().Real().AlmostComplex()).run()
                    sage: Manifolds().Real().AlmostComplex.__module__
                    'sage.categories.manifolds'
                """
                return self._with_axiom('AlmostComplex')

        class Differentiable(CategoryWithAxiom):
            """
            The category of differentiable manifolds over `\RR`.

            A `d`-dimensional differentiable manifold is a manifold whose
            underlying vector space is `\RR^d` and differentiable atlas.
            """

        class Smooth(CategoryWithAxiom):
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
                    sage: Manifolds().Real().Smooth().super_categories() # indirect doctest
                    [Category of differentiable real manifolds]
                """
                return [Manifolds().Real().Differentiable()]

        class Analytic(CategoryWithAxiom):
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
                    sage: Manifolds().Real().Analytic().super_categories() # indirect doctest
                    [Category of smooth real manifolds]
                """
                return [Manifolds().Real().Smooth()]

        class AlmostComplex(CategoryWithAxiom):
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
                    sage: Manifolds().Real().Analytic().super_categories() # indirect doctest
                    [Category of smooth real manifolds]
                """
                return [Manifolds().Real().Smooth()]

    class Complex(CategoryWithAxiom):
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
                    sage: Manifolds().Complex().Real()
                    Traceback (most recent call last):
                    ...
                    TypeError: a complex manifold is not a real manifold
                """
                raise TypeError("a complex manifold is not a real manifold")

