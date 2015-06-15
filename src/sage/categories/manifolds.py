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
        def Differentiable(self):
            """
            Return the subcategory of the differentiable objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds().Differentiable()
                Category of differentiable manifolds

            TESTS::

                sage: TestSuite(Manifolds().Differentiable()).run()
                sage: Manifolds().Differentiable.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Differentiable')

        @cached_method
        def Smooth(self):
            """
            Return the subcategory of the smooth objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds().Smooth()
                Category of smooth manifolds

            TESTS::

                sage: TestSuite(Manifolds().Smooth()).run()
                sage: Manifolds().Smooth.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Smooth')

        @cached_method
        def Complex(self):
            """
            Return the full subcategory of the complex objects of ``self``.

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

    class FiniteDimensional(CategoryWithAxiom):
        """
        Category of finite dimensional manifolds.
        """

    class Connected(CategoryWithAxiom):
        """
        The category of connected manifolds.
        """

    class Differentiable(CategoryWithAxiom):
        """
        The category of differentiable manifolds.
        """

    class Smooth(CategoryWithAxiom):
        """
        The category of smooth manifolds.
        """
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            A smooth manifold is differentiable.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds().Smooth().super_categories() # indirect doctest
                [Category of differentiable manifolds]
            """
            return [Manifolds().Differentiable()]

    class Complex(CategoryWithAxiom):
        r"""
        The category of complex manifolds, i.e., the underlying vector
        space is `\CC^d`.
        """

