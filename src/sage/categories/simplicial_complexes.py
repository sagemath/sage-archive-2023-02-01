"""
Simplicial Complexes
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
#from sage.categories.cw_complexes import CWComplexes
from sage.categories.sets_cat import Sets

class SimplicialComplexes(Category_singleton):
    r"""
    The category of abstract simplicial complexes.

    An abstract simplicial complex `A` is a collection of sets `X`
    such that:

    - `\emptyset \in A`,
    - if `X \subset Y \in A`, then `X \in A`.

    .. TODO::

        Implement the category of simplicial complexes considered
        as :class:`CW complexes <sage.categories.cw_complexes.CWComplexes>`.

    EXAMPLES::

        sage: from sage.categories.simplicial_complexes import SimplicialComplexes
        sage: C = SimplicialComplexes(); C
        Category of simplicial complexes

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.simplicial_complexes import SimplicialComplexes
            sage: SimplicialComplexes().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class Finite(CategoryWithAxiom):
        """
        Category of finite simplicial complexes.
        """
        class ParentMethods:
            @cached_method
            def dimension(self):
                """
                Return the dimension of ``self``.
                """
                return max(c.dimension() for c in self.facets())

    class ParentMethods:
        @abstract_method
        def facets(self):
            """
            Return the facets of ``self``.
            """

        @abstract_method
        def faces(self):
            """
            Return the faces of ``self``.
            """

