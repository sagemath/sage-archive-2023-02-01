"""
Graphs
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.simplicial_complexes import SimplicialComplexes

class Graphs(Category_singleton):
    r"""
    The category of graphs.

    EXAMPLES::

        sage: from sage.categories.graphs import Graphs
        sage: C = Graphs(); C
        Category of graphs

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.graphs import Graphs
            sage: Graphs().super_categories()
            [Category of simplicial complexes]
        """
        return [SimplicialComplexes()]

    class ParentMethods:
        def dimension(self):
            """
            Return the dimension of ``self`` as a CW complex.
            """
            if self.edges():
                return 1
            return 0

        def facets(self):
            """
            Return the facets of ``self``.
            """
            return self.edges()

        def faces(self):
            """
            Return the faces of ``self``.
            """
            return set(self.edges()).union(self.vertices())

