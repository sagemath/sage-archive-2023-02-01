r"""
FiniteGroups
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.finite_monoids import FiniteMonoids
from sage.categories.groups import Groups

class FiniteGroups(Category):
    r"""
    The category of (multiplicative) finite groups.

    EXAMPLES::

        sage: FiniteGroups()
        Category of finite groups
        sage: FiniteGroups().super_categories()
        [Category of groups, Category of finite monoids]
        sage: FiniteGroups().example()
        General Linear Group of degree 2 over Finite Field of size 3

    TESTS::

        sage: C = FiniteGroups()
        sage: TestSuite(C).run(verbose = True)
        running ._test_category() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass

        sage: G = FiniteGroups().example()
        sage: G.category()
        Category of finite groups
        sage: TestSuite(G).run(skip = "_test_enumerated_set_iter_list")

    TODO: define an iterator compatible with GAP's list method
    Does GAP provide such an iterator?
    """
    @cached_method
    def super_categories(self):
        r"""
        Returns a list of the (immediate) super categories of ``self``.

        EXAMPLES::

            sage: FiniteGroups().super_categories()
            [Category of groups, Category of finite monoids]

        """
        return [Groups(), FiniteMonoids()]

    def example(self):
        """
        Returns an example of finite group, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: G = FiniteGroups().example(); G
            General Linear Group of degree 2 over Finite Field of size 3
        """
        from sage.groups.matrix_gps.general_linear import GL
        return GL(2,3)

    class ParentMethods:

        def semigroup_generators(self):
            """
            Returns semigroup generators for self.

            For finite groups, the group generators are also semigroup
            generators. Hence, this default implementation calls
            :meth:`.group_generators``.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.semigroup_generators()
                [(1,2,3), (2,3,4)]
            """
            return self.group_generators()

        def cardinality(self):
            """
            Returns the cardinality of self, as per :meth:`EnumeratedSets.ParentMethods.cardinality`.

            This default implementation calls :meth:`.order` if
            available, and otherwise resorts to
            meth:`._cardinality_from_iterator`. This is for backward
            compatibility only. Finite groups should override this
            method instead of :meth:`.order`.

            EXAMPLES::

                sage: G = DihedralGroup(6)
                sage: G.cardinality.__module__
                'sage.categories.finite_groups'
                sage: G.cardinality()
                12
            """
            try:
                o = self.order
            except AttributeError:
                return self._cardinality_from_iterator()
            else:
                return o()

        def some_elements(self):
            return self.gens()

        # TODO: merge with that of finite semigroups
        def cayley_graph_disabled(self, connecting_set=None):
            """

            AUTHORS:

            - Bobby Moretti (2007-08-10)

            - Robert Miller (2008-05-01): editing
            """
            if connecting_set is None:
                connecting_set = self.gens()
            else:
                try:
                    for g in connecting_set:
                        assert g in self
                except AssertionError:
                    raise RuntimeError("Each element of the connecting set must be in the group!")
                connecting_set = [self(g) for g in connecting_set]
            from sage.graphs.all import DiGraph
            arrows = {}
            for x in self:
                arrows[x] = {}
                for g in connecting_set:
                    xg = x*g # cache the multiplication
                    if not xg == x:
                        arrows[x][xg] = g

            return DiGraph(arrows, implementation='networkx')

    class ElementMethods:
        pass
