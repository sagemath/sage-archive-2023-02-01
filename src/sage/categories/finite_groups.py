r"""
FiniteGroups
"""
#*****************************************************************************
#  Copyright (C) 2010-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom

class FiniteGroups(CategoryWithAxiom):
    r"""
    The category of finite (multiplicative) groups.

    EXAMPLES::

        sage: C = FiniteGroups(); C
        Category of finite groups
        sage: C.super_categories()
        [Category of finite monoids, Category of groups]
        sage: C.example()
        General Linear Group of degree 2 over Finite Field of size 3

    TESTS::

        sage: TestSuite(C).run()
    """

    def example(self):
        """
        Return an example of finite group, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: G = FiniteGroups().example(); G
            General Linear Group of degree 2 over Finite Field of size 3
        """
        from sage.groups.matrix_gps.linear import GL
        return GL(2,3)

    class ParentMethods:

        def semigroup_generators(self):
            """
            Returns semigroup generators for self.

            For finite groups, the group generators are also semigroup
            generators. Hence, this default implementation calls
            :meth:`~sage.categories.groups.Groups.ParentMethods.group_generators`.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.semigroup_generators()
                Family ((2,3,4), (1,2,3))
            """
            return self.group_generators()

        def cardinality(self):
            """
            Returns the cardinality of ``self``, as per
            :meth:`EnumeratedSets.ParentMethods.cardinality`.

            This default implementation calls :meth:`.order` if
            available, and otherwise resorts to
            :meth:`._cardinality_from_iterator`. This is for backward
            compatibility only. Finite groups should override this
            method instead of :meth:`.order`.

            EXAMPLES:

            We need to use a finite group which uses this default
            implementation of cardinality::

                sage: R.<x> = PolynomialRing(QQ)
                sage: f = x^4 - 17*x^3 - 2*x + 1
                sage: G = f.galois_group(pari_group=True); G
                PARI group [24, -1, 5, "S4"] of degree 4
                sage: G.cardinality.__module__
                'sage.categories.finite_groups'
                sage: G.cardinality()
                24
            """
            try:
                o = self.order
            except AttributeError:
                return self._cardinality_from_iterator()
            else:
                return o()

        def some_elements(self):
            """
            Return some elements of ``self``.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.some_elements()
                [(2,3,4), (1,2,3)]
            """
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

        def conjugacy_class(self, g):
            r"""
            Return the conjugacy class of the element ``g``.

            This is a fall-back method for groups not defined over GAP.

            EXAMPLES::

                sage: W = WeylGroup(['C',6])
                sage: c = W.conjugacy_class(W.an_element())
                sage: type(c)
                <class 'sage.groups.conjugacy_classes.ConjugacyClass_with_category'>
            """
            from sage.groups.conjugacy_classes import ConjugacyClass
            return ConjugacyClass(self, g)

        def conjugacy_classes(self):
            r"""
            Return a list with all the conjugacy classes of the group.

            This will eventually be a fall-back method for groups not defined
            over GAP. Right now just raises a ``NotImplementedError``, until
            we include a non-GAP way of listing the conjugacy classes
            representatives.

            EXAMPLES::

                sage: from sage.groups.group import FiniteGroup
                sage: G = FiniteGroup()
                sage: G.conjugacy_classes()
                Traceback (most recent call last):
                ...
                NotImplementedError: Listing the conjugacy classes for
                group <type 'sage.groups.group.FiniteGroup'> is not implemented
            """
            raise NotImplementedError("Listing the conjugacy classes for group %s is not implemented"%self)

        def conjugacy_classes_representatives(self):
            r"""
            Return a list of the conjugacy classes representatives of the group.

            EXAMPLES::

                sage: G = SymmetricGroup(3)
                sage: G.conjugacy_classes_representatives()
                [(), (1,2), (1,2,3)]
           """
            return [C.representative() for C in self.conjugacy_classes()]

    class ElementMethods:
        def conjugacy_class(self):
            r"""
            Return the conjugacy class of ``self``.

            EXAMPLES::

                sage: H = MatrixGroup([matrix(GF(5),2,[1,2, -1, 1]), matrix(GF(5),2, [1,1, 0,1])])
                sage: h = H(matrix(GF(5),2,[1,2, -1, 1]))
                sage: h.conjugacy_class()
                Conjugacy class of [1 2]
                [4 1] in Matrix group over Finite Field of size 5 with 2 generators (
                [1 2]  [1 1]
                [4 1], [0 1]
                )
            """
            return self.parent().conjugacy_class(self)

