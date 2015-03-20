"""
Base class for groups
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random

from sage.structure.parent cimport Parent
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ

def is_Group(x):
    """
    Return whether ``x`` is a group object.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: F.<a,b> = FreeGroup()
        sage: from sage.groups.group import is_Group
        sage: is_Group(F)
        True
        sage: is_Group("a string")
        False
    """
    from sage.groups.old import Group as OldGroup
    return isinstance(x, (Group, OldGroup))


cdef class Group(Parent):
    """
    Base class for all groups

    TESTS::

        sage: from sage.groups.group import Group
        sage: G = Group()
        sage: TestSuite(G).run(skip = ["_test_an_element",\
                                       "_test_associativity",\
                                       "_test_elements",\
                                       "_test_elements_eq_reflexive",\
                                       "_test_elements_eq_symmetric",\
                                       "_test_elements_eq_transitive",\
                                       "_test_elements_neq",\
                                       "_test_inverse",\
                                       "_test_one",\
                                       "_test_pickling",\
                                       "_test_prod",\
                                       "_test_some_elements"])
    """
    def __init__(self, base=None, gens=None, category=None):
        """
        The Python constructor

        TESTS::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.category()
            Category of groups
            sage: G = Group(category=Groups()) # todo: do the same test with some subcategory of Groups when there will exist one
            sage: G.category()
            Category of groups
            sage: G = Group(category = CommutativeAdditiveGroups())
            Traceback (most recent call last):
            ...
            ValueError: (Category of commutative additive groups,) is not a subcategory of Category of groups
            sage: G._repr_option('element_is_atomic')
            False

        Check for #8119::

            sage: G = SymmetricGroup(2)
            sage: h = hash(G)
            sage: G.rename('S2')
            sage: h == hash(G)
            True
        """
        from sage.categories.groups import Groups
        if category is None:
            category = Groups()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Groups()) for cat in category):
                raise ValueError("%s is not a subcategory of %s"%(category, Groups()))
        Parent.__init__(self, base=base, gens=gens, category=category)

    def __contains__(self, x):
        r"""
        Test whether `x` defines a group element.

        INPUT:

        - ``x`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: 4 in G               #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        try:
            self(x)
        except TypeError:
            return False
        return True

    def is_abelian(self):
        """
        Test whether this group is abelian.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_abelian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_commutative(self):
        r"""
        Test whether this group is commutative.

        This is an alias for is_abelian, largely to make groups work
        well with the Factorization class.

        (Note for developers: Derived classes should override is_abelian, not
        is_commutative.)

        EXAMPLE::

            sage: SL(2, 7).is_commutative()
            False
        """
        return self.is_abelian()

    def order(self):
        """
        Returns the number of elements of this group, which is either a
        positive integer or infinity.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.order()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_finite(self):
        """
        Returns True if this group is finite.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.order() != infinity

    def is_multiplicative(self):
        """
        Returns True if the group operation is given by \* (rather than
        +).

        Override for additive groups.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_multiplicative()
            True
        """
        return True

    def _an_element_(self):
        """
        Return an element

        OUTPUT:

        An element of the group.

        EXAMPLES:

            sage: G = AbelianGroup([2,3,4,5])
            sage: G.an_element()
            f0*f1*f2*f3
        """
        from sage.misc.all import prod
        return prod(self.gens())

    def random_element(self, bound=None):
        """
        Return a random element of this group.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def quotient(self, H):
        """
        Return the quotient of this group by the normal subgroup
        `H`.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.quotient(G)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

cdef class AbelianGroup(Group):
    """
    Generic abelian group.
    """
    def is_abelian(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.groups.group import AbelianGroup
            sage: G = AbelianGroup()
            sage: G.is_abelian()
            True
        """
        return True

cdef class FiniteGroup(Group):
    """
    Generic finite group.
    """

    def __init__(self, base=None, gens=None, category=None):
        """
        The Python constructor

        TESTS::

            sage: from sage.groups.group import FiniteGroup
            sage: G = FiniteGroup()
            sage: G.category()
            Category of finite groups
        """
        from sage.categories.finite_groups import FiniteGroups
        if category is None:
            category = FiniteGroups()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(FiniteGroups()) for cat in category):
                raise ValueError("%s is not a subcategory of %s"%(category, FiniteGroups()))
        Parent.__init__(self, base=base, gens=gens, category=category)

    def is_finite(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.groups.group import FiniteGroup
            sage: G = FiniteGroup()
            sage: G.is_finite()
            True
        """
        return True

    def cayley_graph(self, connecting_set=None):
        """
        Return the Cayley graph for this finite group.

        INPUT:

        - ``connecting_set`` -- (optional) list of elements to use for
          edges, default is the stored generators

        OUTPUT:

        The Cayley graph as a Sage DiGraph object. To plot the graph
        with with different colors

        EXAMPLES::

            sage: D4 = DihedralGroup(4); D4
            Dihedral group of order 8 as a permutation group
            sage: G = D4.cayley_graph()
            sage: show(G, color_by_label=True, edge_labels=True)
            sage: A5 = AlternatingGroup(5); A5
            Alternating group of order 5!/2 as a permutation group
            sage: G = A5.cayley_graph()
            sage: G.show3d(color_by_label=True, edge_size=0.01, edge_size2=0.02, vertex_size=0.03)
            sage: G.show3d(vertex_size=0.03, edge_size=0.01, edge_size2=0.02, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, xres=700, yres=700, iterations=200) # long time (less than a minute)
            sage: G.num_edges()
            120
            sage: G = A5.cayley_graph(connecting_set=[A5.gens()[0]])
            sage: G.num_edges()
            60
            sage: g=PermutationGroup([(i+1,j+1) for i in range(5) for j in range(5) if j!=i])
            sage: g.cayley_graph(connecting_set=[(1,2),(2,3)])
            Digraph on 120 vertices

            sage: s1 = SymmetricGroup(1); s = s1.cayley_graph(); s.vertices()
            [()]

        AUTHORS:

        - Bobby Moretti (2007-08-10)
        - Robert Miller (2008-05-01): editing
        """
        if connecting_set is None:
            connecting_set = self.gens()
        else:
            if any(g not in self for g in connecting_set):
                raise ValueError("Each element of the connecting set must be in the group!")
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

