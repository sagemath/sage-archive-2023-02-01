r"""
Finite Coxeter Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#  Copyright (C) 2009    Nicolas Borie <nicolas dot borie at math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category import Category
from sage.categories.coxeter_groups import CoxeterGroups
from sage.categories.finite_groups import FiniteGroups
from sage.combinat.backtrack import search_forest_iterator

class FiniteCoxeterGroups(Category):
    r"""
    The category of finite Coxeter groups.

    EXAMPLES::

        sage: FiniteSemigroups()
        Category of finite semigroups
        sage: FiniteSemigroups().super_categories()
        [Category of semigroups, Category of finite enumerated sets]

        sage: G = FiniteCoxeterGroups().example()
        sage: G.cayley_graph(side = "right").plot()

    Here are some further examples::

        sage: FiniteWeylGroups().example()
        The symmetric group on {0, ..., 3}

        sage: WeylGroup(["B", 3])
        Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)

    Those other examples will eventually be also in this category::

        sage: SymmetricGroup(4)
        Symmetric group of order 4! as a permutation group
        sage: DihedralGroup(5)
        Dihedral group of order 10 as a permutation group
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: FiniteCoxeterGroups().super_categories()
            [Category of coxeter groups, Category of finite groups]
        """
        return [CoxeterGroups(), FiniteGroups()]

    class ParentMethods:
        def __iter__(self):
            r"""
            Returns an iterator over the elements of this finite Coxeter group.

            EXAMPLES::

                sage: D5 = FiniteCoxeterGroups().example(5)
                sage: sorted(list(D5)) # indirect doctest (but see :meth:`._test_enumerated_set_iter_list`)
                [(),
                 (1,),
                 (1, 2),
                 (1, 2, 1),
                 (1, 2, 1, 2),
                 (1, 2, 1, 2, 1),
                 (2,),
                 (2, 1),
                 (2, 1, 2),
                 (2, 1, 2, 1)]
            """
            def succ(u):
                for i in u.descents(positive = True, side = "right"):
                    u1 = u.apply_simple_reflection(i, "right")
                    if i == u1.first_descent():
                        yield u1
                return
            return search_forest_iterator((self.one(),), succ)

        @lazy_attribute
        def w0(self):
            r"""
            Return the longest element of self.

            This attribute is deprecated.

            EXAMPLES::

                sage: D8 = FiniteCoxeterGroups().example(8)
                sage: D8.w0
                (1, 2, 1, 2, 1, 2, 1, 2)
                sage: D3 = FiniteCoxeterGroups().example(3)
                sage: D3.w0
                (1, 2, 1)
            """
            return self.long_element()

        def long_element(self, index_set = None):
            r"""

            INPUT:

            - ``index_set`` - a subset (as a list or iterable) of the
              nodes of the dynkin diagram; (default: all of them)

            Returns the longest element of ``self``, or of the
            parabolic subgroup corresponding to the given index_set.

            Should this method be called maximal_element? longest_element?

            EXAMPLES::

                sage: D10 = FiniteCoxeterGroups().example(10)
                sage: D10.long_element()
                (1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
                sage: D10.long_element([1])
                (1,)
                sage: D10.long_element([2])
                (2,)
                sage: D10.long_element([])
                ()

                sage: D7 = FiniteCoxeterGroups().example(7)
                sage: D7.long_element()
                (1, 2, 1, 2, 1, 2, 1)

            """
            if index_set is None:
                index_set = self.index_set()
            w = self.one()
            while True:
                i = w.first_descent(index_set = index_set, positive = True)
                if i is None:
                    return w
                else:
                    w = w.apply_simple_reflection(i)

        @cached_method
        def bruhat_poset(self):
            """
            Returns the Bruhat poset of ``self``

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.bruhat_poset()
                sage: P
                Finite poset containing 6 elements
                sage: P.show()

            Here are some typical operations on this poset::

                sage: W = WeylGroup(["A", 3])
                sage: P = W.bruhat_poset()
                sage: u = W.from_reduced_word([3,1])
                sage: v = W.from_reduced_word([3,2,1,2,3])
                sage: P(u) <= P(v)
                True
                sage: len(P.interval(P(u), P(v)))
                10
                sage: P.is_join_semilattice()
                False

            See also :class:`Poset`.

            TESTS::

                sage: [len(WeylGroup(["A", n]).bruhat_poset().cover_relations()) for n in [1,2,3]]
                [1, 8, 58]

            TODO:

            - Use the symmetric group in the examples (for nicer
              output), and print the edges for a stronger test.
            - The constructed poset should be lazy, in order to
              handle large / infinite Coxeter groups.
            """
            from sage.combinat.posets.posets import Poset
            covers = tuple([u, v] for v in self for u in v.bruhat_lower_covers() )
            return Poset((self, covers), cover_relations = True)

        @cached_method
        def weak_poset(self, side = "right"):
            """
            INPUT:

            - ``side`` -- "left" or "right" (default: "right")

            Returns the left (resp. right) poset for weak order.  In
            this poset, `u` is smaller than `v` if some reduced word
            of `u` is a right (resp. left) factor of some reduced word
            of `v`.

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.weak_poset()
                sage: P
                Finite poset containing 6 elements
                sage: P.show()

            This poset is in fact a lattice::

                sage: W = WeylGroup(["B", 3])
                sage: P = W.weak_poset(side = "left")
                sage: P.is_join_semilattice(), P.is_meet_semilattice() # todo: implement is_lattice
                (True, True)

            As a bonus feature, one can create the left-right weak
            poset::

                sage: W = WeylGroup(["A",2])
                sage: P = W.weak_poset(side = "twosided")
                sage: P.show()
                sage: len(P.hasse_diagram().edges())
                8

            This is the transitive closure of the union of left and
            right order. In this poset, `u` is smaller than `v` if
            some reduced word of `u` is a factor of some reduced word
            of `v`.

            TESTS::

                sage: [len(WeylGroup(["A", n]).weak_poset(side = "right").cover_relations()) for n in [1,2,3]]
                [1, 6, 36]
                sage: [len(WeylGroup(["A", n]).weak_poset(side = "left" ).cover_relations()) for n in [1,2,3]]
                [1, 6, 36]

            TODO:

            - Use the symmetric group in the examples (for nicer
              output), and print the edges for a stronger test.
            - The constructed poset should be lazy, in order to
              handle large / infinite Coxeter groups.

            """
            from sage.combinat.posets.posets import Poset
            if side == "twosided":
                covers = tuple([u, v] for u in self for v in u.upper_covers(side="left")+u.upper_covers(side="right") )
            else:
                covers = tuple([u, v] for u in self for v in u.upper_covers(side=side) )
            return Poset((self, covers), cover_relations = True)

    class ElementMethods:

        @cached_in_parent_method
        def bruhat_upper_covers(self):
            r"""
            Returns all the elements that cover ``self`` in Bruhat order.

            EXAMPLES::

                sage: W = WeylGroup(["A",4])
                sage: w = W.from_reduced_word([3,2])
                sage: print([v.reduced_word() for v in w.bruhat_upper_covers()])
                [[4, 3, 2], [3, 4, 2], [2, 3, 2], [3, 1, 2], [3, 2, 1]]

                sage: W = WeylGroup(["B",6])
                sage: w = W.from_reduced_word([1,2,1,4,5])
                sage: C = w.bruhat_upper_covers()
                sage: len(C)
                9
                sage: print([v.reduced_word() for v in C])
                [[6, 4, 5, 1, 2, 1], [4, 5, 6, 1, 2, 1], [3, 4, 5, 1, 2, 1], [2, 3, 4, 5, 1, 2],
                [1, 2, 3, 4, 5, 1], [4, 5, 4, 1, 2, 1], [4, 5, 3, 1, 2, 1], [4, 5, 2, 3, 1, 2],
                [4, 5, 1, 2, 3, 1]]
                sage: ww = W.from_reduced_word([5,6,5])
                sage: CC = ww.bruhat_upper_covers()
                sage: print([v.reduced_word() for v in CC])
                [[6, 5, 6, 5], [4, 5, 6, 5], [5, 6, 4, 5], [5, 6, 5, 4], [5, 6, 5, 3], [5, 6, 5, 2],
                [5, 6, 5, 1]]

            Recursive algorithm: write `w` for ``self``. If `i` is a
            non-descent of `w``, then the covers of `w` are exactly
            `\{ws_i, u_1s_i, u_2s_i,..., u_js_i\}', where the 'u_k'
            are those covers of 'ws_i' that have a descent at `i`.
            """

            i = self.first_descent(positive=True)
            if i is not None:
                wsi = self.apply_simple_reflection(i)
                return [u.apply_simple_reflection(i) for u in wsi.bruhat_upper_covers() if u.has_descent(i)] + [wsi]
            else:
                return []
