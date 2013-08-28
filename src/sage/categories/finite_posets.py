r"""
Finite posets

Here is some terminology used in this file:

- An *order filter* (or *upper set*) is a subset `S` such that if `x \leq y` and `x\in S` then `y\in S`.

- An *order ideal* (or *lower set*) is a subset `S` such that if `x \leq y` and `y\in S` then `x\in S`.
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.posets import Posets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class FinitePosets(Category):
    r"""
    The category of finite posets i.e. finite sets with a partial
    order structure.

    EXAMPLES::

        sage: FinitePosets()
        Category of finite posets
        sage: FinitePosets().super_categories()
        [Category of posets, Category of finite enumerated sets]
        sage: FinitePosets().example()
        NotImplemented

    .. seealso:: :class:`~sage.categories.posets.Posets`, :func:`Poset`

    TESTS::

        sage: C = FinitePosets()
        sage: TestSuite(C).run()

    """
    @cached_method
    def super_categories(self):
        r"""
        Returns a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: FinitePosets().super_categories()
            [Category of posets, Category of finite enumerated sets]
        """
        return [Posets(), FiniteEnumeratedSets()]

    class ParentMethods:

        ##########################################################################
        # Properties of this poset

        def is_lattice(self):
            r"""
            Returns whether this poset is both a meet and a join semilattice.

            EXAMPLES::

                sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
                sage: P.is_lattice()
                True

                sage: P = Poset([[1,2],[3],[3],[]])
                sage: P.is_lattice()
                True

                sage: P = Poset({0:[2,3],1:[2,3]})
                sage: P.is_lattice()
                False
            """
            return self.is_meet_semilattice() and self.is_join_semilattice()

        def is_selfdual(self):
            r"""
            Returns whether this poset is *self-dual*, that is
            isomorphic to its dual poset.

            EXAMPLE::

                sage: P = Poset(([1,2,3],[[1,3],[2,3]]),cover_relations=True)
                sage: P.is_selfdual()
                False

                sage: P = Poset(([1,2,3,4],[[1,3],[1,4],[2,3],[2,4]]),cover_relations=True)
                sage: P.is_selfdual()
                True

            """
            return self.is_isomorphic( self.dual() )


        ##########################################################################
        # Properties of morphisms

        def is_poset_isomorphism(self, f, codomain):
            r"""
            INPUT:

            - ``f`` -- a function from ``self`` to ``codomain``
            - ``codomain`` -- a poset

            Returns whether `f` is an isomorphism of posets form
            ``self`` to ``codomain``.

            EXAMPLES:

            We build the poset `D` of divisors of 30, and check that
            it is isomorphic to the boolean lattice `B` of the subsets
            of `\{2,3,5\}` ordered by inclusion, via the reverse
            function `f: B \mapsto D, b \rightarrow \prod_{x\in b} x`::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5])], attrcall("issubset")))
                sage: def f(b): return D(prod(b))
                sage: B.is_poset_isomorphism(f, D)
                True

            On the other hand, `f` is not an isomorphism to the chain
            of divisors of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(prod(b))
                sage: B.is_poset_isomorphism(f, P)
                False

            A non surjective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3])], attrcall("issubset")))
                sage: def f(b): return D(prod(b))
                sage: B.is_poset_isomorphism(f, D)
                False

            A non injective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5,6])], attrcall("issubset")))
                sage: def f(b): return D(gcd(prod(b), 30))
                sage: B.is_poset_isomorphism(f, D)
                False

            .. note:: since ``D`` and ``B`` are not facade posets, ``f`` is
               responsible for the conversions between integers and subsets to
               elements of ``D`` and ``B`` and back.

            .. seealso:: :meth:`FiniteLatticePosets.ParentMethods.is_lattice_morphism`
            """
            image = set(f(x) for x in self)
            if len(image) != self.cardinality():
                # Not injective
                return False
            if len(image) != codomain.cardinality():
                # Not surjective
                return False
            for x in self:
                if set(f(y) for y in self.upper_covers(x)) != set(codomain.upper_covers(f(x))):
                    return False
            return True

        def is_poset_morphism(self, f, codomain):
            r"""
            INPUT:

            - ``f`` -- a function from ``self`` to ``codomain``
            - ``codomain`` -- a poset

            Returns whether `f` is a morphism of posets form ``self``
            to ``codomain``, that is

            .. math::  x\leq y \Longrightarrow f(x) \leq f(y)

            EXAMPLES:

            We build the boolean lattice of the subsets of
            `\{2,3,5,6\}` and the lattice of divisors of `30`, and
            check that the map `b \mapsto \gcd(\prod_{x\in b}, 30)` is a
            morphism of posets::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5,6])], attrcall("issubset")))
                sage: def f(b): return D(gcd(prod(b), 30))
                sage: B.is_poset_morphism(f, D)
                True

            .. note:: since ``D`` and ``B`` are not facade posets, ``f`` is responsible
               for the conversions between integers and subsets to elements of
               ``D`` and ``B`` and back.

            `f` is also a morphism of posets to the chain of divisors
            of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(gcd(prod(b), 30))
                sage: B.is_poset_morphism(f, P)
                True

            FIXME: should this be ``is_order_preserving_morphism``?

            .. seealso:: :meth:`is_poset_isomorphism`

            TESTS:

            Base cases::

                sage: P = Posets.ChainPoset(2)
                sage: Q = Posets.AntichainPoset(2)
                sage: f = lambda x: 1-x
                sage: P.is_poset_morphism(f, P)
                False
                sage: P.is_poset_morphism(f, Q)
                False
                sage: Q.is_poset_morphism(f, Q)
                True
                sage: Q.is_poset_morphism(f, P)
                True

                sage: P = Poset(); P
                Finite poset containing 0 elements
                sage: P.is_poset_morphism(f, P)
                True

            """
            for x in self:
                for y in self.upper_covers(x):
                    if not codomain.is_lequal(f(x),f(y)):
                        return False
            return True

        ##########################################################################
        # About order ideals, order filters and the like

        def order_ideal_generators(self, ideal, direction='down'):
            r"""
            Generators for an order ideal (resp. an order filter)

            INPUT:

            - ``ideal`` -- an order ideal `I` of ``self``, as a list (or iterable)
            - ``direction`` -- 'up' or 'down' (default: 'down')

            Returns the minimal set of generators for the ideal
            `I`. It forms an antichain.

            EXAMPLES:

            We build the boolean lattice of all subsets of `\{1,2,3\}`
            ordered by inclusion, and compute an order ideal there::

                sage: P = Poset((Subsets([1,2,3]), attrcall("issubset")))
                sage: I = P.order_ideal([Set([1,2]), Set([2,3]), Set([1])]); I
                [{}, {1}, {3}, {2}, {1, 2}, {2, 3}]

            Then, we retrieve the generators of this ideal::

                sage: P.order_ideal_generators(I)
                {{1, 2}, {2, 3}}

            If ``direction`` is 'up', then this instead computes
            the minimal generators for an order filter::

                sage: I = P.order_filter([Set([1,2]), Set([2,3]), Set([1])]); I
                [{1}, {1, 3}, {1, 2}, {2, 3}, {1, 2, 3}]
                sage: P.order_ideal_generators(I, direction='up')
                {{2, 3}, {1}}

            Complexity: `O(n+m)` where `n` is the cardinality of `I`,
            and `m` the number of upper covers of elements of `I`.
            """
            if direction == 'down':
                covers = self.upper_covers
            else:
                covers = self.lower_covers
            ideal_as_set = set(ideal)
            from sage.sets.set import Set
            return Set(x for x in ideal if all(y not in ideal_as_set
                                               for y in covers(x)))

        def order_filter_generators(self, filter):
            r"""
            Generators for an order filter

            INPUT:

            - ``filter`` -- an order filter of ``self``, as a list (or iterable)

            EXAMPLES::

                sage: P = Poset((Subsets([1,2,3]), attrcall("issubset")))
                sage: I = P.order_filter([Set([1,2]), Set([2,3]), Set([1])]); I
                [{1}, {1, 3}, {1, 2}, {2, 3}, {1, 2, 3}]
                sage: P.order_filter_generators(I)
                {{2, 3}, {1}}

            .. seealso:: :meth:`order_ideal_generators`
            """
            return self.order_ideal_generators(filter, direction='up')

        def order_ideal_complement_generators(self, antichain, direction='up'):
            r"""
            The generators of the complement of an order ideal (resp. order filter)

            INPUT:

            - ``antichain`` -- an antichain of ``self``, as a list (or
              iterable), or, more generally, generators of an order ideal
              (resp. order filter)
            - ``direction`` -- 'up' or 'down' (default: 'up')

            OUTPUT:

            - the complement order filter (resp. order ideal), represented
              by its generating antichain

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: P.order_ideal_complement_generators([1])
                set([2])
                sage: P.order_ideal_complement_generators([3])
                set([])
                sage: P.order_ideal_complement_generators([1,2])
                set([3])
                sage: P.order_ideal_complement_generators([1,2,3])
                set([])

                sage: P.order_ideal_complement_generators([1], direction="down")
                set([2])
                sage: P.order_ideal_complement_generators([3], direction="down")
                set([1, 2])
                sage: P.order_ideal_complement_generators([1,2], direction="down")
                set([])
                sage: P.order_ideal_complement_generators([1,2,3], direction="down")
                set([])

            .. WARNING::

                This is a brute force implementation, building the
                order ideal generated by the antichain, and searching
                for order filter generators of its complement
            """
            if direction == 'up':
                I = self.order_ideal(antichain)
            else:
                I = self.order_filter(antichain)
            I_comp = set(self).difference(I)
            return set(self.order_ideal_generators(I_comp, direction = direction))

        panyushev_complement = order_ideal_complement_generators

        def rowmotion(self, order_ideal):
            r"""
            The image of the order ideal ``order_ideal`` under rowmotion
            in ``self``.

            INPUT:

            - ``order_ideal`` -- an order ideal of ``self``, as a set

            OUTPUT:

            - the image of ``order_ideal`` under rowmotion, as a set again

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [8], 5: [], 6: [5], 7: [1, 4], 8: []} )
                sage: I = Set({2, 6, 1, 7})
                sage: P.rowmotion(I)
                {1, 3, 4, 5, 6, 7}

                sage: P = Poset( {} )
                sage: I = Set({})
                sage: P.rowmotion(I)
                Set of elements of {}
            """
            result = order_ideal
            for i in reversed(self.linear_extension()):
                result = self.order_ideal_toggle(result, i)
            return result

        def panyushev_orbits(self, element_constructor = set):
            r"""
            Return the Panyushev orbits of antichains in ``self``.

            INPUT:

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the antichains
              before they are returned.

            OUTPUT:

            - the partition of the set of all antichains of ``self`` into
              orbits under Panyushev complementation. This is returned as
              a list of lists ``L`` such that for each ``L`` and ``i``,
              cyclically:
              ``self.order_ideal_complement_generators(L[i]) == L[i+1]``.
              The entries ``L[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: P.panyushev_orbits()
                [[set([2]), set([1])], [set([]), set([1, 2]), set([3])]]
                sage: P.panyushev_orbits(element_constructor=list)
                [[[2], [1]], [[], [1, 2], [3]]]
                sage: P.panyushev_orbits(element_constructor=frozenset)
                [[frozenset([2]), frozenset([1])],
                 [frozenset([]), frozenset([1, 2]), frozenset([3])]]
                sage: P.panyushev_orbits(element_constructor=tuple)
                [[(2,), (1,)], [(), (1, 2), (3,)]]
                sage: P = Poset( {} )
                sage: P.panyushev_orbits()
                [[set([])]]
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning the orbits.
            AC = set(self.antichains(element_constructor = frozenset))
            orbits = []
            while AC:
                A = AC.pop()
                orbit = [ A ]
                while True:
                    A = frozenset(self.order_ideal_complement_generators(A))
                    if A not in AC: break
                    orbit.append( A )
                    AC.remove( A )
                orbits.append(map(element_constructor, orbit))
            return orbits

        def rowmotion_orbits(self, element_constructor = set):
            r"""
            Return the rowmotion orbits of order ideals in ``self``.

            INPUT:

            - ``element_constructor`` (defaults to ``set``) -- a type
              constructor (``set``, ``tuple``, ``list``, ``frozenset``,
              ``iter``, etc.) which is to be applied to the antichains
              before they are returned.

            OUTPUT:

            - the partition of the set of all order ideals of ``self``
              into orbits under rowmotion. This is returned as
              a list of lists ``L`` such that for each ``L`` and ``i``,
              cyclically: ``self.rowmotion(L[i]) == L[i+1]``.
              The entries ``L[i]`` are sets by default, but depending on
              the optional keyword variable ``element_constructors``
              they can also be tuples, lists etc.

            EXAMPLES::

                sage: P = Poset( {1: [2, 3], 2: [], 3: [], 4: [2]} )
                sage: sorted(len(o) for o in P.rowmotion_orbits())
                [3, 5]
                sage: sorted(P.rowmotion_orbits(element_constructor=list))
                [[[1, 3], [4], [1], [4, 1, 3], [4, 1, 2]], [[4, 1], [4, 1, 2, 3], []]]
                sage: sorted(P.rowmotion_orbits(element_constructor=tuple))
                [[(1, 3), (4,), (1,), (4, 1, 3), (4, 1, 2)], [(4, 1), (4, 1, 2, 3), ()]]
                sage: P = Poset({})
                sage: sorted(P.rowmotion_orbits(element_constructor=tuple))
                [[()]]
            """
            pan_orbits = self.panyushev_orbits(element_constructor = list)
            return [[element_constructor(self.order_ideal(oideal)) for oideal in orbit] for orbit in pan_orbits]

        def toggling_orbits(self, vs, element_constructor = set):
            r"""
            Returns the orbits of order ideals in ``self`` under the
            sequence of toggles given by ``vs``.

            INPUT:

            - ``vs``: a list (or other iterable) of elements of ``self``
              (but since the output depends on the order, sets should
              not be used as ``vs``).

            OUTPUT:

            - a partition of the order ideals of ``self``, as a list of
              sets ``L`` such that for each ``L`` and ``i``, cyclically:
              ``self.order_ideal_toggles(L[i], vs) == L[i+1]``.

            EXAMPLES::

                sage: P = Poset( {1: [2, 4], 2: [], 3: [4], 4: []} )
                sage: sorted(len(o) for o in P.toggling_orbits([1, 2]))
                [2, 3, 3]
                sage: P = Poset( {1: [3], 2: [1, 4], 3: [], 4: [3]} )
                sage: sorted(len(o) for o in P.toggling_orbits((1, 2, 4, 3)))
                [3, 3]
            """
            # TODO: implement a generic function taking a set and
            # bijections on this set, and returning the orbits.
            OI = set(self.order_ideals_lattice())
            orbits = []
            while OI:
                A = OI.pop()
                orbit = [ A ]
                while True:
                    A = self.order_ideal_toggles(A, vs)
                    if A not in OI: break
                    orbit.append( A )
                    OI.remove( A )
                orbits.append(map(element_constructor, orbit))
            return orbits

        def order_ideals_lattice(self, as_ideals=True):
            r"""
            Returns the lattice of order ideals of a poset `P`,
            ordered by inclusion. The usual notation is `J(P)`.

            The underlying set is by default the set of order ideals
            of `P`. It can be alternatively chosen to be the set of
            antichains of `P`.

            INPUT:

            - ``as_ideals`` -- Boolean, if ``True`` (default) returns
              a poset on the set of order ideals, otherwise on the set
              of antichains

            EXAMPLES::

                sage: P = Posets.PentagonPoset(facade = True)
                sage: P.cover_relations()
                [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]
                sage: J = P.order_ideals_lattice(); J
                Finite lattice containing 8 elements
                sage: list(J)
                [{}, {0}, {0, 2}, {0, 1}, {0, 1, 2}, {0, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3, 4}]

            As a lattice on antichains::

                sage: J2 = P.order_ideals_lattice(False); J2
                Finite lattice containing 8 elements
                sage: list(J2)
                [(0,), (1, 2), (1, 3), (1,), (2,), (3,), (4,), ()]

            TESTS::

                sage: J = Posets.DiamondPoset(4, facade = True).order_ideals_lattice(); J
                Finite lattice containing 6 elements
                sage: list(J)
                [{}, {0}, {0, 2}, {0, 1}, {0, 1, 2}, {0, 1, 2, 3}]
                sage: J.cover_relations()
                [[{}, {0}], [{0}, {0, 2}], [{0}, {0, 1}], [{0, 2}, {0, 1, 2}], [{0, 1}, {0, 1, 2}], [{0, 1, 2}, {0, 1, 2, 3}]]

            .. NOTE:: we use facade posets in the examples above just
               to ensure a nicer ordering in the output.
            """
            from sage.combinat.posets.lattices import LatticePoset
            if as_ideals:
                from sage.misc.misc import attrcall
                from sage.sets.set import Set
                ideals = [Set( self.order_ideal(antichain) ) for antichain in self.antichains()]
                return LatticePoset((ideals,attrcall("issubset")))
            else:
                from sage.misc.cachefunc import cached_function
                antichains = [tuple(a) for a in self.antichains()]
                @cached_function
                def is_above(a,xb):
                    return any(self.is_lequal(xa,xb) for xa in a)
                def cmp(a,b):
                    return all(is_above(a,xb) for xb in b)
                return LatticePoset((antichains,cmp))

        @abstract_method(optional = True)
        def antichains(self):
            r"""
            Returns all antichains of ``self``.

            EXAMPLES::

                sage: A = Posets.PentagonPoset().antichains(); A
                Set of antichains of Finite lattice containing 5 elements
                sage: list(A)
                [[], [0], [1], [1, 2], [1, 3], [2], [3], [4]]
            """

        def directed_subsets(self, direction):
            r"""
            Return the order filters (resp. order ideals) of ``self``, as lists.

            If ``direction`` is 'up', returns the order filters (upper sets).

            If ``direction`` is 'down', returns the order ideals (lower sets).

            INPUT:

            - ``direction`` -- 'up' or 'down'

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")), facade=True)
                sage: A = P.directed_subsets('up')
                sage: sorted(list(A))
                [[], [1, 2, 4, 3, 6, 12], [2, 4, 3, 6, 12], [2, 4, 6, 12], [3, 6, 12], [4, 3, 6, 12], [4, 6, 12], [4, 12], [6, 12], [12]]

            TESTS::

                sage: list(Poset().directed_subsets('up'))
                [[]]
            """
            return self.antichains().map(lambda elements: self.directed_subset(elements, direction))
