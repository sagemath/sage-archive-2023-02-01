r"""
Finite posets
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
                sage: def f(b): return D(prod(b.element))
                sage: B.is_poset_isomorphism(f, D)
                True

            On the other hand, `f` is not an isomorphism to the chain
            of divisors of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(prod(b.element))
                sage: B.is_poset_isomorphism(f, P)
                False

            A non surjective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3])], attrcall("issubset")))
                sage: def f(b): return D(prod(b.element))
                sage: B.is_poset_isomorphism(f, D)
                False

            A non injective case::

                sage: B = Poset(([frozenset(s) for s in Subsets([2,3,5,6])], attrcall("issubset")))
                sage: def f(b): return D(gcd(prod(b.element), 30))
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
                sage: def f(b): return D(gcd(prod(b.element), 30))
                sage: B.is_poset_morphism(f, D)
                True

            .. note:: since ``D`` and ``B`` are not facade posets, ``f`` is responsible
               for the conversions between integers and subsets to elements of
               ``D`` and ``B`` and back.

            `f` is also a morphism of posets to the chain of divisors
            of 30, ordered by usual comparison::

                sage: P = Poset((divisors(30), operator.le))
                sage: def f(b): return P(gcd(prod(b.element), 30))
                sage: B.is_poset_morphism(f, P)
                True

            FIXME: should this be ``is_order_preserving_morphism``?

            .. seealso:: :meth:`is_poset_isomorphism`

            TESTS:

            Base cases::

                sage: P = Posets.ChainPoset(2)
                sage: Q = Posets.AntichainPoset(2)
                sage: f = lambda x: 1-x.element
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
        # About order ideals and the like

        def order_ideal_generators(self, ideal, direction = 'down'):
            r"""
            Generators for an order ideal (resp. filter)

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
            the minimal generators for an upper ideal::

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
            return Set(x for x in ideal if all(y not in ideal_as_set for y in covers(x)))

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
            The generators of the complement of an order ideal.

            INPUT:

            - ``antichain`` -- an antichain of ``self``, as a list (or iterable)
            - ``direction`` -- 'up' or 'down' (default: 'up')

            Returns the minimal set of generators for the complement. It forms
            an antichain.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: P.order_ideal_complement_generators([1])
                set([2])
                sage: P.order_ideal_complement_generators([1,2])
                set([3])
                sage: P.order_ideal_complement_generators([1,2,3])
                set([])
            """
            I = self.order_ideal( antichain )
            I_comp = [ elem for elem in self if elem not in I ]
            return set(self.order_ideal_generators(I_comp, direction = direction))

        panyushev_complement = order_ideal_complement_generators

        def panyushev_orbits(self, element_constructor = set):
            r"""
            Returns the Panyushev orbits of antichains in ``self``

            OUTPUT:

            - a partition of the antichains of ``self``, as a list of
              lists ``L`` such that for each ``L`` and ``i``, cyclically:
              ``self.order_ideal_complement_generators(L[i]) == L[i+1]``.

            EXAMPLES::

                sage: P = Poset( ( [1,2,3], [ [1,3], [2,3] ] ) )
                sage: P.panyushev_orbits()
                [[set([3]), set([]), set([1, 2])], [set([2]), set([1])]]
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

        def order_ideals_lattice(self):
            r"""
            Returns the lattice of order ideals of a poset `P`,
            ordered by inclusion. The usual notation is `J(P)`.

            EXAMPLES::

                sage: P = Posets.PentagonPoset(facade = True)
                sage: P.cover_relations()
                [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]
                sage: J = P.order_ideals_lattice(); J
                Finite lattice containing 8 elements
                sage: list(J)
                [{}, {0}, {0, 2}, {0, 1}, {0, 1, 2}, {0, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3, 4}]

            TESTS::

                sage: J = Posets.DiamondPoset(4, facade = True).order_ideals_lattice(); J
                Finite lattice containing 6 elements
                sage: list(J)
                [{}, {0}, {0, 2}, {0, 1}, {0, 1, 2}, {0, 1, 2, 3}]
                sage: J.cover_relations()
                [[{}, {0}], [{0}, {0, 2}], [{0}, {0, 1}], [{0, 2}, {0, 1, 2}], [{0, 1}, {0, 1, 2}], [{0, 1, 2}, {0, 1, 2, 3}]]

            .. note:: we use facade posets in the examples above just
               to ensure a nicer ordering in the output.
            """
            from sage.misc.misc import attrcall
            from sage.sets.set import Set
            from sage.combinat.posets.lattices import LatticePoset
            ideals = [Set( self.order_ideal(antichain) ) for antichain in self.antichains()]
            return LatticePoset((ideals,attrcall("issubset")))

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
