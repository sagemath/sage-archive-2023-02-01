r"""
Coxeter Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
# With contributions from Dan Bump, Steve Pon, Qiang Wang, Anne Schilling

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.abstract_method import abstract_method
from sage.misc.constant_function import ConstantFunction
from sage.misc.misc import attrcall
from sage.categories.category_singleton import Category_singleton
from sage.categories.groups import Groups
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.sage_object import have_same_parent
from sage.combinat.finite_class import FiniteCombinatorialClass
from sage.misc.flatten import flatten

class CoxeterGroups(Category_singleton):
    r"""
    The category of Coxeter groups.

    A *Coxeter group* is a group `W` with a distinguished (finite)
    family of involutions `(s_i)_{i\in I}`, called the *simple
    reflections*, subject to relations of the form `(s_is_j)^{m_{i,j}} = 1`.

    `I` is the *index set* of `W` and `|I|` is the *rank* of `W`.

    See http://en.wikipedia.org/wiki/Coxeter_group for details.

    EXAMPLES::

        sage: C = CoxeterGroups()
        sage: C                            # todo: uppercase for Coxeter
        Category of coxeter groups
        sage: C.super_categories()
        [Category of groups, Category of enumerated sets]

        sage: W = C.example(); W
        The symmetric group on {0, ..., 3}

        sage: W.simple_reflections()
        Finite family {0: (1, 0, 2, 3), 1: (0, 2, 1, 3), 2: (0, 1, 3, 2)}

    Here are some further examples::

        sage: FiniteCoxeterGroups().example()
        The 5-th dihedral group of order 10
        sage: FiniteWeylGroups().example()
        The symmetric group on {0, ..., 3}
        sage: WeylGroup(["B", 3])
        Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)

    Those will eventually be also in this category::

        sage: SymmetricGroup(4)
        Symmetric group of order 4! as a permutation group
        sage: DihedralGroup(5)
        Dihedral group of order 10 as a permutation group

    TODO: add a demo of usual computations on Coxeter groups.

    SEE ALSO: :class:`WeylGroups`, :mod:`sage.combinat.root_system`

    TESTS::

        sage: W = CoxeterGroups().example(); TestSuite(W).run(verbose = "True")
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_has_descent() . . . pass
        running ._test_inverse() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_reduced_word() . . . pass
        running ._test_simple_projections() . . . pass
        running ._test_some_elements() . . . pass
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CoxeterGroups().super_categories()
            [Category of groups, Category of enumerated sets]
        """
        return [Groups(), EnumeratedSets()]

    class ParentMethods:

        @abstract_method
        def index_set(self):
            """
            Returns the index set of (the simple reflections of)
            ``self``, as a list (or iterable).

            EXAMPLES::

                sage: W = FiniteCoxeterGroups().example(); W
                The 5-th dihedral group of order 10
                sage: W.index_set()
                [1, 2]
            """
            # return self.simple_reflections().keys()

        def _an_element_(self):
            """
            Implements: :meth:`Sets.ParentMethods.an_element` by
            returning the product of the simple reflections (a Coxeter
            element).

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: W.an_element()               # indirect doctest
                (1, 2, 3, 0)

            """
            return self.prod(self.simple_reflections())

        def some_elements(self):
            """
            Implements :meth:`Sets.ParentMethods.some_elements` by
            returning some typical element of `self`.

            EXAMPLES::

                sage: W=WeylGroup(['A',3])
                sage: W.some_elements()
                [[0 1 0 0]
                [1 0 0 0]
                [0 0 1 0]
                [0 0 0 1],
                 [1 0 0 0]
                [0 0 1 0]
                [0 1 0 0]
                [0 0 0 1],
                 [1 0 0 0]
                [0 1 0 0]
                [0 0 0 1]
                [0 0 1 0],
                 [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1],
                 [0 0 0 1]
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]]
                sage: W.order()
                24
            """
            return list(self.simple_reflections()) + [ self.one(), self.an_element() ]

        def __iter__(self):
            r"""
            Returns an iterator over the elements of this Coxeter group.

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

                sage: W = WeylGroup(["A",2,1])
                sage: g = iter(W)
                sage: g.next()
                [1 0 0]
                [0 1 0]
                [0 0 1]
                sage: g.next()
                [-1  1  1]
                [ 0  1  0]
                [ 0  0  1]
                sage: g.next()
                [ 0 -1  2]
                [ 1 -1  1]
                [ 0  0  1]
            """
            return iter(self.weak_order_ideal(predicate = ConstantFunction(True)))

        def weak_order_ideal(self, predicate, side ="right", category = None):
            """
            Returns a weak order ideal defined by a predicate

            INPUT:

            - ``predicate``: a predicate on the elements of ``self`` defining an
              weak order ideal in ``self``
            - ``side``: "left" or "right" (default: "right")

            OUTPUT: an enumerated set

            EXAMPLES::

                sage: D6 = FiniteCoxeterGroups().example(5)
                sage: I = D6.weak_order_ideal(predicate = lambda w: w.length() <= 3)
                sage: I.cardinality()
                7
                sage: list(I)
                [(), (1,), (1, 2), (1, 2, 1), (2,), (2, 1), (2, 1, 2)]

            We now consider an infinite Coxeter group::

                sage: W = WeylGroup(["A",1,1])
                sage: I = W.weak_order_ideal(predicate = lambda w: w.length() <= 2)
                sage: list(iter(I))
                [[1 0]
                 [0 1],
                 [-1  2]
                 [ 0  1],
                 [ 3 -2]
                 [ 2 -1],
                 [ 1  0]
                 [ 2 -1],
                 [-1  2]
                 [-2  3]]

            Even when the result is finite, some features of
            :class:`FiniteEnumeratedSets` are not available::

                sage: I.cardinality() # todo: not implemented
                5
                sage: list(I)         # todo: not implemented

            unless this finiteness is explicitly specified::

                sage: I = W.weak_order_ideal(predicate = lambda w: w.length() <= 2,
                ...                          category = FiniteEnumeratedSets())
                sage: I.cardinality()
                5
                sage: list(I)
                [[1 0]
                 [0 1],
                 [-1  2]
                 [ 0  1],
                 [ 3 -2]
                 [ 2 -1],
                 [ 1  0]
                 [ 2 -1],
                 [-1  2]
                 [-2  3]]

            .. rubric:: Background

            The weak order is returned as a :class:`SearchForest`.
            This is achieved by assigning to each element `u1` of the
            ideal a single ancestor `u=u1 s_i`, where `i` is the
            smallest descent of `u`.

            This allows for iterating through the elements in
            roughly Constant Amortized Time and constant memory
            (taking the operations and size of the generated objects
            as constants).
            """
            from sage.combinat.backtrack import SearchForest
            def succ(u):
                for i in u.descents(positive = True, side = side):
                    u1 = u.apply_simple_reflection(i, side)
                    if i == u1.first_descent(side = side) and predicate(u1):
                        yield u1
                return
            from sage.categories.finite_coxeter_groups import FiniteCoxeterGroups
            default_category = FiniteEnumeratedSets() if self in FiniteCoxeterGroups() else EnumeratedSets()
            return SearchForest((self.one(),), succ, category = default_category.or_subcategory(category))

        def grassmannian_elements(self, side = "right"):
            """
            INPUT:

            - ``side``: "left" or "right" (default: "right")

            Returns the left or right grassmanian elements of self, as an enumerated set

            EXAMPLES::

                sage: S = CoxeterGroups().example()
                sage: G = S.grassmannian_elements()
                sage: G.cardinality()
                12
                sage: G.list()
                [(0, 1, 2, 3), (1, 0, 2, 3), (2, 0, 1, 3), (3, 0, 1, 2), (0, 2, 1, 3), (1, 2, 0, 3), (0, 3, 1, 2), (1, 3, 0, 2), (2, 3, 0, 1), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)]
                sage: sorted(tuple(w.descents()) for w in G)
                [(), (0,), (0,), (0,), (1,), (1,), (1,), (1,), (1,), (2,), (2,), (2,)]
                sage: G = S.grassmannian_elements(side = "left")
                sage: G.cardinality()
                12
                sage: sorted(tuple(w.descents(side = "left")) for w in G)
                [(), (0,), (0,), (0,), (1,), (1,), (1,), (1,), (1,), (2,), (2,), (2,)]
            """
            order_side = "left" if side == "right" else "right"
            return self.weak_order_ideal(attrcall("is_grassmannian", side = side), side = order_side)

        def from_reduced_word(self, word):
            r"""
            INPUT:

            - ``word`` - a list (or iterable) of elements of ``self.index_set()``

            Returns the group element corresponding to the given
            word. Namely, if ``word`` is `[i_1,i_2,\ldots,i_k]`, then
            this returns the corresponding product of simple
            reflections `s_{i_1} s_{i_2} \cdots s_{i_k}`.

            Note: the main use case is for constructing elements from
            reduced words, hence the name of this method. But actually
            the input word need *not* be reduced.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: s = W.simple_reflections()
                sage: W.from_reduced_word([0,2,0,1])
                (0, 3, 1, 2)
                sage: W.from_reduced_word((0,2,0,1))
                (0, 3, 1, 2)
                sage: s[0]*s[2]*s[0]*s[1]
                (0, 3, 1, 2)

            See also :meth:'._test_reduced_word'::

                sage: W._test_reduced_word()

            TESTS::

                sage: W=WeylGroup(['E',6])
                sage: W.from_reduced_word([2,3,4,2])
                [ 0  1  0  0  0  0  0  0]
                [ 0  0 -1  0  0  0  0  0]
                [-1  0  0  0  0  0  0  0]
                [ 0  0  0  1  0  0  0  0]
                [ 0  0  0  0  1  0  0  0]
                [ 0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  1  0]
                [ 0  0  0  0  0  0  0  1]

            """
            return self.one().apply_simple_reflections(word, side = 'right')

        def _test_reduced_word(self, **options):
            """
            Runs sanity checks on :meth:'CoxeterGroups.ElementMethods.reduced_word' and
            :meth:'.from_reduced_word`.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W._test_reduced_word()

            """
            tester = self._tester(**options)
            s = self.simple_reflections()
            for x in tester.some_elements():
                red = x.reduced_word()
                tester.assertEquals(self.from_reduced_word(red), x)
                tester.assertEquals(self.prod((s[i] for i in red)), x)

        def simple_reflection(self, i):
            """
            INPUT:

            - ``i`` - an element from the index set.

            Returns the simple reflection `s_i`

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: W.simple_reflection(1)
                (0, 2, 1, 3)
                sage: s = W.simple_reflections()
                sage: s[1]
                (0, 2, 1, 3)

            """
            assert(i in self.index_set())
            return self.one().apply_simple_reflection(i) # don't care about left/right

        @cached_method
        def simple_reflections(self):
            r"""
            Returns the simple reflections `(s_i)_{i\in I}`, as a family.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: s = W.simple_reflections()
                sage: s
                Finite family {0: (1, 0, 2, 3), 1: (0, 2, 1, 3), 2: (0, 1, 3, 2)}
                sage: s[0]
                (1, 0, 2, 3)
                sage: s[1]
                (0, 2, 1, 3)
                sage: s[2]
                (0, 1, 3, 2)


            This default implementation uses :meth:`.index_set` and
            :meth:`.simple_reflection`.
            """
            from sage.sets.family import Family
            return Family(self.index_set(), self.simple_reflection)

        def group_generators(self):
            r"""
            Implements :meth:`Groups.ParentMethods.group_generators`
            by returning the simple reflections of ``self``.

            EXAMPLES::

                sage: D10 = FiniteCoxeterGroups().example(10)
                sage: D10.group_generators()
                Finite family {1: (1,), 2: (2,)}
                sage: SymmetricGroup(5).group_generators()
                Finite family {1: (1,2), 2: (2,3), 3: (3,4), 4: (4,5)}

            Those give semigroup generators, even for an infinite group::

                sage: W = WeylGroup(["A",2,1])
                sage: W.semigroup_generators()
                Finite family {0: [-1  1  1]
                                  [ 0  1  0]
                                  [ 0  0  1],
                               1: [ 1  0  0]
                                  [ 1 -1  1]
                                  [ 0  0  1],
                               2: [ 1  0  0]
                                  [ 0  1  0]
                                  [ 1  1 -1]}
            """
            return self.simple_reflections()

        semigroup_generators = group_generators

        def simple_projection(self, i, side = 'right', toward_max = True):
            r"""
            INPUT:

            - ``i`` - an element of the index set of self

            Returns the simple projection `\pi_i` (or `\overline\pi_i` if toward_max is False).

            See :meth:`.simple_projections` for the options. and for
            the definition of the simple projections.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: s = W.simple_reflections()
                sage: sigma=W.an_element()
                sage: sigma
                (1, 2, 3, 0)
                sage: u0=W.simple_projection(0)
                sage: d0=W.simple_projection(0,toward_max=False)
                sage: sigma.length()
                3
                sage: pi=sigma*s[0]
                sage: pi.length()
                4
                sage: u0(sigma)
                (2, 1, 3, 0)
                sage: pi
                (2, 1, 3, 0)
                sage: u0(pi)
                (2, 1, 3, 0)
                sage: d0(sigma)
                (1, 2, 3, 0)
                sage: d0(pi)
                (1, 2, 3, 0)

            """
            assert(i in self.index_set() or i==0) # i==0 allows for affine descents
            return lambda x: x.apply_simple_projection(i, side = side, toward_max = toward_max) # should use default_keyword

        @cached_method
        def simple_projections(self, side = 'right', toward_max = True):
            r"""
            INPUT:

            - ``self`` - a Coxeter group `W`
            - ``side`` - 'left' or 'right' (default: 'right')
            - ``toward_max`` - a boolean (default: True) specifying
              the direction of the projection

            Returns the simple projections of `W`, as a family.

            To each simple reflection `s_i` of `W`, corresponds a
            *simple projection* `\pi_i` from `W` to `W` defined by:

                      `\pi_i(w) = w s_i` if `i` is not a descent of `w`
                      `\pi_i(w) = w` otherwise.

            The simple projections `(\pi_i)_{i\in I}` move elements
            down the right permutohedron, toward the maximal element.
            They satisfy the same relations as the simple reflections,
            except that `\pi_i^2=\pi`, whereas `s_i^2 = s_i`. As such,
            the simple projections generate the `0`-Hecke monoid.

            By symmetry, one can also define the projections
            `(\overline\pi_i)_{i\in I}` (option ``toward_max``):

                      `\overline\pi_i(w) = w s_i` if `i` is a descent of `w`
                      `\overline\pi_i(w) = w` otherwise.

            as well as the analogues acting on the left (option ``side``).

            TODO: the name `toward_max` is not explicit and
            should/will be replaced; suggestions welcome.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: s = W.simple_reflections()
                sage: sigma=W.an_element()
                sage: sigma
                (1, 2, 3, 0)
                sage: pi=W.simple_projections()
                sage: pi
                Finite family {0: <function <lambda> at ...>, 1: <function <lambda> at ...>, 2: <function <lambda> ...>}
                sage: pi[1](sigma)
                (1, 3, 2, 0)
                sage: W.simple_projection(1)(sigma)
                (1, 3, 2, 0)
            """
            from sage.sets.family import Family
            return Family(self.index_set(), lambda i: self.simple_projection(i, side = side, toward_max = toward_max))

        def bruhat_interval(self, x, y):
            """
            Returns the list of t such that x <= t <= y.

            EXAMPLES::

                sage: W = WeylGroup("A3", prefix="s")
                sage: [s1,s2,s3]=W.simple_reflections()
                sage: W.bruhat_interval(s2,s1*s3*s2*s1*s3)
                [s1*s2*s3*s2*s1, s2*s3*s2*s1, s3*s1*s2*s1, s1*s2*s3*s1, s1*s2*s3*s2, s3*s2*s1, s2*s3*s1, s2*s3*s2, s1*s2*s1, s3*s1*s2, s1*s2*s3, s2*s1, s3*s2, s2*s3, s1*s2, s2]
                sage: W = WeylGroup(['A',2,1], prefix="s")
                sage: [s0,s1,s2]=W.simple_reflections()
                sage: W.bruhat_interval(1,s0*s1*s2)
                [s0*s1*s2, s1*s2, s0*s2, s0*s1, s2, s1, s0, 1]
            """
            if x == 1:
                x = self.one()
            if y == 1:
                y = self.one()
            if x == y:
                return [x]
            ret = []
            if not x.bruhat_le(y):
                return ret
            ret.append([y])
            while ret[-1] != []:
                nextlayer = []
                for z in ret[-1]:
                    for t in z.bruhat_lower_covers():
                        if t not in nextlayer:
                            if x.bruhat_le(t):
                                nextlayer.append(t)
                ret.append(nextlayer)
            return flatten(ret)

        # TODO: Groups() should have inverse() call __invert__
        # With strong doc stating that this is just a convenience for the user
        # and links to ~ / __invert__

        # parabolic_subgroup

        def _test_simple_projections(self, **options):
            """
            Runs sanity checks on :meth:`.simple_projections`
            and :meth:`CoxeterGroups.ElementMethods.apply_simple_projection`

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W._test_simple_projections()
            """
            tester = self._tester(**options)
            for side in ['left', 'right']:
                pi  = self.simple_projections(side = side)
                opi = self.simple_projections(side = side, toward_max = False)
                for i in self.index_set():
                    for w in tester.some_elements():
                        tester.assert_( pi[i](w) == w.apply_simple_projection(i, side = side))
                        tester.assert_( pi[i](w) == w.apply_simple_projection(i, side = side, toward_max = True ))
                        tester.assert_(opi[i](w) == w.apply_simple_projection(i, side = side, toward_max = False))
                        tester.assert_( pi[i](w).has_descent(i, side = side))
                        tester.assert_(not opi[i](w).has_descent(i, side = side))
                        tester.assertEquals(set([pi[i](w), opi[i](w)]),
                                            set([w, w.apply_simple_reflection(i, side = side)]))


        def _test_has_descent(self, **options):
            """
            Runs sanity checks on the method
            :meth:`CoxeterGroups.ElementMethods.has_descent` of the
            elements of self.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W._test_has_descent()
            """
            tester = self._tester(**options)
            s = self.simple_reflections()
            for i in self.index_set():
                tester.assert_(not self.one().has_descent(i))
                tester.assert_(not self.one().has_descent(i, side = 'left'))
                tester.assert_(not self.one().has_descent(i, side = 'right'))
                tester.assert_(self.one().has_descent(i, positive = True))
                tester.assert_(self.one().has_descent(i, positive = True, side = 'left'))
                tester.assert_(self.one().has_descent(i, positive = True, side = 'right'))
                for j in self.index_set():
                    tester.assertEquals(s[i].has_descent(j, side = 'left' ), i==j)
                    tester.assertEquals(s[i].has_descent(j, side = 'right'), i==j)
                    tester.assertEquals(s[i].has_descent(j                ), i==j)
                    tester.assertEquals(s[i].has_descent(j, positive = True, side = 'left' ), i!=j)
                    tester.assertEquals(s[i].has_descent(j, positive = True, side = 'right'), i!=j)
                    tester.assertEquals(s[i].has_descent(j, positive = True,               ), i!=j)
                    if i == j:
                        continue
                    u = s[i] * s[j]
                    v = s[j] * s[i]
                    tester.assert_((s[i]*s[j]).has_descent(i, side = 'left' ))
                    tester.assert_((s[i]*s[j]).has_descent(j, side = 'right'))
                    tester.assertEquals((s[i]*s[j]).has_descent(j, side = 'left' ), u == v)
                    tester.assertEquals((s[i]*s[j]).has_descent(i, side = 'right'), u == v)

    class ElementMethods:
        def has_descent(self, i, side = 'right', positive=False):
            """
            Returns whether i is a (left/right) descent of self.

            See :meth:`.descents` for a description of the options.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: s = W.simple_reflections()
                sage: w = s[0] * s[1] * s[2]
                sage: w.has_descent(2)
                True
                sage: [ w.has_descent(i)                  for i in [0,1,2] ]
                [False, False, True]
                sage: [ w.has_descent(i, side = 'left')   for i in [0,1,2] ]
                [True, False, False]
                sage: [ w.has_descent(i, positive = True) for i in [0,1,2] ]
                [True, True, False]

            This default implementation delegates the work to
            :meth:`.has_left_descent` and :meth:`.has_right_descent`.
            """
            assert isinstance(positive, bool)
            if side == 'right':
                return self.has_right_descent(i) != positive
            else:
                assert side == 'left'
                return self.has_left_descent(i)  != positive

#        @abstract_method(optional = True)
        def has_right_descent(self, i):
            """
            Returns whether ``i`` is a right descent of self.

            EXAMPLES::

                sage: W = CoxeterGroups().example(); W
                The symmetric group on {0, ..., 3}
                sage: w = W.an_element(); w
                (1, 2, 3, 0)
                sage: w.has_right_descent(0)
                False
                sage: w.has_right_descent(1)
                False
                sage: w.has_right_descent(2)
                True
            """
            return (~self).has_left_descent(i)

        def has_left_descent(self, i):
            """
            Returns whether `i` is a left descent of self.

            This default implementation uses that a left descent of
            `w` is a right descent of `w^{-1}`.

            EXAMPLES::

                sage: W = CoxeterGroups().example(); W
                The symmetric group on {0, ..., 3}
                sage: w = W.an_element(); w
                (1, 2, 3, 0)
                sage: w.has_left_descent(0)
                True
                sage: w.has_left_descent(1)
                False
                sage: w.has_left_descent(2)
                False

            TESTS::

                sage: w.has_left_descent.__module__
                'sage.categories.coxeter_groups'
            """
            return (~self).has_right_descent(i)

        def first_descent(self, side = 'right', index_set=None, positive=False):
            """
            Returns the first left (resp. right) descent of self, as
            ane element of ``index_set``, or ``None`` if there is none.

            See :meth:`.descents` for a description of the options.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: s = W.simple_reflections()
                sage: w = s[2]*s[0]
                sage: w.first_descent()
                0
                sage: w = s[0]*s[2]
                sage: w.first_descent()
                0
                sage: w = s[0]*s[1]
                sage: w.first_descent()
                1
            """
            if index_set is None:
                index_set = self.parent().index_set()
            for i in index_set:
                if self.has_descent(i, side = side, positive = positive):
                    return i
            return None


        def descents(self, side = 'right', index_set=None, positive=False):
            """
            INPUT:

            - ``index_set`` - a subset (as a list or iterable) of the nodes of the dynkin diagram;
              (default: all of them)
            - ``side`` - 'left' or 'right' (default: 'right')
            - ``positive`` - a boolean (default: ``False``)

            Returns the descents of self, as a list of elements of the
            index_set.

            The ``index_set`` option can be used to restrict to the
            parabolic subgroup indexed by ``index_set``.

            If positive is ``True``, then returns the non-descents
            instead

            TODO: find a better name for ``positive``: complement? non_descent?

            Caveat: the return type may change to some other iterable
            (tuple, ...) in the future. Please use keyword arguments
            also, as the order of the arguments may change as well.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: s=W.simple_reflections()
                sage: w=s[0]*s[1]
                sage: w.descents()
                [1]
                sage: w=s[0]*s[2]
                sage: w.descents()
                [0, 2]

                TODO: side, index_set, positive
            """
            if index_set==None:
                index_set=self.parent().index_set()
            return [ i for i in index_set if self.has_descent(i, side = side, positive = positive) ]

        def is_grassmannian(self, side = "right"):
            """
            INPUT:

            - ``side`` - "left" or "right" (default: "right")

            Tests whether ``self`` is Grassmannian, i.e. it has at
            most one descent on the right (resp. on the left).

            EXAMPLES::

                sage: W = CoxeterGroups().example(); W
                The symmetric group on {0, ..., 3}
                sage: s = W.simple_reflections()
                sage: W.one().is_grassmannian()
                True
                sage: s[1].is_grassmannian()
                True
                sage: (s[1]*s[2]).is_grassmannian()
                True
                sage: (s[0]*s[1]).is_grassmannian()
                True
                sage: (s[1]*s[2]*s[1]).is_grassmannian()
                False

                sage: (s[0]*s[2]*s[1]).is_grassmannian(side = "left")
                False
                sage: (s[0]*s[2]*s[1]).is_grassmannian(side = "right")
                True
                sage: (s[0]*s[2]*s[1]).is_grassmannian()
                True
            """
            return len(self.descents(side = side)) <= 1

        def reduced_word_reverse_iterator(self):
            """
            Returns a reverse iterator on a reduced word for self.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: s = W.simple_reflections()
                sage: sigma = s[0]*s[1]*s[2]
                sage: rI=sigma.reduced_word_reverse_iterator()
                sage: [i for i in rI]
                [2, 1, 0]
                sage: s[0]*s[1]*s[2]==sigma
                True
                sage: sigma.length()
                3

            SEE ALSO :meth:`.reduced_word`

            Default implementation: recursively remove the first right
            descent until the identity is reached (see :meth:`.first_descent` and
            :meth:`apply_simple_reflection`).

            """
            while True:
                i = self.first_descent()
                if i is None:
                    return
                self = self.apply_simple_reflection(i, 'right')
                yield i

        def reduced_word(self):
            r"""
            Returns a reduced word for self.

            This is a word `[i_1,i_2,\ldots,i_k]` of minimal length
            such that `s_{i_1} s_{i_2} \cdots s_{i_k}=self`, where `s`
            are the simple reflections.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: s=W.simple_reflections()
                sage: w=s[0]*s[1]*s[2]
                sage: w.reduced_word()
                [0, 1, 2]
                sage: w=s[0]*s[2]
                sage: w.reduced_word()
                [2, 0]

            SEE ALSO: :meth:`.reduced_words`, :meth:`.reduced_word_reverse_iterator`, :meth:`length`

            """
            result = list(self.reduced_word_reverse_iterator())
            return list(reversed(result))

        #def lex_min_reduced_word(w):
        #    return list(reversed((w.inverse()).reduced_word()))

        def reduced_words(self):
            r"""
            Returns all reduced words for self.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: s=W.simple_reflections()
                sage: w=s[0]*s[2]
                sage: w.reduced_words()
                [[2, 0], [0, 2]]
                sage: W=WeylGroup(['E',6])
                sage: w=W.from_reduced_word([2,3,4,2])
                sage: w.reduced_words()
                [[3, 2, 4, 2], [2, 3, 4, 2], [3, 4, 2, 4]]

            TODO: the result should be full featured finite enumerated
            set (e.g. counting can be done much faster than iterating).
            """
            descents = self.descents()
            if descents == []:
                return [[]]
            else:
                return [ r + [i]
                         for i in self.descents()
                         for r in (self.apply_simple_reflection(i)).reduced_words()
                         ]

        def length(self):
            r"""
            Returns the length of self, that is the minimal length of
            a product of simple reflections giving self.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: s1 = W.simple_reflection(1)
                sage: s2 = W.simple_reflection(2)
                sage: s1.length()
                1
                sage: (s1*s2).length()
                2
                sage: W = CoxeterGroups().example()
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[1]*s[0]
                sage: w.length()
                3
                sage: W = CoxeterGroups().example()
                sage: sum((x^w.length()) for w in W) - expand(prod(sum(x^i for i in range(j+1)) for j in range(4))) # This is scandalously slow!!!
                0

            SEE ALSO: :meth:`.reduced_word`

            TODO: Should use reduced_word_iterator (or reverse_iterator)

            """
            return len(self.reduced_word())

        def coset_representative(self, index_set, side = 'right'):
            r"""
            INPUT:

            - ``index_set`` - a subset (or iterable) of the nodes of the dynkin diagram
            - ``side`` - 'left' or 'right'

            Returns the unique shortest element of the coxeter group
            $W$ which is in the same left (resp. right) coset as
            ``self``, with respect to the parabolic subgroup $W_I$.

            EXAMPLES::
                sage: W = CoxeterGroups().example(5)
                sage: s = W.simple_reflections()
                sage: w = s[2]*s[1]*s[3]
                sage: w.coset_representative([]).reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([1]).reduced_word()
                [2, 3]
                sage: w.coset_representative([1,2]).reduced_word()
                [2, 3]
                sage: w.coset_representative([1,3]                 ).reduced_word()
                [2]
                sage: w.coset_representative([2,3]                 ).reduced_word()
                [2, 1]
                sage: w.coset_representative([1,2,3]               ).reduced_word()
                []
                sage: w.coset_representative([],      side = 'left').reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([1],     side = 'left').reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([1,2],   side = 'left').reduced_word()
                [3]
                sage: w.coset_representative([1,3],   side = 'left').reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([2,3],   side = 'left').reduced_word()
                [1]
                sage: w.coset_representative([1,2,3], side = 'left').reduced_word()
                []

            """
            while True:
                i = self.first_descent(side = side, index_set = index_set)
                if i is None:
                    return self
                self = self.apply_simple_reflection(i, side = side)

        def apply_simple_projection(self, i, side = 'right', toward_max = True):
            r"""
            INPUT:

            - ``i`` - an element of the index set of the Coxeter group
            - ``side`` - 'left' or 'right' (default: 'right')
            - ``toward_max`` - a boolean (default: True) specifying
              the direction of the projection

            Returns the result of the application of the simple
            projection `\pi_i` (resp. `\overline\pi_i`) on self.

            See :meth:`CoxeterGroups.ParentMethods.simple_projections`
            for the definition of the simple projections.

            EXAMPLE::

               sage: W=CoxeterGroups().example()
               sage: w=W.an_element()
               sage: w
               (1, 2, 3, 0)
               sage: w.apply_simple_projection(2)
               (1, 2, 3, 0)
               sage: w.apply_simple_projection(2, toward_max=False)
               (1, 2, 0, 3)

            """
            if self.has_descent(i, side = side) is not toward_max:
                return self.apply_simple_reflection(i, side = side)
            else:
                return self

        def binary_factorizations(self, predicate = ConstantFunction(True)):
            """
            Returns the set of all the factorizations `self = u v` such
            that `l(self) = l(u) + l(v)`.

            Iterating through this set is Constant Amortized Time
            (counting arithmetic operations in the coxeter group as
            constant time) complexity, and memory linear in the length
            of `self`.

            One can pass as optional argument a predicate p such that
            `p(u)` implies `p(u')` for any `u` left factor of `self`
            and `u'` left factor of `u`. Then this returns only the
            factorizations `self = uv` such `p(u)` holds.

            EXAMPLES:

            We construct the set of all factorizations of the maximal
            element of the group::

                sage: W = WeylGroup(['A',3])
                sage: s = W.simple_reflections()
                sage: w0 = W.from_reduced_word([1,2,3,1,2,1])
                sage: w0.binary_factorizations().cardinality()
                24

            The same number of factorizations, by bounded length::

                sage: [w0.binary_factorizations(lambda u: u.length() <= l).cardinality() for l in [-1,0,1,2,3,4,5,6]]
                [0, 1, 4, 9, 15, 20, 23, 24]

            The number of factorizations of the elements just below
            the maximal element::

                sage: [(s[i]*w0).binary_factorizations().cardinality() for i in [1,2,3]]
                [12, 12, 12]
                sage: w0.binary_factorizations(lambda u: False).cardinality()
                0

            TESTS::

                sage: w0.binary_factorizations().category()
                Category of finite enumerated sets
            """
            from sage.combinat.backtrack import SearchForest
            W = self.parent()
            if not predicate(W.one()):
                return FiniteCombinatorialClass([])
            s = W.simple_reflections()
            def succ((u,v)):
                for i in v.descents(side = 'left'):
                    u1 = u * s[i]
                    if i == u1.first_descent() and predicate(u1):
                        yield (u1, s[i]*v)
            return SearchForest(((W.one(), self),), succ, category = FiniteEnumeratedSets())

        # TODO: standardize / cleanup
        def apply_simple_reflections(self, word, side = 'right'):
            """
            INPUT:

            - ``word`` -- A sequence of indices of Coxeter generators
            - ``side`` -- Indicates multiplying from left or right

            Returns the result of the (left/right) multiplication of
            word to self.  ``self`` is not changed.

            EXAMPLES::

               sage: W=CoxeterGroups().example()
               sage: w=W.an_element(); w
               (1, 2, 3, 0)
               sage: w.apply_simple_reflections([0,1])
               (2, 3, 1, 0)
               sage: w
               (1, 2, 3, 0)
               sage: w.apply_simple_reflections([0,1],side='left')
               (0, 1, 3, 2)
            """
            for i in word:
                self = self.apply_simple_reflection(i, side)
            return self


        def apply_simple_reflection_left(self, i):
            """
            Returns ``self`` multiplied by the simple reflection ``s[i]`` on the left

            This low level method is used intensively. Coxeter groups
            are encouraged to override this straightforward
            implementation whenever a faster approach exists.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: w = W.an_element(); w
                (1, 2, 3, 0)
                sage: w.apply_simple_reflection_left(0)
                (0, 2, 3, 1)
                sage: w.apply_simple_reflection_left(1)
                (2, 1, 3, 0)
                sage: w.apply_simple_reflection_left(2)
                (1, 3, 2, 0)

            TESTS::

                sage: w.apply_simple_reflection_left.__module__
                'sage.categories.coxeter_groups'
            """
            s = self.parent().simple_reflections()
            return s[i] * self

        def apply_simple_reflection_right(self, i):
            """
            Returns ``self`` multiplied by the simple reflection ``s[i]`` on the right

            This low level method is used intensively. Coxeter groups
            are encouraged to override this straightforward
            implementation whenever a faster approach exists.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: w = W.an_element(); w
                (1, 2, 3, 0)
                sage: w.apply_simple_reflection_right(0)
                (2, 1, 3, 0)
                sage: w.apply_simple_reflection_right(1)
                (1, 3, 2, 0)
                sage: w.apply_simple_reflection_right(2)
                (1, 2, 0, 3)

            TESTS::

                sage: w.apply_simple_reflection_right.__module__
                'sage.categories.coxeter_groups'
            """
            s = self.parent().simple_reflections()
            return self * s[i]

        def apply_simple_reflection(self, i, side = 'right'):
            """
            Returns ``self`` multiplied by the simple reflection ``s[i]``

            INPUT:

            - ``i`` -- an element of the index set
            - ``side`` -- "left" or "right" (default: "right")

            This default implementation simply calls
            :meth:`apply_simple_reflection_left` or
            :meth:`apply_simple_reflection_right`.

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: w = W.an_element(); w
                (1, 2, 3, 0)
                sage: w.apply_simple_reflection(0, side = "left")
                (0, 2, 3, 1)
                sage: w.apply_simple_reflection(1, side = "left")
                (2, 1, 3, 0)
                sage: w.apply_simple_reflection(2, side = "left")
                (1, 3, 2, 0)

                sage: w.apply_simple_reflection(0, side = "right")
                (2, 1, 3, 0)
                sage: w.apply_simple_reflection(1, side = "right")
                (1, 3, 2, 0)
                sage: w.apply_simple_reflection(2, side = "right")
                (1, 2, 0, 3)

            By default, ``side`` is "right"::

                sage: w.apply_simple_reflection(0)
                (2, 1, 3, 0)

            TESTS::

                sage: w.apply_simple_reflection_right.__module__
                'sage.categories.coxeter_groups'
            """
            if side == 'right':
                return self.apply_simple_reflection_right(i)
            else:
                return self.apply_simple_reflection_left(i)

        def _mul_(self, other):
            r"""
            Returns the product of ``self`` and ``other``

            This default implementation computes a reduced word of
            ``other`` using :meth:`reduced_word`, and applies the
            corresponding simple reflections on ``self`` using
            :meth:`apply_simple_reflections`.

            EXAMPLES::

                sage: W = FiniteCoxeterGroups().example(); W
                The 5-th dihedral group of order 10
                sage: w = W.an_element()
                sage: w
                (1, 2)
                sage: w._mul_(w)
                (1, 2, 1, 2)
                sage: w._mul_(w)._mul_(w)
                (2, 1, 2, 1)

            This method is called when computing ``self*other``::

                sage: w * w
                (1, 2, 1, 2)

            TESTS::

                sage: w._mul_.__module__
                'sage.categories.coxeter_groups'
            """
            return self.apply_simple_reflections(other.reduced_word())

        def inverse(self):
            """
            Returns the inverse of self

            EXAMPLES::

                sage: W=WeylGroup(['B',7])
                sage: w=W.an_element()
                sage: u=w.inverse()
                sage: u==~w
                True
                sage: u*w==w*u
                True
                sage: u*w
                [1 0 0 0 0 0 0]
                [0 1 0 0 0 0 0]
                [0 0 1 0 0 0 0]
                [0 0 0 1 0 0 0]
                [0 0 0 0 1 0 0]
                [0 0 0 0 0 1 0]
                [0 0 0 0 0 0 1]

            """

            return self.parent().one().apply_simple_reflections(self.reduced_word_reverse_iterator())

        __invert__ = inverse

        @cached_in_parent_method
        def bruhat_lower_covers(self):
            """
            Returns all elements that ``self`` covers in (strong) Bruhat order.

            If ``w = self`` has a descent at `i`, then the elements that
            `w` covers are exactly `\{ws_i, u_1s_i, u_2s_i,..., u_js_i\}`,
            where the `u_k` are elements that `ws_i` covers that also
            do not have a descent at `i`.

            EXAMPLES::

                sage: W = WeylGroup(["A",3])
                sage: w = W.from_reduced_word([3,2,3])
                sage: print([v.reduced_word() for v in w.bruhat_lower_covers()])
                [[3, 2], [2, 3]]

                sage: W = WeylGroup(["A",3])
                sage: print([v.reduced_word() for v in W.simple_reflection(1).bruhat_lower_covers()])
                [[]]
                sage: print([v.reduced_word() for v in W.one().bruhat_lower_covers()])
                []
                sage: W = WeylGroup(["B",4,1])
                sage: w = W.from_reduced_word([0,2])
                sage: print([v.reduced_word() for v in w.bruhat_lower_covers()])
                [[2], [0]]

            We now show how to construct the Bruhat poset::

                sage: W = WeylGroup(["A",3])
                sage: covers = tuple([u, v] for v in W for u in v.bruhat_lower_covers() )
                sage: P = Poset((W, covers), cover_relations = True)
                sage: P.show()

            Alternatively, one can just use::

                sage: P = W.bruhat_poset()

            The algorithm is taken from Stembridge's coxeter/weyl package for Maple.
            """
            desc = self.first_descent()
            if desc is not None:
                ww = self.apply_simple_reflection(desc)
                return [u.apply_simple_reflection(desc) for u in ww.bruhat_lower_covers() if not u.has_descent(desc)] + [ww]
            else:
                return []

        @cached_in_parent_method
        def bruhat_upper_covers(self):
            r"""
            Returns all elements that cover ``self`` in (strong) Bruhat order.

            The algorithm works recursively, using the 'inverse' of the method described for
            lower covers :meth:`bruhat_lower_covers`. Namely, it runs through all `i` in the
            index set: if `w`=``self`` has no right descent `i`, then `w s_i` is a cover;
            if `w` has a decent at `i`, then `u_j s_i` is a cover of `w` where `u_j` is a cover
            of `w s_i`.

            EXAMPLES::

                sage: W = WeylGroup(['A',3,1])
                sage: w = W.from_reduced_word([1,2,1])
                sage: sorted([x.reduced_word() for x in w.bruhat_upper_covers()])
                [[0, 1, 2, 1], [1, 2, 0, 1], [1, 2, 1, 0], [1, 2, 3, 1], [2, 3, 1, 2], [3, 1, 2, 1]]

                sage: W = WeylGroup(['A',3])
                sage: w = W.long_element()
                sage: w.bruhat_upper_covers()
                []

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([1,2,1])
                sage: S = [v for v in W if w in v.bruhat_lower_covers()]
                sage: C = w.bruhat_upper_covers()
                sage: set(S) == set(C)
                True
            """
            Covers = []
            for i in self.parent().index_set():
                if i in self.descents():
                    Covers += [ x.apply_simple_reflection(i) for x in self.apply_simple_reflection(i).bruhat_upper_covers()
                                if i not in x.descents() ]
                else:
                    Covers += [ self.apply_simple_reflection(i) ]
            return [x for x in set(Covers)]

        @cached_in_parent_method
        def bruhat_le(self, other):
            """
            Bruhat comparison

            INPUT:

            - other - an element of the same Coxeter group

            OUTPUT: a boolean

            Returns whether ``self`` <= ``other`` in the Bruhat order.

            EXAMPLES::

                sage: W = WeylGroup(["A",3])
                sage: u = W.from_reduced_word([1,2,1])
                sage: v = W.from_reduced_word([1,2,3,2,1])
                sage: u.bruhat_le(u)
                True
                sage: u.bruhat_le(v)
                True
                sage: v.bruhat_le(u)
                False
                sage: v.bruhat_le(v)
                True
                sage: s = W.simple_reflections()
                sage: s[1].bruhat_le(W.one())
                False

            The implementation uses the equivalent condition that any
            reduced word for ``other`` contains a reduced word for
            ``self`` as subword. See Stembridge, A short derivation of
            the Mobius function for the Bruhat order. J. Algebraic
            Combin. 25 (2007), no. 2, 141--148, Proposition 1.1.

            Complexity: `O(l * c)`, where `l` is the minimum of the
            lengths of `u` and of `v`, and `c` is the cost of the low
            level methods :meth:`first_descent`, :meth:`has_descent`,
            :meth:`apply_simple_reflection`, etc. Those are typically
            `O(n)`, where `n` is the rank of the Coxeter group.

            TESTS:

            We now run consistency tests with permutations and
            :meth:`bruhat_lower_covers`::

                sage: W = WeylGroup(["A",3])
                sage: P4 = Permutations(4)
                sage: def P4toW(w): return W.from_reduced_word(w.reduced_word())
                sage: for u in P4:
                ...       for v in P4:
                ...           assert u.bruhat_lequal(v) == P4toW(u).bruhat_le(P4toW(v))

                sage: W = WeylGroup(["B",3])
                sage: P = W.bruhat_poset() # This is built from bruhat_lower_covers
                sage: Q = Poset((W, attrcall("bruhat_le")))                             # long time (10s)
                sage: all( u.bruhat_le(v) == P.is_lequal(u,v) for u in W for v in W ) # long time  (7s)
                True
                sage: all( P.is_lequal(u,v) == Q.is_lequal(u,v) for u in W for v in W)       # long time  (9s)
                True

            """
            assert have_same_parent(self, other)
            # could first compare the length, when that information is cheap
            desc = other.first_descent()
            if desc is not None:
                return self.apply_simple_projection(desc, toward_max = False).bruhat_le(other.apply_simple_reflection(desc))
            else:
                return self == other

        def weak_le(self, other, side = 'right'):
            """
            comparison in weak order

            INPUT:

            - other - an element of the same Coxeter group
            - side - 'left' or 'right'  (default: 'right')

            OUTPUT: a boolean

            Returns whether ``self`` <= ``other`` in left
            (resp. right) weak order, that is if 'v' can be obtained
            from 'v' by length increasing multiplication by simple
            reflections on the left (resp. right).

            EXAMPLES::

                sage: W = WeylGroup(["A",3])
                sage: u = W.from_reduced_word([1,2])
                sage: v = W.from_reduced_word([1,2,3,2])
                sage: u.weak_le(u)
                True
                sage: u.weak_le(v)
                True
                sage: v.weak_le(u)
                False
                sage: v.weak_le(v)
                True

            Comparison for left weak order is achieved with the option ``side``::

                sage: u.weak_le(v, side = 'left')
                False

            The implementation uses the equivalent condition that any
            reduced word for `u` is a right (resp. left) prefix of
            some reduced word for `v`.

            Complexity: `O(l * c)`, where `l` is the minimum of the
            lengths of `u` and of `v`, and `c` is the cost of the low
            level methods :meth:`first_descent`, :meth:`has_descent`,
            :meth:`apply_simple_reflection`. Those are typically
            `O(n)`, where `n` is the rank of the Coxeter group.

            We now run consistency tests with permutations::

                sage: W = WeylGroup(["A",3])
                sage: P4 = Permutations(4)
                sage: def P4toW(w): return W.from_reduced_word(w.reduced_word())
                sage: for u in P4:  # long time (5s on sage.math, 2011)
                ...       for v in P4:
                ...           assert u.permutohedron_lequal(v) == P4toW(u).weak_le(P4toW(v))
                ...           assert u.permutohedron_lequal(v, side='left') == P4toW(u).weak_le(P4toW(v), side='left')
            """
            assert have_same_parent(self, other)
            # could first compare the length, when that information is cheap
            prefix_side = 'left' if side == 'right' else 'right'

            while True:
                desc = self.first_descent(side = prefix_side)
                if desc is None:
                    return True
                if not other.has_descent(desc, side = prefix_side):
                    return False
                self = self.apply_simple_reflection(desc, side = prefix_side)
                other = other.apply_simple_reflection(desc, side = prefix_side)

        def weak_covers(self, side = 'right', index_set = None, positive = False):
            """
            Returns all elements that ``self`` covers in weak order.

            INPUT:

            - side - 'left' or 'right'  (default: 'right')
            - positive - a boolean (default: False)
            - index_set - a list of indices or None

            OUTPUT: a list

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,2,1])
                sage: [x.reduced_word() for x in w.weak_covers()]
                [[3, 2]]

            To obtain instead elements that cover self, set ``positive = True``::

                sage: [x.reduced_word() for x in w.weak_covers(positive = True)]
                [[3, 1, 2, 1], [2, 3, 2, 1]]

            To obtain covers for left weak order, set the option side to 'left'::

                sage: [x.reduced_word() for x in w.weak_covers(side='left')]
                [[2, 1]]
                sage: w = W.from_reduced_word([3,2,3,1])
                sage: [x.reduced_word() for x in w.weak_covers()]
                [[2, 3, 2], [3, 2, 1]]
                sage: [x.reduced_word() for x in w.weak_covers(side='left')]
                [[3, 2, 1], [2, 3, 1]]

            Covers w.r.t. a parabolic subgroup are obtained with the option ``index_set``::

                sage: [x.reduced_word() for x in w.weak_covers(index_set = [1,2])]
                [[2, 3, 2]]
            """
            return [ self.apply_simple_reflection(i, side=side)
                     for i in self.descents(side=side, index_set = index_set, positive = positive) ]


        def lower_covers(self, side = 'right', index_set = None):
            """
            Returns all elements that ``self`` covers in weak order.

            INPUT:

            - side - 'left' or 'right' (default: 'right')
            - index_set - a list of indices or None

            OUTPUT: a list

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,2,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[3, 2]]

            To obtain covers for left weak order, set the option side to 'left'::

                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[2, 1]]
                sage: w = W.from_reduced_word([3,2,3,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[2, 3, 2], [3, 2, 1]]

            Covers w.r.t. a parabolic subgroup are obtained with the option ``index_set``::

                sage: [x.reduced_word() for x in w.lower_covers(index_set = [1,2])]
                [[2, 3, 2]]
                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[3, 2, 1], [2, 3, 1]]
            """
            return self.weak_covers(side = side, index_set = index_set, positive = False)

        def upper_covers(self, side = 'right', index_set = None):
            """
            Returns all elements that cover ``self`` in weak order.

            INPUT:

            - side - 'left' or 'right' (default: 'right')
            - index_set - a list of indices or None

            OUTPUT: a list

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([2,3])
                sage: [x.reduced_word() for x in w.upper_covers()]
                [[2, 3, 1], [2, 3, 2]]

            To obtain covers for left weak order, set the option ``side`` to 'left'::

                sage: [x.reduced_word() for x in w.upper_covers(side = 'left')]
                [[1, 2, 3], [2, 3, 2]]

            Covers w.r.t. a parabolic subgroup are obtained with the option ``index_set``::

                sage: [x.reduced_word() for x in w.upper_covers(index_set = [1])]
                [[2, 3, 1]]
                sage: [x.reduced_word() for x in w.upper_covers(side = 'left', index_set = [1])]
                [[1, 2, 3]]
            """
            return self.weak_covers(side = side, index_set = index_set, positive = True)

