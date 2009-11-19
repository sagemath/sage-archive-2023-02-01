r"""
Coxeter Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
# With contributions from Qiang Wang and Steve Pon

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.abstract_method import abstract_method
from sage.misc.constant_function import ConstantFunction
from sage.categories.category import Category
from sage.categories.groups import Groups
from sage.combinat.backtrack import SearchForest
from sage.combinat.finite_class import FiniteCombinatorialClass

class CoxeterGroups(Category):
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
        [Category of groups]

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
        running ._test_an_element() ... done
        running ._test_associativity() ... done
        running ._test_element_pickling() ... done
        running ._test_inverse() ... done
        running ._test_one() ... done
        running ._test_pickling() ... done
        running ._test_prod() ... done
        running ._test_reduced_word() ... done
        running ._test_some_elements() ... done
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CoxeterGroups().super_categories()
            [Category of groups]
        """
        return [Groups()]

    class ParentMethods:

        @abstract_method
        def index_set(self):
            """
            Returns the index set of (the simple reflections of)
            ``self``, as a list (or iterable).
            """
            # return self.simple_reflections().keys()

        @cached_method
        def an_element(self):
            """
            Implements: :meth:`Sets.ParentMethods.an_element` by
            returning the product of the simple reflections (a Coxeter
            element).

            EXAMPLES::

                sage: W=CoxeterGroups().example()
                sage: W
                The symmetric group on {0, ..., 3}
                sage: W.an_element()
                (1, 2, 3, 0)

            """
            return self.prod(self.simple_reflections())

        an_element_force = an_element

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
                 [0 1 0 0]
                [1 0 0 0]
                [0 0 1 0]
                [0 0 0 1]]
                sage: W.order()
                24
            """
            return list(self.simple_reflections()) + [ self.one(), self.an_element() ]

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
            return self.one().apply_simple_reflections(word, side = "right")

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
            Returns the simple reflections `(s_i)_{i\in I}, as a family.

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

        def semigroup_generators(self):
            r"""
            Implements :meth:`SemiGroups.ParentMethods.semigroup_generators`
            by returning the simple reflections of ``self``.

            EXAMPLES::

                sage: D10 = FiniteCoxeterGroups().example(10)
                sage: D10.semigroup_generators()
                Finite family {1: (1,), 2: (2,)}
                sage: D5 = FiniteCoxeterGroups().example(5)
                sage: D5.semigroup_generators()
                Finite family {1: (1,), 2: (2,)}
            """
            return self.simple_reflections()


        def simple_projection(self, i, side = "right", toward_max = True):
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
        def simple_projections(self, side = "right", toward_max = True):
            r"""
            INPUT:
             - ``self`` - a Coxeter group `W`
             - ``side`` - "left" or "right" (default: "right")
             - ``toward_max`` - a boolean (default: True) specifying the direction of the projection

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

        def _test_simple_projections(self, **options):
            """
            Runs sanity checks on :meth:`.simple_projections`
            and :meth:`CoxeterGroups.ElementMethods.apply_simple_projection`

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W._test_simple_projections()
            """
            tester = self._tester(**options)
            for side in ["left", "right"]:
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
                tester.assert_(not self.one().has_descent(i, side = "left"))
                tester.assert_(not self.one().has_descent(i, side = "right"))
                tester.assert_(self.one().has_descent(i, positive = True))
                tester.assert_(self.one().has_descent(i, positive = True, side = "left"))
                tester.assert_(self.one().has_descent(i, positive = True, side = "right"))
                for j in self.index_set():
                    tester.assertEquals(s[i].has_descent(j, side = "left" ), i==j)
                    tester.assertEquals(s[i].has_descent(j, side = "right"), i==j)
                    tester.assertEquals(s[i].has_descent(j                ), i==j)
                    tester.assertEquals(s[i].has_descent(j, positive = True, side = "left" ), i!=j)
                    tester.assertEquals(s[i].has_descent(j, positive = True, side = "right"), i!=j)
                    tester.assertEquals(s[i].has_descent(j, positive = True,               ), i!=j)
                    if i == j:
                        continue
                    u = s[i] * s[j]
                    v = s[j] * s[i]
                    tester.assert_((s[i]*s[j]).has_descent(i, side = "left" ))
                    tester.assert_((s[i]*s[j]).has_descent(j, side = "right"))
                    tester.assertEquals((s[i]*s[j]).has_descent(j, side = "left" ), u == v)
                    tester.assertEquals((s[i]*s[j]).has_descent(i, side = "right"), u == v)

    class ElementMethods:
        def has_descent(self, i, side = "right", positive=False):
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
                sage: [ w.has_descent(i, side = "left")   for i in [0,1,2] ]
                [True, False, False]
                sage: [ w.has_descent(i, positive = True) for i in [0,1,2] ]
                [True, True, False]

            This default implementation delegates the work to
            :meth:`.has_left_descent` and :meth:`.has_right_descent`.
            """
            assert isinstance(positive, bool)
            if side == "right":
                return self.has_right_descent(i) != positive
            else:
                assert side == "left"
                return self.has_left_descent(i)  != positive

        @abstract_method(optional = True)
        def has_right_descent(self, i):
            """
            Returns whether i is a right descent of self.

            EXAMPLES::

#                sage:
            """

        def has_left_descent(self, i):
            """
            Returns whether `i` is a left descent of self.

            This default implementation uses that a left descent of
            `w` is a right descent of `w^{-1}`.
            """
            return (~self).has_right_descent(i)

        def first_descent(self, side = "right", index_set=None, positive=False):
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


        def descents(self, side = "right", index_set=None, positive=False):
            """
            INPUT:
             - ``index_set`` - a subset (as a list or iterable) of the nodes of the dynkin diagram;
               (default: all of them)
             - ``side`` - "left" or "right" (default: "right")
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

            SEE ALSO :meth:`.reduced_word.

            Default implementation: recursively remove the first right
            descent until the identity is reached (see :meth:`.first_descent` and
            :meth:`apply_simple_reflection`).

            """
            while True:
                i = self.first_descent()
                if i is None:
                    return
                self = self.apply_simple_reflection(i, "right")
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

        def coset_representative(self, index_set, side = "right"):
            r"""
            INPUT:
             - ``index_set`` - a subset (or iterable) of the nodes of the dynkin diagram
             - ``side`` - "left" or "right"

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
                sage: w.coset_representative([],      side = "left").reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([1],     side = "left").reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([1,2],   side = "left").reduced_word()
                [3]
                sage: w.coset_representative([1,3],   side = "left").reduced_word()
                [2, 3, 1]
                sage: w.coset_representative([2,3],   side = "left").reduced_word()
                [1]
                sage: w.coset_representative([1,2,3], side = "left").reduced_word()
                []

            """
            while True:
                i = self.first_descent(side = side, index_set = index_set)
                if i is None:
                    return self
                self = self.apply_simple_reflection(i, side = side)

        def apply_simple_projection(self, i, side = "right", toward_max = True):
            """
            INPUT:
             - ``i`` - an element of the index set of the Coxeter group
             - ``side`` - "left" or "right" (default: "right")
             - ``toward_max`` - a boolean (default: True) specifying the direction of the projection

            Returns the result of the application of the simple
            projection `\pi_i` (resp. `\opi_i`) on self.

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
                sage: W = WeylGroup(['A',3])
                sage: s = W.simple_reflections()
                sage: w0 = W.from_reduced_word([1,2,3,1,2,1])
                sage: w0.binary_factorizations().cardinality()
                24
                sage: [w0.binary_factorizations(lambda u: u.length() <= l).cardinality() for l in [-1,0,1,2,3,4,5,6]]
                [0, 1, 4, 9, 15, 20, 23, 24]
                sage: [(s[i]*w0).binary_factorizations().cardinality() for i in [1,2,3]]
                [12, 12, 12]
                sage: w0.binary_factorizations(lambda u: False).cardinality()
                0

            """
            W = self.parent()
            if not predicate(W.one()):
                return FiniteCombinatorialClass([])
            s = W.simple_reflections()
            def succ((u,v)):
                for i in v.descents(side = "left"):
                    u1 = u * s[i]
                    if i == u1.first_descent() and predicate(u1):
                        yield (u1, s[i]*v)
            return SearchForest(((W.one(), self),), succ)

        # TODO: standardize / cleanup
        def apply_simple_reflections(self, word, side = "right"):
            """
            INPUT:
             - "word": A sequence of indices of Coxeter generators
             - "side": Indicates multiplying from left or right

            Returns the result of the (left/right) multiplication of
            word to self.  self is not changed.

            EXAMPLE::

               sage: W=CoxeterGroups().example()
               sage: w=W.an_element()
               sage: w
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

            ... multiplies s[i] by x...

            """
            s = self.parent().simple_reflections()
            return s[i] * self

        def apply_simple_reflection_right(self, i):
            """

            ... multiplies x by s[i] ...

            """
            s = self.parent().simple_reflections()
            return self * s[i]

        def apply_simple_reflection(self, i, side = "right"):
            """

            ... multiplies x by s[i] ...

            """
            if side == "right":
                return self.apply_simple_reflection_right(i)
            else:
                return self.apply_simple_reflection_left(i)

        def _mul_(self, other):
            """
            Returns the product of self and other, called when computing
            self*other

            EXAMPLES::

                sage: W=WeylGroup(['B',7])
                sage: w=W.an_element()
                sage: u=w.inverse()
                sage: u*w
                [1 0 0 0 0 0 0]
                [0 1 0 0 0 0 0]
                [0 0 1 0 0 0 0]
                [0 0 0 1 0 0 0]
                [0 0 0 0 1 0 0]
                [0 0 0 0 0 1 0]
                [0 0 0 0 0 0 1]

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
            Returns all elements that self covers in Bruhat order.
            If w = self has a descent at i, then the elements that
            w covers are exactly {ws_i, u_1s_i, u_2s_i,..., u_js_i},
            where the u_k are elements that ws_i covers that also
            do not have a descent at i

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

            We now show how to construct the Bruhat order poset.

            EXAMPLES::

                sage: elms = [w for w in WeylGroup(["A",3])]
                sage: fcn = lambda p,q : q in p.bruhat_lower_covers()
                sage: P = Poset((elms, fcn), cover_relations = True) # todo: not implemented (depends on poset_improvements-fs.patch)
                sage: P.show() # todo: not implemented (depends on poset_improvements-fs.patch)
            """
            desc = self.first_descent()
            if desc is not None:
                ww = self.apply_simple_reflection(desc)
                return [u.apply_simple_reflection(desc) for u in ww.bruhat_lower_covers() if not u.has_descent(desc)] + [ww]
            else:
                return []

        def weak_covers(self, side = 'right', index_set = None, positive = False):
            """
            Returns all elements that self covers in weak order.

            To obtain left weak order one can change the option side to 'left'.

            To obtain elements that cover self, set ``positive = False``.

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,2,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[3, 2]]
                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[2, 1]]
                sage: w = W.from_reduced_word([3,2,3,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[2, 3, 2], [3, 2, 1]]
                sage: [x.reduced_word() for x in w.lower_covers(index_set = [1,2])]
                [[2, 3, 2]]
                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[3, 2, 1], [2, 3, 1]]
            """
            return [ self.apply_simple_reflection(i,side=side) for i in self.descents(side=side, index_set = index_set, positive = positive) ]


        def lower_covers(self, side = 'right', index_set = None):
            """
            Returns all elements that self covers in right weak order.
            To obtain left weak order one can change the option side to 'left'.

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,2,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[3, 2]]
                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[2, 1]]
                sage: w = W.from_reduced_word([3,2,3,1])
                sage: [x.reduced_word() for x in w.lower_covers()]
                [[2, 3, 2], [3, 2, 1]]
                sage: [x.reduced_word() for x in w.lower_covers(index_set = [1,2])]
                [[2, 3, 2]]
                sage: [x.reduced_word() for x in w.lower_covers(side='left')]
                [[3, 2, 1], [2, 3, 1]]
            """
            return self.weak_covers(side = side, index_set = index_set, positive = False)

        def upper_covers(self, side = 'right', index_set = None):
            """
            Returns all elements that cover self in right weak order.
            To obtain left weak order one can change the option side to 'left'.

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,2,1])
                sage: [x.reduced_word() for x in w.upper_covers()]
                [[3, 1, 2, 1], [2, 3, 2, 1]]
                sage: [x.reduced_word() for x in w.upper_covers(index_set = [2])]
                [[3, 1, 2, 1]]
                sage: [x.reduced_word() for x in w.upper_covers(side='left')]
                [[3, 1, 2, 1], [2, 3, 2, 1]]
                sage: w = W.from_reduced_word([3,2,3,1])
                sage: [x.reduced_word() for x in w.upper_covers()]
                [[2, 3, 1, 2, 1]]
                sage: [x.reduced_word() for x in w.upper_covers(side='left')]
                [[1, 2, 3, 2, 1]]
            """
            return self.weak_covers(side = side, index_set = index_set, positive = True)

        # TODO: Groups() should have inverse() call __invert__
        # With strong doc stating that this is just a convenience for the user
        # and links to ~ / __invert__

        # bruhat_order()
        # weak_order(side = "left", "right")
        # parabolic_subgroup
