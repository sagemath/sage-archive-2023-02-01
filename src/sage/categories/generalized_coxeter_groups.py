r"""
Generalized Coxeter Groups
"""
#*****************************************************************************
#  Copyright (C) 2016 Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom, axiom
from sage.categories.groups import Groups
from sage.categories.complex_reflection_groups import ComplexReflectionGroups

class GeneralizedCoxeterGroups(Category_singleton):
    r"""
    The category of generalized Coxeter groups.

    A generalized Coxeter group is a group with a presentation of
    the following form:

    .. MATH::

        \langle s_i \mid s_i^{p_i}, s_i s_j \cdots = s_j s_i \cdots \rangle,

    where `p_i > 1`, `i \in I`, and the factors in the braid relation
    occur `m_{ij} = m_{ji}` times for all `i \neq j \in I`.

    EXAMPLES::

        sage: from sage.categories.generalized_coxeter_groups import GeneralizedCoxeterGroups
        sage: C = GeneralizedCoxeterGroups(); C
        Category of generalized coxeter groups

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.generalized_coxeter_groups import GeneralizedCoxeterGroups
            sage: GeneralizedCoxeterGroups().super_categories()
            [Category of finitely generated groups]
        """
        return [Groups().FinitelyGenerated()]

    class SubcategoryMethods:
        Irreducible = axiom("Irreducible")

    class Finite(CategoryWithAxiom):
        """
        The category of finite generalized Coxeter groups.
        """
        def extra_super_categories(self):
            """
            Implement that a finite generalized Coxeter group is a
            well-generated complex reflection group.

            EXAMPLES::

                sage: from sage.categories.generalized_coxeter_groups import GeneralizedCoxeterGroups
                sage: Cat = GeneralizedCoxeterGroups().Finite()
                sage: Cat.extra_super_categories()
                [Category of finite well generated complex reflection groups]
                sage: Cat.is_subcategory(ComplexReflectionGroups())
                True
            """
            return [ComplexReflectionGroups().Finite().WellGenerated()]

    class Irreducible(CategoryWithAxiom):
        """
        The category of irreducible generalized Coxeter groups.
        """

    class ParentMethods:
        @abstract_method
        def index_set(self):
            """
            Return the index set of (the simple reflections of)
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
            Implement: :meth:`Sets.ParentMethods.an_element` by
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
            Implement :meth:`Sets.ParentMethods.some_elements` by
            returning some typical element of ``self`.

            EXAMPLES::

                sage: W=WeylGroup(['A',3])
                sage: W.some_elements()
                [
                [0 1 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [0 0 0 1]
                [1 0 0 0]  [0 0 1 0]  [0 1 0 0]  [0 1 0 0]  [1 0 0 0]
                [0 0 1 0]  [0 1 0 0]  [0 0 0 1]  [0 0 1 0]  [0 1 0 0]
                [0 0 0 1], [0 0 0 1], [0 0 1 0], [0 0 0 1], [0 0 1 0]
                ]
                sage: W.order()
                24
            """
            return list(self.simple_reflections()) + [ self.one(), self.an_element() ]

        def simple_reflection_orders(self):
            """
            Return the orders of the simple reflections.
            """
            one = self.one()
            s = self.simple_reflections()
            return [s[i].order() for i in self.index_set()]

        def simple_reflection(self, i):
            """
            Return the simple reflection `s_i`.

            INPUT:

            - ``i`` -- an element from the index set

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
            if not i in self.index_set():
                raise ValueError("%s is not in the Dynkin node set %s"%(i,self.index_set()))
            return self.one().apply_simple_reflection(i) # don't care about left/right

        @cached_method
        def simple_reflections(self):
            r"""
            Return the simple reflections `(s_i)_{i\in I}` as a family.

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

        @cached_method
        def rank(self):
            r"""
            Return the rank of ``self``.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
                sage: W.rank()
                3
            """
            return len(self.simple_reflections())

        def group_generators(self):
            r"""
            Return the simple reflections of ``self``.

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

    class ElementMethods:
        # TODO: standardize / cleanup
        def apply_simple_reflections(self, word, side='right'):
            """
            Return the result of the (left/right) multiplication of
            word to ``self``.

            INPUT:

            - ``word`` -- a sequence of indices of Coxeter generators
            - ``side`` -- indicates multiplying from left or right

            EXAMPLES::

               sage: W = CoxeterGroups().example()
               sage: w = W.an_element(); w
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
            Return ``self`` multiplied by the simple reflection ``s[i]``
            on the left.

            This low level method is used intensively. Coxeter groups
            are encouraged to override this straightforward
            implementation whenever a faster approach exists.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
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
                'sage.categories.generalized_coxeter_groups'
            """
            s = self.parent().simple_reflections()
            return s[i] * self

        def apply_simple_reflection_right(self, i):
            """
            Return ``self`` multiplied by the simple reflection ``s[i]``
            on the right.

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
                'sage.categories.generalized_coxeter_groups'
            """
            s = self.parent().simple_reflections()
            return self * s[i]

        def apply_simple_reflection(self, i, side='right'):
            """
            Return ``self`` multiplied by the simple reflection ``s[i]``.

            INPUT:

            - ``i`` -- an element of the index set
            - ``side`` -- (default: ``"right"``) ``"left"`` or ``"right"``

            This default implementation simply calls
            :meth:`apply_simple_reflection_left` or
            :meth:`apply_simple_reflection_right`.

            EXAMPLES::

                sage: W = CoxeterGroups().example()
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
                'sage.categories.generalized_coxeter_groups'
            """
            if side == 'right':
                return self.apply_simple_reflection_right(i)
            else:
                return self.apply_simple_reflection_left(i)

        def _mul_(self, other):
            r"""
            Return the product of ``self`` and ``other``

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

            This method is called when computing ``self * other``::

                sage: w * w
                (1, 2, 1, 2)

            TESTS::

                sage: w._mul_.__module__
                'sage.categories.generalized_coxeter_groups'
            """
            return self.apply_simple_reflections(other.reduced_word())

        def inverse(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: W = WeylGroup(['B',7])
                sage: w = W.an_element()
                sage: u = w.inverse()
                sage: u == ~w
                True
                sage: u * w == w * u
                True
                sage: u * w
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

        def apply_conjugation_by_simple_reflection(self, i):
            r"""
            Conjugate ``self`` by the ``i``-th simple reflection.

            EXAMPLES::

                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([3,1,2,1])
                sage: w.apply_conjugation_by_simple_reflection(1).reduced_word()
                [3, 2]
            """
            return (self.apply_simple_reflection(i)).apply_simple_reflection(i, side='left')

