r"""
Finite Complex Reflection Groups
"""
#*****************************************************************************
#       Copyright (C) 2011-2015 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.all import prod
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom, axiom
from sage.categories.groups import Groups
from sage.rings.all import ZZ

class FiniteComplexReflectionGroups(CategoryWithAxiom):
    r"""
    The category of finite complex reflection groups.

    This is the base category for finite subgroups of the general
    linear group over a complex vector space which are generated 
    by complex reflections.

    EXAMPLES::

        sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
        sage: ComplexReflectionGroups().Finite()
        Category of finite complex reflection groups
        sage: ComplexReflectionGroups().Finite().super_categories()
        [Category of finite groups, Category of complex reflection groups]
        sage: ComplexReflectionGroups().Finite().all_super_categories()
        [Category of finite complex reflection groups,
         Category of finite groups,
         Category of finite monoids,
         Category of complex reflection groups,
         Category of groups,
         Category of monoids,
         Category of finite semigroups,
         Category of semigroups,
         Category of inverse unital magmas,
         Category of unital magmas,
         Category of magmas,
         Category of finite sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    An example of a finite reflection group::

        sage: W = ComplexReflectionGroups().Finite().example(); W
        Reducible complex reflection group of rank 4 and type A2 x B2

        sage: W.reflections()
        Finite family {0: (1,8)(2,5)(9,12), 1: (1,5)(2,9)(8,12),
                       2: (3,10)(4,7)(11,14), 3: (3,6)(4,11)(10,13),
                       4: (1,9)(2,8)(5,12), 5: (4,14)(6,13)(7,11),
                       6: (3,13)(6,10)(7,14)}

    ``W`` is in the category of complex reflection groups::

        sage: W in ComplexReflectionGroups().Finite()
        True
    """
    def example(self):
        r"""
        Return an example of a complex reflection group.

        EXAMPLES::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: ComplexReflectionGroups().Finite().example()
            Reducible complex reflection group of rank 4 and type A2 x B2
        """
        from sage.combinat.root_system.reflection_group_real import ReflectionGroup
        return ReflectionGroup((1,1,3), (2,1,2))

    class SubcategoryMethods:
        WellGenerated = axiom("WellGenerated")

    class ParentMethods:
        @abstract_method(optional=True)
        def degrees(self):
            r"""
            Return the degrees of ``self`` in increasing order.

            EXAMPLES::

                sage: W = ColoredPermutations(1,4)
                sage: W.degrees()
                [2, 3, 4]

                sage: W = ColoredPermutations(3,3)
                sage: W.degrees()
                [3, 6, 9]

                sage: W = ReflectionGroup(31)
                sage: W.degrees()
                [8, 12, 20, 24]
            """

        @abstract_method(optional=True)
        def codegrees(self):
            r"""
            Return the codegrees of ``self`` in decreasing order.

            EXAMPLES::

                sage: W = ColoredPermutations(1,4)
                sage: W.codegrees()
                [2, 1, 0]

                sage: W = ColoredPermutations(3,3)
                sage: W.codegrees()
                [6, 3, 0]

                sage: W = ReflectionGroup(31)
                sage: W.codegrees()
                [28, 16, 12, 0]
            """

        def nr_reflecting_hyperplanes(self):
            r"""
            Return the number of reflecting hyperplanes of ``self``.

            It is given by the sum of the codegrees of ``self``
            plus its rank.

            For real groups, this coincides with the number of
            reflections.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.nr_reflecting_hyperplanes()
                3
                sage: W = ColoredPermutations(2,3)
                sage: W.nr_reflecting_hyperplanes()
                9
                sage: W = ColoredPermutations(4,3)
                sage: W.nr_reflecting_hyperplanes()
                15
                sage: W = ReflectionGroup((4,2,3))
                sage: W.nr_reflecting_hyperplanes()
                15
            """
            return sum(self.codegrees()) + self.rank()

        def nr_reflections(self):
            r"""
            Return the number of reflections of ``self``.

            It is given by the sum of the degrees of self minus its
            rank.

            For real groups, this coincides with the number of
            reflecting hyperplanes.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.nr_reflections()
                3
                sage: W = ColoredPermutations(2,3)
                sage: W.nr_reflections()
                9
                sage: W = ColoredPermutations(4,3)
                sage: W.nr_reflections()
                21
                sage: W = ReflectionGroup((4,2,3))
                sage: W.nr_reflections()
                15
            """
            return sum(self.degrees()) - self.rank()

        def rank(self):
            r"""
            Return the rank of ``self``.

            This is the dimension of the underlying vector space.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.rank()
                2
                sage: W = ColoredPermutations(2,3)
                sage: W.rank()
                3
                sage: W = ColoredPermutations(4,3)
                sage: W.rank()
                3
                sage: W = ReflectionGroup((4,2,3))
                sage: W.rank()
                3
            """
            return len(self.degrees())

        def nr_irreducible_components(self):
            r"""
            Return the number of irreducible components of ``self``.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.nr_irreducible_components()
                1

                sage: W = ReflectionGroup((1,1,3),(2,1,3))
                sage: W.nr_irreducible_components()
                2
            """
            return len(self.irreducible_components())

        def cardinality(self):
            r"""
            Return the cardinality of ``self``.

            It is given by the product of the degrees of ``self``.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.cardinality()
                6
                sage: W = ColoredPermutations(2,3)
                sage: W.cardinality()
                48
                sage: W = ColoredPermutations(4,3)
                sage: W.cardinality()
                384
                sage: W = ReflectionGroup((4,2,3))
                sage: W.cardinality()
                192
            """
            return ZZ.prod(self.degrees())

        def is_well_generated(self):
            r"""
            Return ``True`` if ``self`` is well-generated.

            This is, if ``self`` is generated by `n` many
            reflections where `n` is the rank of ``self``.

            .. NOTE::

                - All finite real reflection groups are well generated.
                - The complex reflection groups of type `G(r,1,n)` and
                  of type `G(r,r,n)` are well generated.
                - The complex reflection groups of type `G(r,p,n)`
                  with `1 < p < r` are *not* well generated.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.is_well_generated()
                True

                sage: W = ColoredPermutations(4,3)
                sage: W.is_well_generated()
                True

                sage: W = ReflectionGroup((4,2,3))
                sage: W.is_well_generated()
                False

                sage: W = ReflectionGroup((4,4,3))
                sage: W.is_well_generated()
                True
            """
            return self.nr_simple_reflections() == self.rank()

        def is_real(self):
            r"""
            Return ``True`` if ``self`` is real.

            For irreducible reflection groups, this holds if and
            only if `2` is a degree of ``self``.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3)
                sage: W.is_real()
                True

                sage: W = ColoredPermutations(4,3)
                sage: W.is_real()
                False
            """
            return self.degrees().count(2) == self.nr_irreducible_components()

    class ElementMethods:

        @abstract_method(optional=True)
        def to_matrix(self):
            r"""
            Return the matrix presentation of ``self`` acting on a
            vector space `V`.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: [t.to_matrix() for t in W]
                [
                [1 0]  [-1  0]  [ 1  1]  [-1 -1]  [ 0  1]  [ 0 -1]
                [0 1], [ 1  1], [ 0 -1], [ 1  0], [-1 -1], [-1  0]
                ]

                sage: W = ColoredPermutations(1,3)
                sage: [t.to_matrix() for t in W]
                [
                [1 0 0]  [1 0 0]  [0 1 0]  [0 0 1]  [0 1 0]  [0 0 1]
                [0 1 0]  [0 0 1]  [1 0 0]  [1 0 0]  [0 0 1]  [0 1 0]
                [0 0 1], [0 1 0], [0 0 1], [0 1 0], [1 0 0], [1 0 0]
                ]

            A different representation is given by the
            colored permutations::

                sage: W = ColoredPermutations(3, 1)
                sage: [t.to_matrix() for t in W]
                [[1], [zeta3], [-zeta3 - 1]]
            """

        def _matrix_(self):
            """
            Return ``self`` as a matrix.

            EXAMPLES::

                sage: W = ReflectionGroup((1,1,3))
                sage: [matrix(t) for t in W]
                [
                [1 0]  [-1  0]  [ 1  1]  [-1 -1]  [ 0  1]  [ 0 -1]
                [0 1], [ 1  1], [ 0 -1], [ 1  0], [-1 -1], [-1  0]
                ]
            """
            return self.to_matrix()

        def character_value(self):
            r"""
            Return the value at ``self`` of the character of the
            reflection representation given by :meth:`to_matrix`.

            EXAMPLES::

                sage: W = ColoredPermutations(1,3); W
                1-colored permutations of size 3
                sage: [t.character_value() for t in W]
                [3, 1, 1, 0, 0, 1]

            Note that this could be a different (faithful)
            representation that given by the corresponding
            root system::

                sage: W = ReflectionGroup((1,1,3)); W
                Irreducible complex reflection group of rank 2 and type A2
                sage: [t.character_value() for t in W]
                [2, 0, 0, -1, -1, 0]

                sage: W = ColoredPermutations(2,2); W
                2-colored permutations of size 2
                sage: [t.character_value() for t in W]
                [2, 0, 0, -2, 0, 0, 0, 0]

                sage: W = ColoredPermutations(3,1); W
                3-colored permutations of size 1
                sage: [t.character_value() for t in W]
                [1, zeta3, -zeta3 - 1]
            """
            return self.to_matrix().trace()

    class Irreducible(CategoryWithAxiom):

        def example(self):
            r"""
            Return an example of an irreducible complex
            reflection group.

            EXAMPLES::

                sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
                sage: ComplexReflectionGroups().Finite().Irreducible().example()
                Irreducible complex reflection group of rank 3 and type G(4,2,3)
            """
            from sage.combinat.root_system.reflection_group_real import ReflectionGroup
            return ReflectionGroup((4,2,3))

        class ParentMethods:
            def coxeter_number(self):
                r"""
                Return the Coxeter number of an irreducible
                reflection group.

                This is defined as `\frac{N + N^*}{n}` where
                `N` is the number of reflections, `N^*` is the
                number of reflecting hyperplanes, and `n` is the
                rank of ``self``.

                EXAMPLES::

                    sage: W = ReflectionGroup(31)
                    sage: W.coxeter_number()
                    30
                """
                return (self.nr_reflecting_hyperplanes() + self.nr_reflections()) // self.rank()

    class WellGenerated(CategoryWithAxiom):
        r"""
        The category of finite well-generated finite complex
        reflection groups.

        This is the base category for finite subgroups of the
        special linear group which are generated by `n`
        reflections where `n` is the finite dimension of the
        underlying vector space.

        An example of a finite well-generated complex reflection
        group::

            sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
            sage: W = ComplexReflectionGroups().Finite().WellGenerated().example(); W
            Reducible complex reflection group of rank 4 and type A2 x G(3,1,2)

        ``W`` is in the category of complex reflection groups::

            sage: W in ComplexReflectionGroups().Finite().WellGenerated()
            True

        TESTS::

            sage: TestSuite(W).run()
            sage: TestSuite(ComplexReflectionGroups().Finite().WellGenerated()).run()
        """
        def example(self):
            r"""
            Return an example of a well-generated complex
            reflection group.

            EXAMPLES::

                sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
                sage: ComplexReflectionGroups().Finite().WellGenerated().example()
                Reducible complex reflection group of rank 4 and type A2 x G(3,1,2)
            """
            from sage.combinat.root_system.reflection_group_real import ReflectionGroup
            return ReflectionGroup((1,1,3), (3,1,2))

        class Irreducible(CategoryWithAxiom):
            r"""
            The category of finite irreducible well-generated
            finite complex reflection groups.
            """
            def example(self):
                r"""
                Return an example of an irreducible well-generated
                complex reflection group.

                EXAMPLES::

                    sage: from sage.categories.complex_reflection_groups import ComplexReflectionGroups
                    sage: ComplexReflectionGroups().Finite().WellGenerated().Irreducible().example()
                    4-colored permutations of size 3
                """
                from sage.combinat.colored_permutations import ColoredPermutations
                return ColoredPermutations(4, 3)

            class ParentMethods:
                def coxeter_number(self):
                    r"""
                    Return the Coxeter number of a well-generated,
                    irreducible reflection group. This is defined to be
                    the order of a regular element in ``self``, and is
                    equal to the highest degree of self.

                    .. NOTE::

                        This method overwrites the more general
                        method for complex reflection groups since
                        the expression given here is quicker to
                        compute.

                    EXAMPLES::

                        sage: W = ColoredPermutations(1,3)
                        sage: W.coxeter_number()
                        3

                        sage: W = ColoredPermutations(4,3)
                        sage: W.coxeter_number()
                        12

                        sage: W = ReflectionGroup((4,4,3))
                        sage: W.coxeter_number()
                        8
                    """
                    return max(self.degrees())

                def number_of_reflections_of_full_support(self):
                    r"""
                    Return the number of reflections with full
                    support.

                    EXAMPLES::

                        sage: W = ColoredPermutations(1,4)
                        sage: W.number_of_reflections_of_full_support()
                        1

                        sage: W = ColoredPermutations(3,3)
                        sage: W.number_of_reflections_of_full_support()
                        3
                    """
                    n = self.rank()
                    h = self.coxeter_number()
                    l = self.cardinality()
                    codegrees = self.codegrees()[:-1]
                    return n * h * prod(codegrees) // l

                @cached_method
                def rational_catalan_number(self, p, polynomial=False):
                    r"""
                    Return the ``p``-th rational Catalan number
                    associated to ``self``.

                    It is defined by

                    .. MATH::

                        \prod_{i = 1}^n \frac{p + (p(d_i-1)) \mod h)}{d_i},

                    where `d_1, \ldots, d_n` are the degrees and
                    `h` is the Coxeter number. See [STW2016]_
                    for this formula.

                    REFERENCES:

                    .. [STW2016] C. Stump, H. Thomas, N. Williams.
                       *Cataland II*, in preparation, 2016.

                    EXAMPLES::

                        sage: W = ColoredPermutations(1,3)
                        sage: [ W.rational_catalan_number(p) for p in [5,7,8] ]
                        [7, 12, 15]

                        sage: W = ColoredPermutations(2,2)
                        sage: [ W.rational_catalan_number(p) for p in [7,9,11] ]
                        [10, 15, 21]
                    """
                    from sage.rings.all import ZZ
                    from sage.arith.all import gcd
                    from sage.combinat.q_analogues import q_int

                    h = self.coxeter_number()
                    if not gcd(h,p) == 1:
                        raise ValueError("parameter p = %s is not coprime to the Coxeter number %s"%(p,h))

                    if polynomial:
                        f = q_int
                    else:
                        f = lambda n: n

                    num = prod(f(p + (p*(deg-1))%h) for deg in self.degrees())
                    den = prod(f(deg              ) for deg in self.degrees())
                    ret = num / den
                    if ret in ZZ:
                        ret = ZZ(ret)
                    return ret

                def fuss_catalan_number(self, m, positive=False, polynomial=False):
                    r"""
                    Return the ``m``-th Fuss-Catalan number
                    associated to ``self``.

                    It is defined by

                    .. MATH::

                        \prod_{i = 1}^n \frac{d_i + mh}{d_i},

                    where `d_1, \ldots, d_n` are the degrees and
                    `h` is the Coxeter number.
                    See [Arm2006]_ for further information.

                    .. NOTE::

                        - For the symmetric group `S_n`, it reduces to the
                          Fuss-Catalan number `\frac{1}{mn+1}\binom{(m+1)n}{n}`.
                        - The Fuss-Catalan numbers for `G(r, 1, n)` all
                          coincide for `r > 1`.

                    REFERENCES:

                    .. [Arm2006] D. Armstrong. *Generalized noncrossing
                       partitions and combinatorics of Coxeter groups*.
                       Mem. Amer. Math. Soc., 2006.

                    EXAMPLES::

                        sage: W = ColoredPermutations(1,3)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [5, 12, 22]

                        sage: W = ColoredPermutations(1,4)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [14, 55, 140]

                        sage: W = ColoredPermutations(1,5)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [42, 273, 969]

                        sage: W = ColoredPermutations(2,2)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [6, 15, 28]

                        sage: W = ColoredPermutations(2,3)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [20, 84, 220]

                        sage: W = ColoredPermutations(2,4)
                        sage: [ W.fuss_catalan_number(i) for i in [1,2,3] ]
                        [70, 495, 1820]
                    """
                    h = self.coxeter_number()
                    if positive:
                        p = m*h-1
                    else:
                        p = m*h+1

                    return self.rational_catalan_number(p, polynomial=polynomial)

                def catalan_number(self,positive=False,polynomial=False):
                    r"""
                    Return the Catalan number associated to ``self``.

                    It is defined by

                    .. MATH::

                        \prod_{i = 1}^n \frac{d_i + h}{d_i},

                    where `d_1, \ldots, d_n` are the degrees and where
                    `h` is the Coxeter number.
                    See [Arm2006]_ for further information.

                    .. NOTE::

                        - For the symmetric group `S_n`, it reduces to the
                          Catalan number `\frac{1}{n+1} \binom{2n}{n}`.
                        - The Catalan numbers for `G(r,1,n)` all coincide
                          for `r > 1`.

                    EXAMPLES::

                        sage: [ColoredPermutations(1,n).catalan_number() for n in [3,4,5]]
                        [5, 14, 42]

                        sage: [ColoredPermutations(2,n).catalan_number() for n in [3,4,5]]
                        [20, 70, 252]

                        sage: [ReflectionGroup((2,2,n)).catalan_number() for n in [3,4,5]]
                        [14, 50, 182]
                    """
                    return self.fuss_catalan_number(1, positive=positive,
                                                    polynomial=polynomial)

