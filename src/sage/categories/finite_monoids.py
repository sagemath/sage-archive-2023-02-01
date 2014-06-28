r"""
Finite Monoids
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom

class FiniteMonoids(CategoryWithAxiom):
    """
    The category of finite (multiplicative) :class:`monoids <Monoids>`.

    A finite monoid is a :class:`finite sets <FiniteSets>` endowed
    with an associative unital binary operation `*`.

    EXAMPLES::

        sage: FiniteMonoids()
        Category of finite monoids
        sage: FiniteMonoids().super_categories()
        [Category of monoids, Category of finite semigroups]

    TESTS::

        sage: TestSuite(FiniteMonoids()).run()
    """

    class ElementMethods:
        def pseudo_order(self):
            r"""
            Returns the pair `[k, j]` with `k` minimal and `0\leq j <k` such
            that ``self^k == self^j``.

            Note that `j` is uniquely determined.

            EXAMPLES::

                sage: M = FiniteMonoids().example(); M
                An example of a finite multiplicative monoid: the integers modulo 12

                sage: x = M(2)
                sage: [ x^i for i in range(7) ]
                [1, 2, 4, 8, 4, 8, 4]
                sage: x.pseudo_order()
                [4, 2]

                sage: x = M(3)
                sage: [ x^i for i in range(7) ]
                [1, 3, 9, 3, 9, 3, 9]
                sage: x.pseudo_order()
                [3, 1]

                sage: x = M(4)
                sage: [ x^i for i in range(7) ]
                [1, 4, 4, 4, 4, 4, 4]
                sage: x.pseudo_order()
                [2, 1]

                sage: x = M(5)
                sage: [ x^i for i in range(7) ]
                [1, 5, 1, 5, 1, 5, 1]
                sage: x.pseudo_order()
                [2, 0]

            TODO: more appropriate name? see, for example, Jean-Eric Pin's
            lecture notes on semigroups.
            """
            self_powers = {self.parent().one(): 0}
            k = 1
            self_power_k = self
            while not self_power_k in self_powers:
                self_powers[self_power_k] = k
                k += 1
                self_power_k = self_power_k * self
            return [k, self_powers[self_power_k]]
