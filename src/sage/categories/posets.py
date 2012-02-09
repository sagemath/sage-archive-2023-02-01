r"""
Posets
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.categories.category import Category
from sage.categories.sets_cat import Sets

class Posets(Category):
    r"""
    The category of posets i.e. sets with a partial order structure.

    EXAMPLES::

        sage: Posets()
        Category of posets
        sage: Posets().super_categories()
        [Category of sets]
        sage: P = Posets().example(); P
        An example of a poset: sets ordered by inclusion

    The partial order is implemented by the mandatory method
    :meth:`~Posets.ParentMethods.le`::

        sage: x = P(Set([1,3])); y = P(Set([1,2,3]))
        sage: x, y
        ({1, 3}, {1, 2, 3})
        sage: P.le(x, y)
        True
        sage: P.le(x, x)
        True
        sage: P.le(y, x)
        False

    The other comparison methods are called
    :meth:`~Posets.ParentMethods.lt`, :meth:`~Posets.ParentMethods.ge`,
    :meth:`~Posets.ParentMethods.gt`, following Python's naming
    convention in :mod:`operator`. Default implementations are
    provided::

        sage: P.lt(x, x)
        False
        sage: P.ge(y, x)
        True

    Unless the poset is a facade (see :class:`Sets.Facades`), one can
    compare directly its elements using the usual Python operators::

        sage: D = Poset((divisors(30), attrcall("divides")))
        sage: D(3) <= D(6)
        True
        sage: D(3) <= D(3)
        True
        sage: D(3) <= D(5)
        False
        sage: D(3) < D(3)
        False
        sage: D(10) >= D(5)
        True

    At this point, this has to be implemented by hand. Once #10130
    will be resolved, this will be automatically provided by this category::

        sage: x < y      # todo: not implemented
        True
        sage: x < x      # todo: not implemented
        False
        sage: x <= x     # todo: not implemented
        True
        sage: y >= x     # todo: not implemented
        True

    .. seealso:: :func:`Poset`, :class:`FinitePosets`, :class:`LatticePosets`

    TESTS::

        sage: C = Posets()
        sage: TestSuite(C).run()

    """
    @cached_method
    def super_categories(self):
        r"""
        Returns a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: Posets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def example(self, choice = None):
        r"""
        Returns examples of objects of ``Posets()``, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: Posets().example()
            An example of a poset: sets ordered by inclusion

            sage: Posets().example("facade")
            An example of a facade poset: the positive integers ordered by divisibility
        """
        from sage.categories.examples.posets import FiniteSetsOrderedByInclusion, PositiveIntegersOrderedByDivisibilityFacade
        if choice == "facade":
            return PositiveIntegersOrderedByDivisibilityFacade()
        else:
            return FiniteSetsOrderedByInclusion()

    def __iter__(self):
        r"""
        Iterator over representatives of the isomorphism classes of
        posets with finitely many vertices.

        .. warning:: this feature may become deprecated, since it does
           of course not iterate through all posets.

        EXAMPLES::

            sage: P = Posets()
            sage: it = iter(P)
            sage: for _ in range(10): print it.next();
            Finite poset containing 0 elements
            Finite poset containing 1 elements
            Finite poset containing 2 elements
            Finite poset containing 2 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 4 elements
        """
        from sage.combinat.posets.posets import FinitePosets_n
        n = 0
        while True:
            for P in FinitePosets_n(n):
                yield P
            n += 1


    class ParentMethods:

        @abstract_method
        def le(self, x, y):
            r"""
            Returns whether `x \le y` in this poset

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.le( 3, 6 )
                True
                sage: D.le( 3, 3 )
                True
                sage: D.le( 3, 5 )
                False
            """

        def lt(self, x, y):
            r"""
            Returns whether `x < y` in this poset

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`le`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.lt( 3, 6 )
                True
                sage: D.lt( 3, 3 )
                False
                sage: D.lt( 3, 5 )
                False
            """
            return self.le(x,y) and x != y

        def ge(self, x, y):
            r"""
            Returns whether `x < y` in this poset

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`le`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.ge( 6, 3 )
                True
                sage: D.ge( 3, 3 )
                True
                sage: D.ge( 3, 5 )
                False
            """
            return self.le(y,x)

        def gt(self, x, y):
            r"""
            Returns whether `x < y` in this poset

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`lt`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.gt( 3, 6 )
                False
                sage: D.gt( 3, 3 )
                False
                sage: D.gt( 3, 5 )
                False
            """
            return self.lt(y,x)

        @abstract_method(optional = True)
        def upper_covers(self, x):
            r"""
            Returns the upper covers of `x`, that is the elements `y`
            such that `x<y` and there exists no `z` such that `x<z<y`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.upper_covers(3)
                [6, 15]
            """

        @abstract_method(optional = True)
        def lower_covers(self, x):
            r"""
            Returns the lower covers of `x`, that is the elements `y`
            such that `y<x` and there exists no `z` such that `y<z<x`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.lower_covers(15)
                [3, 5]
            """

        @abstract_method(optional = True)
        def order_ideal(self, gens):
            r"""
            Returns the order ideal in ``self`` generated by ``gens``.

            EXAMPLES::

                sage: B = Posets.BooleanLattice(4)
                sage: B.order_ideal([7,10])
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
            """

        @abstract_method(optional = True)
        def order_filter(self,elements):
            r"""
            Returns the order filter generated by a list of elements.

            `I` is an order filter if, for any `x` in `I` and `y` such that
            `y \ge x`, then `y` is in `I`.

            EXAMPLES::

                sage: B = Posets.BooleanLattice(4)
                sage: B.order_filter([3,8])
                [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            """

        def principal_order_ideal(self, x):
            r"""
            Returns the order ideal generated by an element ``x``.

            EXAMPLES::

                sage: B = Posets.BooleanLattice(4)
                sage: B.principal_order_ideal(6)
                [0, 2, 4, 6]
            """
            return self.order_ideal([x])

        def principal_order_filter(self, x):
            r"""
            Returns the order filter generated by an element ``x``.

            EXAMPLES::

                sage: B = Posets.BooleanLattice(4)
                sage: B.principal_order_filter(2)
                [2, 3, 6, 7, 10, 11, 14, 15]
            """
            return self.order_filter([x])

    class ElementMethods:
        pass
        # TODO: implement x<y, x<=y, x>y, x>=y appropriately once #10130 is resolved
        #
        # def __le__(self, other):
        #     r"""
        #     Returns whether ``self`` is smaller or equal to ``other``
        #     in the poset.
        #
        #     EXAMPLES::
        #
        #         sage: P = Posets().example(); P
        #         An example of poset: sets ordered by inclusion
        #         sage: x = P(Set([1,3])); y = P(Set([1,2,3]))
        #         sage: x.__le__(y)
        #         sage: x <= y
        #     """
        #     return self.parent().le(self, other)
