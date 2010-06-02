r"""
Modules
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.all import Bimodules, HomCategory, Fields
from category_types import Category_module
from sage.misc.cachefunc import cached_method

class Modules(Category_module):
    r"""
    The category of all modules over a base ring `R`

    A `R`-module `M` is a left and right `R`-module over a commutative
    ring `R` such that:

    .. math::  r*(x*s) = (r*x)*s \qquad  \forall r,s \in R \text{ and } x\in M

    INPUT:

      - ``base_ring`` -- a ring `R`
      - ``dispatch`` -- a boolean (for internal use; default: ``True``)

    When the base ring is a field, the category of vector spaces is
    returned instead (unless ``dispatch == False``).

    EXAMPLES::

        sage: Modules(IntegerRing())
        Category of modules over Integer Ring
        sage: Modules(RationalField())
        Category of vector spaces over Rational Field

        sage: Modules(Integers(9))
        Category of modules over Ring of integers modulo 9

        sage: Modules(Integers(9)).super_categories()
        [Category of bimodules over Ring of integers modulo 9 on the left and Ring of integers modulo 9 on the right]
        sage: Modules(Integers(9)).all_super_categories()
        [Category of modules over Ring of integers modulo 9,
         Category of bimodules over Ring of integers modulo 9 on the left and Ring of integers modulo 9 on the right,
         Category of left modules over Ring of integers modulo 9,
         Category of right modules over Ring of integers modulo 9,
         Category of commutative additive groups,
         Category of commutative additive monoids,
         Category of commutative additive semigroups,
         Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

        sage: Modules(ZZ).super_categories()
        [Category of bimodules over Integer Ring on the left and Integer Ring on the right]

        sage: Modules == RingModules
        True

        sage: Modules(ZZ[x]).is_abelian()   # see #6081
        True

    TESTS::

        sage: TestSuite(Modules(ZZ)).run()

    TODO:

     - Implement a FreeModules(R) category, when so prompted by a concrete use case
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, dispatch = True):
        """
        This method implements the default behavior of dispatching
        ``Modules(field)`` to ``VectorSpaces(field)``. This feature
        will later be extended to modules over a principal ideal
        domain/ring or over a semiring.

        TESTS::

            sage: C = Modules(ZZ); C
            Category of modules over Integer Ring
            sage: C is Modules(ZZ, dispatch = False)
            True
            sage: C is Modules(ZZ, dispatch = True)
            True
            sage: C._reduction
            (<class 'sage.categories.modules.Modules'>, (Integer Ring,), {'dispatch': False})
            sage: TestSuite(C).run()

            sage: Modules(QQ) is VectorSpaces(QQ)
            True
            sage: Modules(QQ, dispatch = True) is VectorSpaces(QQ)
            True

            sage: C = Modules(NonNegativeIntegers()); C   # todo: not implemented
            Category of semiring modules over Non negative integers

            sage: C = Modules(QQ, dispatch = False); C
            Category of modules over Rational Field
            sage: C._reduction
            (<class 'sage.categories.modules.Modules'>, (Rational Field,), {'dispatch': False})
            sage: TestSuite(C).run()

        """
        if dispatch and base_ring in Fields():
            from vector_spaces import VectorSpaces
            return VectorSpaces(base_ring)
        result = super(Modules, cls).__classcall__(cls, base_ring)
        result._reduction[2]['dispatch'] = False
        return result

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Modules(ZZ).super_categories()
            [Category of bimodules over Integer Ring on the left and Integer Ring on the right]

        Nota bene::

            sage: Modules(QQ)
            Category of vector spaces over Rational Field
            sage: Modules(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_ring()
        return [Bimodules(R,R)]

    class ParentMethods:
        pass

    class ElementMethods:

        def __mul__(left, right):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
                sage: x = F.monomial("a")
                sage: x * int(2)
                2*B['a']

            TODO: make a better unit test once Modules().example() is implemented
            """
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, right, operator.mul)

        def __rmul__(right, left):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(QQ, ["a", "b"])
                sage: x = F.monomial("a")
                sage: int(2) * x
                2*B['a']

            TODO: make a better unit test once Modules().example() is implemented
            """
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, right, operator.mul)


    class HomCategory(HomCategory):
        """
        The category of homomorphisms sets `\hom(X,Y)` for `X`, `Y` modules
        """

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Modules(ZZ).hom_category().extra_super_categories()
                [Category of modules over Integer Ring]
            """
            return [Modules(self.base_category.base_ring())]

        class ParentMethods:
            @cached_method
            def zero(self):
                """
                EXAMPLES::

                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: H = Hom(E, F)
                    sage: f = H.zero()
                    sage: f
                    Generic morphism:
                      From: Free module generated by {1, 2, 3} over Integer Ring
                      To:   Free module generated by {2, 3, 4} over Integer Ring
                    sage: f(E.monomial(2))
                    0
                    sage: f(E.monomial(3)) == F.zero()
                    True
                """
                return self(lambda x: self.codomain().zero())

    class EndCategory(HomCategory):
        """
        The category of endomorphisms sets `End(X)` for `X` module (this is
        not used yet)
        """

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Hom(ZZ^3, ZZ^3).category().extra_super_categories() # todo: not implemented
                [Category of algebras over Integer Ring]
            """
            from algebras import Algebras
            return [Algebras(self.base_category.base_ring())]
