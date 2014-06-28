r"""
Algebras

AUTHORS:

- David Kohel & William Stein (2005): initial revision
- Nicolas M. Thiery (2008-2011): rewrote for the category framework
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.dual import DualObjectsCategory
from sage.categories.tensor import TensorProductsCategory
from sage.categories.associative_algebras import AssociativeAlgebras

class Algebras(CategoryWithAxiom_over_base_ring):
    r"""
    The category of associative and unital algebras over a given base ring.

    An associative and unital algebra over a ring `R` is a module over
    `R` which is itself a ring.

    .. WARNING::

        :class:`Algebras` will be eventually be replaced by
        :class:`.magmatic_algebras.MagmaticAlgebras`
        for consistency with e.g. :wikipedia:`Algebras` which assumes
        neither associativity nor the existence of a unit (see
        :trac:`15043`).

    .. TODO:: Should `R` be a commutative ring?

    EXAMPLES::

        sage: Algebras(ZZ)
        Category of algebras over Integer Ring
        sage: sorted(Algebras(ZZ).super_categories(), key=str)
        [Category of associative algebras over Integer Ring,
         Category of rings,
         Category of unital algebras over Integer Ring]

    TESTS::

        sage: TestSuite(Algebras(ZZ)).run()
    """
    _base_category_class_and_axiom = (AssociativeAlgebras, 'Unital')

    # For backward compatibility?
    def __contains__(self, x):
        """
        Membership testing

        EXAMPLES::

            sage: QQ[x] in Algebras(QQ)
            True

            sage: QQ^3 in Algebras(QQ)
            False
            sage: QQ[x] in Algebras(CDF)
            False
        """
        if super(Algebras, self).__contains__(x):
            return True
        from sage.rings.ring import Algebra
        return isinstance(x, Algebra) and x.base_ring() == self.base_ring()

    # def extra_super_categories(self):
    #     """
    #     EXAMPLES::

    #         sage: Algebras(ZZ).super_categories()
    #         [Category of associative algebras over Integer Ring, Category of rings]
    #     """
    #     R = self.base_ring()
    #     return [Rings()] # TODO: won't be needed when Rings() will be Rngs().Unital()

    Commutative = LazyImport('sage.categories.commutative_algebras', 'CommutativeAlgebras', at_startup=True)
    Graded      = LazyImport('sage.categories.graded_algebras',      'GradedAlgebras')
    WithBasis   = LazyImport('sage.categories.algebras_with_basis',  'AlgebrasWithBasis')

    class ElementMethods:
        # TODO: move the content of AlgebraElement here or higher in the category hierarchy
        def _div_(self, y):
            """
            Division by invertible elements

            # TODO: move in Monoids

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ).example()
                sage: x = C(2); x
                2*B[word: ]
                sage: y = C.algebra_generators().first(); y
                B[word: a]

                sage: y._div_(x)
                1/2*B[word: a]
                sage: x._div_(y)
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= B[word: a])
            """
            return self.parent().product(self, ~y)

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of algebras constructed as cartesian products of algebras

        This construction gives the direct product of algebras. See
        discussion on:

         - http://groups.google.fr/group/sage-devel/browse_thread/thread/35a72b1d0a2fc77a/348f42ae77a66d16#348f42ae77a66d16
         - http://en.wikipedia.org/wiki/Direct_product
        """
        def extra_super_categories(self):
            """
            A cartesian product of algebras is endowed with a natural
            algebra structure.

            EXAMPLES::

                sage: C = Algebras(QQ).CartesianProducts()
                sage: C.extra_super_categories()
                [Category of algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of Cartesian products of commutative additive groups,
                 Category of Cartesian products of distributive magmas and additive magmas,
                 Category of Cartesian products of semigroups,
                 Category of Cartesian products of unital magmas,
                 Category of algebras over Rational Field]
            """
            return [self.base_category()]


    class TensorProducts(TensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).TensorProducts().extra_super_categories()
                [Category of algebras over Rational Field]
                sage: Algebras(QQ).TensorProducts().super_categories()
                [Category of algebras over Rational Field]

            Meaning: a tensor product of algebras is an algebra
            """
            return [self.base_category()]

        class ParentMethods:
            #def coproduct(self):
            #    tensor products of morphisms are not yet implemented
            #    return tensor(module.coproduct for module in self.modules)
            pass

        class ElementMethods:
            pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = Algebras(QQ)
                sage: C.dual()
                Category of duals of algebras over Rational Field
                sage: C.dual().extra_super_categories()
                [Category of coalgebras over Rational Field]

            .. WARNING::

                This is only correct in certain cases (finite dimension, ...).
                See :trac:`15647`.
            """
            from sage.categories.coalgebras import Coalgebras
            return [Coalgebras(self.base_category().base_ring())]
