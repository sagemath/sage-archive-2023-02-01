"""
Cartesian product functorial construction
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category, CovariantFunctorialConstruction
from sage.misc.cachefunc import cached_method
import sage.structure.parent
import sage.structure.element

class CategoryWithCartesianProduct(Category):
    """
    A category with cartesian product is a category endowed with a cartesian product
    functor (operation on its parents and on its elements).

    Technically, let ``CClass`` be a class inheriting from
    ``CategoryWithCartesianProduct``. An instance `C` of ``CClass`` is a category.

    ``CClass`` must implement a method :meth:`.cartesian_product_category`
    which returns the category of cartesian products of parents in `C`. With
    the default implementation of :meth:`.cartesian_product_category`, it is
    sufficient to provide a class ``CClass.CartesianProductCategory`` whose
    constructor takes as parameter the category `C` and returns the
    desired category.

    If `C` is a subcategory of another category with cartesian product `D`,
    ``C.cartesian_product_category()`` is automatically considered as a
    subcategory of `D.cartesian_product_category()``.

    See also :class:`CovariantFunctorialConstruction`.

    TESTS::

        sage: TestSuite(CategoryWithCartesianProduct()).run() # mostly to silence sage -coverage on this abstract class
    """

    @cached_method
    def cartesian_product_category(self):
        """
        The category of cartesian products of parents in ``self``

        EXAMPLES::

            sage: ModulesWithBasis(QQ).cartesian_product_category()
            Category of cartesian products of modules with basis over Rational Field
        """
        return self.CartesianProductCategory(self)

    class ParentMethods:
        def cartesian_product(*parents):
            """
            Returns the cartesian product of the parents

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example(); A.rename("A")
                sage: A.cartesian_product(A,A)
                A (+) A (+) A
            """

            return parents[0].CartesianProduct(parents, category = cartesian_product.category_from_parents(parents))

    class ElementMethods:
        def cartesian_product(*elements):
            """
            Returns the cartesian product of its arguments, as an element of
            the cartesian product of the parents of those elements.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example()
                sage: (a,b,c) = A.algebra_generators()
                sage: a.cartesian_product(b, c)
                B[(0, word: a)] + B[(1, word: b)] + B[(2, word: c)]

            FIXME: is this a policy that we want to enforce on all parents?
            """
            assert(all(isinstance(element, sage.structure.element.Element) for element in elements))
            parents = [element.parent() for element in elements]
            return cartesian_product(parents)._cartesian_product_of_elements(elements) # good name???


class CartesianProductCategory(CategoryWithCartesianProduct):
    """
    An abstract base class for all CartesianProductCategory's defined in
    CategoryWithCartesianProduct's.
    """

    def __init__(self, category, name=None):
        """
        TESTS::

            sage: C = sage.categories.cartesian_product.CartesianProductCategory(ModulesWithBasis(QQ), "my_category")
            sage: C
            Category of cartesian products of modules with basis over Rational Field
            sage: C.base_category
            Category of modules with basis over Rational Field
            sage: latex(C)
            \mathbf{my_category}
        """
        Category.__init__(self, name)
        self.base_category = category

    def _repr_(self):
        """
        EXAMPLES::

            sage: ModulesWithBasis(QQ).cartesian_product_category() # indirect doctest
            Category of cartesian products of modules with basis over Rational Field

        """
        return "Category of cartesian products of %s"%repr(self.base_category)[12:]

    def cartesian_product_category(self):
        """
        Returns the category of cartesian products of ``self``

        By associativity of cartesian products, this is ``self`` (a
        cartesian_product of cartesian_products of `A`'s is a cartesian product of `A`'s)

        EXAMPLES::

            sage: ModulesWithBasis(QQ).cartesian_product_category().cartesian_product_category()
            Category of cartesian products of modules with basis over Rational Field
        """
        return self


class CartesianProductFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the cartesian_product functor
    """
    functor_name = "cartesian_product"
    functor_category = "cartesian_product_category"
    FunctorialCategory = CategoryWithCartesianProduct
    symbol = " (+) "

cartesian_product = CartesianProductFunctor()
