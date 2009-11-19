"""
Tensor product functorial construction
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

class CategoryWithTensorProduct(Category):
    """
    A category with tensor product is a category endowed with a tensor product
    functor (operation on its parents and on its elements).

    Technically, let ``CClass`` be a class inheriting from
    ``CategoryWithTensorProduct``. An instance `C` of ``CClass`` is a category.

    ``CClass`` must implement a method
    :meth:`.tensor_product_category` which returns the category of
    tensor products of parents in `C`. With the default implementation
    of :meth:`.tensor_product_category`, it is sufficient to provide a
    class ``CClass.TensorProductCategory`` whose constructor takes as
    parameter the category `C` and returns the desired category.

    If `C` is a subcategory of another category with tensor product
    `D`, ``C.tensor_product_category()`` is automatically considered
    as a subcategory of `D.tensor_product_category()``.

    See also :class:`CovariantFunctorialConstruction`.

    TESTS::

        sage: TestSuite(CategoryWithTensorProduct()).run() # mostly to silence sage -coverage on this abstract class
    """

    @cached_method
    def tensor_category(self):
        """
        The category of tensor products of parents in ``self``

        EXAMPLES::

            sage: ModulesWithBasis(QQ).tensor_category()
            Category of tensor products of modules with basis over Rational Field
        """
        return self.TensorCategory(self)

    class ParentMethods:
        def tensor(*parents):
            """
            Returns the tensor product of the parents

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example(); A.rename("A")
                sage: A.tensor(A,A)
                A # A # A
            """
            return parents[0].Tensor(parents, category = tensor.category_from_parents(parents))

    class ElementMethods:
        def tensor(*elements):
            """
            Returns the tensor product of its arguments, as an element of
            the tensor product of the parents of those elements.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example()
                sage: (a,b,c) = A.algebra_generators()
                sage: a.tensor(b, c)
                B[word: a] # B[word: b] # B[word: c]

            FIXME: is this a policy that we want to enforce on all parents?
            """
            assert(all(isinstance(element, sage.structure.element.Element) for element in elements))
            parents = [element.parent() for element in elements]
            return tensor(parents)._tensor_of_elements(elements) # good name???


class TensorCategory(CategoryWithTensorProduct):
    """
    An abstract base class for all TensorCategory's defined in
    CategoryWithTensorProduct's.
    """

    def __init__(self, category, name=None):
        """
        TESTS::

            sage: C = sage.categories.tensor.TensorCategory(ModulesWithBasis(QQ), "my_category")
            sage: C
            Category of tensor products of modules with basis over Rational Field
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

            sage: ModulesWithBasis(QQ).tensor_category() # indirect doctest
            Category of tensor products of modules with basis over Rational Field

        """
        return "Category of tensor products of %s"%repr(self.base_category)[12:]

    def tensor_category(self):
        """
        Returns the category of tensor products of objects of ``self``

        By associativity of tensor products, this is ``self``
        (a tensor product of tensor products of A's is a tensor product of A's)

        EXAMPLES::

            sage: ModulesWithBasis(QQ).tensor_category().tensor_category()
            Category of tensor products of modules with basis over Rational Field

        """
        return self


class TensorFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the tensor functor
    """
    functor_name = "tensor"
    functor_category = "tensor_category"
    FunctorialCategory = CategoryWithTensorProduct
    symbol = " # "

tensor = TensorFunctor()
