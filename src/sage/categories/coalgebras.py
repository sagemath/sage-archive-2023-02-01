r"""
Coalgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Modules, Algebras, CategoryWithTensorProduct, TensorCategory, DualityCategory, tensor
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method

class Coalgebras(Category_over_base_ring, CategoryWithTensorProduct, DualityCategory):
    """
    The category of coalgebras

    EXAMPLES::

        sage: Coalgebras(QQ)
        Category of coalgebras over Rational Field
        sage: Coalgebras(QQ).super_categories()
        [Category of modules over Rational Field]
        sage: Coalgebras(QQ).all_super_categories() # todo: update once the hierarchy is more appropriate
        [Category of coalgebras over Rational Field,
         Category of modules over Rational Field,
         Category of bimodules over Rational Field on the left and Rational Field on the right,
         Category of left modules over Rational Field,
         Category of right modules over Rational Field,
         Category of commutative additive groups,
         Category of commutative additive monoids,
         Category of commutative additive semigroups,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    TESTS::

        sage: TestSuite(Coalgebras(ZZ)).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Coalgebras(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        return [Modules(self.base_ring())]

    def dual(self):
        """
        Returns the dual category

        EXAMPLES:

        The category of coalgebras over the Rational Field is dual
        to the category of algebras over the same field::

            sage: C = Coalgebras(QQ)
            sage: C.dual()
            Category of algebras over Rational Field
        """
        return Algebras(self.base_ring())

    class ParentMethods:
        #def __init_add__(self): # The analogue of initDomainAdd
        #    # Will declare the coproduct of self to the coercion mechanism when it exists
        #    pass

        @cached_method
        def tensor_square(self):
            """
            Returns the tensor square of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: A.tensor_square()
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field # An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
            """
            return tensor([self, self])

        @abstract_method
        def coproduct(self, x):
            """
            Returns the coproduct of x.

            Eventually, there will be a default implementation,
            delegating to the overloading mechanism and forcing the
            conversion back

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, A.coproduct(a)
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, A.coproduct(b)
                (B[(1,3)], B[(1,3)] # B[(1,3)])
            """
            #return self.tensor_square()(overloaded_coproduct(x))

    class ElementMethods:
        def coproduct(self):
            """
            Returns the coproduct of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.coproduct()
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, b.coproduct()
                (B[(1,3)], B[(1,3)] # B[(1,3)])
            """
            return self.parent().coproduct(self)

    class TensorCategory(TensorCategory):
        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: Coalgebras(QQ).tensor_category().super_categories()
                [Category of coalgebras over Rational Field]
            """
            return [Coalgebras(self.base_category.base_ring())]

        class ParentMethods:
            #def coproduct(self):
            #    tensor products of morphisms are not yet implemented
            #    return tensor(module.coproduct for module in self.modules)
            pass

        class ElementMethods:
            pass
