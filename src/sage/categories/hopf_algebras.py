r"""
HopfAlgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Category, CategoryWithTensorProduct, TensorCategory, DualityCategory, Bialgebras
#from sage.categories.tensor import tensor
from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute

class HopfAlgebras(Category_over_base_ring, CategoryWithTensorProduct, DualityCategory):
    """
    The category of Hopf algebras

    EXAMPLES::

        sage: HopfAlgebras(QQ)
        Category of hopf algebras over Rational Field
        sage: HopfAlgebras(QQ).super_categories()
        [Category of bialgebras over Rational Field]

    TESTS::

        sage: TestSuite(HopfAlgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: HopfAlgebras(QQ).super_categories()
            [Category of bialgebras over Rational Field]
        """
        R = self.base_ring()
        return [Bialgebras(R)]

    def dual(self):
        """
        Returns the dual category

        EXAMPLES:

        The category of Hopf algebras over any field is self dual::

            sage: C = HopfAlgebras(QQ)
            sage: C.dual()
            Category of hopf algebras over Rational Field
        """
        return self

    class ElementMethods:
        def antipode(self):
            """
            Returns the antipode of self.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.antipode()
                (B[(1,2,3)], B[(1,3,2)])
                sage: b, b.antipode()
                (B[(1,3)], B[(1,3)])

            TESTS::

                sage: all(x.antipode() * x == A.one() for x in A.basis())
                True
            """
            return self.parent().antipode(self)
            # Variant: delegates to the overloading mechanism
            # result not guaranted to be in self
            # This choice should be done consistently with coproduct, ...
            # return operator.antipode(self)

    class ParentMethods:
        #def __setup__(self): # Check the conventions for _setup_ or __setup__
        #    if self.implements("antipode"):
        #        coercion.declare(operator.antipode, [self], self.antipode)
        #
        #@lazy_attribute
        #def antipode(self):
        #    # delegates to the overloading mechanism but
        #    # guarantees that the result is in self
        #    compose(self, operator.antipode, domain=self)
        pass

    class Morphism(Category):
        """
        The category of Hopf algebra morphisms
        """
        pass


    class TensorCategory(TensorCategory):
        """
        The category of Hopf algebras constructed by tensor product of Hopf algebras
        """
        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: HopfAlgebras(QQ).tensor_category().super_categories()
                [Category of hopf algebras over Rational Field]
            """
            return [HopfAlgebras(self.base_category.base_ring())]

        class ParentMethods:
            # TODO: enable when tensor product of morphisms will be implemented
            #@lazy_attribute
            #def antipode(self):
            #    return tensor([module.antipode for module in self.modules])
            pass

        class ElementMethods:
            pass

    class DualCategory(Category_over_base_ring):
        """
        The category of Hopf algebras constructed as dual of a Hopf algebra
        """

        class ParentMethods:
            #@lazy_attribute
            #def antipode(self):
            #    self.dual().antipode.dual() # Check that this is the correct formula
            pass
