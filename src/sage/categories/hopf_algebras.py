r"""
Hopf algebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from category import Category
from category_types import Category_over_base_ring
from sage.categories.bialgebras import Bialgebras
from sage.categories.tensor import TensorProductsCategory # tensor
from sage.categories.realizations import RealizationsCategory
from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute

class HopfAlgebras(Category_over_base_ring):
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

    WithBasis = LazyImport('sage.categories.hopf_algebras_with_basis',  'HopfAlgebrasWithBasis')

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


    class TensorProducts(TensorProductsCategory):
        """
        The category of Hopf algebras constructed by tensor product of Hopf algebras
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: C = HopfAlgebras(QQ).TensorProducts()
                sage: C.extra_super_categories()
                [Category of hopf algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of hopf algebras over Rational Field,
                 Category of tensor products of algebras over Rational Field,
                 Category of tensor products of coalgebras over Rational Field]
            """
            return [self.base_category()]

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

    class Realizations(RealizationsCategory):

        class ParentMethods:

            # TODO:
            # - Use @conditionally_defined once it's in Sage, for a nicer idiom
            # - Do the right thing (TM): once we will have proper
            #   overloaded operators (as in MuPAD-Combinat; see #8900),
            #   we won't need to specify explicitly to which parent one
            #   should coerce the input to calculate the antipode; so it
            #   will be sufficient to put this default implementation in
            #   HopfAlgebras.ParentMethods.
            def antipode_by_coercion(self, x):
                """
                Returns the image of ``x`` by the antipode

                This default implementation coerces to the default
                realization, computes the antipode there, and coerces the
                result back.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: R.antipode_by_coercion.__module__
                    'sage.categories.hopf_algebras'
                    sage: R.antipode_by_coercion(R[1,3,1])
                    -R[2, 1, 2]
                """
                R = self.realization_of().a_realization()
                return self(R(x).antipode())
