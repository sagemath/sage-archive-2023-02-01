r"""
HopfAlgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.all import CategoryWithTensorProduct, TensorCategory, DualityCategory, HopfAlgebras, BialgebrasWithBasis
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

class HopfAlgebrasWithBasis(Category_over_base_ring, CategoryWithTensorProduct): #, DualityCategory):
    """
    The category of Hopf algebras with a distinguished basis

    EXAMPLES::

        sage: C = HopfAlgebrasWithBasis(QQ)
        sage: C
        Category of hopf algebras with basis over Rational Field
        sage: C.super_categories()
        [Category of bialgebras with basis over Rational Field, Category of hopf algebras over Rational Field]

    We now show how to use a simple hopf algebra, namely the group algebra of the dihedral group
    (see also AlgebrasWithBasis)::

        sage: A = C.example(); A
        An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
        sage: A.__custom_name = "A"
        sage: A.category()
        Category of hopf algebras with basis over Rational Field

        sage: A.one_basis()
        ()
        sage: A.one()
        B[()]

        sage: A.base_ring()
        Rational Field
        sage: A.basis().keys()
        Dihedral group of order 6 as a permutation group

        sage: [a,b] = A.algebra_generators()
        sage: a, b
        (B[(1,2,3)], B[(1,3)])
        sage: a^3, b^2
        (B[()], B[()])
        sage: a*b
        B[(1,2)]

        sage: A.product           # todo: not quite ...
        Generic endomorphism of A
        sage: A.product(b,b)
        B[()]

        sage: A.zero().coproduct()
        0
        sage: A.zero().coproduct().parent()
        A # A
        sage: a.coproduct()
        B[(1,2,3)] # B[(1,2,3)]

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.hopf_algebras_with_basis.MyGroupAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.combinat.free_module.MyGroupAlgebra_with_category.element_class'>

    Let us look at the code for implementing A::

        sage: A??                       # todo: not implemented

    TESTS::

        sage: TestSuite(A).run()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: HopfAlgebrasWithBasis(QQ).super_categories()
            [Category of bialgebras with basis over Rational Field, Category of hopf algebras over Rational Field]
        """
        R = self.base_ring()
        return [BialgebrasWithBasis(R), HopfAlgebras(R)]

    def example(self, G = None):
        """
        Returns an example of algebra with basis::

            sage: HopfAlgebrasWithBasis(QQ[x]).example()
            An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Univariate Polynomial Ring in x over Rational Field

        An other group can be specified as optional argument::

            sage: HopfAlgebrasWithBasis(QQ).example(SymmetricGroup(4))
            An example of Hopf algebra with basis: the group algebra of the SymmetricGroup(4) over Rational Field
        """
        from sage.categories.examples.hopf_algebras_with_basis import MyGroupAlgebra
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        if G is None:
            G = DihedralGroup(3)
        return MyGroupAlgebra(self.base_ring(), G)

#     This is only correct in the finite dimensional / graded case
#     def dual(self):
#         """
#         Returns the dual category

#         EXAMPLES:

#         The category of Hopf algebras over any field is self dual::

#             sage: C = HopfAlgebrasWithBasis(QQ)
#             sage: C.dual()
#             Category of hopf algebras with basis over Rational Field
#         """
#         return self

    class ParentMethods:

        @lazy_attribute
        def antipode(self):
            """
            If :meth:`.antipode_basis` is available, construct the
            antipode morphism from ``self`` to ``self`` by extending
            it by linearity

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(ZZ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Integer Ring
                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: [a,b] = A.algebra_generators()
                sage: a, A.antipode(a)
                (B[(1,2,3)], B[(1,3,2)])
                sage: b, A.antipode(b)
                (B[(1,3)], B[(1,3)])

            TESTS::

                sage: all(A.antipode(x) * x == A.one() for x in A.basis())
                True
            """
            if hasattr(self, "antipode_on_basis"):
                # Should give the information that this is an anti-morphism of algebra
                return self._module_morphism(self.antipode_on_basis, codomain = self)

    class ElementMethods:
        pass

    class TensorCategory(TensorCategory):
        """
        The category of hopf algebras with basis constructed by tensor product of hopf algebras with basis
        """

        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: HopfAlgebrasWithBasis(QQ).tensor_category().super_categories()
                [Category of hopf algebras with basis over Rational Field]
            """
            return [HopfAlgebrasWithBasis(self.base_category.base_ring())]

        class ParentMethods:
            # todo: antipode
            pass

        class ElementMethods:
            pass

