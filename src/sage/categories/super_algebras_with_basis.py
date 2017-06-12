r"""
Super algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.algebras import Algebras
from sage.categories.modules import Modules

class SuperAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super algebras with a distinguished basis

    EXAMPLES::

        sage: C = Algebras(ZZ).WithBasis().Super(); C
        Category of super algebras with basis over Integer Ring

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: C = Algebras(ZZ).WithBasis().Super()
            sage: sorted(C.super_categories(), key=str) # indirect doctest
            [Category of graded algebras with basis over Integer Ring,
             Category of super algebras over Integer Ring,
             Category of super modules with basis over Integer Ring]
        """
        return [self.base_category().Graded()]

    class ParentMethods:
        def graded_algebra(self):
            r"""
            Return the associated graded module to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            .. SEEALSO::

                :meth:`~sage.categories.filtered_modules_with_basis.ParentMethods.graded_algebra`

            EXAMPLES::

                sage: W.<x,y> = algebras.DifferentialWeyl(QQ)
                sage: W.graded_algebra()
                Graded Algebra of Differential Weyl algebra of
                 polynomials in x, y over Rational Field
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)

