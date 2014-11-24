r"""
Graded algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.misc.cachefunc import cached_method

class GradedAlgebrasWithBasis(GradedModulesCategory):
    """
    The category of graded algebras with a distinguished basis

    EXAMPLES::

        sage: C = GradedAlgebrasWithBasis(ZZ); C
        Category of graded algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered algebras with basis over Integer Ring,
         Category of graded algebras over Integer Ring,
         Category of graded modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        # This needs to be copied in GradedAlgebras because we need to have
        #   FilteredAlgebrasWithBasis as an extra super category
        @cached_method
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            This is ``self``, because ``self`` is already graded.
            See :meth:`~sage.categories.filtered_algebras_with_basis.FilteredAlgebrasWithBasis.graded_algebra`
            for the general behavior of this method, and see
            :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and properties of associated graded
            algebras.

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).m()
                sage: m.graded_algebra() is m
                True

            .. TODO::

                Add examples showing that the three methods are
                overridden correctly.
            """
            return self

        # .. TODO::
        #     Possibly override ``to_graded_conversion`` and
        #     ``from_graded_conversion`` with identity morphisms?
        #     I have to admit I don't know the right way to construct
        #     identity morphisms other than using the identity matrix.
        #     Also, ``projection`` could be overridden by, well, a
        #     projection.

    class ElementMethods:
        pass

