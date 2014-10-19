r"""
Graded modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory

class GradedModulesWithBasis(GradedModulesCategory):
    """
    The category of graded modules with a distinguished basis.

    EXAMPLES::

        sage: C = GradedModulesWithBasis(ZZ); C
        Category of graded modules with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered modules with basis over Integer Ring,
         Category of graded modules over Integer Ring]
        sage: C is ModulesWithBasis(ZZ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        r"""
        Adds :class:`FilteredModulesWithBasis` to the super categories
        of ``self`` since every graded module admits a filtraion.

        EXAMPLES::

            sage: GradedModulesWithBasis(ZZ).extra_super_categories()
            [Category of filtered modules with basis over Integer Ring]
        """
        from sage.categories.filtered_modules_with_basis import FilteredModulesWithBasis
        return [FilteredModulesWithBasis(self.base_ring())]

    class ParentMethods:
        pass

    class ElementMethods:
        pass

