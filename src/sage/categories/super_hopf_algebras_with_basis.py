r"""
Super Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory

class SuperHopfAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super Hopf algebras with a distinguished basis.

    EXAMPLES::

        sage: C = HopfAlgebras(ZZ).WithBasis().Super(); C
        Category of super hopf algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of super algebras with basis over Integer Ring,
         Category of super coalgebras with basis over Integer Ring,
         Category of super hopf algebras over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """

