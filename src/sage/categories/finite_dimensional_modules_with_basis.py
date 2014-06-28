r"""
Finite dimensional modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class FiniteDimensionalModulesWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of finite dimensional modules with a distinguished basis

    EXAMPLES::

      sage: C = FiniteDimensionalModulesWithBasis(ZZ); C
      Category of finite dimensional modules with basis over Integer Ring
      sage: sorted(C.super_categories(), key=str)
      [Category of finite dimensional modules over Integer Ring,
       Category of modules with basis over Integer Ring]
      sage: C is Modules(ZZ).WithBasis().FiniteDimensional()
      True

    TESTS::

        sage: TestSuite(C).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass
