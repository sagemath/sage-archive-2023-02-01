r"""
Filtered Algebras
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.filtered_modules import FilteredModulesCategory

class FilteredAlgebras(FilteredModulesCategory):
    r"""
    The category of filtered algebras.

    An algebra `A` over `R` is *filtered* if `A` is a filtered `R`-module
    such that `F_i \cdot F_j \subseteq F_{i+j}` for all `i, j` in the
    filteration group.

    EXAMPLES::

        sage: Algebras(ZZ).Filtered()
        Category of filtered algebras over Integer Ring
        sage: Algebras(ZZ).Filtered().super_categories()
        [Category of algebras over Integer Ring,
         Category of filtered modules over Integer Ring]

    TESTS::

        sage: TestSuite(Algebras(ZZ).Filtered()).run()
    """
    class ParentMethods:
        @abstract_method(optional=True)
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Algebra of An example of a filtered module with basis:
                 the universal enveloping algebra of
                 Lie algebra of RR^3 with cross product over Integer Ring
            """

    class ElementMethods:
        pass

