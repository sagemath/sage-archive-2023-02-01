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

    An algebra `A` over a commutative ring `R` is *filtered* if
    `A` is endowed with a structure of a filtered `R`-module
    (whose underlying `R`-module structure is identical with
    that of the `R`-algebra `A`) such that
    `F_i \cdot F_j \subseteq F_{i+j}` for all `i, j \in \NN`.

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

            .. TODO::

                Should this really be here and not in the
                ``_with_basis`` class? The notion of an associated
                graded algebra of a filtered algebra (without
                basis) exists, but are we ever going to get it into
                Sage, and, more importantly: will it cooperate with
                the with-basis version?

            EXAMPLES::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Algebra of An example of a filtered algebra with basis:
                 the universal enveloping algebra of
                 Lie algebra of RR^3 with cross product over Integer Ring
            """

    class ElementMethods:
        pass

