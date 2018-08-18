r"""
Graded Lie Algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.graded_modules import GradedModulesCategory

class GradedLieAlgebras(GradedModulesCategory):
    r"""
    Category of graded Lie algebras.

    TESTS::

        sage: C = LieAlgebras(QQ).Graded()
        sage: TestSuite(C).run()
    """
    class SubcategoryMethods:
        def Stratified(self):
            r"""
            Return the full subcategory of stratified objects of ``self``.

            A Lie algebra is stratified if it is graded and generated as a
            Lie algebra by its component of degree one.

            EXAMPLES::

                sage: LieAlgebras(QQ).Graded().Stratified()
                Category of stratified graded Lie algebras over Rational Field
            """
            return self._with_axiom("Stratified")

    class Stratified(CategoryWithAxiom_over_base_ring):
        r"""
        Category of stratified Lie algebras.

        TESTS::

            sage: C = LieAlgebras(QQ).Graded().Stratified()
            sage: TestSuite(C).run()
        """
        pass

