r"""
R-Trivial Semigroups
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2009-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from semigroups import Semigroups

class RTrivialSemigroups(CategoryWithAxiom):
    def extra_super_categories(self):
        r"""
        Implement the fact that a `R`-trivial semigroup is `H`-trivial.

        EXAMPLES::

            sage: Semigroups().RTrivial().extra_super_categories()
            [Category of htrivial semigroups]
        """
        return [Semigroups().HTrivial()]

    def Commutative_extra_super_categories(self):
        r"""
        Implement the fact that a commutative `R`-trivial semigroup is `J`-trivial.

        EXAMPLES::

            sage: Semigroups().RTrivial().Commutative_extra_super_categories()
            [Category of jtrivial semigroups]

        TESTS::

            sage: Semigroups().RTrivial().Commutative() is Semigroups().JTrivial().Commutative()
            True
        """
        return [self.JTrivial()]
