r"""
J-Trivial Semigroups
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2009-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.semigroups import Semigroups

class JTrivialSemigroups(CategoryWithAxiom):
    def extra_super_categories(self):
        """
        Implement the fact that a `J`-trivial semigroup is `L` and `R`-trivial.

        EXAMPLES::

            sage: Semigroups().JTrivial().extra_super_categories()
            [Category of ltrivial semigroups, Category of rtrivial semigroups]
        """
        return [Semigroups().LTrivial(), Semigroups().RTrivial()]
