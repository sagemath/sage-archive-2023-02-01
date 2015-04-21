r"""
H-Trivial Semigroups
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

class HTrivialSemigroups(CategoryWithAxiom):
    def Finite_extra_super_categories(self):
        r"""
        Implement the fact that a finite `H`-trivial is aperiodic

        EXAMPLES::

            sage: Semigroups().HTrivial().Finite_extra_super_categories()
            [Category of aperiodic semigroups]
            sage: Semigroups().HTrivial().Finite() is Semigroups().Aperiodic().Finite()
            True
        """
        return [Semigroups().Aperiodic()]
