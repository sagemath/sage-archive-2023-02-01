# -*- coding: utf-8 -*-
r"""
H-trivial semigroups
"""
#*****************************************************************************
#  Copyright (C) 2016 Nicolas M. Thiéry <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

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
