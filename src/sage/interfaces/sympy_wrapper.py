"""
Wrapper Class for Sage Sets as SymPy Sets
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sympy.core.basic import Basic
from sympy.core.decorators import sympify_method_args
from sympy.core.sympify import sympify
from sympy.sets.sets import Set

@sympify_method_args
class SageSet(Set):

    def __new__(cls, sage_set):
        return Basic.__new__(cls, sage_set)

    def _sage_(self):
        return self._args[0]

    @property
    def is_empty(self):
        return self._sage_().is_empty()

    @property
    def is_finite_set(self):
        return self._sage_().is_finite()

    @property
    def is_iterable(self):
        from sage.categories.enumerated_sets import EnumeratedSets
        return self._sage_() in EnumeratedSets()

    def __iter__(self):
        for element in self._sage_():
            yield sympify(element)

    def _contains(self, other):
        return other in self._sage_()

    def __len__(self):
        return len(self._sage_())
