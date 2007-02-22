r"""
Base class for all monoids
"""

import sage.categories.all
from sage.structure.parent_gens import ParentWithGens

def is_Monoid(x):
    return isinstance(x, Monoid_class)

class Monoid_class(ParentWithGens):
    def category(self):
        return sage.categories.all.Monoids()
