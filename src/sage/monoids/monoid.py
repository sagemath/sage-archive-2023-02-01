r"""
Base class for all monoids
"""

import sage.categories.all
from sage.structure.all import Generators

def is_Monoid(x):
    return isinstance(x, Monoid_class)

class Monoid_class(Generators):
    def category(self):
        return sage.categories.all.Monoids()
