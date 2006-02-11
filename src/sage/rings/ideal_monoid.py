"""
Monoid of Ring Ideals
"""

from sage.monoids.monoid import Monoid_class
from commutative_ring import is_CommutativeRing
def IdealMonoid(R):
    return IdealMonoid_c(R)

class IdealMonoid_c(Monoid_class):
    def __init__(self, R):
        if not is_CommutativeRing(R):
            raise TypeError, "R (=%s) must be a commutative ring"%R
        self.__R = R

    def _repr_(self):
        return "Monoid of ideals of %s"%self.__R

    def ring(self):
        return self.__R

    def __call__(self, x):
        return self.__R.ideal(x)
