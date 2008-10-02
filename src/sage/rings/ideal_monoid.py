"""
Monoid of Ring Ideals
"""

from sage.monoids.monoid import Monoid_class
from commutative_ring import is_CommutativeRing
from sage.structure.parent_gens import ParentWithGens
import sage.rings.integer_ring
import ideal

def IdealMonoid(R):
    return IdealMonoid_c(R)

class IdealMonoid_c(Monoid_class):
    def __init__(self, R):
        if not is_CommutativeRing(R):
            raise TypeError, "R must be a commutative ring"
        self.__R = R
        ParentWithGens.__init__(self, sage.rings.integer_ring.ZZ)

    def _repr_(self):
        return "Monoid of ideals of %s"%self.__R

    def ring(self):
        return self.__R

    def __call__(self, x):
        if isinstance(x, ideal.Ideal_generic):
            x = x.gens()
        return self.__R.ideal(x)

    def _coerce_impl(self, x):
        R = self.__R
        if isinstance(x, ideal.Ideal_generic):
            x = [R.coerce(y) for y in x.gens()]
            return R.ideal(x)
        else:
            return R.ideal(R._coerce_(x))
