from sage.rings.all import ZZ
from morphism import Morphism

class HeckeOperator(Morphism):
    def __init__(self, abvar, n):
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n must be positive"
        self._abvar = abvar
        self._n = n

    def _repr_(self):
        return "Hecke operator T_{%s} on %s"%(self._n, self._abvar)

    def index(self):
        return self._n

    def characteristic_polynomial(self, var='x'):
        return self._abvar.rational_homology().hecke_polynomial(self._n, var)

    def charpoly(self, var='x'):
        return self.characteristic_polynomial(var)

