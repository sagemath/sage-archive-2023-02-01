from sage.modular.hecke.all import HeckeModule_free_module


class Homology(HeckeModule_free_module):
    def __init__(self, abvar):
        self._abvar = abvar

    def _repr_(self):
        return "Homology of %s"%self._abvar

    def hecke_matrix(self, n):
        raise NotImplementedError

class IntegralHomology(Homology):
    def _repr_(self):
        return "Integral Homology of %s"%self._abvar

    def hecke_matrix(self, n):
        return self._abvar._integral_hecke_matrix(n)

class RationalHomology(Homology):
    def _repr_(self):
        return "Rational Homology of %s"%self._abvar
    def hecke_matrix(self, n):
        return self._abvar._rational_hecke_matrix(n)


class Homology_over_base(Homology):
    def __init__(self, abvar, base_ring):
        Homology.__init__(self, abvar)
        self._base_ring = base_ring

    def _repr_(self):
        return "Homology with coefficients in %s of %s"%(self._base_ring, self._abvar)

    def hecke_matrix(self, n):
        return self._abvar._integral_hecke_matrix(n).change_ring(self._base_ring)
