from sage.modular.hecke.all import HeckeModule_free_module


class Homology(HeckeModule_free_module):
    def __init__(self, abvar):
        self._abvar = abvar

    def _repr_(self):
        return "Homology associated to %s"%self._abvar

class IntegralHomology(Homology):
    pass

class RationalHomology(Homology):
    pass


class Homology_over_base(Homology):
    def __init__(self, abvar, base_ring):
        Homology.__init__(self, abvar)
        self._base_ring = base_ring
