from sage.structure.sage_object import SageObject

class TorsionPoint(SageObject):
    def __init__(self, abvar, x):
        self._abvar = abvar
        self._x = x

    def _repr_(self):
        return "Torsion point defined by %s on %s"%\
               (self._abvar, self._x)

