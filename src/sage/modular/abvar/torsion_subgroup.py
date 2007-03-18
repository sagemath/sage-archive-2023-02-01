from sage.structure.sage_object import SageObject

class TorsionSubgroup(SageObject):
    def __init__(self, abvar):
        self._abvar = abvar

    def _repr_(self):
        return "Torsion subgroup of %s"%self._abvar

    def multiple_of_order(self, bound=10):
        raise NotImplementedError
