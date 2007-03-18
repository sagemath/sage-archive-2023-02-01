from sage.modules.module        import Module

class FiniteSubgroup(Module):
    def __init__(self, abvar):
        self._abvar = abvar

    def _generators(self):
        """
        Return a list of vectors that define elements of the rational
        homology that generate this finite subgroup.

        Raises a ValueError if no explicit presentation of this finite
        subgroup is known.
        """
        raise ValueError, "no explicit presentation of this finite subgroup is known"

    def abelian_variety(self):
        return self._abvar

    def _repr_(self):
        return "Finite subgroup of %s"%self._abvar

    def order(self):
        raise NotImplementedError, "unable to compute order of this group"

    def gens(self):
        try:
            return self.__gens
        except AttributeError:
            pass
        G = self._generators()


    def gen(self, n):
        return self.gens()[n]


