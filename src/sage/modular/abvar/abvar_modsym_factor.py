"""
Abelian Variety that is defined by a Hecke stable factor of
an ambient modular symbols space.
"""

from abvar             import ModularAbelianVariety_modsym
from sage.rings.all    import QQ

class ModAbVar_modsym_factor(ModularAbelianVariety_modsym):
    def __init__(self, ambient, modsym):
        self._ambient = ambient
        self._modsym = {1:modsym}
        assert modsym.base_ring() == QQ
        assert modsym.sign() == 1
        ModularAbelianVariety_modsym.__init__(self, level = modsym.level(), base_ring = QQ)

    def ambient_variety(self):
        """
        Return the ambient variety of which this modular abelian
        variety is a *quotient*.

        EXAMPLES:
            sage: A = J0(33)[1]; A
            Modular abelian subvariety of dimension 2 and level 33
            sage: A.ambient_variety()
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(33)
        """
        return self._ambient

    def modular_symbols(self, sign=0):
        """
        Return space of modular symbols (with given sign) associated
        to this modular abelian variety.

        INPUT:
            sign -- integer, either -1, 0 or 1 (default: 0)

        EXAMPLES:
            sage: A = J0(33)[1]; A
            Modular abelian subvariety of dimension 2 and level 33
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: A.modular_symbols(1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field
            sage: A.modular_symbols(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(33) of weight 2 with sign -1 over Rational Field
        """
        sign = int(sign)
        try:
            return self._modsym[sign]
        except KeyError:
            pass
        M = self._modsym[1]
        A = M.modular_symbols_of_sign(sign)
        self._modsym[sign] = A
        return A

    def _repr_(self):
        return "Modular abelian quotient variety of dimension %s and level %s"%(\
            self.dimension(), self.level())

##     def is_subvariety(self, other):
##         """
##         Return True if self is a subvariety of other, as they sit in an ambient
##         modular abelian variety.

##         EXAMPLES:

##         """
##         if not isinstance(other, ModularAbelianVariety_modsym):
##             return False
##         A = self.ambient_variety()
##         if A == other:
##             return True
##         B = other.ambient_variety()
##         if A != B:
##             return False
##         return self.modular_symbols(1).is_submodule(other.modular_symbols(1))


## class ModAbVar_modsym_simple_new_factor(ModAbVar_modsym_factor):

##     def _repr_(self):
##         return "New simple modular abelian quotient variety of dimension %s and level %s"%(\
##             self.dimension(), self.level())
