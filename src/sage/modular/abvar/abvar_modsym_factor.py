"""
Abelian Variety that is defined by a Hecke stable factor of
an ambient modular symbols space.
"""

from abvar             import ModularAbelianVariety
from sage.rings.all    import QQ

class ModAbVar_modsym_factor(ModularAbelianVariety):
    def __init__(self, modsym):
        self._modsym = modsym
        assert modsym.base_ring() == QQ
        ModularAbelianVariety.__init__(self, level = modsym.level(), base_ring = QQ)
