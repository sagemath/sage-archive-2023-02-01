from abvar             import ModularAbelianVariety
from sage.rings.all    import QQ

class ModAbVar_ambient_jacobian(ModularAbelianVariety):
    def __init__(self, group):
        self._group = group
        ModularAbelianVariety.__init__(self, level = group.level(), base_ring = QQ)

    def _repr_(self):
        return "Jacobian of the modular curve associated to %s"%self._group

    def group(self):
        return self._group




