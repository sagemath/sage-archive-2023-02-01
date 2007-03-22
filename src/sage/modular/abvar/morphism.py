import sage.categories.all

class Morphism(sage.categories.all.Morphism):
    def __init__(self, parent, x):
        self._parent = parent
        self._x = x

    def _repr_(self):
        return "Morphism defined by %s in %s"%(self._x, self._parent)
