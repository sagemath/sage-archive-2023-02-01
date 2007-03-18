from sage.categories.homset import HomsetWithBase

class Homspace(HomsetWithBase):
    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    def _repr_(self):
        return "Space of homomorphisms from %s to %s"%\
               (self._domain, self._codomain)
