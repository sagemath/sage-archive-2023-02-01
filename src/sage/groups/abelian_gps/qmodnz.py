from sage.groups.group import AbelianGroup
from sage.rings.all import ZZ, QQ
from sage.groups.abelian_gps.qmodnz_element import QmodnZ_Element
from sage.categories.groups import Groups
from sage.arith import srange

class QmodnZ(AbelianGroup):
    Element = QmodnZ_Element
    def __init__(self, n=1):
        self.n = QQ(n).abs()
        category = Groups().Commutative().Topological().Infinite()
        AbelianGroup.__init__(self, base=ZZ, category=category)
        self._populate_coercion_lists_(coerce_list=[QQ])

    def _repr_(self):
        if self.n == 1:
            return "Q/Z"
        elif self.n in ZZ:
            return "Q/%sZ"%(self.n)
        else:
            return "Q/(%s)Z"%(self.n)

    def __eq__(self, other):
        if type(self) == type(other):
            return self.n == other.n
        return NotImplemented

    def __ne__(self, other):
        if type(self) == type(other):
            return self.n != other.n
        return NotImplemented

    def _element_constructor_(self, x):
        return self.element_class(self, QQ(x))

    def random_element(self, *args, **kwds):
        return self(QQ.random_element(*args, **kwds))

    def __iter__(self):
        if self.n == 0:
            for x in QQ:
                yield self(x)
        else:
            yield self(0)
            d = ZZ(1)
            while True:
                for a in d.coprime_integers((d*self.n).floor()):
                    yield self(a/d)
                d += 1
