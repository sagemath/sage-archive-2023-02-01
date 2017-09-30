
from sage.structure.element import AdditiveGroupElement
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.structure.richcmp import richcmp, op_EQ, op_NE

class QmodnZ_Element(AdditiveGroupElement):
    def __init__(self, parent, x, construct=False):
        AdditiveGroupElement.__init__(self, parent)
        # x = (a/b) = q(n/m) + r/mb
        # am = q(nb) + r
        # r < nb so r/mb < n/m
        if construct:
            self._x = x
            return
        n = parent.n.numerator()
        if n == 0:
            self._x = x
        else:
            m = parent.n.denominator()
            a = x.numerator()
            b = x.denominator()
            q, r = (a*m).quo_rem(n*b)
            self._x = r/(m*b)

    def lift(self):
        return self._x

    def __neg__(self):
        if self._x == 0:
            return self
        else:
            QZ = self.parent()
            return QZ.element_class(QZ, QZ.n - self._x, True)

    def _add_(self, other):
        QZ = self.parent()
        ans = self._x + other._x
        if ans > QZ.n:
            ans -= QZ.n
        return QZ.element_class(QZ, ans, True)

    def _sub_(self, other):
        QZ = self.parent()
        ans = self._x - other._x
        if ans < 0:
            ans += QZ.n
        return QZ.element_class(QZ, ans, True)

    def _rmul_(self, c):
        QZ = self.parent()
        return QZ.element_class(QZ, self._x * c)

    def _lmul_(self, c):
        return self._rmul_(c)

    def __truediv__(self, other):
        QZ = self.parent()
        other = ZZ(Other)
        return QZ.element_class(QZ, self._x / other, True)

    def _repr_(self):
        return repr(self._x)

    def __hash__(self):
        return hash(self._x)

    def _richcmp_(self, right, op):
        if op == op_EQ or op == op_NE:
            return richcmp(self._x, right._x, op)
        else:
            return NotImplemented

    def additive_order(self):
        # a/b * k = n/m * r
        QZ = self.parent()
        if QZ.n == 0:
            if self._x == 0:
                return ZZ(1)
            return infinity
        return (self._x / QZ.n).denominator()
