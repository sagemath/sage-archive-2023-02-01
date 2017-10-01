
from sage.structure.element import AdditiveGroupElement
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.structure.richcmp import richcmp, op_EQ, op_NE

class QmodnZ_Element(AdditiveGroupElement):
    r"""
    The ``QmodnZ_Element`` class represents an element of the abelian group Q/nZ.

    INPUT:

    - ``q`` -- a rational number.

    - ``parent`` -- the parent abelian group Q/nZ

    OUTPUT:

    The element `q` of abelian group Q/nZ, in standard form.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.qmodnz import QmodnZ
        sage: G = QQ/(19*ZZ)
        sage: q = G(400/19)
        sage: q
        39/19

    """

    def __init__(self, parent, x, construct=False):
        r"""
        Create an element of Q/nZ.

        EXAMPLES::

            sage: G = QQ/(3*ZZ)
            sage: g = G.random_element()
            sage: g
            3/8
        """

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
        r"""
        Returns the rational number representation of ``self``

        EXAMPLES::
            sage: G = QQ/(5*ZZ)
            sage: g = G(2/4); g
            1/2
            sage: q = lift(g); q
            1/2

        TESTS::
            sage: q.parent() == QQ
            True
        """

        return self._x

    def __neg__(self):
        r"""
        Returns the additive inverse of ``self`` in Q/nZ.

        EXAMPLES::
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: -g
            2/21

        TESTS::
            sage: G = QmodnZ(19/23)
            sage: g = G(15/23)
            sage: -g
            4/23
            sage: g + -g == G(0)
            True
        """

        if self._x == 0:
            return self
        else:
            QZ = self.parent()
            return QZ.element_class(QZ, QZ.n - self._x, True)

    def _add_(self, other):
        r"""
        Returns the sum of ``self`` and ``other`` in Q/nZ

        EXAMPLES::
            sage: G = QmodnZ(9/10)
            sage: g = G(5)
            sage: h = G(1/2)
            sage: g + h
            1/10
            sage: g + h == G(1/10) #indirect doctest
            True

        TESTS::
            sage: h + g == G(1/10)
            True
        """

        QZ = self.parent()
        ans = self._x + other._x
        if ans >= QZ.n:
            ans -= QZ.n
        return QZ.element_class(QZ, ans, True)

    def _sub_(self, other):
        r"""
        Returns the difference of ``self`` and ``other`` in Q/nZ

        EXAMPLES::
            sage: G = QmodnZ(9/10)
            sage: g = G(4)
            sage: h = G(1/2)
            sage: g - h
            4/5
            sage: h - g
            1/10
            sage: g - h == G(4/5) #indirect doctest
            True
            sage: h - g == G(1/10) #indirect doctest
            True
        """

        QZ = self.parent()
        ans = self._x - other._x
        if ans < 0:
            ans += QZ.n
        return QZ.element_class(QZ, ans, True)

    def _rmul_(self, c):
        r"""
        Returns the (right) scalar product of ``self`` by ``c`` in Q/nZ.

        EXAMPLES::
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: g*6
            1/7
        """
        QZ = self.parent()
        return QZ.element_class(QZ, self._x * c)

    def _lmul_(self, c):
        r"""
        Returns the (left) scalar product of ``self`` by ``c`` in Q/nZ.

        EXAMPLES::
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: 6*g
            1/7

        TESTS::
            sage: 6*g == g*6
            True
            sage: 6*g == 5*g
            False
        """
        return self._rmul_(c)

    def __div__(self, other):
        #TODO: This needs to be implemented.
        QZ = self.parent()
        other = ZZ(other)
        return QZ.element_class(QZ, self._x / other, True)

    def _repr_(self):
        r"""
        Display the element

        EXAMPLES::

            sage: G = QmodnZ(8)
            sage: g = G(25/7); g;
            25/7

        """

        return repr(self._x)

    def __hash__(self):
        r"""
        TESTS::

            sage: G = QmodnZ(4)
            sage: g = G(4/5)
            sage: hash(g)
            -7046029254386353128
            sage: hash(G(3/4))
            3938850096065010962
            sage: hash(G(1))
            1
        """

        return hash(self._x)

    def _richcmp_(self, right, op):
        r"""
        Compare ``self`` with ``right``.

        OUTPUT:

        boolean

        EXAMPLES::

            sage: G = QQ/(4*ZZ)
            sage: g = G(4/5)
            sage: h = G(6/7)
            sage: g == h
            False
            sage: g == g
            True
        """
        if op == op_EQ or op == op_NE:
            return richcmp(self._x, right._x, op)
        else:
            return NotImplemented

    def additive_order(self):
        r"""
        Returns the order of ``self`` in the abelian group Q/nZ.

        EXAMPLES::
            sage: G = QmodnZ(12)
            sage: g = G(5/3)
            sage: g.additive_order()
            36
            sage: (-g).additive_order() # indirect doctest
            36
        """

        # a/b * k = n/m * r
        QZ = self.parent()
        if QZ.n == 0:
            if self._x == 0:
                return ZZ(1)
            return infinity
        return (self._x / QZ.n).denominator()
