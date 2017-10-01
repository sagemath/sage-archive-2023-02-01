from sage.groups.group import AbelianGroup
from sage.rings.all import ZZ, QQ
from sage.groups.abelian_gps.qmodnz_element import QmodnZ_Element
from sage.categories.groups import Groups
from sage.arith import srange

class QmodnZ(AbelianGroup):
    r"""
    The ``QmodnZ`` class represents the abelian group Q/nZ.

    INPUT:

    The constructor may be called in any of the following ways.

    #. ``QmodnZ(n)``, where

        - `n` -- a rational number (including 0 or negative rational numbers).

    #. ``QQ/(n*ZZ)``, where

        - `n` -- an integer (including 0 or negative integers).


    OUTPUT:

    The abelian group Q/nZ.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.qmodnz import QmodnZ
        sage: QQ/(19*ZZ)
        Q/19Z

        sage: QmodnZ(19)
        Q/19Z

        sage: QmodnZ(2/3)
        Q/(2/3)Z
    """

    Element = QmodnZ_Element
    def __init__(self, n=1):
        r"""
        EXAMPLES::

            sage: G = QmodnZ(2)
            sage: G
            Q/2Z

        TESTS::

            sage: G = QQ/(19*ZZ)
            sage: TestSuite(G).run()
        """
        self.n = QQ(n).abs()
        category = Groups().Commutative().Topological().Infinite()
        AbelianGroup.__init__(self, base=ZZ, category=category)
        self._populate_coercion_lists_(coerce_list=[QQ])

    def _repr_(self):
        r"""
        Display the group.

        EXAMPLES::

            sage: G = QmodnZ(1); G
            Q/Z

            sage: G = QQ/(3*ZZ); G
            Q/3Z

            sage: G = QmodnZ(1/5); G
            Q/(1/5)Z
        """
        if self.n == 1:
            return "Q/Z"
        elif self.n in ZZ:
            return "Q/%sZ"%(self.n)
        else:
            return "Q/(%s)Z"%(self.n)

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are identical.

        This means they are the same abelian group.

        EXAMPLES::

            sage: G = QQ/(5*ZZ)
            sage: H = QQ/(5*ZZ)
            sage: K = QQ/(3*ZZ)
            sage: L = QmodnZ(25/5)
            sage: G == H #indirect doctest
            True
            sage: G == K #indirect doctest
            False

        TESTS::

            sage: G == G
            True
            sage: H == G
            True
            sage: G == L
            True
            sage: K == 7
            False
        """
        if type(self) == type(other):
            return self.n == other.n
        return NotImplemented

    def __ne__(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are not identical.

        This means they are different abelian groups.

        EXAMPLES::

            sage: G = QQ/(5*ZZ)
            sage: H = QQ/(5*ZZ)
            sage: K = QQ/(3*ZZ)
            sage: L = QmodnZ(25/5)
            sage: G != H #indirect doctest
            False
            sage: G != K #indirect doctest
            True
        """

        if type(self) == type(other):
            return self.n != other.n
        return NotImplemented

    #TODO: Disallow order comparisons between different Q/nZ's
    # e.g., sage: QmodnZ(10/3) > QmodnZ(5/3)
    # returns False.

    def _element_constructor_(self, x):
        r"""
        Construct an element in Q/nZ.

        EXAMPLES::

            sage: G = QmodnZ(2/3)
            sage: G(5/6)
            1/6
        """
        return self.element_class(self, QQ(x))

    def random_element(self, *args, **kwds):
        r"""
        Return a random element of Q/nZ.

        EXAMPLES::

            sage: G = QQ/(6*ZZ)
            sage: G.random_element()
            1
            sage: G.random_element()
            35/6
            sage: G.random_element()
            1/4

        Extra positional or keyword arguments are passed through::

            sage: G = QmodnZ(4/5)
            sage: G.random_element(distribution='1/n')
            1/2
            sage: G.random_element(distribution='1/n')
            3/5
            sage: G.random_element(distribution='1/n')
            11/20
        """
        return self(QQ.random_element(*args, **kwds))

    def __iter__(self):
        r"""
        Creates an iterator that generates the elements of Q/nZ without
        repetition, organized by increasing denominator; for a fixed denominator
        elements are listed by increasing numerator.

        EXAMPLES:

            The first 19 elements of Q/5Z::

            sage: import itertools
            sage: lst = [a for a in itertools.islice(QQ/(5*ZZ),19)]; lst
            [0, 1, 2, 3, 4, 1/2, 3/2, 5/2, 7/2, 9/2, 1/3, 2/3, 4/3, 5/3, 7/3, 8/3, 10/3, 11/3, 13/3]
        """

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
