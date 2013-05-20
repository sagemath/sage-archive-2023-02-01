r"""
Strata of quadratic differentials on Riemann surfaces
"""
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.rational import Rational

class QuadraticStratum(SageObject):
    r"""
    Stratum of quadratic differentials.
    """
    def __init__(self, *l):
        """
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: loads(dumps(a)) == a
            True
        """
        if l == []:
            raise ValueError("The list must be non empty !")

        if isinstance(l[0], list) or isinstance(l[0], tuple):
            self._zeroes = []
            for (i, j) in l.iteritems():
                assert(isinstance(i, (int, Integer)) and
                       isinstance(j, (int, Integer)))
                self._zeroes += [i]*j

        else:
            for i in l:
                assert(isinstance(i, (int, Integer)))
            self._zeroes = sorted(list(l), reverse=True)


        self._genus = sum(l)/4 + 1
        assert(isinstance(self._genus, (int, Integer))  or
               (isinstance(self._genus,Rational) and self._genus.is_integral()))

        self._genus = Integer(self._genus)

    def __repr__(self):
        r"""
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: print a
            Q(-1, -1, -1, -1)
        """
        return "Q(" + str(self._zeroes)[1:-1] + ")"

    def __str__(self):
        r"""
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: print a
            Q(-1, -1, -1, -1)
        """
        return "Q(" + str(self._zeroes)[1:-1] + ")"

    def __eq__(self, other):
        r"""
        TESTS::

            sage: QuadraticStratum(0) == QuadraticStratum(0)
            True
            sage: QuadraticStratum(4) == QuadraticStratum(0)
            False
        """
        return type(self) == type(other) and self._zeroes == other._zeroes

    def __ne__(self, other):
        r"""
        TESTS::

            sage: QuadraticStratum(0) != QuadraticStratum(0)
            False
            sage: QuadraticStratum(4) != QuadraticStratum(0)
            True
        """
        return type(self) != type(other) or self._zeroes != other._zeroes

    def genus(self):
        r"""
        Returns the genus.

        EXAMPLES:

        ::

            sage: QuadraticStratum(-1,-1,-1,-1).genus()
            0
        """
        return self._genus
