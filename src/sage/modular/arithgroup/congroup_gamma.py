r"""
Congruence Subgroup `\Gamma(N)`
"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from congroup_generic import CongruenceSubgroup
from arithgroup_element import ArithmeticSubgroupElement

_gamma_cache = {}
def Gamma_constructor(N):
    r"""
    Return the congruence subgroup `\Gamma(N)`.

    EXAMPLES::

        sage: Gamma(5) # indirect doctest
        Congruence Subgroup Gamma(5)
        sage: G = Gamma(23)
        sage: G is Gamma(23)
        True
        sage: G == loads(dumps(G))
        True
        sage: G is loads(dumps(G))
        True
    """
    from all import SL2Z
    if N == 1: return SL2Z
    try:
        return _gamma_cache[N]
    except KeyError:
        _gamma_cache[N] = Gamma_class(N)
        return _gamma_cache[N]

class Gamma_class(CongruenceSubgroup):
    r"""
    The principal congruence subgroup `\Gamma(N)`.
    """


    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: Gamma(133)._repr_()
            'Congruence Subgroup Gamma(133)'
        """
        return "Congruence Subgroup Gamma(%s)"%self.level()

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES::

            sage: Gamma(5).__reduce__()
            (<function Gamma_constructor at ...>, (5,))
        """
        return Gamma_constructor, (self.level(),)

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: Gamma(3) == SymmetricGroup(8)
            False
            sage: Gamma(3) == Gamma1(3)
            False
            sage: Gamma(5) < Gamma(6)
            True
            sage: Gamma(5) == Gamma(5)
            True
        """
        if not is_Gamma(other):
            return cmp(type(self), type(other))
        else:
            return cmp(self.level(), other.level())

    def __call__(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES::

            sage: G = Gamma(5)
            sage: G([1, 0, -10, 1])
            [ 1   0]
            [-10  1]
            sage: G(matrix(ZZ, 2, [26, 5, 5, 1]))
            [26  5]
            [ 5  1]
            sage: G([1, 1, 6, 7])
            Traceback (most recent call last):
            ...
            TypeError: matrix must have diagonal entries (=1, 7) congruent to 1
            modulo 5, and off-diagonal entries (=1,6) divisible by 5
        """
        if isinstance(x, ArithmeticSubgroupElement) and x.parent() == self:
            return x
        x = ArithmeticSubgroupElement(self, x, check=check)
        if not check:
            return x

        a = x.a()
        b = x.b()
        c = x.c()
        d = x.d()
        N = self.level()
        if (a%N == 1) and (c%N == 0) and (d%N == 1) and (b%N == 0):
            return x
        else:
            raise TypeError, "matrix must have diagonal entries (=%s, %s)\
            congruent to 1 modulo %s, and off-diagonal entries (=%s,%s)\
            divisible by %s" %(a, d, N, b, c, N)


def is_Gamma(x):
    r"""
    Return True if x is a congruence subgroup of type Gamma.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_Gamma
        sage: is_Gamma(Gamma0(13))
        False
        sage: is_Gamma(Gamma(4))
        True
    """

    return isinstance(x, Gamma_class)

