r"""
The modular group `{\rm SL}_2(\ZZ)`
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

from congroup_gamma0 import Gamma0_class
from sage.rings.integer_ring import ZZ
from sage.modular.cusps import Cusp
from sage.rings.arith import gcd
from sage.modular.modsym.p1list import lift_to_sl2z

def is_SL2Z(x):
    r"""
    Return True if x is the modular group `{\rm SL}_2(\ZZ)`.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_SL2Z
        sage: is_SL2Z(SL2Z)
        True
        sage: is_SL2Z(Gamma0(6))
        False
    """
    return isinstance(x, SL2Z_class)

class SL2Z_class(Gamma0_class):
    r"""
    The full modular group `{\rm SL}_2(\ZZ)`, regarded as a congruence
    subgroup of itself.
    """

    def __init__(self):
        r"""
        The modular group ${\rm SL}_2(\Z)$.

        EXAMPLES::

            sage: G = SL2Z; G
            Modular Group SL(2,Z)
            sage: G.gens()
            ([ 0 -1]
            [ 1  0], [1 1]
            [0 1])
            sage: G.0
            [ 0 -1]
            [ 1  0]
            sage: G.1
            [1 1]
            [0 1]
            sage: latex(G)
            \mbox{\rm SL}_2(\Bold{Z})
            sage: G([1,-1,0,1])
            [ 1 -1]
            [ 0  1]
            sage: loads(G.dumps()) == G
            True
            sage: SL2Z.0 * SL2Z.1
            [ 0 -1]
            [ 1  1]

            sage: SL2Z == loads(dumps(SL2Z))
            True
            sage: SL2Z is loads(dumps(SL2Z))
            True
        """
        Gamma0_class.__init__(self, 1)

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES::

            sage: SL2Z.__reduce__()
            (<function _SL2Z_ref at ...>, ())
        """
        return _SL2Z_ref, ()

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: SL2Z._repr_()
            'Modular Group SL(2,Z)'
        """
        return "Modular Group SL(2,Z)"

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES::

            sage: SL2Z._latex_()
            '\\mbox{\\rm SL}_2(\\Bold{Z})'
            sage: latex(SL2Z)
            \mbox{\rm SL}_2(\Bold{Z})
        """
        return "\\mbox{\\rm SL}_2(%s)"%(ZZ._latex_())

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES::

            sage: SL2Z.is_subgroup(SL2Z)
            True
            sage: SL2Z.is_subgroup(Gamma1(1))
            True
            sage: SL2Z.is_subgroup(Gamma0(6))
            False
        """
        return right.level() == 1

    def reduce_cusp(self, c):
        r"""
        Return the unique reduced cusp equivalent to c under the
        action of self. Always returns Infinity, since there is only
        one equivalence class of cusps for $SL_2(Z)$.

        EXAMPLES::

            sage: SL2Z.reduce_cusp(Cusps(-1/4))
            Infinity
        """
        return Cusp(1,0)

    def random_element(self, bound=100):
        r"""
        Return a random element of `{\rm SL}_2(\ZZ)` with entries whose
        absolute value is strictly less than bound (default 100).

        (Algorithm: Generate a random pair of integers at most bound. If they
        are not coprime, throw them away and start again. If they are, find an
        element of `{\rm SL}_2(\ZZ)` whose bottom row is that, and
        left-multiply it by `\begin{pmatrix} 1 & w \\ 0 & 1\end{pmatrix}` for
        an integer `w` randomly chosen from a small enough range that the
        answer still has entries at most bound.)

        It is, unfortunately, not true that all elements of SL2Z with entries <
        bound appear with equal probability; those with larger bottom rows are
        favoured, because there are fewer valid possibilities for w.

        EXAMPLES::

            sage: SL2Z.random_element() # random
            sage: SL2Z.random_element(5) # still random
        """

        if bound <= 1: raise ValueError, "Don't be silly"
        c = ZZ.random_element(1-bound, bound)
        d = ZZ.random_element(1-bound, bound)
        if gcd(c,d) != 1: # try again
            return self.random_element(bound)
        else:
            a,b,c,d = lift_to_sl2z(c,d,0)
            whi = bound
            wlo = bound
            if c > 0:
                whi = min(whi, ((bound - a)/ZZ(c)).ceil())
                wlo = min(wlo, ((bound + a)/ZZ(c)).ceil())
            elif c < 0:
                whi = min(whi, ((bound + a)/ZZ(-c)).ceil())
                wlo = min(wlo, ((bound - a)/ZZ(-c)).ceil())

            if d > 0:
                whi = min(whi, ((bound - b)/ZZ(d)).ceil())
                wlo = min(wlo, ((bound + b)/ZZ(d)).ceil())
            elif d < 0:
                whi = min(whi, ((bound + b)/ZZ(-d)).ceil())
                wlo = min(wlo, ((bound - b)/ZZ(-d)).ceil())

            w = ZZ.random_element(1-wlo, whi)
            a += c*w
            b += d*w
            return self([a,b,c,d])


SL2Z = SL2Z_class()

def _SL2Z_ref():
    """
    Return SL2Z. (Used for pickling SL2Z.)

    EXAMPLES::

        sage: sage.modular.arithgroup.congroup_sl2z._SL2Z_ref()
        Modular Group SL(2,Z)
        sage: sage.modular.arithgroup.congroup_sl2z._SL2Z_ref() is SL2Z
        True
    """
    return SL2Z


