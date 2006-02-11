r"""
Congruence subgroups of $\SL_2(\Z)$
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random

import sage.rings.arith as arith

from sage.groups.group import Group

from sage.rings.integer_ring import IntegerRing

from sage.rings.infinity import infinity

from congroup_element import CongruenceSubgroupElement

class CongruenceSubgroup(Group):
    def __init__(self, level):
        level=int(level)
        if level <= 0:
            raise ArithmeticError, "Congruence groups only defined for positive levels."
        self.__level = level

    def _repr_(self):
        return "Generic congruence subgroup"

    def __hash__(self):
        return hash(str(self))

    def are_equivalent(self, x, y):
        raise NotImplementedError

    def coset_reps(self):
        raise NotImplementedError

    def generators(self):
        raise NotImplementedError

    def level(self):
        return self.__level

    def __cmp__(self, right):
        raise NotImplementedError

    def is_abelian(self):
        return False

    def is_finite(self):
        return False

    def order(self):
        return infinity

    def __call__(self, x, check=True):
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        raise NotImplementedError

class Gamma0(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_0(N)$.

    EXAMPLES:
        sage: G = Gamma0(11); G
        Congruence Subgroup Gamma0(11)
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self, level):
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        return "Congruence Subgroup Gamma0(%s)"%self.level()

    def _latex_(self):
        return "\\Gamma_0(%s)"%self.level()

    def __cmp__(self, right):
        if not isinstance(right, Gamma0):
            return -1
        if self.level() < right.level():
            return -1
        elif self.level() > right.level():
            return 1
        return 0

class SL2Z(Gamma0):
    r"""
    The modular group $\SL_2(\Z)$.

    EXAMPLES:
        sage: G = SL2Z(); G
        Modular Group SL(2,Z)
        sage: loads(G.dumps()) == G
            True
    """
    def __init__(self):
        Gamma0.__init__(self, 1)

    def _repr_(self):
        return "Modular Group SL(2,Z)"

    def _latex_(self):
        return "\\mbox{\\rm SL}_2(%s)"%(IntegerRing()._latex_())

    def gens(self):
        try:
            return self.__gens
        except AttributeError:
            self.__gens = (self([0, -1, 1, 0]), self([1, 1, 0, 1]))
            return self.__gens

    def ngens(self):
        return 2

    def __call__(self, x, check=True):
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        return CongruenceSubgroupElement(self, x, check=check)


class Gamma1(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_1(N)$.

    EXAMPLES:
        sage: G = Gamma1(11); G
        Congruence Subgroup Gamma1(11)
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self, level):
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        return "Congruence Subgroup Gamma1(%s)"%self.level()

    def _latex_(self):
        return "\\Gamma_1(%s)"%self.level()

    def __cmp__(self, right):
        if not isinstance(right, Gamma1):
            return -1
        if self.level() < right.level():
            return -1
        elif self.level() > right.level():
            return 1
        return 0

class GammaH(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_H(N)$.

    TODO: This is NOT really implemented yet.
    """
    def __init__(self, level, H):
        CongruenceSubgroup.__init__(self, level)
        self.__H = H

    def _repr_(self):
        return "Congruence Subgroup Gamma_H(%s,H=%s)"%(self.level(), self.__H)

    def _latex_(self):
        return "\\Gamma_H(%s)"%self.level()



import sage.ext.congroup_pyx
degeneracy_coset_representatives_gamma0 = sage.ext.congroup_pyx.degeneracy_coset_representatives_gamma0
degeneracy_coset_representatives_gamma1 = sage.ext.congroup_pyx.degeneracy_coset_representatives_gamma1



## def xxx_degeneracy_coset_representatives_gamma0(N, M, t):
##     r"""
##     Let $N$ be a positive integer and $M$ a multiple of $N$.  Let $t$ be a
##     divisor of $N/M$, and let $T$ be the $2x2$ matrix $T=[0,0; 0,t]$.
##     This function returns representatives for the orbit set
##     $\Gamma_0(N) \backslash T \Gamma_0(M)$, where $\Gamma_0(N)$
##     acts on the left on $T \Gamma_0(M)$.

##     INPUT:
##         N -- int
##         M -- int (divisor of N)
##         t -- int (divisors of N/M)

##     OUTPUT:
##         list -- list of lists [a,b,c,d], where [a,b,c,d] should
##                 be viewed as a 2x2 matrix.

##     This function is used for computation of degeneracy maps between
##     spaces of modular symbols, hence its name.

##     We use that $T^{-1}*(a,b;c,d)*T = (a,bt,c/t,d),$
##     that the group $T^{-1}Gamma_0(N) T$ is contained in $\Gamma_0(M)$,
##     and that $\Gamma_0(N) T$ is contained in $T \Gamma_0(M)$.

##     ALGORITHM:
##     \begin{enumerate}
##     \item Compute representatives for $\Gamma_0(N/t,t)$ inside of $\Gamma_0(M)$:
##           COSET EQUIVALENCE:
##            Two right cosets represented by $[a,b;c,d]$ and
##            $[a',b';c',d']$ of $\Gamma_0(N/t,t)$ in $\SL_2(\Z)$
##            are equivalent if and only if
##            $(a,b)=(a',b')$ as points of $\P^1(\Z/t\Z)$,
##            i.e., $ab' \con ba' \pmod{t}$,
##            and $(c,d) = (c',d')$ as points of $\P^1(\Z/(N/t)\Z)$.

##         ALGORITHM to list all cosets:
##         \begin{enumerate}
##            \item Compute the number of cosets.
##            \item Compute a random element x of Gamma_0(M).
##            \item Check if x is equivalent to anything generated so
##                far; if not, add x to the list.
##            \item Continue until the list is as long as the bound
##                computed in A.
##         \end{enumerate}

##     \item There is a bijection between $\Gamma_0(N)\backslash T \Gamma_0(M)$
##           and $\Gamma_0(N/t,t) \Gamma_0(M)$ given by $T r$ corresponds to $r$.
##           Consequently we obtain coset representatives for
##           $\Gamma_0(N)\backslash T \Gamma_0(M)$ by left multiplying by $T$ each
##           coset representative of $\Gamma_0(N/t,t) \Gamma_0(M)$ found
##           in step 1.
##     \end{enumerate}
##     """
##     import sage.modular.dims as dims
##     N = int(N); M = int(M); t = int(t)
##     if N % M != 0:
##         raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)
##     n     = dims.idxG0(N) // dims.idxG0(M)          # number of coset representatives.
##     Ndivd = N // t
##     R     = []                        # coset reps found so far.
##     halfmax = 2*(n+10)
##     while len(R) < n:
##         # try to find another coset representative.
##         cc = M*random.randrange(-halfmax, halfmax+1)
##         dd =   random.randrange(-halfmax, halfmax+1)
##         g, bb, aa = arith.xgcd(-cc,dd)
##         if g == 0: continue
##         cc //= g
##         dd //= g
##         if cc % M != 0: continue
##         # Test if we've found a new coset representative.
##         is_new = True
##         for r in R:
##             if ((r[1]*aa - r[0]*bb) % t == 0) and \
##                      ((r[3]*cc - r[2]*dd) % Ndivd == 0):
##                 is_new = False
##                 break
##         # If our matrix is new add it to the list.
##         if is_new: R.append([aa,bb,cc,dd])
##     # Return the list left multiplied by T.
##     return [[r[0], r[1], r[2]*t, r[3]*t] for r in R]
