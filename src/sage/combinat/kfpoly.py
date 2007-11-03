#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2007 John Stembridge
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

from partition import Partitions_n, Partition_class, Partition, Partitions_all
from sage.rings.all import PolynomialRing, ZZ

global_t = PolynomialRing(ZZ, 't').gen()

def KostkaFoulkesPolynomial(mu, nu, t=None):
    r"""
    Returns the Kostka-Foulkes polynomial K_{mu nu}(t).

    EXAMPLES:
        sage: KostkaFoulkesPolynomial([2,2],[2,2])
        1
        sage: KostkaFoulkesPolynomial([2,2],[4])
        0
        sage: KostkaFoulkesPolynomial([2,2],[1,1,1,1])
        t^6 + t^5 + t^4 + t^2
        sage: KostkaFoulkesPolynomial([2,2],[2,1,1])
        t^2 + t
        sage: q = PolynomialRing(QQ,'q').gen()
        sage: KostkaFoulkesPolynomial([2,2],[2,1,1],q)
        q^2 + q
    """
    if mu not in Partitions_all():
        raise ValueError, "mu must be a partition"
    if nu not in Partitions_all():
        raise ValueError, "nu must be a partition"

    if sum(mu) != sum(nu):
        raise ValueError, "mu and nu must be partitions of the same size"

    return kfpoly(mu, nu, t)

def kfpoly(mu, nu, t=None):
    """
     kfpoly(mu,nu,t) computes the Kostka-Foulkes polynomial K[mu,nu](t)
     by generating all rigging sequences for the shape mu, and then
     selecting those of content nu.
    """
    if mu == nu:
        return 1
    elif mu == []:
        return 0

    if t is None:
        t = global_t

    nuc = Partition(nu).conjugate()

    def f(x):
        if x[0] == nuc:
            return weight(x, t)
        else:
            return 0

    res = sum([f(rg) for rg in riggings(mu)])
    return res

def schur_to_hl(mu, t=None):
    """
    Returns a dictionary corresponding to s(mu) in Hall-Littlewood
    P basis.

    EXAMPLES:
        sage: from sage.combinat.kfpoly import *
        sage: schur_to_hl([1,1,1])
        {[1, 1, 1]: 1}
        sage: a = schur_to_hl([2,1])
        sage: for m,c in sorted(a.iteritems()): print m, c
        [1, 1, 1] t^2 + t
        [2, 1] 1
        sage: a = schur_to_hl([3])
        sage: for m,c in sorted(a.iteritems()): print m, c
        [1, 1, 1] t^3 + 1
        [2, 1] t
        [3] 1

    """
    if mu == []:
        return mu
    if t is None:
        t = global_t

    res = {}
    for rg in riggings(mu):
        res[ rg[0].conjugate() ] = res.get(rg[0], 0) + weight(rg, t)

    return res

def riggings(part):
    """
    Generate all possible rigging sequences for a
    fixed partition.

    EXAMPLES:
        sage: from sage.combinat.kfpoly import *
        sage: riggings([3])
        [[[1, 1, 1]], [[2, 1]], [[3]]]
        sage: riggings([2,1])
        [[[2, 1], [1]], [[3], [1]]]
        sage: riggings([1,1,1])
        [[[3], [2], [1]]]
        sage: riggings([2,2])
        [[[2, 2], [1, 1]], [[3, 1], [1, 1]], [[4], [1, 1]], [[4], [2]]]
    """
    l = len(part)
    res = [ [[],[]] ]
    sa = 0
    for i in sorted(part):
        sa += i
        new_res = []
        for nu in res:
            new_res += [ [new]+nu for new in compat(sa, nu[0], nu[1])]
        res = new_res


    return map(lambda x: x[:l], res)

def compat(n, mu, nu):
    """
    Generate all possible partitions of n that can
    precede mu,nu in a rigging sequence.

    EXAMPLES:
        sage: from sage.combinat.kfpoly import *
        sage: compat(4, [1], [2,1])
        [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
        sage: compat(3, [1], [2,1])
        [[1, 1, 1], [2, 1], [3]]
        sage: compat(2, [1], [])
        [[2]]
        sage: compat(3, [1], [])
        [[2, 1], [3]]
        sage: compat(3, [2], [1])
        [[3]]
        sage: compat(4, [1,1], [])
        [[2, 2], [3, 1], [4]]
        sage: compat(4, [2], [])
        [[4]]

    """
    sp = map(lambda p: p.conjugate(),Partitions_n(n))
    l = max( len(mu), len(nu))
    mmu = list(mu) + [0]*(l-len(mu))
    nnu = list(nu) + [0]*(l-len(nu))

    bd = []
    sa = 0
    for i in range(l):
        sa += 2*mmu[i]-nnu[i]
        bd.append(sa)

    for i in range(len(sp)):
        if dom(sp[i],bd):
            break

    if i >= len(sp):
        return Partition_class([])
    else:
        return sp[i:]

def dom(mu, snu):
    l = len(snu)
    mmu = list(mu)+[0]*l
    sa = 0
    for i in range(l):
        sa += mmu[i]
        if sa < snu[i]:
            return False
    return True

def weight(rg, t=None):
    """
    Returns the wieght of a rigging.

    EXAMPLES:
        sage: from sage.combinat.kfpoly import *
        sage: t = PolynomialRing(ZZ, 't').gen()
        sage: weight([[2,1], [1]], t)
        1
        sage: weight([[3],[1]], t)
        t^2 + t
        sage: weight([[2,1],[3]], t)
        t^4
    """
    if t is None:
        t = global_t
    nu = rg + [ [] ]
    l = 1 + max( map(len, nu) )
    nu = [ list(mu) + [0]*l for mu in nu ]

    res = t**int(sum( [ i*(i-1)/2 for i in rg[-1] ] ))

    for k in range(1, len(nu)-1):
        sa = 0
        for i in range( max( len(rg[k]), len(rg[k-1])) ):
            sa += nu[k-1][i]-2*nu[k][i]+nu[k+1][i]
            res *= q_bin(nu[k][i]-nu[k+1][i+1], sa, t)
            mu = nu[k-1][i]-nu[k][i]
            res *= t**int((mu*(mu-1)/2))

    return res


def q_bin(a,b,t=None):
    """
    EXAMPLES:
        sage: from sage.combinat.kfpoly import *
        sage: t = PolynomialRing(ZZ, 't').gen()
        sage: q_bin(4,2, t)
        t^8 + t^7 + 2*t^6 + 2*t^5 + 3*t^4 + 2*t^3 + 2*t^2 + t + 1
        sage: q_bin(4,3, t)
        t^12 + t^11 + 2*t^10 + 3*t^9 + 4*t^8 + 4*t^7 + 5*t^6 + 4*t^5 + 4*t^4 + 3*t^3 + 2*t^2 + t + 1
    """
    if t is None:
        t = global_t
    res = 1
    for i in range(1, min(a,b)+1):
        res *= (1-t**(a+b+1-i))/(1-t**i)
    return res
