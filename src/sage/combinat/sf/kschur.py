"""
k-Schur Functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
from sage.rings.all import Integer, gcd, lcm, QQ, is_PolynomialRing, is_FractionField
import sage.combinat.partition
from sage.misc.misc import prod
import sfa
import copy


def kSchurFunctions(R, k, t=None):
    """
    Returns the k-Schur functions. See the examples below for caveats
    on their use.

    EXAMPLES::

        sage: ks3 = kSchurFunctions(QQ, 3); ks3
        k-Schur Functions at level 3 over Univariate Polynomial Ring in t over Rational Field
        sage: s = SFASchur(ks3.base_ring()) # Allow 't' coefficients for the Schurs.
        sage: t = ks3.t # Allow 't' as input
        sage: s(ks3([3,2,1]))
        s[3, 2, 1] + t*s[4, 1, 1] + t*s[4, 2] + t^2*s[5, 1]

    ::

        sage: ks3(s([3, 2, 1]) + t*s([4, 1, 1]) + t*s([4, 2]) + t^2*s([5, 1]))
        ks3[3, 2, 1]

    ::

        sage: ks3([4,3,2,1])  # k-Schurs are indexed by partitions with first part \le k.
        0

    Attempting to convert a function that is not in the linear span of
    the k-Schur's raises an error.

    ::

        sage: ks3(s([4]))
        Traceback (most recent call last):
        ...
        ValueError: s[4] is not in the space spanned by k-Schur Functions at level 3 over Univariate Polynomial Ring in t over Rational Field.

    Note that the product of k-Schurs is not guaranteed to be in the
    space spanned by the k-Schurs. In general, we only have that a
    k-Schur times a j-Schur is a (k+j)-Schur. This fact is not
    currently incorporated into the Sage design, so the multiplication
    of k-Schur functions may return an error. This example shows how to
    get around this 'manually'.

    ::

        sage: ks2 = kSchurFunctions(QQ, 2)
        sage: ks2([2,1])^2
        Traceback (most recent call last):
        ...
        ValueError: s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + (2*t+2)*s[3, 2, 1] + (t^2+1)*s[3,
        3] + (2*t+1)*s[4, 1, 1] + (t^2+2*t+1)*s[4, 2] + (t^2+2*t)*s[5, 1] + t^2*s[6] is not in the
        space spanned by k-Schur Functions at level 2 over Univariate Polynomial Ring in t over
        Rational Field.

    ::

        sage: f = s(ks2([2,1]))^2; f # Convert to Schur functions first and multiply there.
        s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + (2*t+2)*s[3, 2, 1] + (t^2+1)*s[3,
        3] + (2*t+1)*s[4, 1, 1] + (t^2+2*t+1)*s[4, 2] + (t^2+2*t)*s[5, 1] + t^2*s[6]
        sage: ks4 = kSchurFunctions(QQ, 4)
        sage: ks4(f) # The product of two 'ks2's is a 'ks4'.
        ks4[2, 2, 1, 1] + ks4[2, 2, 2] + ks4[3, 1, 1, 1] + (t+2)*ks4[3, 2, 1] + (t^2+1)*ks4[3, 3] + (t+1)*ks4[4, 1, 1] + ks4[4, 2]

    However, at t=1, the product of k-Schurs is in the span of the
    k-Schurs. Below are some examples at t=1.

    ::

        sage: ks3 = kSchurFunctions(QQ, 3, 1); ks3
        k-Schur Functions at level 3 with t=1 over Rational Field
        sage: s = SFASchur(ks3.base_ring())
        sage: ks3(s([3]))
        ks3[3]
        sage: s(ks3([3,2,1]))
        s[3, 2, 1] + s[4, 1, 1] + s[4, 2] + s[5, 1]
        sage: ks3([2,1])^2
        ks3[2, 2, 1, 1] + ks3[2, 2, 2] + ks3[3, 1, 1, 1]
    """
    return cache_t(R, k, t)


class kSchurFunctions_generic(sfa.SymmetricFunctionAlgebra_generic):
    def _change_by_triangularity(self, el, to_other_cache, unitriang=False):
        """
        Returns self(el) converted by triangularity.

        INPUT:


        -  ``el`` - a symmetric function

        -  ``to_other_cache`` - a dictionary containing the
           change of basis from self to el's basis

        -  ``unitriang`` - a boolean, if True, the coefficient
           of part in self( el.parent()(part) ) is assumed to be 1.


        EXAMPLES::

            sage: ks3 = kSchurFunctions(QQ, 3)
            sage: ks3._s_cache(3)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: d = ks3._self_to_s_cache
            sage: l(d[3])
            [([1, 1, 1], [([1, 1, 1], 1)]), ([2, 1], [([2, 1], 1)]), ([3], [([3], 1)])]
            sage: el = ks3._s([2,1]) + ks3._s([1,1,1])
            sage: ks3._change_by_triangularity(el, d, True)
            ks3[1, 1, 1] + ks3[2, 1]
        """
        orig = el
        P = el.parent()
        zero = self.base_ring()(0)
        out = {}
        while el != 0:
            l = el.support()
            l.sort()
            part2 = l[0]
            n = part2.size()

            if not to_other_cache[n][part2]:
                raise ValueError,"%s is not in the space spanned by %s."%(orig,self)

            c = el.coefficient(part2)
            if not unitriang:
                c /= to_other_cache[n][part2][part2]
            el -= c*P._from_dict(to_other_cache[n][part2])

            out[part2] = out.get(part2,zero) + c

        return self._from_dict(out)


    def _multiply(self, left, right):
        """
        Multiply left and right by coverting to the Schurs, multiplying
        there, and converting back. Note that the product of k-Schurs with
        t is only guaranteed to be a sum of k-Schurs when t = 1.

        EXAMPLES::

            sage: ks3 = kSchurFunctions(QQ, 3)

        ::

            sage: ks3([1])^2 # indirect doctest
            ks3[1, 1] + ks3[2]

        ::

            sage: ks3([2,1])^2
            Traceback (most recent call last):
            ...
            ValueError: s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] +
            s[3, 3] + s[4, 1, 1] + s[4, 2] is not in the space spanned by k-Schur
            Functions at level 3 over Univariate Polynomial Ring in t over Rational
            Field.

        ::

            sage: ks3 = kSchurFunctions(QQ,3,1)
            sage: ks3([2,1])^2
            ks3[2, 2, 1, 1] + ks3[2, 2, 2] + ks3[3, 1, 1, 1]
        """
        return self( self._s(left) * self._s(right) )


class kSchurFunction_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    pass

s_to_k_cache = {}
k_to_s_cache = {}
class kSchurFunction_t(kSchurFunction_generic):
    pass
class kSchurFunctions_t(kSchurFunctions_generic):
    def __init__(self, R, k, t=None):
        """
        EXAMPLES::

            sage: kSchurFunctions(QQ, 3).base_ring()
            Univariate Polynomial Ring in t over Rational Field
            sage: kSchurFunctions(QQ, 3, t=1).base_ring()
            Rational Field
            sage: ks3 = kSchurFunctions(QQ, 3)
            sage: ks3 == loads(dumps(ks3))
            True
        """
        self.k = k
        self._name = "k-Schur Functions at level %s"%k
        self._prefix = "ks%s"%k
        self._element_class = kSchurFunction_t
        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])

        if t is None:
            R = R['t']
            self.t = R.gen()
        elif t not in R:
            raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
        else:
            self.t = R(t)
            if str(t) != 't':
                self._name += " with t=%s"%self.t


        CombinatorialAlgebra.__init__(self, R)

        self._s = sfa.SFASchur(self.base_ring())
        self._s_to_self_cache = s_to_k_cache.get(k, {})
        self._self_to_s_cache = k_to_s_cache.get(k, {})

    def _coerce_start(self, x):
        """
        Coerce things into the k-Schurs through the Schurs.

        EXAMPLES::

            sage: ks3 = kSchurFunctions(QQ, 3)
            sage: s = SFASchur(QQ)
            sage: ks3([4,3,2,1]) # indirect doctest
            0
            sage: ks3(s([2,1]))
            ks3[2, 1]
            sage: ks3(s([4]))
            Traceback (most recent call last):
            ...
            ValueError: s[4] is not in the space spanned by k-Schur Functions at level 3 over Univariate Polynomial Ring in t over Rational Field.
        """
        if x in sage.combinat.partition.Partitions():
            if len(x) > 0 and max(x) > self.k:
                return self(0)
            x = sage.combinat.partition.Partition(x)
            return self._from_dict({x:self.base_ring()(1)})

        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._s(x)
            for p in x.support():
                self._s_cache(p.size())
            return self._change_by_triangularity(x, self._self_to_s_cache, True)
        else:
            raise TypeError

    def _s_cache(self, n):
        """
        Computes the change of basis from the kSchurs to the Schurs for
        partitions of size n.

        EXAMPLES::

            sage: ks3 = kSchurFunctions(QQ, 3)
            sage: ks3._s_cache(3)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(ks3._self_to_s_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]), ([2, 1], [([2, 1], 1)]), ([3], [([3], 1)])]
        """
        if n in self._self_to_s_cache:
            return

        R = self.base_ring()
        t = self.t
        s = self._s
        zero = s(0)

        if n == 0:
            p = sage.combinat.partition.Partition_class([])
            self._self_to_s_cache[0] = {p: {p:R(1)}}
            return
        else:
            self._self_to_s_cache[n] = {}


        #Fill in the cache from the k-Schurs to the Schurs
        for p in sage.combinat.partition.Partitions_n(n):
            if max(p) > self.k:
                self._self_to_s_cache[n][p] = {}
                continue
            katom = p.k_atom(self.k)
            res = sum( [t**tab.charge()*s(tab.shape()) for tab in katom], zero)
            self._self_to_s_cache[n][p] = res.monomial_coefficients()



#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_t = Cache(kSchurFunctions_t)
