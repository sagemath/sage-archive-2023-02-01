"""
LLT Polynomials
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
import sfa
import sage.combinat.ribbon_tableau as ribbon_tableau
import sage.combinat.skew_partition
from sage.rings.all import QQ, ZZ
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra
import sage.combinat.partition


def LLT(R, k, t=None):
    """
    Returns a class for working with LLT polynomials.

    EXAMPLES::

        sage: L3 = LLT(QQ,3); L3
        LLT polynomials at level 3 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: L3.cospin([3,2,1])
        (t+1)*m[1, 1] + m[2]
        sage: L3.hcospin()
        LLT polynomials in the HCosp basis at level 3 with t=t over Fraction Field of Univariate Polynomial Ring in t over Rational Field
    """
    return cache_llt(R, k, t)

class LLT_class:
    def __init__(self, R, k, t=None):
        """
        EXAMPLES::

            sage: L3 = LLT(QQ,3)
            sage: L3 == loads(dumps(L3))
            True
        """
        self._k =k
        self._name = "LLT polynomials at level %s"%self._k

        if t is None:
            R = R['t'].fraction_field()
            self._t = R.gen()
        elif t not in R:
            raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
        else:
            self._t = R(t)
            self._name += " with t=%s"%self._t
        self._name += " over %s"%R

        self._m = sfa.SFAMonomial(R)

    def __repr__(self):
        """
        EXAMPLES::

            sage: LLT(QQ,3)
            LLT polynomials at level 3 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: LLT(QQ,3,2)
            LLT polynomials at level 3 with t=2 over Rational Field
        """
        return self._name

    def base_ring(self):
        """
        Returns the base ring of self.

        EXAMPLES::

            sage: LLT(QQ,3).base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._m.base_ring()

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: L3Q = LLT(QQ,3)
            sage: L4Q = LLT(QQ,4)
            sage: L3Z = LLT(ZZ, 4)
            sage: cmp(L3Q, L3Q)
            0
            sage: cmp(L3Q, L3Z)
            -1
            sage: cmp(L3Q, L4Q)
            -1
            sage: cmp(L3Q, QQ)
            1
            sage: cmp(LLT(QQ,3,2), LLT(QQ,3,3))
            -1
        """
        if other.__class__ is not self.__class__:
            return cmp(other.__class__, self.__class__)
        if self._k !=  other._k:
            return cmp(self._k, other._k)
        if self._m != other._m:
            return cmp(self._m, other._m)
        if self._t != other._t:
            return cmp(self._t, other._t)
        return 0


    def level(self):
        """
        Returns the level of self.

        EXAMPLES::

            sage: LLT(QQ,3).level()
            3
        """
        return self._k

    def _llt_generic(self, skp, stat):
        """
        Takes in partition, list of partitions, or a list of skew
        partitions as well as a statistic which takes in two partitions and
        a level and spits out a coefficient.

        EXAMPLES::

            sage: L3 = LLT(QQ,3)
            sage: f = lambda skp,mu,level: QQ(1)
            sage: L3._llt_generic([3,2,1],f)
            m[1, 1] + m[2]
        """
        if skp in sage.combinat.partition.Partitions():
            m = (sum(skp) / self.level()).floor()
            if m == 0:
                raise ValueError, "level (%=) must divide %s "%(sum(skp), self.level())
            mu = sage.combinat.partition.Partitions( ZZ(sum(skp) / self.level()) )


        elif isinstance(skp, list) and skp[0] in sage.combinat.skew_partition.SkewPartitions():
            #skp is a list of skew partitions
            skp =  [sage.combinat.partition.Partition(core_and_quotient=([], skp[i][0])) for i in range(len(skp))]
            skp += [sage.combinat.partition.Partition(core_and_quotient=([], skp[i][1])) for i in range(len(skp))]
            mu = sage.combinat.partition.Partitions(ZZ(sum( [ s.size() for s in skp] ) / self.level()))


        elif isinstance(skp, list) and skp[0] in sage.combinat.partition.Partitions():
            #skp is a list of partitions
            skp = sage.combinat.partition.Partition(core_and_quotient=([], skp))
            mu = sage.combinat.partition.Partitions( ZZ(sum(skp) / self.level() ))
        else:
            raise ValueError, "LLT polynomials not defined for %s"%skp


        BR = self.base_ring()
        return sum([ stat(skp,nu,self.level()).subs(t=self._t)*self._m(nu) for nu in mu])

    def spin_square(self, skp):
         r"""
         Returns the spin polynomial associated with skp with the
         substitution `t \rightarrow t^2` made.

         EXAMPLES::

             sage: L3 = LLT(QQ,3)
             sage: L3.spin_square([2,1])
             t*m[1]
             sage: L3.spin_square([3,2,1])
             (t^3+t)*m[1, 1] + t^3*m[2]
             sage: L3.spin_square([[1],[1],[1]])
             (t^6+2*t^4+2*t^2+1)*m[1, 1, 1] + (t^6+t^4+t^2)*m[2, 1] + t^6*m[3]
         """
         return self._llt_generic(skp, ribbon_tableau.spin_polynomial_square)

    def cospin(self, skp):
        """
        EXAMPLES::

            sage: L3 = LLT(QQ,3)
            sage: L3.cospin([2,1])
            m[1]
            sage: L3.cospin([3,2,1])
            (t+1)*m[1, 1] + m[2]
            sage: s = SFASchur(L3.base_ring())
            sage: s(L3.cospin([[2],[1],[2]]))
            t^4*s[2, 2, 1] + t^3*s[3, 1, 1] + (t^3+t^2)*s[3, 2] + (t^2+t)*s[4, 1] + s[5]
        """
        return self._llt_generic(skp, ribbon_tableau.cospin_polynomial)

##     def llt_inv(self, skp):
##         """
##         """
##         l = sage.combinat.partitions( sum( [ p.size() for p in skp ] ) ).list()
##         res = m(0)
##         for p in l:
##             inv_p = [ ktuple.inversions() for ktuple in kTupleTableaux(skp, p) ]
##             res += sum([t**x for x in inv_p])*m(p)
##         return res

    def hcospin(self):
        """
        Returns the HCopsin basis.

        EXAMPLES::

            sage: LLT(QQ,3).hcospin()
            LLT polynomials in the HCosp basis at level 3 with t=t over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return cache_llt_cospin(self._m.base_ring(), self._k, self._t)

    def hspin(self):
        """
        Returns the HSpin basis.

        EXAMPLES::

            sage: LLT(QQ,3).hspin()
            LLT polynomials in the HSp basis at level 3 with t=t over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return cache_llt_spin(self._m.base_ring(), self._k, self._t)



class LLT_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, level, t=None):
        """
        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: LLT_spin(QQ, 3)
            LLT polynomials in the HSp basis at level 3 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: LLT_spin(QQ, 3, t=2)
            LLT polynomials in the HSp basis at level 3 with t=2 over Rational Field
        """
        if t is None:
            R = R['t'].fraction_field()
            self.t = R.gen()
        elif t not in R:
            raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
        else:
            self.t = R(t)
            self._name += " with t=%s"%self.t

        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])
        self._level = level

        CombinatorialAlgebra.__init__(self, R)

    def level(self):
        """
        Returns the level of self.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: HSp3.level()
            3
        """
        return self._level

class LLTElement_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    pass



def LLTHSpin(R, level, t=None):
    """
    Returns the LLT polynomials in the HSpin basis at level level.

    EXAMPLES::

        sage: HSp3 = LLTHSpin(QQ,3)
        sage: HSp3([1])^2
        HSp[1, 1] + (-t+1)*HSp[2]
    """
    return cache_llt_spin(R, level, t)

hsp_to_m_cache={}
m_to_hsp_cache={}
class LLT_spin(LLT_generic):
    def __init__(self, R, level, t=None):
        """
        TESTS::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: HSp3 == loads(dumps(HSp3))
            True
        """
        self._name = "LLT polynomials in the HSp basis at level %s"%level
        self._prefix = "HSp"
        self._element_class = LLTElement_spin

        LLT_generic.__init__(self, R, level, t=t)

        self._m = sfa.SFAMonomial(self.base_ring())

        if level not in hsp_to_m_cache:
            hsp_to_m_cache[level] = {}
            m_to_hsp_cache[level] = {}
        self._self_to_m_cache = hsp_to_m_cache[level]
        self._m_to_self_cache = m_to_hsp_cache[level]


    def _multiply(self, left, right):
        """
        Convert to the monomial basis, do the multiplication there, and
        convert back to the HSp basis.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: HSp3([1])^2 #indirect doctest
            HSp[1, 1] + (-t+1)*HSp[2]
        """
        return self( self._m(left) * self._m(right) )

    def _to_m(self, part):
        """
        Returns a function which gives the coefficient of a partition
        in the monomial expansion of self(part).

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: f21 = HSp3._to_m(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [t, t + 1, t + 2]
        """
        BR = self.base_ring()
        level = self.level()
        f = lambda part2: BR(ribbon_tableau.spin_polynomial([level*i for i in part], part2, level))
        return f

    def _m_cache(self, n):
        """
        Compute the change of basis from the monomial symmetric functions
        to self.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: HSp3._m_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( HSp3._self_to_m_cache[2] )
            [([1, 1], [([1, 1], t + 1), ([2], t)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l( HSp3._m_to_self_cache[2] )
            [([1, 1], [([1, 1], 1), ([2], -t)]), ([2], [([1, 1], -1), ([2], t + 1)])]
        """
        self._invert_morphism(n, self.base_ring(), self._self_to_m_cache, \
                              self._m_to_self_cache, to_other_function = self._to_m)

    def _coerce_start(self, x):
        """
        Coerce things into the LLT HSp basis through the monomials.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HSp3 = LLT_spin(QQ, 3)
            sage: s = SFASchur(HSp3.base_ring())
            sage: HSp3._coerce_start(s([2]))
            HSp[2]
            sage: HSp3._coerce_start(s([1,1]))
            HSp[1, 1] - t*HSp[2]
        """
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._m(x)
            return self._from_cache(x, self._m_cache, self._m_to_self_cache, t=self.t)
        else:
            raise TypeError


class LLTElement_spin(LLTElement_generic):
    pass



def LLTHCospin(R, level, t=None):
    """
    Returns the LLT polynomials in the HCospin basis at level level.

    EXAMPLES::

        sage: HCosp3 = LLTHCospin(QQ,3)
        sage: HCosp3([1])^2
        1/t*HCosp3[1, 1] + ((t-1)/t)*HCosp3[2]
    """
    return cache_llt_cospin(R, level, t)

hcosp_to_m_cache={}
m_to_hcosp_cache={}
class LLT_cospin(LLT_generic):
    def __init__(self, R, level, t=None):
        """
        TESTS::

            sage: from sage.combinat.sf.llt import *
            sage: HCosp3 = LLT_cospin(QQ, 3)
            sage: HCosp3 == loads(dumps(HCosp3))
            True
        """
        self._name = "LLT polynomials in the HCosp basis at level %s"%level
        self._prefix = "HCosp"+str(level)
        self._element_class = LLTElement_cospin

        LLT_generic.__init__(self, R, level, t=t)

        self._m = sfa.SFAMonomial(self.base_ring())

        if level not in hcosp_to_m_cache:
            hcosp_to_m_cache[level] = {}
            m_to_hcosp_cache[level] = {}
        self._self_to_m_cache = hcosp_to_m_cache[level]
        self._m_to_self_cache = m_to_hcosp_cache[level]


    def _multiply(self, left, right):
        """
        Convert to the monomial basis, do the multiplication there, and
        convert back to the Hcosp basis.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HCosp3 = LLT_cospin(QQ, 3)
            sage: HCosp3([1])^2 #indirect doctest
            1/t*HCosp3[1, 1] + ((t-1)/t)*HCosp3[2]
        """
        return self( self._m(left) * self._m(right) )

    def _to_m(self, part):
        """
        Returns a function which gives the coefficient of part2 in the
        monomial expansion of self(part).

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HCosp3 = LLT_cospin(QQ, 3)
            sage: f21 = HCosp3._to_m(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [1, t + 1, 2*t + 1]
        """
        BR = self.base_ring()
        level = self.level()
        f = lambda part2: BR(ribbon_tableau.cospin_polynomial([level*i for i in part], part2, level))
        return f

    def _m_cache(self, n):
        """
        Compute the change of basis from the monomial symmetric functions
        to self.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HCosp3 = LLT_cospin(QQ, 3)
            sage: HCosp3._m_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( HCosp3._self_to_m_cache[2] )
            [([1, 1], [([1, 1], t + 1), ([2], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l( HCosp3._m_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/t), ([2], -1/t)]),
             ([2], [([1, 1], -1/t), ([2], (t + 1)/t)])]
        """
        self._invert_morphism(n, self.base_ring(), self._self_to_m_cache, \
                              self._m_to_self_cache, to_other_function = self._to_m)

    def _coerce_start(self, x):
        """
        Coerce things into the LLT HCosp basis through the monomials.

        EXAMPLES::

            sage: from sage.combinat.sf.llt import *
            sage: HCosp3 = LLT_cospin(QQ, 3)
            sage: s = SFASchur(HCosp3.base_ring())
            sage: HCosp3._coerce_start(s([2]))
            HCosp3[2]
            sage: HCosp3._coerce_start(s([1,1]))
            1/t*HCosp3[1, 1] - 1/t*HCosp3[2]
        """
        BR = self.base_ring()
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._m(x)
            return self._from_cache(x, self._m_cache, self._m_to_self_cache, t=self.t)
        else:
            raise TypeError

class LLTElement_cospin(LLTElement_generic):
    pass


from sage.misc.cache import Cache
cache_llt = Cache(LLT_class)
cache_llt_cospin = Cache(LLT_cospin)
cache_llt_spin = Cache(LLT_spin)
