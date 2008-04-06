"""
Jack Polynomials
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

def JackPolynomialsP(R, t=None):
    """
    Returns the algebra of Jack polynomials in the P basis.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: JackPolynomialsP(QQ)
        Jack polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsP(QQ,t=-1)
        Jack polynomials in the P basis with t=-1 over Rational Field

      At t = 1, the Jack polynomials on the P basis are the Schur symmetric functions.
        sage: P = JackPolynomialsP(QQ,1)
        sage: s = SFASchur(QQ)
        sage: P([2,1])^2
        JackP[2, 2, 1, 1] + JackP[2, 2, 2] + JackP[3, 1, 1, 1] + 2*JackP[3, 2, 1] + JackP[3, 3] + JackP[4, 1, 1] + JackP[4, 2]
        sage: s([2,1])^2
        s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

      At t = 2, the Jack polynomials on the P basis are the zonal polynomials.
        sage: P = JackPolynomialsP(QQ,2)
        sage: Z = ZonalPolynomials(QQ)
        sage: P([2])^2
        64/45*JackP[2, 2] + 16/21*JackP[3, 1] + JackP[4]
        sage: Z([2])^2
        64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]
        sage: Z(P([2,1]))
        Z[2, 1]
        sage: P(Z([2,1]))
        JackP[2, 1]
    """
    return cache_p(R,t)

def JackPolynomialsQ(R, t=None):
    """
    Returns the algebra of Jack polynomials in the Q basis.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: JackPolynomialsQ(QQ)
        Jack polynomials in the Q basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsQ(QQ,t=-1)
        Jack polynomials in the Q basis with t=-1 over Rational Field
    """
    return cache_q(R,t)

def JackPolynomialsJ(R, t=None):
    """
    Returns the algebra of Jack polynomials in the J basis.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: JackPolynomialsJ(QQ)
        Jack polynomials in the J basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsJ(QQ,t=-1)
        Jack polynomials in the J basis with t=-1 over Rational Field

      At t = 1, the Jack polynomials in the J basis are scalar multiples of the Schur
      functions with the scalar given by a Partition's hook_product method at 1.
        sage: J = JackPolynomialsJ(QQ, t=1)
        sage: s = SFASchur(J.base_ring())
        sage: p = Partition([3,2,1,1])
        sage: s(J(p)) == p.hook_product(1)*s(p)
        True

      At t = 2, the Jack polynomials on the J basis are scalar multiples of the
      zonal polynomials with the scalar given by a Partition's hook_product method
      at 1.
        sage: t = 2
        sage: J = JackPolynomialsJ(QQ,t=t)
        sage: Z = ZonalPolynomials(J.base_ring())
        sage: p = Partition([2,2,1])
        sage: Z(J(p)) == p.hook_product(t)*Z(p)
        True

    """
    return cache_j(R,t)

def JackPolynomialsQp(R, t=None):
    """
    Returns the algebra of Jack polynomials in the Qp, which is dual to the P
    basis with respect to the standard scalar product.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: Qp = JackPolynomialsQp(QQ)
        sage: P = JackPolynomialsP(QQ)
        sage: a = Qp([2])
        sage: a.scalar(P([2]))
        1
        sage: a.scalar(P([1,1]))
        0
        sage: P(Qp([2]))
        ((t-1)/(t+1))*JackP[1, 1] + JackP[2]

    """
    res = JackPolynomialsP(R, t=None).dual_basis(prefix="JackQp")
    res._name = "Jack polynomials in the Qp basis"
    if t is not None:
        res._name += "with t=%s"%t
    return res

def ZonalPolynomials(R):
    """
    Returns the algebra of zonal polynomials.

    EXAMPLES:
        sage: Z = ZonalPolynomials(QQ)
        sage: a = Z([2])
        sage: Z([2])^2
        64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]
    """
    return cache_z(R)

def _makeZonalPolynomials(R):
    """
    A utility function that is used to make the zonal polynomials
    from the Jack polynomials on the P basis.  This routine is
    called by the cache.

    EXAMPLES:
        sage: from sage.combinat.sf.jack import _makeZonalPolynomials
        sage: _makeZonalPolynomials(QQ)
        Zonal polynomials over Rational Field
    """
    res = JackPolynomialsP(R, t=R(2))
    res = copy.copy(res)
    res._name = "Zonal polynomials"
    res._prefix = "Z"
    return res


###################################################################
def c1(part, t):
    """
    EXAMPLES:
        sage: from sage.combinat.sf.jack import c1
        sage: t = QQ['t'].gen()
        sage: [c1(p,t) for p in Partitions(3)]
        [2*t^2 + 3*t + 1, t + 2, 6]
    """
    return prod([1+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))])

def c2(part, t):
    """
    EXAMPLES:
        sage: from sage.combinat.sf.jack import c2
        sage: t = QQ['t'].gen()
        sage: [c2(p,t) for p in Partitions(3)]
        [6*t^3, 2*t^3 + t^2, t^3 + 3*t^2 + 2*t]
    """
    return prod([t+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))])

####################################################################

class JackPolynomials_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: JackPolynomialsJ(QQ).base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: JackPolynomialsJ(QQ, t=2).base_ring()
            Rational Field
        """
        if t is None:
            R = R['t'].fraction_field()
            self.t = R.gen()
        elif t not in R:
            raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
        else:
            self.t = R(t)
            if str(t) != 't':
                self._name += " with t=%s"%self.t

        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])

        CombinatorialAlgebra.__init__(self, R)


    def _normalize_coefficients(self, c):
        """
        If our coefficient ring is the field of fractions over a univariate
        polynomial ring over the rationals, then we should clear both the
        numerator and denominator of the denominators of their coefficients.

        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: t = P.base_ring().gen()
            sage: a = 2/(1/2*t+1/2)
            sage: P._normalize_coefficients(a)
            4/(t + 1)
            sage: a = 1/(1/3+1/6*t)
            sage: P._normalize_coefficients(a)
            6/(t + 2)
            sage: a = 24/(4*t^2 + 12*t + 8)
            sage: P._normalize_coefficients(a)
            6/(t^2 + 3*t + 2)
        """
        BR = self.base_ring()
        if is_FractionField(BR) and BR.base_ring() == QQ:
            denom = c.denominator()
            numer = c.numerator()

            #Clear the denominators
            a = lcm([i.denominator() for i in denom.coeffs()])
            b = lcm([i.denominator() for i in numer.coeffs()])
            l = Integer(a).lcm(Integer(b))
            denom *= l
            numer *= l

            #Divide through by the gcd of the numerators
            a = gcd([i.numerator() for i in denom.coeffs()])
            b = gcd([i.numerator() for i in numer.coeffs()])
            l = Integer(a).gcd(Integer(b))

            denom = denom / l
            numer = numer / l

            return c.parent()(numer, denom)
        else:
            return c

class JackPolynomial_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    def scalar_jack(self, x):
        """
        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: Q = JackPolynomialsQ(QQ)
            sage: p = Partitions(3).list()
            sage: matrix([[P(a).scalar_jack(Q(b)) for a in p] for b in p])
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        parent = self.parent()

        #Make the power-sum basis
        p = sfa.SFAPower(parent.base_ring())

        p_self = p(self)
        p_x    = p(x)
        res = 0
        for m,c in p_self:
            if m in p_x:
                res += scalar_jack(m,m,self.parent().t)*(c*p_x.coefficient(m))

        return parent._normalize_coefficients(res)


def scalar_jack(part1, part2, t):
    """
    Returns the Jack scalar product between p(part1) and p(part2)
    where p is the power-sum basis.

    EXAMPLES:
        sage: Q.<t> = QQ[]
        sage: from sage.combinat.sf.jack import scalar_jack
        sage: matrix([[scalar_jack(p1,p2,t) for p1 in Partitions(4)] for p2 in Partitions(4)])
        [   4*t      0      0      0      0]
        [     0  3*t^2      0      0      0]
        [     0      0  8*t^2      0      0]
        [     0      0      0  4*t^3      0]
        [     0      0      0      0 24*t^4]
    """
    if part1 != part2:
        return 0
    else:
        return part1.centralizer_size()*t**len(part1)

#P basis

class JackPolynomial_p(JackPolynomial_generic):
    def scalar_jack(self, x):
        """
        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: l = [P(p) for p in Partitions(3)]
            sage: matrix([[a.scalar_jack(b) for a in l] for b in l])
            [  6*t^3/(2*t^2 + 3*t + 1)                         0                         0]
            [                        0     (2*t^3 + t^2)/(t + 2)                         0]
            [                        0                         0 1/6*t^3 + 1/2*t^2 + 1/3*t]
        """
        if isinstance(x, JackPolynomials_p):
            P = self.parent()
            f = lambda p1, p2: c2(p1, P.t)/c1(p1, P.t) if p1 == p2 else 0
            return P._apply_multi_module_morphism(self, x, f, orthogonal=True)
        else:
            return JackPolynomial_generic.scalar_jack(self, x)

class JackPolynomials_p(JackPolynomials_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: P == loads(dumps(P))
            True
        """
        self._name = "Jack polynomials in the P basis"
        self._prefix = "JackP"
        self._element_class = JackPolynomial_p
        JackPolynomials_generic.__init__(self, R, t=t)

        self._m = sfa.SFAMonomial(self.base_ring())
        self._m_to_self_cache = {}
        self._self_to_m_cache = {}

    def _coerce_start(self, x):
        """
        Coerce things into the Jack polynomials P basis.

        EXAMPLES:
            sage: Q = JackPolynomialsQ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: J = JackPolynomialsJ(QQ)

            sage: P(Q([2,1])) # indirect doctest
            ((t+2)/(2*t^3+t^2))*JackP[2, 1]
            sage: P(Q([3]))
            ((2*t^2+3*t+1)/(6*t^3))*JackP[3]
            sage: P(Q([1,1,1]))
            (6/(t^3+3*t^2+2*t))*JackP[1, 1, 1]

            sage: P(J([3]))
            (2*t^2+3*t+1)*JackP[3]
            sage: P(J([2,1]))
            (t+2)*JackP[2, 1]
            sage: P(J([1,1,1]))
            6*JackP[1, 1, 1]

            sage: s = SFASchur(QQ)
            sage: P(s([2,1]))
            ((2*t-2)/(t+2))*JackP[1, 1, 1] + JackP[2, 1]
            sage: s(_)
            s[2, 1]
        """
        BR = self.base_ring()
        if isinstance(x, JackPolynomial_q):
            f = lambda part: (1/self(part).scalar_jack(self(part)))*self(part)
            return self._apply_module_morphism(x, f).map_coefficients(self._normalize_coefficients)
        elif isinstance(x, JackPolynomial_j):
            f = lambda part: c1(part, self.t)*self(part)
            return self._apply_module_morphism(x,f).map_coefficients(self._normalize_coefficients)
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._m(x)
            return self._from_cache(x, self._m_cache, self._m_to_self_cache)
        else:
            raise TypeError

    def _m_cache(self, n):
        """
        Computes the change of basis between the Jack polynomials in the P basis
        and the monomial symmetric functions.  This uses Gram-Schmidt to go
        to the monomials, and then that matrix is simply inverted.

        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: P._m_cache(2)
            sage: l(P._self_to_m_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 2/(t + 1)), ([2], 1)])]
            sage: l(P._m_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -2/(t + 1)), ([2], 1)])]
            sage: P._m_cache(3)
            sage: l(P._m_to_self_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], -6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 6/(t^2 + 3*t + 2)), ([2, 1], -3/(2*t + 1)), ([3], 1)])]
            sage: l(P._self_to_m_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 6/(2*t^2 + 3*t + 1)), ([2, 1], 3/(2*t + 1)), ([3], 1)])]
        """
        if n in self._self_to_m_cache:
            return
        else:
            self._self_to_m_cache[n] = {}

        self._gram_schmidt(n, self._m, lambda p: scalar_jack(p,p,self.t), \
                           self._self_to_m_cache[n], upper_triangular=True)
        self._invert_morphism(n, self.base_ring(), self._self_to_m_cache, \
                              self._m_to_self_cache, to_other_function = self._to_m)

    def _to_m(self, part):
        """
        Return a function that takes in a partition lambda that returns the
        coefficient of lambda in the expansion of self(part) in the monomial
        basis.

        This assumes that the cache from the Jack polynomials in the P basis
        to the monomial symmetric functions has already been computed.

        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: P._m_cache(2)
            sage: p2 = Partition([2])
            sage: p11 = Partition([1,1])
            sage: f = P._to_m(p2)
            sage: f(p11)
            2/(t + 1)
            sage: f(p2)
            1
            sage: f = P._to_m(p11)
            sage: f(p2)
            0
            sage: f(p11)
            1
        """
        f = lambda part2: self._self_to_m_cache[sum(part)][part].get(part2, 0)
        return f

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: P([1])^2 # indirect doctest
            (2*t/(t+1))*JackP[1, 1] + JackP[2]
            sage: P._m(_)
            2*m[1, 1] + m[2]
            sage: P = JackPolynomialsP(QQ, 2)
            sage: P([2,1])^2
            125/63*JackP[2, 2, 1, 1] + 25/12*JackP[2, 2, 2] + 25/18*JackP[3, 1, 1, 1] + 12/5*JackP[3, 2, 1] + 4/3*JackP[3, 3] + 4/3*JackP[4, 1, 1] + JackP[4, 2]
            sage: P._m(_)
            45*m[1, 1, 1, 1, 1, 1] + 51/2*m[2, 1, 1, 1, 1] + 29/2*m[2, 2, 1, 1] + 33/4*m[2, 2, 2] + 9*m[3, 1, 1, 1] + 5*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]
        """
        return self( self._m(left)*self._m(right) )



#J basis


class JackPolynomial_j(JackPolynomial_generic):
    pass

class JackPolynomials_j(JackPolynomials_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: J = JackPolynomialsJ(QQ)
            sage: J == loads(dumps(J))
            True
        """
        self._name = "Jack polynomials in the J basis"
        self._prefix = "JackJ"
        self._element_class = JackPolynomial_j
        JackPolynomials_generic.__init__(self, R, t=t)

        self._P = JackPolynomialsP(self.base_ring(),t=self.t)

    def _coerce_start(self, x):
        """
        EXAMPLES:
            sage: J = JackPolynomialsJ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: J(sum(P(p) for p in Partitions(3)))
            1/6*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]

            sage: s = SFASchur(J.base_ring())
            sage: J(s([3])) # indirect doctest
            ((t^2-3*t+2)/(6*t^2+18*t+12))*JackJ[1, 1, 1] + ((2*t-2)/(2*t^2+5*t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]
            sage: J(s([2,1]))
            ((t-1)/(3*t+6))*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1]
            sage: J(s([1,1,1]))
            1/6*JackJ[1, 1, 1]
        """
        BR = self.base_ring()
        if isinstance(x, JackPolynomial_p):
            f = lambda m: BR(1/c1(m, self.t))
            return self._change_by_proportionality(x, f).map_coefficients(self._normalize_coefficients)
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._P( x ) )
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: J = JackPolynomialsJ(QQ)
            sage: J([1])^2 #indirect doctest
            (t/(t+1))*JackJ[1, 1] + (1/(t+1))*JackJ[2]
            sage: J([2])^2
            (2*t^2/(2*t^2+3*t+1))*JackJ[2, 2] + (4*t/(3*t^2+4*t+1))*JackJ[3, 1] + ((t+1)/(6*t^2+5*t+1))*JackJ[4]
        """
        return self( self._P(left) * self._P(right) )



#Q basis
class JackPolynomial_q(JackPolynomial_generic):
    pass

class JackPolynomials_q(JackPolynomials_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: Q = JackPolynomialsQ(QQ)
            sage: Q == loads(dumps(Q))
            True
        """
        self._name = "Jack polynomials in the Q basis"
        self._prefix = "JackQ"
        self._element_class = JackPolynomial_q
        JackPolynomials_generic.__init__(self, R, t=t)

        self._P = JackPolynomialsP(self.base_ring(), t=self.t)

    def _coerce_start(self, x):
        """
        EXAMPLES:
            sage: Q = JackPolynomialsQ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: Q(sum(P(p) for p in Partitions(3)))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]

            sage: s = SFASchur(Q.base_ring())
            sage: Q(s([3])) # indirect doctest
            (1/6*t^3-1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3-2*t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]
            sage: Q(s([2,1]))
            (1/3*t^3-1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1]
            sage: Q(s([1,1,1]))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1]
        """
        BR = self.base_ring()
        if isinstance(x, JackPolynomial_p):
            f = lambda m: BR(self._P(m).scalar_jack(self._P(m)))
            return self._change_by_proportionality(x, f).map_coefficients(self._normalize_coefficients)
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._P( x ) )
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: Q = JackPolynomialsQ(QQ)
            sage: Q([1])^2 # indirect doctest
            JackQ[1, 1] + (2/(t+1))*JackQ[2]
            sage: Q([2])^2
            JackQ[2, 2] + (2/(t+1))*JackQ[3, 1] + ((6*t+6)/(6*t^2+5*t+1))*JackQ[4]

        """
        return self( self._P(left) * self._P(right) )




#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_p = Cache(JackPolynomials_p)
cache_j = Cache(JackPolynomials_j)
cache_q = Cache(JackPolynomials_q)
cache_z = Cache(_makeZonalPolynomials)

