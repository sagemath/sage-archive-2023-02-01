"""
Macdonald Polynomials -- under development.
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

from sage.combinat.combinatorial_algebra import CombinatorialAlgebra
import sfa
import sage.combinat.partition
from sage.matrix.all import MatrixSpace
from sage.rings.all import QQ
from sage.misc.misc import prod
import functools
QQqt = QQ['q,t'].fraction_field()

def MacdonaldPolynomialsP(R, q=None, t=None):
    """
    Returns the Macdonald polynomials on the P basis.  These are
    upper triangularly related to the monomial symmetric functions
    and are orthogonal with respect to the qt-Hall scalar product.

    EXAMPLES:
        sage: P = MacdonaldPolynomialsP(QQ); P
        Macdonald polynomials in the P basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: m = SFAMonomial(P.base_ring())
        sage: P.transition_matrix(m,2)
        [                          1 (q*t - q + t - 1)/(q*t - 1)]
        [                          0                           1]
        sage: P([1,1]).scalar_qt(P([2]))
        0
        sage: P([2]).scalar_qt(P([2]))
        (q^3 - q^2 - q + 1)/(q*t^2 - q*t - t + 1)
        sage: P([1,1]).scalar_qt(P([1,1]))
        (q^2*t - q*t - q + 1)/(t^3 - t^2 - t + 1)


      When q = 0, the Macdonald polynomials on the P basis are the same
      as the Hall-Littlewood polynomials on the P basis.
        sage: P = MacdonaldPolynomialsP(QQ,q=0)
        sage: P([2])^2
        (t+1)*McdP[2, 2] + (-t+1)*McdP[3, 1] + McdP[4]
        sage: HLP = HallLittlewoodP(QQ)
        sage: HLP([2])^2
        (t+1)*P[2, 2] + (-t+1)*P[3, 1] + P[4]

    """
    return cache_p(R, q, t)

def MacdonaldPolynomialsQ(R, q=None, t=None):
    """
    Returns the Macdonald polynomials on the Q basis.  These are
    dual to the Macdonald polynomials on the P basis with
    respect to the qt-Hall scalar product.

    EXAMPLES:
        sage: Q = MacdonaldPolynomialsQ(QQ); Q
        Macdonald polynomials in the Q basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: P = MacdonaldPolynomialsP(QQ)
        sage: Q([2]).scalar_qt(P([2]))
        1
        sage: Q([2]).scalar_qt(P([1,1]))
        0
        sage: Q([1,1]).scalar_qt(P([2]))
        0
        sage: Q([1,1]).scalar_qt(P([1,1]))
        1
        sage: Q(P([2]))
        ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
        sage: Q(P([1,1]))
        ((q^2*t-q*t-q+1)/(t^3-t^2-t+1))*McdQ[1, 1]

    """
    return cache_q(R, q, t)

def MacdonaldPolynomialsJ(R, q=None, t=None):
    """
    Returns the Macdonald polynomials on the J basis also known as
    the integral form of the Macdonald polynomials.  These are
    scalar multiples of both the P and Q bases. When expressed in
    the P or Q basis, the scaling coefficients are polynomials in
    q and t rather than rational functions.

    EXAMPLES:
        sage: J = MacdonaldPolynomialsJ(QQ); J
        Macdonald polynomials in the J basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: P = MacdonaldPolynomialsP(QQ)
        sage: Q = MacdonaldPolynomialsQ(QQ)
        sage: P(J([2]))
        (q*t^2-q*t-t+1)*McdP[2]
        sage: P(J([1,1]))
        (t^3-t^2-t+1)*McdP[1, 1]
        sage: Q(J([2]))
        (q^3-q^2-q+1)*McdQ[2]
        sage: Q(J([1,1]))
        (q^2*t-q*t-q+1)*McdQ[1, 1]

    """
    return cache_j(R, q, t)

def MacdonaldPolynomialsH(R, q=None, t=None):
    """
    Returns the Macdonald polynomials on the H basis. When the H
    basis is expanded on the Schur basis, the coefficients are
    the qt-Kostka numbers.

    EXAMPLES:
        sage: H = MacdonaldPolynomialsH(QQ)
        sage: s = SFASchur(H.base_ring())
        sage: s(H([2]))
        q*s[1, 1] + s[2]
        sage: s(H([1,1]))
        s[1, 1] + t*s[2]
    """
    return cache_h(R, q, t)


def MacdonaldPolynomialsHt(R, q=None, t=None):
    """
    Returns the Macdonald polynomials on the Ht basis. The elements
    of the Ht basis are eigenvectors of the nabla operator. When
    expanded on the Schur basis, the coefficients are the modified
    qt-Kostka numbers.

    EXAMPLES:
        sage: Ht = MacdonaldPolynomialsHt(QQ)
        sage: [Ht(p).nabla() for p in Partitions(3)]
        [q^3*McdHt[3], q*t*McdHt[2, 1], t^3*McdHt[1, 1, 1]]

        sage: s = SFASchur(Ht.base_ring())
        sage: from sage.combinat.sf.macdonald import qt_kostka
        sage: q,t = Ht.base_ring().gens()
        sage: s(Ht([2,1]))
        q*t*s[1, 1, 1] + (q+t)*s[2, 1] + s[3]
        sage: qt_kostka([1,1,1],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
        q*t
        sage: qt_kostka([2,1],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
        q + t
        sage: qt_kostka([3],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
        1

    """
    return cache_ht(R, q, t)

def MacdonaldPolynomialsS(R, q=None, t=None):
    """
    Returns the modified Schur functions defined by the
    plethystic subsitution  $S_{\mu} = s_{\mu}[X(1-t)]$.
    When the Macdonald polynomials in the J basis are expressed
    in terms of the modified Schur functions, the coefficients
    are qt-Kostka numbers.

    EXAMPLES:
        sage: S = MacdonaldPolynomialsS(QQ)
        sage: J = MacdonaldPolynomialsJ(QQ)
        sage: S(J([2]))
        q*McdS[1, 1] + ((-1)/(-1))*McdS[2]
        sage: S(J([1,1]))
        McdS[1, 1] + t*McdS[2]
        sage: from sage.combinat.sf.macdonald import qt_kostka
        sage: qt_kostka([2],[1,1])
        t
        sage: qt_kostka([1,1],[2])
        q

    """
    return cache_s(R, q, t)









##############################################

def c1(part, q, t):
    """
    This function returns the qt-Hall scalar product between
    J(part) and P(part).

    EXAMPLES:
        sage: from sage.combinat.sf.macdonald import c1
        sage: R.<q,t> = QQ[]
        sage: c1(Partition([2,1]),q,t)
        -q^4*t + 2*q^3*t - q^2*t + q^2 - 2*q + 1
        sage: c1(Partition([1,1]),q,t)
        q^2*t - q*t - q + 1
    """
    res = 1
    for i in range(part.size()):
        res *= 1-q**(sum(part.arm_lengths(),[])[i]+1)*t**(sum(part.leg_lengths(),[])[i])
    return res

def c2(part, q, t):
    """
    This function returns the qt-Hall scalar product between
    J(part) and Q(part).

    EXAMPLES:
        sage: from sage.combinat.sf.macdonald import c2
        sage: R.<q,t> = QQ[]
        sage: c2(Partition([1,1]),q,t)
        t^3 - t^2 - t + 1
        sage: c2(Partition([2,1]),q,t)
        -q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1
     """
    res = 1
    for i in range(part.size()):
        res *= 1-q**(sum(part.arm_lengths(),[])[i])*t**(sum(part.leg_lengths(),[])[i]+1)
    return res


#Generic MacdonaldPolynomials
class MacdonaldPolynomials_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, q=None, t=None):
        """
        EXAMPLES:
            sage: MacdonaldPolynomialsP(QQ)
            Macdonald polynomials in the P basis over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
            sage: MacdonaldPolynomialsP(QQ,t=2)
            Macdonald polynomials in the P basis with t=2 over Fraction Field of Univariate Polynomial Ring in q over Rational Field
            sage: MacdonaldPolynomialsP(QQ,q=2)
            Macdonald polynomials in the P basis with q=2 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: MacdonaldPolynomialsP(QQ,q=2,t=2)
            Macdonald polynomials in the P basis with q=2 and t=2 over Rational Field
        """
        #Neither q nor t are specified
        if t is None and q is None:
            R = R['q,t'].fraction_field()
            self.q, self.t = R.gens()
        elif t is not None and q is None:
            if t not in R:
                raise ValueError, "t (=%s) must be in R (=%s)"%(t,R)
            self.t = R(t)
            self._name += " with t=%s"%self.t

            R = R['q'].fraction_field()
            self.q = R.gen()

        elif t is None and q is not None:
            if q not in R:
                raise ValueError, "q (=%s) must be in R (=%s)"%(q,R)
            self.q = R(q)
            self._name += " with q=%s"%self.q

            R = R['t'].fraction_field()
            self.t = R.gen()
        else:
            if t not in R or q not in R:
                raise ValueError
            self.q = R(q)
            self.t = R(t)
            self._name += " with q=%s and t=%s"%(self.q, self.t)

        self._combinatorial_class = sage.combinat.partition.Partitions()
        self._one = sage.combinat.partition.Partition([])
        CombinatorialAlgebra.__init__(self, R)

class MacdonaldPolynomial_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    def scalar_qt(self, x):
        """
        Returns the qt-Hall scalar product of self and x by converting
        both to the power-sum basis.

        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: H([1]).scalar_qt(H([1]))
            (-q + 1)/(-t + 1)

        """
        P = self.parent()
        p = sfa.SFAPower(self.parent().base_ring())
        p_self = p(self)
        p_x = p(x)
        f = lambda part1, part2: part1.centralizer_size(t=P.t, q=P.q)
        return P._apply_multi_module_morphism(p_self, p_x, f, orthogonal=True)


    def omega_qt(self):
        """
        Returns the image of self under the $\omega_{qt}$ automorphism.

        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: H([1,1]).omega_qt()
            ((2*q^2-2*q*t-2*q+2*t)/(t^3-t^2-t+1))*McdH[1, 1] + ((q-1)/(t-1))*McdH[2]

        """
        P = self.parent()
        p = sfa.SFAPower(self.parent().base_ring())
        p_self = p(self)
        f = lambda part: (-1)**(part.size()-part.length())*prod([(1-P.q**i)/(1-P.t**i) for i in part])*p(part)
        return P( p._apply_module_morphism(p_self, f) )

    def nabla(self):
        """
        Returns the value of the nabla operator applied to self.  The
        eigenvectors of the nabla operator are the Macdonald polynomials
        in the Ht basis.

        EXAMPLES:
            sage: from sage.combinat.sf.macdonald import *
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: P([1,1]).nabla()
            ((q^2*t+q*t^2-2*t)/(q*t-1))*McdP[1, 1] + McdP[2]
        """
        parent = self.parent()
        Ht = MacdonaldPolynomials_ht(parent.base_ring(), q=parent.q, t=parent.t)
        return parent( Ht(self).nabla() )



#P basis
class MacdonaldPolynomials_p(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: P == loads(dumps(P))
            True
        """
        self._name = "Macdonald polynomials in the P basis"
        self._prefix = "McdP"
        self._element_class = MacdonaldPolynomial_p
        MacdonaldPolynomials_generic.__init__(self, R, q, t)
        self._J = MacdonaldPolynomialsJ(self.base_ring(), self.q, self.t)

        _set_cache(cache_p, self)


    def _coerce_start(self, x):
        """
        Coerce things into the P basis

        1) Q basis (proportional)
        2) J basis (proportional)
        3) Other symmetric functions

        EXAMPLES:
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: s = SFASchur(P.base_ring())

            sage: P._coerce_start(Q([2]))
            ((q*t^2-q*t-t+1)/(q^3-q^2-q+1))*McdP[2]
            sage: P._coerce_start(Q([2,1]))
            ((-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1)/(-q^4*t+2*q^3*t-q^2*t+q^2-2*q+1))*McdP[2, 1]

            sage: P(J([2]))
            (q*t^2-q*t-t+1)*McdP[2]
            sage: P(J([2,1]))
            (-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1)*McdP[2, 1]

            sage: P(s([2]))
            ((q-t)/(q*t-1))*McdP[1, 1] + McdP[2]
            sage: P(s([2,1]))
            ((q*t-t^2+q-t)/(q*t^2-1))*McdP[1, 1, 1] + McdP[2, 1]


        """
        BR = self.base_ring()

        #Macdonald Q-basis
        if isinstance(x, MacdonaldPolynomial_q):
            f = lambda m: 1/self(m).scalar_qt(self(m))
            return self._change_by_proportionality(x, f)
        elif isinstance(x, MacdonaldPolynomial_j):
            f = lambda m:  c2(m, self.q, self.t)
            return self._change_by_proportionality(x, f)
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._J( x ) )
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: P([1])^2 #indirect doctest
            ((q*t+q-t-1)/(q*t-1))*McdP[1, 1] + McdP[2]
        """
        return self( self._J(left)*self._J(right) )


class MacdonaldPolynomial_p(MacdonaldPolynomial_generic):
    def scalar_qt(self, x):
        """
        Returns the qt-Hall scalar product of self and x.  If x is
        in the Macdonald P or Q basis, then specialized code is used;
        otherwise, both are converted to the power-sums and the scalar
        product is carried out there.

        EXAMPLES:
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: a = P([2])
            sage: b = Q([2])
            sage: a.scalar_qt(a)
            (q^3 - q^2 - q + 1)/(q*t^2 - q*t - t + 1)
            sage: a.scalar_qt(b)
            1

        """
        if isinstance(x, MacdonaldPolynomial_p):
            P = self.parent()
            f = lambda part1, part2: c1(part1, P.q, P.t)/c2(part1, P.q, P.t)
            return P._apply_multi_module_morphism(self, x, f, orthogonal=True)
        elif isinstance(x, MacdonaldPolynomial_q):
            return x.scalar_qt(self)
        else:
            return MacdonaldPolynomial_generic.scalar_qt(self, x)


#Q basis
class MacdonaldPolynomials_q(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: Q == loads(dumps(Q))
            True
        """
        self._name = "Macdonald polynomials in the Q basis"
        self._prefix = "McdQ"
        self._element_class = MacdonaldPolynomial_q
        MacdonaldPolynomials_generic.__init__(self, R, q, t)

        self._J = MacdonaldPolynomialsJ(self.base_ring(), self.q, self.t)
        self._P = MacdonaldPolynomialsP(self.base_ring(), self.q, self.t)

        _set_cache(cache_q, self)


    def _coerce_start(self, x):
        """
        Coerce things into the Q basis

        1) P basis (proportional)
        2) J basis (proportional)
        3) Other symmetric functions

        EXAMPLES:
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: s = SFASchur(P.base_ring())

            sage: Q._coerce_start(J([2]))
            (q^3-q^2-q+1)*McdQ[2]

            sage: Q(P([2]))
            ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
            sage: P(Q(P([2])))
            McdP[2]
            sage: Q(P(Q([2])))
            McdQ[2]

            sage: Q(s([2]))
            ((q^2-q*t-q+t)/(t^3-t^2-t+1))*McdQ[1, 1] + ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
        """
        BR = self.base_ring()

        #Macdonald P-basis
        if isinstance(x, MacdonaldPolynomial_p):
            f = lambda m: self._P(m).scalar_qt(self._P(m))
            return self._change_by_proportionality(x, f)
        elif isinstance(x, MacdonaldPolynomial_j):
            return self._coerce_start( self._P(x) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._J( x ) )
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: J([1])^2 #indirect doctest
            ((q-1)/(q*t-1))*McdJ[1, 1] + ((t-1)/(q*t-1))*McdJ[2]
        """
        return self( self._J(left)*self._J(right) )

class MacdonaldPolynomial_q(MacdonaldPolynomial_generic):
    def scalar_qt(self, x):
        """
        Returns the qt-Hall scalar product of self and x.  If x is
        in the Macdonald P basis, then specialized code is used;
        otherwise, both are converted to the power-sums and the scalar
        product is carried out there.

        EXAMPLES:
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: a = Q([2])
            sage: a.scalar_qt(a)
            (-q*t^2 + q*t + t - 1)/(-q^3 + q^2 + q - 1)
            sage: a.scalar_qt(H([1,1]))
            t

        """
        if isinstance(x, MacdonaldPolynomial_p):
            P = self.parent()
            f = lambda part1, part2: 1
            return P._apply_multi_module_morphism(self, x, f, orthogonal=True)
        else:
            return MacdonaldPolynomial_generic.scalar_qt(self, x)


#J basis
j_to_s_cache = {}
s_to_j_cache = {}

class MacdonaldPolynomials_j(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: J == loads(dumps(J))
            True
        """
        self._name = "Macdonald polynomials in the J basis"
        self._prefix = "McdJ"
        self._element_class = MacdonaldPolynomial_j
        MacdonaldPolynomials_generic.__init__(self, R, q, t)

        self._s = sfa.SFASchur(self.base_ring())
        self._j_to_s_cache = j_to_s_cache
        self._s_to_j_cache = s_to_j_cache

        self._S = MacdonaldPolynomialsS(QQ)

        _set_cache(cache_j, self)


    def _coerce_start(self, x):
        """
        Coerce things into the J basis

        1) P basis (proportional)
        2) Q basis (proportional)
        4) Schur functions ( Lascoux, Lapointe and Morse creation operators on Macdonald polynomials )
        3) Everything else through s

        EXAMPLES:
            sage: P = MacdonaldPolynomialsP(QQ)
            sage: Q = MacdonaldPolynomialsQ(QQ)
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: s = SFASchur(P.base_ring())

            sage: J._coerce_start(P([2]))
            (1/(q*t^2-q*t-t+1))*McdJ[2]

            sage: J(Q([2]))
            (1/(q^3-q^2-q+1))*McdJ[2]


            sage: s(J([2]))
            (-q*t+t^2+q-t)*s[1, 1] + (q*t^2-q*t-t+1)*s[2]
            sage: J(s([2]))
            ((q-t)/(q*t^4-q*t^3-q*t^2-t^3+q*t+t^2+t-1))*McdJ[1, 1] + (1/(q*t^2-q*t-t+1))*McdJ[2]

        """
        BR = self.base_ring()

        #Macdonald P-basis
        if isinstance(x, MacdonaldPolynomial_p):
            f = lambda m: 1/c2(m, self.q, self.t)
            return self._change_by_proportionality(x, f)
        #Macdonald Q-basis
        elif isinstance(x, MacdonaldPolynomial_q):
            return self._coerce_start(x.parent()._P(x))
        #Everything else through s
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._s(x)
            return self._from_cache(x, self._s_cache, self._s_to_j_cache, q=self.q, t=self.t)
        else:
            raise TypeError

    def _s_cache(self, n):
        """
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the J basis and the Schur functions.

        EXAMPLES:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: J._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( J._s_to_j_cache[2] )
            [([1, 1], [([1, 1], 1/(t^3 - t^2 - t + 1))]),
             ([2],
              [([1, 1], (q - t)/(q*t^4 - q*t^3 - q*t^2 - t^3 + q*t + t^2 + t - 1)),
               ([2], 1/(q*t^2 - q*t - t + 1))])]
            sage: l( J._j_to_s_cache[2] )
            [([1, 1], [([1, 1], t^3 - t^2 - t + 1)]),
             ([2], [([1, 1], -q*t + t^2 + q - t), ([2], q*t^2 - q*t - t + 1)])]

        """
        self._invert_morphism(n, QQqt, self._j_to_s_cache, \
                              self._s_to_j_cache, to_other_function = self._to_s, \
                              upper_triangular=False)

    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of a partition in the Schur
        expansion of self(part).

        EXAMPLES:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: f21 = J._to_s(Partition([2,1]))
            sage: [f21(part) for part in Partitions(3)]
            [0,
             -q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1,
             q*t^3 - t^4 - q*t^2 + t^3 - q*t + t^2 + q - t]
        """
        q,t = QQqt.gens()
        res = self._S(1)
        for k in reversed(part):
            res = res.creation(k)
        res = res._omega_qt_in_schurs()
        res = res.map_coefficients(lambda c: c(t,q))
        f = lambda part2: res.coefficient(part2)
        return f


    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: J([1])^2 #indirect doctest
            ((q-1)/(q*t-1))*McdJ[1, 1] + ((t-1)/(q*t-1))*McdJ[2]
            sage: J([2])^2
            ((q^3-q^2-q+1)/(q^3*t^2-q^2*t-q*t+1))*McdJ[2, 2] + ((q^3*t-q^3+q^2*t-q^2-q*t+q-t+1)/(q^4*t^2-q^3*t-q*t+1))*McdJ[3, 1] + ((q*t^2-q*t-t+1)/(q^5*t^2-q^3*t-q^2*t+1))*McdJ[4]

        """
        return self( self._s(left)*self._s(right) )


class MacdonaldPolynomial_j(MacdonaldPolynomial_generic):
    def scalar_qt(self, x):
        """
        Returns the qt-Hall scalar product of self and x.  If x is
        in the Macdonald J basis, then specialized code is used;
        otherwise, both are converted to the power-sums and the scalar
        product is carried out there.

        EXAMPLES:
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: J([1,1]).scalar_qt(J([1,1]))
            q^2*t^4 - q^2*t^3 - q*t^4 - q^2*t^2 + q^2*t + 2*q*t^2 + t^3 - t^2 - q - t + 1
            sage: J([1,1]).scalar_qt(J([2]))
            0
            sage: J([2]).scalar_qt(J([2]))
            q^4*t^2 - q^4*t - q^3*t^2 - q^2*t^2 + q^3 + 2*q^2*t + q*t^2 - q^2 - q - t + 1

        """
        if isinstance(x, MacdonaldPolynomial_j):
            J = self.parent()
            f = lambda part1, part2: c1(part1, J.q, J.t)*c2(part2, J.q, J.t)
            return J._apply_multi_module_morphism(self, x, f, orthogonal=True)
        else:
            return MacdonaldPolynomial_generic.scalar_qt(self, x)



#H basis
h_to_s_cache = {}
s_to_h_cache = {}
class MacdonaldPolynomials_h(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: H == loads(dumps(H))
            True
        """
        self._name = "Macdonald polynomials in the H basis"
        self._prefix = "McdH"
        self._element_class = MacdonaldPolynomial_h
        MacdonaldPolynomials_generic.__init__(self, R, q, t)

        self._s = sfa.SFASchur(self.base_ring())
        self._self_to_s_cache = h_to_s_cache
        self._s_to_self_cache = s_to_h_cache

        self._Ht = MacdonaldPolynomialsHt(QQ)

        _set_cache(cache_h, self)


    def _coerce_start(self, x):
        """
        Coerce things into the H basis through the Schur functions.

        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: s = SFASchur(H.base_ring())
            sage: H._coerce_start(s([2]))
            ((-q)/(-q*t+1))*McdH[1, 1] + (1/(-q*t+1))*McdH[2]
        """
        BR = self.base_ring()
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._s(x)
            return self._from_cache(x, self._s_cache, self._s_to_self_cache, q=self.q, t=self.t)
        else:
            raise TypeError

    def _s_cache(self, n):
        """
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the H basis and the Schur functions.

        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: H._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( H._s_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/(-q*t + 1)), ([2], (-t)/(-q*t + 1))]),
             ([2], [([1, 1], (-q)/(-q*t + 1)), ([2], 1/(-q*t + 1))])]
            sage: l( H._self_to_s_cache[2] )
            [([1, 1], [([1, 1], 1), ([2], t)]), ([2], [([1, 1], q), ([2], 1)])]
        """
        self._invert_morphism(n, QQqt, self._self_to_s_cache, \
                              self._s_to_self_cache, to_other_function = self._to_s)

    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of a partition in the Schur
        expansion of self(part).

        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: f21 = H._to_s(Partition([2,1]))
            sage: [f21(part) for part in Partitions(3)]
            [t, q*t + 1, q]
        """
        q,t = QQqt.gens()
        s = sfa.SFASchur(QQqt)
        res = s(self._Ht(part)).map_coefficients(lambda c: c.subs(t=1/t))
        res *= t**part.weighted_size()
        f = lambda part2: res.coefficient(part2)
        return f

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: H = MacdonaldPolynomialsH(QQ)
            sage: H([1])^2 #indirect doctest
            ((q-1)/(q*t-1))*McdH[1, 1] + ((t-1)/(q*t-1))*McdH[2]
        """
        return self( self._s(left)*self._s(right) )

class MacdonaldPolynomial_h(MacdonaldPolynomial_generic):
    pass

#HTt basis
ht_to_s_cache = {}
s_to_ht_cache = {}
class MacdonaldPolynomials_ht(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: Ht == loads(dumps(Ht))
            True
        """
        self._name = "Macdonald polynomials in the Ht basis"
        self._prefix = "McdHt"
        self._element_class = MacdonaldPolynomial_ht
        MacdonaldPolynomials_generic.__init__(self, R, q, t)

        self._s = sfa.SFASchur(self.base_ring())
        self._self_to_s_cache = ht_to_s_cache
        self._s_to_self_cache = s_to_ht_cache
        self._J = MacdonaldPolynomialsJ(self.base_ring(), self.q, self.t)

        _set_cache(cache_ht, self)


    def _coerce_start(self, x):
        """
        Coerce things into the Ht basis though the Schur functions.

        EXAMPLES:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: s = SFASchur(Ht.base_ring())
            sage: Ht._coerce_start(s([2,1]))
            ((-q)/(-q*t^2+t^3+q^2-q*t))*McdHt[1, 1, 1] + ((q^2+q*t+t^2)/(-q^2*t^2+q^3+t^3-q*t))*McdHt[2, 1] + ((-t)/(q^3-q^2*t-q*t+t^2))*McdHt[3]
            sage: Ht._coerce_start(s([2]))
            ((-q)/(-q+t))*McdHt[1, 1] + (t/(-q+t))*McdHt[2]

        """
        BR = self.base_ring()
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._s(x)
            return self._from_cache(x, self._s_cache, self._s_to_self_cache, q=self.q, t=self.t)
        else:
            raise TypeError

    def _s_cache(self, n):
        """
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the Ht basis and the Schur functions.

        EXAMPLES:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: Ht._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( Ht._s_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/(-q + t)), ([2], (-1)/(-q + t))]),
             ([2], [([1, 1], (-q)/(-q + t)), ([2], t/(-q + t))])]
            sage: l( Ht._self_to_s_cache[2] )
            [([1, 1], [([1, 1], t), ([2], 1)]), ([2], [([1, 1], q), ([2], 1)])]
        """
        self._invert_morphism(n, QQqt, self._self_to_s_cache, \
                              self._s_to_self_cache, to_other_function = self._to_s)

    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of a partition in the Schur
        expansion of self(part).

        EXAMPLES:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: f21 = Ht._to_s(Partition([2,1]))
            sage: [f21(part) for part in Partitions(3)]
            [1, q + t, q*t]
            sage: f22 = Ht._to_s(Partition([2,2]))
            sage: [f22(part) for part in Partitions(4)]
            [1, q*t + q + t, q^2 + t^2, q^2*t + q*t^2 + q*t, q^2*t^2]

        """
        q,t = QQqt.gens()
        s = sfa.SFASchur(QQqt)
        J = MacdonaldPolynomialsJ(QQ)
        res = 0
        for p in sage.combinat.partition.Partitions(sum(part)):
            res += (J(part).scalar_t(s(p), t))*s(p)
        res = res.map_coefficients(lambda c: c.subs(t=1/t))
        res *= t**part.weighted_size()
        f = lambda part2: res.coefficient(part2)
        return f

    def _multiply(self, left, right):
        """
        EXAMPLES:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: Ht([1])^2 #indirect doctest
            ((-q+1)/(-q+t))*McdHt[1, 1] + ((t-1)/(-q+t))*McdHt[2]
        """
        return self( self._s(left)*self._s(right) )

class MacdonaldPolynomial_ht(MacdonaldPolynomial_generic):
    def nabla(self):
        """
        Returns the value of the nabla operator applied to self.  The
        eigenvectors of the nabla operator are the Macdonald polynomials
        in the Ht basis.

        EXAMPLES:
            sage: Ht = MacdonaldPolynomialsHt(QQ)
            sage: t = Ht.t; q = Ht.q;
            sage: a = sum(Ht(p) for p in Partitions(3))
            sage: a.nabla() == t^3*Ht([1,1,1])+q*t*Ht([2,1]) + q^3*Ht([3])
            True
        """
        Ht = self.parent()
        f = lambda part: (Ht.t)**(part.weighted_size())*(Ht.q)**(part.conjugate().weighted_size())*Ht(part)
        return Ht._apply_module_morphism(self, f)

#S basis
S_to_s_cache = {}
s_to_S_cache = {}
class MacdonaldPolynomials_s(MacdonaldPolynomials_generic):
    def __init__(self, R, q=None, t=None):
        """
        TESTS:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: S == loads(dumps(S))
            True
        """
        self._name = "Macdonald polynomials in the S basis"
        self._prefix = "McdS"
        self._element_class = MacdonaldPolynomial_s
        MacdonaldPolynomials_generic.__init__(self, R, q, t)

        self._s = sfa.SFASchur(self.base_ring())
        self._self_to_s_cache = S_to_s_cache
        self._s_to_self_cache = s_to_S_cache

        _set_cache(cache_s, self)



    def _multiply(self, left, right):
        """
        The multiplication of the modified Schur functions behaves the
        same as the multiplication of the Schur functions.

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: S([2])^2 #indirect doctest
            McdS[2, 2] + McdS[3, 1] + McdS[4]
        """
        s_left = self._s._from_element(left)
        s_right = self._s._from_element(right)
        product = s_left*s_right
        return self._from_element(product)


    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of a partition in the Schur
        expansion of self(part).

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: S2 = S._to_s(Partition([2]))
            sage: S2(Partition([2]))
            -t + 1
            sage: S2(Partition([1,1]))
            t^2 - t

        """
        #Covert to the power sum
        p = sfa.SFAPower(QQqt)
        s = sfa.SFASchur(QQqt)
        p_x = p(s(part))
        t = self.t
        f = lambda m, c: (m, c*prod([(1-t**k) for k in m]))
        res = s(p_x.map_mc(f))
        f = lambda part2: res.coefficient(part2)
        return f

    def _s_cache(self, n):
        """
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the S basis and the Schur functions.

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: S._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( S._s_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/(t^3 - t^2 - t + 1)), ([2], t/(t^3 - t^2 - t + 1))]), ([2], [([1, 1], t/(t^3 - t^2 - t + 1)), ([2], (-1)/(-t^3 + t^2 + t - 1))])]
            sage: l( S._self_to_s_cache[2] )
            [([1, 1], [([1, 1], -t + 1), ([2], t^2 - t)]), ([2], [([1, 1], t^2 - t), ([2], -t + 1)])]

        """
        self._invert_morphism(n, QQqt, self._self_to_s_cache, \
                              self._s_to_self_cache, to_other_function = self._to_s)


    def _coerce_start(self, x):
        """
        Coerce things into the S basis though the Schur functions.

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: J = MacdonaldPolynomialsJ(QQ)
            sage: S._coerce_start(J([2]))
            q*McdS[1, 1] + ((-1)/(-1))*McdS[2]
            sage: S._coerce_start(J([1,1]))
            McdS[1, 1] + t*McdS[2]
        """
        BR = self.base_ring()
        if isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            x = self._s(x)
            return self._from_cache(x, self._s_cache, self._s_to_self_cache, q=self.q, t=self.t)
        else:

            raise TypeError

class MacdonaldPolynomial_s(MacdonaldPolynomial_generic):
    def omega_qt(self):
        """
        Returns the image of self under the Frobenius / omega automorphism.

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: S([1,1]).omega_qt()
            (t/(t^3-t^2-t+1))*McdS[1, 1] + ((-1)/(-t^3+t^2+t-1))*McdS[2]

        """
        S = self.parent()
        return S( self._omega_qt_in_schurs() )

    def _creation_by_determinant_helper(self, k, part):
        """
        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: a = S([2,1])
            sage: a._creation_by_determinant_helper(2,[1])
            (q^3*t-q^2*t-q+1)*McdS[2, 1] + (q^3-q^2*t-q+t)*McdS[3]

        """
        S = self.parent()
        q,t = S.q, S.t

        part += [0]*(k-len(part))

        if len(part) > k:
            raise ValueError, "the column to add is too small"

        #Create the matrix over the homogeneous symmetric
        #functions and take its determinant
        MS = MatrixSpace(sfa.SFAHomogeneous(S.base_ring()), k, k)
        h  = MS.base_ring()
        m = []
        for i in range(k):
            row = [0]*max(0, (i+1)-2-part[i])
            for j in range(max(0, (i+1)-2-part[i]),k):
                value = part[i]+j-i+1
                p = [value] if value > 0 else []
                row.append( (1-q**(part[i]+j-i+1)*t**(k-(j+1)))*h(p) )
            m.append(row)
        M = MS(m)
        res = M.det()

        #Convert to the Schurs
        res = S._s( res )
        return S._from_element(res)

    def _creation_by_determinant(self, k):
        """
        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: a = S(1)
            sage: a._creation_by_determinant(1)
            (-q+1)*McdS[1]
            sage: a._creation_by_determinant(2)
            (q^2*t-q*t-q+1)*McdS[1, 1] + (q^2-q*t-q+t)*McdS[2]
        """
        S = self.parent()
        f = functools.partial(self._creation_by_determinant_helper,k)
        return S._apply_module_morphism(self, f)



    def _omega_qt_in_schurs(self):
        """
        Returns the image of self under the omega_qt automorphism in
        the Schur basis.

        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: a = S([2,1]) + S([1,1,1])
            sage: a._omega_qt_in_schurs()
            s[2, 1] + s[3]

        """
        S = self.parent()
        f = lambda part: S._s(part.conjugate())
        return S._s._apply_module_morphism(self, f)

    def creation(self, k):
        """
        EXAMPLES:
            sage: S = MacdonaldPolynomialsS(QQ)
            sage: a = S(1)
            sage: a.creation(1)
            (-q+1)*McdS[1]
            sage: a.creation(2)
            (q^2*t-q*t-q+1)*McdS[1, 1] + (q^2-q*t-q+t)*McdS[2]
        """
        return self._creation_by_determinant(k)


###############
_qt_kostka_cache = {}
def qt_kostka(lam, mu):
    """
    Returns the $K_{\lambda\mu}(q,t)$ by computing the change of
    basis from the Macdonald H basis to the Schurs.

    EXAMPLES:
        sage: from sage.combinat.sf.macdonald import qt_kostka
        sage: qt_kostka([2,1,1],[1,1,1,1])
        t^3 + t^2 + t
        sage: qt_kostka([1,1,1,1],[2,1,1])
        q
        sage: qt_kostka([1,1,1,1],[3,1])
        q^3
        sage: qt_kostka([1,1,1,1],[1,1,1,1])
        1
        sage: qt_kostka([2,1,1],[2,2])
        q^2*t + q*t + q
        sage: qt_kostka([2,2],[2,2])
        q^2*t^2 + 1
        sage: qt_kostka([4],[3,1])
        t
        sage: qt_kostka([2,2],[3,1])
        q^2*t + q
        sage: qt_kostka([3,1],[2,1,1])
        q*t^3 + t^2 + t
        sage: qt_kostka([2,1,1],[2,1,1])
        q*t^2 + q*t + 1
        sage: qt_kostka([2,1],[1,1,1,1])
        0

    """
    lam = sage.combinat.partition.Partition(lam)
    mu = sage.combinat.partition.Partition(mu)
    R = QQ['q,t']

    if lam.size() != mu.size():
        return R(0)

    if (lam,mu) in _qt_kostka_cache:
        return _qt_kostka_cache[(lam,mu)]

    H = MacdonaldPolynomialsH(QQ)
    s = sfa.SFASchur(H.base_ring())

    parts = sage.combinat.partition.Partitions(mu.size())

    for p2 in parts:
        res = s(H(p2))
        for p1 in parts:
            _qt_kostka_cache[(p1,p2)] = R(res.coefficient(p1).numerator())

    return _qt_kostka_cache[(lam,mu)]



def _set_cache(c, self):
    """
    Since the Macdonald polynomials could be called many ways,
    this is a utility routine that adds those other ways to
    the cache.

    EXAMPLES:
        sage: S = MacdonaldPolynomialsS(QQ)
        sage: R = QQ['q,t'].fraction_field()
        sage: q,t = R.gens()
        sage: S2 = MacdonaldPolynomialsS(R,q,t)
        sage: S2 is S #indirect doctest
        True

    """
    br = self.base_ring()
    q,t = self.q, self.t
    c[c.key(br, q=q, t=t)] = self
    c[c.key(br, q, t)] = self
    c[c.key(br, q, t)] = self
    c[c.key(br, q, t=t)] = self



#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_p = Cache(MacdonaldPolynomials_p)
cache_j = Cache(MacdonaldPolynomials_j)
cache_q = Cache(MacdonaldPolynomials_q)
cache_h = Cache(MacdonaldPolynomials_h)
cache_ht = Cache(MacdonaldPolynomials_ht)
cache_s = Cache(MacdonaldPolynomials_s)
