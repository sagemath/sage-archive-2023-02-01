"""
Hall-Littlewood Polynomials
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

from sage.combinat.combinat import CombinatorialClass
from sage.libs.symmetrica.all import hall_littlewood
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
import sfa
import sage.combinat.partition
import kfpoly
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ, QQ
from sage.misc.misc import prod


##################################
#Still under major development!!!#
##################################

def HallLittlewoodP(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    P basis.  This is the same as the HL basis in John Stembridge's
    SF examples file.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodP(QQ)
        Hall-Littlewood polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodP(QQ, t=-1)
        Hall-Littlewood polynomials in the P basis with t=-1 over Rational Field
        sage: HLP = HallLittlewoodP(QQ)
        sage: s = SFASchur(HLP.base_ring())
        sage: s(HLP([2,1]))
        (-t^2-t)*s[1, 1, 1] + s[2, 1]

      The Hall-Littlewood polynomials in the P basis at t = 0 are the Schur
      functions.
        sage: HLP = HallLittlewoodP(QQ,t=0)
        sage: s = SFASchur(HLP.base_ring())
        sage: s(HLP([2,1])) == s([2,1])
        True

      The Hall-Littlewood polynomials in the P basis at t = 1 are the monomial
      symmetric functions.
        sage: HLP = HallLittlewoodP(QQ,t=1)
        sage: m = SFAMonomial(HLP.base_ring())
        sage: m(HLP([2,2,1])) == m([2,2,1])
        True
    """
    return cache_p(R,t)

def HallLittlewoodQ(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    Q basis.  This is the same as the Q basis in John Stembridge's
    SF examples file.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodQ(QQ)
        Hall-Littlewood polynomials in the Q basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodQ(QQ, t=-1)
        Hall-Littlewood polynomials in the Q basis with t=-1 over Rational Field

    """
    return cache_q(R,t)

def HallLittlewoodQp(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood
    Qp basis.  This is dual to the Hall-Littlewood P basis with
    respect to the standard scalar product.

    If t is not specified, then the base ring will be obtained by making the univariate
    polynomial ring over R with the variable t and taking its fraction field.

    EXAMPLES:
        sage: HallLittlewoodQp(QQ)
        Hall-Littlewood polynomials in the Qp basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodQp(QQ, t=-1)
        Hall-Littlewood polynomials in the Qp basis with t=-1 over Rational Field

    """
    return cache_qp(R,t)


##################################

class HallLittlewood_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, t=None):
        """
        TESTS:
            sage: HallLittlewoodP(QQ)
            Hall-Littlewood polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: HallLittlewoodP(QQ,t=2)
            Hall-Littlewood polynomials in the P basis with t=2 over Rational Field
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

        CombinatorialAlgebra.__init__(self, R)


    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix between self and basis
        for the homogenous component of degree n.

        EXAMPLES:
            sage: HLP = HallLittlewoodP(QQ)
            sage: s   = SFASchur(HLP.base_ring())
            sage: HLP.transition_matrix(s, 4)
            [             1             -t              0            t^2           -t^3]
            [             0              1             -t             -t      t^3 + t^2]
            [             0              0              1             -t            t^3]
            [             0              0              0              1 -t^3 - t^2 - t]
            [             0              0              0              0              1]
            sage: HLQ = HallLittlewoodQ(QQ)
            sage: HLQ.transition_matrix(s,3)
            [                        -t + 1                        t^2 - t                     -t^3 + t^2]
            [                             0                  t^2 - 2*t + 1           -t^4 + t^3 + t^2 - t]
            [                             0                              0 -t^6 + t^5 + t^4 - t^2 - t + 1]
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp.transition_matrix(s,3)
            [      1       0       0]
            [      t       1       0]
            [    t^3 t^2 + t       1]


        """
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)

class HallLittlewoodElement_generic(sfa.SymmetricFunctionAlgebraElement_generic):
    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in n variables.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).expand(2)
            x0^2 + (-t + 1)*x0*x1 + x1^2
            sage: HLQ([2]).expand(2)
            (-t + 1)*x0^2 + (t^2 - 2*t + 1)*x0*x1 + (-t + 1)*x1^2
            sage: HLQp([2]).expand(2)
            x0^2 + x0*x1 + x1^2
        """
        sp = self.parent()
        BR = sp.base_ring()
        s = sfa.SFASchur(BR)
        return s(self).expand(n, alphabet=alphabet)


    def scalar(self, x):
        """
        Returns standard scalar product between self and s.

        This is the default implementation that converts both self
        and x into Schur functions and performs the scalar product
        that basis.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).scalar(HLQp([2]))
            1
            sage: HLP([2]).scalar(HLQp([1,1]))
            0
        """
        sp = self.parent()
        xp = x.parent()
        BR = sp.base_ring()

        s = sfa.SFASchur(BR)
        s_self = s(self)
        s_x = s(x)
        return s_self.scalar(s_x)


    def scalar_hl(self, x, t=None):
        """
        Returns the standard Hall-Littlewood scalar product of self and x.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLP([2]).scalar_hl(HLQ([2]))
            1
            sage: HLP([2]).scalar_hl(HLQ([1,1]))
            0
        """
        parent = self.parent()
        p = sfa.SFAPower(parent.base_ring())
        f = lambda part1, part2: part1.centralizer_size(t=parent.t)
        return parent._apply_multi_module_morphism(p(self),p(x),f,orthogonal=True)



###########
# P basis #
###########
p_to_s_cache = {}
s_to_p_cache = {}

class HallLittlewoodElement_p(HallLittlewoodElement_generic):
    pass

class HallLittlewood_p(HallLittlewood_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: P = HallLittlewoodP(QQ)
            sage: P == loads(dumps(P))
            True
        """
        self._name = "Hall-Littlewood polynomials in the P basis"
        self._prefix = "P"
        self._element_class = HallLittlewoodElement_p

        HallLittlewood_generic.__init__(self, R, t=t)

        self._s = sfa.SFASchur(self.base_ring())

        self._self_to_s_cache = p_to_s_cache
        self._s_to_self_cache = s_to_p_cache



    def _multiply(self, left, right):
        """
        Convert to the Schur basis, do the multiplication there, and
        convert back to the P basis.

        EXAMPLES:
            sage: HLP = HallLittlewoodP(QQ)
            sage: HLP([2])^2
            (t+1)*P[2, 2] + (-t+1)*P[3, 1] + P[4]

        """
        return self( self._s(left) * self._s(right) )


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood P basis.

        1) Hall-Littlewood polynomials in the Q basis
        2) Hall-Littlewood polynomials in the Qp basis (via the Schurs)
        3) Classical symmetric functions

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLP(HLQ([2]))
            (-t+1)*P[2]
            sage: HLP(HLQp([2]))
            t*P[1, 1] + P[2]
            sage: HLP(s([2]))
            t*P[1, 1] + P[2]
            sage: HLP(p([2]))
            (t-1)*P[1, 1] + P[2]

        """
        BR = self.base_ring()
        if isinstance(x, HallLittlewoodElement_q):
            return self._change_by_proportionality(x, self._q_to_self)
        elif isinstance(x, HallLittlewoodElement_qp):
            return self( self._s(x) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            #Convert x to the Schur basis
            x = self._s(x)
            return self._from_cache(x, self._s_cache, self._s_to_self_cache,t=self.t)
        else:
            raise TypeError

    def _q_to_self(self, m):
        """
        Returns the scalar coefficient on self(m) when converting from the
        Q basis to the P basis.  Note that this assumes that m is a
        Partition object.

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLP._q_to_self(Partition([2,1]))
            t^2 - 2*t + 1

        """
        t = self.t
        coeff = (1-t)**len(m)
        for i in m.to_exp():
            for j in range(1,i+1):
                coeff *= (1-t**j)/(1-t)
        return coeff

    def _s_to_self(self, part):
        """
        Returns a function which gives the coefficient of part2 in the expansion
        of the Schur functions s(part) in self.

        EXAMPLES:
            sage: P = HallLittlewoodP(QQ)
            sage: f21 = P._s_to_self(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [0, 1, t^2 + t]
        """
        Zt = ZZ['t']
        t = Zt.gen()
        zero = Zt(0)
        res_dict = kfpoly.schur_to_hl(part, t)
        f = lambda part2: res_dict.get(part2,zero)
        return f


    def _s_cache(self, n):
        """
        Computes the change of basis between the P polynomials and
        the Schur functions for partitions of size n.

        Uses the fact that the transformation matrix is upper-triangular
        in order to obtain the inverse transformation.

        EXAMPLES:
            sage: P = HallLittlewoodP(QQ)
            sage: P._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(P._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], t), ([2], 1)])]
            sage: l(P._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -t), ([2], 1)])]

        """

        self._invert_morphism(n, ZZ['t'], self._self_to_s_cache, \
                              self._s_to_self_cache, to_self_function = self._s_to_self, \
                              upper_triangular=True, ones_on_diagonal=True)





###########
# Q basis #
###########
class HallLittlewoodElement_q(HallLittlewoodElement_generic):
    pass

class HallLittlewood_q(HallLittlewood_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: Q = HallLittlewoodQ(QQ)
            sage: Q == loads(dumps(Q))
            True
        """
        self._name = "Hall-Littlewood polynomials in the Q basis"
        self._prefix = "Q"
        self._element_class = HallLittlewoodElement_q

        HallLittlewood_generic.__init__(self, R, t=t)

        self._P = HallLittlewood_p(R, t=t)


    def _multiply(self, left, right):
        """
        Converts to the P basis, does the multiplication there, and
        converts back to the Q basis.

        EXAMPLES:
            sage: HLQ = HallLittlewoodQ(QQ)
            sage: HLQ([2])^2
            Q[2, 2] + (-t+1)*Q[3, 1] + (-t+1)*Q[4]

        """
        return self( self._P(left) * self._P(right) )


    def _coerce_start(self, x):
        """
        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQ( HLP([2,1]) + HLP([3]) )
            (1/(t^2-2*t+1))*Q[2, 1] + (1/(-t+1))*Q[3]
            sage: HLQ(HLQp([2]))
            (-t/(-t^3+t^2+t-1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(s([2]))
            (-t/(-t^3+t^2+t-1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(p([2]))
            (-1/(-t^2+1))*Q[1, 1] + (1/(-t+1))*Q[2]
        """
        if isinstance(x, HallLittlewoodElement_p):
            return self._change_by_proportionality(x, self._p_to_self)
        elif isinstance(x, HallLittlewoodElement_qp):
            return self( self._P(x) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            return self( self._P(x) )
        else:
            raise TypeError

    def _p_to_self(self, m):
        """
        Returns the scalar coefficient on self(m) when converting from the
        Q basis to the P basis.  Note that this assumes that m is a
        Partition object.

        EXAMPLES:
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQ._p_to_self(Partition([2,1]))
            1/(t^2 - 2*t + 1)

        """
        t = self.t
        coeff = 1/(1-t)**len(m)
        for i in m.to_exp():
            for j in range(1,i+1):
                coeff *= (1-t)/(1-t**j)
        return coeff

############
# Qp basis #
############
qp_to_s_cache = {}
s_to_qp_cache = {}

class HallLittlewoodElement_qp(HallLittlewoodElement_generic):
    pass


class HallLittlewood_qp(HallLittlewood_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES:
            sage: Qp = HallLittlewoodQp(QQ)
            sage: Qp == loads(dumps(Qp))
            True
        """
        self._name = "Hall-Littlewood polynomials in the Qp basis"
        self._prefix = "Qp"
        self._element_class = HallLittlewoodElement_qp

        HallLittlewood_generic.__init__(self,R, t=t)

        self._s = sfa.SFASchur(self.base_ring())

        self._self_to_s_cache = qp_to_s_cache
        self._s_to_self_cache = s_to_qp_cache


    def _coerce_start(self, x):
        """
        Coerce things into the Hall-Littlewood Qp basis.

        1) Hall-Littlewood polynomials in the Q basis (via the Schurs)
        2) Hall-Littlewood polynomials in the P basis (via the Schurs)
        3) Symmetric Functions

        EXAMPLES:
            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQp(HLP([2]))
            -t*Qp[1, 1] + (t^2+1)*Qp[2]
            sage: HLQp(HLQ([2]))
            (t^2-t)*Qp[1, 1] + (-t^3+t^2-t+1)*Qp[2]
            sage: HLQp(s([2]))
            Qp[2]
            sage: HLQp(p([2]))
            -Qp[1, 1] + (t+1)*Qp[2]
        """
        if isinstance(x, HallLittlewoodElement_q):
            return self( self._s( x ) )
        elif isinstance(x, HallLittlewoodElement_p):
            return self( self._s( x ) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebraElement_generic):
            sx = self._s( x )
            return self._from_cache(sx, self._s_cache, self._s_to_self_cache,t=self.t)
        else:
            raise TypeError

    def _multiply(self, left, right):
        """
        Converts the Hall-Littlewood polynomial in the Qp basis
        to a Schur function, performs the multiplication there,
        and converts it back to the Qp basis.

        EXAMPLES:
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp([2])^2
            Qp[2, 2] + (-t+1)*Qp[3, 1] + (-t+1)*Qp[4]

        """
        return self( self._s(left)*self._s(right) )

    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of part2 in the Schur
        expansion of self(part).

        EXAMPLES:
            sage: Qp = HallLittlewoodQp(QQ)
            sage: f21 = Qp._to_s(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [t, 1, 0]
        """
        Zt = ZZ['t']
        t = Zt.gen()
        zero = Zt(0)

        if part == []:
            return lambda part2: Zt(1)

        res = hall_littlewood(part)
        f = lambda part2: res.coefficient(part2).subs(x=t)
        return f


    def _s_cache(self, n):
        """
        Computes the change of basis between the Qp polynomials and
        the Schur functions for partitions of size n.

        Uses the fact that the transformation matrix is lower-triangular
        in order to obtain the inverse transformation.

        EXAMPLES:
            sage: Qp = HallLittlewoodQp(QQ)
            sage: Qp._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(Qp._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1), ([2], -t)]), ([2], [([2], 1)])]
            sage: l(Qp._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1), ([2], t)]), ([2], [([2], 1)])]

        """
        self._invert_morphism(n, ZZ['t'], self._self_to_s_cache, \
                              self._s_to_self_cache, to_other_function = self._to_s, \
                              lower_triangular=True, ones_on_diagonal=True)



#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_p = Cache(HallLittlewood_p)
cache_q = Cache(HallLittlewood_q)
cache_qp = Cache(HallLittlewood_qp)

