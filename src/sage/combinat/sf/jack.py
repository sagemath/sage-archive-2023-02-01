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
import sage.categories.all
from sage.rings.all import Integer, gcd, lcm, QQ, is_FractionField
from sage.misc.misc import prod
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom, End
import sfa

def JackPolynomialsP(R, t=None):
    """
    Returns the algebra of Jack polynomials in the P basis.

    If t is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over R with the variable t
    and taking its fraction field.

    EXAMPLES::

        sage: JackPolynomialsP(QQ)
        Jack polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsP(QQ,t=-1)
        Jack polynomials in the P basis with t=-1 over Rational Field

    At t = 1, the Jack polynomials on the P basis are the Schur
    symmetric functions.

    ::

        sage: P = JackPolynomialsP(QQ,1)
        sage: s = SFASchur(QQ)
        sage: P([2,1])^2
        JackP[2, 2, 1, 1] + JackP[2, 2, 2] + JackP[3, 1, 1, 1] + 2*JackP[3, 2, 1] + JackP[3, 3] + JackP[4, 1, 1] + JackP[4, 2]
        sage: s([2,1])^2
        s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2]

    At t = 2, the Jack polynomials on the P basis are the zonal
    polynomials.

    ::

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

    If t is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over R with the variable t
    and taking its fraction field.

    EXAMPLES::

        sage: JackPolynomialsQ(QQ)
        Jack polynomials in the Q basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsQ(QQ,t=-1)
        Jack polynomials in the Q basis with t=-1 over Rational Field
    """
    return cache_q(R,t)

def JackPolynomialsJ(R, t=None):
    """
    Returns the algebra of Jack polynomials in the J basis.

    If t is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over R with the variable t
    and taking its fraction field.

    EXAMPLES::

        sage: JackPolynomialsJ(QQ)
        Jack polynomials in the J basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: JackPolynomialsJ(QQ,t=-1)
        Jack polynomials in the J basis with t=-1 over Rational Field

    At t = 1, the Jack polynomials in the J basis are scalar multiples
    of the Schur functions with the scalar given by a Partition's
    hook_product method at 1.

    ::

        sage: J = JackPolynomialsJ(QQ, t=1)
        sage: s = SFASchur(J.base_ring())
        sage: p = Partition([3,2,1,1])
        sage: s(J(p)) == p.hook_product(1)*s(p)
        True

    At t = 2, the Jack polynomials on the J basis are scalar multiples
    of the zonal polynomials with the scalar given by a Partition's
    hook_product method at 1.

    ::

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
    Returns the algebra of Jack polynomials in the Qp, which is dual to
    the P basis with respect to the standard scalar product.

    If t is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over R with the variable t
    and taking its fraction field.

    EXAMPLES::

        sage: P = JackPolynomialsP(QQ)
        sage: Qp = JackPolynomialsQp(QQ)
        sage: a = Qp([2])
        sage: a.scalar(P([2]))
        1
        sage: a.scalar(P([1,1]))
        0
        sage: P(Qp([2]))                        # todo: missing auto normalization
        ((2*t-2)/(2*t+2))*JackP[1, 1] + JackP[2]
        sage: P._normalize(P(Qp([2])))
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

    EXAMPLES::

        sage: Z = ZonalPolynomials(QQ)
        sage: a = Z([2])
        sage: Z([2])^2
        64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]
    """
    return cache_z(R)

###################################################################
def c1(part, t):
    """
    EXAMPLES::

        sage: from sage.combinat.sf.jack import c1
        sage: t = QQ['t'].gen()
        sage: [c1(p,t) for p in Partitions(3)]
        [2*t^2 + 3*t + 1, t + 2, 6]
    """
    return prod([1+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))],
                t.parent()(1)) # FIXME: use .one()

def c2(part, t):
    """
    EXAMPLES::

        sage: from sage.combinat.sf.jack import c2
        sage: t = QQ['t'].gen()
        sage: [c2(p,t) for p in Partitions(3)]
        [6*t^3, 2*t^3 + t^2, t^3 + 3*t^2 + 2*t]
    """
    return prod([t+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))],
                t.parent()(1)) # FIXME: use .one()

####################################################################

class JackPolynomials_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES::

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

        sfa.SymmetricFunctionAlgebra_generic.__init__(self, R)

        # Bases defined by orthotriangularity should inherit from some
        # common category BasesByOrthotriangularity (shared with Jack, HL, orthotriang, Mcdo)
        if hasattr(self, "_m_cache"):
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = sage.categories.all.ModulesWithBasis(R)
            self._m = sfa.SFAMonomial(self.base_ring())
            self   .register_coercion(SetMorphism(Hom(self._m, self, category), self._m_to_self))
            self._m.register_coercion(SetMorphism(Hom(self, self._m, category), self._self_to_m))

    def _m_to_self(self, x):
        """
        Isomorphism from the monomial basis into self

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ,t=2)
            sage: m = SFAMonomial(P.base_ring())
            sage: P._m_to_self(m[2,1])
            -3/2*JackP[1, 1, 1] + JackP[2, 1]

        This is for internal use only. Please use instead::

            sage: P(m[2,1])
            -3/2*JackP[1, 1, 1] + JackP[2, 1]
        """
        return self._from_cache(x, self._m_cache, self._m_to_self_cache)

    def _self_to_m(self, x):
        r"""
        Isomorphism from self to the monomial basis

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ,t=2)
            sage: m = SFAMonomial(P.base_ring())
            sage: P._self_to_m(P[2,1])
            3/2*m[1, 1, 1] + m[2, 1]

        This is for internal use only. Please use instead::

            sage: m(P[2,1])
            3/2*m[1, 1, 1] + m[2, 1]
        """
        return self._m._from_cache(x, self._m_cache, self._self_to_m_cache)

    def c1(self, part):
        """
        Returns the t-Jack scalar product between J(part) and P(part).

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: P.c1(Partition([2,1]))
            t + 2
        """
        return c1(part, self.t)

    def c2(self, part):
        """
        Returns the t-Jack scalar product between J(part) and Q(part).

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: P.c2(Partition([2,1]))
            2*t^3 + t^2
        """
        return c2(part, self.t)

    def _normalize_coefficients(self, c):
        """
        If our coefficient ring is the field of fractions over a univariate
        polynomial ring over the rationals, then we should clear both the
        numerator and denominator of the denominators of their
        coefficients.

        EXAMPLES::

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

    def _normalize(self, x):
        """
        Normalize the coefficients of `x`

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: t = P.base_ring().gen()
            sage: a = 2/(1/2*t+1/2)
            sage: b = 1/(1/3+1/6*t)
            sage: c = 24/(4*t^2 + 12*t + 8)
            sage: P._normalize( a*P[1] + b*P[2] + c*P[2,1] )
            (4/(t+1))*JackP[1] + (6/(t+2))*JackP[2] + (6/(t^2+3*t+2))*JackP[2, 1]

        TODO: this should be a method on the elements (what's the
        standard name for such methods?)
        """
        return x.map_coefficients(self._normalize_coefficients)

    def _normalize_morphism(self, category):
        """
        Returns the normalize morphism

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ); P.rename("JackP")
            sage: normal = P._normalize_morphism(AlgebrasWithBasis(P.base_ring()))
            sage: normal.parent()
            Set of Homomorphisms from JackP to JackP
            sage: normal.category_for()
            Category of algebras with basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field

            sage: t = P.t
            sage: a = 2/(1/2*t+1/2)
            sage: b = 1/(1/3+1/6*t)
            sage: c = 24/(4*t^2 + 12*t + 8)
            sage: normal( a*P[1] + b*P[2] + c*P[2,1] )
            (4/(t+1))*JackP[1] + (6/(t+2))*JackP[2] + (6/(t^2+3*t+2))*JackP[2, 1]

        TODO: this method should not be needed once short idioms to
        construct morphisms will be available
        """
        return SetMorphism(End(self, category), self._normalize)

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
      def scalar_jack(self, x):
        """
        EXAMPLES::

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
    Returns the Jack scalar product between p(part1) and p(part2) where
    p is the power-sum basis.

    EXAMPLES::

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

class JackPolynomials_p(JackPolynomials_generic):

    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: P == loads(dumps(P))
            True
        """
        self._name = "Jack polynomials in the P basis"
        self._prefix = "JackP"

        self._m_to_self_cache = {}
        self._self_to_m_cache = {}
        JackPolynomials_generic.__init__(self, R, t=t)

    def _coerce_start_disabled(self, x):
        """
        Coerce things into the Jack polynomials P basis.

        EXAMPLES::

            sage: Q = JackPolynomialsQ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: J = JackPolynomialsJ(QQ)

        ::

            sage: P(Q([2,1])) # indirect doctest
            ((t+2)/(2*t^3+t^2))*JackP[2, 1]
            sage: P(Q([3]))
            ((2*t^2+3*t+1)/(6*t^3))*JackP[3]
            sage: P(Q([1,1,1]))
            (6/(t^3+3*t^2+2*t))*JackP[1, 1, 1]

        ::

            sage: P(J([3]))
            (2*t^2+3*t+1)*JackP[3]
            sage: P(J([2,1]))
            (t+2)*JackP[2, 1]
            sage: P(J([1,1,1]))
            6*JackP[1, 1, 1]

        ::

            sage: s = SFASchur(QQ) # todo: not implemented
            sage: s = SFASchur(P.base_ring())
            sage: P(s([2,1]))
            ((2*t-2)/(t+2))*JackP[1, 1, 1] + JackP[2, 1]
            sage: s(_)
            s[2, 1]
        """

    def _m_cache(self, n):
        """
        Computes the change of basis between the Jack polynomials in the P
        basis and the monomial symmetric functions. This uses Gram-Schmidt
        to go to the monomials, and then that matrix is simply inverted.

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: P._m_cache(2)
            sage: l(P._self_to_m_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -2/(-t - 1)), ([2], 1)])]
            sage: l(P._m_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 2/(-t - 1)), ([2], 1)])]
            sage: P._m_cache(3)
            sage: l(P._m_to_self_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], -6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], -6/(-t^2 - 3*t - 2)), ([2, 1], -3/(2*t + 1)), ([3], 1)])]
            sage: l(P._self_to_m_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], -6/(-2*t^2 - 3*t - 1)), ([2, 1], 3/(2*t + 1)), ([3], 1)])]
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
        coefficient of lambda in the expansion of self(part) in the
        monomial basis.

        This assumes that the cache from the Jack polynomials in the P
        basis to the monomial symmetric functions has already been
        computed.

        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: P._m_cache(2)
            sage: p2 = Partition([2])
            sage: p11 = Partition([1,1])
            sage: f = P._to_m(p2)
            sage: f(p11)
            -2/(-t - 1)
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
        EXAMPLES::

            sage: P = JackPolynomialsP(QQ)
            sage: P([1])^2 # indirect doctest
            (-2*t/(-t-1))*JackP[1, 1] + JackP[2]
            sage: P._m(_)
            2*m[1, 1] + m[2]
            sage: P = JackPolynomialsP(QQ, 2)
            sage: P([2,1])^2
            125/63*JackP[2, 2, 1, 1] + 25/12*JackP[2, 2, 2] + 25/18*JackP[3, 1, 1, 1] + 12/5*JackP[3, 2, 1] + 4/3*JackP[3, 3] + 4/3*JackP[4, 1, 1] + JackP[4, 2]
            sage: P._m(_)
            45*m[1, 1, 1, 1, 1, 1] + 51/2*m[2, 1, 1, 1, 1] + 29/2*m[2, 2, 1, 1] + 33/4*m[2, 2, 2] + 9*m[3, 1, 1, 1] + 5*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]
        """
        return self( self._m(left)*self._m(right) )


    def scalar_jack_basis(self, part1, part2 = None):
        """
        Returns the scalar product of P[part1] and P[part2].

        Todo: check all the results!

        EXAMPLES:
            sage: P = JackPolynomialsP(QQ)
            sage: P.scalar_jack_basis(Partition([2,1]), Partition([1,1,1]))
            0
            sage: P._normalize_coefficients(P.scalar_jack_basis(Partition([3,2,1]), Partition([3,2,1])))
            (12*t^6 + 20*t^5 + 11*t^4 + 2*t^3)/(2*t^3 + 11*t^2 + 20*t + 12)

        With a single argument, takes part2 = part1
            sage: P.scalar_jack_basis(Partition([2,1]), Partition([2,1]))
            (2*t^3 + t^2)/(t + 2)

        NT: those results do not quite with Macdonald Symmetric
        Function and Orthogonal Polynomials p.12 (11.3). Is this P
        basis a normalization variant of that of Macdo?
        """

        if part2 is not None and part1 != part2:
            return self.base_ring().zero()
        return self.c2(part1) / self.c1(part1)


    class Element(JackPolynomials_generic.Element):
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
                return P._apply_multi_module_morphism(self, x, P.scalar_jack_basis, orthogonal=True)
            else:
                return JackPolynomials_generic.Element.scalar_jack(self, x)

#Zonal polynomials ( =P(2) )

class ZonalPolynomials(JackPolynomials_p):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: Z = ZonalPolynomials(QQ)
            sage: P = Z._P; P
            Jack polynomials in the P basis with t=2 over Rational Field
            sage: Z(P[2,1] + 2*P[3,1])
            Z[2, 1] + 2*Z[3, 1]
            sage: P(Z[2,1] + 2*Z[3,1])
            JackP[2, 1] + 2*JackP[3, 1]

        TESTS::

            sage: TestSuite(Z).run(elements = [Z[1], Z[1,1]])
            sage: TestSuite(Z).run(skip = ["_test_associativity", "_test_prod"])  # long time (7s on sage.math, 2011)

        Note: ``Z.an_element()`` is of degree 4; so we skip the
        ``_test_associativity`` and ``_test_prod`` which involve
        (currently?) expensive calculations up to degree 12.
        """
        JackPolynomials_p.__init__(self, R, t=R(2))
        self._name = "Zonal polynomials"
        self._prefix = "Z"
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        self._P = JackPolynomialsP(self.base_ring(), t=R(2))
        self   .register_coercion(SetMorphism(Hom(self._P, self, category), self.sum_of_terms))
        self._P.register_coercion(SetMorphism(Hom(self, self._P, category), self._P.sum_of_terms))

    class Element(JackPolynomials_p.Element):
        pass


#J basis

class JackPolynomials_j(JackPolynomials_generic):

    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: J = JackPolynomialsJ(QQ)
            sage: J == loads(dumps(J))
            True
        """
        self._name = "Jack polynomials in the J basis"
        self._prefix = "JackJ"
        JackPolynomials_generic.__init__(self, R, t=t)

        # Should be shared with _q (and possibly other bases in Macdo/HL) as BasesByRenormalization
        self._P = JackPolynomialsP(R, t)
        # temporary until Hom(GradedHopfAlgebrasWithBasis) works better
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        phi = self.module_morphism(diagonal = self.c1, codomain = self._P, category = category)
        # should use module_morphism(on_coeffs = ...) once it exists
        self._P.register_coercion(self._P._normalize_morphism(category) * phi)
        self   .register_coercion(self   ._normalize_morphism(category) *~phi)

    def _coerce_start_disabled(self, x):
        """
        EXAMPLES::

            sage: J = JackPolynomialsJ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: J(sum(P(p) for p in Partitions(3)))
            1/6*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]

        ::

            sage: s = SFASchur(J.base_ring())
            sage: J(s([3])) # indirect doctest
            ((-t^2+3*t-2)/(-6*t^2-18*t-12))*JackJ[1, 1, 1] + ((2*t-2)/(2*t^2+5*t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]
            sage: J(s([2,1]))
            ((t-1)/(3*t+6))*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1]
            sage: J(s([1,1,1]))
            1/6*JackJ[1, 1, 1]
        """

    def _multiply(self, left, right):
        """
        EXAMPLES::

            sage: J = JackPolynomialsJ(QQ)
            sage: J([1])^2 #indirect doctest
            (-t/(-t-1))*JackJ[1, 1] + (1/(t+1))*JackJ[2]
            sage: J([2])^2
            (-2*t^2/(-2*t^2-3*t-1))*JackJ[2, 2] + (-4*t/(-3*t^2-4*t-1))*JackJ[3, 1] + ((t+1)/(6*t^2+5*t+1))*JackJ[4]
        """
        return self( self._P(left) * self._P(right) )

    class Element(JackPolynomials_generic.Element):
        pass



#Q basis
class JackPolynomials_q(JackPolynomials_generic):
    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: Q = JackPolynomialsQ(QQ)
            sage: Q == loads(dumps(Q))
            True
        """
        self._name = "Jack polynomials in the Q basis"
        self._prefix = "JackQ"
        JackPolynomials_generic.__init__(self, R, t=t)

        # Should be shared with _j (and possibly other bases in Macdo/HL) as BasesByRenormalization
        self._P = JackPolynomialsP(R, t)
        # temporary until Hom(GradedHopfAlgebrasWithBasis) works better
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        phi = self._P.module_morphism(diagonal = self._P.scalar_jack_basis, codomain = self, category = category)
        self   .register_coercion(self   ._normalize_morphism(category) *  phi)
        self._P.register_coercion(self._P._normalize_morphism(category) * ~phi)

    def _coerce_start_disabled(self, x):
        """
        EXAMPLES::

            sage: Q = JackPolynomialsQ(QQ)
            sage: P = JackPolynomialsP(QQ)
            sage: Q(sum(P(p) for p in Partitions(3)))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]

        ::

            sage: s = SFASchur(Q.base_ring())
            sage: Q(s([3])) # indirect doctest
            (1/6*t^3-1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3-2*t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]
            sage: Q(s([2,1]))
            (1/3*t^3-1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1]
            sage: Q(s([1,1,1]))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1]
        """

    def _multiply(self, left, right):
        """
        EXAMPLES::

            sage: Q = JackPolynomialsQ(QQ)
            sage: Q([1])^2 # indirect doctest
            JackQ[1, 1] + (2/(t+1))*JackQ[2]
            sage: Q([2])^2
            JackQ[2, 2] + (2/(t+1))*JackQ[3, 1] + ((6*t+6)/(6*t^2+5*t+1))*JackQ[4]
        """
        return self( self._P(left) * self._P(right) )

    class Element(JackPolynomials_generic.Element):
        pass



#############
#   Cache   #
#############
from sage.misc.cache import Cache
cache_p = Cache(JackPolynomials_p)
cache_j = Cache(JackPolynomials_j)
cache_q = Cache(JackPolynomials_q)
cache_z = Cache(ZonalPolynomials)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_j', JackPolynomials_j.Element)
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_p', JackPolynomials_p.Element)
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_q', JackPolynomials_q.Element)
