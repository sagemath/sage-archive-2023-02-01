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

from sage.libs.symmetrica.all import hall_littlewood
import sfa
import sage.combinat.partition
import kfpoly
from sage.matrix.all import matrix
from sage.rings.all import ZZ
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom


##################################
#Still under major development!!!#
##################################

def HallLittlewoodP(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood `P`
    basis. This is the same as the `HL` basis in John Stembridge's SF
    examples file.

    If `t` is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over `R` with the variable `t`
    and taking its fraction field.

    EXAMPLES::

        sage: HallLittlewoodP(QQ)
        Hall-Littlewood polynomials in the P basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodP(QQ, t=-1)
        Hall-Littlewood polynomials in the P basis with t=-1 over Rational Field
        sage: HLP = HallLittlewoodP(QQ)
        sage: s = SFASchur(HLP.base_ring())
        sage: s(HLP([2,1]))
        (-t^2-t)*s[1, 1, 1] + s[2, 1]

    The Hall-Littlewood polynomials in the `P` basis at `t = 0` are the
    Schur functions.

    ::

        sage: HLP = HallLittlewoodP(QQ,t=0)
        sage: s = SFASchur(HLP.base_ring())
        sage: s(HLP([2,1])) == s([2,1])
        True

    The Hall-Littlewood polynomials in the `P` basis at `t = 1` are the
    monomial symmetric functions.

    ::

        sage: HLP = HallLittlewoodP(QQ,t=1)
        sage: m = SFAMonomial(HLP.base_ring())
        sage: m(HLP([2,2,1])) == m([2,2,1])
        True


    We end with some examples of coercions between:

        1. Hall-Littlewood `P` basis.

        2. Hall-Littlewood polynomials in the `Q` basis

        3.  Hall-Littlewood polynomials in the `Q^\prime` basis (via the Schurs)

        4. Classical symmetric functions

    EXAMPLES::

        sage: HLP  = HallLittlewoodP(QQ)
        sage: HLQ  = HallLittlewoodQ(QQ)
        sage: HLQp = HallLittlewoodQp(QQ)
        sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
        sage: HLP(HLQ([2])) # indirect doctest
        (-t+1)*P[2]
        sage: HLP(HLQp([2]))
        t*P[1, 1] + P[2]
        sage: HLP(s([2]))
        t*P[1, 1] + P[2]
        sage: HLP(p([2]))
        (t-1)*P[1, 1] + P[2]

    TESTS::

        sage: HLP(s[[]])
        P[]
        sage: HLQ(s[[]])
        Q[]
        sage: HLQp(s[[]])
        Qp[]
    """
    return cache_p(R,t)

def HallLittlewoodQ(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood `Q`
    basis. This is the same as the `Q` basis in John Stembridge's SF
    examples file.

    If `t` is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over `R` with the variable `t`
    and taking its fraction field.

    EXAMPLES::

        sage: HallLittlewoodQ(QQ)
        Hall-Littlewood polynomials in the Q basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: HallLittlewoodQ(QQ, t=-1)
        Hall-Littlewood polynomials in the Q basis with t=-1 over Rational Field
    """
    return cache_q(R,t)

def HallLittlewoodQp(R,t=None):
    """
    Returns the algebra of symmetric functions in Hall-Littlewood `Q^\prime` (Qp)
    basis. This is dual to the Hall-Littlewood `P` basis with respect to
    the standard scalar product.

    If `t` is not specified, then the base ring will be obtained by
    making the univariate polynomial ring over `R` with the variable `t`
    and taking its fraction field.

    EXAMPLES::

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
        TESTS::

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

        sfa.SymmetricFunctionAlgebra_generic.__init__(self, R)

#        print self, self._s
        # This coercion is broken: HLP = HallLittlewoodP(QQ); HLP(HLP._s[1])

        # Bases defined by orthotriangularity should inherit from some
        # common category BasesByOrthotriangularity (shared with Jack, HL, orthotriang, Mcdo)
        if hasattr(self, "_s_cache"):
            self._s = sfa.SFASchur(R)
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = sage.categories.all.ModulesWithBasis(self.base_ring())
            self   .register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
            self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))

    def _s_to_self(self, x):
        """
        Isomorphism from the Schur basis into self

        EXAMPLES::

            sage: P = HallLittlewoodP(QQ,t=2)
            sage: s = SFASchur(P.base_ring())
            sage: P._s_to_self(s[2,1])
            6*P[1, 1, 1] + P[2, 1]

        This is for internal use only. Please use instead::

            sage: P(s[2,1])
            6*P[1, 1, 1] + P[2, 1]
        """
        return self._from_cache(x, self._s_cache, self._s_to_self_cache, t = self.t)

    def _self_to_s(self, x):
        r"""
        Isomorphism from self to the Schur basis

        EXAMPLES::

            sage: P = HallLittlewoodP(QQ,t=2)
            sage: s = SFASchur(P.base_ring())
            sage: P._self_to_s(P[2,1])
            -6*s[1, 1, 1] + s[2, 1]

        This is for internal use only. Please use instead::

            sage: s(P[2,1])
            -6*s[1, 1, 1] + s[2, 1]
        """
        return self._s._from_cache(x, self._s_cache, self._self_to_s_cache, t = self.t) # do we want this t = self.t?

    def transition_matrix(self, basis, n):
        """
        Returns the transitions matrix between self and basis for the
        homogenous component of degree n.

        EXAMPLES::

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

# TODO: move this as ElementHallLittlewood_generic.Element
class ElementHallLittlewood_generic_Element(sfa.SymmetricFunctionAlgebra_generic.Element):
    def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in `n`
        variables.

        EXAMPLES::

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

        This is the default implementation that converts both self and x
        into Schur functions and performs the scalar product that basis.

        EXAMPLES::

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

        EXAMPLES::

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

HallLittlewood_generic.Element = ElementHallLittlewood_generic_Element


###########
# P basis #
###########
p_to_s_cache = {}
s_to_p_cache = {}

class HallLittlewood_p(HallLittlewood_generic):

    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: P = HallLittlewoodP(QQ)
            sage: TestSuite(P).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
        """
        self._name = "Hall-Littlewood polynomials in the P basis"
        self._prefix = "P"

        self._self_to_s_cache = p_to_s_cache
        self._s_to_self_cache = s_to_p_cache

        HallLittlewood_generic.__init__(self, R, t=t)



    # We probably want to get rid of this one
    def _multiply(self, left, right):
        """
        Convert to the Schur basis, do the multiplication there, and
        convert back to the `P` basis.

        EXAMPLES::

            sage: HLP = HallLittlewoodP(QQ)
            sage: HLP([2])^2 # indirect doctest
            (t+1)*P[2, 2] + (-t+1)*P[3, 1] + P[4]
        """
        return self( self._s(left) * self._s(right) )

    def _q_to_self(self, m):
        """
        Returns the scalar coefficient on self(m) when converting from the
        `Q` basis to the `P` basis. Note that this assumes that m is a
        Partition object.

        Todo: find a better name

        EXAMPLES::

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

    def _s_to_self_base(self, part):
        """
        Returns a function which gives the coefficient of a partition
        in the expansion of the Schur functions ``s(part)`` in self.

        EXAMPLES::

            sage: P = HallLittlewoodP(QQ)
            sage: f21 = P._s_to_self_base(Partition([2,1]))
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
        Computes the change of basis between the `P` polynomials and the
        Schur functions for partitions of size `n`.

        Uses the fact that the transformation matrix is upper-triangular in
        order to obtain the inverse transformation.

        EXAMPLES::

            sage: P = HallLittlewoodP(QQ)
            sage: P._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(P._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], t), ([2], 1)])]
            sage: l(P._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -t), ([2], 1)])]
        """

        self._invert_morphism(n, ZZ['t'], self._self_to_s_cache, \
                              self._s_to_self_cache, to_self_function = self._s_to_self_base, \
                              upper_triangular=True, ones_on_diagonal=True)





###########
# Q basis #
###########

class HallLittlewood_q(HallLittlewood_generic):
    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: Q = HallLittlewoodQ(QQ)
            sage: Q == loads(dumps(Q))
            True
        """
        self._name = "Hall-Littlewood polynomials in the Q basis"
        self._prefix = "Q"

        HallLittlewood_generic.__init__(self, R, t=t)

        self._P = HallLittlewood_p(R, t=t)
        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = sage.categories.all.ModulesWithBasis(self.base_ring())

        phi = self.module_morphism(diagonal = self._P._q_to_self, codomain = self._P, category = category)
        self._P.register_coercion(phi)
        self   .register_coercion(~phi)

    # TODO: discard (except doctests)
    def _multiply(self, left, right):
        """
        Converts to the `P` basis, does the multiplication there, and
        converts back to the `Q` basis.

        EXAMPLES::

            sage: HLQ = HallLittlewoodQ(QQ)
            sage: HLQ([2])^2 # indirect doctest
            Q[2, 2] + (-t+1)*Q[3, 1] + (-t+1)*Q[4]
        """
        return self( self._P(left) * self._P(right) )


    # TODO: discard (except doctests)
    def _coerce_start_disabled(self, x):
        """
        EXAMPLES::

            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQ( HLP([2,1]) + HLP([3]) )
            (1/(t^2-2*t+1))*Q[2, 1] + (1/(-t+1))*Q[3]
            sage: HLQ(HLQp([2])) # indirect doctest
            (t/(t^3-t^2-t+1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(s([2]))
            (t/(t^3-t^2-t+1))*Q[1, 1] + (1/(-t+1))*Q[2]
            sage: HLQ(p([2]))
            (1/(t^2-1))*Q[1, 1] + (1/(-t+1))*Q[2]
        """
        pass

    def _p_to_self(self, m):
        """
        Returns the scalar coefficient on self(m) when converting from the
        `Q` basis to the `P` basis. Note that this assumes that m is a
        Partition object.

        Note: this is not used anymore!

        EXAMPLES::
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



class HallLittlewood_qp(HallLittlewood_generic):

    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, R, t=None):
        """
        EXAMPLES::

            sage: Qp = HallLittlewoodQp(QQ)
            sage: Qp == loads(dumps(Qp))
            True
        """
        self._name = "Hall-Littlewood polynomials in the Qp basis"
        self._prefix = "Qp"

        self._self_to_s_cache = qp_to_s_cache
        self._s_to_self_cache = s_to_qp_cache

        HallLittlewood_generic.__init__(self,R, t=t)


    # TODO: discard (except doctests)
    def _coerce_start_disabled(self, x):
        """
        Coerce things into the Hall-Littlewood `Q^\prime` basis.

        1. Hall-Littlewood polynomials in the `Q` basis (via the Schurs)

        2. Hall-Littlewood polynomials in the `P` basis (via the Schurs)

        3. Classical symmetric functions

        EXAMPLES::

            sage: HLP  = HallLittlewoodP(QQ)
            sage: HLQ  = HallLittlewoodQ(QQ)
            sage: HLQp = HallLittlewoodQp(QQ)
            sage: s = SFASchur(HLP.base_ring()); p = SFAPower(HLP.base_ring())
            sage: HLQp(HLP([2])) # indirect doctest
            -t*Qp[1, 1] + (t^2+1)*Qp[2]
            sage: HLQp(HLQ([2]))
            (t^2-t)*Qp[1, 1] + (-t^3+t^2-t+1)*Qp[2]
            sage: HLQp(s([2]))
            Qp[2]
            sage: HLQp(p([2]))
            -Qp[1, 1] + (t+1)*Qp[2]
        """
        if isinstance(x, HallLittlewood_q.Element):
            return self( self._s( x ) )
        elif isinstance(x, HallLittlewood_p.Element):
            return self( self._s( x ) )
        elif isinstance(x, sfa.SymmetricFunctionAlgebra_generic.Element):
            sx = self._s( x )
            return self._from_cache(sx, self._s_cache, self._s_to_self_cache,t=self.t)
        else:
            raise TypeError

    # TODO: discard (except doctests)
    def _multiply(self, left, right):
        """
        Converts the Hall-Littlewood polynomial in the `Q^\prime` basis to a Schur
        function, performs the multiplication there, and converts it back
        to the `Q^\prime` basis.

        EXAMPLES::

            sage: HLQp = HallLittlewoodQp(QQ)
            sage: HLQp([2])^2 # indirect doctest
            Qp[2, 2] + (-t+1)*Qp[3, 1] + (-t+1)*Qp[4]
        """
        return self( self._s(left)*self._s(right) )

    def _to_s(self, part):
        """
        Returns a function which gives the coefficient of a partition
        in the Schur expansion of ``self(part)``.

        EXAMPLES::

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
        Computes the change of basis between the `Q^\prime` polynomials and the
        Schur functions for partitions of size `n`.

        Uses the fact that the transformation matrix is lower-triangular in
        order to obtain the inverse transformation.

        EXAMPLES::

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

# Unpickling backward compatibility
sage.structure.sage_object.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_p', HallLittlewood_p.Element)
sage.structure.sage_object.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_q', HallLittlewood_q.Element)
sage.structure.sage_object.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_qp', HallLittlewood_qp.Element)
