"""
Hyperelliptic curves over a padic field.
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hyperelliptic_generic

from sage.rings.all import PowerSeriesRing, PolynomialRing, ZZ, QQ
from sage.misc.functional import ceil, log

class HyperellipticCurve_padic_field(hyperelliptic_generic.HyperellipticCurve_generic):


# The functions below were prototyped at the 2007 Arizona Winter School by
# Robert Bradshaw and Ralf Gerkmann, working with Miljan Brakovevic and
# Kiran Kedlaya
# All of the below is with respect to the Monsky Washnitzer cohomology.

    def local_analytic_interpolation(self, P, Q):
        """
        Construct a linear interpolation from P to Q in a power series
        in the local parameter t, with precision equal to the p-adic
        precision of the underlying ring.

        Returns a point $X(t) = ( x(t) : y(t) : z(t) )$ such that
        $X(0) = P and X(1) = Q$

        This is implemented by linearly interpolating x, solving formally for y,
        and letting z(t) = 1.

        For this to make sense, P and Q must be in the same residue series, neither equal to infinity.

        TODO: remove the last condition?
        """
        prec = self.base_ring().precision_cap()
        t = PowerSeriesRing(self.base_ring(), 't', prec).gen(0)
        x = P[0]+t*(Q[0]-P[0])
        pts = self.lift_x(x)
        if (pts[0][1] - P[1]).valuation() > 0:
            return pts[0]
        else:
            return pts[1]

    def tiny_integrals(self, F, P, Q):
        """
        Evaluate the integrals of $f_i dx/y$ from P to Q for each f_i in F
        by formally integrating a power series in a local parameter $t$

        P and Q MUST be in the same residue disk for this result to make sense.
        """
        x, y, z = self.local_analytic_interpolation(P, Q)
        dt = x.derivative() / y
        integrals = []
        for f in F:
            f = f(x,y)
            I = (f*dt).integral()
            integrals.append(I(1))
        return integrals

    def tiny_integrals_on_basis(self, P, Q):
        """
        Evaluate the integrals of $dx/y$ and $x dx/y$
        by formally integrating a power series in a local parameter $t$

        P and Q MUST be in the same residue disk for this result to make sense.
        """
        R = PolynomialRing(self.base_ring(), ['x', 'y'])
        x, y = R.gens()
        return self.tiny_integrals([x**i for i in range(2*self.genus())], P, Q)


    def teichmuller(self, P):
        """
        Find a Teichm\:uller point in the same residue class of P.

        Because this lift of frobenius acts as $x \mapsto x^p$,
        take the Teichmuler lift of $x$ and then find a matching y
        from that.

        EXAMPLES:
            sage: K = pAdicField(7, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: E.frobenius(P) == P
            False
            sage: TP = E.teichmuller(P); TP
            (0 : 2 + 3*7 + 3*7^2 + 3*7^4 + O(7^5) : 1 + O(7^5))
            sage: E.frobenius(TP) == TP
            True
            sage: (TP[0] - P[0]).valuation() > 0, (TP[1] - P[1]).valuation() > 0
            (True, True)
        """
        x = padic_teichmuller(P[0])
        pts = self.lift_x(x)
        if (pts[0][1] - P[1]).valuation() > 0:
            return pts[0]
        else:
            return pts[1]

    def coleman_integrals_on_basis(self, P, Q):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        from sage.misc.profiler import Profiler
        prof = Profiler()
        prof("setup")
        K = self.base_ring()
        p = K.prime()
        from sage.modules.free_module import VectorSpace
        from sage.matrix.constructor import matrix
        dim = 2*self.genus()
        V = VectorSpace(K, dim)

        prof("tiny integrals")
        TP = self.teichmuller(P)
#        print "TP", TP
        if P == TP:
            P_to_TP = V(0)
        else:
            P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
#        print " P to TP:", P_to_TP[0]

        TQ = self.teichmuller(Q)
#        print "TQ", TQ
        if Q == TQ:
            TQ_to_Q = V(0)
        else:
            TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
#        print "TQ to  Q:", TQ_to_Q[0]

        prof("mw calc")
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)

        prof("eval f")
        R = forms[0].base_ring()
        L = [f(R(TP[0]), R(TP[1])) - f(R(TQ[0]), R(TQ[1])) for f in forms]
#        L = [f(TP[0], TP[1]) - f(TQ[0], TQ[1]) for f in forms]
        b = 2*V(L)

        prof("lin alg")
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys**(-1) * b

#        print "TP to TQ: ", TP_to_TQ[0]
#        print "\n"
        prof("done")
        print prof
        return P_to_TP + TP_to_TQ + TQ_to_Q

    coleman_integrals_on_basis_hyperelliptic = coleman_integrals_on_basis


    def monsky_washnitzer_gens(self):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        return S.gens()

    def invariant_differential(self):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        return MW.invariant_differential()

    def coleman_integral(self, w, P, Q):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, QQ)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        print "w", w
        w = MW(w)
        f, rem = w.reduce_fast()
        vec = rem.H1_vector()
        basis_values = self.coleman_integrals_on_basis(P, Q)
        dim = len(basis_values)
        return f(P[0], P[1]) - f(Q[0], Q[1]) + sum([vec[i] * basis_values[i] for i in range(dim)]) # this is just a dot product...



# TODO: add this to padics (if it isn't there in the new version already).
def padic_teichmuller(a):
    K = a.parent()
    p = K.prime()
    p_less_1_inverse = ~K(p-1)
    one = K(1)
    x = K(ZZ(a) % p)
    if x == 0:
        return x
    for _ in range(ceil(log(K.precision_cap(),2))):
        x = ((p-2)*x + x**(2-p)) * p_less_1_inverse
    return x





