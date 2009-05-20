"""
Hyperelliptic curves over a padic field.
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hyperelliptic_generic

from sage.rings.all import PowerSeriesRing, PolynomialRing, ZZ, QQ, Integers
from sage.misc.functional import ceil, log
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix

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
        pts = self.lift_x(x, all=True)
        if pts[0][1][0] == P[1]:
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
#        print "dt", dt
        integrals = []
        for f in F:
#            print "f", f
            f_dt = f(x,y) * dt
#            print "f_dt", f_dt
            I = sum([f_dt[n]/(n+1) for n in xrange(f_dt.degree()+1)]) # \int_0^1 f dt
#            print "I", I
            integrals.append(I)
#            integrals.append(f_dt.integral()(1))
        return integrals

    def tiny_integrals_on_basis(self, P, Q):
        """
        Evaluate the integrals of $dx/y$ and $x dx/y$
        by formally integrating a power series in a local parameter $t$

        P and Q MUST be in the same residue disk for this result to make sense.

        TEST:
            sage: K = pAdicField(17, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: TP = E.teichmuller(P);
            sage: E.tiny_integrals_on_basis(P, TP)
            [2*17 + 11*17^2 + 3*17^3 + 16*17^4 + O(17^5), 15*17 + 11*17^2 + 16*17^3 + 11*17^4 + O(17^5)]
        """
        if P == Q:
            V = VectorSpace(self.base_ring(), 2*self.genus())
            return V(0)
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
        K = P[0].parent()
        x = K.teichmuller(P[0])
        pts = self.lift_x(x, all=True)
        p = K.prime()
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
        dim = 2*self.genus()
        V = VectorSpace(K, dim)

        prof("tiny integrals")
        TP = self.teichmuller(P)
#        print "TP", TP
        P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
#        print " P to TP:", P_to_TP[0]

        TQ = self.teichmuller(Q)
#        print "TQ", TQ
        TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
#        print "TQ to  Q:", TQ_to_Q[0]

        prof("mw calc")
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)

        prof("eval f")
        # another hack due to slow padics
        R = forms[0].base_ring()
        try:
            prof("eval f %s"%R)
            L = [f(R(TP[0]), R(TP[1])) - f(R(TQ[0]), R(TQ[1])) for f in forms]
        except ValueError:
            prof("changing rings")
            forms = [f.change_ring(self.base_ring()) for f in forms]
            prof("eval f %s"%self.base_ring())
            L = [f(TP[0], TP[1]) - f(TQ[0], TQ[1]) for f in forms]
        b = 2*V(L)
#        print "b =", b

        prof("lin alg")
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys**(-1) * b

#        print "TP to TQ: ", TP_to_TQ[0]
#        print "\n"
        prof("done")
#        print prof
        return P_to_TP + TP_to_TQ + TQ_to_Q

    coleman_integrals_on_basis_hyperelliptic = coleman_integrals_on_basis


    def invariant_differential(self):
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        return MW.invariant_differential()

    def coleman_integral(self, w, P, Q):
        """
        Example of Leprevost from Kiran Kedlaya
        The first two should be zero as $(P-Q) = 30(P-Q)$ in the Jacobian
        and $dx/y$ and $x dx/y$ are holomorphic.

        sage: K = pAdicField(11, 6)
        sage: x = polygen(K)
        sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
        sage: P = C(-1, 1); P1 = C(-1, -1)
        sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
        sage: x, y = C.monsky_washnitzer_gens()
        sage: w = C.invariant_differential()
        sage: w.coleman_integral(P, Q)
        O(11^6)
        sage: C.coleman_integral(x*w, P, Q)
        O(11^6)
        sage: C.coleman_integral(x^2*w, P, Q)
        3*11 + 2*11^2 + 7*11^3 + 2*11^4 + 10*11^5 + O(11^6)

        sage: p = 71; m = 4
        sage: K = pAdicField(p, m)
        sage: x = polygen(K)
        sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
        sage: P = C(-1, 1); P1 = C(-1, -1)
        sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
        sage: x, y = C.monsky_washnitzer_gens()
        sage: w = C.invariant_differential()
        sage: w.integrate(P, Q), (x*w).integrate(P, Q)
        (O(71^4), O(71^4))
        sage: R, R1 = C.lift_x(4, all=True)
        sage: w.integrate(P, R)
        42*71 + 63*71^2 + 55*71^3 + O(71^4)
        sage: w.integrate(P, R) + w.integrate(P1, R1)
        O(71^4)

    A simple example, integrating dx::

        sage: R.<x> = QQ['x']
        sage: E= HyperellipticCurve(x^3-4*x+4)
        sage: K = Qp(5,10)
        sage: EK = E.change_ring(K)
        sage: P = EK(2, 2)
        sage: Q = EK.teichmuller(P)
        sage: x, y = EK.monsky_washnitzer_gens()
        sage: EK.coleman_integral(x.diff(), P, Q)
        5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
        sage: Q[0] - P[0]
        5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
        """
        # TODO: implement jacobians and show the relationship directly
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, QQ)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        basis_values = self.coleman_integrals_on_basis(P, Q)
        dim = len(basis_values)
        return f(Q[0], Q[1]) - f(P[0], P[1]) + sum([vec[i] * basis_values[i] for i in range(dim)]) # this is just a dot product...



