"""
Hyperelliptic curves over a padic field.
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import hyperelliptic_generic

from sage.rings.all import PowerSeriesRing, PolynomialRing, ZZ, QQ, Integers, Integer, O, pAdicField, mod, LaurentSeriesRing, GF, RR
from sage.misc.functional import ceil, log, sqrt
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix, identity_matrix
import sage.rings.power_series_poly
from sage.modules.all import vector, FreeModule


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
        take the Teichmuller lift of $x$ and then find a matching y
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

    Yet another example::

        sage: R.<x> = QQ['x']
        sage: H = HyperellipticCurve(x*(x-1)*(x+9))
        sage: K = Qp(7,10)
        sage: HK = H.change_ring(K)
        sage: import sage.schemes.elliptic_curves.monsky_washnitzer as mw
        sage: M_frob, forms = mw.matrix_of_frobenius_hyperelliptic(HK)
        sage: w = HK.invariant_differential()
        sage: x,y = HK.monsky_washnitzer_gens()
        sage: f = forms[0]
        sage: S= HK(9,36)
        sage: Q = HK.teichmuller(S)
        sage: P = HK(-1,4)
        sage: b = x*w*w._coeff.parent()(f)
        sage: HK.coleman_integral(b,P,Q)
        7 + 7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^6 + 5*7^7 + 3*7^8 + 4*7^9 + 4*7^10 + O(7^11)

        """
        # TODO: implement Jacobians and show the relationship directly
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        basis_values = self.coleman_integrals_on_basis(P, Q)
        dim = len(basis_values)
        return f(Q[0], Q[1]) - f(P[0], P[1]) + sum([vec[i] * basis_values[i] for i in range(dim)]) # this is just a dot product...


    def frobenius(self, P):
        """
        For a point P = (x,y), this lift of Frobenius maps x to x^p
        and solves for the appropriate y-coordinate that is
        congruent to y modulo p.

        If sqrt(x^p) is not defined over the same field as x, returns an error.

        INPUT:
            - P a point on self (can be over an extension of Q_p)

        OUTPUT:
        the Frobenius of P

        EXAMPLES:
            sage: R.<x> = Qp(11,5)['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = H(5,30)
            sage: H.frobenius(P)
            (5 + 2*11 + 3*11^2 + 2*11^4 + O(11^5) : 8 + 11 + 8*11^2 + 8*11^3 + 8*11^4 + O(11^5) : 1 + O(11^5))
            sage: Q = H(1,6)
            sage: H.frobenius(Q)
            (1 + O(11^5) : 6 + O(11^5) : 1 + O(11^5))
            sage: H.frobenius(Q) == H.teichmuller(Q)
            True
            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: Q = H(0,0)
            sage: u,v = H.local_coord(Q,prec=100)
            sage: K = Qp(11,5)
            sage: L.<a> = K.extension(x^20-11)
            sage: HL = H.change_ring(L)
            sage: S = HL(u(a),v(a))
            sage: HL.frobenius(S)
            (8*a^22 + 10*a^42 + 4*a^44 + 2*a^46 + 9*a^48 + 8*a^50 + a^52 + 7*a^54 +
            7*a^56 + 5*a^58 + 9*a^62 + 5*a^64 + a^66 + 6*a^68 + a^70 + 6*a^74 +
            2*a^76 + 2*a^78 + 4*a^82 + 5*a^84 + 2*a^86 + 7*a^88 + a^90 + 6*a^92 +
            a^96 + 5*a^98 + 2*a^102 + 2*a^106 + 6*a^108 + 8*a^110 + 3*a^112 +
            a^114 + 8*a^116 + 10*a^118 + 3*a^120 + O(a^122) :
            a^11 + 7*a^33 + 7*a^35 + 4*a^37 + 6*a^39 + 9*a^41 + 8*a^43 + 8*a^45 +
            a^47 + 7*a^51 + 4*a^53 + 5*a^55 + a^57 + 7*a^59 + 5*a^61 + 9*a^63 +
            4*a^65 + 10*a^69 + 3*a^71 + 2*a^73 + 9*a^75 + 10*a^77 + 6*a^79 +
            10*a^81 + 7*a^85 + a^87 + 4*a^89 + 8*a^91 + a^93 + 8*a^95 + 2*a^97 +
            7*a^99 + a^101 + 3*a^103 + 6*a^105 + 7*a^107 + 4*a^109 + O(a^111) :
            1 + O(a^100))

        AUTHOR:
            - Jennifer Balakrishnan (2008-02)
        """
        K = self.base_ring()
        p = K.prime()
        x = (P[0])**p
        if x == P[0]:
            return P
        try:
            y = (self.hyperelliptic_polynomials()[0](x)).sqrt()
            if y == 0:
                return self.lift_x(x)
            pts = self.lift_x(x, all=True)
            if (pts[0][1] - P[1]).valuation() > 0:
                return pts[0]
            else:
                return pts[1]
        except (TypeError, NotImplementedError):
            pol = self.hyperelliptic_polynomials()[0]
            yfrob2 = pol(x)
            c = yfrob2.list()[0]
            v = yfrob2.valuation()
            a = yfrob2.parent().gen()
            y = self.newton_sqrt(yfrob2,c.sqrt()*a**(v//2),K.precision_cap())
            if y**2 != yfrob2:
                print "Need more precision for the square root of frobenius!"
            try:
                return self(x,y)
            except ValueError:
                return self._curve_over_ram_extn(x,y)

    def newton_sqrt(self,f,x0, prec):
        """
        NOTE: this function should eventually be moved to p-adic power series ring

        takes the square root of the power series f by Newton's method

        INPUT:
            - f  power series wtih coefficients in Q_p or an extension
            - x0 seeds the Newton iteration
            - prec precision

        OUTPUT:
        the square root of f

        EXAMPLES:
        sage: R.<x> = QQ['x']
        sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
        sage: Q = H(0,0)
        sage: u,v = H.local_coord(Q,prec=100)
        sage: K = Qp(11,5)
        sage: HK = H.change_ring(K)
        sage: L.<a> = K.extension(x^20-11)
        sage: HL = H.change_ring(L)
        sage: S = HL(u(a),v(a))
        sage: f = H.hyperelliptic_polynomials()[0]
        sage: y = HK.newton_sqrt( f(u(a)^11), a^11,5)
        sage: y^2 - f(u(a)^11)
        O(a^122)

        AUTHOR:
            - Jennifer Balakrishnan

        """
        z = x0
        try:
            x = f.parent().variable_name()
            if x!='a' :  #this is to distinguish between extensions of Qp that are finite vs. not
                S = f.base_ring()[[x]]
                x = S.gen()
        except ValueError:
            pass
        z = x0
        loop_prec = (log(RR(prec))/log(RR(2))).ceil()
        for i in range(loop_prec):
            z = (z+f/z)/2
        try:
            return z + O(x**prec)
        except (NameError,ArithmeticError,TypeError):
            return z

    def coleman_integrals_on_basis_no_teichmuller(self, P, Q):
        """
        This is an alternate implementation of Coleman integration
        that does not use Teichm\"{u}ller points.

        INPUT:
            - P  point on self
            - Q  point on self

        OUTPUT:
        the Coleman integrals \int_P^Q w_i

        EXAMPLES:
            We check against the other implementation:
            sage: R.<x> = PolynomialRing(pAdicField(11,5))
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = H(1,6)
            sage: Q = H(-2,12)
            sage: H.coleman_integrals_on_basis(P,Q)
            (9 + 7*11 + 11^2 + 4*11^3 + O(11^4),
             1 + 11 + 2*11^2 + 9*11^3 + O(11^4),
             10*11^-1 + 5 + 7*11 + 3*11^2 + O(11^3),
             8*11^-1 + 9 + 8*11 + 3*11^2 + O(11^3))
            sage: H.coleman_integrals_on_basis_no_teichmuller(P,Q)
            (9 + 7*11 + 11^2 + 4*11^3 + O(11^4),
             1 + 11 + 2*11^2 + 9*11^3 + O(11^4),
             10*11^-1 + 5 + 7*11 + 3*11^2 + O(11^3),
             8*11^-1 + 9 + 8*11 + 3*11^2 + O(11^3))

        AUTHOR:
            - Jennifer Balakrishnan
        """
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        K = self.base_ring()
        p = K.prime()
        dim = 2*self.genus()
        V = VectorSpace(K,dim)
        FP = self.frobenius(P)
        if P == FP:
            P_to_FP = V(dim*[0])
        else:
            P_to_FP = V(self.tiny_integrals_on_basis(P, FP))
        FQ = self.frobenius(Q)
        if Q == FQ:
            FQ_to_Q = V(dim*[0])
        else:
            FQ_to_Q = V(self.tiny_integrals_on_basis(FQ, Q))
        try:
            L = [f(K(P[0]), K(P[1])) - f(K(Q[0]),K(Q[1])) for f in forms]
        except ValueError:
            forms = [f.change_ring(R) for f in forms]
            L = [f(K(P[0]), K(P[1])) - f(K(Q[0]),K(Q[1])) for f in forms]
        b = 2*V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        B = (~M_sys)
        v = [B.list()[i].valuation() for i in range(len(B.list()))]
        vv= min(v)
        B = (p**(-vv)*B).change_ring(K)
        B = p**(vv)*B
        return B*(b-P_to_FP-FQ_to_Q)

    def coleman_integral_no_teichmuller(self,w,P,Q):
        """
        Computes Coleman integrals without using Teichm\"{u}ller points

        INPUT:
            - w differential
            - P point on self
            - Q point on self

        OUTPUT:
        the Coleman integral \int_P^Q w

        EXAMPLES:
        We check against the other implementation of Coleman integrals

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
        sage: C.coleman_integral_no_teichmuller(w,P,Q)
        O(11^6)
        sage: C.coleman_integral(w,P,Q)
        O(11^6)
        sage: C.coleman_integral_no_teichmuller(x*w,P, Q)
        O(11^6)
        sage: C.coleman_integral(x*w,P,Q)
        O(11^6)
        sage: C.coleman_integral_no_teichmuller(x^2*w,P,Q)
        3*11 + 2*11^2 + 7*11^3 + 2*11^4 + 10*11^5 + O(11^6)
        sage: C.coleman_integral(x^2*w,P,Q)
        3*11 + 2*11^2 + 7*11^3 + 2*11^4 + 10*11^5 + O(11^6)

        AUTHOR:
            - Jennifer Balakrishnan
        """
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        g = self.genus()
        const = f(Q[0],Q[1])-f(P[0],P[1])
        if vec == vector(2*g*[0]):
            return const
        else:
            basis_values = self.coleman_integrals_on_basis_no_teichmuller(P, Q)
            dim = len(basis_values)
            dot = sum([vec[i] * basis_values[i] for i in range(dim)])
            return const + dot

    def coleman_integral_from_weierstrass(self,w,Q):
        """
        the Coleman integral \int_P^Q w, where P is a Weierstrass point and w
        is an odd differential

        INPUT:
            - w: differential such that under the hyperelliptic involution i, i*w = -w
            - Q: non-Weierstrass point

        OUTPUT:
        the Coleman integral \int_P^Q w

        EXAMPLE:
            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurve(x^3-4*x+4)
            sage: K = Qp(5,6)
            sage: EK = E.change_ring(K)
            sage: P = EK(2,2)
            sage: Q = EK(1,1)
            sage: w = EK.invariant_differential()
            sage: EK.coleman_integral_from_weierstrass(w,Q)
            4*5 + 4*5^2 + 3*5^4 + 2*5^5 + O(5^6)
            sage: EK.coleman_integral_from_weierstrass(w,P)
            4*5 + 3*5^2 + 4*5^3 + 2*5^4 + O(5^6)
            sage: EK.coleman_integral(w,P,Q)
            5^2 + 5^3 + 2*5^5 + O(5^6)
            sage: EK.coleman_integral_from_weierstrass(w,Q) - EK.coleman_integral_from_weierstrass(w,P)
            5^2 + 5^3 + 2*5^5 + O(5^6)

        AUTHOR:
            - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        x,y = self.local_coordinates_at_infinity(2*prec)
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        if w._coeff(x,-y)*x.derivative()/(-2*y)+w._coeff(x,y)*x.derivative()/(2*y) == 0:
            return self.coleman_integral(w,self(Q[0],-Q[1]), self(Q[0],Q[1]))/2
        else :
            return "sorry, the differential is not odd! we can't do this case"

    def curve_over_ram_extn(self,deg):
        """
        self over Qp(p^(1/deg))

        INPUT:
            - deg: the degree of the ramified extension

        OUTPUT:
        self over Qp(p^(1/deg))

        EXAMPLES:
        sage: R.<x> = QQ['x']
        sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
        sage: K = Qp(11,5)
        sage: HK = H.change_ring(K)
        sage: HL = HK.curve_over_ram_extn(2)
        sage: HL
        Hyperelliptic Curve over Eisenstein Extension of 11-adic Field with capped relative precision 5 in a defined by (1 + O(11^5))*x^2 + (O(11^6))*x + (10*11 + 10*11^2 + 10*11^3 + 10*11^4 + 10*11^5 + O(11^6)) defined by (1 + O(a^10))*y^2 = (1 + O(a^10))*x^5 + (10 + 8*a^2 + 10*a^4 + 10*a^6 + 10*a^8 + O(a^10))*x^3 + (7 + a^2 + O(a^10))*x^2 + (7 + 3*a^2 + O(a^10))*x

        AUTHOR:
            - Jennifer Balakrishnan

        """
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
        K = self.base_ring()
        p = K.prime()
        A = PolynomialRing(QQ,'x')
        x = A.gen()
        J = K.extension(x**deg-p,names='a')
        pol = self.hyperelliptic_polynomials()[0]
        H = HyperellipticCurve(A(pol))
        HJ = H.change_ring(J)
        self._curve_over_ram_extn = HJ
        self._curve_over_ram_extn._curve_over_Qp = self
        return HJ

    def get_boundary_point(self, curve_over_extn, P):
        """
        given self over an extension field, find a point in the disc of P near the boundary

        INPUT:
             - curve_over_extn: self over a totally ramified extension
             - P: Weierstrass point

        OUTPUT:
        a point in the disc of P near the boundary

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(3,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: J.<a> = K.extension(x^30-3)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: S
            (1 + 2*a^2 + 2*a^6 + 2*a^18 + a^32 + a^34 + a^36 + 2*a^38 + 2*a^40 + a^42 + 2*a^44 + a^48 + 2*a^50 + 2*a^52 + a^54 + a^56 + 2*a^60 + 2*a^62 + a^70 + 2*a^72 + a^76 + 2*a^78 + a^82 + a^88 + a^96 + 2*a^98 + 2*a^102 + a^104 + 2*a^106 + a^108 + 2*a^110 + a^112 + 2*a^116 + a^126 + 2*a^130 + 2*a^132 + a^144 + 2*a^148 + 2*a^150 + a^152 + 2*a^154 + a^162 + a^164 + a^166 + a^168 + a^170 + a^176 + a^178 + O(a^180) : a + O(a^181) : 1 + O(a^180))

        AUTHOR:
            - Jennifer Balakrishnan

        """
        J = curve_over_extn.base_ring()
        a = J.gen()
        prec2 = J.precision_cap()
        x,y = self.local_coord(P,prec2)
        return curve_over_extn(x(a),y(a))

    def P_to_S(self, P, S):
        """
        given a finite Weierstrass point P and a point S
        in the same disc, computes the Coleman integrals \int_P^S w_i


        INPUT:
            - P: finite Weierstrass point
            - S: point in disc of P

        OUTPUT:
        Coleman integrals \int_P^S w_i

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: HJ = HK.curve_over_ram_extn(10)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: HK.P_to_S(P, S)
            (4*a + 3*a^3 + 4*a^11 + 4*a^13 + 4*a^17 + 4*a^19 + 2*a^21 + 4*a^23 + 2*a^25 + 4*a^27 + 4*a^29 + a^31 + 4*a^33 + O(a^35), 2*a^-5 + 4*a + 4*a^3 + 2*a^7 + a^11 + 2*a^13 + a^15 + a^17 + 4*a^19 + 4*a^21 + 3*a^23 + 4*a^25 + 2*a^29 + 3*a^31 + 2*a^33 + O(a^35))

        AUTHOR:
            - Jennifer Balakrishnan

        """
        prec = self.base_ring().precision_cap()
        deg = (S[0]).parent().defining_polynomial().degree()
        prec2= prec*deg
        x,y = self.local_coord(P,prec2)
        g = self.genus()
        integrals = [((x**k*x.derivative()/y).integral()) for k in range(2*g)]
        val = [I(S[1]) for I in integrals]
        return vector(val)

    def coleman_integral_P_to_S(self,w,P,S):
        """
        given a finite Weierstrass point P and a point S
        in the same disc, computes the Coleman integral \int_P^S w

        INPUT:
            - w: differential
            - P: Weierstrass point
            - S: point in the same disc of P (S is defined over an extension of \Q_p; coordinates
                 of S are given in terms of uniformizer a)

        OUTPUT:
        Coleman integral \int_P^S w in terms of a

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: J.<a> = K.extension(x^10-5)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: S[0]-P[0] == HK.coleman_integral_P_to_S(x.diff(),P,S)
            True
            sage: HK.coleman_integral_P_to_S(HK.invariant_differential(),P,S) == HK.P_to_S(P,S)[0]
            True

        AUTHOR:
            - Jennifer Balakrishnan

        """
        prec = self.base_ring().precision_cap()
        deg = S[0].parent().defining_polynomial().degree()
        prec2= prec*deg
        x,y = self.local_coord(P,prec2)
        g = self.genus()
        if w.reduce_fast()[1] == vector(2*g*[0]):  #this relates to an earlier normalization
            c = 2
        else:
            c = 1
        int_sing = (w.coeff()(x,y)*x.derivative()/(c*y)).integral()  ##/y or /(2*y)
        int_sing_a = int_sing(S[1])
        return int_sing_a

    def S_to_Q(self,S,Q):
        """
        given S a point on self over an extension field, computes the
        Coleman integrals \int_S^Q w_i for all i

        **one should be able to feed S,Q into coleman_integral,
        but currently that segfaults

        INPUT:
            - S: a point with coordinates in an extension of \Q_p (with unif. a)
            - Q: a non-Weierstrass point defined over \Q_p

        OUTPUT:
        the Coleman integrals \int_S^Q w_i in terms of a

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: w = HK.invariant_differential()
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: P_to_S = HK.P_to_S(P,S)
            sage: S_to_Q = HJ.S_to_Q(S,Q)
            sage: P_to_S + S_to_Q
            (4*a^40 + 2*a^80 + 2*a^100 + O(a^105), 2*a^20 + 4*a^40 + 3*a^60 + a^100 + O(a^105))
            sage: HK.coleman_integral_from_weierstrass(w,Q)
            4*5^2 + 2*5^4 + 2*5^5 + 5^6 + O(5^7)
            sage: HK.coleman_integral_from_weierstrass(x*w,Q)
            2*5 + 4*5^2 + 3*5^3 + 5^5 + 2*5^6 + O(5^7)

        AUTHOR:
            - Jennifer Balakrishnan

        """
        FS = self.frobenius(S)
        FS = (FS[0],FS[1])
        FQ = self.frobenius(Q)
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        try:
            HJ = self._curve_over_ram_extn
            K = HJ.base_ring()
        except AttributeError:
            K = self.base_ring()
        g = self.genus()
        prec2 = K.precision_cap()
        p = K.prime()
        dim = 2*g
        V = VectorSpace(K,dim)
        if S == FS:
            S_to_FS = V(dim*[0])
        else:
            P = self(ZZ(FS[0][0]),ZZ(FS[1][0]))
            x,y = self.local_coord(P,prec2)
            integrals = [(x**i*x.derivative()/y).integral() for i in range(dim)]
            S_to_FS = vector([I(FS[1])-I(S[1]) for I in integrals])
        if Q == FQ:
            FQ_to_Q = V(dim*[0])
        else:
            FQ_to_Q = V(self.tiny_integrals_on_basis(FQ, Q))
        try:
            L = [f(K(S[0]), K(S[1])) - f(K(Q[0]), K(Q[1])) for f in forms]
        except ValueError:
            forms = [f.change_ring(K) for f in forms]
            L = [f(S[0], S[1]) - f(Q[0], Q[1]) for f in forms]
        b = 2*V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        B = (~M_sys)
        v = [B.list()[i].valuation() for i in range(len(B.list()))]
        vv= min(v)
        B = (p**(-vv)*B).change_ring(K)
        B = p**(vv)*B
        return B*(b-S_to_FS-FQ_to_Q)

    def coleman_integral_S_to_Q(self,w,S,Q):
        """
        computes the Coleman integral \int_S^Q w

        **one should be able to feed S,Q into coleman_integral,
        but currently that segfaults

        INPUT:
            - w: a differential
            - S: a point with coordinates in an extension of \Q_p
            - Q: a non-Weierstrass point defined over \Q_p

        OUTPUT:
        the Coleman integrals \int_S^Q w_i

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ,P)
            sage: P_to_S = HK.coleman_integral_P_to_S(y.diff(),P,S)
            sage: S_to_Q = HJ.coleman_integral_S_to_Q(y.diff(),S,Q)
            sage: P_to_S  + S_to_Q
            3 + O(a^120)
            sage: HK.coleman_integral_from_weierstrass(y.diff(),Q)
            3 + O(5^6)

        AUTHOR:
            - Jennifer Balakrishnan

        """
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        R = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(R)
        w = MW(w)
        f, vec = w.reduce_fast()
        g = self.genus()
        const = f(Q[0],Q[1])-f(S[0],S[1])
        if vec == vector(2*g*[0]):
            return const
        else:
            basis_values = self.S_to_Q(S, Q)
            dim = len(basis_values)
            dot = sum([vec[i] * basis_values[i] for i in range(dim)])
            return const + dot

    def coleman_integral_from_weierstrass_via_boundary(self, w, P, Q, d):
        """
        computes the Coleman integral \int_P^Q w via a boundary point
        in the disc of P, defined over a degree d extension

        INPUT:
            - w: a differential
            - P: a Weierstrass point
            - Q: a non-Weierstrass point
            - d: degree of extension where coordinates of boundary point lie

        OUTPUT:
        the Coleman integral \int_P^Q w, written in terms of the uniformizer
        a of the degree d extension

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1,0)
            sage: Q = HK(0,3)
            sage: x,y = HK.monsky_washnitzer_gens()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(y.diff(),P,Q,20)
            3 + O(a^120)
            sage: HK.coleman_integral_from_weierstrass(y.diff(),Q)
            3 + O(5^6)
            sage: w = HK.invariant_differential()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(w,P,Q,20)
            4*a^40 + 2*a^80 + 2*a^100 + O(a^105)
            sage: HK.coleman_integral_from_weierstrass(w,Q)
            4*5^2 + 2*5^4 + 2*5^5 + 5^6 + O(5^7)

        AUTHOR:
            - Jennifer Balakrishnan

        """
        HJ = self.curve_over_ram_extn(deg)
        S = self.get_boundary_point(HJ,P)
        P_to_S = self.coleman_integral_P_to_S(w,P,S)
        S_to_Q = HJ.coleman_integral_S_to_Q(w,S,Q)
        return P_to_S + S_to_Q
