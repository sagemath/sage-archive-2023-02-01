"""
Elliptic curves over padic fields
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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


import sage.rings.ring as ring
from ell_field import EllipticCurve_field
import ell_point
from sage.rings.all import PowerSeriesRing, PolynomialRing, IntegerModRing

import sage.databases.cremona


class EllipticCurve_padic_field(EllipticCurve_field):
    """
    Elliptic curve over a padic field.
    """
    def __init__(self, x, y=None):
        if y is None:
            if isinstance(x, list):
                ainvs = x
                field = ainvs[0].parent()
        else:
            if isinstance(y, str):
                field = x
                X = sage.databases.cremona.CremonaDatabase()[label]
                ainvs = [field(a) for a in X.a_invariants()]
            else:
                field = x
                ainvs = y
        if not (isinstance(field, ring.Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(self, [field(x) for x in ainvs])

        self._point_class = ell_point.EllipticCurvePoint_field

    def _pari_(self):
        try:
            return self.__pari
        except AttributeError:
            pass
        F = self.base_ring()
        self.__pari = pari('ellinit(%s,%s,%s,%s,%s)'%(F.characteristic(), [b._pari_() for b in self.ainvs()]))
        return self.__pari


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
        Kt = PowerSeriesRing(self.base_ring(), 't', prec)
        t = Kt.gen(0)
        x = P[0]+t*(Q[0]-P[0])
        a1, a2, a3, a4, a6 = self.a_invariants()
        if a1 != 0 or a2 != 0:
            raise NotImplementedError, "Curve must be in weierstrass normal form."
        y = (x*x*x + a2*x*x + a4*x + a6).sqrt()
        if (y(0) - P[1]).valuation() == 0:
            y = -y
        return self.point([x, y, Kt(0)], check=False)

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
        return self.tiny_integrals([R(1), R.gen(0)], P, Q)

    def p_power_frobenious(self):
        K = self.base_field()
        p = K.prime()
        x = PolynomialRing(K, 'x').gen(0)

        a1, a2, a3, a4, a6 = self.a_invariants()
        if a1 != 0 or a2 != 0:
            raise NotImplementedError, "Curve must be in weierstrass normal form."

        f = x*x*x + a2*x*x + a4*x + a6
        h = (f(x**p) - f**p)

        def frob(P):
            x0 = P[0]
            y0 = P[1]
            uN = (1 + h(x0)/y0**(2*p)).sqrt()
            yres=y0**p * uN
            xres=x0**p
            if (yres-y0).valuation() == 0:
                yres=-yres
            return self.point((xres,yres))

        return frob

    def teichmuller(self, P, frob=None):
        """
        Find a Teichm\:uller point in the same residue class of P.

        TODO: what kind of convergence am I guerenteed here?
        """
        if frob is None:
            frob = self.p_power_frobenious()
        for i in xrange(self.base_field().precision_cap()):
            P = frob(P)
        return P


    def coleman_integrals_on_basis(self, P, Q):
        """
        Return the coleman integral of dx/y and x dx/y from P to Q

        sage: K = pAdicField(p, prec)
        sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
        sage: P = E(K(14/3), K(11/2))
        sage: res = E.coleman_integrals_on_basis(); res
        (7*13^6 + O(13^7), 2 + 7*13 + 2*13^2 + 5*13^3 + 10*13^4 + 7*13^5 + 8*13^6 + O(13^7))

        As the Coleman integral of dx/y is in invariant under translation, it should
        evaluate to zero between a torsion point and its multiples.

        sage: res[0].valuation() >= 6
        True
        """
        K = self.base_field()
        p = K.prime()
        from sage.modules.free_module import VectorSpace
        from sage.matrix.constructor import matrix
        V = VectorSpace(K, 2)

        frob = self.p_power_frobenious()
        TP = self.teichmuller(P, frob)
#        print "TP", TP
        P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
#        print " P to TP:", P_to_TP[0]

        TQ = self.teichmuller(Q, frob)
#        print "TQ", TQ
        TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
#        print "TQ to  Q:", TQ_to_Q[0]

        # TODO fis matrix_of_frobenius code to use real padics
        prec = K.precision_cap()
        pseudo_Qp = IntegerModRing(p**prec)
        x = PolynomialRing(pseudo_Qp,'x').gen(0)
        q = x*x*x + pseudo_Qp(self.a2())*x*x + pseudo_Qp(self.a4())*x + pseudo_Qp(self.a6())

        import monsky_washnitzer
        M_frob, f0, f1 = monsky_washnitzer.matrix_of_frobenius(q, p, prec, None, True)
        f0 *= 2
        f1 *= 2

        TPx = pseudo_Qp(TP[0])
        TPy = pseudo_Qp(TP[1])
        TQx = pseudo_Qp(TQ[0])
        TQy = pseudo_Qp(TQ[1])

        L = [K(f0(TPx)(TPy) - f0(TQx)(TQy)), K(f1(TPx)(TPy) - f1(TQx)(TQy))]
        from sage.rings.all import ZZ
        L = [ZZ(t) for t in L] # pass through ZZ due to bug in p-adics
        b = V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys**(-1) * b

#        print "TP to TQ: ", TP_to_TQ[0]
#        print "\n"
        return P_to_TP + TP_to_TQ + TQ_to_Q
