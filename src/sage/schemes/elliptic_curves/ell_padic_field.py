"""
Elliptic curves over padic fields
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein   <wstein@gmail.com>
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
from sage.rings.all import PowerSeriesRing, PolynomialRing, IntegerModRing, ZZ, QQ
from sage.misc.functional import ceil, log
from sage.libs.pari.all import pari

# Elliptic curves are very different than genus > 1 hyperelliptic curves,
# there is an "is a" relationship here, and common implementation with regard
# Coleman integration.
from sage.schemes.hyperelliptic_curves.hyperelliptic_padic_field import HyperellipticCurve_padic_field

import sage.databases.cremona


class EllipticCurve_padic_field(EllipticCurve_field, HyperellipticCurve_padic_field):
    """
    Elliptic curve over a padic field.
    """
    def __init__(self, x, y=None):
        """
        Constructor from [a1,a2,a3,a4,a6] or [a4,a6].
        EXAMPLES:
        sage: Qp=pAdicField(17)
        sage: E=EllipticCurve(Qp,[2,3]); E
        Elliptic Curve defined by y^2  = x^3 + (2+O(17^20))*x + (3+O(17^20)) over 17-adic Field with capped relative precision 20
        sage: E == loads(dumps(E))
        True
        """
        if y is None:
            if isinstance(x, list):
                ainvs = x
                field = ainvs[0].parent()
        else:
            if isinstance(y, str):
                field = x
                X = sage.databases.cremona.CremonaDatabase()[y]
                ainvs = [field(a) for a in X.a_invariants()]
            else:
                field = x
                ainvs = y
        if not (isinstance(field, ring.Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(self, [field(x) for x in ainvs])

        self._point_class = ell_point.EllipticCurvePoint_field
        self._genus = 1

    def _pari_(self):
        """
        Convert the elliptic curve to Pari. It requires the p-adic j-valuation
        to be negative, i.e., the j-valuation must not be a p-adic integer.
        EXAMPLES:
        sage: Qp=pAdicField(5, prec=3)
        sage: E=EllipticCurve(Qp,[3, 4])
        sage: E._pari_()
        [O(5^3), O(5^3), O(5^3), 3 + O(5^3), 4 + O(5^3), O(5^3), 1 + 5 + O(5^3), 1 + 3*5 + O(5^3), 1 + 3*5 + 4*5^2 + O(5^3), 1 + 5 + 4*5^2 + O(5^3), 4 + 3*5 + 5^2 + O(5^3), 2*5 + 4*5^2 + O(5^3), 3*5^-1 + O(5), [4 + 4*5 + 4*5^2 + O(5^3)], 1 + 2*5 + 4*5^2 + O(5^3), 1 + 5 + 4*5^2 + O(5^3), 2*5 + 4*5^2 + O(5^3), 3 + 3*5 + 3*5^2 + O(5^3), 0]
        sage: E.j_invariant()
        3*5^-1 + O(5)
        sage: E=EllipticCurve(Qp,[1, 1])
        sage: E.j_invariant() # the j-invariant is a p-adic integer
        2 + 4*5^2 + O(5^3)
        sage: E._pari_()
        Traceback (most recent call last):
        ...
        PariError:  (8)
        """
        try:
            return self.__pari
        except AttributeError:
            pass
        F = self.base_ring()
        self.__pari = pari('ellinit([%s,%s,%s,%s,%s])'%tuple([b._pari_() for b in self.ainvs()]))
        return self.__pari


    def frobenius(self, P=None):
        """
        Returns the Frobenius as a function on the group of points of
        this elliptic curve.

        EXAMPLE:
            sage: Qp=pAdicField(13)
            sage: E=EllipticCurve(Qp,[1,1])
            sage: type(E.frobenius())
            <type 'function'>
            sage: point=E(0,1)
            sage: E.frobenius(point)
            (0 : 1 + O(13^20) : 1 + O(13^20))
        """
        try:
            _frob = self._frob
        except AttributeError:
            K = self.base_field()
            p = K.prime()
            x = PolynomialRing(K, 'x').gen(0)

            a1, a2, a3, a4, a6 = self.a_invariants()
            if a1 != 0 or a2 != 0:
                raise NotImplementedError, "Curve must be in weierstrass normal form."

            f = x*x*x + a2*x*x + a4*x + a6
            h = (f(x**p) - f**p)

            # internal function: I don't know how to doctest it...
            def _frob(P):
                x0 = P[0]
                y0 = P[1]
                uN = (1 + h(x0)/y0**(2*p)).sqrt()
                yres=y0**p * uN
                xres=x0**p
                if (yres-y0).valuation() == 0:
                    yres=-yres
                return self.point([xres,yres, K(1)])

            self._frob = _frob

        if P is None:
            return _frob
        else:
            return _frob(P)

    def coleman_integrals_on_basis(self, P, Q):
        """
        Return the coleman integral of dx/y and x dx/y from P to Q.

        EXAMPLES:
            sage: K = pAdicField(13, 7)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: res = E.coleman_integrals_on_basis(P, 2*P); res
            (O(13^7), 2 + 7*13 + 2*13^2 + 5*13^3 + 10*13^4 + 7*13^5 + 2*13^6 + O(13^7))

        As the Coleman integral of dx/y is in invariant under
        translation, it should evaluate to zero between a torsion
        point and its multiples.

            sage: w = E.invariant_differential()
            sage: w.coleman_integral(P, 2*P)
            O(13^7)
        """
        from sage.misc.profiler import Profiler
        prof = Profiler()
        prof("setup")
        K = self.base_field()
        p = K.prime()
        from sage.modules.free_module import VectorSpace
        from sage.matrix.constructor import matrix
        V = VectorSpace(K, 2)

        prof("tiny integrals")
        TP = self.teichmuller(P)
#        print "TP", TP
        P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
#        print " P to TP:", P_to_TP[0]

        TQ = self.teichmuller(Q)
#        print "TQ", TQ
        TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
#        print "TQ to  Q:", TQ_to_Q[0]

        prof("mw setup")
        import monsky_washnitzer
        # TODO fis matrix_of_frobenius code to use real padics
        prec = K.precision_cap()
        extra_prec = monsky_washnitzer.adjusted_prec(p, prec)
        pseudo_Qp = IntegerModRing(p**extra_prec)
        x = PolynomialRing(pseudo_Qp,'x').gen(0)
        q = x*x*x + pseudo_Qp(self.a2())*x*x + pseudo_Qp(self.a4())*x + pseudo_Qp(self.a6())

        prof("mw calc")
        M_frob, f0, f1 = monsky_washnitzer.matrix_of_frobenius(q, p, extra_prec, None, True)
        prof("eval")
        f0 *= 2
        f1 *= 2

        TPx = pseudo_Qp(TP[0])
        TPy = pseudo_Qp(TP[1])
        TQx = pseudo_Qp(TQ[0])
        TQy = pseudo_Qp(TQ[1])

        prof("linalg")
        L = [K(f0(TPx)(TPy) - f0(TQx)(TQy)), K(f1(TPx)(TPy) - f1(TQx)(TQy))]
        from sage.rings.all import ZZ
        L = [ZZ(t) for t in L] # pass through ZZ due to bug in p-adics
        b = V(L)
#        print "b =", b
        M_sys = matrix(K, M_frob).transpose() - 1
#        print M_sys
        TP_to_TQ = M_sys**(-1) * b

#        print "TP to TQ: ", TP_to_TQ[0]
#        print "\n"
        prof("done")
#        print prof
        return P_to_TP + TP_to_TQ + TQ_to_Q


