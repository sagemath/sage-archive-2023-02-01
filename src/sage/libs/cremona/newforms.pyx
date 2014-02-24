from sage.rings.rational import Rational
from sage.rings.integer cimport Integer
from sage.schemes import elliptic_curves
from sage.modular.all import Cusp

#*****************************************************************************
#       Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef class ECModularSymbol:
    """
    Modular symbol associated with an elliptic curve,
    using John Cremona's newforms class.

    EXAMPLES::

        sage: from sage.libs.cremona.newforms import ECModularSymbol
        sage: E = EllipticCurve('11a')
        sage: M = ECModularSymbol(E,1); M
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i) for i in range(1,11)]
        [0, 2, 1, -1, -2, -2, -1, 1, 2, 0]

    """
    def __init__(self, E, sign=1):
        """
        Construct the modular symbol.

        EXAMPLES::

            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E)
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E)

        TESTS:

        This one is from :trac: `8042` (note that Cremona's code is more limited
        when run on a 32-bit platform, but works fine on 64-bit in this
        case; see :trac: `8114`).

        ::
            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('858k2')
            sage: ECModularSymbol(E)
            Traceback (most recent call last):                           # 32-bit
            ...                                                          # 32-bit
            OverflowError: Python int too large to convert to C long     # 32-bit
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + x*y = x^3 + 16353089*x - 335543012233 over Rational Field         # 64-bit
        """
        cdef ZZ_c a1, a2, a3, a4, a6, N
        cdef Curve *C
        cdef Curvedata *CD
        cdef CurveRed *CR
        cdef int n, t

        a1 = new_bigint(int(E.a1()))
        a2 = new_bigint(int(E.a2()))
        a3 = new_bigint(int(E.a3()))
        a4 = new_bigint(int(E.a4()))
        a6 = new_bigint(int(E.a6()))

        sig_on()
        C = new_Curve(a1,a2,a3,a4,a6)
        CD = new_Curvedata(C[0],0)
        CR = new_CurveRed(CD[0])
        N = getconductor(CR[0])
        n = I2int(N)
        self.n = n
        self.sign = sign

        self.nfs = new_newforms(n,0)
        self.nfs.createfromcurve(sign,CR[0])
        self._E = E
        sig_off()

    def __repr__(self):
        """
        TESTS::

            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return "Modular symbol with sign %s over Rational Field attached to %s"%(self.sign, self._E)

    def __call__(self, r):
        """
        Computes the (rational) value of self at a rational number r.

        EXAMPLES::

            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,10)]
            [0, 2, 1, -1, -2, -2, -1, 1, 2]
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,10)]
            [0, 0, 0, 0, 1, 0, 1, 1, 0]
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E)
            sage: [M(1/i) for i in range(1,10)]
            [0, 0, 0, 0, 4, 0, 2, 0, -2]

        TESTS (see :trac: `11211`)::

            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: M(oo)
            2/5
            sage: M(7/5)
            3
            sage: M("garbage")
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert garbage to a Cusp
            sage: M(7/5)
            3
        """
        cdef rational _r
        cdef rational _s
        cdef mpz_t *z_n, *z_d
        cdef ZZ_c *Z_n, *Z_d
        cdef long n, d

        r = Cusp(r)
        d = r.denominator()
        n = r.numerator()
        if d != 0:
            n = n % d
        sig_on()
        _r = new_rational(n,d)
        _s = self.nfs.plus_modular_symbol(_r)
        sig_off()
        return Rational((rational_num(_s), rational_den(_s)))
