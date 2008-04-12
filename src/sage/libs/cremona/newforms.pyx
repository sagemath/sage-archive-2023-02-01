from sage.rings.rational import Rational
from sage.rings.integer cimport Integer
from sage.schemes import elliptic_curves

cdef class ECModularSymbol:
    """
    Modular symbol associated with an elliptic curve,
    using John Cremona's newforms class.

    EXAMPLES:
        sage: from sage.libs.cremona.newforms import ECModularSymbol
        sage: E = EllipticCurve('11a')
        sage: M = ECModularSymbol(E); M
        Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: [M(1/i) for i in range(1,11)]
        [0, 2, 1, -1, -2, -2, -1, 1, 2, 0]

    """
    def __init__(self,E):
        """
        Construct the modular symbol.

        EXAMPLES:
            sage: from sage.libs.cremona.newforms import ECModularSymbol
            sage: E = EllipticCurve('11a')
            sage: M = ECModularSymbol(E)
            sage: E = EllipticCurve('37a')
            sage: M = ECModularSymbol(E)
            sage: E = EllipticCurve('389a')
            sage: M = ECModularSymbol(E)
        """
        cdef ZZ_c a1, a2, a3, a4, a6, N
        cdef Curve *C
        cdef Curvedata *CD
        cdef CurveRed *CR
        cdef int n, t

        _sig_on
        a1 = new_bigint(int(E.a1()))
        a2 = new_bigint(int(E.a2()))
        a3 = new_bigint(int(E.a3()))
        a4 = new_bigint(int(E.a4()))
        a6 = new_bigint(int(E.a6()))
        C = new_Curve(a1,a2,a3,a4,a6)
        CD = new_Curvedata(C[0],0)
        CR = new_CurveRed(CD[0])
        N = getconductor(CR[0])
        n = I2int(N)
        self.n = n

        self.nfs = new_newforms(n,1,0,0)
        self.nfs.createfromcurve(CR[0])
        self._E = E
        _sig_off

    def __repr__(self):
        """
        TESTS:
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
        return "Modular symbol with sign 1 over Rational Field attached to %s"%self._E

    def __call__(self, r):
        """
        Computes the (rational) value of self at a rational number r.

        EXAMPLES:
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
        """
        cdef rational _r
        cdef rational _s
        cdef mpz_t *z_n, *z_d
        cdef ZZ_c *Z_n, *Z_d
        cdef long n, d

        _sig_on
        r = Rational(r)
        d = r.denom()
        n = r.numer() % d
        _r = new_rational(n,d)
        _s = self.nfs.plus_modular_symbol(_r)
        r = Rational((rational_num(_s), rational_den(_s)))
        _sig_off
        return r


#cdef Integer integer_from_ZZ(ZZ_c a):
#    cdef Integer r = Integer(None)
#    ZZ_to_mpz(&r.value, &a)
#    return r
