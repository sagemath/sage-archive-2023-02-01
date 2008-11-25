from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.lfunctions.dokchitser import Dokchitser
from l_series_coeffs import gross_zagier_L_series


class GrossZagierLseries(SageObject):

    def __init__(self, E, A, prec=53):
        """
        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
        """
        Q = A.ideal().quadratic_form().reduced_form()
        self._A = A
        self._Q = Q
        self._E = E
        self._N = N = E.conductor()
        a, b, c = self._Q = Q
        D = b ** 2 - 4 * a * c
        self._dokchister = Dokchitser(N ** 2 * D ** 2 / Integer(4),
                                      [0, 0, 1, 1],
                                      weight=2, eps=-1, prec=prec)
        self._nterms = nterms = Integer(self._dokchister.gp()('cflength()'))
        if nterms > 1e6:
            # just takes way to long
            raise ValueError("Too many terms: %s" % nterms)
        an_list = list(gross_zagier_L_series(E.anlist(nterms + 1), Q, N, A.ideal().number_field().zeta_order()))
        self._dokchister.gp().set('a', an_list[1:])
        self._dokchister.init_coeffs('a[k]', 1)

    def __call__(self, s, der=0):
        r"""
        Return the value at `s`

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G(3)
            -0.272819809220447
        """
        return self._dokchister(s, der)

    def taylor_series(self, s, nterms):
        r"""
        Return the Taylor series at `s`

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G.taylor_series(2,3)
            -0.598849375222341 + 0.431769334217468*z - 0.0128395302793324*z^2 + O(z^3)
        """
        return self._dokchister.taylor_series(s, nterms)

    def _repr_(self):
        """
        Return the string representation

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: GrossZagierLseries(e, A)
            Gross Zagier L-series attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field with ideal class Fractional ideal class (2, 1/2*a)
        """
        return "Gross Zagier L-series attached to %s with ideal class %s" % (self._E, self._A)
