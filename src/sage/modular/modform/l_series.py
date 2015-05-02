from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.lfunctions.dokchitser import Dokchitser
from l_series_coeffs import gross_zagier_L_series


class GrossZagierLseries(SageObject):

    def __init__(self, E, A, prec=53):
        r"""
        Class for the Gross-Zagier L-series.

        INPUT:

        - ``E`` -- an elliptic curve over `\QQ`

        - ``A`` -- an ideal class in a quadratic number field

        - ``prec`` -- an integer (default 53) giving the required precision

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
        """
        self._E = E
        self._N = N = E.conductor()
        self._A = A
        ideal = A.ideal()
        Q = ideal.quadratic_form().reduced_form()
        a, b, c = self._Q = Q
        D = b ** 2 - 4 * a * c
        self._dokchister = Dokchitser(N ** 2 * D ** 2,
                                      [0, 0, 1, 1],
                                      weight=2, eps=-1, prec=prec)
        self._nterms = nterms = Integer(self._dokchister.gp()('cflength()'))
        if nterms > 1e6:
            # just takes way to long
            raise ValueError("Too many terms: {}".format(nterms))
        zeta_ord = ideal.number_field().zeta_order()
        an_list = gross_zagier_L_series(E.anlist(nterms + 1), Q, N, zeta_ord)
        self._dokchister.gp().set('a', an_list[1:])
        self._dokchister.init_coeffs('a[k]', 1)

    def __call__(self, s, der=0):
        r"""
        Return the value at `s`.

        INPUT:

        - ``s` -- complex number

        - ``der`` -- ? (default 0)

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G(3)
            -0.272946890617590
        """
        return self._dokchister(s, der)

    def taylor_series(self, s, nterms):
        r"""
        Return the Taylor series at `s`.

        INPUT:

        - ``s` -- complex number

        - ``nterms`` -- number of terms in the Taylor series

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: G = GrossZagierLseries(e, A)
            sage: G.taylor_series(2,3)
            -0.613002046122894 + 0.490374999263514*z - 0.122903033710382*z^2 + O(z^3)
        """
        return self._dokchister.taylor_series(s, nterms)

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: e = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-40)
            sage: A = K.class_group().gen(0)
            sage: from sage.modular.modform.l_series import GrossZagierLseries
            sage: GrossZagierLseries(e, A)
            Gross Zagier L-series attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field with ideal class Fractional ideal class (2, 1/2*a)
        """
        msg = "Gross Zagier L-series attached to {} with ideal class {}"
        return msg.format(self._E, self._A)
