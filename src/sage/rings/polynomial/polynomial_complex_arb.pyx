# -*- coding: utf-8
r"""
Univariate polynomials over `\CC` with interval coefficients using Arb.

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls using Arb <sage.rings.complex_arb>`

TESTS:

    sage: type(polygen(ComplexBallField(140)))
    <type 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>
    sage: Pol.<x> = CBF[]
    sage: (x+1/2)^3
    x^3 + 1.500000000000000*x^2 + 0.7500000000000000*x + 0.1250000000000000

"""

include "cysignals/signals.pxi"

from sage.rings.integer cimport Integer, smallInteger
from sage.rings.complex_arb cimport ComplexBall

from sage.rings.complex_arb import ComplexBallField

cdef inline long prec(Polynomial_complex_arb pol):
    return pol._parent._base._prec

cdef class Polynomial_complex_arb(Polynomial):
    r"""
    Wrapper for `Arb <http://fredrikj.net/arb/>`_ polynomials of type
    ``acb_poly_t``

    EXAMPLES::

        sage: Pol.<x> = CBF[]
        sage: type(x)
        <type 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>

        sage: Pol(), Pol(1), Pol([0,1,2]), Pol({1: pi, 3: i})
        (0,
         1.000000000000000,
         2.000000000000000*x^2 + x,
         I*x^3 + ([3.141592653589793 +/- 5.61e-16])*x)

        sage: Pol("x - 2/3")
        x + [-0.666666666666667 +/- 4.82e-16]
        sage: Pol(polygen(QQ))
        x

        sage: [Pol.has_coerce_map_from(P) for P in
        ....: QQ['x'], QuadraticField(-1), RealBallField(100)]
        [True, True, True]
        sage: [Pol.has_coerce_map_from(P) for P in
        ....: QQ['y'], RR, CC, RDF, CDF, RIF, CIF, RealBallField(20)]
        [False, False, False, False, False, False, False, False]
    """

    # Memory management and initialization

    def __cinit__(self):
        r"""
        TESTS::

            sage: ComplexBallField(2)['y']()
            0
        """
        acb_poly_init(self.__poly)

    def __dealloc__(self):
        r"""
        TESTS::

            sage: pol = CBF['x']()
            sage: del pol
        """
        acb_poly_clear(self.__poly)

    cdef Polynomial_complex_arb _new(self):
        r"""
        Return a new polynomial with the same parent as this one.
        """
        cdef Polynomial_complex_arb res = Polynomial_complex_arb.__new__(Polynomial_complex_arb)
        res._parent = self._parent
        res._is_gen = 0
        return res

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Initialize this polynomial to the specified value.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
            sage: Pol = CBF['x']
            sage: Polynomial_complex_arb(Pol)
            0
            sage: Polynomial_complex_arb(Pol, is_gen=True)
            x
            sage: Polynomial_complex_arb(Pol, 42, is_gen=True)
            x
            sage: Polynomial_complex_arb(Pol, CBF(1))
            1.000000000000000
            sage: Polynomial_complex_arb(Pol, [])
            0
            sage: Polynomial_complex_arb(Pol, [0])
            0
            sage: Polynomial_complex_arb(Pol, [0, 2, 0])
            2.000000000000000*x
            sage: Polynomial_complex_arb(Pol, (1,))
            1.000000000000000
            sage: Polynomial_complex_arb(Pol, (CBF(i), 1))
            x + I
            sage: Polynomial_complex_arb(Pol, {10: pi})
            ([3.141592653589793 +/- 5.61e-16])*x^10
            sage: Polynomial_complex_arb(Pol, pi)
            [3.141592653589793 +/- 5.61e-16]
        """
        cdef ComplexBall ball
        cdef long length, i

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            acb_poly_set_coeff_si(self.__poly, 1, 1)
        elif x is None:
            acb_poly_zero(self.__poly)
        elif isinstance(x, Polynomial_complex_arb):
            acb_poly_set(self.__poly, (<Polynomial_complex_arb> x).__poly)
        elif isinstance(x, ComplexBall):
            acb_poly_set_coeff_acb(self.__poly, 0, (<ComplexBall> x).value)
        else:
            Coeff = parent.base_ring()
            if isinstance(x, (list, tuple)):
                length = len(x)
                sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                for i in range(length):
                    ball = Coeff(x[i])
                    acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            elif isinstance(x, dict):
                if len(x) == 0:
                    acb_poly_zero(self.__poly)
                else:
                    length = max(int(i) for i in x) + 1
                    sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                    for i, c in x.iteritems():
                        ball = Coeff(x[i])
                        acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            else:
                ball = Coeff(x)
                acb_poly_set_coeff_acb(self.__poly, 0, ball.value)

    # Access

    def degree(self):
        r"""
        Return the (apparent) degree of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2 + 1).degree()
            2
            sage: pol = (x/3 + 1) - x/3; pol
            ([+/- 1.12e-16])*x + 1.000000000000000
            sage: pol.degree()
            1
            sage: Pol([1, 0, 0, 0]).degree()
            0
        """
        return smallInteger(acb_poly_degree(self.__poly))

    cdef get_unsafe(self, Py_ssize_t n):
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._parent._base
        acb_poly_get_coeff_acb(res.value, self.__poly, n)
        return res

    def list(self):
        r"""
        Return the coefficient list of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2/3).list()
            [0, 0, [0.3333333333333333 +/- 7.04e-17]]
            sage: Pol(0).list()
            []
            sage: Pol([0, 1, RBF(0, rad=.1), 0]).list()
            [0, 1.000000000000000, [+/- 0.101]]
        """
        cdef unsigned long length = acb_poly_length(self.__poly)
        return [self.get_unsafe(n) for n in range(length)]
