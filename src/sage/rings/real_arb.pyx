r"""
Arbitrary precision real intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-21): Initial version.

This is a rudimentary binding to the optional `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

You may have to run ``sage -i arb`` to use the arb library.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

Comparison
==========

.. WARNING::

    Identical :class:`RealBall` objects are understood to give
    permission for algebraic simplification. This assumption is made
    to improve performance.  For example, setting ``z = x*x`` sets `z`
    to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: from sage.rings.real_arb import RealBallField # optional - arb
    sage: RBF = RealBallField() # optional - arb
    sage: a = RBF(1) # optional - arb
    sage: b = RBF(1) # optional - arb
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    True
    sage: a = RBF(1/3) # optional - arb
    sage: b = RBF(1/3) # optional - arb
    sage: a.is_exact() # optional - arb
    False
    sage: b.is_exact() # optional - arb
    False
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    False

A ball is non-zero if and only if it does not contain zero. ::

    sage: a = RBF(RIF(-0.5, 0.5)) # optional - arb
    sage: bool(a) # optional - arb
    False
    sage: a != 0 # optional - arb
    False
    sage: b = RBF(1/3) # optional - arb
    sage: bool(b) # optional - arb
    True
    sage: b != 0 # optional - arb
    True

A ball ``left`` is less than a ball ``right`` if all elements of
``left`` are less than all elements of ``right``. ::

    sage: a = RBF(RIF(1, 2)) # optional - arb
    sage: b = RBF(RIF(3, 4)) # optional - arb
    sage: a < b # optional - arb
    True
    sage: a <= b # optional - arb
    True
    sage: a > b # optional - arb
    False
    sage: a >= b # optional - arb
    False
    sage: a = RBF(RIF(1, 3)) # optional - arb
    sage: b = RBF(RIF(2, 4)) # optional - arb
    sage: a < b # optional - arb
    False
    sage: a <= b # optional - arb
    False
    sage: a > b # optional - arb
    False
    sage: a >= b # optional - arb
    False

TESTS::

    sage: from sage.rings.real_arb import RBF
    sage: (RBF(pi) * identity_matrix(QQ, 3)).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Real ball field with 53 bits precision

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include "sage/ext/python.pxi"
include "sage/ext/stdsage.pxi"

import operator

import sage.symbolic.constants

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.structure.unique_representation import UniqueRepresentation

cimport sage.rings.integer
cimport sage.rings.rational
cimport sage.structure.element

from sage.libs.arb.arb cimport *
from sage.libs.arb.arf cimport arf_t, arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag
from sage.libs.arb.arf cimport arf_equal, arf_is_nan, arf_is_neg_inf, arf_is_pos_inf
from sage.libs.arb.mag cimport mag_t, mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS, mag_is_inf, mag_is_finite
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_fits_slong_p, mpz_get_ui, mpz_get_si
from sage.libs.mpfi cimport mpfi_get_left, mpfi_get_right, mpfi_interv_fr
from sage.libs.mpfr cimport mpfr_t, mpfr_init2, mpfr_clear, mpfr_sgn, MPFR_PREC_MIN
from sage.libs.mpfr cimport GMP_RNDN, GMP_RNDU, GMP_RNDD, GMP_RNDZ
from sage.rings.real_double cimport RealDoubleElement
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision):
    """
    Convert an MPFI interval to an Arb ball.

    INPUT:

    - ``target`` -- an ``arb_t``.

    - ``source`` -- an ``mpfi_t``.

    - ``precision`` -- an integer `\ge 2`.
    """
    cdef mpfr_t left
    cdef mpfr_t right

    if _do_sig(precision): sig_on()

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    mpfi_get_left(left, source)
    mpfi_get_right(right, source)

    arb_set_interval_mpfr(target,
                          left,
                          right,
                          precision)

    mpfr_clear(left)
    mpfr_clear(right)

cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1:
    """
    Convert an Arb ball to an MPFI interval.

    INPUT:

    - ``target`` -- an ``mpfi_t``.

    - ``source`` -- an ``arb_t``.

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: cython("\n".join([ # optional - arb
        ....:     '#cinclude $SAGE_ROOT/local/include/flint',
        ....:     '#clib arb',
        ....:     'from sage.rings.real_mpfi cimport RealIntervalFieldElement',
        ....:     'from sage.libs.arb.arb cimport *',
        ....:     'from sage.rings.real_arb cimport arb_to_mpfi',
        ....:     'from sage.rings.real_mpfi import RIF',
        ....:     '',
        ....:     'cdef extern from "arb.h":',
        ....:     '    void arb_pow_ui(arb_t y, const arb_t b, unsigned long e, long prec)',
        ....:     '',
        ....:     'cdef RealIntervalFieldElement result',
        ....:     'cdef arb_t arb',
        ....:     'arb_init(arb)',
        ....:     'result = RIF(0)',
        ....:     'arb_set_ui(arb, 65536)',
        ....:     'arb_pow_ui(arb, arb, 65536**3 * 65535, 53)',
        ....:     'arb_to_mpfi(result.value, arb, 53)',
        ....:     'arb_clear(arb)'
        ....: ]))
        Traceback (most recent call last):
        ...
        ArithmeticError: Error converting arb to mpfi. Overflow?
    """
    cdef mpfr_t left
    cdef mpfr_t right

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    try:
        sig_on()
        arb_get_interval_mpfr(left, right, source)
        mpfi_interv_fr(target, left, right)
        sig_off()
    except RuntimeError:
        raise ArithmeticError("Error converting arb to mpfi. Overflow?")
    finally:
        mpfr_clear(left)
        mpfr_clear(right)

class RealBallField(UniqueRepresentation, Parent):
    r"""
    An approximation of the field of real numbers using mid-rad intervals, also
    known as balls.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: RBF = RealBallField() # optional - arb; indirect doctest
        sage: RBF(1) # optional - arb
        1.000000000000000

    ::

        sage: from sage.rings.real_arb import RBF
        sage: (1/2*RBF(1)) + AA(sqrt(2)) - 1 + polygen(QQ, x)
        x + [0.914213562373095 +/- 4.10e-16]
    """
    Element = RealBall

    @staticmethod
    def __classcall__(cls, long precision=53, category=None):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField(53) is RealBallField() # optional - arb
            True
        """
        return super(RealBallField, cls).__classcall__(cls, precision, category)

    def __init__(self, precision, category):
        r"""
        Initialize the real ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(1) # optional - arb
            1.000000000000000
            sage: RealBallField(0) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.
            sage: RealBallField(1) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.
        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        super(RealBallField, self).__init__(
                #category=category or sage.categories.magmas_and_additive_magmas.MagmasAndAdditiveMagmas().Infinite(),
                # FIXME: RBF is not even associative, but CompletionFunctor only works with rings.
                category=category or sage.categories.rings.Rings().Infinite())
        self._prec = precision

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField() # optional - arb
            Real ball field with 53 bits precision
            sage: RealBallField(106) # optional - arb
            Real ball field with 106 bits precision
        """
        return "Real ball field with {} bits precision".format(self._prec)

    def _coerce_map_from_(self, other):
        r"""
        Parents that canonically coerce into real ball fields include:

        - some exact or lazy parents representing subsets of the reals, such as
          ``ZZ``, ``QQ``, ``AA``, and ``RLF``;

        - real ball fields with a larger precision.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().has_coerce_map_from(RealBallField(54))
            True
            sage: RealBallField().has_coerce_map_from(RealBallField(52))
            False
            sage: RealBallField().has_coerce_map_from(RIF) # optional - arb
            False
            sage: RealBallField().has_coerce_map_from(SR) # optional - arb
            False
            sage: RealBallField().has_coerce_map_from(RR) # optional - arb
            False
        """
        from sage.rings.qqbar import AA
        from sage.rings.real_lazy import RLF
        if isinstance(other, RealBallField):
            return (other._prec >= self._prec)
        elif (other is ZZ) or (other is QQ) or (other is AA) or (other is RLF):
            return True
        else:
            return False

    def _element_constructor_(self, mid=None, rad=None):
        """
        Convert ``mid`` to an element of this real ball field, perhaps
        non-canonically.

        In addition to the inputs supported by
        :meth:`ElementConstructor.__init__`,
        anything that is convertible to a real interval can also be used to
        construct a real ball::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(RIF(0, 1))                  # optional - arb; indirect doctest
            [+/- 1.01]
            sage: RBF(1)                          # optional - arb
            1.000000000000000
            sage: RBF(x)                          # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a RealIntervalFieldElement

        Various symbolic constants can be converted without going through real
        intervals. (This is faster and yields tighter error bounds.) ::

            sage: RBF(e)
            [2.718281828459045 +/- 5.49e-16]
            sage: RBF(pi)
            [3.141592653589793 +/- 5.62e-16]
        """
        try:
            return self.element_class(self, mid, rad)
        except TypeError:
            pass

        try:
            return self.element_class(self, mid.pyobject(), rad)
        except (AttributeError, TypeError):
            pass

        try:
            mid = RealIntervalField(self._prec)(mid)
        except TypeError:
            raise TypeError("unable to convert {} to a RealIntervalFieldElement".format(mid))
        return self.element_class(self, mid, rad)

    def gens(self):
        r"""
        EXAMPLE::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().gens()
            (1.000000000000000,)
        """
        return (self.one(),)

    def construction(self):
        """
        Return the construction of a real ball field as a completion of the
        rationals.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField(42) # optional - arb
            sage: functor, base = RBF.construction()
            sage: functor, base
            (Completion[+Infinity], Rational Field)
            sage: functor(base) is RBF
            True
        """
        from sage.categories.pushout import CompletionFunctor
        functor = CompletionFunctor(sage.rings.infinity.Infinity,
                                    self._prec,
                                    {'type': 'Ball'})
        return functor, QQ

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().precision() # optional - arb
            53
        """
        return self._prec

    def is_exact(self):
        """
        Real ball fields are not exact.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().is_exact() # optional - arb
            False
        """
        return False

    def characteristic(self):
        """
        Real ball fields have characteristic zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().characteristic() # optional - arb
            0
        """
        return 0

    # Ball functions of non-ball arguments

    def gamma(self, x):
        """
        Return a ball enclosing the gamma function of ``x``.

        This works even if ``x`` itself is not a ball, and may be more
        efficient in the case where ``x`` is an integer or a rational number.

        .. seealso: :meth`~sage.rings.real_arb.RealBall.gamma`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF.gamma(5)
            24.00000000000000
            sage: RBF.gamma(10**20)
            [+/- 5.92e+1956570551809674821757]
            sage: RBF.gamma(1/3)
            [2.678938534707747 +/- 9.06e-16]
            sage: RBF.gamma(-5)
            nan

        TESTS::

            sage: RBF.gamma(RLF(pi))
            [2.2880377953400 +/- 4.29e-14]
        """
        cdef RealBall res
        cdef sage.rings.integer.Integer x_as_Integer
        cdef sage.rings.rational.Rational x_as_Rational
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        res = self.element_class(self)
        try:
            x_as_Integer = ZZ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, x_as_Integer.value)
                arb_gamma_fmpz(res.value, tmpz, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpz_clear(tmpz)
            return res
        except TypeError:
            pass
        try:
            x_as_Rational = QQ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, x_as_Rational.value)
                arb_gamma_fmpq(res.value, tmpq, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpq_clear(tmpq)
            return res
        except TypeError:
            pass
        return self.coerce(x).gamma()

    def zeta(self, s):
        """
        Return a ball enclosing the Riemann zeta function of ``s``.

        This works even if ``s`` itself is not a ball, and may be more
        efficient in the case where ``s`` is an integer.

        .. seealso: :meth`~sage.rings.real_arb.RealBall.zeta`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF.zeta(3) # abs tol 5e-16
            [1.202056903159594 +/- 2.87e-16]
            sage: RBF.zeta(1)
            nan
            sage: RBF.zeta(1/2)
            [-1.460354508809587 +/- 1.94e-16]
        """
        cdef RealBall res
        cdef sage.rings.integer.Integer s_as_Integer
        try:
            s_as_Integer = ZZ.coerce(s)
            if mpz_fits_ulong_p(s_as_Integer.value):
                res = self.element_class(self)
                if _do_sig(self._prec): sig_on()
                arb_zeta_ui(res.value, mpz_get_ui(s_as_Integer.value), self._prec)
                if _do_sig(self._prec): sig_off()
                return res
        except TypeError:
            pass
        return self.coerce(s).zeta()

    def bernoulli(self, n):
        """
        Return a ball enclosing the ``n``-th Bernoulli number.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: [RBF.bernoulli(n) for n in range(4)]
            [1.000000000000000, -0.5000000000000000, [0.1666666666666667 +/- 7.04e-17], 0]
            sage: RBF.bernoulli(2**20)
            [-1.823002872104961e+5020717 +/- 7.16e+5020701]
            sage: RBF.bernoulli(2**1000)
            Traceback (most recent call last):
            ...
            ValueError: argument too large

        TESTS::

            sage: RBF.bernoulli(2r)
            [0.1666666666666667 +/- 7.04e-17]
            sage: RBF.bernoulli(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Integer Ring
            sage: RBF.bernoulli(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index
        """
        cdef RealBall res
        cdef sage.rings.integer.Integer n_as_Integer = ZZ.coerce(n)
        if mpz_fits_ulong_p(n_as_Integer.value):
            res = self.element_class(self)
            if _do_sig(self._prec): sig_on()
            arb_bernoulli_ui(res.value, mpz_get_ui(n_as_Integer.value), self._prec)
            if _do_sig(self._prec): sig_off()
            return res
        elif n < 0:
            raise ValueError("expected a nonnegative index")
        else:
            # TODO: Fall back to a Sage implementation in this case?
            raise ValueError("argument too large")


cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: _ = RealBallField()(1).psi() # optional - arb; indirect doctest
        sage: _ = RealBallField(1500)(1).psi() # optional - arb
    """
    return (prec > 1000)

cdef inline long prec(RealBall ball):
    return ball._parent._prec

cdef class RealBall(RingElement):
    """
    Hold one ``arb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: a = RealBallField()(RIF(1))                     # optional - arb; indirect doctest
        sage: b = a.psi()                         # optional - arb
        sage: b                                   # optional - arb
        [-0.577215664901533 +/- 3.85e-16]
        sage: b._interval()        # optional - arb
        -0.577215664901533?
    """

    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()(RIF(1)) # optional - arb; indirect doctest
            1.000000000000000
        """
        arb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(1)) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        arb_clear(self.value)

    def __init__(self, parent, mid=None, rad=None):
        """
        Initialize the :class:`RealBall`.

        INPUT:

        - ``parent`` -- a :class:`RealBallField`.

        - ``mid`` (optional) --  ball midpoint, see examples below. If omitted,
          initialize the ball to zero, ignoring the ``rad`` argument.

        - ``rad`` (optional) -- a :class:`RealDoubleElement`, ball radius. If
          the midpoint is not exactly representable in floating-point, the
          radius is adjusted to account for the roundoff error.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField()
            sage: RBF()
            0

        One can create exact real balls using elements of various exact parents,
        or using floating-point numbers::

            sage: RBF(3)
            3.000000000000000
            sage: RBF(3r)
            3.000000000000000
            sage: RBF(1/3)
            [0.3333333333333333 +/- 7.04e-17]
            sage: RBF(3.14)
            [3.140000000000000 +/- 1.25e-16]

        ::

            sage: RBF(3, 0.125r)
            [3e+0 +/- 0.126]
            sage: RBF(pi, 0.125r)
            [3e+0 +/- 0.267]

        Note that integers and floating-point numbers are ''not'' rounded to
        the parent's precision::

            sage: b = RBF(11111111111111111111111111111111111111111111111); b
            [1.111111111111111e+46 +/- 1.12e+30]
            sage: b.mid().exact_rational()
            11111111111111111111111111111111111111111111111

        Similarly, converting a real ball from one real ball field to another
        (with a different precision) only changes the way it is displayed and
        the precision of operations involving it, not the actual representation
        of its center::

            sage: RBF100 = RealBallField(100)
            sage: b100 = RBF100(1/3); b100
            [0.333333333333333333333333333333 +/- 4.65e-31]
            sage: b53 = RBF(b100); b53
            [0.3333333333333333 +/- 3.34e-17]
            sage: RBF100(b53)
            [0.333333333333333333333333333333 +/- 4.65e-31]

        Special values are supported::

            sage: RBF(oo).mid(), RBF(-oo).mid(), RBF(unsigned_infinity).mid()
            (+infinity, -infinity, 0.000000000000000)
            sage: RBF(NaN)
            nan

        TESTS::

            sage: from sage.rings.real_arb import RealBall
            sage: RealBall(RBF, sage.symbolic.constants.Pi()) # abs tol 1e-16
            [3.141592653589793 +/- 5.62e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Log2()) # abs tol 1e-16
            [0.693147180559945 +/- 4.06e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Catalan())
            [0.915965594177219 +/- 9.43e-17]
            sage: RealBall(RBF, sage.symbolic.constants.Khinchin())
            [2.685452001065306 +/- 6.82e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Glaisher())
            [1.282427129100623 +/- 6.02e-16]
            sage: RealBall(RBF, sage.symbolic.constants.e)
            [2.718281828459045 +/- 5.49e-16]
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef arf_t  tmpr
        cdef mag_t  tmpm

        super(RealBall, self).__init__(parent)

        if mid is None:
            return

        elif isinstance(mid, RealBall):
            arb_set(self.value, (<RealBall> mid).value) # no rounding!
        elif isinstance(mid, int):
            arb_set_si(self.value, PyInt_AS_LONG(mid)) # no rounding!
        elif isinstance(mid, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> mid).value)
            arb_set_fmpz(self.value, tmpz) # no rounding!
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.rings.rational.Rational):
            if _do_sig(prec(self)): sig_on()
            fmpq_init(tmpq)
            fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> mid).value)
            arb_set_fmpq(self.value, tmpq, prec(self))
            fmpq_clear(tmpq)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, RealNumber):
            if _do_sig(prec(self)): sig_on()
            arf_init(tmpr)
            arf_set_mpfr(tmpr, (<RealNumber> mid).value)
            arb_set_arf(self.value, tmpr) # no rounding!
            arf_clear(tmpr)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.rings.infinity.AnInfinity):
            if isinstance(mid, sage.rings.infinity.PlusInfinity):
                arb_pos_inf(self.value)
            elif isinstance(mid, sage.rings.infinity.MinusInfinity):
                arb_neg_inf(self.value)
            else:
                arb_zero_pm_inf(self.value)
        elif isinstance(mid, sage.symbolic.constants.Constant):
            if _do_sig(prec(self)): sig_on()
            try:
                if isinstance(mid, sage.symbolic.constants.NotANumber):
                    arb_indeterminate(self.value)
                elif isinstance(mid, sage.symbolic.constants.Pi):
                    arb_const_pi(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Log2):
                    arb_const_log2(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Catalan):
                    arb_const_catalan(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Khinchin):
                    arb_const_khinchin(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Glaisher):
                    arb_const_glaisher(self.value, prec(self))
                else:
                    raise TypeError("unsupported constant")
            finally:
                if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.symbolic.constants_c.E):
            if _do_sig(prec(self)): sig_on()
            arb_const_e(self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, RealIntervalFieldElement):
            mpfi_to_arb(self.value,
                (<RealIntervalFieldElement> mid).value,
                prec(self))
        else:
            raise TypeError("unsupported midpoint type")

        if rad is not None:
            if isinstance(rad, float):
                mag_init(tmpm)
                mag_set_d(tmpm, PyFloat_AS_DOUBLE(rad))
                mag_add(arb_radref(self.value), arb_radref(self.value), tmpm)
                mag_clear(tmpm)
            else:
                raise TypeError("rad should be a Python float")


    cdef RealBall _new(self):
        """
        Return a new real ball element with the same parent as ``self``.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField
            sage: RealBallField()(2)**2 # indirect doctest
            4.000000000000000

        """
        cdef RealBall x
        x = RealBall.__new__(RealBall)
        x._parent = self._parent
        return x

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: from sage.rings.real_arb import RealBallField # optional - arb
           sage: RealBallField()(RIF(1.9, 2))  # optional - arb
           [2e+0 +/- 0.101]
        """
        cdef char* c_result
        cdef bytes py_string

        c_result = arb_get_str(self.value, (prec(self) * 31) // 100, 0)
        try:
            py_string = c_result
        finally:
            flint_free(c_result)

        return py_string

    # Conversions

    cpdef RealIntervalFieldElement _interval(self):
        """
        Return a :mod:`real interval <sage.rings.real_mpfr>` containing this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(2))                   # optional - arb
            sage: RIF(a)                                        # optional - arb, indirect doctest
            2
        """
        cdef RealIntervalFieldElement result
        result = RealIntervalField(prec(self))(0)
        arb_to_mpfi(result.value, self.value, prec(self))
        return result

    def _integer_(self, _):
        """
        Check that this ball contains a single integer and return that integer.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: ZZ(RBF(1, rad=0.1r))
            1
            sage: ZZ(RBF(1, rad=1.0r))
            Traceback (most recent call last):
            ...
            ValueError: [+/- 2.01] does not contain a unique integer
            sage: ZZ(RBF(pi))
            Traceback (most recent call last):
            ...
            ValueError: [3.141592653589793 +/- 5.61e-16] does not contain a unique integer

        """
        cdef sage.rings.integer.Integer res
        cdef fmpz_t tmp
        fmpz_init(tmp)
        try:
            if arb_get_unique_fmpz(tmp, self.value):
                res = sage.rings.integer.Integer.__new__(sage.rings.integer.Integer)
                fmpz_get_mpz(res.value, tmp)
            else:
                raise ValueError("{} does not contain a unique integer".format(self))
        finally:
            fmpz_clear(tmp)
        return res

    def _mpfr_(self, RealField_class field):
        """
        Convert this real ball to a real number.

        This attempts to do something sensible for all rounding modes, as
        illustrated below.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: mypi = RBF(pi)
            sage: RR(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDU')(mypi)
            3.14159265358980
            sage: Reals(rnd='RNDD')(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDZ')(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDZ')(-mypi)
            -3.14159265358979
            sage: Reals(rnd='RNDU')(-mypi)
            -3.14159265358979

        ::

            sage: b = RBF(RIF(-1/2, 1))
            sage: RR(b)
            0.250000000000000
            sage: Reals(rnd='RNDU')(b)
            1.00000000093133
            sage: Reals(rnd='RNDD')(b)
            -0.500000000931323
            sage: Reals(rnd='RNDZ')(b)
            0.250000000000000
        """
        cdef RealNumber left, mid, right
        cdef long prec = field.precision()
        cdef int sl, sr
        if (field.rnd == GMP_RNDN or
                field.rnd == GMP_RNDZ and arb_contains_zero(self.value)):
            mid = RealNumber(field, None)
            sig_str("unable to convert to MPFR (exponent out of range?)")
            arf_get_mpfr(mid.value, arb_midref(self.value), field.rnd)
            sig_off()
            return mid
        else:
            left = RealNumber(field, None)
            right = RealNumber(field, None)
            sig_str("unable to convert to MPFR (exponent out of range?)")
            arb_get_interval_mpfr(left.value, right.value, self.value)
            sig_off()
            if field.rnd == GMP_RNDD:
                return left
            elif field.rnd == GMP_RNDU:
                return right
            elif field.rnd == GMP_RNDZ:
                sl, sr = mpfr_sgn(left.value), mpfr_sgn(left.value)
                if sr > 0 and sl > 0:
                    return left
                elif sr < 0 and sl < 0:
                    return right
                else:
                    return field(0)
        raise ValueError("unknown rounding mode")

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(0).is_zero() # optional - arb
            True
            sage: RBF(RIF(-0.5, 0.5)).is_zero() # optional - arb
            False
        """
        return arb_is_zero(self.value)

    def __nonzero__(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: bool(RBF(pi)) # optional - arb
            True
            sage: bool(RBF(RIF(-0.5, 0.5))) # optional - arb
            False
        """
        return arb_is_nonzero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(1).is_exact() # optional - arb
            True
            sage: RBF(RIF(0.1, 0.2)).is_exact() # optional - arb
            False
        """
        return arb_is_exact(self.value)

    def __richcmp__(left, right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.real_arb`.

        EXAMPLES::

                sage: from sage.rings.real_arb import RealBallField # optional - arb
                sage: RBF = RealBallField() # optional - arb
                sage: a = RBF(1) # optional - arb
                sage: b = RBF(1) # optional - arb
                sage: a is b # optional - arb
                False
                sage: a == b # optional - arb
                True
                sage: a = RBF(1/3) # optional - arb
                sage: a.is_exact() # optional - arb
                False
                sage: b = RBF(1/3) # optional - arb
                sage: b.is_exact() # optional - arb
                False
                sage: a == b # optional - arb
                False
        """
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.real_arb`.

        EXAMPLES::

                sage: from sage.rings.real_arb import RealBallField # optional - arb
                sage: RBF = RealBallField() # optional - arb
                sage: a = RBF(1) # optional - arb
                sage: b = RBF(1) # optional - arb
                sage: a is b # optional - arb
                False
                sage: a == b # optional - arb
                True

        TESTS:

            Balls whose intersection consists of one point::

                sage: a = RBF(RIF(1, 2)) # optional - arb
                sage: b = RBF(RIF(2, 4)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Balls with non-trivial intersection::

                sage: a = RBF(RIF(1, 4)) # optional - arb
                sage: a = RBF(RIF(2, 5)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            One ball contained in another::

                sage: a = RBF(RIF(1, 4)) # optional - arb
                sage: b = RBF(RIF(2, 3)) # optional - arb
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                False
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Disjoint balls::

                sage: a = RBF(1/3) # optional - arb
                sage: b = RBF(1/2) # optional - arb
                sage: a < b # optional - arb
                True
                sage: a <= b # optional - arb
                True
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                False
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                True

            Exact elements::

                sage: a = RBF(2) # optional - arb
                sage: b = RBF(2) # optional - arb
                sage: a.is_exact() # optional - arb
                True
                sage: b.is_exact() # optional - arb
                True
                sage: a < b # optional - arb
                False
                sage: a <= b # optional - arb
                True
                sage: a > b # optional - arb
                False
                sage: a >= b # optional - arb
                True
                sage: a == b # optional - arb
                True
                sage: a != b # optional - arb
                False

            Special values::

                sage: inf = RBF(+infinity)
                sage: other_inf = RBF(+infinity, 42.r)
                sage: neg_inf = RBF(-infinity)
                sage: extended_line = 1/RBF(0)
                sage: exact_nan = inf - inf
                sage: exact_nan.mid(), exact_nan.rad()
                (NaN, 0.00000000)
                sage: other_exact_nan = inf - inf

            ::

                sage: exact_nan == exact_nan, exact_nan <= exact_nan, exact_nan >= exact_nan
                (True, True, True)
                sage: exact_nan != exact_nan, exact_nan < exact_nan, exact_nan > exact_nan
                (False, False, False)
                sage: from operator import eq, ne, le, lt, ge, gt
                sage: ops = [eq, ne, le, lt, ge, gt]
                sage: any(op(exact_nan, other_exact_nan) for op in ops)
                False
                sage: any(op(exact_nan, b) for op in ops for b in [RBF(1), extended_line, inf, neg_inf])
                False

            ::

                sage: neg_inf < a < inf and inf > a > neg_inf
                True
                sage: neg_inf <= b <= inf and inf >= b >= neg_inf
                True
                sage: neg_inf <= extended_line <= inf and inf >= extended_line >= neg_inf
                True
                sage: neg_inf < extended_line or extended_line < inf
                False
                sage: inf > extended_line or extended_line > neg_inf
                False

            ::

                sage: all(b <= b == b >= b and not (b < b or b != b or b > b)
                ....:     for b in [inf, neg_inf, other_inf])
                True
                sage: any(b1 == b2 for b1 in [inf, neg_inf, a, extended_line]
                ....:              for b2 in [inf, neg_inf, a, extended_line]
                ....:              if not b1 is b2)
                False
                sage: all(b1 != b2 and not b1 == b2
                ....:     for b1 in [inf, neg_inf, a]
                ....:     for b2 in [inf, neg_inf, a]
                ....:     if not b1 is b2)
                True
                sage: neg_inf <= -other_inf == neg_inf == -other_inf < other_inf == inf <= other_inf
                True
                sage: any(inf < b or b > inf
                ....:     for b in [inf, other_inf,  a, extended_line])
                False
                sage: any(inf <= b or b >= inf for b in [a, extended_line])
                False
        """
        cdef RealBall lt, rt
        cdef arb_t difference

        lt = left
        rt = right

        if lt is rt:
            return op == Py_EQ or op == Py_GE or op == Py_LE
        elif arb_is_finite(lt.value) or arb_is_finite(rt.value):
            if op == Py_EQ:
                return arb_is_exact(lt.value) and arb_equal(lt.value, rt.value)
            arb_init(difference)
            arb_sub(difference, lt.value, rt.value, prec(lt))
            if op == Py_NE:
                result = arb_is_nonzero(difference)
            elif op == Py_GT:
                result = arb_is_positive(difference)
            elif op == Py_GE:
                result = arb_is_nonnegative(difference)
            elif op == Py_LT:
                result = arb_is_negative(difference)
            elif op == Py_LE:
                result = arb_is_nonpositive(difference)
            arb_clear(difference)
            return result
        elif arf_is_nan(arb_midref(lt.value)) or arf_is_nan(arb_midref(rt.value)):
            return False
        elif mag_is_inf(arb_radref(lt.value)):
            # left is the whole extended real line
            if op == Py_GE:
                return arf_is_neg_inf(arb_midref(rt.value)) and mag_is_finite(arb_radref(rt.value))
            elif op == Py_LE:
                return arf_is_pos_inf(arb_midref(rt.value)) and mag_is_finite(arb_radref(rt.value))
            else:
                return False
        elif mag_is_inf(arb_radref(rt.value)):
            # right is the whole extended real line
            if op == Py_GE:
                return arf_is_pos_inf(arb_midref(lt.value)) and mag_is_finite(arb_radref(lt.value))
            elif op == Py_LE:
                return arf_is_neg_inf(arb_midref(lt.value)) and mag_is_finite(arb_radref(lt.value))
            else:
                return False
        else:
            # both left and right are special, neither is nan, and neither is
            # [-∞,∞], so they are both points at infinity
            if op == Py_EQ:
                return arf_equal(arb_midref(lt.value), arb_midref(rt.value))
            elif op == Py_NE:
                return not arf_equal(arb_midref(lt.value), arb_midref(rt.value))
            elif op == Py_GT:
                return (arf_is_pos_inf(arb_midref(lt.value))
                        and arf_is_neg_inf(arb_midref(rt.value)))
            elif op == Py_LT:
                return (arf_is_neg_inf(arb_midref(lt.value))
                        and arf_is_pos_inf(arb_midref(rt.value)))
            elif op == Py_GE or op == Py_LE:
                return True
        assert False, "not reached"

    # Center and radius

    def mid(self):
        """
        Return the center of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField, RBF # optional - arb
            sage: RealBallField(16)(1/3).mid()
            0.3333
            sage: RealBallField(16)(1/3).mid().parent()
            Real Field with 15 bits of precision

        ::

            sage: b = RBF(2)^(2^1000)
            sage: b.mid()
            Traceback (most recent call last):
            ...
            RuntimeError: unable to convert to MPFR (exponent out of range?)
        """
        cdef long mid_prec = arb_bits(self.value) or prec(self)
        if mid_prec < MPFR_PREC_MIN:
            mid_prec = MPFR_PREC_MIN
        cdef RealField_class mid_field = RealField(mid_prec)
        return self._mpfr_(mid_field)

    def rad(self):
        """
        Return the radius of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()(1/3).rad()
            5.5511151e-17
            sage: RealBallField()(1/3).rad().parent()
            Real Field with 30 bits of precision
        """
        # Should we return a real number with rounding towards +∞ (or away from
        # zero if/when implemented)?
        cdef RealField_class rad_field = RealField(MAG_BITS)
        cdef RealNumber rad = RealNumber(rad_field, None)
        cdef arf_t tmp
        arf_init(tmp)
        arf_set_mag(tmp, arb_radref(self.value))
        cdef int rnd = arf_get_mpfr(rad.value, tmp, GMP_RNDN)
        arf_clear(tmp)
        if rnd != 0:
            raise OverflowError("Unable to represent the radius of this ball within the exponent range of RealNumbers")
        return rad

    # Precision

    def round(self):
        """
        Return a copy of this ball with center rounded to the precision of the
        parent.

        .. SEEALSO:: :meth:`trim`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: b = RBF(pi.n(100))
            sage: b.mid()
            3.141592653589793238462643383
            sage: b.round().mid()
            3.1415926535898
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_set_round(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def accuracy(self):
        """
        Return the effective relative accuracy of this ball measured in bits.

        The accuracy is defined as the difference between the position of the
        top bit in the midpoint and the top bit in the radius and , minus one.
        The result is clamped between plus/minus ``ARF_PREC_EXACT``.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi).accuracy()
            51
            sage: RBF(1).accuracy()
            9223372036854775807
            sage: RBF(NaN).accuracy()
            -9223372036854775807
        """
        return arb_rel_accuracy_bits(self.value)

    def trim(self):
        """
        Return a trimmed copy of this ball.

        Round ``self`` to a number of bits equal to the :meth:`accuracy` of
        ``self`` (as indicated by its radius), plus a few guard bits. The
        resulting ball is guaranteed to contain ``self``, but is more economical
        if ``self`` has less than full accuracy.

        .. SEEALSO:: :meth:`round`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: b = RBF(RIF(3.1415,3.1416))
            sage: b.mid()
            3.14155000000000
            sage: b.trim().mid()
            3.14155000
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_trim(res.value, self.value)
        if _do_sig(prec(self)): sig_off()
        return res

    # Comparisons and predicates

    def is_zero(self):
        """
        Return True iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(0).is_zero()
            True
            sage: RBF(0, rad=0.25r).is_zero()
            False

        """
        return bool(arb_is_zero(self.value))

    def is_nonzero(self):
        """
        Return True iff zero is not contained in the interval represented
        by this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi).is_nonzero()
            True
            sage: RBF(1, rad=2.r).is_nonzero()
            False
        """
        return bool(arb_is_nonzero(self.value))

    def is_finite(self):
        """
        Return True iff the midpoint and radius of this ball are both
        finite floating-point numbers, i.e. not infinities or NaN.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: (RBF(2)^(2^1000)).is_finite()
            True
            sage: RBF(oo).is_finite()
            False
        """
        return bool(arb_is_finite(self.value))

    def is_exact(self):
        """
        Return True iff the radius of this ball is zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).is_exact()
            True
            sage: RBF(pi).is_exact()
            False
        """
        return bool(arb_is_exact(self.value))

    def identical(self, RealBall other):
        """
        Return True iff ``self`` and ``other`` are equal as balls, i.e.
        have both the same midpoint and radius.

        Note that this is not the same thing as testing whether both ``self``
        and ``other`` certainly represent the same real number, unless either
        ``self`` or ``other`` is exact (and neither contains NaN). To test
        whether both operands might represent the same mathematical quantity,
        use :meth:`overlaps` or :meth:`contains`, depending on the
        circumstance.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).identical(RBF(3)-RBF(2))
            True
            sage: RBF(1, rad=0.25r).identical(RBF(1, rad=0.25r))
            True
            sage: RBF(1).identical(RBF(1, rad=0.25r))
            False
        """
        return bool(arb_equal(self.value, other.value))

    def overlaps(self, RealBall other):
        """
        Return True iff ``self`` and ``other`` have some point in common.

        If either ``self`` or ``other`` contains NaN, this method always
        returns nonzero (as a NaN could be anything, it could in particular
        contain any number that is included in the other operand).

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi).overlaps(RBF(pi) + 2**(-100))
            True
            sage: RBF(pi).overlaps(RBF(3))
            False
        """
        return bool(arb_overlaps(self.value, other.value))

    def contains_exact(self, other):
        """
        Returns nonzero *iff* the given number (or ball) ``other`` is contained
        in the interval represented by ``self``.

        If ``self`` contains NaN, this function always returns nonzero (as
        it could represent anything, and in particular could represent all the
        points included in ``other``). If ``other`` contains NaN and ``self``
        does not, it always returns zero.

        Use ``other in self`` for a test that works for a wider range of inputs
        but may return false negatives.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: b = RBF(1)
            sage: b.contains_exact(1)
            True
            sage: b.contains_exact(QQ(1))
            True
            sage: b.contains_exact(1.)
            True
            sage: b.contains_exact(b)
            True

        ::

            sage: RBF(1/3).contains_exact(1/3)
            True
            sage: RBF(sqrt(2)).contains_exact(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError

        TESTS::

            sage: b.contains_exact(1r)
            True

        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        if _do_sig(prec(self)): sig_on()
        try:
            if isinstance(other, RealBall):
                res = arb_contains(self.value, (<RealBall> other).value)
            elif isinstance(other, int):
                res = arb_contains_si(self.value, PyInt_AS_LONG(other))
            elif isinstance(other, sage.rings.integer.Integer):
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> other).value)
                res = arb_contains_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
            elif isinstance(other, sage.rings.rational.Rational):
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> other).value)
                res = arb_contains_fmpq(self.value, tmpq)
                fmpq_clear(tmpq)
            elif isinstance(other, RealNumber):
                res = arb_contains_mpfr(self.value, (<RealNumber> other).value)
            else:
                raise TypeError
        finally:
            if _do_sig(prec(self)): sig_off()
        return bool(res)

    def __contains__(self, other):
        """
        Return True if ``other`` can be verified to be contained in ``self``.

        The test is done using interval arithmetic with a precision determined
        by the parent of ``self`` and may return false negatives.

        .. SEEALSO:: :meth:`contains_exact`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF, RealBallField
            sage: sqrt(2) in RBF(sqrt(2))
            True

        A false negative::

            sage: sqrt(2) in RBF(RealBallField(100)(sqrt(2)))
            False
        """
        return self.contains_exact(self._parent(other))

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: -RBF(1/3)
            [-0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall res = self._new()
        arb_neg(res.value, self.value)
        return res

    def __invert__(self):
        """
        Return the inverse of this ball.

        The result is guaranteed to contain the inverse of any point of the
        input ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: ~RBF(5)
            [0.2000000000000000 +/- 4.45e-17]
            sage: ~RBF(0)
            [+/- inf]
            sage: RBF(RIF(-0.1,0.1))
            [+/- 0.101]

        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_inv(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef ModuleElement _add_(self, ModuleElement other):
        """
        Return the sum of two balls, rounded to the ambient field's precision.

        The resulting ball is guaranteed to contain the sums of any two points
        of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1) + RBF(1/3)
            [1.333333333333333 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_add(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef ModuleElement _sub_(self, ModuleElement other):
        """
        Return the difference of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the differences of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1) - RBF(1/3)
            [0.666666666666667 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sub(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _mul_(self, RingElement other):
        """
        Return the product of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the products of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(-2) * RBF(1/3)
            [-0.666666666666667 +/- 4.82e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_mul(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _div_(self, RingElement other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi)/RBF(e)
            [1.155727349790922 +/- 8.49e-16]
            sage: RBF(2)/RBF(0)
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_div(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __pow__(base, expo, _):
        """
        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(e)^17
            [24154952.753575 +/- 3.47e-7]
            sage: RBF(e)^(-1)
            [0.367879441171442 +/- 4.52e-16]
            sage: RBF(e)^(1/2)
            [1.648721270700128 +/- 5.00e-16]
            sage: RBF(e)^RBF(pi)
            [23.1406926327793 +/- 9.20e-14]

        ::

            sage: RBF(-1)^(1/3)
            nan
            sage: RBF(0)^(-1)
            [+/- inf]
            sage: RBF(-e)**RBF(pi)
            nan

        TESTS::

            sage: RBF(e)**(2r)
            [7.38905609893065 +/- 4.75e-15]
            sage: RBF(e)**(-1r)
            [0.367879441171442 +/- 4.52e-16]
        """
        cdef fmpz_t tmpz
        if not isinstance(base, RealBall):
            return sage.structure.element.bin_op(base, expo, operator.pow)
        cdef RealBall self = base
        cdef RealBall res = self._new()
        if isinstance(expo, int) and expo > 0:
            if _do_sig(prec(self)): sig_on()
            arb_pow_ui(res.value, self.value, PyInt_AS_LONG(expo), prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> expo).value)
            arb_pow_fmpz(res.value, self.value, tmpz, prec(self))
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, RealBall):
            if _do_sig(prec(self)): sig_on()
            arb_pow(res.value, self.value, (<RealBall> expo).value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            return sage.structure.element.bin_op(base, expo, operator.pow)
        return res

    def sqrt(self):
        """
        Return the square root of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(2).sqrt()
            [1.414213562373095 +/- 2.99e-16]
            sage: RBF(-1/3).sqrt()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sqrtpos(self):
        """
        Return the square root of this ball, assuming that it represents a
        nonnegative number.

        Any negative numbers in the input interval are discarded, and the
        output ball does not contain any negative numbers (unless the radius is
        infinite).

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(2).sqrtpos()
            [1.414213562373095 +/- 2.99e-16]
            sage: RBF(-1/3).sqrtpos()
            0
            sage: RBF(0, rad=2.r).sqrtpos()
            [+/- 1.42]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrtpos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rsqrt(self):
        """
        Return the reciprocal square root of ``self``.

        At high precision, this is faster than computing a square root.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(2).rsqrt()
            [0.707106781186547 +/- 5.73e-16]
            sage: RBF(0).rsqrt()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_rsqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sqrt1pm1(self):
        """
        Return `\sqrt{1+\mathrm{self}}-1`, computed accurately when ``self`` is
        close to zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: eps = RBF(10^(-20))
            sage: (1 + eps).sqrt() - 1
            [+/- 1.12e-16]
            sage: eps.sqrt1pm1()
            [5.00000000000000e-21 +/- 2.54e-36]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrt1pm1(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Elementary functions

    def log(self):
        """
        Return the natural logarithm of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(3).log()
            [1.098612288668110 +/- 6.63e-16]
            sage: RBF(-1/3).log()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_log(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def log1p(self):
        """
        Return ``log(1 + self)``, computed accurately when ``self`` is close to
        zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: eps = RBF(1e-30)
            sage: (1 + eps).log()
            [+/- 2.23e-16]
            sage: eps.log1p()
            [1.00000000000000e-30 +/- 2.68e-46]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_log1p(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def exp(self):
        """
        Return the exponential of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).exp()
            [2.718281828459045 +/- 5.41e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_exp(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def expm1(self):
        """
        Return ``exp(self) - 1``, computed accurately when ``self`` is close to
        zero.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: eps = RBF(1e-30)
            sage: exp(eps) - 1
            [+/- 3.16e-30]
            sage: eps.expm1()
            [1.000000000000000e-30 +/- 8.34e-47]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_expm1(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sin(self):
        """
        Return the sine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi).sin() # abs tol 1e-16
            [+/- 5.69e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cos(self):
        """
        Return the cosine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(pi).cos() # abs tol 1e-16
            [-1.00000000000000 +/- 6.69e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tan(self):
        """
        Return the tangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).tan()
            [1.557407724654902 +/- 3.26e-16]
            sage: RBF(pi/2).tan()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_tan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cot(self):
        """
        Return the cotangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).cot()
            [0.642092615934331 +/- 4.79e-16]
            sage: RBF(pi).cot()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cot(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsin(self):
        """
        Return the arcsine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).arcsin()
            [1.570796326794897 +/- 6.65e-16]
            sage: RBF(1, rad=.125r).arcsin()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_asin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccos(self):
        """
        Return the arccosine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).arccos()
            0
            sage: RBF(1, rad=.125r).arccos()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_acos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctan(self):
        """
        Return the arctangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).arctan()
            [0.785398163397448 +/- 3.91e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_atan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sinh(self):
        """
        Return the hyperbolic sine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).sinh()
            [1.175201193643801 +/- 6.18e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cosh(self):
        """
        Return the hyperbolic cosine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).cosh()
            [1.543080634815244 +/- 5.28e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tanh(self):
        """
        Return the hyperbolic tangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).tanh()
            [0.761594155955765 +/- 2.81e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_tanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def coth(self):
        """
        Return the hyperbolic cotangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).coth()
            [1.313035285499331 +/- 4.97e-16]
            sage: RBF(0).coth()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_coth(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).arcsinh()
            [0.881373587019543 +/- 1.87e-16]
            sage: RBF(0).arcsinh()
            0
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_asinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(2).arccosh()
            [1.316957896924817 +/- 6.61e-16]
            sage: RBF(1).arccosh()
            0
            sage: RBF(0).arccosh()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_acosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of this ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(0).arctanh()
            0
            sage: RBF(1/2).arctanh()
            [0.549306144334055 +/- 3.32e-16]
            sage: RBF(1).arctanh()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_atanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Special functions

    def gamma(self):
        """
        Return the image of this ball by the Euler Gamma function.

        For integer and rational arguments,
        :meth`~sage.rings.real_arb.RealBall.gamma` may be faster.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1/2).gamma()
            [1.772453850905516 +/- 3.72e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_gamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def log_gamma(self):
        """
        Return the image of this ball by the logarithmic Gamma function.

        The complex branch structure is assumed, so if ``self`` ≤ 0, the result
        is an indeterminate interval.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1/2).log_gamma()
            [0.572364942924700 +/- 4.87e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_lgamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rgamma(self):
        """
        Return the image of this ball by the function 1/Γ, avoiding division by
        zero at the poles of the gamma function.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(-1).rgamma()
            0
            sage: RBF(3).rgamma()
            0.5000000000000000
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_rgamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RealBall psi(self):
        """
        Compute the digamma function with argument self.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(1).psi()  # optional - arb
            [-0.577215664901533 +/- 3.85e-16]
        """

        cdef RealBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_digamma(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def zeta(self, a=None):
        """
        Return the image of this ball by the Hurwitz zeta function.

        For ``a = 1`` (or ``a = None``), this computes the Riemann zeta function.

        Use :meth`~sage.rings.real_arb.RealBallField.zeta` to compute the
        Riemann zeta function of a small integer without first converting it to
        a real ball.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: RBF(-1).zeta()
            [-0.0833333333333333 +/- 4.36e-17]
            sage: RBF(-1).zeta(1)
            [-0.0833333333333333 +/- 6.81e-17]
            sage: RBF(-1).zeta(2) # abs tol 1e-16
            [-1.083333333333333 +/- 4.09e-16]
        """
        cdef RealBall a_ball
        cdef RealBall res = self._new()
        if a is None:
            if _do_sig(prec(self)): sig_on()
            arb_zeta(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            a_ball = self._parent.coerce(a)
            if _do_sig(prec(self)): sig_on()
            arb_hurwitz_zeta(res.value, self.value, a_ball.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def polylog(self, s):
        """
        Return the polylogarithm `\operatorname{Li}_s(\mathrm{self})`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: polylog(0, -1)
            -1/2
            sage: RBF(-1).polylog(0)
            [-0.50000000000000 +/- 1.78e-15]
            sage: polylog(1, 1/2)
            -log(1/2)
            sage: RBF(1/2).polylog(1)
            [0.6931471805599 +/- 5.08e-14]
            sage: RBF(1/3).polylog(1/2)
            [0.44210883528067 +/- 6.75e-15]
            sage: RBF(1/3).polylog(RLF(pi))
            [0.34728895057225 +/- 5.51e-15]

        TESTS::

            sage: RBF(1/3).polylog(2r)
            [0.36621322997706 +/- 4.62e-15]
        """
        cdef RealBall s_as_ball
        cdef sage.rings.integer.Integer s_as_Integer
        cdef RealBall res = self._new()
        try:
            s_as_Integer = ZZ.coerce(s)
            if mpz_fits_slong_p(s_as_Integer.value):
                if _do_sig(prec(self)): sig_on()
                arb_polylog_si(res.value, mpz_get_si(s_as_Integer.value), self.value, prec(self))
                if _do_sig(prec(self)): sig_off()
                return res
        except TypeError:
            pass
        s_as_ball = self._parent.coerce(s)
        if _do_sig(prec(self)): sig_on()
        arb_polylog(res.value, s_as_ball.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

RBF = RealBallField()
