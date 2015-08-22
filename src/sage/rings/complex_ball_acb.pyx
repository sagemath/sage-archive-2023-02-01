r"""
Arbitrary Precision Complex Intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`sage.rings.real_arb`
    - :mod:`sage.rings.complex_interval_field`
    - :mod:`sage.rings.complex_interval`

Data Structure
==============

A :class:`ComplexBall` represents a complex number with error bounds. It wraps
an Arb object of type ``acb_t``, which  consists of a pair of real number balls
representing the real and imaginary part with separate error bounds.

A :class:`ComplexBall` thus represents a rectangle `[m_1-r_1, m_1+r_1] +
[m_2-r_2, m_2+r_2] i` in the complex plane. This is used in Arb instead of a
disk or square representation (consisting of a complex floating-point midpoint
with a single radius), since it allows implementing many operations more
conveniently by splitting into ball operations on the real and imaginary parts.
It also allows tracking when complex numbers have an exact (for example exactly
zero) real part and an inexact imaginary part, or vice versa.

Comparison
==========

.. WARNING::

    Identical :class:`ComplexBall` objects are understood to give
    permission for algebraic simplification. This assumption is made
    to improve performance. For example, setting ``z = x*x`` sets `z`
    to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: from sage.rings.complex_ball_acb import CBF
    doctest:...: FutureWarning: This class/method/function is marked as experimental.
    It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/17218 for details.
    sage: a = CBF(1, 2)
    sage: b = CBF(1, 2)
    sage: a is b
    False
    sage: a == b
    True
    sage: a = CBF(1/3, 1/5)
    sage: b = CBF(1/3, 1/5)
    sage: a.is_exact()
    False
    sage: b.is_exact()
    False
    sage: a is b
    False
    sage: a == b
    False

A ball is non-zero if and only if it does not contain zero. ::

    sage: a = CBF(RIF(-0.5, 0.5))
    sage: bool(a)
    False
    sage: a != 0
    False
    sage: b = CBF(1/3, 1/5)
    sage: bool(b)
    True
    sage: b != 0
    True

Coercion
========

Automatic coercions work as expected::

    sage: from sage.rings.real_arb import RealBallField
    sage: bpol = 1/3*CBF(i) + AA(sqrt(2)) + (polygen(RealBallField(20), 'x') + QQbar(i))
    sage: bpol
    x + [1.41421 +/- 5.09e-6] + [1.33333 +/- 3.97e-6]*I
    sage: bpol.parent()
    Univariate Polynomial Ring in x over Complex ball field with 20 bits precision

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

import sage.categories.fields

cimport sage.rings.integer
cimport sage.rings.rational

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_hypgeom cimport *
from sage.libs.arb.arf cimport arf_t, arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag, arf_set
from sage.libs.arb.mag cimport mag_t, mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS, mag_is_inf, mag_is_finite, mag_zero
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.misc.superseded import experimental
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi
from sage.rings.real_arb import RealBallField
from sage.structure.element cimport Element, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation

cdef inline bint acb_is_nonzero(const acb_t z):
    return arb_is_nonzero(&z.real) or arb_is_nonzero(&z.imag)

cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source):
    """
    Convert a :class:`ComplexIntervalFieldElement` to an ``acb``.

    INPUT:

    - ``target`` -- an ``acb_t``

    - ``source`` -- a :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef long precision
    precision = source.parent().precision()
    mpfi_to_arb(&target.real, source.__re, precision)
    mpfi_to_arb(&target.imag, source.__im, precision)

cdef int acb_to_ComplexIntervalFieldElement(
    ComplexIntervalFieldElement target,
    const acb_t source) except -1:
    """
    Convert an ``acb`` to a :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``target`` -- a :class:`ComplexIntervalFieldElement`

    - ``source`` -- an ``acb_t``

    OUTPUT:

    A :class:`ComplexIntervalFieldElement`.
    """
    cdef long precision = target._prec

    arb_to_mpfi(target.__re, &source.real, precision)
    arb_to_mpfi(target.__im, &source.imag, precision)
    return 0

class ComplexBallField(UniqueRepresentation, Parent):
    r"""
    An approximation of the field of complex numbers using pairs of mid-rad
    intervals.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: from sage.rings.complex_ball_acb import ComplexBallField
        sage: CBF = ComplexBallField() # indirect doctest
        sage: CBF(1)
        1.000000000000000

    TESTS::

        sage: ComplexBallField(0)
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
        sage: ComplexBallField(1)
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
    """
    Element = ComplexBall

    @staticmethod
    def __classcall__(cls, long precision=53, category=None):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField(53) is ComplexBallField()
            True
        """
        return super(ComplexBallField, cls).__classcall__(cls, precision, category)

    @experimental(17218)
    def __init__(self, precision, category):
        r"""
        Initialize the complex ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: CBF(1)
            1.000000000000000

        TESTS::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.base()
            Real ball field with 53 bits precision
            sage: CBF.base_ring()
            Real ball field with 53 bits precision

        There are direct coercions from ZZ and QQ (for which arb provides
        construction functions)::

            sage: CBF.coerce_map_from(ZZ)
            Conversion map:
            From: Integer Ring
            To:   Complex ball field with 53 bits precision
            sage: CBF.coerce_map_from(QQ)
            Conversion map:
            From: Rational Field
            To:   Complex ball field with 53 bits precision

        Various other coercions are available through real ball fields or CLF::

            sage: CBF.coerce_map_from(RLF)
            Composite map:
            From: Real Lazy Field
            To:   Complex ball field with 53 bits precision
            Defn:   Conversion map:
                    From: Real Lazy Field
                    To:   Real ball field with 53 bits precision
                    then
                    Conversion map:
                    From: Real ball field with 53 bits precision
                    To:   Complex ball field with 53 bits precision
            sage: CBF.has_coerce_map_from(AA)
            True
            sage: CBF.has_coerce_map_from(QuadraticField(-1))
            True
            sage: CBF.has_coerce_map_from(QQbar)
            True
            sage: CBF.has_coerce_map_from(CLF)
            True
        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        real_field = RealBallField(precision)
        super(ComplexBallField, self).__init__(
                base=real_field,
                category=category or [sage.categories.fields.Fields()])
        self._prec = precision
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.real_lazy import CLF
        self._populate_coercion_lists_([ZZ, QQ, real_field, CLF])

    def _real_field(self):
        """
        TESTS::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF._real_field()
            Real ball field with 53 bits precision
        """
        return self._base

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField()
            Complex ball field with 53 bits precision
            sage: ComplexBallField(106)
            Complex ball field with 106 bits precision
        """
        return "Complex ball field with {} bits precision".format(self._prec)

    def construction(self):
        """
        Returns the construction of a complex ball field as the algebraic
        closure of the real ball field with the same precision.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: functor, base = CBF.construction()
            sage: functor, base
            (AlgebraicClosureFunctor, Real ball field with 53 bits precision)
            sage: functor(base) is CBF
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._real_field())

    def ngens(self):
        r"""
        EXAMPLE::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.ngens()
            1
        """
        return 1

    def gen(self, i):
        r"""
        EXAMPLE::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.0
            1.000000000000000*I
            sage: CBF.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: only one generator
        """
        if i == 0:
            return self(0, 1)
        else:
            raise ValueError("only one generator")

    def gens(self):
        r"""
        EXAMPLE::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.gens()
            (1.000000000000000*I,)
            sage: CBF.gens_dict()
            {'1.000000000000000*I': 1.000000000000000*I}
        """
        return (self(0, 1),)

    def _coerce_map_from_(self, other):
        r"""
        Parents that canonically coerce into complex ball fields include:

        - anything that coerces into the corresponding real ball field;

        - real and complex ball fields with a larger precision;

        - various exact or lazy parents representing subsets of the complex
          numbers, such as ``QQbar``, ``CLF``, and number fields equipped
          with complex embeddings.

        TESTS::

            sage: from sage.rings.complex_ball_acb import CBF, ComplexBallField
            sage: from sage.rings.real_arb import RBF, RealBallField
            sage: CBF.coerce_map_from(CBF)
            Identity endomorphism of Complex ball field with 53 bits precision
            sage: CBF.coerce_map_from(ComplexBallField(100))
            Conversion map:
            From: Complex ball field with 100 bits precision
            To:   Complex ball field with 53 bits precision
            sage: CBF.has_coerce_map_from(ComplexBallField(42))
            False
            sage: CBF.has_coerce_map_from(RealBallField(54))
            True
            sage: CBF.has_coerce_map_from(RealBallField(52))
            False

        Check that there are no coercions from interval or floating-point parents::

            sage: CBF.has_coerce_map_from(RIF)
            False
            sage: CBF.has_coerce_map_from(CIF)
            False
            sage: CBF.has_coerce_map_from(RR)
            False
            sage: CBF.has_coerce_map_from(CC)
            False
        """
        if isinstance(other, (RealBallField, ComplexBallField)):
            return (other._prec >= self._prec)

    def _element_constructor_(self, x=None, y=None):
        r"""
        Convert (x, y) to an element of this complex ball field, perhaps
        non-canonically.

        INPUT:

        - ``x``, ``y`` (optional) -- either a complex number, interval or ball,
          or two real ones (see examples below for more information on accepted
          number types).

        .. SEEALSO:: :meth:`sage.rings.real_arb.RealBallField._element_constructor_`

        EXAMPLES::

            sage: from sage.rings.real_arb import RBF
            sage: from sage.rings.complex_ball_acb import CBF, ComplexBallField
            sage: CBF()
            0
            sage: CBF(1) # indirect doctest
            1.000000000000000
            sage: CBF(1, 1)
            1.000000000000000 + 1.000000000000000*I
            sage: CBF(pi, sqrt(2))
            [3.141592653589793 +/- 5.61e-16] + [1.414213562373095 +/- 4.10e-16]*I
            sage: CBF(I)
            1.000000000000000*I
            sage: CBF(pi+I/3)
            [3.141592653589793 +/- 5.61e-16] + [0.3333333333333333 +/- 7.04e-17]*I
            sage: CBF(QQbar(i/7))
            [0.1428571428571428 +/- 9.09e-17]*I
            sage: CBF(AA(sqrt(2)))
            [1.414213562373095 +/- 4.10e-16]
            sage: CBF(CIF(0, 1))
            1.000000000000000*I
            sage: CBF(RBF(1/3))
            [0.3333333333333333 +/- 7.04e-17]
            sage: CBF(RBF(1/3), RBF(1/6))
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: CBF(1/3)
            [0.3333333333333333 +/- 7.04e-17]
            sage: CBF(1/3, 1/6)
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: ComplexBallField(106)(1/3, 1/6)
            [0.33333333333333333333333333333333 +/- 6.94e-33] + [0.16666666666666666666666666666666 +/- 7.70e-33]*I
            sage: CBF(infinity, NaN)
            [+/- inf] + nan*I
            sage: CBF(x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a ComplexBall

        TESTS::

            sage: CBF(1+I, 2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert I + 1 to a RealBall
        """
        try:
            return self.element_class(self, x, y)
        except TypeError:
            pass

        if y is None:
            try:
                x = self._real_field()(x)
                return self.element_class(self, x)
            except (TypeError, ValueError):
                pass
            try:
                y = self._real_field()(x.imag())
                x = self._real_field()(x.real())
                return self.element_class(self, x, y)
            except (AttributeError, TypeError):
                pass
            try:
                x = ComplexIntervalField(self._prec)(x)
                return self.element_class(self, x)
            except TypeError:
                pass
            raise TypeError("unable to convert {} to a ComplexBall".format(x))
        else:
            x = self._real_field()(x)
            y = self._real_field()(y)
            return self.element_class(self, x, y)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.an_element() # indirect doctest
            [0.3333333333333333 +/- 1.49e-17] - [0.1666666666666667 +/- 4.26e-17]*I
        """
        return self(1.0/3, -1.0/6)

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField().precision()
            53
        """
        return self._prec

    def is_exact(self):
        """
        Complex ball fields are not exact.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField().is_exact()
            False
        """
        return False

    def is_finite(self):
        """
        Complex ball fields are infinite.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField().is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Complex ball fields have characteristic zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField().characteristic()
            0
        """
        return 0

    def some_elements(self):
        """
        Complex ball fields contain elements with exact, inexact, infinite, or
        undefined real and imaginary parts.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF.some_elements()
                [1.000000000000000,
                -0.5000000000000000*I,
                1.000000000000000 + [0.3333333333333333 +/- 1.49e-17]*I,
                [-0.3333333333333333 +/- 1.49e-17] + 0.2500000000000000*I,
                [-2.175556475109056e+181961467118333366510562 +/- 1.29e+181961467118333366510545],
                [+/- inf],
                [+/- inf]*I,
                [+/- inf] + [+/- inf]*I,
                nan,
                nan + nan*I,
                [+/- inf] + nan*I]
        """
        return [self(1), self(0, -1./2), self(1, 1./3), self(-1./3, 1./4),
                -self(1, 1)**(sage.rings.integer.Integer(2)**80),
                self('inf'), self(1/3, 'inf'), self('inf', 'inf'),
                self('nan'), self('nan', 'nan'), self('inf', 'nan')]

cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: from sage.rings.complex_ball_acb import ComplexBallField
    """
    return (prec > 1000)

cdef inline long prec(ComplexBall ball):
    return ball._parent._prec

cdef inline Parent real_ball_field(ComplexBall ball):
    return ball._parent._real_field()

cdef class ComplexBall(RingElement):
    """
    Hold one ``acb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: from sage.rings.complex_ball_acb import ComplexBallField
        sage: a = ComplexBallField()(1, 1)
        sage: a
        1.000000000000000 + 1.000000000000000*I
        sage: a._interval()
        1 + 1*I
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: ComplexBallField(2)(0) # indirect doctest
            0
        """
        acb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: a = ComplexBallField(2)(0) # indirect doctest
            sage: del a
        """
        acb_clear(self.value)

    def __init__(self, parent, x=None, y=None):
        """
        Initialize the :class:`ComplexBall`.

        INPUT:

        - ``parent`` -- a :class:`ComplexBallField`.

        - ``x``, ``y`` (optional) -- either a complex number, interval or ball,
          or two real ones.

        .. SEEALSO:: :meth:`ComplexBallField._element_constructor_`

        TESTS::

            sage: from sage.rings.complex_ball_acb import ComplexBallField, ComplexBall
            sage: from sage.rings.real_arb import RBF
            sage: CBF53, CBF100 = ComplexBallField(53), ComplexBallField(100)
            sage: ComplexBall(CBF100)
            0
            sage: ComplexBall(CBF100, ComplexBall(CBF53, ComplexBall(CBF100, 1/3)))
            [0.333333333333333333333333333333 +/- 4.65e-31]
            sage: ComplexBall(CBF100, RBF(pi))
            [3.141592653589793 +/- 5.61e-16]
            sage: ComplexBall(CBF100, -3r)
            -3.000000000000000000000000000000
            sage: ComplexBall(CBF100, 10^100)
            1.000000000000000000000000000000e+100
            sage: ComplexBall(CBF100, CIF(1, 2))
            1.000000000000000000000000000000 + 2.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, RBF(1/3), RBF(1))
            [0.3333333333333333 +/- 7.04e-17] + 1.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, 1, 2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported initializer
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq

        Element.__init__(self, parent)

        if x is None:
            return
        elif y is None:
            if isinstance(x, ComplexBall):
                acb_set(self.value, (<ComplexBall> x).value)
            elif isinstance(x, RealBall):
                acb_set_arb(self.value, (<RealBall> x).value)
            elif isinstance(x, int):
                acb_set_si(self.value, PyInt_AS_LONG(x))
            elif isinstance(x, sage.rings.integer.Integer):
                if _do_sig(prec(self)): sig_on()
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> x).value)
                acb_set_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
                if _do_sig(prec(self)): sig_off()
            elif isinstance(x, sage.rings.rational.Rational):
                if _do_sig(prec(self)): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> x).value)
                acb_set_fmpq(self.value, tmpq, prec(self))
                fmpq_clear(tmpq)
                if _do_sig(prec(self)): sig_off()
            elif isinstance(x, ComplexIntervalFieldElement):
                ComplexIntervalFieldElement_to_acb(self.value,
                                                   <ComplexIntervalFieldElement> x)
            else:
                raise TypeError("unsupported initializer")
        elif isinstance(x, RealBall) and isinstance(y, RealBall):
            arb_set(acb_realref(self.value), (<RealBall> x).value)
            arb_set(acb_imagref(self.value), (<RealBall> y).value)
        else:
            raise TypeError("unsupported initializer")

    cdef ComplexBall _new(self):
        """
        Return a new complex ball element with the same parent as ``self``.
        """
        cdef ComplexBall x
        x = ComplexBall.__new__(ComplexBall)
        x._parent = self._parent
        return x

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField
           sage: CBF = ComplexBallField()
           sage: CBF(1/3)
           [0.3333333333333333 +/- 7.04e-17]
           sage: CBF(0, 1/3)
           [0.3333333333333333 +/- 7.04e-17]*I
           sage: ComplexBallField()(1/3, 1/6)
           [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I

        TESTS::

           sage: CBF(1-I/2)
           1.000000000000000 - 0.5000000000000000*I
        """
        cdef arb_t real = acb_realref(self.value)
        cdef arb_t imag = acb_imagref(self.value)
        if arb_is_zero(imag):
            return self.real()._repr_()
        elif arb_is_zero(real):
            return "{}*I".format(self.imag()._repr_())
        elif arb_is_exact(imag) and arb_is_negative(imag):
            return "{} - {}*I".format(self.real()._repr_(),
                                        (-self.imag())._repr_())
        else:
            return "{} + {}*I".format(self.real()._repr_(),
                                        self.imag()._repr_())

    cpdef ComplexIntervalFieldElement _interval(self):
        """
        Return :class:`ComplexIntervalFieldElement` of the same value.

        OUTPUT:

        A :class:`ComplexIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: a = CBF(CIF(2, 2))
            sage: a._interval()
            2 + 2*I
        """
        cdef ComplexIntervalFieldElement target = ComplexIntervalField(prec(self))(0)
        acb_to_ComplexIntervalFieldElement(target, self.value)
        return target

    # Real and imaginary part, midpoint

    cpdef RealBall real(self):
        """
        Return the real part of this ball.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField
           sage: CBF = ComplexBallField()
           sage: a = CBF(1/3, 1/5)
           sage: a.real()
           [0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall r
        r = real_ball_field(self)(0)
        arb_set(r.value, &self.value.real)
        return r

    cpdef RealBall imag(self):
        """
        Return the imaginary part of this ball.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField
           sage: CBF = ComplexBallField()
           sage: a = CBF(1/3, 1/5)
           sage: a.imag()
           [0.2000000000000000 +/- 4.45e-17]
        """
        cdef RealBall r
        r = real_ball_field(self)(0)
        arb_set(r.value, &self.value.imag)
        return r

    def mid(self):
        """
        Return the floating-point complex number formed by the centers of the
        real and imaginary parts of this ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF(1/3, 1).mid()
            0.333333333333333 + 1.00000000000000*I
            sage: CBF(1/3, 1).mid().parent()
            Complex Field with 53 bits of precision
            sage: CBF('inf', 'nan').mid()
            +infinity - NaN*I
            sage: CBF('nan', 'inf').mid()
            NaN + +infinity*I
            sage: CBF('nan').mid()
            NaN
            sage: CBF('inf').mid()
            +infinity
            sage: CBF(0, 'inf').mid()
            +infinity*I
        """
        re, im = self.real().mid(), self.imag().mid()
        field = ComplexField(max(prec(self), re.prec(), im.prec()))
        return field(re, im)

    def squash(self):
        """
        Return an exact ball with the same midpoint as this ball.

        .. SEEALSO:: :meth:`mid`

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: mid = CBF(1/3, 1/10).squash()
            sage: mid
            [0.3333333333333333 +/- 1.49e-17] + [0.09999999999999999 +/- 1.68e-18]*I
            sage: mid.parent()
            Complex ball field with 53 bits precision
            sage: mid.is_exact()
            True
        """
        cdef ComplexBall res = self._new()
        arf_set(arb_midref(acb_realref(res.value)), arb_midref(acb_realref(self.value)))
        arf_set(arb_midref(acb_imagref(res.value)), arb_midref(acb_imagref(self.value)))
        mag_zero(arb_radref(acb_realref(res.value)))
        mag_zero(arb_radref(acb_imagref(res.value)))
        return res

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: CBF(0).is_zero()
            True
            sage: CBF(RIF(-0.5, 0.5)).is_zero()
            False
        """
        return acb_is_zero(self.value)

    def __nonzero__(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: bool(CBF(pi, 1/3))
            True
            sage: bool(CBF(RIF(-0.5, 0.5), 1/3))
            True
            sage: bool(CBF(1/3, RIF(-0.5, 0.5)))
            True
            sage: bool(CBF(RIF(-0.5, 0.5), RIF(-0.5, 0.5)))
            False
        """
        return acb_is_nonzero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: CBF(1).is_exact()
            True
            sage: CBF(1/3, 1/3).is_exact()
            False
        """
        return acb_is_exact(self.value)

    cpdef _richcmp_(left, Element right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.complex_ball_acb`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField
            sage: CBF = ComplexBallField()
            sage: a = CBF(1)
            sage: b = CBF(1)
            sage: a is b
            False
            sage: a == b
            True
            sage: a = CBF(1/3)
            sage: a.is_exact()
            False
            sage: b = CBF(1/3)
            sage: b.is_exact()
            False
            sage: a == b
            False
            sage: a = CBF(1, 2)
            sage: b = CBF(1, 2)
            sage: a is b
            False
            sage: a == b
            True

        TESTS:

        Balls whose intersection consists of one point::

            sage: a = CBF(RIF(1, 2), RIF(1, 2))
            sage: b = CBF(RIF(2, 4), RIF(2, 4))
            sage: a < b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a > b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a <= b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a >= b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a == b
            False
            sage: a != b
            False

        Balls with non-trivial intersection::

            sage: a = CBF(RIF(1, 4), RIF(1, 4))
            sage: a = CBF(RIF(2, 5), RIF(2, 5))
            sage: a == b
            False
            sage: a != b
            False

        One ball contained in another::

            sage: a = CBF(RIF(1, 4), RIF(1, 4))
            sage: b = CBF(RIF(2, 3), RIF(2, 3))
            sage: a == b
            False
            sage: a != b
            False

        Disjoint balls::

            sage: a = CBF(1/3, 1/3)
            sage: b = CBF(1/5, 1/5)
            sage: a == b
            False
            sage: a != b
            True

        Exact elements::

            sage: a = CBF(2, 2)
            sage: b = CBF(2, 2)
            sage: a.is_exact()
            True
            sage: b.is_exact()
            True
            sage: a == b
            True
            sage: a != b
            False
        """
        cdef ComplexBall lt, rt
        cdef acb_t difference

        lt = left
        rt = right

        if op == Py_EQ:
            return lt is rt or (
                acb_is_exact(lt.value) and acb_is_exact(rt.value)
                and acb_equal(lt.value, rt.value))

        if op == Py_NE:
            return not acb_overlaps(lt.value, rt.value)

        elif op == Py_GT or op == Py_GE or op == Py_LT or op == Py_LE:
            raise TypeError("No order is defined for ComplexBalls.")

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: -CBF(1/3 + I)
            [-0.3333333333333333 +/- 7.04e-17] - 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        acb_neg(res.value, self.value)
        return res

    def conjugate(self):
        """
        Return the complex conjugate of this ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF(-2 + I/3).conjugate()
            -2.000000000000000 + [-0.3333333333333333 +/- 7.04e-17]*I
        """
        cdef ComplexBall res = self._new()
        acb_conj(res.value, self.value)
        return res

    cpdef ModuleElement _add_(self, ModuleElement other):
        """
        Return the sum of two balls, rounded to the ambient field's precision.

        The resulting ball is guaranteed to contain the sums of any two points
        of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF(1) + CBF(I)
            1.000000000000000 + 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_add(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef ModuleElement _sub_(self, ModuleElement other):
        """
        Return the difference of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the differences of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF(1) - CBF(I)
            1.000000000000000 - 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sub(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __invert__(self):
        """
        Return the inverse of this ball.

        The result is guaranteed to contain the inverse of any point of the
        input ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: ~CBF(i/3)
            [-3.00000000000000 +/- 9.44e-16]*I
            sage: ~CBF(0)
            [+/- inf]
            sage: ~CBF(RIF(10,11))
            [0.1 +/- 9.53e-3]
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_inv(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _mul_(self, RingElement other):
        """
        Return the product of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the products of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: CBF(-2, 1)*CBF(1, 1/3)
            [-2.333333333333333 +/- 5.37e-16] + [0.333333333333333 +/- 4.82e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_mul(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef RingElement _div_(self, RingElement other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import CBF
            sage: from sage.rings.real_arb import RBF

            sage: CBF(-2, 1)/CBF(1, 1/3)
            [-1.50000000000000 +/- 1.27e-15] + [1.500000000000000 +/- 8.94e-16]*I

            sage: CBF(2+I)/CBF(0)
            [+/- inf] + [+/- inf]*I
            sage: CBF(1)/CBF(0)
            [+/- inf]
            sage: CBF(1)/CBF(RBF(0, 1.))
            nan
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_div(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

CBF = ComplexBallField()
