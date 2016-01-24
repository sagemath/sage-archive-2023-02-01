# -*- coding: utf-8
r"""
Arbitrary precision complex balls using Arb

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Real balls using Arb <sage.rings.real_arb>`
    - :mod:`Complex interval field (using MPFI) <sage.rings.complex_interval_field>`
    - :mod:`Complex intervals (using MPFI) <sage.rings.complex_interval>`

Data Structure
==============

A :class:`ComplexBall` represents a complex number with error bounds. It wraps
an Arb object of type ``acb_t``, which  consists of a pair of real number balls
representing the real and imaginary part with separate error bounds. (See the
documentation of :mod:`sage.rings.real_arb` for more information.)

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

    In accordance with the semantics of Arb, identical :class:`ComplexBall`
    objects are understood to give permission for algebraic simplification.
    This assumption is made to improve performance. For example, setting ``z =
    x*x`` sets `z` to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

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

A ball is non-zero in the sense of usual comparison if and only if it does not
contain zero::

    sage: a = CBF(RIF(-0.5, 0.5))
    sage: a != 0
    False
    sage: b = CBF(1/3, 1/5)
    sage: b != 0
    True

However, ``bool(b)`` returns ``False`` for a ball ``b`` only if ``b`` is exactly
zero::

    sage: bool(a)
    True
    sage: bool(b)
    True
    sage: bool(CBF.zero())
    False

Coercion
========

Automatic coercions work as expected::

    sage: bpol = 1/3*CBF(i) + AA(sqrt(2)) + (polygen(RealBallField(20), 'x') + QQbar(i))
    sage: bpol
    x + [1.41421 +/- 5.09e-6] + [1.33333 +/- 3.97e-6]*I
    sage: bpol.parent()
    Univariate Polynomial Ring in x over Complex ball field with 20 bits precision
    sage: bpol/3
    ([0.333333 +/- 4.93e-7])*x + [0.47140 +/- 5.39e-6] + [0.44444 +/- 4.98e-6]*I

TESTS::

    sage: polygen(CBF, x)^3
    x^3

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
include "sage/ext/stdsage.pxi"

import sage.categories.fields

cimport sage.rings.integer
cimport sage.rings.rational

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.int cimport PyInt_AS_LONG
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

from sage.libs.mpfr cimport MPFR_RNDU
from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_hypgeom cimport *
from sage.libs.arb.arf cimport arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag, arf_set
from sage.libs.arb.mag cimport mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS, mag_is_inf, mag_is_finite, mag_zero
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi
from sage.rings.real_arb import RealBallField
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.rings.ring import Field
from sage.structure.element cimport Element, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation

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
    mpfi_to_arb(acb_realref(target), source.__re, precision)
    mpfi_to_arb(acb_imagref(target), source.__im, precision)

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

    arb_to_mpfi(target.__re, acb_realref(source), precision)
    arb_to_mpfi(target.__im, acb_imagref(source), precision)
    return 0

class ComplexBallField(UniqueRepresentation, Field):
    r"""
    An approximation of the field of complex numbers using pairs of mid-rad
    intervals.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

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

            sage: ComplexBallField(53) is ComplexBallField()
            True
        """
        return super(ComplexBallField, cls).__classcall__(cls, precision, category)

    def __init__(self, precision, category):
        r"""
        Initialize the complex ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: CBF = ComplexBallField()
            sage: CBF(1)
            1.000000000000000

        TESTS::

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
                base_ring=real_field,
                category=category or sage.categories.fields.Fields().Infinite())
        self._prec = precision
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.real_lazy import CLF
        self._populate_coercion_lists_([ZZ, QQ, real_field, CLF])

    def _real_field(self):
        """
        TESTS::

            sage: CBF._real_field()
            Real ball field with 53 bits precision
        """
        return self._base

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: ComplexBallField()
            Complex ball field with 53 bits precision
            sage: ComplexBallField(106)
            Complex ball field with 106 bits precision
        """
        return "Complex ball field with {} bits precision".format(self._prec)

    def construction(self):
        """
        Return the construction of a complex ball field as the algebraic
        closure of the real ball field with the same precision.

        EXAMPLES::

            sage: functor, base = CBF.construction()
            sage: functor, base
            (AlgebraicClosureFunctor, Real ball field with 53 bits precision)
            sage: functor(base) is CBF
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._base)

    def complex_field(self):
        """
        Return the complex ball field with the same precision, i.e. ``self``

        EXAMPLES::

            sage: CBF.complex_field() is CBF
            True
        """
        return ComplexBallField(self._prec)

    def ngens(self):
        r"""
        Return 1 as the only generator is the imaginary unit.

        EXAMPLE::

            sage: CBF.ngens()
            1
        """
        return 1

    def gen(self, i):
        r"""
        For i = 0, return the imaginary unit in this complex ball field.

        EXAMPLE::

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
        Return the tuple of generators of this complex ball field, i.e.
        ``(i,)``.

        EXAMPLE::

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

        EXAMPLES::

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

        .. SEEALSO::

            :meth:`sage.rings.real_arb.RealBallField._element_constructor_`

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
                x = self._base(x)
                return self.element_class(self, x)
            except (TypeError, ValueError):
                pass
            try:
                y = self._base(x.imag())
                x = self._base(x.real())
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
            x = self._base(x)
            y = self._base(y)
            return self.element_class(self, x, y)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: CBF.an_element() # indirect doctest
            [0.3333333333333333 +/- 1.49e-17] - [0.1666666666666667 +/- 4.26e-17]*I
        """
        return self(1.0/3, -1.0/6)

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: ComplexBallField().precision()
            53
        """
        return self._prec

    def is_exact(self):
        """
        Complex ball fields are not exact.

        EXAMPLES::

            sage: ComplexBallField().is_exact()
            False
        """
        return False

    def is_finite(self):
        """
        Complex ball fields are infinite.

        They already specify it via their category, but we currently need to
        re-implement this method due to the legacy implementation in
        :class:`sage.rings.ring.Ring`.

        EXAMPLES::

            sage: ComplexBallField().is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Complex ball fields have characteristic zero.

        EXAMPLES::

            sage: ComplexBallField().characteristic()
            0
        """
        return 0

    def some_elements(self):
        """
        Complex ball fields contain elements with exact, inexact, infinite, or
        undefined real and imaginary parts.

        EXAMPLES::

            sage: CBF.some_elements()
            [1.000000000000000,
             -0.5000000000000000*I,
             1.000000000000000 + [0.3333333333333333 +/- 1.49e-17]*I,
             [-0.3333333333333333 +/- 1.49e-17] + 0.2500000000000000*I,
             [-2.175556475109056e+181961467118333366510562 +/- 1.29e+181961467118333366510545],
             [+/- inf],
             [0.3333333333333333 +/- 1.49e-17] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I,
             nan,
             nan + nan*I,
             [+/- inf] + nan*I]
        """
        return [self(1), self(0, -1./2), self(1, 1./3), self(-1./3, 1./4),
                -self(1, 1)**(sage.rings.integer.Integer(2)**80),
                self('inf'), self(1./3, 'inf'), self('inf', 'inf'),
                self('nan'), self('nan', 'nan'), self('inf', 'nan')]

cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.
    """
    return (prec > 1000)

cdef inline long prec(ComplexBall ball):
    return ball._parent._prec

cdef inline Parent real_ball_field(ComplexBall ball):
    return ball._parent._base

cdef class ComplexBall(RingElement):
    """
    Hold one ``acb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: a = ComplexBallField()(1, 1)
        sage: a
        1.000000000000000 + 1.000000000000000*I
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: ComplexBallField(2)(0) # indirect doctest
            0
        """
        acb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

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

            sage: from sage.rings.complex_arb import ComplexBall
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

        RingElement.__init__(self, parent)

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

    def __hash__(self):
        """
        TESTS::

            sage: hash(CBF(1/3)) == hash(RBF(1/3))
            True
            sage: hash(CBF(1/3 + 2*i)) != hash(CBF(1/3 + i))
            True
        """
        if self.is_real():
            return hash(self.real())
        else:
            return (hash(self.real()) // 3) ^ hash(self.imag())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: CBF(1/3)
           [0.3333333333333333 +/- 7.04e-17]
           sage: CBF(0, 1/3)
           [0.3333333333333333 +/- 7.04e-17]*I
           sage: CBF(1/3, 1/6)
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

    def _is_atomic(self):
        r"""
        Declare that complex balls print atomically in some cases.

        TESTS::

            sage: CBF(-1/3)._is_atomic()
            True

        This method should in principle ensure that ``CBF['x']([1, -1/3])``
        is printed as::

            sage: CBF['x']([1, -1/3]) # todo - not tested
            [-0.3333333333333333 +/- 7.04e-17]*x + 1.000000000000000

        However, this facility is not really used in Sage at this point, and we
        still get::

            sage: CBF['x']([1, -1/3])
            ([-0.3333333333333333 +/- 7.04e-17])*x + 1.000000000000000
        """
        return self.is_real() or self.real().is_zero()

    # Conversions

    cpdef ComplexIntervalFieldElement _complex_mpfi_(self, parent):
        """
        Return :class:`ComplexIntervalFieldElement` of the same value.

        EXAMPLES::

            sage: CIF(CBF(1/3, 1/3)) # indirect doctest
            0.3333333333333333? + 0.3333333333333333?*I
        """
        cdef ComplexIntervalFieldElement res = parent.zero()
        res = res._new() # FIXME after modernizing CIF
        acb_to_ComplexIntervalFieldElement(res, self.value)
        return res

    def _integer_(self, _):
        """
        Check that this ball contains a single integer and return that integer.

        EXAMPLES::

            sage: ZZ(CBF(-42, RBF(.1, rad=.2))) # indirect doctest
            -42
            sage: ZZ(CBF(i))
            Traceback (most recent call last):
            ...
            ValueError: 1.000000000000000*I does not contain a unique integer
        """
        cdef sage.rings.integer.Integer res
        cdef fmpz_t tmp
        fmpz_init(tmp)
        try:
            if acb_get_unique_fmpz(tmp, self.value):
                res = sage.rings.integer.Integer.__new__(sage.rings.integer.Integer)
                fmpz_get_mpz(res.value, tmp)
            else:
                raise ValueError("{} does not contain a unique integer".format(self))
        finally:
            fmpz_clear(tmp)
        return res

    def _rational_(self):
        """
        Check that this ball contains a single rational number and return that
        number.

        EXAMPLES::

            sage: QQ(CBF(12345/2^5))
            12345/32
            sage: QQ(CBF(i))
            Traceback (most recent call last):
            ...
            ValueError: 1.000000000000000*I does not contain a unique rational number
        """
        if acb_is_real(self.value) and acb_is_exact(self.value):
            return self.real().mid().exact_rational()
        else:
            raise ValueError("{} does not contain a unique rational number".format(self))

    def _complex_mpfr_field_(self, parent):
        r"""
        Convert this complex ball to a complex number.

        INPUT:

        - ``parent`` - :class:`~sage.rings.complex_field.ComplexField_class`,
          target parent.

        EXAMPLES::

            sage: CC(CBF(1/3, 1/3))
            0.333333333333333 + 0.333333333333333*I
            sage: ComplexField(100)(CBF(1/3, 1/3))
            0.33333333333333331482961625625 + 0.33333333333333331482961625625*I
        """
        real_field = parent._base
        return parent(real_field(self.real()), real_field(self.imag()))

    def _real_mpfi_(self, parent):
        r"""
        Try to convert this complex ball to a real interval.

        Fail if the imaginary part is not exactly zero.

        INPUT:

        - ``parent`` - :class:`~sage.rings.real_mpfi.RealIntervalField_class`,
          target parent.

        EXAMPLES::

            sage: RIF(CBF(RBF(1/3, rad=1e-5)))
            0.3334?
            sage: RIF(CBF(RBF(1/3, rad=1e-5), 1e-10))
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
        """
        if acb_is_real(self.value):
            return parent(self.real())
        else:
            raise ValueError("nonzero imaginary part")

    def _mpfr_(self, parent):
        r"""
        Try to convert this complex ball to a real number.

        Fail if the imaginary part is not exactly zero.

        INPUT:

        - ``parent`` - :class:`~sage.rings.real_mpfr.RealField_class`,
          target parent.

        EXAMPLES::

            sage: RR(CBF(1/3))
            0.333333333333333
            sage: RR(CBF(1, 1/3) - CBF(0, 1/3))
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
        """
        if acb_is_real(self.value):
            return parent(self.real())
        else:
            raise ValueError("nonzero imaginary part")

    # Real and imaginary part, midpoint, radius

    cpdef RealBall real(self):
        """
        Return the real part of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_arb.RealBall`.

        EXAMPLES::

           sage: CBF = ComplexBallField()
           sage: a = CBF(1/3, 1/5)
           sage: a.real()
           [0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        arb_set(r.value, acb_realref(self.value))
        return r

    cpdef RealBall imag(self):
        """
        Return the imaginary part of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_arb.RealBall`.

        EXAMPLES::

           sage: CBF = ComplexBallField()
           sage: a = CBF(1/3, 1/5)
           sage: a.imag()
           [0.2000000000000000 +/- 4.45e-17]
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        arb_set(r.value, acb_imagref(self.value))
        return r

    def __abs__(self):
        """
        Return the absolute value of this complex ball.

        EXAMPLES::

            sage: CBF(1 + i).abs() # indirect doctest
            [1.414213562373095 +/- 2.99e-16]
            sage: abs(CBF(i))
            1.000000000000000

            sage: CBF(1 + i).abs().parent()
            Real ball field with 53 bits precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_abs(r.value, self.value, prec(self))
        return r

    def below_abs(self, test_zero=False):
        """
        Return a lower bound for the absolute value of this complex ball.

        INPUT:

        - ``test_zero`` (boolean, default ``False``) -- if ``True``,
          make sure that the returned lower bound is positive, raising
          an error if the ball contains zero.

        OUTPUT:

        A ball with zero radius

        EXAMPLES::

            sage: b = ComplexBallField(8)(1+i).below_abs()
            sage: b
            [1.4 +/- 0.0141]
            sage: b.is_exact()
            True
            sage: QQ(b)*128
            181
            sage: (CBF(1/3) - 1/3).below_abs()
            0
            sage: (CBF(1/3) - 1/3).below_abs(test_zero=True)
            Traceback (most recent call last):
            ...
            ValueError: ball contains zero

        .. SEEALSO:: :meth:`above_abs`
        """
        cdef RealBall res = RealBall(real_ball_field(self))
        acb_get_abs_lbound_arf(arb_midref(res.value), self.value, prec(self))
        if test_zero and arb_contains_zero(res.value):
            assert acb_contains_zero(self.value)
            raise ValueError("ball contains zero")
        return res

    def above_abs(self):
        """
        Return an upper bound for the absolute value of this complex ball.

        OUTPUT:

        A ball with zero radius

        EXAMPLES::

            sage: b = ComplexBallField(8)(1+i).above_abs()
            sage: b
            [1.4 +/- 0.0219]
            sage: b.is_exact()
            True
            sage: QQ(b)*128
            182

        .. SEEALSO:: :meth:`below_abs`
        """
        cdef RealBall res = RealBall(real_ball_field(self))
        acb_get_abs_ubound_arf(arb_midref(res.value), self.value, prec(self))
        return res

    def arg(self):
        """
        Return the argument of this complex ball.

        EXAMPLES::

            sage: CBF(1 + i).arg()
            [0.785398163397448 +/- 3.91e-16]
            sage: CBF(-1).arg()
            [3.141592653589793 +/- 5.61e-16]
            sage: CBF(-1).arg().parent()
            Real ball field with 53 bits precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_arg(r.value, self.value, prec(self))
        return r

    def mid(self):
        """
        Return the midpoint of this ball.

        OUTPUT:

        :class:`~sage.rings.complex_number.ComplexNumber`, floating-point
        complex number formed by the centers of the real and imaginary parts of
        this ball.

        EXAMPLES::

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

        .. SEEALSO:: :meth:`squash`
        """
        re, im = self.real().mid(), self.imag().mid()
        field = sage.rings.complex_field.ComplexField(
                max(prec(self), re.prec(), im.prec()))
        return field(re, im)

    def squash(self):
        """
        Return an exact ball with the same midpoint as this ball.

        OUTPUT:

        A :class:`ComplexBall`.

        EXAMPLES::

            sage: mid = CBF(1/3, 1/10).squash()
            sage: mid
            [0.3333333333333333 +/- 1.49e-17] + [0.09999999999999999 +/- 1.68e-18]*I
            sage: mid.parent()
            Complex ball field with 53 bits precision
            sage: mid.is_exact()
            True

        .. SEEALSO:: :meth:`mid`
        """
        cdef ComplexBall res = self._new()
        arf_set(arb_midref(acb_realref(res.value)), arb_midref(acb_realref(self.value)))
        arf_set(arb_midref(acb_imagref(res.value)), arb_midref(acb_imagref(self.value)))
        mag_zero(arb_radref(acb_realref(res.value)))
        mag_zero(arb_radref(acb_imagref(res.value)))
        return res

    def rad(self):
        """
        Return an upper bound for the error radius of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfr.RealNumber` of the same precision as
        the radii of real balls.

        .. WARNING::

            Unlike a :class:`~sage.rings.real_arb.RealBall`,
            a :class:`ComplexBall` is *not* defined
            by its midpoint and radius. (Instances of :class:`ComplexBall` are
            actually rectangles, not balls.)

        EXAMPLES::

            sage: CBF(1 + i).rad()
            0.00000000
            sage: CBF(i/3).rad()
            1.1102230e-16
            sage: CBF(i/3).rad().parent()
            Real Field with 30 bits of precision

        TESTS::

            sage: (CBF(0, 1/3) << (2^64)).rad()
            Traceback (most recent call last):
            ...
            RuntimeError: unable to convert the radius to MPFR (exponent out of range?)
        """
        # Should we return a real number with rounding towards +∞ (or away from
        # zero if/when implemented)?
        cdef RealField_class rad_field = RealField(MAG_BITS)
        cdef RealNumber rad = RealNumber(rad_field, None)
        cdef arf_t tmp
        arf_init(tmp)
        acb_get_rad_ubound_arf(tmp, self.value, MAG_BITS)
        sig_str("unable to convert the radius to MPFR (exponent out of range?)")
        if arf_get_mpfr(rad.value, tmp, MPFR_RNDU):
            sig_error()
        sig_off()
        arf_clear(tmp)
        return rad

    # Should we implement rad_as_ball? If we do, should it return an enclosure
    # of the radius (which radius?), or an upper bound?

    # Precision and accuracy

    def round(self):
        """
        Return a copy of this ball rounded to the precision of the parent.

        EXAMPLES:

        It is possible to create balls whose midpoint is more precise that
        their parent's nominal precision (see :mod:`~sage.rings.real_arb` for
        more information)::

            sage: b = CBF(exp(I*pi/3).n(100))
            sage: b.mid()
            0.50000000000000000000000000000 + 0.86602540378443864676372317075*I

        The ``round()`` method rounds such a ball to its parent's precision::

            sage: b.round().mid()
            0.500000000000000 + 0.866025403784439*I

        .. SEEALSO:: :meth:`trim`
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_set_round(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def accuracy(self):
        """
        Return the effective relative accuracy of this ball measured in bits.

        This is computed as if calling
        :meth:`~sage.rings.real_arb.RealBall.accuracy()`
        on the real ball whose midpoint is the larger out of the real and
        imaginary midpoints of this complex ball, and whose radius is the
        larger out of the real and imaginary radii of this complex ball.

        EXAMPLES::

            sage: CBF(exp(I*pi/3)).accuracy()
            51
            sage: CBF(I/2).accuracy() == CBF.base().maximal_accuracy()
            True
            sage: CBF('nan', 'inf').accuracy() == -CBF.base().maximal_accuracy()
            True

        .. SEEALSO::

            :meth:`~sage.rings.real_arb.RealBallField.maximal_accuracy`
        """
        return acb_rel_accuracy_bits(self.value)

    def trim(self):
        """
        Return a trimmed copy of this ball.

        Return a copy of this ball with both the real and imaginary parts
        trimmed (see :meth:`~sage.rings.real_arb.RealBall.trim()`).

        EXAMPLES::

            sage: b = CBF(1/3, RBF(1/3, rad=.01))
            sage: b.mid()
            0.333333333333333 + 0.333333333333333*I
            sage: b.trim().mid()
            0.333333333333333 + 0.333333015441895*I

        .. SEEALSO:: :meth:`round`
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_trim(res.value, self.value)
        if _do_sig(prec(self)): sig_off()
        return res

    def add_error(self, ampl):
        """
        Increase the radii of the real and imaginary parts by (an upper bound
        on) ``ampl``.

        If ``ampl`` is negative, the radii remain unchanged.

        INPUT:

        - ``ampl`` - A **real** ball (or an object that can be coerced to a
          real ball).

        OUTPUT:

        A new complex ball.

        EXAMPLES::

            sage: CBF(1+i).add_error(10^-16)
            [1.000000000000000 +/- 1.01e-16] + [1.000000000000000 +/- 1.01e-16]*I
        """
        return ComplexBall(self._parent, self.real().add_error(ampl), self.imag().add_error(ampl))

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: CBF = ComplexBallField()
            sage: CBF(0).is_zero()
            True
            sage: CBF(RIF(-0.5, 0.5)).is_zero()
            False

        .. SEEALSO:: :meth:`is_nonzero`
        """
        return acb_is_zero(self.value)

    def is_nonzero(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        .. NOTE::

            This method is not the negation of :meth:`is_zero`: it only
            returns ``True`` if zero is known not to be contained in the ball.

            Use ``bool(b)`` (or, equivalently, ``not b.is_zero()``) to check if
            a ball ``b`` **may** represent a nonzero number (for instance, to
            determine the “degree” of a polynomial with ball coefficients).

        EXAMPLES::

            sage: CBF = ComplexBallField()
            sage: CBF(pi, 1/3).is_nonzero()
            True
            sage: CBF(RIF(-0.5, 0.5), 1/3).is_nonzero()
            True
            sage: CBF(1/3, RIF(-0.5, 0.5)).is_nonzero()
            True
            sage: CBF(RIF(-0.5, 0.5), RIF(-0.5, 0.5)).is_nonzero()
            False

        .. SEEALSO:: :meth:`is_zero`
        """
        return (arb_is_nonzero(acb_realref(self.value))
                or arb_is_nonzero(acb_imagref(self.value)))

    def __nonzero__(self):
        """
        Return ``True`` iff this complex ball is not the zero ball, i.e. if the
        midpoint and radius of its real and imaginary parts are not all zero.

        This is the preferred way, for instance, to determine the “degree” of a
        polynomial with ball coefficients.

        .. WARNING::

            A “nonzero” ball in the sense of this method may represent the
            value zero. Use :meth:`is_nonzero` to check that a complex number
            represented by a ``ComplexBall`` object is known to be nonzero.

        EXAMPLES::

            sage: bool(CBF(0)) # indirect doctest
            False
            sage: bool(CBF(i))
            True
            sage: bool(CBF(RIF(-0.5, 0.5)))
            True
        """
        return not acb_is_zero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: CBF = ComplexBallField()
            sage: CBF(1).is_exact()
            True
            sage: CBF(1/3, 1/3).is_exact()
            False
        """
        return acb_is_exact(self.value)

    def is_real(self):
        """
        Return ``True`` iff the imaginary part of this ball is exactly zero.

        EXAMPLES::

            sage: CBF(1/3, 0).is_real()
            True
            sage: (CBF(i/3) - CBF(1, 1/3)).is_real()
            False
            sage: CBF('inf').is_real()
            True
        """
        return acb_is_real(self.value)

    cpdef _richcmp_(left, Element right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.complex_arb`.

        EXAMPLES::

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

    def identical(self, ComplexBall other):
        """
        Return whether ``self`` and ``other`` represent the same ball.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`.

        OUTPUT:

        Return True iff ``self`` and ``other`` are equal as sets, i.e. if their
        real and imaginary parts each have the same midpoint and radius.

        Note that this is not the same thing as testing whether both ``self``
        and ``other`` certainly represent the complex real number, unless
        either ``self`` or ``other`` is exact (and neither contains NaN). To
        test whether both operands might represent the same mathematical
        quantity, use :meth:`overlaps` or ``in``, depending on the
        circumstance.

        EXAMPLES::

            sage: CBF(1, 1/3).identical(1 + CBF(0, 1)/3)
            True
            sage: CBF(1, 1).identical(1 + CBF(0, 1/3)*3)
            False
        """
        return acb_equal(self.value, other.value)

    def overlaps(self, ComplexBall other):
        """
        Return True iff ``self`` and ``other`` have some point in common.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`.

        EXAMPLES::

            sage: CBF(1, 1).overlaps(1 + CBF(0, 1/3)*3)
            True
            sage: CBF(1, 1).overlaps(CBF(1, 'nan'))
            True
            sage: CBF(1, 1).overlaps(CBF(0, 'nan'))
            False
        """
        return acb_overlaps(self.value, other.value)

    def contains_exact(self, other):
        """
        Return ``True`` *iff* ``other`` is contained in ``self``.

        Use ``other in self`` for a test that works for a wider range of inputs
        but may return false negatives.

        INPUT:

        - ``other`` -- :class:`ComplexBall`,
          :class:`~sage.rings.integer.Integer`,
          or :class:`~sage.rings.rational.Rational`

        EXAMPLES::

            sage: CBF(RealBallField(100)(1/3), 0).contains_exact(1/3)
            True
            sage: CBF(1).contains_exact(1)
            True
            sage: CBF(1).contains_exact(CBF(1))
            True

            sage: CBF(sqrt(2)).contains_exact(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError: unsupported type: <type 'sage.symbolic.expression.Expression'>
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        if _do_sig(prec(self)): sig_on()
        try:
            if isinstance(other, ComplexBall):
                res = acb_contains(self.value, (<ComplexBall> other).value)
            elif isinstance(other, sage.rings.integer.Integer):
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> other).value)
                res = acb_contains_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
            elif isinstance(other, sage.rings.rational.Rational):
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> other).value)
                res = acb_contains_fmpq(self.value, tmpq)
                fmpq_clear(tmpq)
            else:
                raise TypeError("unsupported type: " + str(type(other)))
        finally:
            if _do_sig(prec(self)): sig_off()
        return res

    def __contains__(self, other):
        """
        Return True if ``other`` can be verified to be contained in ``self``.

        Depending on the type of ``other``, the test may use interval
        arithmetic with a precision determined by the parent of ``self`` and
        may return false negatives.

        EXAMPLES::

            sage: 1/3*i in CBF(0, 1/3)
            True

        A false negative::

            sage: RLF(1/3) in CBF(RealBallField(100)(1/3), 0)
            False

        .. SEEALSO:: :meth:`contains_exact`
        """
        if not isinstance(other, (
                ComplexBall,
                sage.rings.integer.Integer,
                sage.rings.rational.Rational)):
            other = self._parent(other)
        return self.contains_exact(other)

    def contains_zero(self):
        """
        Return ``True`` iff this ball contains zero.

        EXAMPLES::

            sage: CBF(0).contains_zero()
            True
            sage: CBF(RIF(-1,1)).contains_zero()
            True
            sage: CBF(i).contains_zero()
            False
        """
        return acb_contains_zero(self.value)

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

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

            sage: CBF(-2, 1)*CBF(1, 1/3)
            [-2.333333333333333 +/- 5.37e-16] + [0.333333333333333 +/- 4.82e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_mul(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __lshift__(val, shift):
        r"""
        If ``val`` is a ``ComplexBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the left
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer, may be negative.

        EXAMPLES::

            sage: CBF(i/3) << 2
            [1.333333333333333 +/- 4.82e-16]*I
            sage: CBF(i) << -2
            0.2500000000000000*I

        TESTS::

            sage: CBF(i) << (2^65)
            [3.636549880934858e+11106046577046714264 +/- 1.91e+11106046577046714248]*I
            sage: 'a' << CBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for <<: 'str' and
            'ComplexBall'
            sage: CBF(1) << 1/2
            Traceback (most recent call last):
            ...
            TypeError: shift should be an integer
        """
        cdef fmpz_t tmpz
        # the ComplexBall might be shift, not val
        if not isinstance(val, ComplexBall):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))
        cdef ComplexBall self = val
        cdef ComplexBall res = self._new()
        if isinstance(shift, int):
             acb_mul_2exp_si(res.value, self.value, PyInt_AS_LONG(shift))
        elif isinstance(shift, sage.rings.integer.Integer):
            sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> shift).value)
            acb_mul_2exp_fmpz(res.value, self.value, tmpz)
            fmpz_clear(tmpz)
            sig_off()
        else:
            raise TypeError("shift should be an integer")
        return res

    def __rshift__(val, shift):
        r"""
        If ``val`` is a ``ComplexBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the right
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer, may be negative.

        EXAMPLES::

            sage: CBF(1+I) >> 2
            0.2500000000000000 + 0.2500000000000000*I

        TESTS::

            sage: 'a' >> CBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for >>: 'str' and
            'ComplexBall'
        """
        # the ComplexBall might be shift, not val
        if isinstance(val, ComplexBall):
            return val << (-shift)
        else:
            raise TypeError("unsupported operand type(s) for >>: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))

    cpdef RingElement _div_(self, RingElement other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

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
