r"""
Arbitrary precision real intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-21): Initial version.

This is a rudimentary binding to the optional `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

You may have to run ``sage -i arb`` to use the arb library.

.. WARNING::

    Identical :class:`RealBall` objects are understood to give
    permission for algebraic simplification. This assumption is made
    to improve performance.  For example, setting ``z = x*x`` sets `z`
    to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Comparison
==========

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

import sage.categories.sets_cat
from sage.libs.arb.arb cimport *
from sage.libs.arb.arf cimport arf_t, arf_get_mpfr
from sage.libs.arb.mag cimport mag_t
from sage.libs.flint.flint cimport flint_free
from sage.libs.mpfi cimport mpfi_get_left, mpfi_get_right, mpfi_interv_fr
from sage.libs.mpfr cimport mpfr_t, mpfr_init2, mpfr_clear, GMP_RNDN
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.structure.unique_representation import UniqueRepresentation

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

cdef void arb_to_mpfi(mpfi_t target, arb_t source, const long precision):
    """
    Convert an Arb ball to an MPFI interval.

    INPUT:

    - ``target`` -- an ``mpfi_t``.

    - ``source`` -- an ``arb_t``.

    - ``precision`` -- an integer `\ge 2`.

    OUTPUT:

    None.
    """
    cdef mpfr_t left
    cdef mpfr_t right

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    arb_get_interval_mpfr(left, right, source)
    mpfi_interv_fr(target, left, right)

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

    TESTS::

        sage: RealBallField(0) # optional - arb
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
        sage: RealBallField(1) # optional - arb
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
    """
    Element = RealBall

    @staticmethod
    def __classcall__(cls, long precision=53):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField(53) is RealBallField() # optional - arb
            True
        """
        return super(RealBallField, cls).__classcall__(cls, precision)

    def __init__(self, precision):
        r"""
        Initialize the real ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RBF = RealBallField() # optional - arb
            sage: RBF(1) # optional - arb
            1.000000000000000
        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        super(RealBallField, self).__init__(categories=[sage.categories.sets_cat.Sets])
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

    def _coerce_map_from_(self, S):
        r"""
        Currently, there is no coercion.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()._coerce_map_from_(RIF) # optional - arb
            False
            sage: RealBallField()._coerce_map_from_(SR) # optional - arb
            False
        """
        return False

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construct a :class:`RealBall`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional: arb
            sage: RBF = RealBallField() # optional: arb
            sage: RBF(RIF(1)) # optional: arb; indirect doctest
            1.000000000000000
        """
        return self.element_class(self, *args, **kwds)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional: arb
            sage: RBF = RealBallField() # optional: arb
            sage: RBF._an_element_() # optional: arb; indirect doctest
            [0.3333333333333333 +/- 1.49e-17]
        """
        return self(1.0/3)

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

    def is_finite(self):
        """
        Real ball fields are infinite.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField().is_finite() # optional - arb
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

cdef class RealBall(Element):
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

    def __init__(self, parent, x=None):
        """
        Initialize the :class:`RealBall` using ``x``.

        INPUT:

        - ``parent`` -- a :class:`RealBallField`.

        - ``x`` -- (default: ``None``) ``None`` or a
          :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: RealBallField()(RIF(0, 1))                  # optional - arb; indirect doctest
            [+/- 1.01]
            sage: RealBallField()(1)                          # optional - arb
            1.000000000000000
            sage: RealBallField()(x)                          # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: unable to convert to a RealIntervalFieldElement
        """
        super(RealBall, self).__init__(parent)

        if x is None:
            return

        if not isinstance(x, RealIntervalFieldElement):
            try:
                x = RealIntervalField(prec(self))(x)
            except TypeError:
                raise TypeError("unable to convert to a RealIntervalFieldElement")

        mpfi_to_arb(self.value,
                    (<RealIntervalFieldElement> x).value,
                    prec(self))

    cdef RealBall _new(self):
        """
        Return a new real ball element with the same parent as ``self``.
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

    cpdef RealIntervalFieldElement _interval(self):
        """
        Return a :mod:`real interval <sage.rings.real_mpfr>` containing this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(2))                     # optional - arb
            sage: a._interval()        # optional - arb
            2
        """

        cdef RealIntervalFieldElement result

        result = RealIntervalField(prec(self))(0)
        arb_to_mpfi(result.value, self.value, prec(self))

        return result

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
        """
        cdef RealBall lt, rt
        cdef arb_t difference

        lt = left
        rt = right

        if op == Py_EQ:
            return (lt is rt) or (
                arb_is_exact(lt.value) and arb_is_exact(rt.value)
                and arb_equal(lt.value, rt.value))

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

    cpdef RealBall psi(self):
        """
        Compute the digamma function with argument self.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(1))                     # optional - arb
            sage: a.psi()  # optional - arb
            [-0.577215664901533 +/- 3.85e-16]
        """

        cdef RealBall result

        result = self._new()

        if _do_sig(prec(self)): sig_on()
        arb_digamma(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()

        return result

