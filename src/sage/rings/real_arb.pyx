r"""
Arbitrary precision real intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-21): Initial version.

This is a rudimentary binding to the optional `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

You may have to run ``sage -i arb`` to use the arb library.
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

import sage.categories.sets_cat
from sage.libs.arb.arb cimport *
from sage.libs.arb.arf cimport arf_t, arf_get_mpfr
from sage.libs.arb.fmpr cimport fmpr_t, fmpr_init, fmpr_clear, fmpr_get_mpfr
from sage.libs.arb.mag cimport mag_t, mag_get_fmpr
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
        self.precision = precision

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
        return "Real ball field with {} bits precision".format(self.precision)

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
        return self._element_constructor_(self)


cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: from sage.rings.real_arb import RealBallField # optional - arb
        sage: _ = RealBallField()(1).psi() # optional - arb; indirect doctest
        sage: _ = RealBallField(1500)(1).psi() # optional - arb
    """
    return (prec > 1000)

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
          :class:`RealIntervalFieldElement`.

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
                x = RealIntervalField(self._parent.precision)(x)
            except TypeError:
                raise TypeError("unable to convert to a RealIntervalFieldElement")

        mpfi_to_arb(self.value,
                    (<RealIntervalFieldElement> x).value,
                    self._parent.precision)

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

        c_result = arb_get_str(self.value,
                               (self._parent.precision * 31) // 100,
                               0)
        try:
            py_string = c_result
        finally:
            flint_free(c_result)

        return py_string

    cpdef RealIntervalFieldElement _interval(self):
        """
        Return a :mod:`real interval <sage.rings.real_mpfr>` containing this ball.

        OUTPUT:

        A :class:`RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: a = RealBallField()(RIF(2))                     # optional - arb
            sage: a._interval()        # optional - arb
            2
        """

        cdef RealIntervalFieldElement result

        result = RealIntervalField(self._parent.precision)(0)
        arb_to_mpfi(result.value,
                    self.value,
                    self._parent.precision)

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

        if _do_sig(self._parent.precision): sig_on()
        arb_digamma(result.value, self.value, self._parent.precision)
        if _do_sig(self._parent.precision): sig_off()

        return result

