r"""
Arbitrary Precision Real Intervals using Arb

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
from sage.libs.arb.arb cimport *
from sage.libs.mpfr cimport mpfr_inits2, mpfr_clears, mpfr_t
from sage.rings.real_mpfi import RealIntervalField

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const unsigned long precision):
    """
    Convert an ``mpfi`` to an ``arb``.

    INPUT:

    - ``target`` -- an ``arb_t``.

    - ``source`` -- an ``mpfi_t``.

    - ``precision`` -- a positive integer.

    OUTPUT:

    None.
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

cdef void arb_to_mpfi(mpfi_t target, arb_t source, const unsigned long precision):
    """
    Convert an ``arb`` to an ``mpfi``.

    INPUT:

    - ``target`` -- an ``mpfi_t``.

    - ``source`` -- an ``arb_t``.

    - ``precision`` -- a positive integer.

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




cdef class Arb(SageObject):
    """
    Hold one ``arb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    INPUT:

    None.

    OUTPUT:

    None.

    EXAMPLES::

        sage: from sage.rings.real_arb import Arb # optional - arb
        sage: a = Arb(RIF(1))                     # optional - arb
        sage: b = a.psi()                         # optional - arb
        sage: b.RealIntervalFieldElement()        # optional - arb
        -0.577215664901533?
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb # optional - arb
            sage: Arb(precision=2) # optional - arb; indirect doctest
            <type 'sage.rings.real_arb.Arb'>
        """
        arb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb # optional - arb
            sage: a = Arb(precision=2) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        arb_clear(self.value)

    def __init__(self, value=None, unsigned long precision=0):
        """
        Initialize Arb using value.

        INPUT:

        - ``value`` -- (default: ``None``) ``None`` or a
          :class:`RealIntervalFieldElement`.

        - ``precision`` -- (default: ``0``) a non-negative
          integer. Must be given unless ``value`` is not ``None``.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb # optional - arb
            sage: a = Arb(RIF(0, 1))                  # optional - arb
            sage: b = Arb(1)                          # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: value must be None or a RealIntervalFieldElement.
            sage: c = Arb() # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: precision must be given.
        """
        cdef RealIntervalFieldElement element
        cdef int prec

        if value is None:
            if precision > 0:
                self.precision = precision
            else:
                raise TypeError("precision must be given.")
        elif isinstance(value, RealIntervalFieldElement):
            element = <RealIntervalFieldElement> value
            self.precision = value.parent().precision()
            mpfi_to_arb(self.value, element.value, self.precision)

        else:
            raise TypeError("value must be None or a "
                            "RealIntervalFieldElement.")

    cpdef RealIntervalFieldElement RealIntervalFieldElement(self):
        """
        Return :class:`RealIntervalFieldElement` of the same value.

        INPUT:

        None.

        OUTPUT:

        A :class:`RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb # optional - arb
            sage: a = Arb(RIF(2))                     # optional - arb
            sage: a.RealIntervalFieldElement()        # optional - arb
            2
        """

        cdef RealIntervalFieldElement result

        result = RealIntervalField(self.precision)(0)
        arb_to_mpfi(result.value, self.value, self.precision)

        return result

    cpdef Arb psi(self):
        """
        Compute the digamma function with argument self.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`Arb`.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb # optional - arb
            sage: a = Arb(RIF(1))                     # optional - arb
            sage: a.psi().RealIntervalFieldElement()  # optional - arb
            -0.577215664901533?
        """

        cdef Arb result
        cdef int prec

        result = Arb(precision=self.precision)

        arb_digamma(result.value, self.value, self.precision)
        return result
