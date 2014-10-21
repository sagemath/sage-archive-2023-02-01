r"""
Arbitrary Precision Real Intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-21): Initial version.

This is a rudimentary binding to the `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************
"""
include 'mpfi.pxi'

from sage.libs.arb.arb cimport *
from sage.libs.mpfr cimport mpfr_inits2, mpfr_clears, mpfr_t
from sage.rings.real_mpfi cimport RealIntervalFieldElement
from sage.rings.real_mpfi import RealIntervalField
from sage.structure.sage_object cimport SageObject

cdef class Arb(SageObject):
    """
    Hold one ``arb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    INPUT:

    None.

    OUTPUT:

    None.

    EXAMPLES::
    """
    cdef arb_t value
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: Arb() # indirect doctest
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

            sage: from sage.rings.real_arb import Arb
            sage: a = Arb() # indirect doctest
            sage: del a
        """
        arb_clear(self.value)

    def __init__(self, value=None):
        """
        Initialize Arb using value.

        INPUT:

        - `value` -- (default: ``None``) ``None`` or a
          :class:`RealIntervalFieldElement`.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: a = Arb(RIF(0, 1))
            sage: b = Arb(1)
            Traceback (most recent call last):
            ...
            TypeError: value must be None or a RealIntervalFieldElement.
        """
        cdef RealIntervalFieldElement element
        cdef mpfr_t left
        cdef mpfr_t right
        cdef int prec

        if value is None:
            pass
        elif isinstance(value, RealIntervalFieldElement):
            element = <RealIntervalFieldElement> value
            prec = value.parent().precision()
            mpfr_init2(left, prec)
            mpfr_init2(right, prec)
            mpfi_get_left(left, element.value)
            mpfi_get_right(right, element.value)
            arb_set_interval_mpfr(self.value,
                                  left,
                                  right,
                                  prec)
            mpfr_clear(left)
            mpfr_clear(right)

        else:
            raise TypeError("value must be None or a "
                            "RealIntervalFieldElement.")

    cpdef RealIntervalFieldElement(self):
        """
        Return :class:`RealIntervalFieldElement` of the same value.

        INPUT:

        None.

        OUTPUT:

        A :class:`RealIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: a = Arb(RIF(2))
            sage: a.RealIntervalFieldElement()
            2
        """

        cdef mpfr_t left
        cdef mpfr_t right
        cdef int prec
        cdef RealIntervalFieldElement result

        prec = max(arb_bits(self.value), 2)
        mpfr_init2(left, prec)
        mpfr_init2(right, prec)
        arb_get_interval_mpfr(left, right, self.value)
        result = RealIntervalField(prec)(0)
        mpfi_interv_fr(result.value, left, right)
        mpfr_clear(left)
        mpfr_clear(right)

        return result

    cpdef psi(self):
        """
        Compute the digamma function with argument self.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`Arb`.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: a = Arb(RIF(1))
            sage: a.psi().RealIntervalFieldElement()
            -1.?
        """

        cdef Arb result
        cdef int prec

        result = Arb()

        prec = arb_bits(self.value)
        arb_digamma(result.value, self.value, prec)
        return result



def _test_arb_():
    """
    EXAMPLES::

        sage: from sage.rings.real_arb import _test_arb_
        sage: _test_arb_()
        (948239929664858513 * 2^-59) +/- (692955552 * 2^-58)
    """

    cdef arb_t x
    cdef arb_t y
    arb_init(x)
    arb_init(y)
    arb_set_ui(x, 2)
    arb_zeta(y, x, 53)
    arb_print(y)
    arb_clear(x)
    arb_clear(y)
