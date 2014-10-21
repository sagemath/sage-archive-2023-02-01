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

from sage.libs.arb.arb cimport *
from sage.libs.mpfr cimport mpfr_inits2, mpfr_clears, mpfr_t
from sage.libs.mpfi cimport mpfi_interv_fr
from sage.rings.real_mpfi cimport RealIntervalFieldElement, RealIntervalField

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
    cdef __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: Arb() # indirect doctest
            Memory initialized.
        """
        arb_init(self.value)
        print "Memory initialized."

    cdef __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.real_arb import Arb
            sage: a = Arb() # indirect doctest
            Memory initialized.
            sage: del a
            Memory deallocated.
        """
        arb_clear(self.value)
        print "Memory deallocated."

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
            sage: a.print_d()
            sage: b = Arb(1)
            Traceback (most recent call last):
            ...
            TypeError: value must be None or a RealIntervalFieldElement.
        """
        if value is None:
            pass
        elif isinstance(value, RealIntervalFieldElement):
            arb_set_interval_mpfr(self.value,
                                  value.left().value,
                                  value.right().value,
                                  value.parent().precision())
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
        """

        cdef mpfr_t left, right
        cdef int prec

        prec = arb_bits(self.value)
        mpfr_inits2(prec, left, right)
        arb_get_interval_mpfr(left, right, self.value)
        result = RealIntervalField(prec)(0)
        mpfi_interv_fr(result.value, left, right)
        mpfr_clears(left, right)

        return result

    cpdef psi(self):
        """
        Compute the digamma function with argument self.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`Arb`.

        EXAMPLES::

            sage: a = Arb(RIF(1))
            sage: a.psi().RealIntervalFieldElement()
        """

        cdef Arb result
        arb_digamma(result.value, self.value)
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
