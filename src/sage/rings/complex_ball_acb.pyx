r"""
Arbitrary Precision Complex Intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

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
from sage.libs.arb.acb cimport *
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi

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
    cdef unsigned long precision
    precision = source.parent().precision()
    mpfi_to_arb(&target.real, source.__re, precision)
    mpfi_to_arb(&target.imag, source.__im, precision)

cdef ComplexIntervalFieldElement acb_to_ComplexIntervalFieldElement(
    const acb_t source,
    const unsigned long precision):
    """
    Convert an ``acb`` to a :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``source`` -- an ``acb_t``

     - ``precision`` -- a positive integer

    OUTPUT:

    A :class:`ComplexIntervalFieldElement`.
    """
    cdef ComplexIntervalFieldElement result

    result = ComplexIntervalField(precision)(0)
    arb_to_mpfi(result.__re, &source.real, precision)
    arb_to_mpfi(result.__im, &source.imag, precision)
    return result

cdef class Acb(SageObject):
    """
    Hold one ``acb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    INPUT:

    - ``value`` -- (default: ``None``) ``None`` or a
      :class:`ComplexIntervalFieldElement`.

    - ``precision`` -- (default: ``0``) a non-negative
      integer. Must be given unless ``value`` is not ``None``.

    OUTPUT:

    None.

    EXAMPLES::

        sage: from sage.rings.complex_interval_acb import Acb # optional - arb
        sage: a = Acb(CIF(1, 1))              # optional - arb
        sage: a.ComplexIntervalFieldElement() # optional - arb
        1 + 1*I
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.complex_interval_acb import Acb # optional - arb
            sage: Acb(precision=2) # optional - arb; indirect doctest
            <type 'sage.rings.complex_interval_acb.Acb'>
        """
        acb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.complex_interval_acb import Acb # optional - arb
            sage: a = Acb(precision=2) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        acb_clear(self.value)

    def __init__(self, value=None, unsigned long precision=0):
        """
        Initialize Acb using value.

        INPUT:

        - ``value`` -- (default: ``None``) ``None`` or a
          :class:`ComplexIntervalFieldElement`.

        - ``precision`` -- (default: ``0``) a non-negative
          integer. Must be given unless ``value`` is not ``None``.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.complex_interval_acb import Acb # optional - arb
            sage: a = Acb(CIF(0, 1))                  # optional - arb
            sage: b = Acb(1)                          # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: value must be None or a ComplexIntervalFieldElement.
            sage: c = Acb() # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: precision must be given.
        """
        cdef ComplexIntervalFieldElement element

        if value is None:
            if precision > 0:
                self._precision_ = precision
            else:
                raise TypeError("precision must be given.")
        elif isinstance(value, ComplexIntervalFieldElement):
            element = <ComplexIntervalFieldElement> value
            self._precision_ = value.parent().precision()
            ComplexIntervalFieldElement_to_acb(self.value, element)

        else:
            raise TypeError("value must be None or a "
                            "ComplexIntervalFieldElement.")

    cpdef ComplexIntervalFieldElement ComplexIntervalFieldElement(self):
        """
        Return :class:`ComplexIntervalFieldElement` of the same value.

        INPUT:

        None.

        OUTPUT:

        A :class:`ComplexIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.complex_interval_acb import Acb # optional - arb
            sage: a = Acb(CIF(2, 2))                # optional - arb
            sage: a.ComplexIntervalFieldElement()   # optional - arb
            2 + 2*I
        """

        return acb_to_ComplexIntervalFieldElement(self.value, self._precision_)
