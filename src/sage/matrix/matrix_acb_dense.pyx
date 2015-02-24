r"""
Arbitrary Precision Complex Interval Matrices using Arb

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

from sage.libs.arb.acb_mat cimport *
from sage.matrix.constructor import matrix
from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse
from sage.rings.complex_interval_field import ComplexIntervalField_class, ComplexIntervalField
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_ball_acb cimport (
    ComplexIntervalFieldElement_to_acb,
    acb_to_ComplexIntervalFieldElement)


cdef void matrix_to_acb_mat(acb_mat_t target, source):
    """
    Convert a matrix containing :class:`ComplexIntervalFieldElement` to an ``acb_mat_t``.

    INPUT:

    - ``target`` -- an ``acb_mat_t``

    - ``source`` -- a matrix consisting of :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef unsigned long nrows, ncols, r, c, precision

    nrows = acb_mat_nrows(target)
    ncols = acb_mat_ncols(target)

    for r in range(nrows):
        for c in range(ncols):
            ComplexIntervalFieldElement_to_acb(acb_mat_entry(target, r, c),
                                               source[r][c])

cdef ComplexIntervalFieldElement _to_CIF(acb_t source, ComplexIntervalFieldElement template):
    cdef ComplexIntervalFieldElement result
    result = template._new()
    acb_to_ComplexIntervalFieldElement(
        result, source)
    return result

cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF):
    """
    Convert an ``acb_mat_t`` to a matrix containing :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``source`` -- an ``acb_mat_t``

    - ``precision`` -- a positive integer.

    OUTPUT:

    A :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`
    containing :class:`ComplexIntervalFieldElement`.
    """
    cdef unsigned long nrows, ncols, r, c
    cdef ComplexIntervalFieldElement template

    nrows = acb_mat_nrows(source)
    ncols = acb_mat_ncols(source)
    template = CIF(0)

    return matrix(
                  [[_to_CIF(acb_mat_entry(source, r, c), template)
                    for c in range(ncols)]
                   for r in range(nrows)])


cdef class Acb_mat(SageObject):
    """
    Class holding an ``acb_mat``, a matrix of complex balls implemented in Arb.

    INPUT:

    - ``nrows`` -- number of rows, a non-negative integer.

    - ``ncols`` -- number of columns, a non-negative integer.

    - ``value`` -- (default: ``None``) ``None`` or a
      :class:`ComplexIntervalFieldElement`.

    - ``precision`` -- (default: ``0``) a non-negative
      integer. Must be given unless ``value`` is not ``None``.

    OUTPUT:

    None.

    EXAMPLES::

        sage: from sage.matrix.matrix_acb_dense import Acb_mat # optional - arb
        sage: Acb_mat(nrows=3, ncols=2, precision=2) # optional - arb; indirect doctest
        <type 'sage.matrix.matrix_acb_dense.Acb_mat'>
    """
    def __cinit__(self,
                  unsigned long nrows,
                  unsigned long ncols,
                  value=None,
                  unsigned long precision=0):
        """
        Allocate memory for the encapsulated value.

        INPUT:

        - ``nrows`` -- number of rows, a non-negative integer.

        - ``ncols`` -- number of columns, a non-negative integer.

        - ``value`` -- (default: ``None``) ``None`` or a matrix.

        - ``precision`` -- (default: ``0``) a non-negative
          integer.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.matrix.matrix_acb_dense import Acb_mat # optional - arb
            sage: Acb_mat(nrows=3, ncols=2, precision=2) # optional - arb; indirect doctest
            <type 'sage.matrix.matrix_acb_dense.Acb_mat'>
        """
        acb_mat_init(self.value, nrows, ncols)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        INPUT:

        None.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.matrix.matrix_acb_dense import Acb_mat # optional - arb
            sage: a = Acb_mat(nrows=3, ncols=2, precision=2) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        acb_mat_clear(self.value)


    def __init__(self,
                 unsigned long nrows,
                 unsigned long ncols,
                 value=None,
                 unsigned long precision=0):
        """
        Initialize Acb using value.

        INPUT:

        - ``nrows`` -- number of rows, a non-negative integer.

        - ``ncols`` -- number of columns, a non-negative integer.

        - ``value`` -- (default: ``None``) ``None`` or a matrix.

        - ``precision`` -- (default: ``0``) a non-negative
          integer. Must be given unless ``value`` is not ``None``.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.matrix.matrix_acb_dense import Acb_mat # optional - arb
            sage: a = Acb_mat(2, 2, matrix([[CIF(1), CIF(1)], [CIF(0), CIF(1)]]))                 # optional - arb
            sage: c = Acb_mat(3, 2) # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: precision must be given.
        """
        if value is None:
            if precision > 0:
                self._precision_ = precision
            else:
                raise TypeError("precision must be given.")
        elif isinstance(value, Matrix_generic_dense) or \
                isinstance(value, Matrix_generic_sparse):
            self._precision_ = value.parent().base_ring().precision()
            matrix_to_acb_mat(self.value, value)

        else:
            raise TypeError("value must be None or a "
                            "matrix.")

    cpdef  Matrix_generic_dense _matrix_(self):
        """
        Return :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense` of the same value.

        INPUT:

        None.

        OUTPUT:

        A :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`.

        EXAMPLES::

            sage: from sage.matrix.matrix_acb_dense import Acb_mat # optional - arb
            sage: a = Acb_mat(2, 2, matrix([[CIF(1), CIF(1)], [CIF(0), CIF(1)]]))                 # optional - arb
            sage: matrix(a)   # optional - arb
            [1 1]
            [0 1]
        """

        return acb_mat_to_matrix(self.value, ComplexIntervalField(self._precision_))
