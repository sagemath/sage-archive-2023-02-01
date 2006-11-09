"""
Generic matrices
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport sage.structure.element
import  sage.structure.element
cimport sage.structure.mutability

cdef class Matrix(sage.structure.element.ModuleElement):
    # Properties of any matrix  (plus _parent, inherited from base class)
    cdef Py_ssize_t _nrows
    cdef Py_ssize_t _ncols
    cdef object _cache
    cdef public object _base_ring
    cdef sage.structure.mutability.Mutability _mutability

    cdef int _will_use_strassen(self, Matrix right) except -2
    cdef int _will_use_strassen_echelon(self) except -2
    cdef int _strassen_default_cutoff(self, Matrix right) except -2
    cdef int _strassen_default_echelon_cutoff(self) except -2

    cdef _mul_c_impl(self, Matrix right)
    cdef int _cmp_c_impl(self, sage.structure.element.Element right) except -2

    cdef long _hash(self) except -1

    # Pivots.
    cdef _set_pivots(self, X)

    # Cache
    cdef clear_cache(self)
    cdef fetch(self, key)
    cdef cache(self, key, x)

    # Mutability and bounds checking
    cdef check_bounds(self, Py_ssize_t i, Py_ssize_t j)
    cdef check_mutability(self)
    cdef check_bounds_and_mutability(self, Py_ssize_t i, Py_ssize_t j)

    # Unsafe entry access
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x)
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j)
    cdef _coerce_element(self, x)

    # Strassen
    cdef subtract_strassen_product(self, result, A, B, int cutoff)

    # Row and column operations
    cdef check_row_bounds_and_mutability(self, Py_ssize_t r1, Py_ssize_t r2)
    cdef check_column_bounds_and_mutability(self, Py_ssize_t c1, Py_ssize_t c2)
    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2)
    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2)
    cdef add_multiple_of_row_c(self, Py_ssize_t i, Py_ssize_t j,    s,   Py_ssize_t col_start)
    cdef add_multiple_of_column_c(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t row_start)
    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col)
    cdef rescale_col_c(self, Py_ssize_t i, s, Py_ssize_t start_row)





cdef class MatrixWindow:
    cdef Py_ssize_t _row, _col, _nrows, _ncols
    cdef Matrix _matrix
