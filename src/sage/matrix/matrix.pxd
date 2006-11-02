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

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t

cimport sage.structure.element
import  sage.structure.element
cimport sage.structure.mutability

cdef class Matrix(sage.structure.element.ModuleElement):
    # Properties of any matrix  (plus _parent, inherited from base class)
    cdef size_t _nrows
    cdef size_t _ncols
    cdef object _cache
    cdef sage.structure.mutability.Mutability _mutability

    cdef int _will_use_strassen(self, Matrix right) except -1
    cdef int _strassen_default_cutoff(self, Matrix right) except -1

    cdef _mul_cousin_cdef(self, Matrix right)
    cdef _cmp_sibling_cdef(self, Matrix right)

    # Pivots.
    cdef _set_pivots(self, X)

    # Cache
    cdef clear_cache(self)
    cdef fetch(self, key)
    cdef cache(self, key, x)

    # Mutability and bounds checking
    cdef check_bounds(self, size_t i, size_t j)
    cdef check_mutability(self)
    cdef check_bounds_and_mutability(self, size_t i, size_t j)

    # Unsafe entry access
    cdef set_unsafe(self, size_t i, size_t j, object x)
    cdef get_unsafe(self, size_t i, size_t j)

    # Pickling:
    #cdef pickle(self)
    #cdef unpickle(self, data, int version)


    # Strassen
    cdef subtract_strassen_product(result, A, B, int cutoff)

cdef class MatrixWindow:
    cdef size_t _row, _col, _nrows, _ncols
    cdef Matrix _matrix
