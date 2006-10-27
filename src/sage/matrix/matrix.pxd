"""
Declaration file for generic matrices
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
    cdef int _strassen_default_cutoff(self) except -1

    cdef _mul_cousins_cdef(self, Matrix right)
    cdef classical_multiply_cdef(self, Matrix right)

    # Cache
    cdef clear_cache_cdef(self)
    cdef fetch(self, key)
    cdef cache(self, key, x)

    # Mutability and bounds checking
    cdef check_bounds(self, size_t i, size_t j)
    cdef check_mutability(self)
    cdef check_bounds_and_mutability(self, size_t i, size_t j)

