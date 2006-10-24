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

cdef class Matrix(sage.structure.element.ModuleElement):
    cdef object _mutability

    # todo -- delete these two
    cdef object __nrows
    cdef object __ncols

    cdef object __dict
    cdef object __determinant
    cdef object __charpoly
    cdef object __sparse_columns
    cdef object __sparse_rows
    cdef object __eigenvectors
    cdef object __rank
    cdef object __echelon_form

    cdef size_t _nrows
    cdef size_t _ncols
