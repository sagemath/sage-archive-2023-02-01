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

cimport matrix1

cdef class Matrix(matrix1.Matrix):

    cdef matrix_window_c(self, Py_ssize_t row, Py_ssize_t col, Py_ssize_t nrows, Py_ssize_t ncols)
