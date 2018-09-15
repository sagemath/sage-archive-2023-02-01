r"""
Inline conversions between LinBox and Sage

Each LinBox type has a corresponding Sage types and we use the following
conventions for conversion functions

- ``new_linbox_XXX`` : create a new linbox object
- ``new_sage_XXX``   : create a new Sage object
- ``set_linbox_XXX`` : set the entries of the linbox object
- ``set_sage_XXX``   : set the entries of the Sage object

For matrices that uses a flint datastructure, see the lower level conversions
in the module ``linbox_flint_interface``.
"""
#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#       Copyright (C) 2008 Clement Pernet
#       Copyright (C) 2018 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .givaro cimport Modular_int64
from .linbox cimport SparseMatrix_Modular_int64

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse

from sage.modules.vector_modn_sparse cimport c_vector_modint

################################################
# matrix_modn_sparse (sparse matrix over Z/nZ) #
################################################

cdef inline void set_linbox_matrix_modn_sparse(SparseMatrix_Modular_int64& A, Matrix_modn_sparse m):
    r"""
    set the entries of a LinBox matrix from a Sage matrix.

    INPUT:

    - A -- LinBox matrix
    - m -- Sage matrix
    """
    cdef c_vector_modint * row
    cdef size_t i, j
    for i in range(m._nrows):
        row = m.rows + i
        for j in range(row.num_nonzero):
            A.setEntry(i, row.positions[j], row.entries[j])

cdef inline SparseMatrix_Modular_int64 * new_linbox_matrix_modn_sparse(Modular_int64 &F, Matrix_modn_sparse m):
    r"""
    Return a new LinBox matrix from a Sage matrix.

    Such matrix has to be deallocated with a "del" statement.

    INPUT:

    - F -- LinBox field
    - m -- Sage matrix
    """
    cdef SparseMatrix_Modular_int64 * A = new SparseMatrix_Modular_int64(F, m._nrows, m._ncols)
    set_linbox_matrix_modn_sparse(A[0], m)
    return A


