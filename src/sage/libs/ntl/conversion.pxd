r"""
Inline conversions between LinBox and Sage

Each LinBox type has a corresponding Sage types and we use the following
conventions for conversion functions

- ``new_ntl_XXX`` : create a new ntl object
- ``new_sage_XXX``   : create a new Sage object
- ``set_ntl_XXX`` : set the entries of the ntl object
- ``set_sage_XXX``   : set the entries of the Sage object

For matrices that uses a flint datastructure, see the lower level conversions
in the module ``ntl_flint_interface``.
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

from .ntl cimport mat_ZZ_p_c, ZZ_pContext

from sage.matrix.matrix_modn_dense cimport Matrix_modn_dense


################################################
# matrix_modn_sparse (sparse matrix over Z/nZ) #
################################################

cdef inline void set_ntl_matrix_modn_sparse(mat_ZZ_p_c& A, Matrix_modn_dense m):
    r"""
    set the entries of a LinBox matrix from a Sage matrix.

    INPUT:

    - A -- LinBox matrix
    - m -- Sage matrix
    """
    cdef size_t i, j
    for i in range(m._nrows):
        for j in range(m._ncols):
            A.put(i, j, m[i,j])

cdef inline mat_ZZ_p_c * new_ntl_matrix_modn_sparse(ZZ_pContext &F, Matrix_modn_dense m):
    r"""
    Return a new LinBox matrix from a Sage matrix.

    Such matrix has to be deallocated with a "del" statement.

    INPUT:

    - F -- LinBox field
    - m -- Sage matrix
    """
    cdef mat_ZZ_p_c * A = new mat_ZZ_p_c(m._nrows, m._ncols)
    set_ntl_matrix_modn_sparse(A[0], m)
    return A
