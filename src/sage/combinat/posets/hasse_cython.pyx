# -*- coding: utf-8 -*-
r"""
Some fast computations for finite posets
"""
# ****************************************************************************
#       Copyright (C) 2020 Frédéric Chapoton <chapoton@unistra.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_str, sig_off

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mat cimport *

from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.constructor import Matrix


cpdef Matrix_integer_dense moebius_matrix_fast(list positions):
    """
    Compute the Möbius matrix of a poset by a specific triangular inversion.

    INPUT:

    a list of sets describing the poset, as given by the
    lazy attribute ``_leq_storage`` of Hasse diagrams.

    OUTPUT:

    a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import moebius_matrix_fast
        sage: D = [{0,1},{1}]
        sage: moebius_matrix_fast(D)
        [ 1 -1]
        [ 0  1]
        sage: P = posets.TamariLattice(5)
        sage: H = P._hasse_diagram
        sage: D = H._leq_storage
        sage: moebius_matrix_fast(D)
        42 x 42 dense matrix over Integer Ring (...)
    """
    cdef fmpz_mat_t M
    cdef Matrix_integer_dense A
    cdef Py_ssize_t n = len(positions)
    cdef Py_ssize_t i
    cdef int j, k
    sig_str("FLINT exception")
    fmpz_mat_init(M, n, n)
    sig_off()
    fmpz_mat_one(M)
    for i in range(n - 1, -1, -1):
        for j in positions[i]:
            if j != i:
                for k in positions[j]:
                    fmpz_sub(fmpz_mat_entry(M, i, k),
                             fmpz_mat_entry(M, i, k),
                             fmpz_mat_entry(M, j, k))
    A = Matrix(ZZ, n, n)
    fmpz_mat_set(A._matrix, M)
    fmpz_mat_clear(M)
    return A


cpdef Matrix_integer_dense coxeter_matrix_fast(list positions):
    """
    Compute the Coxeter matrix of a poset by a specific algorithm.

    INPUT:

    a list of sets describing the poset, as given by the
    lazy attribute ``_leq_storage`` of Hasse diagrams.

    OUTPUT:

    a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import coxeter_matrix_fast
        sage: D = [{0,1},{1}]
        sage: coxeter_matrix_fast(D)
        [ 0 -1]
        [ 1 -1]
        sage: P = posets.TamariLattice(5)
        sage: H = P._hasse_diagram
        sage: D = H._leq_storage
        sage: coxeter_matrix_fast(D)
        42 x 42 dense matrix over Integer Ring (...)
    """
    cdef fmpz_mat_t M
    cdef Matrix_integer_dense A
    cdef Py_ssize_t n = len(positions)
    cdef Py_ssize_t i
    cdef int j, k
    sig_str("FLINT exception")
    fmpz_mat_init(M, n, n)
    sig_off()
    fmpz_mat_one(M)
    for i in range(n - 1, -1, -1):
        for j in positions[i]:
            if j != i:
                for k in positions[j]:
                    fmpz_sub(fmpz_mat_entry(M, k, i),
                             fmpz_mat_entry(M, k, i),
                             fmpz_mat_entry(M, k, j))
    for i in range(n):
        for j in positions[i]:
            if j != i:
                for k in range(n):
                    fmpz_add(fmpz_mat_entry(M, i, k),
                             fmpz_mat_entry(M, i, k),
                             fmpz_mat_entry(M, j, k))
    fmpz_mat_scalar_mul_si(M, M, -1)
    A = Matrix(ZZ, n, n)
    fmpz_mat_set(A._matrix, M)
    fmpz_mat_clear(M)
    return A
