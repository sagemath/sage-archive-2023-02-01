r"""
A binary matrix datatype in Cython

It's almost a copy of the bitset datatype, but allows a differentiation of the
rows of the matrix. That's the only advantage compared to storing all the rows
of a matrix in a loooooonng bitset.

a ``binary_matrix_t`` structure contains :

- ``long n_cols`` -- number of columns

- ``long n_rows`` -- number of rows

- ``long width`` -- number of ``unsigned long`` per row

- ``bitset_s * rows`` -- ``rows+i`` points toward a block of type ``bitset_t``
  containing the bits of row `i`.

"""
from sage.misc.binary_matrix cimport *
include 'sage/data_structures/bitset.pxi'

cdef inline binary_matrix_init(binary_matrix_t m, long n_rows, long n_cols):
    r"""
    Allocates the binary matrix.
    """
    cdef int i

    m.n_cols = n_cols
    m.n_rows = n_rows
    m.width = (n_cols - 1)/(8*sizeof(unsigned long)) + 1
    m.rows = <bitset_s *>sage_malloc(n_rows * sizeof(bitset_s))
    if m.rows == NULL:
        raise MemoryError

    for i from 0 <= i < n_rows:
        bitset_init(m.rows+i, n_cols)

cdef inline binary_matrix_free(binary_matrix_t m):
    r"""
    Frees the memory allocated by the matrix
    """
    cdef int i

    for i from 0 <= i < m.n_rows:
        bitset_free(m.rows+i)
    sage_free(m.rows)

cdef inline binary_matrix_fill(binary_matrix_t m, bint bit):
    r"""
    Fill the whole matrix with a bit
    """
    cdef int i

    if bit: # set the matrix to 1
        for i from 0 <= i < m.n_rows:
            bitset_set_first_n(m.rows+i, m.n_cols)
    else:
        for i from 0 <= i < m.n_rows:
            bitset_clear(m.rows+i)

cdef inline binary_matrix_complement(binary_matrix_t m):
    r"""
    Complements all of the matrix' bits.
    """
    cdef int i
    for i from 0 <= i < m.n_rows:
        bitset_complement(m.rows+i, m.rows+i)

cdef inline binary_matrix_set1(binary_matrix_t m, long row, long col):
    r"""
    Sets an entry to 1
    """
    bitset_add(m.rows+row, col)

cdef inline binary_matrix_set0(binary_matrix_t m, long row, long col):
    r"""
    Sets an entry to 0
    """
    bitset_discard(m.rows+row, col)

cdef inline binary_matrix_set(binary_matrix_t m, long row, long col, bint value):
    r"""
    Sets an entry
    """
    if value:
        binary_matrix_set1(m,row,col)
    else:
        binary_matrix_set0(m,row,col)

cdef inline bint binary_matrix_get(binary_matrix_t m, long row, long col):
    r"""
    Returns the value of a given entry
    """
    return bitset_in(m.rows+row, col)

cdef inline binary_matrix_print(binary_matrix_t m):
    r"""
    Prints the binary matrix
    """
    cdef int i,j
    import sys
    for i from 0 <= i < m.n_rows:
        # If you want to print the *whole* matrix, including the useless bits,
        # use the following line instead
        #
        # for j in (m.width*8*sizeof(unsigned long)):
        for j from 0 <= j < m.n_cols:
            sys.stdout.write("1" if binary_matrix_get(m, i, j) else ".",)
        print ""
