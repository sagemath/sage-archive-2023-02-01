r"""
A binary matrix datatype in Cython

It's almost a copy of the bitset datatype, but allows a differentiation of the
rows of the matrix. That's the only advantage compared to storing all the rows
of a matrix in a loooooonng bitset.

a ``binary_matrix_t`` structure contains :

- ``long n_cols`` -- number of columns

- ``long n_rows`` -- number of rows

- ``long width`` -- number of ``unsigned long`` per row

- ``unsigned long ** rows`` -- ``rows[i]`` points toward the ``width`` blocks of
  type ``unsigned long`` containing the bits of row `i`.

.. NOTE::

    The rows are stored contiguously in memory, i.e. ``row[i] = row[i-1] +
    width``.
"""
include "binary_matrix_pxd.pxi"

cdef inline binary_matrix_init(binary_matrix_t m, long n_rows, long n_cols):
    r"""
    Allocates the binary matrix.
    """
    cdef int i

    m.n_cols = n_cols
    m.n_rows = n_rows
    m.width = (n_cols - 1)/(8*sizeof(unsigned long)) + 1
    m.rows = <unsigned long**>sage_malloc(n_rows * sizeof(unsigned long *))
    if m.rows == NULL:
        raise MemoryError

    m.rows[0] = <unsigned long*>sage_malloc(n_rows * m.width * sizeof(unsigned long))
    if m.rows[0] == NULL:
        sage_free(m.rows)
        raise MemoryError

    for i in range(1,n_rows):
        m.rows[i] = m.rows[i-1] + m.width
        m.rows[i][m.width-1] = 0

cdef inline binary_matrix_free(binary_matrix_t m):
    r"""
    Frees the memory allocated by the matrix
    """
    sage_free(m.rows[0])
    sage_free(m.rows)

cdef inline binary_matrix_fill(binary_matrix_t m, bint bit):
    r"""
    Fill the whole matrix with a bit
    """
    memset(m.rows[0],-(<char> bit),m.width * m.n_rows * sizeof(unsigned long))

cdef inline binary_matrix_set1(binary_matrix_t m, long row, long col):
    r"""
    Sets an entry to 1
    """
    m.rows[row][col >> index_shift] |= (<unsigned long>1) << (col & offset_mask)

cdef inline binary_matrix_set0(binary_matrix_t m, long row, long col):
    r"""
    Sets an entry to 0
    """
    m.rows[row][col >> index_shift] &= ~((<unsigned long>1) << (col & offset_mask))

cdef inline binary_matrix_set(binary_matrix_t m, long row, long col, bint value):
    r"""
    Sets an entry
    """
    binary_matrix_set0(m,row,col)
    m.rows[row][col >> index_shift] |= (<unsigned long>1) << (value & offset_mask)

cdef inline bint binary_matrix_get(binary_matrix_t m, long row, long col):
    r"""
    Returns the value of a given entry
    """
    return (m.rows[row][col >> index_shift] >> (col & offset_mask)) & 1

cdef inline binary_matrix_print(binary_matrix_t m):
    r"""
    Prints the binary matrix
    """
    cdef int i,j
    import sys
    for i in range(m.n_rows):
        for j in range(m.n_cols):
            sys.stdout.write("1" if binary_matrix_get(m, i, j) else ".",)
        print ""
