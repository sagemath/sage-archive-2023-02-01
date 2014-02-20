include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

cdef extern from *:
    int __builtin_popcountl(unsigned long)

# Constants from bitset.pxd
cdef extern from *:
    int index_shift "(sizeof(unsigned long)==8 ? 6 : 5)"
    unsigned long offset_mask "(sizeof(unsigned long)==8 ? 0x3F : 0x1F)"


cdef struct binary_matrix_s:
    long n_cols
    long n_rows

    # Number of "unsigned long" per row
    long width

    unsigned long ** rows

ctypedef binary_matrix_s[1] binary_matrix_t

