from sage.data_structures.bitset cimport bitset_t

cdef struct binary_matrix_s:
    long n_cols
    long n_rows
    bitset_t* rows
    
ctypedef binary_matrix_s[1] binary_matrix_t
