#include "sage/misc/binary_matrix_pxd.pxi"
cimport sage.misc.binary_matrix
from sage.misc.binary_matrix cimport binary_matrix_t


cdef dict dense_graph_init(binary_matrix_t m, g, translation = ?)
