from sage.ext.sparse_matrix_pyx import *


def make_sparse_rational_matrix(nrows, ncols, entries):
    return Matrix_mpq(nrows, ncols, entries = entries, init=False)

