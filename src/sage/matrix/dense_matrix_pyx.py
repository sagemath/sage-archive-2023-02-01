from sage.ext.dense_matrix_pyx import *

def make_rational_matrix(nrows, ncols, entries):
    return Matrix_rational(nrows, ncols, entries = entries, construct = True)
