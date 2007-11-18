from sage.matrix.all import MatrixSpace
from sage.rings.all import QQ

from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.rings.rational cimport Rational


cdef class Matrix:
    def __init__(self):
        self.M = NULL

    def __repr__(self):
        return "%s x %s Cremona matrix over Rational Field"%(self.nrows(), self.ncols())

    def str(self):
        return self.sage_matrix_over_QQ(sparse=False).str()

    cdef set(self, mat*  M):
        if self.M:
            raise RuntimeError, "self.M is already set."
        self.M = M

    def __dealloc__(self):
        if self.M:
            delete_mat(self.M)

    def nrows(self):
        return nrows(self.M[0])

    def ncols(self):
        return ncols(self.M[0])

    def rank(self):
        return rank(self.M[0])

    def add_scalar(self, scalar s):
        return new_Matrix(addscalar(self.M[0], s))

    def charpoly(self, var='x'):
        return self.sage_matrix_over_QQ(sparse=False).charpoly(var)

    def sage_matrix_over_QQ(self, sparse=True):
        """
        Return corresponding Sage matrix over the rational numbers.

        INPUTS:
            sparse -- (default: True) whether the return matrix has a sparse representation
        """
        cdef long n = self.nrows()
        cdef long i, j, k
        cdef scalar* v = <scalar*> self.M.get_entries()   # coercion needed to deal with const

        cdef Matrix_rational_sparse T
        cdef Matrix_rational_dense Td

        # Ugly code...
        if sparse:
            T = MatrixSpace(QQ, n, sparse=sparse).zero_matrix()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        T.set_unsafe(i, j, Rational(v[k]))
                    k += 1
            return T
        else:
            Td = MatrixSpace(QQ, n, sparse=sparse).zero_matrix()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        Td.set_unsafe(i, j, Rational(v[k]))
                    k += 1
            return Td

cdef class MatrixFactory:
    cdef new_matrix(self, mat M):
        return new_Matrix(M)


cdef Matrix new_Matrix(mat M):
    cdef Matrix A = Matrix()
    A.set(new_mat(M))
    return A
