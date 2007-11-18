from sage.matrix.all import MatrixSpace
from sage.rings.all import QQ

from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.rings.rational cimport Rational


cdef class Matrix:
    """
    A Cremona Matrix.

    EXAMPLES:
        sage: M = CremonaModularSymbols(225)
        sage: t = M.hecke_matrix(2)
        sage: type(t)
        <type 'sage.libs.cremona.mat.Matrix'>
        sage: t
        61 x 61 Cremona matrix over Rational Field
    """
    def __init__(self):
        """
        Called when the matrix is being created.

        EXAMPLES:
            sage: t = CremonaModularSymbols(11).hecke_matrix(2); t
            3 x 3 Cremona matrix over Rational Field
            sage: type(t)
            <type 'sage.libs.cremona.mat.Matrix'>
        """
        self.M = NULL

    def __repr__(self):
        """
        String representation of this matrix.  Use print self.str() to
        print out the matrix entries on the screen.

        EXAMPLES:
            sage: M = CremonaModularSymbols(23)
            sage: t = M.hecke_matrix(2); t
            5 x 5 Cremona matrix over Rational Field
            sage: print t.str()
            [ 3  0  0  0  0]
            [ 0  0  0  1  1]
            [ 0 -1 -1 -1  0]
            [ 1  0 -1 -1 -1]
            [ 0  1  1  0  0]
        """
        return "%s x %s Cremona matrix over Rational Field"%(self.nrows(), self.ncols())

    def str(self):
        """
        Return full string representation of this matrix, never in compact form.

        EXAMPLES:
            sage: M = CremonaModularSymbols(22, sign=1)
            sage: t = M.hecke_matrix(13)
            sage: t.str()
            '[14  0  0  0  0]\n[-4 12  0  8  4]\n[ 0 -6  4 -6  0]\n[ 4  2  0  6 -4]\n[ 0  0  0  0 14]'
        """
        return self.sage_matrix_over_QQ(sparse=False).str()

    cdef set(self, mat*  M):
        if self.M:
            raise RuntimeError, "self.M is already set."
        self.M = M

    def __dealloc__(self):
        if self.M:
            delete_mat(self.M)

    def nrows(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: M = CremonaModularSymbols(19, sign=1)
            sage: t = M.hecke_matrix(13); t
            2 x 2 Cremona matrix over Rational Field
            sage: t.nrows()
            2
        """
        return nrows(self.M[0])

    def ncols(self):
        """
        Return the number of columns of this matrix.

        EXAMPLES:
            sage: M = CremonaModularSymbols(1234, sign=1)
            sage: t = M.hecke_matrix(3); t.ncols()
            156
            sage: M.dimension()
            156
        """
        return ncols(self.M[0])

    # Commented out since it gives very weird
    # results when sign != 0.
##     def rank(self):
##         """
##         Return the rank of this matrix.

##         EXAMPLES:
##             sage: M = CremonaModularSymbols(389)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()
##             65
##             sage: M = CremonaModularSymbols(389, cuspidal=True)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()
##             64

##             sage: M = CremonaModularSymbols(389,sign=1)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()   # known bug.
##             16
##         """
##         return rank(self.M[0])

    def add_scalar(self, scalar s):
        """
        Return new matrix obtained by adding s to each diagonal entry of self.

        EXAMPLES:
            sage: M = CremonaModularSymbols(23, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2); print t.str()
            [ 0  1]
            [ 1 -1]
            sage: w = t.add_scalar(3); print w.str()
            [3 1]
            [1 2]
        """
        return new_Matrix(addscalar(self.M[0], s))

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this matrix, viewed as
        as a matrix over the rational numbers.

        ALGORITHM: Note that currently, this function converts this
        matrix into a dense matrix over the rational numbers, then
        calls the charpoly algorithm on that, which I think is
        Linbox's.

        EXAMPLES:
            sage: M = CremonaModularSymbols(33, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2)
            sage: t.charpoly()
            x^3 + 3*x^2 - 4
            sage: t.charpoly().factor()
            (x - 1) * (x + 2)^2
        """
        return self.sage_matrix_over_QQ(sparse=False).charpoly(var)

    def sage_matrix_over_QQ(self, sparse=True):
        """
        Return corresponding Sage matrix over the rational numbers.

        INPUTS:
            sparse -- (default: True) whether the return matrix has a
                      sparse representation

        EXAMPLES:
            sage: M = CremonaModularSymbols(23, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2)
            sage: s = t.sage_matrix_over_QQ(); s
            [ 0  1]
            [ 1 -1]
            sage: type(s)
            <type 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
            sage: s = t.sage_matrix_over_QQ(sparse=False); s
            [ 0  1]
            [ 1 -1]
            sage: type(s)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
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
