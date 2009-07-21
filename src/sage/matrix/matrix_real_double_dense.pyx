"""
Dense matrices over the Real Double Field using NumPy

EXAMPLES::

    sage: b=Mat(RDF,2,3).basis()
    sage: b[0]
    [1.0 0.0 0.0]
    [0.0 0.0 0.0]

We deal with the case of zero rows or zero columns::

    sage: m = MatrixSpace(RDF,0,3)
    sage: m.zero_matrix()
    []

TESTS::

    sage: a = matrix(RDF,2,range(4), sparse=False)
    sage: loads(dumps(a)) == a
    True
    sage: MatrixSpace(RDF,0,0).zero_matrix().inverse()
    []

AUTHORS:

- Jason Grout (2008-09): switch to NumPy backend, factored out the
  Matrix_double_dense class

- Josh Kantor

- William Stein: many bug fixes and touch ups.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################
from sage.rings.real_double import RDF

cimport numpy as cnumpy

numpy=None
scipy=None

cdef class Matrix_real_double_dense(matrix_double_dense.Matrix_double_dense):
    """
    Class that implements matrices over the real double field. These
    are supposed to be fast matrix operations using C doubles. Most
    operations are implemented using numpy which will call the
    underlying BLAS on the system.

    EXAMPLES::

        sage: m = Matrix(RDF, [[1,2],[3,4]])
        sage: m**2
        [ 7.0 10.0]
        [15.0 22.0]
        sage: n= m^(-1); n
        [-2.0  1.0]
        [ 1.5 -0.5]

    To compute eigenvalues the use the functions left_eigenvectors or
    right_eigenvectors

    ::

        sage: p,e = m.right_eigenvectors()

    the result of eigen is a pair (p,e), where p is a list of
    eigenvalues and the e is a matrix whose columns are the
    eigenvectors.

    To solve a linear system Ax = b where A = [[1,2],[3,4]] and
    b = [5,6].

    ::

        sage: b = vector(RDF,[5,6])
        sage: m.solve_left(b)
        (-4.0, 4.5)

    See the commands qr, lu, and svd for QR, LU, and singular value
    decomposition.
    """


    ########################################################################
    # LEVEL 1 functionality
    #   * __cinit__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __richcmp__    -- always the same
    #   * __hash__       -- always simple
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce):
        global numpy
        if numpy is None:
            import numpy
        self._numpy_dtype = numpy.dtype('float64')
        self._numpy_dtypeint = cnumpy.NPY_DOUBLE
        self._python_dtype = float
        # TODO: Make RealDoubleElement instead of RDF for speed
        self._sage_dtype = RDF
        self.__create_matrix__()
        return

    cdef set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value):
        """
        Set the (i,j) entry to value without any type checking or
        bound checking.

        This currently isn't faster than calling self.set_unsafe; should
        we speed it up or is it just a convenience function that has the
        right headers?
        """
        self.set_unsafe(i,j,value)

    cdef double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the (i,j) entry without any type checking or bound checking.

        This currently isn't faster than calling self.get_unsafe; should
        we speed it up or is it just a convenience function that has the
        right headers?
        """
        return self.get_unsafe(i,j)

    def _cholesky_decomposition_(self):
        r"""
        Return the Cholesky factorization of this matrix.

        The input matrix must be symmetric and positive definite or
        ``ValueError`` exception will be raised.

        EXAMPLES:
            sage: M = MatrixSpace(RDF,5)
            sage: r = matrix(RDF,[[   0.,    0.,    0.,    0.,    1.],[   1.,    1.,    1.,    1.,    1.],[  16.,    8.,    4.,    2.,    1.],[  81.,   27.,    9.,    3.,    1.],[ 256.,   64.,   16.,    4.,    1.]])

            sage: m = r*M.identity_matrix()*r.transpose()
            sage: L = m.cholesky_decomposition() # indirect doctest
            sage: L*L.transpose()
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        if self._nrows == 0:   # special case
            return self.copy()

        # cdef matrix_double_dense

        cdef Matrix_real_double_dense M = self._new()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        from numpy.linalg import LinAlgError
        try:
            M._matrix_numpy = scipy.linalg.cholesky(self._matrix_numpy, lower=1)
        except LinAlgError:
            raise ValueError, "The input matrix was not symmetric and positive definite"

        return M
