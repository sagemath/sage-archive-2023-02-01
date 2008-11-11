"""
Dense matrices using a numpy backend.  This serves as a base class for
dense matrices over Real Double Field and Complex Double Field.

EXAMPLES:
    sage: b=Mat(RDF,2,3).basis()
    sage: b[0]
    [1.0 0.0 0.0]
    [0.0 0.0 0.0]


We deal with the case of zero rows or zero columns:
    sage: m = MatrixSpace(RDF,0,3)
    sage: m.zero_matrix()
    []

TESTS:
    sage: a = matrix(RDF,2,range(4), sparse=False)
    sage: loads(dumps(a)) == a
    True
    sage: a = matrix(CDF,2,range(4), sparse=False)
    sage: loads(dumps(a)) == a
    True

AUTHORS:
    -- Jason Grout, Sep 2008: switch to numpy backend, factored out the Matrix_double_dense class
    -- Josh Kantor
    -- William Stein: many bug fixes and touch ups.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################
import math

import sage.rings.real_double
import sage.rings.complex_double
from matrix cimport Matrix
from sage.structure.element cimport ModuleElement,Vector
from constructor import matrix
from sage.modules.free_module_element import vector
cimport sage.structure.element
from matrix_space import MatrixSpace

cimport sage.ext.numpy as cnumpy

numpy=None
scipy=None

# This is for the Numpy C API to work
cnumpy.import_array()


from matrix_complex_double_dense import Matrix_complex_double_dense

cdef class Matrix_double_dense(matrix_dense.Matrix_dense):
    """
    Base class for matrices over the Real Double Field and the Complex
    Double Field.  These are supposed to be fast matrix operations
    using C doubles. Most operations are implemented using numpy which
    will call the underlying BLAS on the system.

    This class cannot be instantiated on its own.  The numpy matrix
    creation depends on several variables that are set in the
    subclasses.

    EXAMPLES:
        sage: m = Matrix(RDF, [[1,2],[3,4]])
        sage: m**2
        [ 7.0 10.0]
        [15.0 22.0]
        sage: n= m^(-1); n
        [-2.0  1.0]
        [ 1.5 -0.5]
    """


    ########################################################################
    # LEVEL 1 functionality
    #   * __new__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __richcmp__    -- always the same
    #   * __hash__       -- alway simple
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        """
        Set up a new matrix
        """
        matrix_dense.Matrix_dense.__init__(self,parent)
        return

    def __create_matrix__(self):
        """
        Create a new uninitialized numpy matrix to hold the data for the class.

        This function assumes that self._numpy_dtypeint and
        self._nrows and self._ncols have already been initialized.

        EXAMPLE:
        In this example, we throw away the current matrix and make a
        new unitialized matrix representing the data for the class.
            sage: a=matrix(RDF, 3, range(9))
            sage: a.__create_matrix__()
        """
        cdef cnumpy.npy_intp dims[2]
        dims[0] = self._nrows
        dims[1] = self._ncols
        self._matrix_numpy = cnumpy.PyArray_SimpleNew(2, dims, self._numpy_dtypeint)
        return

    def __dealloc__(self):
        """ Deallocate any memory that was initialized."""
        return

    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    def __hash__(self):
        """
        Hash this matrix if it's immutable.

        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            88
            sage: A = matrix(CDF,3,range(1,10))
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            88
        """
        return self._hash()

    def LU_valid(self):
        r"""
        Returns \code{True} if the LU form of this matrix has
        already been computed.

        EXAMPLES:
            sage: A= random_matrix(RDF,3) ; A.LU_valid()
            False
            sage: L,U,P = A.LU()
            sage: A.LU_valid()
            True
        """
        return self.fetch('LU_valid') == True

    def __init__(self, parent, entries, copy, coerce):
        """
        Fill the matrix with entries.

        The numpy matrix must have already been allocated.

        EXAMPLES:
            sage: matrix(RDF,3,range(9))
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: matrix(CDF,3,3,2)
            [2.0   0   0]
            [  0 2.0   0]
            [  0   0 2.0]

        TESTS:
            sage: matrix(RDF,3,0)
            []
            sage: matrix(RDF,3,3,0)
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            sage: matrix(RDF,3,3,1)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]
            sage: matrix(RDF,3,3,2)
            [2.0 0.0 0.0]
            [0.0 2.0 0.0]
            [0.0 0.0 2.0]
            sage: matrix(CDF,3,0)
            []
            sage: matrix(CDF,3,3,0)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: matrix(CDF,3,3,1)
            [1.0   0   0]
            [  0 1.0   0]
            [  0   0 1.0]
            sage: matrix(CDF,3,3,range(9))
            [  0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: matrix(CDF,2,2,[CDF(1+I)*j for j in range(4)])
            [          0 1.0 + 1.0*I]
            [2.0 + 2.0*I 3.0 + 3.0*I]
        """
        cdef Py_ssize_t i,j
        cdef cnumpy.npy_intp dims[2]
        dims[0] = self._nrows
        dims[1] = self._ncols
        if isinstance(entries,(tuple, list)):
            if len(entries)!=self._nrows*self._ncols:
                    raise TypeError, "entries has wrong length"

            if coerce:
                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        self.set_unsafe(i,j,self._python_dtype(entries[i*self._ncols+j]))
            else:
                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        self.set_unsafe(i,j,entries[i*self._ncols+j])

        else:
            cnumpy.PyArray_FILLWBYTE(self._matrix_numpy, 0)

            if entries is None:
                z = self._python_dtype(0.0)
            else:
                try:
                    z = self._python_dtype(entries)
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or float"
            if z != 0:
                if self._nrows != self._ncols:
                    raise TypeError, "scalar matrix must be square"
                for i from 0<=i<self._ncols:
                    self.set_unsafe(i,i,z)



    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object value):
        """
        Set the (i,j) entry to value without any bounds checking,
        mutability checking, etc.
        """
        # We assume that Py_ssize_t is the same as cnumpy.npy_intp

        # We must patch the ndarrayobject.h file so that the SETITEM
        # macro does not have a semicolon at the end for this to work.
        # Cython wraps the macro in a function that converts the
        # returned int to a python object, which leads to compilation
        # errors because after preprocessing you get something that
        # looks like "););".  This is bug
        # http://scipy.org/scipy/numpy/ticket/918

        # We call the self._python_dtype function on the value since
        # numpy does not know how to deal with complex numbers other
        # than the builtin complex number type.
        cdef int status
        status = cnumpy.PyArray_SETITEM(self._matrix_numpy,
                        cnumpy.PyArray_GETPTR2(self._matrix_numpy, i, j),
                        self._python_dtype(value))
        #TODO: Throw an error if status == -1


    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the (i,j) entry without any bounds checking, etc.
        """
        # We assume that Py_ssize_t is the same as cnumpy.npy_intp
        return self._sage_dtype(cnumpy.PyArray_GETITEM(self._matrix_numpy,
                                                cnumpy.PyArray_GETPTR2(self._matrix_numpy, i, j)))

    cdef Matrix_double_dense _new(self, int nrows=-1, int ncols=-1):
        """
        Return a new uninitialized matrix with same parent as self.

        INPUT
            nrows -- (default self._nrows) number of rows in returned matrix
            ncols -- (default self._ncols) number of columns in returned matrix

        """
        cdef Matrix_double_dense m
        if nrows == -1:
            nrows = self._nrows
        if ncols == -1:
            ncols = self._ncols
        parent = self.matrix_space(nrows, ncols)
        m = self.__class__.__new__(self.__class__,parent,None,None,None)
        return m




    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two matrices together.

        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
            sage: A+A
            [ 2.0  4.0  6.0]
            [ 8.0 10.0 12.0]
            [14.0 16.0 18.0]
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.copy()

        cdef Matrix_double_dense M, _right, _left
        _right = right
        _left = self

        M = self._new()
        M._matrix_numpy = _left._matrix_numpy + _right._matrix_numpy
        return M

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Return self - right


        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
            sage: (A-A).is_zero()
            True
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.copy()

        cdef Matrix_double_dense M,_right,_left
        _right = right
        _left = self

        M = self._new()
        M._matrix_numpy = _left._matrix_numpy - _right._matrix_numpy
        return M

    def __neg__(self):
        """
        Negate this matrix

        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
            sage: -A
            [-1.0 -2.0 -3.0]
            [-4.0 -5.0 -6.0]
            [-7.0 -8.0 -9.0]
            sage: B = -A ; (A+B).is_zero()
            True
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.copy()

        cdef Matrix_double_dense M
        M = self._new()
        M._matrix_numpy = -self._matrix_numpy
        return M


    #   * cdef _cmp_c_impl
    #   * __copy__
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):                        #unsure how to implement
    # def _unpickle(self, data, int version):   # use version >= 0 #unsure how to implement
    ######################################################################
    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):
        """
        Multiply self*right as matrices.

        EXAMPLES
            sage: A = matrix(RDF,3,range(1,10))
            sage: B = matrix(RDF,3,range(1,13))
            sage: A*B
            [ 38.0  44.0  50.0  56.0]
            [ 83.0  98.0 113.0 128.0]
            [128.0 152.0 176.0 200.0]
        """
        if self._ncols!=right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right"

        if self._nrows == 0 or self._ncols == 0 or right._nrows == 0 or right._ncols == 0:
            return self.matrix_space(self._nrows, right._ncols).zero_matrix()

        cdef Matrix_double_dense M,_right,_left
        M = self._new(self._nrows, right._ncols)
        _right = right
        _left = self
        global numpy
        if numpy is None:
            import numpy

        M._matrix_numpy = numpy.dot(_left._matrix_numpy, _right._matrix_numpy)
        return M

    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    def __invert__(self):
        """
        Invert this matrix.

        EXAMPLES:
            sage: A = Matrix(RDF, [[10, 0], [0, 100]])
            sage: (~A).det()
            0.001

            sage: A = matrix(RDF,3,[2,3,5,7,8,9,11,13,17]);A
            [ 2.0  3.0  5.0]
            [ 7.0  8.0  9.0]
            [11.0 13.0 17.0]
            sage: ~A
            [ -2.71428571429            -2.0   1.85714285714]
            [  2.85714285714             3.0  -2.42857142857]
            [-0.428571428571            -1.0  0.714285714286]

            Note that if this matrix is (nearly) singular, finding
            its inverse will not help much and will give slightly different
            answers on similar platforms depending on the hardware
            and tuning options given to ATLAS:
            sage: A = Matrix(RDF, [[1, 0], [0, 0]])




            sage: A = matrix(RDF,3,range(1,10));A
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            [7.0 8.0 9.0]

            sage: A.determinant() < 10e-12
            True




        """
# see trac ticket 4502 --- there is an issue with the "#random" pragma that needs to be fixed
#                          as for the mathematical side, scipy v0.7 is expected to fix the invertibility failures
#
#            sage: A = Matrix(RDF, [[1, 0], [0, 0]])
#            sage: A.inverse().det()        # random - on some computers, this will be invertible due to numerical error.
#            Traceback (most recent call last):
#            ...
#            LinAlgError: singular matrix
#            sage: A = matrix(RDF,3,range(1,10));A
#            [1.0 2.0 3.0]
#            [4.0 5.0 6.0]
#            [7.0 8.0 9.0]
#
#            sage: A.determinant() < 10e-12
#            True
#            sage: ~A                       # random - on some computers, this will be invertible due to numerical error.
#            Traceback (most recent call last):
#            ...
#            LinAlgError: singular matrix
#
        if self._nrows == 0 or self._ncols == 0:
            return self.copy()

        # Maybe we should cache the (P)LU decomposition and use scipy.lu_solve?
        cdef Matrix_double_dense M
        M = self._new()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        M._matrix_numpy = scipy.linalg.inv(self._matrix_numpy)
        return M

    # def __copy__(self):
    # def _list(self):
    # def _dict(self):


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #
    #    compute_LU(self)
    #    get_LU  #add
    #
    ########################################################################
    cdef _c_compute_LU(self):
        """
        Compute the LU factorization and store it in self._L, self._P,
        and self._U.
        """
        if not hasattr(self, '_L'):
            # It's an awful waste to store a huge matrix for P, which
            # is just a simple permutation, really.
            global scipy
            if scipy is None:
                import scipy
            import scipy.linalg
            self._P, self._L, self._U = scipy.linalg.lu(self._matrix_numpy)
            # Numpy has a different convention than we had with GSL
            # So we invert (transpose) the P to match our prior behavior
            self._P = self._P.T
            self.cache('LU_valid',True)

    def LU(self):
        """
        Computes the LU decomposition of a matrix.

        For and square matrix A we can find matrices P,L, and U. s.t.
           P*A = L*U
        where P is a permutation matrix, L is lower triangular and U
        is upper triangular.

        OUTPUT:
            P, L, U -- as a tuple

        EXAMPLES:
            sage: m = matrix(RDF,4,range(16))
            sage: P,L,U = m.LU()
            sage: P*m
            [12.0 13.0 14.0 15.0]
            [ 0.0  1.0  2.0  3.0]
            [ 8.0  9.0 10.0 11.0]
            [ 4.0  5.0  6.0  7.0]
            sage: L*U
            [12.0 13.0 14.0 15.0]
            [ 0.0  1.0  2.0  3.0]
            [ 8.0  9.0 10.0 11.0]
            [ 4.0  5.0  6.0  7.0]
        """

        if self._ncols!=self._nrows:
            raise TypeError,"LU decomposition only works for square matrix"

        if self._ncols == 0:
            return self.copy(), self.copy(), self.copy()

        if self.fetch('LU_valid')!=True:
            self._c_compute_LU()
        cdef Py_ssize_t i,j,k,l,copy_result
        cdef Matrix_double_dense P, L,U
        P = self._new()
        L = self._new()
        U = self._new()
        P._matrix_numpy = self._P.copy()
        L._matrix_numpy = self._L.copy()
        U._matrix_numpy = self._U.copy()
        return P, L, U

    def eigenspaces(self, var='a'):
        r"""
        Return a list of pairs (e, V) where e runs through all complex
        eigenvalues of this matrix, and V is the corresponding
        left eigenspace (always a 1-dimensional complex vector space).

        EXAMPLES:
            sage: m = matrix(RDF, 3, range(9)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: es = m.eigenspaces()
            sage: es # random
            [(13.3484692283, Vector space of degree 3 and dimension 1 over Real Double Field
            User basis matrix:
            [-0.440242867236 -0.567868371314 -0.695493875393]),
            (-1.34846922835, Vector space of degree 3 and dimension 1 over Real Double Field
            User basis matrix:
            [-0.897878732262 -0.278434036822  0.341010658618]),
            (-9.10854412047e-16, Vector space of degree 3 and dimension 1 over Real Double Field
            User basis matrix:
            [ 0.408248290464 -0.816496580928  0.408248290464])]

            sage: e, v = es[0]
            sage: v = v.basis()[0]
            sage: a = v * m
            sage: b = e * v
            sage: diff = a.change_ring(CDF) - b
            sage: abs(abs(diff)) < 1e-10
            True
            sage: diff # random -- very small numbers
            (-2.6645352591e-15, -7.1054273576e-15, -3.5527136788e-15)
        """
        e, v = self.left_eigenvectors()
        v = v.rows()
        pairs = []
        for l from 0<=l<len(e):
            c = v[l]
            pairs.append((e[l], c.parent().span_of_basis([c], check=False)))
        return pairs

    def left_eigenvectors(self):
        """
        Computes the eigenvalues and *left* eigenvectors of this
        matrix m acting *from the right*.  I.e., vectors v such that
        v*m = lambda*v.

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as an RDF matrix whose rows
                           are the eigenvectors.

        EXAMPLES:
            sage: m = Matrix(RDF, 3, range(9))
            sage: es = m.left_eigenvectors()
            sage: es # random-ish platform-dependent output (low order digits)
            ([13.3484692283, -1.34846922835, -9.10854412047e-16],
            [-0.440242867236 -0.567868371314 -0.695493875393]
            [-0.897878732262 -0.278434036822  0.341010658618]
            [ 0.408248290464 -0.816496580928  0.408248290464])
            sage: e, v = es; e = e[0]; v = v[0]
            sage: abs(abs(e*v - v*m)) < 1e-10
            True
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self.copy()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v,eig = scipy.linalg.eig(self._matrix_numpy, right=False, left=True)
        eig = matrix(eig.T)
        return ([sage.rings.complex_double.CDF(x) for x in v], eig)

    eigenvectors_left = left_eigenvectors

    def right_eigenvectors(self):
        """
        Computes the eigenvalues and *right* eigenvectors of this
        matrix m acting *from the left*.  I.e., vectors v such that
        m * v = lambda*v, where v is viewed as a column vector.

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as an RDF matrix whose columns
                           are the eigenvectors.

        EXAMPLES:
            sage: m = Matrix(RDF, 3, range(9))
            sage: evals,evecs = m.right_eigenvectors()
            sage: evals # random low-order bits
            [13.3484692283, -1.34846922835, -1.1327908706e-15]
            sage: evecs # random low-order bits
            [ 0.164763817282  0.799699663112  0.408248290464]
            [ 0.505774475901  0.104205787719 -0.816496580928]
            [ 0.846785134519 -0.591288087674  0.408248290464]
            sage: max([max(m*evec - evec*eval) for eval,evec in zip(evals, evecs.columns())]) < 1e-14
            True
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self.copy()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v,eig = scipy.linalg.eig(self._matrix_numpy, right=True, left=False)
        return ([sage.rings.complex_double.CDF(x) for x in v], matrix(eig))

    eigenvectors_right = right_eigenvectors

    def solve_left_LU(self, b):
        """
        Solve the equation A*x = b.

        INPUT:
           self -- an invertible matrix
           b -- a vector

        NOTES: This method precomputes and stores the LU decomposition
        before solving. If many equations of the form Ax=b need to be
        solved for a singe matrix A, then this method should be used
        instead of solve. The first time this method is called it will
        compute the LU decomposition.  If the matrix has not changed
        then subsequent calls will be very fast as the precomputed LU
        decomposition will be reused.

        EXAMPLES:
            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF,[1,2,3])
            sage: x = A.solve_left(b); x
            (-0.113695090439, 1.39018087855, -0.333333333333)
            sage: A*x
            (1.0, 2.0, 3.0)

        TESTS:
        We test two degenerate cases:
            sage: A = matrix(RDF, 0, 3, [])
            sage: A.solve_left_LU(vector(RDF,[]))
            (0.0, 0.0, 0.0)
            sage: A = matrix(RDF, 3, 0, [])
            sage: A.solve_left_LU(vector(RDF,3, [1,2,3]))
            ()

        """
        if self._nrows != b.degree():
            raise ValueError, "number of rows of self must equal degree of b"
        if self._nrows == 0 or self._ncols == 0:
            return self._row_ambient_module().zero_vector()

        cdef Matrix_double_dense M = self._new()
        lu = self._L*self._U
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        M._matrix_numpy = scipy.linalg.lu_solve((lu, self._P), b)
        return M

    def solve_left(self,vec):
        """
        Solve the equation A*x = b, where

        EXAMPLES:
            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF,[1,2,3])
            sage: x = A.solve_left(b); x
            (-0.113695090439, 1.39018087855, -0.333333333333)
            sage: A*x
            (1.0, 2.0, 3.0)

        TESTS:
        We test two degenerate cases:
            sage: A = matrix(RDF, 0, 3, [])
            sage: A.solve_left(vector(RDF,[]))
            ()
            sage: A = matrix(RDF, 3, 0, [])
            sage: A.solve_left(vector(RDF,3, [1,2,3]))
            (0.0, 0.0, 0.0)
        """
        M = self._column_ambient_module()
        if self._nrows == 0 or self._ncols == 0:
            return M.zero_vector()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        return M(scipy.linalg.solve(self._matrix_numpy,vec.numpy()))

    def determinant(self):
        """
        Return the determinant of self.

        ALGORITHM: Use numpy

        EXAMPLES:
            sage: m = matrix(RDF,2,range(4)); m.det()
            -2.0
            sage: m = matrix(RDF,0,[]); m.det()
            1.0
            sage: m = matrix(RDF, 2, range(6)); m.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be square
        """
        if self._nrows == 0 or self._ncols == 0:
            return self._sage_dtype(1)
        if not self.is_square():
            raise ValueError, "self must be square"
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        return self._sage_dtype(scipy.linalg.det(self._matrix_numpy))


    def log_determinant(self):
        """
        Compute the log of the absolute value of the determinant
        using LU decomposition.

        NOTE: This is useful if the usual determinant overlows.

        EXAMPLES:
            sage: m = matrix(RDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: RDF(log(abs(m.determinant())))
            0.69314718056
            sage: m.log_determinant()
            0.69314718056
            sage: m = matrix(RDF,0,0,[]); m
            []
            sage: m.log_determinant()
            0.0
            sage: m = matrix(CDF,2,2,range(4)); m
            [  0 1.0]
            [2.0 3.0]
            sage: RDF(log(abs(m.determinant())))
            0.69314718056
            sage: m.log_determinant()
            0.69314718056
            sage: m = matrix(CDF,0,0,[]); m
            []
            sage: m.log_determinant()
            0.0

        """
        if self._nrows == 0 or self._ncols == 0:
            return sage.rings.real_double.RDF(0)

        if not self.is_square():
            raise ValueError, "self must be square"

        if(self.fetch('LU_valid') !=True):
            self._c_compute_LU()
        global numpy
        if numpy is None:
            import numpy

        return sage.rings.real_double.RDF(sum(numpy.log(abs(numpy.diag(self._U)))))

    def transpose(self):
        """
        Return the transpose of this matrix.

        EXAMPLES:
            sage: m = matrix(RDF,2,3,range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m.transpose()
            [0.0 3.0]
            [1.0 4.0]
            [2.0 5.0]
            sage: m = matrix(RDF,0,3); m
            []
            sage: m.transpose()
            []
            sage: m.transpose().parent()
            Full MatrixSpace of 3 by 0 dense matrices over Real Double Field
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._ncols, self._nrows)

        cdef Matrix_double_dense trans
        trans = self._new(self._ncols, self._nrows)
        trans._matrix_numpy = self._matrix_numpy.transpose()
        if self.subdivisions is not None:
            row_divs, col_divs = self.get_subdivisions()
            trans.subdivide(col_divs, row_divs)
        return trans

    def SVD(self, *args, **kwds):
        r"""
        Return the singular value decomposition of this matrix.

        INPUT:
            A -- a matrix
        OUTPUT:
            U, S, V -- matrices such that $A = U*S*V.conj().transpose()$, where
                       U and V are orthogonal and S is zero off of the diagonal.

        Note that if self is m-by-n, then the dimensions of the
        matrices that this returns are (m,m), (m,n), and (n, n).

        EXAMPLES:
            sage: m = matrix(RDF,4,range(16))
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()    # slightly random output (due to computer architecture)
            [3.45569519412e-16               1.0               2.0               3.0]
            [4.0               5.0               6.0               7.0]
            [8.0               9.0              10.0              11.0]
            [12.0              13.0              14.0              15.0]

        A non-square example:
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: U, S, V = m.SVD()
            sage: U*S*V.transpose()           # random low bits
            [7.62194127257e-17               1.0               2.0]
            [              3.0               4.0               5.0]
            sage: U
            [-0.274721127897 -0.961523947641]
            [-0.961523947641  0.274721127897]
            sage: S
            [7.34846922835           0.0           0.0]
            [          0.0           1.0           0.0]
            sage: V
            [-0.392540507864  0.824163383692  0.408248290464]
            [-0.560772154092  0.137360563949 -0.816496580928]
            [ -0.72900380032 -0.549442255795  0.408248290464]
            sage: U*U.transpose()            # random low bits
            [              1.0 2.13506512817e-16]
            [2.13506512817e-16               1.0]
            sage: max((U*U.transpose()-identity_matrix(2)).list())<1e-15
            True
            sage: V*V.transpose()            # random low bits
            [               1.0  2.02230810223e-16 -2.11947972194e-16]
            [ 2.02230810223e-16                1.0  7.09339271349e-17]
            [-2.11947972194e-16  7.09339271349e-17                1.0]
            sage: max((V*V.transpose()-identity_matrix(3)).list())<1e-15
            True
            sage: max((U*S*V.transpose()-m).list())<1e-15 # check
            True
            sage: m = matrix(RDF,3,2,range(6)); m
            [0.0 1.0]
            [2.0 3.0]
            [4.0 5.0]
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()   # random low order bits
            [-8.13151629364e-19                1.0]
            [               2.0                3.0]
            [               4.0                5.0]

	Due to numerical noise issues on Intel Macs, the following fails if 1e-14
	is changed to 1e-15:
            sage: max((U*S*V.transpose()-m).list())<1e-14 # check
            True

        TESTS:
            sage: m = matrix(RDF, 3, 0, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: m = matrix(RDF, 0, 3, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: def shape(x): return (x.nrows(), x.ncols())
            sage: m = matrix(RDF, 2, 3, range(6))
            sage: map(shape, m.SVD())
            [(2, 2), (2, 3), (3, 3)]
        """
        if len(args)>0 or len(kwds)>0:
            from sage.misc.misc import deprecation
            deprecation("Arguments passed to SVD, but SVD no longer supports different methods (it only uses numpy now).")

        if self._nrows == 0 or self._ncols == 0:
            U_t = self.new_matrix(self._nrows, self._ncols)
            S_t = self.new_matrix(self._nrows, self._ncols)
            V_t = self.new_matrix(self._ncols, self._nrows)
            return U_t, S_t, V_t


        cdef int i, s_dim
        cdef Matrix_double_dense S = self.new_matrix()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        U,_S,V = scipy.linalg.svd(self._matrix_numpy)

        s_dim = len(_S)
        for i from 0 <= i < s_dim:
            S[(i,i)] = _S[i]

        return (matrix(U),S,matrix(V.conj().T))

    def QR(self):
        """
        Return the Q,R factorization of a real matrix A.

        INPUT:
           self -- a real matrix A

        OUTPUT:
           Q, R -- matrices such that A = Q*R such that the columns of Q are
                   orthogonal (i.e., $Q^t Q = I$), and R is upper triangular.

        EXAMPLES:
            sage: m = matrix(RDF,3,range(12)); m
            [ 0.0  1.0  2.0  3.0]
            [ 4.0  5.0  6.0  7.0]
            [ 8.0  9.0 10.0 11.0]
            sage: Q,R = m.QR()
            sage: Q*R
            [ 0.0  1.0  2.0  3.0]
            [ 4.0  5.0  6.0  7.0]
            [ 8.0  9.0 10.0 11.0]

        Note that the columns of Q will be an orthogonal

            sage: Q*Q.transpose()           # slightly random output.
            [1.0                   5.55111512313e-17 -1.11022302463e-16]
            [ 5.55111512313e-17    1.0               -5.55111512313e-17]
            [-1.11022302463e-16    -5.55111512313e-17               1.0]
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._nrows, self._nrows), self.new_matrix()

        cdef Matrix_double_dense Q,R
        Q = self._new(self._nrows, self._nrows)
        R = self._new(self._nrows, self._ncols)
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        Q._matrix_numpy, R._matrix_numpy = scipy.linalg.qr(self._matrix_numpy)
        return Q,R

    def is_symmetric(self, tol = 1e-12):
        """
        Return whether this matrix is symmetric, to the given tolerance.

        EXAMPLES:
            sage: m = matrix(RDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: m.is_symmetric()
            False
            sage: m[1,0] = 1.1; m
            [0.0 1.0]
            [1.1 3.0]
            sage: m.is_symmetric()
            False

        The tolerance inequality is strict:
            sage: m.is_symmetric(tol=0.1)
            False
            sage: m.is_symmetric(tol=0.11)
            True
        """
        cdef Py_ssize_t i, j
        tol = float(tol)
        key = 'symmetric_%s'%tol
        b = self.fetch(key)
        if not b is None:
            return b
        if self._nrows != self._ncols:
            self.cache(key, False)
            return False
        b = True
        for i from 0 < i < self._nrows:
            for j from 0 <= j < i:
                if math.fabs(self.get_unsafe(i,j) - self.get_unsafe(j,i)) > tol:
                    b = False
                    break
        self.cache(key, b)
        return b


    def cholesky(self):
        r"""
        Return the cholesky factorization of this matrix.  The input
        matrix must be symmetric and positive definite or an exception
        will be raised.

        EXAMPLES:
            sage: M = MatrixSpace(RDF,5)
            sage: r = matrix(RDF,[[   0.,    0.,    0.,    0.,    1.],[   1.,    1.,    1.,    1.,    1.],[  16.,    8.,    4.,    2.,    1.],[  81.,   27.,    9.,    3.,    1.],[ 256.,   64.,   16.,    4.,    1.]])

            sage: m = r*M.identity_matrix()*r.transpose()
            sage: L = m.cholesky()
            sage: L*L.transpose()
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:   # special case
            return self.copy()

        cdef matrix_double_dense


        cdef Matrix_double_dense M = self._new()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        M._matrix_numpy = scipy.linalg.cholesky(self._matrix_numpy, lower=1)
        return M

    cdef Vector _vector_times_matrix_(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._row_ambient_module().zero_vector()
        global numpy
        if numpy is None:
            import numpy

        v_numpy = numpy.array([self._python_dtype(i) for i in v])

        M=self._row_ambient_module()
        ans = numpy.dot(v_numpy,self._matrix_numpy)
        return M(ans)

    cdef Vector _matrix_times_vector_(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._column_ambient_module().zero_vector()

        global numpy
        if numpy is None:
            import numpy

        v_numpy = numpy.array([self._python_dtype(i) for i in v], dtype=self._numpy_dtype)

        M=self._column_ambient_module()
        ans = numpy.dot(self._matrix_numpy, v_numpy)
        return M(ans)

    def numpy(self):
        """
        This method returns a copy of the matrix as a numpy array. It
        uses the numpy C/api so is very fast.

        EXAMPLES:
            sage: m = matrix(RDF,[[1,2],[3,4]])
            sage: n = m.numpy()
            sage: import numpy
            sage: numpy.linalg.eig(n)
            (array([-0.37228132,  5.37228132]), array([[-0.82456484, -0.41597356],
                   [ 0.56576746, -0.90937671]]))
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m.numpy()
            array([[ 0.,  1.,  2.],
                   [ 3.,  4.,  5.]])

        TESTS:
            sage: m = matrix(RDF,0,5,[]); m
            []
            sage: m.numpy()
            array([], shape=(0, 5), dtype=float64)
            sage: m = matrix(RDF,5,0,[]); m
            []
            sage: m.numpy()
            array([], shape=(5, 0), dtype=float64)
        """
        return self._matrix_numpy.copy()

    def _replace_self_with_numpy(self,numpy_matrix):
        """
        EXAMPLES:
            sage: import numpy
            sage: a = numpy.array([[1,2],[3,4]], 'float64')
            sage: m = matrix(RDF,2,2,0)
            sage: m._replace_self_with_numpy(a)
            sage: m
            [1.0 2.0]
            [3.0 4.0]
        """
        if self._matrix_numpy.shape != numpy_matrix.shape:
            raise ValueError, "matrix shapes are not the same"
        self._matrix_numpy = numpy_matrix.astype(self._numpy_dtype)


    def _replace_self_with_numpy32(self,numpy_matrix):
        """
        EXAMPLES:
            sage: import numpy
            sage: a = numpy.array([[1,2],[3,4]], 'float32')
            sage: m = matrix(RDF,2,2,0)
            sage: m._replace_self_with_numpy32(a)
            sage: m
            [1.0 2.0]
            [3.0 4.0]
        """
        #TODO find where this is used and change it
        self._replace_self_with_numpy(numpy_matrix)


    def _hadamard_row_bound(self):
        r"""
        Return an integer n such that the absolute value of the
        determinant of this matrix is at most $10^n$.

        EXAMPLES:
            sage: a = matrix(RDF, 3, [1,2,5,7,-3,4,2,1,123])
            sage: a._hadamard_row_bound()
            4
            sage: a.det()
            -2014.0
            sage: 10^4
            10000
        """
        cdef double d = 0, s
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            s = 0
            for j from 0 <= j < self._ncols:
                s += self.get_unsafe(i,j)**2
            d += math.log(s)
        d /= 2
        return int(math.ceil(d / math.log(10)))
