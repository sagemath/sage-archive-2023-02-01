"""
Dense matrices using a NumPy backend.  This serves as a base class for
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
    sage: TestSuite(a).run()
    sage: a = matrix(CDF,2,range(4), sparse=False)
    sage: TestSuite(a).run()

AUTHORS:
    -- Jason Grout, Sep 2008: switch to NumPy backend, factored out the Matrix_double_dense class
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

cimport numpy as cnumpy

numpy=None
scipy=None

# This is for the Numpy C API to work
cnumpy.import_array()


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
    #   * __cinit__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __richcmp__    -- always the same
    #   * __hash__       -- always simple
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce):
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
        new uninitialized matrix representing the data for the class.
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
            sage: A = random_matrix(RDF,3) ; A.LU_valid()
            False
            sage: P, L, U = A.LU()
            sage: A.LU_valid()
            True
        """
        return self.fetch('PLU_factors') is not None

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
        # than the built-in complex number type.
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
            return self.__copy__()

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
            return self.__copy__()

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
            return self.__copy__()

        cdef Matrix_double_dense M
        M = self._new()
        M._matrix_numpy = -self._matrix_numpy
        return M


    #   * cdef _cmp_c_impl
    # x * __copy__
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


        TESTS:
            sage: ~Matrix(RDF, 0,0)
            []
            sage: ~Matrix(RDF, 0,3)
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix
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
#            ZeroDivisionError: singular matrix
#
        if self._nrows != self._ncols:
            raise ArithmeticError, "self must be a square matrix"
        if self._nrows == 0 and self._ncols == 0:
            return self.__copy__()

        # Maybe we should cache the (P)LU decomposition and use scipy.lu_solve?
        cdef Matrix_double_dense M
        M = self._new()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        from numpy.linalg import LinAlgError
        try: ##  Standard error reporting for Sage.
            M._matrix_numpy = scipy.linalg.inv(self._matrix_numpy)
        except LinAlgError:
            raise ZeroDivisionError, "singular matrix"
        return M

    def __copy__(self):
        r"""
        Returns a new copy of this matrix.

        EXAMPLES::

            sage: a = matrix(RDF,1,3, [1,2,-3])
            sage: a
            [ 1.0  2.0 -3.0]
            sage: b = a.__copy__()
            sage: b
            [ 1.0  2.0 -3.0]
            sage: b is a
            False
            sage: b == a
            True
            sage: b[0,0] = 3
            sage: a[0,0] # note that a hasn't changed
            1.0

        ::

            sage: copy(MatrixSpace(RDF,0,0,sparse=False).zero_matrix())
            []
        """
        if self._nrows == 0 or self._ncols == 0:
            # Create a brand new empy matrix. This is needed to prevent a
            # recursive loop: a copy of zero_matrix is asked otherwise.
            return self.__class__(self.parent(), [], self._nrows, self._ncols)

        cdef Matrix_double_dense A
        A = self._new(self._nrows, self._ncols)
        A._matrix_numpy = self._matrix_numpy.copy()
        if self.subdivisions is not None:
            A.subdivide(*self.get_subdivisions())
        return A


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
    #
    ########################################################################

    def LU(self):
        """
        Computes the LU decomposition of a matrix.

        For and square matrix A we can find matrices P, L, and U such
        that ``P*A = L*U``, where P is a permutation matrix, L is lower
        triangular and U is upper triangular.

        The computed decomposition is cached and returned on subsequent calls.

        OUTPUT:
            P, L, U -- as a tuple of immutable matrices

        EXAMPLES::

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

        The result is immutable...::

            sage: L[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: P[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: U[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
        """
        global scipy, numpy
        cdef Matrix_double_dense P, L, U

        if self._ncols!=self._nrows:
            raise TypeError,"LU decomposition only works for square matrix"

        if self._ncols == 0:
            return self.__copy__(), self.__copy__(), self.__copy__()

        PLU = self.fetch('PLU_factors')
        if PLU is None:
            if scipy is None:
                import scipy
            import scipy.linalg
            if numpy is None:
                import numpy
            P = self._new()
            L = self._new()
            U = self._new()
            PM, LM, UM = scipy.linalg.lu(self._matrix_numpy)
            # Numpy has a different convention than we had with GSL
            # So we invert (transpose) the P to match our prior behavior
            # TODO: It's an awful waste to store a huge matrix for P, which
            # is just a simple permutation, really.
            P._matrix_numpy = PM.T.copy()
            L._matrix_numpy = numpy.ascontiguousarray(LM)
            U._matrix_numpy = numpy.ascontiguousarray(UM)
            PLU = (P, L, U)
            for M in PLU: M.set_immutable()
            self.cache('PLU_factors', PLU)
        return PLU

    def eigenspaces_left(self, var='a', algebraic_multiplicity=False):
        r"""
        Computes the left eigenspaces of a matrix of double precision
        real or complex numbers (i.e. RDF or CDF).

        INPUT:

        - ``var`` - ignored for numerical matrices
        - ``algebraic_multiplicity`` - must be set to ``False``
          for numerical matrices, and will raise an error otherwise.

        OUTPUT:

        Return a list of pairs ``(e, V)`` where ``e`` is a (complex)
        eigenvalue and ``V`` is the associated left eigenspace as a
        vector space.

        No attempt is made to determine if an eigenvalue has multiplicity
        greater than one, so all the eigenspaces returned have dimension one.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a sign. So for doctests we have normalized output by creating an
        eigenspace with a canonical basis.

        EXAMPLES::

            sage: m = matrix(RDF, [[-5, 3, 2, 8],[10, 2, 4, -2],[-1, -10, -10, -17],[-2, 7, 6, 13]])
            sage: spectrum = m.eigenspaces_left()
            sage: spectrum[0][0]
            2.0
            sage: (RDF^4).subspace(spectrum[0][1].basis())
            Vector space of degree 4 and dimension 1 over Real Double Field
            Basis matrix:
            [1.0 1.0 1.0 1.0]

            sage: e, V = spectrum[2]
            sage: v = V.basis()[0]
            sage: diff = (v*m).change_ring(CDF) - e*v
            sage: abs(abs(diff)) < 1e-14
            True

        TESTS::

            sage: m.eigenspaces_left(algebraic_multiplicity=True)
            Traceback (most recent call last):
            ...
            ValueError: algebraic_multiplicity must be set to False for double precision matrices
        """
        # For numerical values we leave decisions about
        # multiplicity to the calling routine
        if algebraic_multiplicity:
            raise ValueError, "algebraic_multiplicity must be set to False for double precision matrices"
        spectrum = self.left_eigenvectors()
        pairs = []
        for evalue,evectors,_ in spectrum:
            evector = evectors[0]
            espace = evector.parent().span_of_basis([evector],check=False)
            pairs.append((evalue, espace))
        return pairs

    left_eigenspaces = eigenspaces_left

    def eigenspaces_right(self, var='a', algebraic_multiplicity=False):
        r"""
        Computes the right eigenspaces of a matrix of double precision
        real or complex numbers (i.e. RDF or CDF).

        INPUT:

        - ``var`` - ignored for numerical matrices
        - ``algebraic_multiplicity`` - must be set to ``False``
          for numerical matrices, and will raise an error otherwise.

        OUTPUT:

        Return a list of pairs ``(e, V)`` where ``e`` is a (complex)
        eigenvalue and ``V`` is the associated right eigenspace as a
        vector space.

        No attempt is made to determine if an eigenvalue has multiplicity
        greater than one, so all the eigenspaces returned have dimension one.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a sign. So for doctests we have normalized output by creating an
        eigenspace with a canonical basis.


        EXAMPLES::

            sage: m = matrix(RDF, [[-9, -14, 19, -74],[-1, 2, 4, -11],[-4, -12, 6, -32],[0, -2, -1, 1]])
            sage: m
            [ -9.0 -14.0  19.0 -74.0]
            [ -1.0   2.0   4.0 -11.0]
            [ -4.0 -12.0   6.0 -32.0]
            [  0.0  -2.0  -1.0   1.0]
            sage: spectrum = m.eigenspaces_right()
            sage: spectrum[0][0]
            2.0
            sage: (RDF^4).subspace(spectrum[0][1].basis())
            Vector space of degree 4 and dimension 1 over Real Double Field
            Basis matrix:
            [ 1.0 -2.0  3.0  1.0]

            sage: e, V = spectrum[2]
            sage: v = V.basis()[0]
            sage: diff = (m*v).change_ring(CDF) - e*v
            sage: abs(abs(diff)) < 1e-14
            True

        TESTS::

            sage: m.eigenspaces_right(algebraic_multiplicity=True)
            Traceback (most recent call last):
            ...
            ValueError: algebraic_multiplicity must be set to False for double precision matrices
        """
        # For numerical values we leave decisions about
        # multiplicity to the calling routine
        if algebraic_multiplicity:
            raise ValueError, "algebraic_multiplicity must be set to False for double precision matrices"
        spectrum = self.right_eigenvectors()
        pairs = []
        for evalue,evectors,_ in spectrum:
            evector = evectors[0]
            espace = evector.parent().span_of_basis([evector],check=False)
            pairs.append((evalue, espace))
        return pairs

    right_eigenspaces = eigenspaces_right

    def eigenvalues(self):
        r"""
        Returns a list of the eigenvalues (with multiplicity)
        of this matrix.  The returned eigenvalues are elements of CDF.

        EXAMPLES:
            sage: m = matrix(RDF, 2, 2, [1,2,3,4])
            sage: m.eigenvalues()
            [-0.372281323269, 5.37228132327]
            sage: parent(m.eigenvalues()[0])
            Complex Double Field

            sage: m = matrix(RDF, 2, 2, [0,1,-1,0])
            sage: m.eigenvalues()
            [1.0*I, -1.0*I]

            sage: m = matrix(CDF, 2, 2, [I,1,-I,0])
            sage: m.eigenvalues()
            [-0.624810533844 + 1.30024259022*I, 0.624810533844 - 0.30024259022*I]

            sage: matrix(CDF,0,0).eigenvalues()
            []
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        if self._nrows == 0:
            return []
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        return [sage.rings.complex_double.CDF(x) for x in scipy.linalg.eigvals(self._matrix_numpy)]

    def left_eigenvectors(self):
        r"""
        Compute the left eigenvectors of a matrix of double precision
        real or complex numbers (i.e. RDF or CDF).

        OUTPUT:
        Returns a list of triples, each of the form ``(e,[v],1)``,
        where ``e`` is the eigenvalue, and ``v`` is an associated
        left eigenvector.  If the matrix is of size `n`, then there are
        `n` triples.  Values are computed with the SciPy library.

        The format of this output is designed to match the format
        for exact results.  However, since matrices here have numerical
        entries, the resulting eigenvalues will also be numerical.  No
        attempt is made to determine if two eigenvalues are equal, or if
        eigenvalues might actually be zero.  So the algebraic multiplicity
        of each eigenvalue is reported as 1.  Decisions about equal
        eigenvalues or zero eigenvalues should be addressed in the
        calling routine.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a sign. So for doctests we have normalized output by forcing their
        eigenvectors to have their first non-zero entry equal to one.

        EXAMPLES::

            sage: m = matrix(RDF, [[-5, 3, 2, 8],[10, 2, 4, -2],[-1, -10, -10, -17],[-2, 7, 6, 13]])
            sage: m
            [ -5.0   3.0   2.0   8.0]
            [ 10.0   2.0   4.0  -2.0]
            [ -1.0 -10.0 -10.0 -17.0]
            [ -2.0   7.0   6.0  13.0]
            sage: spectrum = m.left_eigenvectors()
            sage: for i in range(len(spectrum)):
            ...     spectrum[i][1][0]=matrix(RDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum[0]
            (2.0, [(1.0, 1.0, 1.0, 1.0)], 1)
            sage: spectrum[1]
            (1.0, [(1.0, 0.8, 0.8, 0.6)], 1)
            sage: spectrum[2]
            (-2.0, [(1.0, 0.4, 0.6, 0.2)], 1)
            sage: spectrum[3]
            (-1.0, [(1.0, 1.0, 2.0, 2.0)], 1)
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        if self._nrows == 0:
            return [], self.__copy__()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v,eig = scipy.linalg.eig(self._matrix_numpy, right=False, left=True)
        # scipy puts eigenvectors in columns, we will extract from rows
        eig = matrix(eig.T)
        return [(sage.rings.complex_double.CDF(v[i]), [eig[i]], 1) for i in range(len(v))]

    eigenvectors_left = left_eigenvectors

    def right_eigenvectors(self):
        r"""
        Compute the right eigenvectors of a matrix of double precision
        real or complex numbers (i.e. RDF or CDF).

        OUTPUT:
        Returns a list of triples, each of the form ``(e,[v],1)``,
        where ``e`` is the eigenvalue, and ``v`` is an associated
        right eigenvector.  If the matrix is of size `n`, then there
        are `n` triples.  Values are computed with the SciPy library.

        The format of this output is designed to match the format
        for exact results.  However, since matrices here have numerical
        entries, the resulting eigenvalues will also be numerical.  No
        attempt is made to determine if two eigenvalues are equal, or if
        eigenvalues might actually be zero.  So the algebraic multiplicity
        of each eigenvalue is reported as 1.  Decisions about equal
        eigenvalues or zero eigenvalues should be addressed in the
        calling routine.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a sign. So for doctests we have normalized output by forcing their
        eigenvectors to have their first non-zero entry equal to one.

        EXAMPLES::

            sage: m = matrix(RDF, [[-9, -14, 19, -74],[-1, 2, 4, -11],[-4, -12, 6, -32],[0, -2, -1, 1]])
            sage: m
            [ -9.0 -14.0  19.0 -74.0]
            [ -1.0   2.0   4.0 -11.0]
            [ -4.0 -12.0   6.0 -32.0]
            [  0.0  -2.0  -1.0   1.0]
            sage: spectrum = m.right_eigenvectors()
            sage: for i in range(len(spectrum)):
            ...     spectrum[i][1][0]=matrix(RDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum[0]
            (2.0, [(1.0, -2.0, 3.0, 1.0)], 1)
            sage: spectrum[1]
            (1.0, [(1.0, -0.666666666667, 1.33333333333, 0.333333333333)], 1)
            sage: spectrum[2]
            (-2.0, [(1.0, -0.2, 1.0, 0.2)], 1)
            sage: spectrum[3]
            (-1.0, [(1.0, -0.5, 2.0, 0.5)], 1)
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        if self._nrows == 0:
            return [], self.__copy__()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v,eig = scipy.linalg.eig(self._matrix_numpy, right=True, left=False)
        # scipy puts eigenvectors in columns, we will extract from rows
        eig = matrix(eig.T)
        return [(sage.rings.complex_double.CDF(v[i]), [eig[i]], 1) for i in range(len(v))]

    eigenvectors_right = right_eigenvectors

    def solve_left_LU(self, b):
        """
        Solve the equation A*x = b using LU decomposition.

        WARNING: This function is broken. See trac 4932.

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
            sage: x = A.solve_left_LU(b); x
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is not finished (see trac 4932)


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

        raise NotImplementedError, "this function is not finished (see trac 4932)"
        self._c_compute_LU()  # so self._L_M and self._U_M are defined below.
        cdef Matrix_double_dense M = self._new()
        lu = self._L_M*self._U_M
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        M._matrix_numpy = scipy.linalg.lu_solve((lu, self._P_M), b)
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
            ValueError: self must be a square matrix
        """
        if not self.is_square():
            raise ValueError, "self must be a square matrix"
        if self._nrows == 0 or self._ncols == 0:
            return self._sage_dtype(1)
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        return self._sage_dtype(scipy.linalg.det(self._matrix_numpy))


    def log_determinant(self):
        """
        Compute the log of the absolute value of the determinant
        using LU decomposition.

        NOTE: This is useful if the usual determinant overflows.

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
        global numpy
        cdef Matrix_double_dense U

        if self._nrows == 0 or self._ncols == 0:
            return sage.rings.real_double.RDF(0)

        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"

        P, L, U = self.LU()
        if numpy is None:
            import numpy

        return sage.rings.real_double.RDF(sum(numpy.log(abs(numpy.diag(U._matrix_numpy)))))

    def transpose(self):
        """
        Return the transpose of this matrix, without changing self.

        EXAMPLES::

            sage: m = matrix(RDF,2,3,range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m2 = m.transpose()
            sage: m[0,0] = 2
            sage: m2           #note that m2 hasn't changed
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
        trans._matrix_numpy = self._matrix_numpy.transpose().copy()
        if self.subdivisions is not None:
            row_divs, col_divs = self.get_subdivisions()
            trans.subdivide(col_divs, row_divs)
        return trans

    def SVD(self, *args, **kwds):
        r"""
        Return the singular value decomposition of this matrix.

        The U and V matrices are not unique and may be returned with different
        values in the future or on different systems. The S matrix is unique
        and contains the singular values in descending order.

        The computed decomposition is cached and returned on subsequent calls.

        INPUT:
            A -- a matrix

        OUTPUT:
            U, S, V -- immutable matrices such that $A = U*S*V.conj().transpose()$
                       where U and V are orthogonal and S is zero off of the diagonal.

        Note that if self is m-by-n, then the dimensions of the
        matrices that this returns are (m,m), (m,n), and (n, n).

        EXAMPLES::

            sage: m = matrix(RDF,4,range(1,17))
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()
            [ 1.0  2.0  3.0  4.0]
            [ 5.0  6.0  7.0  8.0]
            [ 9.0 10.0 11.0 12.0]
            [13.0 14.0 15.0 16.0]

        A non-square example::

            sage: m = matrix(RDF, 2, range(1,7)); m
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            sage: U, S, V = m.SVD()
            sage: U*S*V.transpose()
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]

        S contains the singular values::

            sage: S.round(4)
            [ 9.508    0.0    0.0]
            [   0.0 0.7729    0.0]
            sage: [round(sqrt(abs(x)),4) for x in (S*S.transpose()).eigenvalues()]
            [9.508, 0.7729]

        U and V are orthogonal matrices::

            sage: U # random, SVD is not unique
            [-0.386317703119 -0.922365780077]
            [-0.922365780077  0.386317703119]
            [-0.274721127897 -0.961523947641]
            [-0.961523947641  0.274721127897]
            sage: (U*U.transpose()).zero_at(1e-15)
            [1.0 0.0]
            [0.0 1.0]
            sage: V # random, SVD is not unique
            [-0.428667133549  0.805963908589  0.408248290464]
            [-0.566306918848  0.112382414097 -0.816496580928]
            [-0.703946704147 -0.581199080396  0.408248290464]
            sage: (V*V.transpose()).zero_at(1e-15)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]

        TESTS::

            sage: m = matrix(RDF,3,2,range(1, 7)); m
            [1.0 2.0]
            [3.0 4.0]
            [5.0 6.0]
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()
            [1.0 2.0]
            [3.0 4.0]
            [5.0 6.0]

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
            sage: for x in m.SVD(): x.is_immutable()
            True
            True
            True
        """
        global scipy, numpy
        cdef Py_ssize_t i
        cdef Matrix_double_dense U, S, V

        if len(args)>0 or len(kwds)>0:
            from sage.misc.misc import deprecation
            deprecation("Arguments passed to SVD, but SVD no longer supports different methods (it only uses numpy now).")

        if self._nrows == 0 or self._ncols == 0:
            U_t = self.new_matrix(self._nrows, self._ncols)
            S_t = self.new_matrix(self._nrows, self._ncols)
            V_t = self.new_matrix(self._ncols, self._nrows)
            return U_t, S_t, V_t

        USV = self.fetch('SVD_factors')
        if USV is None:
            # TODO: More efficient representation of non-square diagonal matrix S
            if scipy is None:
                import scipy
            import scipy.linalg
            if numpy is None:
                import numpy
            U_mat, S_diagonal, V_mat = scipy.linalg.svd(self._matrix_numpy)

            U = self._new(self._nrows, self._nrows)
            S = self._new(self._nrows, self._ncols)
            V = self._new(self._ncols, self._ncols)

            S_mat = numpy.zeros((self._nrows, self._ncols), dtype=self._numpy_dtype)
            for i in range(S_diagonal.shape[0]):
                S_mat[i,i] = S_diagonal[i]

            U._matrix_numpy = numpy.ascontiguousarray(U_mat)
            S._matrix_numpy = S_mat
            V._matrix_numpy = numpy.ascontiguousarray(V_mat.conj().T)
            USV = U, S, V
            for M in USV: M.set_immutable()
            self.cache('SVD_factors', USV)

        return USV

    def QR(self):
        """
        Return the Q,R factorization of a real matrix A.

        The computed decomposition is cached and returned on subsequent calls.

        INPUT:
           self -- a real matrix A

        OUTPUT:
           Q, R -- immutable matrices such that A = Q*R such that the columns of Q are
                   orthogonal (i.e., $Q^t Q = I$), and R is upper triangular.

        EXAMPLES::

            sage: m = matrix(RDF,3,range(0, 12)); m
            [ 0.0  1.0  2.0  3.0]
            [ 4.0  5.0  6.0  7.0]
            [ 8.0  9.0 10.0 11.0]
            sage: Q,R = m.QR()
            sage: Q*R
            [ 0.0  1.0  2.0  3.0]
            [ 4.0  5.0  6.0  7.0]
            [ 8.0  9.0 10.0 11.0]

        Note that Q is an orthogonal matrix::

            sage: (Q*Q.transpose()).zero_at(1e-10)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]

        The result is immutable::

            sage: Q[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: R.is_immutable()
            True

        """
        global scipy
        cdef Matrix_double_dense Q,R

        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._nrows, self._nrows), self.new_matrix()

        QR = self.fetch('QR_factors')
        if QR is None:
            Q = self._new(self._nrows, self._nrows)
            R = self._new(self._nrows, self._ncols)
            if scipy is None:
                import scipy
            import scipy.linalg
            Q._matrix_numpy, R._matrix_numpy = scipy.linalg.qr(self._matrix_numpy)
            Q.set_immutable()
            R.set_immutable()
            QR = (Q, R)
            self.cache('QR_factors', QR)
        return QR

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

        WARNING::

        .. This is deprecated!  Use ``cholesky_decomposition'' instead.

        EXAMPLES::

            sage: M = MatrixSpace(RDF,5)
            sage: r = matrix(RDF,[[   0.,    0.,    0.,    0.,    1.],[   1.,    1.,    1.,    1.,    1.],[  16.,    8.,    4.,    2.,    1.],[  81.,   27.,    9.,    3.,    1.],[ 256.,   64.,   16.,    4.,    1.]])

            sage: m = r*M.identity_matrix()*r.transpose()
            sage: L = m.cholesky()
            doctest... DeprecationWarning: cholesky is deprecated for matrices over RDF; use cholesky_decomposition instead.
            sage: L*L.transpose()
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]
        """
        # deprecation added 2009-05
        from sage.misc.misc import deprecation
        deprecation("cholesky is deprecated for matrices over RDF; use cholesky_decomposition instead.")
        return self.cholesky_decomposition()

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

        Alternatively, numpy automatically calls this function (via
        the magic :meth:`__array__` method) to convert Sage matrices
        to numpy arrays::

            sage: import numpy
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: numpy.array(m)
            array([[ 0.,  1.,  2.],
                   [ 3.,  4.,  5.]])
            sage: numpy.array(m).dtype
            dtype('float64')
            sage: m = matrix(CDF, 2, range(6)); m
            [  0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: numpy.array(m)
            array([[ 0.+0.j,  1.+0.j,  2.+0.j],
                   [ 3.+0.j,  4.+0.j,  5.+0.j]])
            sage: numpy.array(m).dtype
            dtype('complex128')

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
        if (<object>self._matrix_numpy).shape != (<object>numpy_matrix).shape:
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

    def exp(self, method='pade', order=None):
        r"""
        Calculate the exponential of this matrix X, which is the matrix

        e^X = sum_{k=0}^{infty} frac{X^k}{k!}.

        INPUT:

            method -- 'pade', 'eig', or 'taylor'; the method used to
            compute the exponential.

            order -- for Pade or Taylor series methods, order is the
            order of the Pade approximation or the order of the Taylor
            series used.  The current defaults (from scipy) are 7 for
            'pade' and 20 for 'taylor'.

        EXAMPLES:
            sage: A=matrix(RDF, 2, [1,2,3,4]); A
            [1.0 2.0]
            [3.0 4.0]
            sage: A.exp()
            [51.9689561987  74.736564567]
            [112.104846851 164.073803049]
            sage: A.exp(method='eig')
            [51.9689561987  74.736564567]
            [112.104846851 164.073803049]
            sage: A.exp(order=2)
            [51.8888631634 74.6198348038]
            [111.929752206 163.818615369]
            sage: A.exp(method='taylor', order=5)
            [19.9583333333 28.0833333333]
            [       42.125 62.0833333333]
            sage: A.exp(method='taylor')
            [51.9689035511 74.7364878369]
            [112.104731755 164.073635306]

            sage: A=matrix(CDF, 2, [1,2+I,3*I,4]); A
            [        1.0 2.0 + 1.0*I]
            [      3.0*I         4.0]
            sage: A.exp()
            [-19.6146029538 + 12.5177438468*I  3.79496364496 + 28.8837993066*I]
            [-32.3835809809 + 21.8842359579*I  2.26963300409 + 44.9013248277*I]
            sage: A.exp(method='eig')
            [-19.6146029538 + 12.5177438468*I  3.79496364496 + 28.8837993066*I]
            [-32.3835809809 + 21.8842359579*I  2.26963300409 + 44.9013248277*I]
            sage: A.exp(order=2)
            [-19.6130852955 + 12.5327938535*I   3.81156364812 + 28.891438232*I]
            [-32.3827876895 + 21.9087393169*I   2.29565402142 + 44.915581543*I]
            sage: A.exp(method='taylor', order=5)
            [       -6.29166666667 + 14.25*I 14.0833333333 + 15.7916666667*I]
            [               -10.5 + 26.375*I         20.0833333333 + 24.75*I]
            sage: A.exp(method='taylor')
            [-19.6146006163 + 12.5177432169*I  3.79496442472 + 28.8837964828*I]
            [-32.3835771246 + 21.8842351994*I  2.26963458304 + 44.9013203415*I]
        """
        if method not in ('pade', 'eig', 'taylor'):
            raise ValueError, "method must be 'pade', 'eig', or 'taylor'"

        if scipy is None:
            import scipy
        import scipy.linalg

        cdef Matrix_double_dense M
        M = self._new()

        if method=='pade':
            if order is None:
                M._matrix_numpy = scipy.linalg.expm(self._matrix_numpy)
            else:
                M._matrix_numpy = scipy.linalg.expm(self._matrix_numpy, q=order)
        elif method=='eig':
            M._matrix_numpy = scipy.linalg.expm2(self._matrix_numpy)
        elif method=='taylor':
            if order is None:
                M._matrix_numpy = scipy.linalg.expm3(self._matrix_numpy)
            else:
                M._matrix_numpy = scipy.linalg.expm3(self._matrix_numpy, q=order)

        return M

    def zero_at(self, eps):
        """
        Returns a copy of the matrix where elements smaller than or
        equal to ``eps`` are replaced with zeroes. For complex matrices,
        the real and imaginary parts are considered individually.

        This is useful for modifying output from algorithms which have large
        relative errors when producing zero elements, e.g. to create reliable
        doctests.

        INPUT:

          - ``eps`` - Cutoff value

        OUTPUT:

        A modified copy of the matrix.

        EXAMPLES::

            sage: a=matrix([[1, 1e-4r, 1+1e-100jr], [1e-8+3j, 0, 1e-58r]])
            sage: a
            [           1.0         0.0001 1.0 + 1e-100*I]
            [ 1e-08 + 3.0*I              0          1e-58]
            sage: a.zero_at(1e-50)
            [          1.0        0.0001           1.0]
            [1e-08 + 3.0*I             0             0]
            sage: a.zero_at(1e-4)
            [  1.0     0   1.0]
            [3.0*I     0     0]



        """
        global numpy
        cdef Matrix_double_dense M
        if numpy is None:
            import numpy
        eps = float(eps)
        out = self._matrix_numpy.copy()
        if self._sage_dtype is sage.rings.complex_double.CDF:
            out.real[numpy.abs(out.real) <= eps] = 0
            out.imag[numpy.abs(out.imag) <= eps] = 0
        else:
            out[numpy.abs(out) <= eps] = 0
        M = self._new()
        M._matrix_numpy = out
        return M

    def round(self, ndigits=0):
        """
        Returns a copy of the matrix where all entries have been rounded
        to a given precision in decimal digits (default 0 digits).

        INPUT:

         - ``ndigits`` - The precision in number of decimal digits

        OUTPUT:

        A modified copy of the matrix

        EXAMPLES:

            sage: M=matrix(CDF, [[10.234r + 34.2343jr, 34e10r]])
            sage: M
            [10.234 + 34.2343*I     340000000000.0]
            sage: M.round(2)
            [10.23 + 34.23*I  340000000000.0]
            sage: M.round()
            [ 10.0 + 34.0*I 340000000000.0]
        """
        global numpy
        cdef Matrix_double_dense M
        if numpy is None:
            import numpy
        ndigits = int(ndigits)
        M = self._new()
        M._matrix_numpy = numpy.round(self._matrix_numpy, ndigits)
        return M
