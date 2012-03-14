"""
Dense matrices using a NumPy backend.


This serves as a base class for dense matrices over
Real Double Field and Complex Double Field.

AUTHORS:

- Jason Grout, Sep 2008: switch to NumPy backend, factored out the Matrix_double_dense class

- Josh Kantor

- William Stein: many bug fixes and touch ups.


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
    sage: TestSuite(a).run()
    sage: a = matrix(CDF,2,range(4), sparse=False)
    sage: TestSuite(a).run()
"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math

import sage.rings.real_double
import sage.rings.complex_double

from matrix cimport Matrix
from sage.structure.element cimport ModuleElement,Vector
from constructor import matrix
from sage.modules.free_module_element import vector
cimport sage.structure.element
from matrix_space import MatrixSpace
from sage.misc.decorators import rename_keyword

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

    EXAMPLES::

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
        new uninitialized matrix representing the data for the class.::

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

        EXAMPLES::

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
        Returns ``True`` if the LU form of this matrix has
        already been computed.

        EXAMPLES::

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

        EXAMPLES::

            sage: matrix(RDF,3,range(9))
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            [6.0 7.0 8.0]
            sage: matrix(CDF,3,3,2)
            [2.0   0   0]
            [  0 2.0   0]
            [  0   0 2.0]

        TESTS::

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
                    raise TypeError("entries has wrong length")

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
                    raise TypeError("entries must be coercible to a list or float")
            if z != 0:
                if self._nrows != self._ncols:
                    raise TypeError("scalar matrix must be square")
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

        INPUT:

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

        EXAMPLES::

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


        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

            sage: A = matrix(RDF,3,range(1,10))
            sage: B = matrix(RDF,3,range(1,13))
            sage: A*B
            [ 38.0  44.0  50.0  56.0]
            [ 83.0  98.0 113.0 128.0]
            [128.0 152.0 176.0 200.0]
        """
        if self._ncols!=right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right")

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

        EXAMPLES::

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
        and tuning options given to ATLAS::

            sage: A = matrix(RDF,3,range(1,10));A
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            [7.0 8.0 9.0]

            sage: A.determinant() < 10e-12
            True


        TESTS::

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
            raise ArithmeticError("self must be a square matrix")
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
            raise ZeroDivisionError("input matrix must be nonsingular")
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
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
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


    def condition(self, p='frob'):
        r"""
        Returns the condition number of a square nonsingular matrix.

        Roughly speaking, this is a measure of how sensitive
        the matrix is to round-off errors in numerical computations.
        The minimum possible value is 1.0, and larger numbers indicate
        greater sensitivity.

        INPUT:

        - ``p`` - default: 'frob' - controls which norm is used
          to compute the condition number, allowable values are
          'frob' (for the Frobenius norm), integers -2, -1, 1, 2,
          positive and negative infinity. See output discussion
          for specifics.

        OUTPUT:

        The condition number of a matrix is the product of a norm
        of the matrix times the norm of the inverse of the matrix.
        This requires that the matrix be square and invertible
        (nonsingular, full rank).

        Returned value is a double precision floating point value
        in ``RDF``, or ``Infinity``.  Row and column sums described below are
        sums of the absolute values of the entries, where the
        absolute value of the complex number `a+bi` is `\sqrt{a^2+b^2}`.
        Singular values are the "diagonal" entries of the "S" matrix in
        the singular value decomposition.

        - ``p = 'frob'``: the default norm employed in computing
          the condition number, the Frobenius norm, which for a
          matrix `A=(a_{ij})` computes

          .. math::

                \left(\sum_{i,j}\left\lvert{a_{i,j}}\right\rvert^2\right)^{1/2}

        - ``p = 'sv'``: the quotient of the maximal and minimal singular value.
        - ``p = Infinity`` or ``p = oo``: the maximum row sum.
        - ``p = -Infinity`` or ``p = -oo``: the minimum column sum.
        - ``p = 1``: the maximum column sum.
        - ``p = -1``: the minimum column sum.
        - ``p = 2``: the 2-norm, equal to the maximum singular value.
        - ``p = -2``: the minimum singular value.

        ALGORITHM:

        Computation is performed by the ``cond()`` function of
        the SciPy/NumPy library.

        EXAMPLES:

        First over the reals.  ::

            sage: A = matrix(RDF, 4, [(1/4)*x^3 for x in range(16)]); A
            [   0.0   0.25    2.0   6.75]
            [  16.0  31.25   54.0  85.75]
            [ 128.0 182.25  250.0 332.75]
            [ 432.0 549.25  686.0 843.75]
            sage: A.condition()
            9923.88955...
            sage: A.condition(p='frob')
            9923.88955...
            sage: A.condition(p=Infinity)
            22738.5
            sage: A.condition(p=-Infinity)
            17.5
            sage: A.condition(p=1)
            12139.21...
            sage: A.condition(p=-1)
            550.0
            sage: A.condition(p=2)
            9897.8088...
            sage: A.condition(p=-2)
            0.000101032462...

        And over the complex numbers.  ::

            sage: B = matrix(CDF, 3, [x + x^2*I for x in range(9)]); B
            [           0  1.0 + 1.0*I  2.0 + 4.0*I]
            [ 3.0 + 9.0*I 4.0 + 16.0*I 5.0 + 25.0*I]
            [6.0 + 36.0*I 7.0 + 49.0*I 8.0 + 64.0*I]
            sage: B.condition()
            203.851798...
            sage: B.condition(p='frob')
            203.851798...
            sage: B.condition(p=Infinity)
            369.55630...
            sage: B.condition(p=-Infinity)
            5.46112969...
            sage: B.condition(p=1)
            289.251481...
            sage: B.condition(p=-1)
            20.4566639...
            sage: B.condition(p=2)
            202.653543...
            sage: B.condition(p=-2)
            0.00493453005...

        Hilbert matrices are famously ill-conditioned, while
        an identity matrix can hit the minimum with the right norm.  ::

            sage: A = matrix(RDF, 10, [1/(i+j+1) for i in range(10) for j in range(10)])
            sage: A.condition()
            1.633...e+13
            sage: id = identity_matrix(CDF, 10)
            sage: id.condition(p=1)
            1.0

        Return values are in `RDF`.  ::

            sage: A = matrix(CDF, 2, range(1,5))
            sage: A.condition() in RDF
            True

        Rectangular and singular matrices raise errors if p is not 'sv'.  ::

            sage: A = matrix(RDF, 2, 3, range(6))
            sage: A.condition()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square if p is not 'sv', not 2 x 3

            sage: A.condition('sv')
            7.34...

            sage: A = matrix(QQ, 5, range(25))
            sage: A.is_singular()
            True
            sage: B = A.change_ring(CDF)
            sage: B.condition()
            Traceback (most recent call last):
            ...
            LinAlgError: Singular matrix

        Improper values of ``p`` are caught.  ::

            sage: A = matrix(CDF, 2, range(1,5))
            sage: A.condition(p='bogus')
            Traceback (most recent call last):
            ...
            ValueError: condition number 'p' must be +/- infinity, 'frob', 'sv' or an integer, not bogus
            sage: A.condition(p=632)
            Traceback (most recent call last):
            ...
            ValueError: condition number integer values of 'p' must be -2, -1, 1 or 2, not 632

        TESTS:

        Some condition numbers, first by the definition which also exercises
        :meth:`norm`, then by this method.  ::

            sage: A = matrix(CDF, [[1,2,4],[5,3,9],[7,8,6]])
            sage: c = A.norm(2)*A.inverse().norm(2)
            sage: d = A.condition(2)
            sage: abs(c-d) < 1.0e-12
            True
            sage: c = A.norm(1)*A.inverse().norm(1)
            sage: d = A.condition(1)
            sage: abs(c-d) < 1.0e-12
            True
        """
        if not self.is_square() and p != 'sv':
            raise TypeError("matrix must be square if p is not 'sv', not %s x %s" % (self.nrows(), self.ncols()))
        global numpy
        if numpy is None:
            import numpy
        import sage.rings.infinity
        import sage.rings.integer
        import sage.rings.real_double
        if p == sage.rings.infinity.Infinity:
            p = numpy.inf
        elif p == -sage.rings.infinity.Infinity:
            p = -numpy.inf
        elif p == 'frob':
            p = 'fro'
        elif p == 'sv' :
            p = None
        else:
            try:
                p = sage.rings.integer.Integer(p)
            except:
                raise ValueError("condition number 'p' must be +/- infinity, 'frob', 'sv' or an integer, not %s" % p)
            if p not in [-2,-1,1,2]:
                raise ValueError("condition number integer values of 'p' must be -2, -1, 1 or 2, not %s" % p)
        # may raise a LinAlgError if matrix is singular
        c = numpy.linalg.cond(self._matrix_numpy, p=p)
        if c == numpy.inf:
            return sage.rings.infinity.Infinity
        else:
            return sage.rings.real_double.RDF(c)

    def norm(self, p='frob'):
        r"""
        Returns the norm of the matrix.

        INPUT:

        - ``p`` - default: 'frob' - controls which norm is computed,
          allowable values are 'frob' (for the Frobenius norm),
          integers -2, -1, 1, 2, positive and negative infinity.
          See output discussion for specifics.

        OUTPUT:

        Returned value is a double precision floating point value
        in ``RDF``.  Row and column sums described below are
        sums of the absolute values of the entries, where the
        absolute value of the complex number `a+bi` is `\sqrt{a^2+b^2}`.
        Singular values are the "diagonal" entries of the "S" matrix in
        the singular value decomposition.

        - ``p = 'frob'``: the default, the Frobenius norm, which for
          a matrix `A=(a_{ij})` computes

          .. math::

                \left(\sum_{i,j}\left\lvert{a_{i,j}}\right\rvert^2\right)^{1/2}

        - ``p = Infinity`` or ``p = oo``: the maximum row sum.
        - ``p = -Infinity`` or ``p = -oo``: the minimum column sum.
        - ``p = 1``: the maximum column sum.
        - ``p = -1``: the minimum column sum.
        - ``p = 2``: the 2-norm, equal to the maximum singular value.
        - ``p = -2``: the minimum singular value.

        ALGORITHM:

        Computation is performed by the ``norm()`` function of
        the SciPy/NumPy library.

        EXAMPLES:

        First over the reals.  ::

            sage: A = matrix(RDF, 3, range(-3, 6)); A
            [-3.0 -2.0 -1.0]
            [ 0.0  1.0  2.0]
            [ 3.0  4.0  5.0]
            sage: A.norm()
            8.30662386...
            sage: A.norm(p='frob')
            8.30662386...
            sage: A.norm(p=Infinity)
            12.0
            sage: A.norm(p=-Infinity)
            3.0
            sage: A.norm(p=1)
            8.0
            sage: A.norm(p=-1)
            6.0
            sage: A.norm(p=2)
            7.99575670...
            sage: A.norm(p=-2) < 10^-15
            True

        And over the complex numbers.  ::

            sage: B = matrix(CDF, 2, [[1+I, 2+3*I],[3+4*I,3*I]]); B
            [1.0 + 1.0*I 2.0 + 3.0*I]
            [3.0 + 4.0*I       3.0*I]
            sage: B.norm()
            7.0
            sage: B.norm(p='frob')
            7.0
            sage: B.norm(p=Infinity)
            8.0
            sage: B.norm(p=-Infinity)
            5.01976483...
            sage: B.norm(p=1)
            6.60555127...
            sage: B.norm(p=-1)
            6.41421356...
            sage: B.norm(p=2)
            6.66189877...
            sage: B.norm(p=-2)
            2.14921023...

        Since it is invariant under unitary multiplication, the
        Frobenius norm is equal to the square root of the sum of
        squares of the singular values.  ::

            sage: A = matrix(RDF, 5, range(1,26))
            sage: f = A.norm(p='frob')
            sage: U, S, V = A.SVD()
            sage: s = sqrt(sum([S[i,i]^2 for i in range(5)]))
            sage: abs(f-s) < 1.0e-12
            True

        Return values are in `RDF`. ::

            sage: A = matrix(CDF, 2, range(4))
            sage: A.norm() in RDF
            True

        Improper values of ``p`` are caught.  ::

            sage: A.norm(p='bogus')
            Traceback (most recent call last):
            ...
            ValueError: matrix norm 'p' must be +/- infinity, 'frob' or an integer, not bogus
            sage: A.norm(p=632)
            Traceback (most recent call last):
            ...
            ValueError: matrix norm integer values of 'p' must be -2, -1, 1 or 2, not 632
        """
        global numpy
        if numpy is None:
            import numpy
        import sage.rings.infinity
        import sage.rings.integer
        import sage.rings.real_double
        if p == sage.rings.infinity.Infinity:
            p = numpy.inf
        elif p == -sage.rings.infinity.Infinity:
            p = -numpy.inf
        elif p == 'frob':
            p = 'fro'
        else:
            try:
                p = sage.rings.integer.Integer(p)
            except:
                raise ValueError("matrix norm 'p' must be +/- infinity, 'frob' or an integer, not %s" % p)
            if not p in [-2,-1,1,2]:
                raise ValueError("matrix norm integer values of 'p' must be -2, -1, 1 or 2, not %s" % p)
        return sage.rings.real_double.RDF(numpy.linalg.norm(self._matrix_numpy, ord=p))

    def singular_values(self, eps=None):
        r"""
        Returns a sorted list of the singular values of the matrix.

        INPUT:

        - ``eps`` - default: ``None`` - the largest number which
          will be considered to be zero.  May also be set to the
          string 'auto'.  See the discussion below.

        OUTPUT:

        A sorted list of the singular values of the matrix, which are the
        diagonal entries of the "S" matrix in the SVD decomposition.  As such,
        the values are real and are returned as elements of ``RDF``.  The
        list is sorted with larger values first, and since theory predicts
        these values are always positive, for a rank-deficient matrix the
        list should end in zeros (but in practice may not).  The length of
        the list is the minimum of the row count and column count for the
        matrix.

        The number of non-zero singular values will be the rank of the
        matrix.  However, as a numerical matrix, it is impossible to
        control the difference between zero entries and very small
        non-zero entries.  As an informed consumer it is up to you
        to use the output responsibly.  We will do our best, and give
        you the tools to work with the output, but we cannot
        give you a guarantee.

        With ``eps`` set to ``None`` you will get the raw singular
        values and can manage them as you see fit.  You may also set
        ``eps`` to any positive floating point value you wish.  If you
        set ``eps`` to 'auto' this routine will compute a reasonable
        cutoff value, based on the size of the matrix, the largest
        singular value and the smallest nonzero value representable
        by the 53-bit precision values used.  See the discussion
        at page 268 of [WATKINS]_.

        See the examples for a way to use the "verbose" facility
        to easily watch the zero cutoffs in action.

        ALGORITHM:

        The singular values come from the SVD decomposition
        computed by SciPy/NumPy.

        EXAMPLES:

        Singular values close to zero have trailing digits that may vary
        on different hardware.  For exact matrices, the number of non-zero
        singular values will equal the rank of the matrix.  So for some of
        the doctests we round the small singular values that ideally would
        be zero, to control the variability across hardware.

        This matrix has a determinant of one.  A chain of two or
        three theorems implies the product of the singular values
        must also be one.  ::

            sage: A = matrix(QQ, [[ 1,  0,  0,  0,  0,  1,  3],
            ...                   [-2,  1,  1, -2,  0, -4,  0],
            ...                   [ 1,  0,  1, -4, -6, -3,  7],
            ...                   [-2,  2,  1,  1,  7,  1, -1],
            ...                   [-1,  0, -1,  5,  8,  4, -6],
            ...                   [ 4, -2, -2,  1, -3,  0,  8],
            ...                   [-2,  1,  0,  2,  7,  3, -4]])
            sage: A.determinant()
            1
            sage: B = A.change_ring(RDF)
            sage: sv = B.singular_values(); sv
            [20.5239806589, 8.48683702854, 5.86168134845, 2.44291658993,
              0.583197014472, 0.269332872866, 0.00255244880761]
            sage: prod(sv)
            1.0

        An exact matrix that is obviously not of full rank, and then
        a computation of the singular values after conversion
        to an approximate matrix. ::

            sage: A = matrix(QQ, [[1/3, 2/3, 11/3],
            ...                   [2/3, 1/3,  7/3],
            ...                   [2/3, 5/3, 27/3]])
            sage: A.rank()
            2
            sage: B = A.change_ring(CDF)
            sage: sv = B.singular_values()
            sage: sv[0:2]
            [10.1973039..., 0.487045871...]
            sage: sv[2] < 1e-14
            True

        A matrix of rank 3 over the complex numbers.  ::

            sage: A = matrix(CDF, [[46*I - 28, -47*I - 50, 21*I + 51, -62*I - 782, 13*I + 22],
            ...                    [35*I - 20, -32*I - 46, 18*I + 43, -57*I - 670, 7*I + 3],
            ...                    [22*I - 13, -23*I - 23, 9*I + 24, -26*I - 347, 7*I + 13],
            ...                    [-44*I + 23, 41*I + 57, -19*I - 54, 60*I + 757, -11*I - 9],
            ...                    [30*I - 18, -30*I - 34, 14*I + 34, -42*I - 522, 8*I + 12]])
            sage: sv = A.singular_values()
            sage: sv[0:3]
            [1440.733666, 18.4044034134, 6.83970779714]
            sage: (10^-15 < sv[3]) and (sv[3] < 10^-13)
            True
            sage: (10^-16 < sv[4]) and (sv[4] < 10^-14)
            True

        A full-rank matrix that is ill-conditioned.  We use this to
        illustrate ways of using the various possibilities for ``eps``,
        including one that is ill-advised. Notice that the automatically
        computed cutoff gets this (difficult) example slightly wrong.
        This illustrates the impossibility of any automated process always
        getting this right.  Use with caution and judgement.  ::

            sage: entries = [1/(i+j+1) for i in range(12) for j in range(12)]
            sage: B = matrix(QQ, 12, 12, entries)
            sage: B.rank()
            12
            sage: A = B.change_ring(RDF)
            sage: A.condition() > 1.7 * 10^16
            True

            sage: sv = A.singular_values(eps=None)
            sage: [round(sv[i],15) for i in range(12)]
            [1.79537205956, 0.380275245955, 0.0447385487522, 0.00372231223789,
             0.000233089089022, 1.1163357483e-05, 4.08237611e-07,
             1.1228611e-08, 2.25196e-10, 3.111e-12, 2.6e-14, 0.0]
            sage: (10^-17 < sv[11]) and (sv[11] < 10^-15)
            True

            sage: sv = A.singular_values(eps='auto')
            sage: [round(sv[i],15) for i in range(12)]
            [1.79537205956, 0.380275245955, 0.0447385487522, 0.00372231223789,
             0.000233089089022, 1.1163357483e-05, 4.08237611e-07,
             1.1228611e-08, 2.25196e-10, 3.111e-12, 2.6e-14, 0.0]
            sage: sv[11] == 0.0
            True

            sage: sv = A.singular_values(eps=1e-4)
            sage: [round(sv[i],15) for i in range(12)]
            [1.79537205956, 0.380275245955, 0.0447385487522, 0.00372231223789,
             0.000233089089022, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sage: all([sv[i] == 0.0 for i in range(5, 12)])
            True

        With Sage's "verbose" facility, you can compactly see the cutoff
        at work.  In any application of this routine, or those that build upon
        it, it would be a good idea to conduct this exercise on samples.
        We also test here that all the  values are returned in `RDF` since
        singular values are always real. ::

            sage: A = matrix(CDF, 4, range(16))
            sage: set_verbose(1)
            sage: sv = A.singular_values(eps='auto'); sv
            verbose 1 (<module>) singular values,
            smallest-non-zero:cutoff:largest-zero,
            2.2766...:6.2421...e-14:...
            [35.139963659, 2.27661020871, 0.0, 0.0]
            sage: set_verbose(0)

            sage: all([s in RDF for s in sv])
            True

        TESTS:

        Bogus values of the ``eps`` keyword will be caught.  ::

            sage: A.singular_values(eps='junk')
            Traceback (most recent call last):
            ...
            ValueError: could not convert string to float: junk

        REFERENCES:

        .. [WATKINS] Watkins, David S. Fundamentals of Matrix Computations,
           Third Edition.  Wiley, Hoboken, New Jersey, 2010.

        AUTHOR:

        - Rob Beezer - (2011-02-18)
        """
        from sage.misc.misc import verbose
        from sage.rings.real_double import RDF
        global scipy
        # get SVD decomposition, which is a cached quantity
        _, S, _ = self.SVD()
        diag = min(self._nrows, self._ncols)
        sv = [RDF(S[i,i]) for i in range(diag)]
        # no cutoff, send raw data back
        if eps is None:
            verbose("singular values, no zero cutoff specified", level=1)
            return sv
        # set cutoff as RDF element
        if eps == 'auto':
            if scipy is None: import scipy
            eps = 2*max(self._nrows, self._ncols)*scipy.finfo(float).eps*sv[0]
        eps = RDF(eps)
        # locate non-zero entries
        rank = 0
        while rank < diag and sv[rank] > eps:
            rank = rank + 1
        # capture info for watching zero cutoff behavior at verbose level 1
        if rank == 0:
            small_nonzero = None
        else:
            small_nonzero = sv[rank-1]
        if rank < diag:
            large_zero = sv[rank]
        else:
            large_zero = None
        # convert small values to zero, then done
        for i in range(rank, diag):
            sv[i] = RDF(0)
        verbose("singular values, smallest-non-zero:cutoff:largest-zero, %s:%s:%s" % (small_nonzero, eps, large_zero), level=1)
        return sv

    def LU(self):
        r"""
        Returns a decomposition of the (row-permuted) matrix as a product of
        a lower-triangular matrix ("L") and an upper-triangular matrix ("U").

        OUTPUT:

        For an `m\times n` matrix ``A`` this method returns a triple of
        immutable matrices ``P, L, U`` such that

        - ``P*A = L*U``
        - ``P`` is a square permutation matrix, of size `m\times m`,
          so is all zeroes, but with exactly a single one in each
          row and each column.
        - ``L`` is lower-triangular, square of size `m\times m`,
          with every diagonal entry equal to one.
        - ``U`` is upper-triangular with size `m\times n`, i.e.
          entries below the "diagonal" are all zero.

        The computed decomposition is cached and returned on
        subsequent calls, thus requiring the results to be immutable.

        Effectively, ``P`` permutes the rows of ``A``.  Then ``L``
        can be viewed as a sequence of row operations on this matrix,
        where each operation is adding a multiple of a row to a
        subsequent row.  There is no scaling (thus 1's on the diagonal
        of ``L``) and no row-swapping (``P`` does that).  As a result
        ``U`` is close to being the result of Gaussian-elimination.
        However, round-off errors can make it hard to determine
        the zero entries of ``U``.

        .. NOTE::

            Sometimes this decomposition is written as ``A=P*L*U``,
            where ``P`` represents the inverse permutation and is
            the matrix inverse of the ``P`` returned by this method.
            The computation of this matrix inverse can be accomplished
            quickly with just a transpose as the matrix is orthogonal/unitary.

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

        Trac 10839 made this routine available for rectangular matrices.  ::

            sage: A = matrix(RDF, 5, 6, range(30)); A
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 6.0  7.0  8.0  9.0 10.0 11.0]
            [12.0 13.0 14.0 15.0 16.0 17.0]
            [18.0 19.0 20.0 21.0 22.0 23.0]
            [24.0 25.0 26.0 27.0 28.0 29.0]
            sage: P, L, U = A.LU()
            sage: P
            [0.0 0.0 0.0 0.0 1.0]
            [1.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 1.0 0.0 0.0]
            [0.0 0.0 0.0 1.0 0.0]
            [0.0 1.0 0.0 0.0 0.0]
            sage: L.zero_at(0)   # Use zero_at(0) to get rid of signed zeros
            [ 1.0  0.0  0.0  0.0  0.0]
            [ 0.0  1.0  0.0  0.0  0.0]
            [ 0.5  0.5  1.0  0.0  0.0]
            [0.75 0.25  0.0  1.0  0.0]
            [0.25 0.75  0.0  0.0  1.0]
            sage: U.zero_at(0)   # Use zero_at(0) to get rid of signed zeros
            [24.0 25.0 26.0 27.0 28.0 29.0]
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            sage: P*A-L*U
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            sage: P.transpose()*L*U
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 6.0  7.0  8.0  9.0 10.0 11.0]
            [12.0 13.0 14.0 15.0 16.0 17.0]
            [18.0 19.0 20.0 21.0 22.0 23.0]
            [24.0 25.0 26.0 27.0 28.0 29.0]

        Trivial cases return matrices of the right size and
        characteristics.  ::

            sage: A = matrix(RDF, 5, 0, entries=0)
            sage: P, L, U = A.LU()
            sage: P.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Double Field
            sage: L.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Double Field
            sage: U.parent()
            Full MatrixSpace of 5 by 0 dense matrices over Real Double Field
            sage: P*A-L*U
            []

        The results are immutable since they are cached.  ::

            sage: P, L, U = matrix(RDF, 2, 2, range(4)).LU()
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
        m = self._nrows
        n = self._ncols

        # scipy fails on trivial cases
        if m == 0 or n == 0:
            P = self._new(m, m)
            for i in range(m):
                P[i,i]=1
            P.set_immutable()
            L = P
            U = self._new(m,n)
            U.set_immutable()
            return P, L, U

        PLU = self.fetch('PLU_factors')
        if not PLU is None:
            return PLU
        if scipy is None:
            import scipy
        import scipy.linalg
        if numpy is None:
            import numpy
        PM, LM, UM = scipy.linalg.lu(self._matrix_numpy)
        # Numpy has a different convention than we had with GSL
        # So we invert (transpose) the P to match our prior behavior
        # TODO: It's an awful waste to store a huge matrix for P, which
        # is just a simple permutation, really.
        P = self._new(m, m)
        L = self._new(m, m)
        U = self._new(m, n)
        P._matrix_numpy = PM.T.copy()
        L._matrix_numpy = numpy.ascontiguousarray(LM)
        U._matrix_numpy = numpy.ascontiguousarray(UM)
        PLU = (P, L, U)
        for M in PLU:
            M.set_immutable()
        self.cache('PLU_factors', PLU)
        return PLU

    def eigenspaces_left(self, var='a', algebraic_multiplicity=False):
        r"""
        Computes the left eigenspaces of a matrix of double precision
        real or complex numbers (i.e. RDF or CDF).

        .. warning::

            This method returns eigenspaces that are all of
            dimension one, since it is impossible to ascertain
            if the numerical results belong to the same eigenspace.
            So this is deprecated in favor of the eigenmatrix routines,
            such as :meth:`sage.matrix.matrix2.Matrix.eigenmatrix_right`.

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

        EXAMPLES:

        This first test simply raises the deprecation warning.  ::

            sage: A = identity_matrix(RDF, 2)
            sage: es = A.eigenspaces_left()
            doctest:...: DeprecationWarning: Eigenspaces of RDF/CDF matrices are
            deprecated as of Sage version 5.0,
            please use "eigenmatrix_left" instead

        ::

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
        from sage.misc.misc import deprecation
        msg = ('Eigenspaces of RDF/CDF matrices are deprecated as of ',
               'Sage version 5.0',
               ', please use "eigenmatrix_left" instead')
        deprecation(''.join(msg))
        # For numerical values we leave decisions about
        # multiplicity to the calling routine
        if algebraic_multiplicity:
            raise ValueError("algebraic_multiplicity must be set to False for double precision matrices")
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

        .. warning::

            This method returns eigenspaces that are all of
            dimension one, since it is impossible to ascertain
            if the numerical results belong to the same eigenspace.
            So this is deprecated in favor of the eigenmatrix routines,
            such as :meth:`sage.matrix.matrix2.Matrix.eigenmatrix_right`.

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


        EXAMPLES:

        This first test simply raises the deprecation warning.  ::

            sage: A = identity_matrix(RDF, 2)
            sage: es = A.eigenspaces_right()
            doctest:...: DeprecationWarning: Eigenspaces of RDF/CDF matrices are
            deprecated as of Sage version 5.0,
            please use "eigenmatrix_right" instead

        ::

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
            sage: abs(abs(diff)) < 3e-14
            True

        TESTS::

            sage: m.eigenspaces_right(algebraic_multiplicity=True)
            Traceback (most recent call last):
            ...
            ValueError: algebraic_multiplicity must be set to False for double precision matrices
        """
        from sage.misc.misc import deprecation
        msg = ('Eigenspaces of RDF/CDF matrices are deprecated as of ',
               'Sage version 5.0',
               ', please use "eigenmatrix_right" instead')
        deprecation(''.join(msg))
        # For numerical values we leave decisions about
        # multiplicity to the calling routine
        if algebraic_multiplicity:
            raise ValueError("algebraic_multiplicity must be set to False for double precision matrices")
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

        EXAMPLES::

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
            raise ArithmeticError("self must be a square matrix")
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
            raise ArithmeticError("self must be a square matrix")
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
            raise ArithmeticError("self must be a square matrix")
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
        Solve the equation `A x = b` using LU decomposition.

        .. WARNING::

            This function is broken. See trac 4932.

        INPUT:

        - self -- an invertible matrix
        - b -- a vector

        .. NOTE::

            This method precomputes and stores the LU decomposition
            before solving. If many equations of the form Ax=b need to be
            solved for a singe matrix A, then this method should be used
            instead of solve. The first time this method is called it will
            compute the LU decomposition.  If the matrix has not changed
            then subsequent calls will be very fast as the precomputed LU
            decomposition will be reused.

        EXAMPLES::

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

        We test two degenerate cases::

            sage: A = matrix(RDF, 0, 3, [])
            sage: A.solve_left_LU(vector(RDF,[]))
            (0.0, 0.0, 0.0)
            sage: A = matrix(RDF, 3, 0, [])
            sage: A.solve_left_LU(vector(RDF,3, [1,2,3]))
            ()

        """
        if self._nrows != b.degree():
            raise ValueError("number of rows of self must equal degree of b")
        if self._nrows == 0 or self._ncols == 0:
            return self._row_ambient_module().zero_vector()

        raise NotImplementedError("this function is not finished (see trac 4932)")
        self._c_compute_LU()  # so self._L_M and self._U_M are defined below.
        cdef Matrix_double_dense M = self._new()
        lu = self._L_M*self._U_M
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        M._matrix_numpy = scipy.linalg.lu_solve((lu, self._P_M), b)
        return M

    def solve_right(self, b):
        r"""
        Solve the vector equation ``A*x = b`` for a nonsingular ``A``.

        INPUT:

        - ``self`` - a square matrix that is nonsigular (of full rank).
        - ``b`` - a vector of the correct size.  Elements of the vector
          must coerce into the base ring of the coefficient matrix.  In
          particular, if ``b`` has entries from ``CDF`` then ``self``
          must have ``CDF`` as its base ring.

        OUTPUT:

        The unique solution ``x`` to the matrix equation ``A*x = b``,
        as a vector over the same base ring as ``self``.

        ALGORITHM:

        Uses the ``solve()`` routine from the SciPy ``scipy.linalg`` module.

        EXAMPLES:

        Over the reals. ::

            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF,[1,2,3])
            sage: x = A.solve_right(b); x
            (-0.113695090439, 1.39018087855, -0.333333333333)
            sage: x.parent()
            Vector space of dimension 3 over Real Double Field
            sage: A*x
            (1.0, 2.0, 3.0)

        Over the complex numbers.  ::

            sage: A = matrix(CDF, [[      0, -1 + 2*I,  1 - 3*I,        I],
            ...                    [2 + 4*I, -2 + 3*I, -1 + 2*I,   -1 - I],
            ...                    [  2 + I,    1 - I,       -1,        5],
            ...                    [    3*I,   -1 - I,   -1 + I,   -3 + I]])
            sage: b = vector(CDF, [2 -3*I, 3, -2 + 3*I, 8])
            sage: x = A.solve_right(b); x
            (1.96841637... - 1.07606761...*I, -0.614323843... + 1.68416370...*I, 0.0733985765... + 1.73487544...*I, -1.6018683... + 0.524021352...*I)
            sage: x.parent()
            Vector space of dimension 4 over Complex Double Field
            sage: abs(A*x - b) < 1e-14
            True

        The vector of constants, ``b``, can be given in a
        variety of forms, so long as it coerces to a vector
        over the same base ring as the coefficient matrix.  ::

            sage: A=matrix(CDF, 5, [1/(i+j+1) for i in range(5) for j in range(5)])
            sage: A.solve_right([1]*5)
            (5.0..., -120.0..., 630.0..., -1120.0..., 630.0...)

        TESTS:

        A degenerate case. ::

            sage: A = matrix(RDF, 0, 0, [])
            sage: A.solve_right(vector(RDF,[]))
            ()

        The coefficent matrix must be square. ::

            sage: A = matrix(RDF, 2, 3, range(6))
            sage: b = vector(RDF, [1,2,3])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            ValueError: coefficient matrix of a system over RDF/CDF must be square, not 2 x 3

        The coefficient matrix must be nonsingular.  ::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4,5])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            LinAlgError: singular matrix

        The vector of constants needs the correct degree.  ::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Real Double Field incompatible with matrix over Real Double Field

        The vector of constants needs to be compatible with
        the base ring of the coefficient matrix.  ::

            sage: F.<a> = FiniteField(27)
            sage: b = vector(F, [a,a,a,a,a])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Finite Field in a of size 3^3 incompatible with matrix over Real Double Field

        With a coefficient matrix over ``RDF``, a vector of constants
        over ``CDF`` can be accomodated by converting the base ring
        of the coefficient matrix.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: b = vector(CDF, [1+I,2])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Complex Double Field incompatible with matrix over Real Double Field

            sage: B = A.change_ring(CDF)
            sage: B.solve_right(b)
            (-0.5 - 1.5*I, 1.0 + 1.0*I)
        """
        if not self.is_square():
            raise ValueError("coefficient matrix of a system over RDF/CDF must be square, not %s x %s " % (self.nrows(), self.ncols()))
        M = self._column_ambient_module()
        try:
            vec = M(b)
        except TypeError:
            raise TypeError("vector of constants over %s incompatible with matrix over %s" % (b.base_ring(), self.base_ring()))
        if vec.degree() != self.ncols():
            raise ValueError("vector of constants in linear system over RDF/CDF must have degree equal to the number of columns for the coefficient matrix, not %s" % vec.degree() )

        if self._ncols == 0:
            return M.zero_vector()

        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        # may raise a LinAlgError for a singular matrix
        return M(scipy.linalg.solve(self._matrix_numpy, vec.numpy()))

    def solve_left(self, b):
        r"""
        Solve the vector equation ``x*A = b`` for a nonsingular ``A``.

        INPUT:

        - ``self`` - a square matrix that is nonsigular (of full rank).
        - ``b`` - a vector of the correct size.  Elements of the vector
          must coerce into the base ring of the coefficient matrix.  In
          particular, if ``b`` has entries from ``CDF`` then ``self``
          must have ``CDF`` as its base ring.

        OUTPUT:

        The unique solution ``x`` to the matrix equation ``x*A = b``,
        as a vector over the same base ring as ``self``.

        ALGORITHM:

        Uses the ``solve()`` routine from the SciPy ``scipy.linalg`` module,
        after taking the tranpose of the coefficient matrix.

        EXAMPLES:

        Over the reals. ::

            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF,[1,2,3])
            sage: x = A.solve_left(b); x.zero_at(1e-18) # fix noisy zeroes
            (0.666666666..., 0.0, 0.333333333...)
            sage: x.parent()
            Vector space of dimension 3 over Real Double Field
            sage: x*A
            (1.0, 2.0, 3.0)

        Over the complex numbers.  ::

            sage: A = matrix(CDF, [[      0, -1 + 2*I,  1 - 3*I,        I],
            ...                    [2 + 4*I, -2 + 3*I, -1 + 2*I,   -1 - I],
            ...                    [  2 + I,    1 - I,       -1,        5],
            ...                    [    3*I,   -1 - I,   -1 + I,   -3 + I]])
            sage: b = vector(CDF, [2 -3*I, 3, -2 + 3*I, 8])
            sage: x = A.solve_left(b); x
            (-1.55765124... - 0.644483985...*I, 0.183274021... + 0.286476868...*I, 0.270818505... + 0.246619217...*I, -1.69003558... - 0.828113879...*I)
            sage: x.parent()
            Vector space of dimension 4 over Complex Double Field
            sage: abs(x*A - b) < 1e-14
            True

        The vector of constants, ``b``, can be given in a
        variety of forms, so long as it coerces to a vector
        over the same base ring as the coefficient matrix.  ::

            sage: A=matrix(CDF, 5, [1/(i+j+1) for i in range(5) for j in range(5)])
            sage: A.solve_left([1]*5)
            (5.0..., -120.0..., 630.0..., -1120.0..., 630.0...)

        TESTS:

        A degenerate case. ::

            sage: A = matrix(RDF, 0, 0, [])
            sage: A.solve_left(vector(RDF,[]))
            ()

        The coefficent matrix must be square. ::

            sage: A = matrix(RDF, 2, 3, range(6))
            sage: b = vector(RDF, [1,2,3])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            ValueError: coefficient matrix of a system over RDF/CDF must be square, not 2 x 3

        The coefficient matrix must be nonsingular.  ::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4,5])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            LinAlgError: singular matrix

        The vector of constants needs the correct degree.  ::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Real Double Field incompatible with matrix over Real Double Field

        The vector of constants needs to be compatible with
        the base ring of the coefficient matrix.  ::

            sage: F.<a> = FiniteField(27)
            sage: b = vector(F, [a,a,a,a,a])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Finite Field in a of size 3^3 incompatible with matrix over Real Double Field

        With a coefficient matrix over ``RDF``, a vector of constants
        over ``CDF`` can be accomodated by converting the base ring
        of the coefficient matrix.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: b = vector(CDF, [1+I,2])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            TypeError: vector of constants over Complex Double Field incompatible with matrix over Real Double Field

            sage: B = A.change_ring(CDF)
            sage: B.solve_left(b)
            (0.5 - 1.5*I, 0.5 + 0.5*I)
        """
        if not self.is_square():
            raise ValueError("coefficient matrix of a system over RDF/CDF must be square, not %s x %s " % (self.nrows(), self.ncols()))
        M = self._row_ambient_module()
        try:
            vec = M(b)
        except TypeError:
            raise TypeError("vector of constants over %s incompatible with matrix over %s" % (b.base_ring(), self.base_ring()))
        if vec.degree() != self.nrows():
            raise ValueError("vector of constants in linear system over RDF/CDF must have degree equal to the number of rows for the coefficient matrix, not %s" % vec.degree() )

        if self._nrows == 0:
            return M.zero_vector()

        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        # may raise a LinAlgError for a singular matrix
        # call "right solve" routine with the transpose
        return M(scipy.linalg.solve(self._matrix_numpy.T, vec.numpy()))

    def determinant(self):
        """
        Return the determinant of self.

        ALGORITHM:

        Use numpy

        EXAMPLES::

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
            raise ValueError("self must be a square matrix")
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

        .. NOTE::

            This is useful if the usual determinant overflows.

        EXAMPLES::

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
            raise ArithmeticError("self must be a square matrix")

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

        ``.T`` is a convenient shortcut for the transpose::

            sage: m.T
            [2.0 3.0]
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
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
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

        - A -- a matrix

        OUTPUT:

        - U, S, V -- immutable matrices such that $A = U*S*V.conj().transpose()$
          where U and V are orthogonal and S is zero off of the diagonal.

        Note that if self is m-by-n, then the dimensions of the
        matrices that this returns are (m,m), (m,n), and (n, n).

        .. NOTE::

            If all you need is the singular values of the matrix, see
            the more convenient :meth:`singular_values`.

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

        - self -- a real matrix A

        OUTPUT:

        - Q, R -- immutable matrices such that A = Q*R such that the columns of Q are
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

        EXAMPLES::

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

    def is_unitary(self, tol=1e-12):
        r"""
        Returns ``True`` if the columns of the matrix are an orthonormal basis.

        For a matrix with real entries this determines if a matrix is
        "orthogonal" and for a matrix with complex entries this determines
        if the matrix is "unitary."

        INPUT:

        - ``tol`` - default: ``1e-12`` - the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        OUTPUT:

        ``True`` if the matrix is square and its conjugate-transpose is
        its inverse, and ``False`` otherwise.  In other words, a matrix
        is orthogonal or unitary if the product of its conjugate-transpose
        times the matrix is the identity matrix.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance basis.

        EXAMPLES:

        A matrix that is far from unitary. ::

            sage: A = matrix(RDF, 4, range(16))
            sage: A.conjugate().transpose()*A
            [224.0 248.0 272.0 296.0]
            [248.0 276.0 304.0 332.0]
            [272.0 304.0 336.0 368.0]
            [296.0 332.0 368.0 404.0]
            sage: A.is_unitary()
            False

        The QR decoposition will produce a unitary matrix as Q and the
        SVD decomposition will create two unitary matrices, U and V. ::

            sage: A = matrix(CDF, [[   1 - I,   -3*I,  -2 + I,        1, -2 + 3*I],
            ...                    [   1 - I, -2 + I, 1 + 4*I,        0,    2 + I],
            ...                    [      -1, -5 + I,  -2 + I,    1 + I, -5 - 4*I],
            ...                    [-2 + 4*I,  2 - I, 8 - 4*I,  1 - 8*I,  3 - 2*I]])
            sage: Q, R = A.QR()
            sage: Q.is_unitary()
            True
            sage: U, S, V = A.SVD()
            sage: U.is_unitary()
            True
            sage: V.is_unitary()  # not tested - known bug (trac #11248)
            True

        If we make the tolerance too strict we can get misleading results.  ::

            sage: A = matrix(RDF, 10, 10, [1/(i+j+1) for i in range(10) for j in range(10)])
            sage: Q, R = A.QR()
            sage: Q.is_unitary(tol=1e-16)
            False

        Rectangular matrices are not unitary/orthogonal, even if their
        columns form an orthonormal set.  ::

            sage: A = matrix(CDF, [[1,0], [0,0], [0,1]])
            sage: A.is_unitary()
            False

        A trivial case.  ::

            sage: P = matrix(CDF,0,0)
            sage: P.is_unitary()
            True
        """
        global numpy
        tol = float(tol)
        key = 'unitary_%s'%tol
        b = self.fetch(key)
        if not b is None:
            return b
        if not self.is_square():
            self.cache(key, False)
            return False
        if numpy is None:
            import numpy
        cdef Matrix_double_dense P
        P = self.conjugate().transpose()*self
        cdef Py_ssize_t i, j
        unitary = True
        for i from 0 < i < self._nrows:
            # off-diagonal
            for j from 0 <= j < i:
                if numpy.absolute(P.get_unsafe(i,j)) > tol:
                    unitary = False
                    break
            # at diagonal
            if numpy.absolute(P.get_unsafe(i,i)-1) > tol:
                unitary = False
                break
        self.cache(key, unitary)
        return unitary

    def _is_lower_triangular(self, tol):
        r"""
        Returns ``True`` if the entries above the diagonal are all zero.

        INPUT:

        - ``tol`` -  the largest value of the absolute value of the
          difference between two matrix entries for which they will
          still be considered equal.

        OUTPUT:

        Returns ``True`` if each entry above the diagonal (entries
        with a row index less than the column index) is zero.

        EXAMPLES::

            sage: A = matrix(RDF, [[ 2.0, 0.0,  0.0],
            ...                    [ 1.0, 3.0,  0.0],
            ...                    [-4.0, 2.0, -1.0]])
            sage: A._is_lower_triangular(1.0e-17)
            True
            sage: A[1,2] = 10^-13
            sage: A._is_lower_triangular(1.0e-14)
            False
        """
        global numpy
        if numpy is None:
            import numpy
        cdef Py_ssize_t i, j
        for i in range(self._nrows):
            for j in range(i+1, self._ncols):
                if numpy.absolute(self.get_unsafe(i,j)) > tol:
                    return False
        return True

    def is_hermitian(self, tol = 1e-12, algorithm='orthonormal'):
        r"""
        Returns ``True`` if the matrix is equal to its conjugate-transpose.

        INPUT:

        - ``tol`` - default: ``1e-12`` - the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        - ``algorithm`` - default: 'orthonormal' - set to 'orthonormal'
          for a stable procedure and set to 'naive' for a fast
          procedure.

        OUTPUT:

        ``True`` if the matrix is square and equal to the transpose
        with every entry conjugated, and ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_symmetric`
        method is equivalent and faster.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply compares corresponing entries on either
        side of the diagonal (and on the diagonal itself) to see if they are
        conjugates, with equality controlled by the tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix with real entries.

        So the naive algorithm can finish quickly for a matrix that is not
        Hermitian, while the orthonormal algorithm will always compute a
        Schur decomposition before going through a similar check of the matrix
        entry-by-entry.

        EXAMPLES::

            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ...                    [-3 - I,     -4*I,     -2],
            ...                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: A.is_hermitian(algorithm='orthonormal')
            False
            sage: A.is_hermitian(algorithm='naive')
            False
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian(algorithm='orthonormal')
            True
            sage: B.is_hermitian(algorithm='naive')
            True

        A matrix that is nearly Hermitian, but for one non-real
        diagonal entry. ::

            sage: A = matrix(CDF, [[    2,   2-I, 1+4*I],
            ...                    [  2+I,   3+I, 2-6*I],
            ...                    [1-4*I, 2+6*I,     5]])
            sage: A.is_hermitian(algorithm='orthonormal')
            False
            sage: A[1,1] = 132
            sage: A.is_hermitian(algorithm='orthonormal')
            True

        We get a unitary matrix from the SVD routine and use this
        numerical matrix to create a matrix that should be Hermitian
        (indeed it should be the identity matrix), but with some
        imprecision.  We use this to illustrate that if the tolerance
        is set too small, then we can be too strict about the equality
        of entries and may achieve the wrong result (depending on
        the system)::

            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ...                    [-3 - I,     -4*I,     -2],
            ...                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: U, _, _ = A.SVD()
            sage: B=U*U.conjugate_transpose()
            sage: B.is_hermitian(algorithm='naive')
            True
            sage: B.is_hermitian(algorithm='naive', tol=1.0e-17)  # random
            False
            sage: B.is_hermitian(algorithm='naive', tol=1.0e-15)
            True

        A square, empty matrix is trivially Hermitian.  ::

            sage: A = matrix(RDF, 0, 0)
            sage: A.is_hermitian()
            True

        Rectangular matrices are never Hermitian, no matter which
        algorithm is requested.  ::

            sage: A = matrix(CDF, 3, 4)
            sage: A.is_hermitian()
            False

        TESTS:

        The tolerance must be strictly positive.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_hermitian(tol = -3.1)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive, not -3.1

        The ``algorithm`` keyword gets checked.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_hermitian(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk

        AUTHOR:

        - Rob Beezer (2011-03-30)
        """
        import sage.rings.complex_double
        global numpy
        tol = float(tol)
        if tol <= 0:
            raise ValueError('tolerance must be positive, not {0}'.format(tol))
        if not algorithm in ['naive', 'orthonormal']:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))

        key = 'hermitian_{0}_{1}'.format(algorithm, tol)
        h = self.fetch(key)
        if not h is None:
            return h
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True
        if numpy is None:
            import numpy
        cdef Py_ssize_t i, j
        cdef Matrix_double_dense T
        if algorithm == 'orthonormal':
            # Schur decomposition over CDF will be diagonal and real iff Hermitian
            _, T = self.schur(base_ring=sage.rings.complex_double.CDF)
            hermitian = T._is_lower_triangular(tol)
            if hermitian:
                for i in range(T._nrows):
                    if numpy.absolute(numpy.imag(T.get_unsafe(i,i))) > tol:
                        hermitian = False
                        break
        elif algorithm == 'naive':
            hermitian = True
            for i in range(self._nrows):
                for j in range(i+1):
                    if numpy.absolute(self.get_unsafe(i,j) - self.get_unsafe(j,i).conjugate()) > tol:
                        hermitian = False
                        break
                if not hermitian:
                    break
        self.cache(key, hermitian)
        return hermitian

    def is_normal(self, tol=1e-12, algorithm='orthonormal'):
        r"""
        Returns ``True`` if the matrix commutes with its conjugate-transpose.

        INPUT:

        - ``tol`` - default: ``1e-12`` - the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        - ``algorithm`` - default: 'orthonormal' - set to 'orthonormal'
          for a stable procedure and set to 'naive' for a fast
          procedure.

        OUTPUT:

        ``True`` if the matrix is square and commutes with its
        conjugate-transpose, and ``False`` otherwise.

        Normal matrices are precisely those that can be diagonalized
        by a unitary matrix.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply compares entries of the two possible
        products of the matrix with its conjugate-transpose, with equality
        controlled by the tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix.  An orthonormal diagonalization
        is equivalent to being normal.

        So the naive algorithm can finish fairly quickly for a matrix
        that is not normal, once the products have been computed.
        However, the orthonormal algorithm will compute a Schur
        decomposition before going through a similar check of a
        matrix entry-by-entry.

        EXAMPLES:

        First over the complexes.  ``B`` is Hermitian, hence normal.  ::

            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ...                    [-3 - I,     -4*I,     -2],
            ...                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian()
            True
            sage: B.is_normal(algorithm='orthonormal')
            True
            sage: B.is_normal(algorithm='naive')
            True
            sage: B[0,0] = I
            sage: B.is_normal(algorithm='orthonormal')
            False
            sage: B.is_normal(algorithm='naive')
            False

        Now over the reals.  Circulant matrices are normal.  ::

            sage: G = graphs.CirculantGraph(20, [3, 7])
            sage: D = digraphs.Circuit(20)
            sage: A = 3*D.adjacency_matrix() - 5*G.adjacency_matrix()
            sage: A = A.change_ring(RDF)
            sage: A.is_normal()
            True
            sage: A.is_normal(algorithm = 'naive')
            True
            sage: A[19,0] = 4.0
            sage: A.is_normal()
            False
            sage: A.is_normal(algorithm = 'naive')
            False

        Skew-Hermitian matrices are normal.  ::

            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ...                    [-3 - I,     -4*I,     -2],
            ...                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: B = A - A.conjugate_transpose()
            sage: B.is_hermitian()
            False
            sage: B.is_normal()
            True
            sage: B.is_normal(algorithm='naive')
            True

        A small matrix that does not fit into any of the usual categories
        of normal matrices.  ::

            sage: A = matrix(RDF, [[1, -1],
            ...                    [1,  1]])
            sage: A.is_normal()
            True
            sage: not A.is_hermitian() and not A.is_skew_symmetric()
            True

        Sage has several fields besides the entire complex numbers
        where conjugation is non-trivial. ::

            sage: F.<b> = QuadraticField(-7)
            sage: C = matrix(F, [[-2*b - 3,  7*b - 6, -b + 3],
            ...                  [-2*b - 3, -3*b + 2,   -2*b],
            ...                  [   b + 1,        0,     -2]])
            sage: C = C*C.conjugate_transpose()
            sage: C.is_normal()
            True

        We get a unitary matrix from the SVD routine and use this
        numerical matrix to create a matrix that should be normal
        (indeed it should be the identity matrix), but with some
        imprecision.  We use this to illustrate that if the tolerance
        is set too small, then we can be too strict about the equality
        of entries and achieve the wrong result.  ::

            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ...                    [-3 - I,     -4*I,     -2],
            ...                    [-1 + I, -2 - 8*I,  2 + I]])
            sage: U, _, _ = A.SVD()
            sage: B=U*U.conjugate_transpose()
            sage: B.is_normal(algorithm='naive')
            True
            sage: B.is_normal(algorithm='naive', tol=1.0e-34)
            False
            sage: B.is_normal(algorithm='naive', tol=1.0e-31)
            True

        A square, empty matrix is trivially normal.  ::

            sage: A = matrix(CDF, 0, 0)
            sage: A.is_normal()
            True

        Rectangular matrices are never normal, no matter which
        algorithm is requested.  ::

            sage: A = matrix(CDF, 3, 4)
            sage: A.is_normal()
            False

        TESTS:

        The tolerance must be strictly positive.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_normal(tol = -3.1)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive, not -3.1

        The ``algorithm`` keyword gets checked.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_normal(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk

        AUTHOR:

         - Rob Beezer (2011-03-31)
        """
        import sage.rings.complex_double
        global numpy
        tol = float(tol)
        if tol <= 0:
            raise ValueError('tolerance must be positive, not {0}'.format(tol))
        if not algorithm in ['naive', 'orthonormal']:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))

        key = 'normal_{0}_{1}'.format(algorithm, tol)
        b = self.fetch(key)
        if not b is None:
            return b
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True
        cdef Py_ssize_t i, j
        cdef Matrix_double_dense T, left, right
        if algorithm == 'orthonormal':
            # Schur decomposition over CDF will be diagonal iff normal
            _, T = self.schur(base_ring=sage.rings.complex_double.CDF)
            normal = T._is_lower_triangular(tol)
        elif algorithm == 'naive':
            if numpy is None:
                import numpy
            CT = self.conjugate_transpose()
            left = self*CT
            right = CT*self
            normal = True
            # two products are Hermitian, need only check lower triangle
            for i in range(self._nrows):
                for j in range(i+1):
                    if numpy.absolute(left.get_unsafe(i,j) - right.get_unsafe(i,j)) > tol:
                        normal = False
                        break
                if not normal:
                    break
        self.cache(key, normal)
        return normal

    def schur(self, base_ring=None):
        r"""
        Returns the Schur decomposition of the matrix.

        INPUT:

        - ``base_ring`` - optional, defaults to the base ring of ``self``.
          Use this to request the base ring of the returned matrices, which
          will affect the format of the results.

        OUTPUT:

        A pair of immutable matrices.  The first is a unitary matrix `Q`.
        The second, `T`, is upper-triangular when returned over the complex
        numbers, while it is almost upper-triangular over the reals.  In the
        latter case, there can be some `2\times 2` blocks on the diagonal
        which represent a pair of conjugate complex eigenvalues of ``self``.

        If ``self`` is the matrix `A`, then

        .. math::

            A = QT({\overline Q})^t

        where the latter matrix is the conjugate-transpose of ``Q``, which
        is also the inverse of ``Q``, since ``Q`` is unitary.

        Note that in the case of a normal matrix (Hermitian, symmetric, and
        others), the upper-triangular matrix is  a diagonal matrix with
        eigenvalues of ``self`` on the diagonal, and the unitary matrix
        has columns that form an orthonormal basis composed of eigenvectors
        of ``self``.  This is known as "orthonormal diagonalization".

        .. WARNING::

            The Schur decomposition is not unique, as there may be numerous
            choices for the vectors of the orthonormal basis, and consequently
            different possibilities for the upper-triangular matrix.  However,
            the diagonal of the upper-triangular matrix will always contain the
            eigenvalues of the matrix (in the complex version), or `2\times 2`
            block matrices in the real version representing pairs of conjugate
            complex eigenvalues.

            In particular, results may vary across systems and processors.

        EXAMPLES:

        First over the complexes.  The similar matrix is always
        upper-triangular in this case.  ::

            sage: A = matrix(CDF, 4, 4, range(16)) + matrix(CDF, 4, 4, [x^3*I for x in range(0, 16)])
            sage: Q, T = A.schur()
            sage: (Q*Q.conjugate().transpose()).zero_at(1.0e-12)
            [1.0   0   0   0]
            [  0 1.0   0   0]
            [  0   0 1.0   0]
            [  0   0   0 1.0]
            sage: all([T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i)])
            True
            sage: (Q*T*Q.conjugate().transpose()-A).zero_at(1.0e-11)
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: eigenvalues = [T[i,i] for i in range(4)]; eigenvalues
            [30.733... + 4648.541...*I, -0.184... - 159.057...*I, -0.523... + 11.158...*I, -0.025... - 0.642...*I]
            sage: A.eigenvalues()
            [30.733... + 4648.541...*I, -0.184... - 159.057...*I, -0.523... + 11.158...*I, -0.025... - 0.642...*I]
            sage: abs(A.norm()-T.norm()) < 1e-10
            True

        We begin with a real matrix but ask for a decomposition over the
        complexes.  The result will yield an upper-triangular matrix over
        the complex numbers for ``T``. ::

            sage: A = matrix(RDF, 4, 4, [x^3 for x in range(16)])
            sage: Q, T = A.schur(base_ring=CDF)
            sage: (Q*Q.conjugate().transpose()).zero_at(1.0e-12)
            [1.0   0   0   0]
            [  0 1.0   0   0]
            [  0   0 1.0   0]
            [  0   0   0 1.0]
            sage: T.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Complex Double Field
            sage: all([T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i)])
            True
            sage: (Q*T*Q.conjugate().transpose()-A).zero_at(1.0e-11)
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        Now totally over the reals.  But with complex eigenvalues, the
        similar matrix may not be upper-triangular. But "at worst" there
        may be some `2\times 2` blocks on the diagonal which represent
        a pair of conjugate complex eigenvalues. These blocks will then
        just interrupt the zeros below the main diagonal.  This example
        has a pair of these of the blocks. ::

            sage: A = matrix(RDF, 4, 4, [[1, 0, -3, -1],
            ...                          [4, -16, -7, 0],
            ...                          [1, 21, 1, -2],
            ...                          [26, -1, -2, 1]])
            sage: Q, T = A.schur()
            sage: (Q*Q.conjugate().transpose()).zero_at(1.0e-12)
            [1.0 0.0 0.0 0.0]
            [0.0 1.0 0.0 0.0]
            [0.0 0.0 1.0 0.0]
            [0.0 0.0 0.0 1.0]
            sage: all([T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i)])
            False
            sage: all([T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i-1)])
            True
            sage: (Q*T*Q.conjugate().transpose()-A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            sage: sorted(T[0:2,0:2].eigenvalues() + T[2:4,2:4].eigenvalues())
            [-5.710... - 8.382...*I, -5.710... + 8.382...*I, -0.789... - 2.336...*I, -0.789... + 2.336...*I]
            sage: sorted(A.eigenvalues())
            [-5.710... - 8.382...*I, -5.710... + 8.382...*I, -0.789... - 2.336...*I, -0.789... + 2.336...*I]
            sage: abs(A.norm()-T.norm()) < 1e-12
            True

        Starting with complex numbers and requesting a result over the reals
        will never happen.  ::

            sage: A = matrix(CDF, 2, 2, [[2+I, -1+3*I], [5-4*I, 2-7*I]])
            sage: A.schur(base_ring=RDF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert input matrix over CDF to a matrix over RDF

        If theory predicts your matrix is real, but it contains some
        very small imaginary parts, you can specify the cutoff for "small"
        imaginary parts, then request the output as real matrices, and let
        the routine do the rest. ::

            sage: A = matrix(RDF, 2, 2, [1, 1, -1, 0]) + matrix(CDF, 2, 2, [1.0e-14*I]*4)
            sage: B = A.zero_at(1.0e-12)
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
            sage: Q, T = B.schur(RDF)
            sage: Q.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Double Field
            sage: T.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Double Field
            sage: Q.round(6)
            [ 0.707107  0.707107]
            [-0.707107  0.707107]
            sage: T.round(6)
            [ 0.5  1.5]
            [-0.5  0.5]
            sage: (Q*T*Q.conjugate().transpose()-B).zero_at(1.0e-11)
            [0 0]
            [0 0]

        A Hermitian matrix has real eigenvalues, so the similar matrix
        will be upper-triangular.  Furthermore, a Hermitian matrix is
        diagonalizable with respect to an orthonormal basis, composed
        of eigenvectors of the matrix.  Here that basis is the set of
        columns of the unitary matrix.  ::

            sage: A = matrix(CDF, [[        52,   -9*I - 8,    6*I - 187,  -188*I + 2],
            ...                    [   9*I - 8,         12,   -58*I + 59,   30*I + 42],
            ...                    [-6*I - 187,  58*I + 59,         2677, 2264*I + 65],
            ...                    [ 188*I + 2, -30*I + 42, -2264*I + 65,       2080]])
            sage: Q, T = A.schur()
            sage: T = T.zero_at(1.0e-12).change_ring(RDF)
            sage: T.round(6)
            [4680.13301        0.0        0.0        0.0]
            [       0.0 102.715967        0.0        0.0]
            [       0.0        0.0  35.039344        0.0]
            [       0.0        0.0        0.0    3.11168]
            sage: (Q*Q.conjugate().transpose()).zero_at(1.0e-12)
            [1.0   0   0   0]
            [  0 1.0   0   0]
            [  0   0 1.0   0]
            [  0   0   0 1.0]
            sage: (Q*T*Q.conjugate().transpose()-A).zero_at(1.0e-11)
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        Similarly, a real symmetric matrix has only real eigenvalues,
        and there is an orthonormal basis composed of eigenvectors of
        the matrix.  ::

            sage: A = matrix(RDF, [[ 1, -2, 5, -3],
            ...                    [-2,  9, 1,  5],
            ...                    [ 5,  1, 3 , 7],
            ...                    [-3,  5, 7, -8]])
            sage: Q, T = A.schur()
            sage: Q.round(4)
            [-0.3027  -0.751   0.576 -0.1121]
            [  0.139 -0.3892 -0.2648  0.8713]
            [ 0.4361   0.359  0.7599  0.3217]
            [ -0.836  0.3945  0.1438  0.3533]
            sage: T = T.zero_at(10^-12)
            sage: all(abs(e) < 10^-4 for e in (T - diagonal_matrix(RDF, [-13.5698, -0.8508, 7.7664, 11.6542])).list())
            True
            sage: (Q*Q.transpose()).zero_at(1.0e-12)
            [1.0 0.0 0.0 0.0]
            [0.0 1.0 0.0 0.0]
            [0.0 0.0 1.0 0.0]
            [0.0 0.0 0.0 1.0]
            sage: (Q*T*Q.transpose()-A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]

        The results are cached, both as a real factorization and also as a
        complex factorization.  This means the returned matrices are
        immutable.  ::

            sage: A = matrix(RDF, 2, 2, [[0, -1], [1, 0]])
            sage: Qr, Tr = A.schur(base_ring=RDF)
            sage: Qc, Tc = A.schur(base_ring=CDF)
            sage: all([M.is_immutable() for M in [Qr, Tr, Qc, Tc]])
            True
            sage: Tr.round(6) != Tc.round(6)
            True

        TESTS:

        The Schur factorization is only defined for square matrices.  ::

            sage: A = matrix(RDF, 4, 5, range(20))
            sage: A.schur()
            Traceback (most recent call last):
            ...
            ValueError: Schur decomposition requires a square matrix, not a 4 x 5 matrix

        A base ring request is checked. ::

            sage: A = matrix(RDF, 3, range(9))
            sage: A.schur(base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: base ring of Schur decomposition matrices must be RDF or CDF, not Rational Field

        AUTHOR:

        - Rob Beezer (2011-03-31)
        """
        global scipy
        from sage.rings.real_double import RDF
        from sage.rings.complex_double import CDF

        cdef Matrix_double_dense Q, T

        if not self.is_square():
            raise ValueError('Schur decomposition requires a square matrix, not a {0} x {1} matrix'.format(self.nrows(), self.ncols()))
        if base_ring == None:
            base_ring = self.base_ring()
        if not base_ring in [RDF, CDF]:
            raise ValueError('base ring of Schur decomposition matrices must be RDF or CDF, not {0}'.format(base_ring))

        if self.base_ring() != base_ring:
            try:
                self = self.change_ring(base_ring)
            except TypeError:
                raise TypeError('unable to convert input matrix over CDF to a matrix over RDF')
        if base_ring == CDF:
            format = 'complex'
        else:
            format = 'real'

        schur = self.fetch('schur_factors_' + format)
        if not schur is None:
            return schur
        if scipy is None:
            import scipy
        import scipy.linalg
        Q = self._new(self._nrows, self._nrows)
        T = self._new(self._nrows, self._nrows)
        T._matrix_numpy, Q._matrix_numpy = scipy.linalg.schur(self._matrix_numpy, output=format)
        Q.set_immutable()
        T.set_immutable()
        # Our return order is the reverse of NumPy's
        schur = (Q, T)
        self.cache('schur_factors_' + format, schur)
        return schur

    def cholesky(self):
        r"""
        Return the cholesky factorization of this matrix.  The input
        matrix must be symmetric and positive definite or an exception
        will be raised.

        .. WARNING::

            This is deprecated!  Use :func:`~sage.matrix.matrix2.Matrix.cholesky_decomposition` instead.

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

        EXAMPLES::

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

        TESTS::

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

        EXAMPLES::

            sage: import numpy
            sage: a = numpy.array([[1,2],[3,4]], 'float64')
            sage: m = matrix(RDF,2,2,0)
            sage: m._replace_self_with_numpy(a)
            sage: m
            [1.0 2.0]
            [3.0 4.0]
        """
        if (<object>self._matrix_numpy).shape != (<object>numpy_matrix).shape:
            raise ValueError("matrix shapes are not the same")
        self._matrix_numpy = numpy_matrix.astype(self._numpy_dtype)


    def _replace_self_with_numpy32(self,numpy_matrix):
        """

        EXAMPLES::

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

        EXAMPLES::

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

    @rename_keyword(deprecated='Sage version 4.6', method="algorithm")
    def exp(self, algorithm='pade', order=None):
        r"""
        Calculate the exponential of this matrix X, which is the matrix

        .. math::

           e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}.

        INPUT:

        - algorithm -- 'pade', 'eig', or 'taylor'; the algorithm used to
          compute the exponential.

        - order -- for Pade or Taylor series algorithms, order is the
          order of the Pade approximation or the order of the Taylor
          series used.  The current defaults (from scipy) are 7 for
          'pade' and 20 for 'taylor'.

        EXAMPLES::

            sage: A=matrix(RDF, 2, [1,2,3,4]); A
            [1.0 2.0]
            [3.0 4.0]
            sage: A.exp()
            [51.9689561987  74.736564567]
            [112.104846851 164.073803049]
            sage: A.exp(algorithm='eig')
            [51.9689561987  74.736564567]
            [112.104846851 164.073803049]
            sage: A.exp(order=2)
            [51.8888631634 74.6198348038]
            [111.929752206 163.818615369]
            sage: A.exp(algorithm='taylor', order=5)
            [19.9583333333 28.0833333333]
            [       42.125 62.0833333333]
            sage: A.exp(algorithm='taylor')
            [51.9689035511 74.7364878369]
            [112.104731755 164.073635306]

            sage: A=matrix(CDF, 2, [1,2+I,3*I,4]); A
            [        1.0 2.0 + 1.0*I]
            [      3.0*I         4.0]
            sage: A.exp()
            [-19.6146029538 + 12.5177438468*I  3.79496364496 + 28.8837993066*I]
            [-32.3835809809 + 21.8842359579*I  2.26963300409 + 44.9013248277*I]
            sage: A.exp(algorithm='eig')
            [-19.6146029538 + 12.5177438468*I  3.79496364496 + 28.8837993066*I]
            [-32.3835809809 + 21.8842359579*I  2.26963300409 + 44.9013248277*I]
            sage: A.exp(order=2)
            [-19.6130852955 + 12.5327938535*I   3.81156364812 + 28.891438232*I]
            [-32.3827876895 + 21.9087393169*I   2.29565402142 + 44.915581543*I]
            sage: A.exp(algorithm='taylor', order=5)
            [       -6.29166666667 + 14.25*I 14.0833333333 + 15.7916666667*I]
            [               -10.5 + 26.375*I         20.0833333333 + 24.75*I]
            sage: A.exp(algorithm='taylor')
            [-19.6146006163 + 12.5177432169*I  3.79496442472 + 28.8837964828*I]
            [-32.3835771246 + 21.8842351994*I  2.26963458304 + 44.9013203415*I]
        """
        if algorithm not in ('pade', 'eig', 'taylor'):
            raise ValueError("algorithm must be 'pade', 'eig', or 'taylor'")

        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        cdef Matrix_double_dense M
        M = self._new()

        if algorithm=='pade':
            if order is None:
                M._matrix_numpy = scipy.linalg.expm(self._matrix_numpy)
            else:
                M._matrix_numpy = scipy.linalg.expm(self._matrix_numpy, q=order)
        elif algorithm=='eig':
            M._matrix_numpy = scipy.linalg.expm2(self._matrix_numpy)
        elif algorithm=='taylor':
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

        EXAMPLES::

            sage: M=matrix(CDF, [[10.234r + 34.2343jr, 34e10r]])
            sage: M
            [10.234 + 34.2343*I     3.4e+11]
            sage: M.round(2)
            [10.23 + 34.23*I        3.4e+11]
            sage: M.round()
            [10.0 + 34.0*I       3.4e+11]
        """
        global numpy
        cdef Matrix_double_dense M
        if numpy is None:
            import numpy
        ndigits = int(ndigits)
        M = self._new()
        M._matrix_numpy = numpy.round(self._matrix_numpy, ndigits)
        return M
