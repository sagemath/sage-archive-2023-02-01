"""
Dense matrices over the real double field.

Matrix operations use GSl and numpy.

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

AUTHORS:
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

include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'
include '../ext/cdefs.pxi'
include '../ext/python.pxi'

from sage.rings.real_double cimport RealDoubleElement
import sage.rings.real_double
import sage.rings.complex_double
from matrix cimport Matrix
from sage.structure.element cimport ModuleElement,Vector
from sage.modules.real_double_vector cimport RealDoubleVectorSpaceElement
from constructor import matrix
cimport sage.structure.element

cdef extern from "arrayobject.h":
#The following exposes the internal C structure of the numpy python object
# extern class [object PyArrayObject]  tells pyrex that this is
# a compiled python class defined by the C struct PyArrayObject
    cdef enum:
        NPY_OWNDATA = 0x0004 #bit mask so numpy does not free array contents when its destroyed

    ctypedef int intp

##     ctypedef extern class numpy.dtype [object PyArray_Descr]:
##         cdef int type_num, elsize, alignment
##         cdef char type, kind, byteorder, hasobject
##         cdef object fields, typeobj

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
#        cdef object base
#        cdef dtype descr
        cdef int flags

    object PyArray_FromDims(int,int *,int)
    object PyArray_FromDimsAndData(int,int*,int,double *)
    void import_array()


cdef class Matrix_real_double_dense(matrix_dense.Matrix_dense):   # dense
    """
    Class that implements matrices over the real double field. These
    are supposed to be fast matrix operations using C doubles. Most
    operations are implemented using GSl or numpy libraries which will
    call the underlying BLAS on the system.

    EXAMPLES:
        sage: m = Matrix(RDF, [[1,2],[3,4]])
        sage: m**2
        [ 7.0 10.0]
        [15.0 22.0]
        sage: n= m^(-1); n
        [-2.0  1.0]
        [ 1.5 -0.5]

    To compute eigenvalues the use the functions left eigen_vectors or right_eigenvectors

        sage: p,e = m.right_eigenvectors()

    the result of eigen is a pair p,e, where p is a list
    of eigenvalues and the e is a matrix whose columns are the eigenvectors

    To solve a linear system Ax = b
    for A = [[1,2]  and b = [5,6]
             [3,4]]

        sage: b = vector(RDF,[5,6])
        sage: m.solve_left(b)
        (-4.0, 4.5)

    See the commands QR,LU,SVD for QR, LU, and singular value
    decomposition.
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
        matrix_dense.Matrix_dense.__init__(self,parent)
        if self._nrows == 0 or self._ncols == 0:
            self._matrix = NULL
            return
        _sig_on
        self._matrix= <gsl_matrix *> gsl_matrix_calloc(self._nrows, self._ncols)
        _sig_off
        if self._matrix == NULL:
            raise MemoryError, "unable to allocate memory for matrix "
        self._LU = <gsl_matrix *> NULL
        self._p = <gsl_permutation *> NULL


    def __dealloc__(self):
        if self._matrix == NULL:
            return
        gsl_matrix_free(self._matrix)
        if self._LU != NULL:
            gsl_matrix_free(self._LU)
        if self._p !=NULL:
            gsl_permutation_free(self._p)


    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    def __hash__(self):
        """
        Hash this matrix, if it's immutable.

        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
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
        if self._matrix == NULL:
            return
        cdef double z
        cdef Py_ssize_t i,j
        if isinstance(entries,list):
            if len(entries)!=self._nrows*self._ncols:
                    raise TypeError, "entries has wrong length"

            if coerce:

                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        z= float(entries[i*self._ncols+j])
                        gsl_matrix_set(self._matrix, i,j,z)

            else:

                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        gsl_matrix_set(self._matrix, i,j,entries[i*self._ncols +j])


        else:
            if entries is None:
                z = 0.0
            else:
                try:
                    z = float(entries)
                    gsl_matrix_set_zero(self._matrix)
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or float"
            if z != 0:
                if self._nrows != self._ncols:
                    raise TypeError, "scalar matrix must be square"
                for i from 0<=i<self._ncols:
                    gsl_matrix_set(self._matrix,i,i,z)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        cdef double z
        z = float(value)
        gsl_matrix_set(self._matrix,i,j,z)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return RealDoubleElement(gsl_matrix_get(self._matrix,i,j))


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        if self._nrows == 0 or self._ncols == 0: return self
        cdef Matrix_real_double_dense M,_right,_left
        _right = right
        _left = self
        cdef int result_add,result_copy
        if (self._matrix.size1 != _right._matrix.size1 and self._matrix.size2 != _right._matrix.size2):
            raise TypeError, "Cannot add matrices if they have different dimensions"
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_memcpy(M._matrix,_left._matrix)
        result_add = gsl_matrix_add(M._matrix,_right._matrix)
        if result_copy!=GSL_SUCCESS or result_add !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
        # todo -- check error code
        return M


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        EXAMPLES:
            sage: A = matrix(RDF,3,range(1,10))
            sage: (A-A).is_zero()
            True
        """
        if self._nrows == 0 or self._ncols == 0: return self

        cdef Matrix_real_double_dense M,_right,_left
        _right = right
        _left = self
        cdef int result_sub,result_copy
        if (self._matrix.size1 != _right._matrix.size1 and self._matrix.size2 != _right._matrix.size2):
            raise TypeError, "Cannot subtract matrices if they have different dimensions"
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        # todo -- check error code
        result_copy = gsl_matrix_memcpy(M._matrix,_left._matrix)
        result_sub = gsl_matrix_sub(M._matrix,_right._matrix)
        if result_copy!=GSL_SUCCESS or result_sub !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
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
        if self._nrows == 0 or self._ncols == 0: return self

        cdef Matrix_real_double_dense M
        cdef int result_neg, result_copy
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_memcpy(M._matrix,self._matrix)
        result_neg = gsl_matrix_scale(M._matrix,-1.0)
        if result_copy!=GSL_SUCCESS or result_neg !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
        return M



    #   * cdef _cmp_c_impl
    #   * __copy__
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):                        #unsure how to implement
    # def _unpickle(self, data, int version):   # use version >= 0 #unsure how to implement
    ######################################################################
    cdef sage.structure.element.Matrix _matrix_times_matrix_c_impl(self, sage.structure.element.Matrix right):
        cdef int result
        if self._ncols!=right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right"

        parent = self.matrix_space(self._nrows,right._ncols)
        if self._nrows == 0 or self._ncols == 0: return parent.zero_matrix()

        cdef Matrix_real_double_dense M,_right,_left
        _right = right
        _left = self

        M = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result  = gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,_left._matrix,_right._matrix,0,M._matrix)
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
            sage: A.inverse().det()
            nan
            sage: A = matrix(RDF,3,range(1,10));A
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            [7.0 8.0 9.0]

            sage: A.determinant() < 10e-12
            True
            sage: ~A              # slightly random
            [-4.50359962737e+15  9.00719925474e+15 -4.50359962737e+15]
            [ 9.00719925474e+15 -1.80143985095e+16  9.00719925474e+15]
            [-4.50359962737e+15  9.00719925474e+15 -4.50359962737e+15]
        """
        if self._nrows == 0 or self._ncols == 0: return self

        cdef int result_LU, result_invert
        if self.fetch('LU_valid') != True:
            self._c_compute_LU()
        cdef Matrix_real_double_dense M
        parent = self.matrix_space(self._nrows,self._ncols)
        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result_invert = gsl_linalg_LU_invert(self._LU,self._p,M._matrix)
        return M

    # def __copy__(self):
    # def _list(self):
    # def _dict(self):


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
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
        cdef int result_LU
        if self._LU == NULL:
            self._LU = <gsl_matrix *> gsl_matrix_alloc(self._nrows,self._ncols)
        if self._LU == NULL:
            raise MemoryError, "allocation error"
        if self._p ==NULL:
            self._p =<gsl_permutation *> gsl_permutation_alloc(self._nrows)
        if self._p == NULL:
            raise MemoryError, "allocation error"
        gsl_matrix_memcpy(self._LU,self._matrix)
        _sig_on
        result_LU = gsl_linalg_LU_decomp(self._LU,self._p,&self._signum)
        _sig_off
        if result_LU == GSL_SUCCESS:
            self.cache('LU_valid',True)
        else:
            raise ValueError,"Error computing LU decomposition"


    def LU(self):
        """
        Computes the LU decomposition of a matrix.

        For and square matrix A we can find matrices P,L, and U. s.t.
           P*A = L*U
        where P is a permutation matrix, L is lower triangular and U is upper triangular.

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
            return self, self, self

        if self.fetch('LU_valid')!=True:
            self._c_compute_LU()
        cdef Py_ssize_t i,j,k,l,copy_result
        cdef Matrix_real_double_dense P, L,U
        parent = self.matrix_space(self._nrows,self._ncols)
        P=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        L = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        U = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        for i from 0<=i<self._ncols:
            j = gsl_permutation_get(self._p,i)
            P.set_unsafe(i,j,1)
            L.set_unsafe(i,i,1)
            U.set_unsafe(i,i,gsl_matrix_get(self._LU,i,i))
            for l from 0<=l<i:
                L.set_unsafe(i,l,gsl_matrix_get(self._LU,i,l))
                U.set_unsafe(l,i,gsl_matrix_get(self._LU,l,i))


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

        IMPLEMENTATION:
            Uses numpy.
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self

        import numpy
        v, m = numpy.linalg.eig(numpy.transpose(self.numpy()))
        return ([sage.rings.complex_double.CDF(x) for x in v],matrix(m).transpose())

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
            sage: m.right_eigenvectors()           # random-ish platform-dependent output (low order digits)

        IMPLEMENTATION:
            Uses numpy.
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self

        import numpy
        v, m = numpy.linalg.eig(self.numpy())
        return ([sage.rings.complex_double.CDF(x) for x in v],matrix(m))

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
        import solve
        return solve.solve_matrix_real_double_dense(self, b)

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
        cdef RealDoubleVectorSpaceElement _vec, ans
        cdef double *p
        cdef ndarray _result

        M = self._column_ambient_module()
        ans = M.zero_vector()
        if self._nrows == 0 or self._ncols == 0:
            return ans


        _vec=vec

        import numpy
        _result=numpy.linalg.solve(self.numpy(),_vec.numpy())

        p = <double *>_result.data

        memcpy(ans.v.data, _result.data, _result.dimensions[0]*sizeof(double))

        return ans

    def determinant(self):
        """
        Return the determinant of self.

        ALGORITHM: Use GSL (LU decompositon)

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
            return RealDoubleElement(1)
        if not self.is_square():
            raise ValueError, "self must be square"
        if self.fetch('LU_valid')!=True:
            self._c_compute_LU()
        return RealDoubleElement(gsl_linalg_LU_det(self._LU, self._signum))

    def log_determinant(self):
        """
        Compute the log of the absolute value of the determinant
        using GSL (LU decomposition).

        NOTE: This is useful if the usual determinant overlows.

        EXAMPLES:
            sage: m = matrix(RDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: RDF(log(abs(m.determinant())))
            0.69314718056
            sage: m.log_determinant()
            0.69314718056
            sage: m = matrix(RDF,0,0,range(4)); m
            []
            sage: m.log_determinant()
            0.0
        """
        if self._nrows == 0 or self._ncols == 0:
            return RealDoubleElement(0)
        if self.fetch('LU_valid')!=True:
            self._c_compute_LU()
        return RealDoubleElement(gsl_linalg_LU_lndet(self._LU))

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
        cdef Matrix_real_double_dense trans
        cdef int result_copy
        parent  = self.matrix_space(self._ncols,self._nrows)
        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._ncols, self._nrows)
        trans = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_transpose_memcpy(trans._matrix,self._matrix)
        if result_copy !=GSL_SUCCESS:
            raise ValueError, "Error copy matrix"
        if self.subdivisions is not None:
            row_divs, col_divs = self.get_subdivisions()
            trans.subdivide(col_divs, row_divs)
        return trans

    def SVD(self, algorithm='gsl'):
        r"""
        Return the singular value decomposition of this matrix.

        INPUT:
            A -- a matrix
            algorithm -- 'numpy' or 'gsl'
        OUTPUT:
            U, S, V -- matrices such that $A = U S V^t$, where
                       U and V are orthogonal and S is diagonal.


        EXAMPLES:
            sage: m = matrix(RDF,4,range(16))
            sage: U,S,V = m.SVD(algorithm='gsl')
            sage: U*S*V.transpose()    # slightly random output (due to computer architecture)
            [3.45569519412e-16               1.0               2.0               3.0]
            [4.0               5.0               6.0               7.0]
            [8.0               9.0              10.0              11.0]
            [12.0              13.0              14.0              15.0]

        A non-square example:
            sage: m = matrix(RDF, 2, range(6)); m
            [0.0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: U, S, V = m.SVD(algorithm='numpy')
            sage: U*S*V.transpose()           # random low bits
            [7.62194127257e-17               1.0               2.0]
            [              3.0               4.0               5.0]
            sage: U, S, V = m.SVD()
            sage: U
            [-0.274721127897 -0.961523947641]
            [-0.961523947641  0.274721127897]
            sage: S
            [7.34846922835           0.0]
            [          0.0           1.0]
            sage: V
            [-0.392540507864  0.824163383692]
            [-0.560772154092  0.137360563949]
            [ -0.72900380032 -0.549442255795]
            sage: U*S*V.transpose()           # random low bits
            [7.62194127257e-17               1.0               2.0]
            [              3.0               4.0               5.0]
            sage: m = matrix(RDF,3,2,range(6)); m
            [0.0 1.0]
            [2.0 3.0]
            [4.0 5.0]
            sage: U,S,V = m.SVD(algorithm='numpy')
            sage: U*S*V.transpose()   # random low order bits
            [-8.13151629364e-19                1.0]
            [               2.0                3.0]
            [               4.0                5.0]
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()   # random low order bits
            [-8.13151629364e-19                1.0]
            [               2.0                3.0]
            [               4.0                5.0]


        TESTS:
            sage: m = matrix(RDF, 3, 0, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: m = matrix(RDF, 0, 3, []); m
            []
            sage: m.SVD()
            ([], [], [])
        """
        if algorithm == 'numpy':
            return self._SVD_numpy()
        elif algorithm == 'gsl':
            return self._SVD_gsl()
        else:
            raise ValueError, "unknown algorithm"

    def _SVD_gsl(self):
        """
        Return the singular value decomposition of this matrix
        using GSL.  Note that the matrices the dimensions of the
        matrices that this returns are (m,p), (p,p), and (n, p)
        where p = min(m,n).

        EXAMPLES:
            sage: def shape(x): return (x.nrows(), x.ncols())
            sage: m = matrix(RDF, 2, 3, range(6))
            sage: map(shape, m._SVD_gsl())
            [(2, 2), (2, 2), (3, 2)]

        """
        if self._ncols > self._nrows:
            m = self.transpose()
            V_t,S_t,U_t = m.SVD()
            return U_t,S_t,V_t

        if self._nrows == 0 or self._ncols == 0:
            U_t = self.new_matrix(self._nrows, self._ncols)
            S_t = self.new_matrix(self._ncols, self._ncols)
            V_t = self.new_matrix(self._ncols, self._nrows)
            return U_t, S_t, V_t

        cdef Matrix_real_double_dense A,V,_S
        cdef gsl_vector* S
        cdef gsl_vector* work_space
        cdef int result_copy, result_svd, i
        parent_A = self.matrix_space(self._nrows,self._ncols)
        A=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent_A,None,None,None)
        parent_V = self.matrix_space(self._ncols,self._ncols)
        V = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent_V,None,None,None)
        result_copy = gsl_matrix_memcpy(A._matrix,self._matrix)
        S = <gsl_vector *> gsl_vector_alloc(self._ncols)
        work_space = <gsl_vector *> gsl_vector_alloc(self._ncols)
        _sig_on
        result_svd  = gsl_linalg_SV_decomp(A._matrix, V._matrix, S, work_space)
        _sig_off
        parent_S = self.matrix_space(self._ncols,self._ncols)
        _S = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent_S,None,None,None)
        for i from 0<=i<self._ncols:
            _S.set_unsafe(i,i,gsl_vector_get(S,i))
        gsl_vector_free(S)
        gsl_vector_free(work_space)
        return A,_S,V

    def _SVD_numpy(self):
        """
        Return the singular value decomposition of this matrix
        using GSL.  Note that the matrices the dimensions of the
        matrices that this returns are (m,m), (m,n), and (n, n).

        EXAMPLES:
            sage: def shape(x): return (x.nrows(), x.ncols())
            sage: m = matrix(RDF, 2, 3, range(6))
            sage: map(shape, m._SVD_numpy())
            [(2, 2), (2, 3), (3, 3)]

        """
        if self._nrows == 0 or self._ncols == 0:
            U_t = self.new_matrix(self._nrows, self._ncols)
            S_t = self.new_matrix(self._nrows, self._ncols)
            V_t = self.new_matrix(self._ncols, self._nrows)
            return U_t, S_t, V_t

        import numpy.linalg
        cdef int i, s_dim
        P = self.parent()
        CDF = P.base_ring()

        U,_S,V = numpy.linalg.svd(self.numpy())

        #Create the inner diagonal matrix
        s_dim = len(_S)
        S = matrix(CDF, self._nrows, self._ncols, 0)
        for i from 0 <= i < s_dim:
            S[(i,i)] = _S[i]

        return (matrix(U),S,matrix(V).transpose())





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

        cdef gsl_matrix* A
        cdef gsl_vector* v
        cdef Matrix_real_double_dense Q,R

        A = <gsl_matrix *> gsl_matrix_alloc(self._nrows,self._ncols)
        v = <gsl_vector *> gsl_vector_alloc(min(self._nrows,self._ncols))
        cdef int result
        result = gsl_matrix_memcpy(A, self._matrix)
        if result !=GSL_SUCCESS:
            gsl_matrix_free(A)
            gsl_vector_free(v)
            raise ValueError,"Error copying"
        _sig_on
        result = gsl_linalg_QR_decomp(A,v)
        _sig_off
        parent_Q = self.matrix_space(self._nrows,self._nrows)
        parent_R = self.matrix_space(self._nrows,self._ncols)
        Q = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent_Q,None,None,None)
        R = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent_R,None,None,None)
        _sig_on
        result = gsl_linalg_QR_unpack(A,v,Q._matrix,R._matrix)
        _sig_off
        gsl_matrix_free(A)
        gsl_vector_free(v)
        if result!=GSL_SUCCESS:
            raise ValueError,"Error computing QR factorization"
        return Q,R

    def _cholesky_gsl(self,check_symmetry=True):
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:   # special case
            return self

        cdef gsl_matrix *A
        cdef Matrix_real_double_dense result
        cdef int i,j
        if check_symmetry:
            if not self.is_symmetric():
                raise TypeError,"cannot take cholesky decomposition of non-symmetric matrix"
        A=<gsl_matrix*>gsl_matrix_alloc(self._nrows,self._ncols)
        gsl_matrix_memcpy(A,self._matrix)
        _sig_on
        gsl_linalg_cholesky_decomp(A)
        _sig_off
        parent = self.matrix_space(self._nrows,self._ncols)
        result=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        for i from 0<i<self._nrows:
            for j from 0<=j<i:
                gsl_matrix_set(A,j,i,0)
        gsl_matrix_memcpy(result._matrix,A)
        gsl_matrix_free(A)
        return result


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
                if fabs(self.get_unsafe(i,j) - self.get_unsafe(j,i)) > tol:
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

            sage: m = r*M(1)*r.transpose()
            sage: c = m.cholesky()
            sage: c*c.transpose()
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:   # special case
            return self

        cdef Matrix_real_double_dense _M,_result_matrix
        _M=self
        cdef int i
        cdef double *p
        cdef ndarray _n,_m
        import numpy
        _m = numpy.linalg.cholesky(self.numpy())
        parent = self.matrix_space(self._nrows,self._ncols)
        _result_matrix = Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        p = <double *> _m.data
        for i from 0<=i<_M._matrix.size1*_M._matrix.size2:
            _result_matrix._matrix.data[i] = p[i]
        if numpy.isfortran(_m)==True:
            return _result_matrix.transpose()
        return _result_matrix

    cdef Vector _vector_times_matrix_c_impl(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._row_ambient_module().zero_vector()
        cdef RealDoubleVectorSpaceElement v_,ans
        cdef gsl_vector *vec
        cdef Py_ssize_t i,j
        M=self._row_ambient_module()
        v_ = v
        ans=M.zero_vector()
        gsl_blas_dgemv(CblasTrans,1.0,self._matrix, v_.v,0.0,ans.v)
        return ans



    cdef Vector _matrix_times_vector_c_impl(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._column_ambient_module().zero_vector()

        cdef RealDoubleVectorSpaceElement v_,ans
        cdef gsl_vector *vec
        cdef Py_ssize_t i,j
        M=self._column_ambient_module()
        v_ = v
        ans=M.zero_vector()
        gsl_blas_dgemv(CblasNoTrans,1.0,self._matrix, v_.v,0.0,ans.v)
        return ans


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
        if self._ncols == 0:
            import numpy
            return numpy.array([[]]*self._nrows)
        elif self._nrows == 0:
            import numpy
            return numpy.array([[]]*self._ncols).transpose()
        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef Matrix_real_double_dense _M,_result_matrix
        _M=self
        cdef int dims[2]
        cdef double * data
        cdef int i
        cdef object temp
        cdef double *p
        cdef ndarray _n,_m
        dims[0] = _M._matrix.size1
        dims[1] = _M._matrix.size2
        data = <double*> malloc(sizeof(double)*dims[0]*dims[1])
        memcpy(data,_M._matrix.data,sizeof(double)*dims[0]*dims[1])
        temp = PyArray_FromDimsAndData(2, dims, 12,data)
        _n = temp
        _n.flags = _n.flags|(NPY_OWNDATA) # this sets the ownership bit
        return _n

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
        if self._nrows == 0 or self._ncols == 0:
            return
        cdef ndarray n
        cdef double *p
        n = numpy_matrix
        p=<double *>n.data
        memcpy(self._matrix.data,p,sizeof(double)*self._nrows*self._ncols)


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
        if self._nrows == 0 or self._ncols == 0:
            return
        cdef ndarray n
        cdef float *nd
        cdef double *md
        cdef int i

        n = numpy_matrix
        nd = <float *>n.data
        md = <double *>self._matrix.data
        for i from 0 <= i < self._nrows*self._ncols:
            md[i] = <double> nd[i]



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
                s += gsl_matrix_get(self._matrix,i,j) * gsl_matrix_get(self._matrix,i,j)
            d += gsl_sf_log(s)
        d /= 2
        return int(math.ceil(d / gsl_sf_log(10)))
