"""
Dense matrices over the Real Double Field. Matrix operations use GSl and numpy.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################
include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'
include '../ext/cdefs.pxi'
include '../ext/python.pxi'
from sage.rings.real_double cimport RealDoubleElement
import sage.rings.real_double
import numpy
from matrix cimport Matrix
from sage.structure.element cimport ModuleElement
cdef extern from "arrayobject.h":
#The following exposes the internal C structure of the numpy python object
# extern class [object PyArrayObject]  tells pyrex that this is
# a compiled python class defined by the C struct PyArrayObject
    cdef enum:
        NPY_OWNDATA = 0x0004 #bit mask so numpy does not free array contents when its destroyed

    ctypedef int intp

    ctypedef extern class numpy.dtype [object PyArray_Descr]:
        cdef int type_num, elsize, alignment
        cdef char type, kind, byteorder, hasobject
        cdef object fields, typeobj

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef object base
        cdef dtype descr
        cdef int flags

    object PyArray_FromDims(int,int *,int)
    void import_array()


cdef class Matrix_real_double_dense(matrix_dense.Matrix_dense):   # dense
    """Class that implements matrices over the real double field. These are
    supposed to be fast matrix operations using C doubles. Most operations
    are implemented using GSl or numpy libraries which will call the underlying
    BLAS on the system.

    Examples:

    sage: m = Matrix(RDF, [[1,2],[3,4]])
    sage: m**2
    [ 7.0 10.0]
    [15.0 22.0]
    sage: n= m^(-1); n
    [-2.0  1.0]
    [ 1.5 -0.5]

    To compute eigenvalues the use the function eigen

    sage: p,e = m.eigen()

    the result of eigen is a pair p,e . p is a list
    of eigenvalues and the e is a matrix whose columns are the eigenvectors.
    currently e is a numpy matrix and p is a numpy vector. To make it a RDF
    matrix and vector

    sage: p = vector(RDF,list(p))
    sage: e = Matrix(RDF,list(e))

    To solve a linear system Ax = b
    for A = [[1,2]  and b = [5,6]
             [3,4]]

    sage: b = vector(RDF,[5,6])
    sage: m.solve_left(b)
    (-4.0, 4.5)

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
        self._matrix= <gsl_matrix *> gsl_matrix_calloc(self._nrows, self._ncols)
        if self._matrix == NULL:
            raise MemoryError, "unable to allocate memory for matrix "
        self._LU = <gsl_matrix *> NULL
        self._p = <gsl_permutation *> NULL
        self._LU_valid = 0

    def __dealloc__(self):
        gsl_matrix_free(self._matrix)
        if self._LU != NULL:
            gsl_matrix_free(self._LU)
        if self._p !=NULL:
            gsl_permutation_free(self._p)


    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()


    def __init__(self, parent, entries, copy, coerce):
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
            try:
                z=float(entries)
            except TypeError:
                raise TypeError, "entries must to coercible to list or real double "
            if self._nrows != self._ncols:
                raise TypeError, "scalar matrix must be square"
            for i from 0<=i<self._ncols:
                gsl_matrix_set(self._matrix,i,i,z)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        cdef double z
        z = float(value)
        gsl_matrix_set(self._matrix,i,j,z) #sig on here ?
        self._LU_valid  = 0
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return sage.rings.real_double.RDF(gsl_matrix_get(self._matrix,i,j)) #sig on here?


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    cdef ModuleElement _add_c_impl(self, ModuleElement right):
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


    cdef ModuleElement _sub_c_impl(self, ModuleElement right): #matrix.Matrix right):
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
        result_add = gsl_matrix_sub(M._matrix,_right._matrix)
        if result_copy!=GSL_SUCCESS or result_sub !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
        return M

    def __neg__(self):
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
    def _multiply_classical(self, matrix.Matrix right):
        cdef int result
        if self._ncols!=right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right"

        parent = self.matrix_space(self._nrows,right._ncols)
        cdef Matrix_real_double_dense M,_right,_left
        _right = right
        _left = self

        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result  = gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,_left._matrix,_right._matrix,0,M._matrix)
        return M

    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    def __invert__(self):
        cdef int result_LU, result_invert
        if(self._LU_valid != 1):
            self._c_compute_LU()
        cdef Matrix_real_double_dense M
        parent = self.matrix_space(self._nrows,self._ncols)
        M=Matrix_real_double_dense.__new__(Matrix_real_double_dense,parent,None,None,None)
        result_invert = gsl_linalg_LU_invert(self._LU,self._p,M._matrix)
        self._LU_valid = 1
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
        result_LU = gsl_linalg_LU_decomp(self._LU,self._p,&self._signum)
        if result_LU == GSL_SUCCESS:
            self._LU_valid = 1
        else:
            raise ValueError,"Error computing LU decomposition"


    def eigen(self):
        """
        Computes the eigenvalues and eigenvectors of this matrix:

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as a list

        These are still formated via numpy, but this will change.

        EXAMPLES:
            sage: m = Matrix(RDF, 3, range(9))
            sage: m.eigen()           # random-ish platform-dependent output (low order digits)
            (array([  1.33484692e+01,  -1.34846923e+00,  -4.69022010e-16]),
             array([[-0.16476382, -0.79969966,  0.40824829],
                   [-0.50577448, -0.10420579, -0.81649658],
                   [-0.84678513,  0.59128809,  0.40824829]]))

        IMPLEMENTATION:
            Uses numpy.
        """
        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef Matrix_real_double_dense _M
        _M=self
        cdef int dims[2]
        cdef object temp
        cdef ndarray _n
        dims[0] = _M._matrix.size1
        dims[1] = _M._matrix.size2
        temp = PyArray_FromDims(2, dims, 12)#, char_pointer)
        _n = temp                           # this isn't quite working yet so we invalidate the pointer at the end
        _n.flags = _n.flags&(~NPY_OWNDATA) # this perform as a logical AND on NOT(NPY_OWNDATA), which sets that bit to 0
        _n.data = <char *> _M._matrix.data #numpy arrays store their data as char *
        result = numpy.linalg.eig(_n)
        _n.data = <char *> NULL    #keep numpy from deallocating memory
        return result   #todo: make the result a real double matrix

    def solve_left(self, vec):
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
        """
        import solve
        return solve.solve_matrix_real_double_dense(self, vec)



