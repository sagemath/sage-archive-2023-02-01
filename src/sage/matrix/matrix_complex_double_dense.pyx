"""
Dense matrices over the Complex Double Field.

IMPLEMENTATION:
   Specialized matrix operations use GSL and numpy.


EXAMPLES:
    sage: b = Mat(CDF,2,3).basis()
    sage: b[0]
    [1.0   0   0]
    [  0   0   0]

TESTS:
    sage: a = matrix(CDF,2,range(4), sparse=False)
    sage: loads(dumps(a)) == a
    True
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
#include '../gsl/gsl.pxi'
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.modules.complex_double_vector cimport ComplexDoubleVectorSpaceElement
import sage.rings.complex_double
import sage.rings.real_double
from matrix cimport Matrix
from sage.structure.element cimport ModuleElement, Vector
from constructor import matrix
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


cdef class Matrix_complex_double_dense(matrix_dense.Matrix_dense):   # dense
    """Class that implements matrices over the complex double field. These are
    supposed to be fast matrix operations using C doubles. Most operations
    are implemented using GSl or numpy libraries which will call the underlying
    BLAS on the system.

    Examples:

    sage: m = Matrix(CDF, [[1,2*I],[3+I,4]])
    sage: m**2
    [-1.0 + 6.0*I       10.0*I]
    [15.0 + 5.0*I 14.0 + 6.0*I]

    sage: n= m^(-1); n
    [  0.333333333333 + 0.333333333333*I   0.166666666667 - 0.166666666667*I]
    [ -0.166666666667 - 0.333333333333*I 0.0833333333333 + 0.0833333333333*I]

    To compute eigenvalues the use the functions left_eigenvectors or right_eigenvectors

    sage: p,e = m.right_eigenvectors()

    the result of eigen is a pair p,e . p is a list
    of eigenvalues and the e is a matrix whose columns are the eigenvectors.



    To solve a linear system Ax = b
    for A = [[1.0,2*I]  and b = [1,I]
             [3+I,4]]

    sage: b = vector(CDF,[1,I])
    sage: m.solve_left(b)
    (0.5 + 0.5*I, -0.25 - 0.25*I)

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
        self._matrix= <gsl_matrix_complex*> gsl_matrix_complex_calloc(self._nrows, self._ncols)
        _sig_off
        if self._matrix == NULL:
            raise MemoryError, "unable to allocate memory for matrix "
        self._LU = <gsl_matrix_complex *> NULL
        self._p = <gsl_permutation *> NULL


    def __dealloc__(self):
        if self._matrix == NULL:
            return
        gsl_matrix_complex_free(self._matrix)
        if self._LU != NULL:
            gsl_matrix_complex_free(self._LU)
        if self._p !=NULL:
            gsl_permutation_free(self._p)


    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()


    def __init__(self, parent, entries, copy, coerce):
        cdef ComplexDoubleElement z
        cdef Py_ssize_t i,j
        if self._matrix == NULL:
            return
        if isinstance(entries,list):
            if len(entries)!=self._nrows*self._ncols:
                    raise TypeError, "entries has wrong length"

            if coerce:

                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        z= sage.rings.complex_double.CDF(entries[i*self._ncols+j]) #better way to do this?
                        gsl_matrix_complex_set(self._matrix, i,j,z._complex)

            else: # Do I need to coerce here as well or assume CDF already.

                for i from 0<=i<self._nrows:
                    for j from 0<=j<self._ncols:
                        z = sage.rings.complex_double.CDF(entries[i*self._ncols+j])
                        gsl_matrix_complex_set(self._matrix, i,j,z._complex)


        else:
            try:
                z = sage.rings.complex_double.CDF(entries)
            except TypeError:
                raise TypeError, "entries must to coercible to list or complex double "
            if z != 0:
                if self._nrows != self._ncols:
                    raise TypeError, "scalar matrix must be square"
                for i from 0<=i<self._ncols:
                    gsl_matrix_complex_set(self._matrix,i,i,z._complex)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        cdef ComplexDoubleElement z
        z = sage.rings.complex_double.CDF(value)   # do I assume value is already CDF
        gsl_matrix_complex_set(self._matrix,i,j,z._complex) #sig on here ?

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef gsl_complex z
        z= gsl_matrix_complex_get(self._matrix,i,j)
        return sage.rings.complex_double.CDF(GSL_REAL(z),GSL_IMAG(z))


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    cpdef ModuleElement _add_(self, ModuleElement right):
        if self._nrows == 0 or self._ncols == 0: return self
        cdef Matrix_complex_double_dense M,_right,_left
        _right = right
        _left = self
        cdef int result_add,result_copy
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_complex_memcpy(M._matrix,_left._matrix)
        result_add = gsl_matrix_complex_add(M._matrix,_right._matrix)
        if result_copy!=GSL_SUCCESS or result_add !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
        # todo -- check error code
        return M


    cpdef ModuleElement _sub_(self, ModuleElement right): #matrix.Matrix right):
        if self._nrows == 0 or self._ncols == 0: return self

        cdef Matrix_complex_double_dense M,_right,_left
        _right = right
        _left = self
        cdef int result_sub,result_copy
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        # todo -- check error code
        result_copy = gsl_matrix_complex_memcpy(M._matrix,_left._matrix)
        result_add = gsl_matrix_complex_sub(M._matrix,_right._matrix)
        if result_copy!=GSL_SUCCESS or result_sub !=GSL_SUCCESS:
            raise ValueError, "GSL routine had an error"
        return M

    def __neg__(self):
        if self._nrows == 0 or self._ncols == 0: return self

        cdef Matrix_complex_double_dense M
        cdef int result_neg, result_copy
        cdef gsl_complex z
        GSL_SET_COMPLEX(&z,-1.0,0)
        parent = self.matrix_space(self._matrix.size1,self._matrix.size2)
        M=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_complex_memcpy(M._matrix,self._matrix)
        result_neg = gsl_matrix_complex_scale(M._matrix,z)
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
        cdef gsl_complex a, b
        if self._ncols!=right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right"

        parent = self.matrix_space(self._nrows,right._ncols)
        if self._nrows == 0 or self._ncols == 0: return parent.zero_matrix()

        cdef Matrix_complex_double_dense M,_right,_left
        _right = right
        _left = self
        GSL_SET_COMPLEX(&a,1.0,0)
        GSL_SET_COMPLEX(&b,0,0)
        M=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        result  = gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,a,_left._matrix,_right._matrix,b,M._matrix)
        return M

    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    def __invert__(self):
        if self._nrows == 0 or self._ncols == 0: return self

        cdef int result_LU, result_invert
        if self.fetch('LU_valid') != True:
            self._c_compute_LU()
        cdef Matrix_complex_double_dense M
        parent = self.matrix_space(self._nrows,self._ncols)
        M=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        result_invert = gsl_linalg_complex_LU_invert(self._LU,self._p,M._matrix)
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
        cdef int result_LU
        if self._LU == NULL:
            self._LU = <gsl_matrix_complex *> gsl_matrix_complex_alloc(self._nrows,self._ncols)
        if self._LU == NULL:
            raise MemoryError, "allocation error"
        if self._p ==NULL:
            self._p =<gsl_permutation *> gsl_permutation_alloc(self._nrows)
        if self._p == NULL:
            raise MemoryError, "allocation error"
        if self._nrows == 0 or self._ncols == 0:
            raise RuntimeError, "no LU in 0 rows or 0 columns case"
        gsl_matrix_complex_memcpy(self._LU,self._matrix)
        _sig_on
        result_LU = gsl_linalg_complex_LU_decomp(self._LU,self._p,&self._signum)
        _sig_off
        if result_LU == GSL_SUCCESS:
            self.cache('LU_valid',True)
        else:
            raise ValueError,"Error computing LU decomposition"


    def left_eigenvectors(self):
        """
        Computes the eigenvalues and *left* eigenvectors of this
        matrix m acting *from the right*.  I.e., vectors v such that
        v*m = lambda*v.

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as an CDF matrix whose rows
                           are the eigenvectors.

        EXAMPLES:
            sage: I = CDF.gen()
            sage: m = I*Matrix(CDF, 3, range(9))
            sage: vals, vecs = m.left_eigenvectors()
            sage: vecs*m      # random precision
            [      1.83356748856e-16 + 5.8765683663*I      3.11202346481e-16 + 7.58017348024*I      4.39047944106e-16 + 9.28377859418*I]
            [     1.66533453694e-15 - 1.21076184124*I     1.72084568817e-15 - 0.375459730779*I      1.7763568394e-15 + 0.459842379686*I]
            [-3.33066907388e-16 - 6.66133814775e-16*I -5.55111512313e-16 - 7.21644966006e-16*I -7.77156117238e-16 - 7.77156117238e-16*I]

            sage: vals[0] * vecs.row(0)  # random precision
            (1.67272748078e-15 + 5.8765683663*I, 1.82458083857e-15 + 7.58017348024*I, 1.23544530231e-15 + 9.28377859418*I)
            sage: vals[1] * vecs.row(1)  # random precision
            (-2.00045384696e-16 - 1.21076184124*I, 5.36806149949e-16 - 0.375459730779*I, -5.97719297047e-16 + 0.459842379686*I)
            sage: vals[2] * vecs.row(2)  # random precision
            (-2.0580841304e-16 + 1.41000047373e-16*I, 4.11616826079e-16 - 2.82000094746e-16*I, -2.0580841304e-16 + 1.41000047373e-16*I)
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self

        import numpy
        vals,vecs = numpy.linalg.eig(numpy.transpose(self.numpy()))
        return ([sage.rings.complex_double.CDF(x) for x in vals],matrix(vecs).transpose())
#        return vals, vecs


    def right_eigenvectors(self):
        """
        Computes the eigenvalues and *right* eigenvectors of this
        matrix m acting *from the left*.  I.e., vectors v such that
        m * v = lambda*v, where v is viewed as a column vector.

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as an CDF matrix whose columns
                           are the eigenvectors.

        EXAMPLES:
            sage: I = CDF.gen()
            sage: m = I*Matrix(CDF, 3, range(9))
            sage: vals, vecs = m.right_eigenvectors()
            sage: m*vecs      # random precision

            [    -2.18763583539e-17 + 2.19934474494*I      2.22220787778e-16 - 1.07837038763*I                       -1.28830323146e-16]
            [    -5.79030412929e-16 + 6.75131502804*I     2.46412048752e-16 - 0.140518298155*I -3.67002435386e-16 + 1.66533453694e-16*I]
            [     -1.1361844675e-15 + 11.3032853111*I     2.70603309725e-16 + 0.797333791317*I -6.05174547627e-16 + 3.33066907388e-16*I]
            sage: vals[0] * vecs.column(0)  # random precision
            (-2.45441964831e-15 + 2.19934474494*I, -1.11280381694e-15 + 6.75131502804*I, -1.37419154761e-15 + 11.3032853111*I)
            sage: vals[1] * vecs.column(1)  # random precision
            (3.94205510736e-18 - 1.07837038763*I, 2.78424120914e-16 - 0.140518298155*I, -2.91698877571e-16 + 0.797333791317*I)
            sage: vals[2] * vecs.column(2)  # random precision
            (-5.48020775807e-17 - 3.30603810196e-17*I, 1.09604155161e-16 + 6.61207620392e-17*I, -5.48020775807e-17 - 3.30603810196e-17*I)
        IMPLEMENTATION:
            Uses numpy.
        """
        if not self.is_square():
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self

        import numpy
        vals,vecs = numpy.linalg.eig(self.numpy())
        return ([sage.rings.complex_double.CDF(x) for x in vals],matrix(vecs))


    def eigenspaces(self, var='a'):
        r"""
        Return a list of pairs (e, V) where e runs through all complex
        eigenvalues of this matrix, and V is the corresponding
        eigenspace (always a 1-dimensional complex vector space).

        EXAMPLES:
            sage: m = matrix(CDF, 3, 3, [[0, 0, 1], [1, 0, 1], [0, 1, 0]]); m
            [  0   0 1.0]
            [1.0   0 1.0]
            [  0 1.0   0]
            sage: es = m.eigenspaces()
            sage: es # random
            [(1.32471795724 - 5.55111512313e-17*I, Vector space of degree 3 and dimension 1 over Complex Double Field
            User basis matrix:
            [0.413998885523 - 5.55111512313e-17*I 0.548431757932 - 2.22044604925e-16*I                       0.726517398056]), (-0.662358978622 + 0.562279512062*I, Vector space of degree 3 and dimension 1 over Complex Double Field
            User basis matrix:
            [                    0.655865618097 -0.434418480916 + 0.368779799722*I 0.0803836554683 - 0.488529222961*I]), (-0.662358978622 - 0.562279512062*I, Vector space of degree 3 and dimension 1 over Complex Double Field
            User basis matrix:
            [                    0.655865618097 -0.434418480916 - 0.368779799722*I 0.0803836554683 + 0.488529222961*I])]

            sage: e, v = es[0]
            sage: v = v.basis()[0]
            sage: a = v * m
            sage: b = e * v
            sage: diff = a.change_ring(CDF) - b
            sage: abs(abs(diff)) < 1e-10
            True
            sage: diff # random -- very small numbers
            (6.66133814775e-16 - 1.25526431318e-16*I, 3.24590553708e-16*I, 9.99200722163e-16 - 2.37225939001e-16*I)
        """
        e, v = self.left_eigenvectors()
        v = v.rows()
        pairs = []
        for l from 0<=l<len(e):
            c = v[l]
            pairs.append((e[l], c.parent().span_of_basis([c], check=False)))
        return pairs


    def __eigen_numpy(self):
        """
        Computes the eigenvalues and eigenvectors of this matrix acting
        from the left on column vectors.

        OUTPUT:
             eigenvalues -- as a list
             corresponding eigenvectors -- as a matrix whose ** ROWS ** are the eigenvectors of
                       self acting from the ** LEFT ** (which is really weird!)

        IMPLEMENTATION:
            Uses numpy.
        """
        if self._nrows != self._ncols:
            raise ValueError, "self must be square"
        if self._nrows == 0:
            return [], self

        import numpy
        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef Matrix_complex_double_dense _M, _result_matrix
        cdef int dims[2]
        cdef double *p
        cdef object temp
        cdef ndarray _n,_m
        parent = self.matrix_space(self._nrows,self._ncols)
        if self._nrows == 0 or self._ncols == 0:
            return [], parent.zero_matrix()
        _result_matrix = Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        _M=self
        dims[0] = _M._matrix.size1
        dims[1] = _M._matrix.size2
        temp = PyArray_FromDims(2, dims, 15)# 15 is a type code in ndarray.h here 15 is complex double
        _n = temp
        _n.flags = _n.flags&(~NPY_OWNDATA) # this performs a logical AND on NOT(NPY_OWNDATA), which sets that bit to 0
        _n.data = <char *> _M._matrix.data #numpy arrays store their data as char *
        v,_m = numpy.linalg.eig(_n)
        #It may be worthfile to check the flags that it is contiguous.
        p = <double *> _m.data
        for i from 0<=i<2*_M._nrows*_M._ncols:
            _result_matrix._matrix.data[i] = p[i]

        return ( [sage.rings.complex_double.CDF(x) for x in v],_result_matrix)   #todo: make the result a complex double matrix


    def solve_left(self,vec):
        """
        Solve the equation A*x = b, where

        EXAMPLES:
            sage: I = CDF.gen()
            sage: A =I*matrix(CDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1])
            sage: A   # slightly random output
            [1.0*I             2.0*I                5.0*I]
	    [7.59999990463*I   2.29999995232*I      1.0*I]
 	    [1.0*I             2.0*I               -1.0*I]
            sage: b = vector(CDF,[1,2,3])+I*vector(CDF,[1,2,3])
            sage: x = A.solve_left(b); x
             (-0.113695090439 + 0.113695090439*I, 1.39018087855 - 1.39018087855*I, -0.333333333333 + 0.333333333333*I)
            sage: A*x
            (1.0 + 1.0*I, 2.0 + 2.0*I, 3.0 + 3.0*I)


        TESTS:
        We test two degenerate cases:
            sage: A = matrix(CDF, 0, 3, [])
            sage: A.solve_left(vector(CDF,[]))
            ()
            sage: A = matrix(CDF, 3, 0, [])
            sage: A.solve_left(vector(CDF,3, [1,2,3]))
            (0, 0, 0)
        """
        cdef double *p
        cdef ComplexDoubleVectorSpaceElement _vec,ans
        cdef ndarray _result

        M = self._column_ambient_module()
        ans=M.zero_vector()
        if self._nrows == 0 or self._ncols == 0:
            return ans

        import numpy
        _vec=vec
        _result=numpy.linalg.solve(self.numpy(),_vec.numpy())
        p = <double *>_result.data

        memcpy(ans.v.data,_result.data,_result.dimensions[0]*sizeof(double)*2)
        return ans


    def solve_left_LU(self, vec):
        """
        Solve the equation A*x = b, where

        EXAMPLES:
            sage: I = CDF.gen()
            sage: A = I*matrix(CDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1])
            sage: A   # slightly random output
            [1.0*I             2.0*I                5.0*I]
	    [7.59999990463*I   2.29999995232*I      1.0*I]
 	    [1.0*I             2.0*I               -1.0*I]
            sage: b = vector(CDF,[1,2,3])+I*vector(CDF,[1,2,3])
            sage: x = A.solve_left(b); x
            (-0.113695090439 + 0.113695090439*I, 1.39018087855 - 1.39018087855*I, -0.333333333333 + 0.333333333333*I)
            sage: A*x
            (1.0 + 1.0*I, 2.0 + 2.0*I, 3.0 + 3.0*I)

        TESTS:
        We test two degenerate cases:
            sage: A = matrix(CDF, 0, 3, [])
            sage: A.solve_left_LU(vector(CDF,[]))
            ()
            sage: A = matrix(CDF, 3, 0, [])
            sage: A.solve_left_LU(vector(CDF,3, [1,2,3]))
            (0, 0, 0)

        This method precomputes and stores the LU decomposition before
        solving. If many equations of the form Ax=b need to be solved
        for a singe matrix A, then this method should be used instead
        of solve. The first time this method is called it will compute
        the LU decomposition.  If the matrix has not changed then
        subsequent calls will be very fast as the precomputed LU
        decomposition will be used.
        """
        if self._nrows == 0 or self._ncols == 0:
            M=self._column_ambient_module()
            return M.zero_vector()
        if self._nrows == 0 or self._ncols == 0:
            return self.row_module().zero_vector()
        import solve
        return solve.solve_matrix_complex_double_dense(self, vec)

    def determinant(self):
        """
        Return the determinant of self.

        ALGORITHM: Use GSL (LU decompositon)

        EXAMPLES:
            sage: m = matrix(CDF,2,range(4)); m.det()
            -2.0
            sage: m = matrix(CDF,0,[]); m.det()
            1.0
            sage: m = matrix(CDF, 2, range(6)); m.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be square
        """
        if self._nrows == 0 or self._ncols == 0:
            return sage.rings.complex_double.CDF(1)
        if not self.is_square():
            raise ValueError, "self must be square"
        cdef gsl_complex z
        if(self.fetch('LU_valid') !=True):
            self._c_compute_LU()
        z = gsl_linalg_complex_LU_det(self._LU, self._signum)
        return sage.rings.complex_double.CDF(GSL_REAL(z),GSL_IMAG(z))

    def log_determinant(self):
        """
        Compute the log of the absolute value of the determinant
        using GSL (LU decomposition).

        NOTE: This is useful if the usual determinant overlows.

        EXAMPLES:
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

        cdef double z
        if(self.fetch('LU_valid') !=True):
            self._c_compute_LU()
        z = gsl_linalg_complex_LU_lndet(self._LU)
        return sage.rings.real_double.RDF(z)

    def transpose(self):
        """
        Return the transpose of this matrix.

        EXAMPLES:
            sage: m = matrix(CDF,2,3,range(6)); m
            [  0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m.transpose()
            [  0 3.0]
            [1.0 4.0]
            [2.0 5.0]
            sage: m = matrix(CDF,0,3); m
            []
            sage: m.transpose().parent()
            Full MatrixSpace of 3 by 0 dense matrices over Complex Double Field
        """
        cdef Matrix_complex_double_dense trans
        cdef int result_copy
        parent  = self.matrix_space(self._ncols,self._nrows)
        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._ncols, self._nrows)
        trans = Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        result_copy = gsl_matrix_complex_transpose_memcpy(trans._matrix,self._matrix)
        if result_copy !=GSL_SUCCESS:
            raise ValueError, "Error copy matrix"
        return trans

    def LU(self):
        """
        Computes the LU decomposition of a matrix.

        For and square matrix A we can find matrices P,L, and U. s.t.
           P*A = L*U
        where P is a permutation matrix, L is lower triangular and U is upper triangular.

        OUTPUT:
            P, L, U -- as a tuple

        EXAMPLES:
            sage: m=matrix(CDF,4,range(16))
            sage: P,L,U = m.LU()
            sage: P*m
            [12.0 13.0 14.0 15.0]
            [   0  1.0  2.0  3.0]
            [ 8.0  9.0 10.0 11.0]
            [ 4.0  5.0  6.0  7.0]
            sage: L*U
            [12.0 13.0 14.0 15.0]
            [   0  1.0  2.0  3.0]
            [ 8.0  9.0 10.0 11.0]
            [ 4.0  5.0  6.0  7.0]
        """
        if not self.is_square():
            raise TypeError,"LU decomposition only works for square matrix"
        if self._ncols == 0:
            return self, self, self
        if self.fetch('LU_valid') != True:
            self._c_compute_LU()
        cdef Py_ssize_t i,j,k,l,copy_result
        cdef Matrix_complex_double_dense P, L,U
        cdef gsl_complex z
        parent = self.matrix_space(self._nrows,self._ncols)
        P=Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        L = Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        U = Matrix_complex_double_dense.__new__(Matrix_complex_double_dense,parent,None,None,None)
        for i from 0<=i<self._ncols:
            j = gsl_permutation_get(self._p,i)
            GSL_SET_COMPLEX(&z,1,0)
            gsl_matrix_complex_set(P._matrix,i,j,z)
            gsl_matrix_complex_set(L._matrix,i,i,z)
            z = gsl_matrix_complex_get(self._LU,i,i)
            gsl_matrix_complex_set(U._matrix,i,i,z)
            for l from 0<=l<i:
                z = gsl_matrix_complex_get(self._LU,i,l)
                gsl_matrix_complex_set(L._matrix,i,l,z)
                z = gsl_matrix_complex_get(self._LU,l,i)
                gsl_matrix_complex_set(U._matrix,l,i,z)

        return P, L, U

    def numpy(self):
        r"""
        This method returns a copy of the matrix as a numpy array. It
        is fast as the copy is done using the numpy C/api.

        EXAMPLES:
            sage: I = CDF.gen()
            sage: m = matrix(CDF,[[1,2],[3,4]])
            sage: m = I*m
            sage: n = m.numpy()
            sage: import numpy
            sage: numpy.linalg.eig(n)       # random low bits
            (array([  0.00000000e+00-0.37228132j,  -8.03393810e-17+5.37228132j]), array([[ 0.82456484 +0.00000000e+00j,  0.41597356 +8.77187320e-17j],
                   [-0.56576746 -0.00000000e+00j,  0.90937671 +0.00000000e+00j]]))

        TESTS:
            sage: m = matrix(CDF, 2, range(6)); m
            [  0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: m.numpy()
            array([[ 0.+0.j,  1.+0.j,  2.+0.j],
                   [ 3.+0.j,  4.+0.j,  5.+0.j]])
            sage: m = matrix(CDF,0,5,[]); m
            []
            sage: m.numpy()
            array([], shape=(0, 5), dtype=complex128)
            sage: m = matrix(CDF,5,0,[]); m
            []
        """
        if self._ncols == 0:
            import numpy
            return numpy.array([[]]*self._nrows, dtype='complex128')
        elif self._nrows == 0:
            import numpy
            return numpy.array([[]]*self._ncols, dtype='complex128').transpose()

        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef Matrix_complex_double_dense _M,_result_matrix
        _M=self
        cdef int dims[2]
        cdef double * data
        cdef int i
        cdef object temp
        cdef double *p
        cdef ndarray _n,_m
        dims[0] = _M._matrix.size1
        dims[1] = _M._matrix.size2
        data = <double*> malloc(sizeof(double)*dims[0]*dims[1]*2)
        memcpy(data,_M._matrix.data,sizeof(double)*dims[0]*dims[1]*2)
        temp = PyArray_FromDimsAndData(2, dims, 15,data)
        _n = temp
        _n.flags = _n.flags|(NPY_OWNDATA) # this sets the ownership flag
        return _n


    def _replace_self_with_numpy(self,numpy_matrix):
        """
        Replaces the entries of self with the elements of a numpy array with
        dtype = 'complex128'.

        EXAMPLES:
            sage: import numpy
            sage: j = numpy.complex(0,1)
            sage: a = numpy.array([[1*j, 2*j],[3*j, 4*j]], 'complex128')
            sage: m = matrix(CDF, 2, 2, 0)
            sage: m._replace_self_with_numpy(a)
            sage: m
            [1.0*I 2.0*I]
            [3.0*I 4.0*I]
        """
        if self._nrows == 0 or self._ncols == 0:
            return self

        cdef ndarray n
        cdef double *p
        n=numpy_matrix
        p=<double *>n.data
        memcpy(self._matrix.data,p,sizeof(double)*self._nrows*self._ncols*2)

    def _replace_self_with_numpy32(self,numpy_matrix):
        """
        Replaces the entries of self with the elements of a numpy array with
        dtype = 'complex64'.

        EXAMPLES:
            sage: import numpy
            sage: j = numpy.complex(0,1)
            sage: a = numpy.array([[1*j, 2*j],[3*j, 4*j]], 'complex64')
            sage: m = matrix(CDF, 2, 2, 0)
            sage: m._replace_self_with_numpy32(a)
            sage: m
            [1.0*I 2.0*I]
            [3.0*I 4.0*I]
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
        for i from 0 <= i < self._nrows*self._ncols*2:
            md[i] = <double> nd[i]

    cdef Vector _matrix_times_vector_(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._column_ambient_module().zero_vector()

        cdef ComplexDoubleVectorSpaceElement v_,ans
        cdef gsl_complex a,b
        cdef gsl_vector *vec
        cdef Py_ssize_t i,j
        M=self._column_ambient_module()
        v_ = v
        ans=M.zero_vector()
        if self._nrows == 0 or self._ncols == 0:
            return ans
        GSL_SET_COMPLEX(&a,1.0,0)
        GSL_SET_COMPLEX(&b,0,0)
        gsl_blas_zgemv(CblasNoTrans,a,self._matrix, v_.v,b,ans.v)
        return ans

    cdef Vector _vector_times_matrix_(self,Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self._row_ambient_module().zero_vector()

        cdef ComplexDoubleVectorSpaceElement v_,ans
        cdef gsl_complex a,b
        cdef gsl_vector *vec
        cdef Py_ssize_t i,j
        M=self._row_ambient_module()
        v_ = v
        ans=M.zero_vector()
        if self._nrows == 0 or self._ncols == 0:
            return ans
        GSL_SET_COMPLEX(&a,1.0,0)
        GSL_SET_COMPLEX(&b,0,0)
        gsl_blas_zgemv(CblasTrans,a,self._matrix, v_.v,b,ans.v)
        return ans

    def SVD(self):
        r"""
        Return the singular value decomposition of this matrix.

        INPUT:
            A -- a matrix
        OUTPUT:
            U, S, V -- matrices such that $A = U S V^t$, where
                       U and V are orthogonal and S is diagonal.


        EXAMPLES:
            sage: m = matrix(CDF,4,range(16))
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()    # slightly random output (due to computer architecture)
            [3.45569519412e-16               1.0               2.0               3.0]
            [4.0               5.0               6.0               7.0]
            [8.0               9.0              10.0              11.0]
            [12.0              13.0              14.0              15.0]

        A non-square example:
            sage: m = matrix(CDF, 2, range(6)); m
            [  0 1.0 2.0]
            [3.0 4.0 5.0]
            sage: U, S, V = m.SVD()
            sage: U
            [-0.274721127897 -0.961523947641]
            [-0.961523947641  0.274721127897]
            sage: S
            [7.34846922835             0             0]
            [            0           1.0             0]
            sage: V
            [-0.392540507864  0.824163383692  0.408248290464]
            [-0.560772154092  0.137360563949 -0.816496580928]
            [ -0.72900380032 -0.549442255795  0.408248290464]
            sage: U*S*V.transpose()           # random low bits
            [7.62194127257e-17               1.0               2.0]
            [              3.0               4.0               5.0]
            sage: m = matrix(CDF,3,2,range(6)); m
            [  0 1.0]
            [2.0 3.0]
            [4.0 5.0]
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()   # random low order bits
            [-8.13151629364e-19                1.0]
            [               2.0                3.0]
            [               4.0                5.0]

        TESTS:
            sage: m = matrix(CDF, 3, 0, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: m = matrix(CDF, 0, 3, []); m
            []
            sage: m.SVD()
            ([], [], [])
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


