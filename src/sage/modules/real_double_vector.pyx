"""
Real double vectors

AUTHOR:
    -- Josh Kantor (2006-10)
    -- William Stein (2006) rewrite

TESTS:
    sage: v = vector(RDF, [1,2,3,4])
    sage: loads(dumps(v)) == v
"""
###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################


from sage.structure.element cimport ModuleElement, RingElement, Vector

from sage.rings.real_double cimport RealDoubleElement
from sage.rings.real_double import new_RealDoubleElement

cimport free_module_element
import  free_module_element

cimport complex_double_vector
import  complex_double_vector

import sage.rings.complex_double

include '../ext/stdsage.pxi'

cdef extern from "arrayobject.h":
#The following exposes the internal C structure of the numpy python object
# extern class [object PyArrayObject]  tells pyrex that this is
# a compiled python class defined by the C struct PyArrayObject
    cdef enum:
        NPY_OWNDATA = 0x0004 #bit mask so numpy does not free array contents when its destroyed

    ctypedef int intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

    object PyArray_FromDims(int,int *,int)
    object PyArray_FromDimsAndData(int,int*,int,double *)
    void import_array()




cdef class RealDoubleVectorSpaceElement(free_module_element.FreeModuleElement):
    cdef _new_c(self, gsl_vector* v):
        cdef RealDoubleVectorSpaceElement y
        y = PY_NEW(RealDoubleVectorSpaceElement)
        y._parent = self._parent
        y._degree = self._degree
        y.v = v
        return y

    cdef int is_dense_c(self):
        return True
    cdef int is_sparse_c(self):
        return False

    cdef gsl_vector* gsl_vector_copy(self) except NULL:
        """
        Return a copy of the underlying GSL vector of self.
        """
        cdef gsl_vector* v
        v = <gsl_vector*> gsl_vector_alloc(self.v.size)
        if v is NULL:
            raise MemoryError, "error allocating real double vector"
        if gsl_vector_memcpy(v,self.v):
            raise RuntimeError, "error copying real double vector"
        return v

    def __copy__(self, copy=True):
        return self._new_c(self.gsl_vector_copy())

    def __new__(self, parent=None,x=None,coerce=True,copy=True):
        self.v = NULL

    def __init__(self,parent,x,coerce=True,copy=True):
        self._parent = parent
        self._degree = parent.degree()
        cdef int n
        cdef int i
        n = parent.rank()

        try:
            length = len(x)
        except TypeError:
            if x == 0:
                gsl_set_error_handler_off()
                self.v = <gsl_vector*> gsl_vector_calloc(n)
                if self.v is not NULL:
                    return
                else:
                    self.v = NULL
                    raise MemoryError, "Error allocating memory"


        gsl_set_error_handler_off()
        self.v = <gsl_vector*> gsl_vector_calloc(n)

        if self.v is not NULL and length ==n:
            for i from 0 <=i < n:
                gsl_vector_set(self.v, i, x[i] )

        elif self.v is NULL:
            raise MemoryError, "Error allocating memory for vector"

        elif length !=n:
            raise TypeError,"user supplied bector must be same length as rank of vector space"

    def __dealloc__(self):
        if self.v:
            gsl_vector_free(self.v)

    def __len__(self):
        return self.v.size

    def __setitem__(self,size_t i,x):
        if i < 0 or i >=self.v.size:
            raise IndexError
        else:
            gsl_vector_set(self.v,i,x)

    def __getitem__(self,size_t i):
        """
        Return the ith entry of self.

        EXAMPLES:
            sage: v = vector(RDF, [1, sqrt(2), -1]); v
            (1.0, 1.41421356237, -1.0)
            sage: a = v[1]; a
            1.41421356237
            sage: parent(a)
            Real Double Field
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        cdef RealDoubleElement x
        if i < 0 or i >=self.v.size:
            raise IndexError, 'index out of range'
        else:
            x = new_RealDoubleElement()
            x._value = gsl_vector_get(self.v,i)
            return x

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_add(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error adding real double vectors"
        return self._new_c(v)


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_sub(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error subtracting real double vectors"
        return self._new_c(v)

    cdef Vector _vector_times_vector_c_impl(self, Vector right):
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_mul(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error multiplying real double vectors"
        return self._new_c(v)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        cdef gsl_vector* v
        v = self.gsl_vector_copy()
        gsl_vector_scale(v, (<RealDoubleElement>left)._value)
        return self._new_c(v)

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        cdef gsl_vector* v
        v = self.gsl_vector_copy()
        gsl_vector_scale(v, (<RealDoubleElement>right)._value)
        return self._new_c(v)

    def change_ring(self, R):
        """
        EXAMPLES:

        """
        if isinstance(R, sage.rings.complex_double.ComplexDoubleField_class):
            return self.complex_vector()
        else:
            return free_module_element.FreeModuleElement.change_ring(self, R)

    def complex_vector(self):
        """
        Return the associated complex vector, i.e., this vector but with
        coefficients viewed as complex numbers.
        """
        cdef complex_double_vector.ComplexDoubleVectorSpaceElement result
        cdef gsl_vector_complex* v
        v = <gsl_vector_complex*> gsl_vector_complex_calloc(self.v.size)
        cdef Py_ssize_t i
        P = self.parent().change_ring( sage.rings.complex_double.CDF )
        result = complex_double_vector.ComplexDoubleVectorSpaceElement.__new__(
            complex_double_vector.ComplexDoubleVectorSpaceElement, None, None)
        result._parent = P
        result._degree = self._degree
        result.v = v
        for i from 0 <=i< self.v.size:
            (result.v).data[2*i]= self.v.data[i]
        return result

    def fft(self):
        v = self.complex_vector()
        v.fft(inplace=True)
        return v

    def numpy(self):
        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef RealDoubleVectorSpaceElement _V,_result_vector
        _V=self
        cdef int dims[1]
        cdef double * data
        cdef int i
        cdef object temp
        cdef double *p
        cdef ndarray _n,_m
        dims[0] = _V.v.size
#        dims[1] = _M._matrix.size2
        data = <double*> malloc(sizeof(double)*dims[0])
        memcpy(data,_V.v.data,sizeof(double)*dims[0])
        temp = PyArray_FromDimsAndData(1, dims, 12,data)
        _n = temp
        _n.flags = _n.flags|(NPY_OWNDATA) # this sets the ownership flag
        return _n

    def _replace_self_with_numpy(self,numpy_array):
        cdef double *p
        cdef ndarray n
        n=numpy_array
        p=<double *>n.data
        memcpy(self.v.data,p,self.v.size*sizeof(double))

