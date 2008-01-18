"""
Real double vectors

AUTHOR:
    -- Josh Kantor (2006-10)
    -- William Stein (2006) rewrite

TESTS:
    sage: v = vector(RDF, [1,2,3,4])
    sage: loads(dumps(v)) == v
    True
"""
###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################


from sage.structure.element cimport ModuleElement, RingElement, Vector, Element

cimport free_module_element
import  free_module_element

cimport complex_double_vector
import  complex_double_vector

from sage.rings.real_double cimport RealDoubleElement

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
        y._is_mutable = 1
        y._parent = self._parent
        y._degree = self._degree
        y.v = v
        return y

    cdef bint is_dense_c(self):
        return 1
    cdef bint is_sparse_c(self):
        return 0

    cdef gsl_vector* gsl_vector_copy(self) except NULL:
        """
        Return a copy of the underlying GSL vector of self.
        """
        if not self.v:
            return NULL
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

    def __reduce__(self):
        return (unpickle_v0, (self._parent, self.list(), self._degree))

    def __init__(self, parent, x, coerce=True, copy=True):
        self._parent = parent
        self._is_mutable = 1

        cdef int i
        cdef int n = parent.degree()
        self._degree = n
        if n == 0:
            self.v = NULL
            return

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

        if self.v is not NULL and length == n:
            for i from 0 <= i < n:
                gsl_vector_set(self.v, i, x[i])

        elif self.v is NULL:
            raise MemoryError, "Error allocating memory for vector"

        elif length != n:
            raise TypeError,"user supplied vector must be same length as rank of vector space"

    def __dealloc__(self):
        if self.v:
            gsl_vector_free(self.v)

    def __len__(self):
        if not self.v:
            return 0
        return self.v.size

    def __setitem__(self,size_t i,x):
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        if not self.v or i < 0 or i >=self.v.size:
            raise IndexError, "index out of bounds"
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
        if not self.v or i < 0 or i >=self.v.size:
            raise IndexError, 'index out of range'
        else:
            x = <RealDoubleElement>PY_NEW(RealDoubleElement)
            x._value = gsl_vector_get(self.v,i)
            return x

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        if not self.v:
            return self
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_add(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error adding real double vectors"
        return self._new_c(v)


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        if not self.v:
            return self
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_sub(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error subtracting real double vectors"
        return self._new_c(v)

    cdef Element _dot_product_c_impl(self, Vector right):
        """
        Dot product of self and right.

        EXAMPLES:
            sage: v = vector(RDF, [1,2,3]); w = vector(RDF, [2, 4, -3])
            sage: v*w
            1.0
            sage: w*v
            1.0
        """
        cdef RealDoubleElement y = <RealDoubleElement>PY_NEW(RealDoubleElement)
        y._value = 0
        cdef unsigned int i
        for i from 0 <= i < self.v.size:
            y._value += gsl_vector_get(self.v, i) * gsl_vector_get((<RealDoubleVectorSpaceElement>right).v, i)
        return y

    cdef Vector _pairwise_product_c_impl(self, Vector right):
        """
        Return the component-wise product of self and right.

        EXAMPLES:
            sage: v = vector(RDF, [1,2,3]); w = vector(RDF, [2, 4, -3])
            sage: v.pairwise_product(w)
            (2.0, 8.0, -9.0)
        """
        if not right.parent() == self.parent():
            right = self.parent().ambient_module()(right)
        if not self.v:
            return self
        cdef gsl_vector* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_copy()
        if gsl_vector_mul(v, (<RealDoubleVectorSpaceElement> right).v):
            raise RuntimeError, "error multiplying real double vectors"
        return self._new_c(v)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        if not self.v:
            return self
        cdef gsl_vector* v
        v = self.gsl_vector_copy()
        gsl_vector_scale(v, (<RealDoubleElement>left)._value)
        return self._new_c(v)

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        if not self.v:
            return self
        cdef gsl_vector* v
        v = self.gsl_vector_copy()
        gsl_vector_scale(v, (<RealDoubleElement>right)._value)
        return self._new_c(v)

    def change_ring(self, R):
        """
        EXAMPLES:
            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.change_ring(CC)
            (0, 1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: v.change_ring(CDF)
            (0, 1.0, 2.0, 3.0)
            sage: v.change_ring(RR)
            (0.000000000000000, 1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: v = vector(RDF,0)
            sage: v.change_ring(CC)
            ()
        """
        if isinstance(R, sage.rings.complex_double.ComplexDoubleField_class):
            return self.complex_vector()
        else:
            return free_module_element.FreeModuleElement.change_ring(self, R)

    def complex_vector(self):
        """
        Return the associated complex vector, i.e., this vector but with
        coefficients viewed as complex numbers.

        EXAMPLES:
            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.complex_vector()
            (0, 1.0, 2.0, 3.0)
            sage: v = vector(RDF,0)
            sage: v.complex_vector()
            ()
        """
        if not self.v:
            return self.parent().change_ring(sage.rings.complex_double.CDF).zero_vector()

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

    def fft(self, direction = "forward"):
        """
        Return the fast Fourier transform of this vector over the complex numbers.

        INPUT:
           direction -- 'forward' (default) or 'backward'

        EXAMPLES:
            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.fft()
            (6.0, -2.0 + 2.0*I, -2.0, -2.0 - 2.0*I)
            sage: v.fft(direction='backward')
            (1.5, -0.5 - 0.5*I, -0.5, -0.5 + 0.5*I)
            sage: v.fft(direction='backward').fft()     # random low order bits
            (0, 1.0 - 5.74627151417e-18*I, 2.0, 3.0 + 5.74627151417e-18*I)
        """
        if not self.v:
            return self
        v = self.complex_vector()
        v.fft(direction = direction, inplace=True)
        return v

    def numpy(self):
        """
        Return numpy array corresponding to this vector.

        EXAMPLES:
            sage: v = vector(RDF,4,range(4))
            sage: v.numpy()
            array([ 0.,  1.,  2.,  3.])
            sage: v = vector(RDF,0)
            sage: v.numpy()
            array([], shape=(1, 0), dtype=float64)
        """
        if not self.v:
            import numpy
            return numpy.array([[]])

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
        if not self.v:
            return
        cdef double *p
        cdef ndarray n
        n=numpy_array
        p=<double *>n.data
        memcpy(self.v.data,p,self.v.size*sizeof(double))

    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(RDF, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return free_module_element.vector( [e.n(*args, **kwargs) for e in self] )

    #############################
    # statistics
    #############################
    def mean(self):
        return gsl_stats_mean(self.v.data, self.v.stride, self.v.size)

    def variance(self):
        return gsl_stats_variance(self.v.data, self.v.stride, self.v.size)

    #def covariance(self):
    #    return gsl_stats_covariance(self.v.data, self.v.stride, self.v.size)

    def standard_deviation(self):
        """
        EXAMPLES:
        sage: v = vector(RDF, 5, [1,2,3,4,5])
        sage: v.standard_deviation()
        1.5811388300841898
        """
        return gsl_stats_sd(self.v.data, self.v.stride, self.v.size)

    def stats_skew(self):
        return gsl_stats_skew(self.v.data, self.v.stride, self.v.size)

    def stats_kurtosis(self):
        return gsl_stats_kurtosis(self.v.data, self.v.stride, self.v.size)

    def stats_lag1_autocorrelation(self):
        return gsl_stats_lag1_autocorrelation(self.v.data, self.v.stride, self.v.size)



def unpickle_v0(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    return parent(entries)
