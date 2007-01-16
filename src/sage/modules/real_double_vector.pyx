"""
Real double vectors

AUTHOR:
    -- Josh Kantor (2006-10)
    -- William Stein (2006) rewrite
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

cdef class RealDoubleVectorSpaceElement(free_module_element.FreeModuleElement):
    cdef _new_c(self, gsl_vector* v):
        cdef RealDoubleVectorSpaceElement y
        y = PY_NEW(RealDoubleVectorSpaceElement)
        y._parent = self._parent
        y._degree = self._degree
        y.v = v
        return y

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

    def fft(self):
        cdef complex_double_vector.ComplexDoubleVectorSpaceElement result
        cdef int i
        P = self.parent().change_ring( sage.rings.complex_double.CDF )
        result = complex_double_vector.ComplexDoubleVectorSpaceElement.__new__(
            complex_double_vector.ComplexDoubleVectorSpaceElement, P, None)
        for i from 0 <=i< self.v.size:
            (result.v).data[2*i]= self.v.data[i]
        result.fft(inplace=True)
        return result



