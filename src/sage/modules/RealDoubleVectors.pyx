#Modify free_module.py to check if base_ring is RealDoubleField and if so call this constructor#

cimport free_module_element
import  free_module_element

cimport ComplexDoubleVectors
import  ComplexDoubleVectors


import sage.rings.complex_double





cdef class RealDoubleVectorSpace_element(free_module_element.FreeModuleElement):

    def __init__(self,parent,x,coerce_entries=True,copy=True):
        free_module_element.FreeModuleElement.__init__(self,sage.rings.real_double.RDF)
        self._parent=parent
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
        gsl_vector_free(self.v)

    def __len__(self):
        return self.v.size

    def __setitem__(self,size_t i,x):
        if i < 0 or i >=self.v.size:
            raise IndexError
        else:
            gsl_vector_set(self.v,i,x)

 #       self.vec[i] = x

    def __getitem__(self,size_t i):
        if i < 0 or i >=self.v.size:
            raise IndexError
        else:
            return gsl_vector_get(self.v,i)
#        return self.vec[i]


#    def __add__(left,right):
#       return (<RealDoubleVectorSpace_element>left)._add_(<RealDoubleVectorSpace_element>right)

    cdef _add_(RealDoubleVectorSpace_element self, RealDoubleVectorSpace_element right):
        cdef int i
        cdef int n
        cdef int err_1
        cdef int err_2
        n = self.v.size
        cdef RealDoubleVectorSpace_element result
        if self.parent() is right.parent():
            gsl_set_error_handler_off()
            result = <RealDoubleVectorSpace_element> RealDoubleVectorSpace_element.__new__(RealDoubleVectorSpace_element)
            (<RealDoubleVectorSpace_element> result).v = <gsl_vector*> gsl_vector_alloc(n)
            if (<RealDoubleVectorSpace_element> result).v is not NULL:
                err_1 = gsl_vector_memcpy((<RealDoubleVectorSpace_element>result).v,self.v)
                err_2 = gsl_vector_add( (<RealDoubleVectorSpace_element>result).v,(<RealDoubleVectorSpace_element> right).v)

                (<RealDoubleVectorSpace_element>result)._parent= (<RealDoubleVectorSpace_element>right)._parent
                status = (not err_1 and not err_2)
                if status:
                    return result
                else:
                    raise TypeError, "Memory allocation error"


    def fft(self):
        cdef ComplexDoubleVectors.ComplexDoubleVectorSpace_element result
        cdef int i
        P = self.parent().change_ring( sage.rings.complex_double.CDF )
        result = ComplexDoubleVectors.ComplexDoubleVectorSpace_element.__new__(ComplexDoubleVectors.ComplexDoubleVectorSpace_element,
                                                                               P, None)
        for i from 0 <=i< self.v.size:
            (result.v).data[2*i]= self.v.data[i]
        result.fft(inplace=True)
        return result


    def test_add(self,right):
        return self._add_(right)
#    def __getslice__(


#    def free_module(self):
#        return self


