cimport free_module_element
import  free_module_element

from sage.structure.element cimport Element, ModuleElement, RingElement

from sage.rings.complex_double cimport ComplexDoubleElement
from sage.rings.complex_double import CDF


cdef class ComplexDoubleVectorSpace_element(free_module_element.FreeModuleElement):
    cdef _new_c(self, gsl_vector_complex* v):
        cdef ComplexDoubleVectorSpace_element y
        y = PY_NEW(ComplexDoubleVectorSpace_element)
        y._parent = self._parent
        y.v = v
        return y

    cdef gsl_vector_complex* gsl_vector_complex_copy(self):
        """
        Return a copy of the underlying GSL vector of self.
        """
        cdef gsl_vector_complex* v
        v = <gsl_vector_complex*> gsl_vector_complex_alloc(self.v.size)
        if v is NULL:
            raise MemoryError, "error allocating real double vector"
        if gsl_blas_zcopy(self.v , (<ComplexDoubleVectorSpace_element> result).v):
            raise RuntimeError, "error copying real double vector"
        return v

    def __new__(self, parent,x,coerce=True, copy=True):
        self.v = gsl_vector_complex_calloc(parent.rank())
        self._parent = parent

    def __init__(self,parent,x,coerce = True, copy = True):
        #Think about whether coerce = False or copy = False should be behave differently
        free_module_element.FreeModuleElement.__init__(self,sage.rings.complex_double.CDF)
        cdef int i
        cdef int n
        cdef gsl_complex z_temp
        n = parent.rank()
        cdef int length



        try:
            length=len(x)
        except TypeError:
            if x == 0:
                gsl_set_error_handler_off()#This keeps failed gsl routines from calling abort
                if self.v is not NULL:
                    return
                else:
                    raise MemoryError, "error allocating memory"
            else:
                self.v = NULL
                raise TypeError, "must be a list, tuple, vector or 0"

        gsl_set_error_handler_off()#This keeps failed gsl routines from calling abort
        self.v = <gsl_vector_complex *> gsl_vector_complex_calloc(n)
    #    if x == 0 and self.v is not NULL:
    #        return

        if self.v is not NULL and length == n:
            for i from 0 <=i <n:
                z = sage.rings.complex_double.CDF(x[i])
                z_temp = <gsl_complex> (<ComplexDoubleElement> z)._complex
                gsl_vector_complex_set(self.v,i, z_temp)


        elif self.v is NULL:
            raise MemoryError,"Error allocating memory"

        elif length != n:
            raise TypeError, "User supplied vector not same length as rank of vector space"


    def __dealloc__(self):
        if not (self.v is NULL):
            gsl_vector_complex_free(self.v)


    def __len__(self):
        return self.v.size

    def __setitem__(self,size_t i,x):
        cdef gsl_complex z_temp

        if i < 0 or i >= self.v.size:
            raise IndexError
        else:
            z_temp =<gsl_complex> (<ComplexDoubleElement> CDF(x))._complex
            gsl_vector_complex_set(self.v,i,z_temp)

    def __getitem__(self,size_t i):
        cdef gsl_complex z_temp
        if i < 0 or i >= self.v.size:
            raise IndexError
        else:
            z_temp = <gsl_complex> gsl_vector_complex_get(self.v,i)
            return sage.rings.complex_double.CDF(GSL_REAL(z_temp),GSL_IMAG (z_temp))


    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef gsl_vector_complex* v
        gsl_set_error_handler_off()
        v = self.gsl_vector_complex_copy()
        _sig_on
        if gsl_blas_zaxpy(gsl_complex_rect(1,0), (<ComplexDoubleVectorSpace_element> right).v, v):
            _sig_off
            raise RuntimeError
        _sig_off
        return self._new_c(v)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        print "todo: compelx _sub_c_impl may be backwards!"
        gsl_set_error_handler_off()
        cdef gsl_vector_complex* v
        v = self.gsl_vector_complex_copy()
        _sig_on
        if gsl_blas_zaxpy(gsl_complex_rect(-1,0), (<ComplexDoubleVectorSpace_element> right).v, v):
            _sig_off
            raise RuntimeError
        _sig_off
        return self._new_c(v)

    def fft(self,direction = "forward",algorithm = "radix2",inplace=False):
        gsl_set_error_handler_off()
        cdef gsl_fft_complex_wavetable* wavetable
        cdef gsl_fft_complex_workspace* workspace
        cdef int err,err_1
        cdef ComplexDoubleVectorSpace_element result
        cdef double* data
        if inplace:
            data = self.v.data
        else:
            # TODO -- this all needs to be fixed; it's broken.
            raise NotImplementedError
            result = <ComplexDoubleVectorSpace_element>ComplexDoubleVectorSpace_element.__new__(ComplexDoubleVectorSpace_element)
            (<ComplexDoubleVectorSpace_element>result).v=<gsl_vector_complex *> gsl_vector_complex_alloc(self.v.size)# add error check for mem alloc
            data = (<ComplexDoubleVectorSpace_element>result).v.data
            gsl_blas_zcopy(self.v,(<ComplexDoubleVectorSpace_element>result).v)
            (<ComplexDoubleVectorSpace_element> result)._parent = self._parent


        if direction == "forward":
            if ispow(self.v.size):
                err = gsl_fft_complex_radix2_forward(data,self.v.stride,self.v.size)
            else:
                wavetable = gsl_fft_complex_wavetable_alloc(self.v.size)
                workspace = gsl_fft_complex_workspace_alloc(self.v.size)
                if wavetable !=NULL and workspace!= NULL:
                    err = gsl_fft_complex_forward(data,1,self.v.size,wavetable,workspace)

                gsl_fft_complex_wavetable_free(wavetable)
                gsl_fft_complex_workspace_free(workspace)
            if not err:
                pass
            else:
                raise ValueError, "error with transform"



        elif direction == "backward":
            if ispow(self.v.size):
                err  = gsl_fft_complex_radix2_inverse(data,self.v.stride,self.v.size)
            else:
                wavetable = gsl_fft_complex_wavetable_alloc(self.v.size)
                workspace = gsl_fft_complex_workspace_alloc(self.v.size)
                if wavetable !=NULL and workspace!= NULL:
                    err = gsl_fft_complex_inverse(data,1,self.v.size,wavetable,workspace)

                gsl_fft_complex_wavetable_free(wavetable)
                gsl_fft_complex_workspace_free(workspace)
            if not err:
                pass
            else:
                raise ValueError, "error with transform"

        if inplace == True:
            return
        else:
            return result



cdef int ispow(int n):
    while n and n%2==0:
        n = n>>1
    return n == 1

