#Modify free_module.py to check if base_ring is RealDoubleField and if so call this constructor#

#include '../gsl/gsl.pxi'
#ctypedef int size_t

cimport sage.modules.free_module_element
import sage.modules.free_module_element

import sage.modules.free_module

#cdef extern from "Python.h":
#    int PyObject_TypeCheck(object o,object t)



class RealDoubleVectorSpace_class(sage.modules.free_module.FreeModule_ambient_field):
    def __init__(self,n):
        sage.modules.free_module.FreeModule_ambient_field.__init__(self,sage.rings.real_double.RealDoubleField_class )


    def __call__(self,x,coerce_entries=True,copy=True,check_elements=True):

        return RealDoubleVectorSpace_element(n,x)







cdef class RealDoubleVectorSpace_element(sage.modules.free_module_element.FreeModuleElement):

    def __init__(self,n,x,stride=1):
        self.vec = <double*> malloc(sizeof(double)*n)
        self.n = n
        self.stride =1
        cdef int i
        for i from 0 <=i < n:
            self.vec[i]=x[i]


    def __dealloc__(self):
        PyMem_Free(self.vec)

    def __len__(self):
        return self.n

    def __setitem__(self,size_t i,x):
        if i < 0 or i >=self.n:
            raise IndexError
        self.vec[i] = x

    def __getitem__(self,size_t i):
        if i < 0 or i >=self.n:
            raise IndexError
        self.vec[i]
#    def __getslice__(
def free_module():
    pass

    def free_module(self):
        return self

