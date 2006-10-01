import sage.ext.sage_object
cimport sage.ext.sage_object

cdef class Generators(sage.ext.sage_object.SageObject):
    cdef public object __gens
    cdef public object __gens_dict
    cdef public object __list
    cdef public object __names
    cdef public object __latex_names

cdef class MultiplicativeAbelianGenerators(Generators):
    cdef public object __generator_orders

cdef class AdditiveAbelianGenerators(Generators):
    cdef public object __generator_orders
