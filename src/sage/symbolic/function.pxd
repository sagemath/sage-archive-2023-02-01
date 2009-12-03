from sage.structure.sage_object cimport SageObject

cdef class Function(SageObject):
    cdef unsigned int _serial
    cdef int _nargs
    cdef object _name
    cdef object _latex_name
    cdef object _conversions
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class GinacFunction(Function):
    cdef object _ginac_name
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class CustomizableFunction(Function):
    cdef _register_function(self)

cdef class BuiltinFunction(CustomizableFunction):
    cdef _is_registered(self)

cdef class SymbolicFunction(CustomizableFunction):
    # cache hash value
    cdef long _hash_(self)
    cdef bint __hinit
    cdef long __hcache
    cdef _is_registered(self)
