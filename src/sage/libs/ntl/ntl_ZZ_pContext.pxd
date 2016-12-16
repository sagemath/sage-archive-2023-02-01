from .types cimport ZZ_pContext_c
from .ntl_ZZ cimport ntl_ZZ

cdef class ntl_ZZ_pContext_class(object):
    cdef ZZ_pContext_c x
    cdef void restore_c(self)
    cdef ntl_ZZ p
    cdef double p_bits
    cdef object __weakref__

cdef class ntl_ZZ_pContext_factory(object):
    cdef object context_dict
    cdef ntl_ZZ_pContext_class make_c(self, ntl_ZZ v)
