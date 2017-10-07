from .types cimport vec_GF2E_c

cdef extern from "ccobject.h":
    void vec_GF2E_from_str "_from_str<vec_GF2E>"(vec_GF2E_c* dest, char* s)
    object vec_GF2E_to_PyString "_to_PyString<vec_GF2E>"(vec_GF2E_c *x)
