from .types cimport vec_GF2_c

cdef extern from "ccobject.h":
    void vec_GF2_from_str "_from_str<vec_GF2>"(vec_GF2_c* dest, char* s)
    object vec_GF2_to_PyString "_to_PyString<vec_GF2>"(vec_GF2_c *x)
