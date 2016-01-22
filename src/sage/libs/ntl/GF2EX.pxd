from .types cimport GF2EX_c

cdef extern from "ccobject.h":
    void GF2EX_from_str "_from_str<GF2EX>"(GF2EX_c* dest, char* s)
    object GF2EX_to_PyString "_to_PyString<GF2EX>"(GF2EX_c *x)

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void GF2EX_add "add"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_sub "sub"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_mul "mul"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_negate "NTL::negate"(GF2EX_c x, GF2EX_c a)
    void GF2EX_power "NTL::power"(GF2EX_c t, GF2EX_c x, long e)
    int GF2EX_IsOne "IsOne"(GF2EX_c x)
    int GF2EX_IsZero "IsZero"(GF2EX_c x)
