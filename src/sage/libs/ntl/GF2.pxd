from .types cimport GF2_c

cdef extern from "ccobject.h":
    void GF2_from_str "_from_str<GF2>"(GF2_c* dest, char* s)
    object GF2_to_PyString "_to_PyString<GF2>"(GF2_c *x)

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    int GF2_IsOne "IsOne"(GF2_c x)
    int GF2_IsZero "IsZero"(GF2_c x)

    void GF2_add "add"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_sub "sub"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_mul "mul"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_div "div"( GF2_c x, GF2_c a, GF2_c b)
    void GF2_negate "NTL::negate"(GF2_c x, GF2_c a)
    void GF2_power "NTL::power"(GF2_c t, GF2_c x, long e)
    long GF2_deg "deg"(GF2_c x)

    void GF2_conv_long "conv" (GF2_c x, long i)
    long GF2_conv_to_long "rep" (GF2_c x)
