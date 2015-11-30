# distutils: depends = NTL/ZZ.h

from .types cimport ZZ_pEContext_c

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void ZZ_pEContext_restore(ZZ_pEContext_c *c)
