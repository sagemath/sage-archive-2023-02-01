# distutils: depends = NTL/ZZ.h

from .types cimport zz_pContext_c

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    void zz_pContext_restore(zz_pContext_c *c)

