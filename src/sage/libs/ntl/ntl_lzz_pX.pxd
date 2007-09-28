include "decl.pxi"

from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

cdef extern from "ntl_wrap.h":
    struct zz_pX

cdef class ntl_zz_pX:
    cdef zz_pX_c x
    cdef ntl_zz_pContext_class c
    cdef ntl_zz_pX _new(self)
    cdef ntl_zz_pX _add_c_impl(ntl_zz_pX self, ntl_zz_pX other)
    cdef ntl_zz_pX _sub_c_impl(ntl_zz_pX self, ntl_zz_pX other)
    cdef ntl_zz_pX _mul_c_impl(ntl_zz_pX self, ntl_zz_pX other)
    cdef ntl_zz_pX _div_c_impl(ntl_zz_pX self, ntl_zz_pX other)
    cdef ntl_zz_pX _mod_c_impl(ntl_zz_pX self, ntl_zz_pX other)

