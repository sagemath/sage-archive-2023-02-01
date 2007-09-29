include "decl.pxi"

from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

##cdef extern from "ntl_wrap.h":
##    struct zz_p_c

cdef class ntl_zz_p:
    cdef zz_p_c x
    cdef ntl_zz_pContext_class c
    cdef ntl_zz_p _new(ntl_zz_p self)
    cdef ntl_zz_p _add_c_impl(ntl_zz_p self, ntl_zz_p other)
    cdef ntl_zz_p _sub_c_impl(ntl_zz_p self, ntl_zz_p other)
    cdef ntl_zz_p _mul_c_impl(ntl_zz_p self, ntl_zz_p other)
    cdef ntl_zz_p _div_c_impl(ntl_zz_p self, ntl_zz_p other)

