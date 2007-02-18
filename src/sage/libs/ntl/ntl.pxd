cdef extern from "gmp.h":
    struct mpz_t

cdef extern from "ntl_wrap.h":
    struct ZZ
    cdef int ZZ_to_int(ZZ* x)
    cdef ZZ* int_to_ZZ(int value)
    cdef void ZZ_set_from_int(ZZ* x, int value)

cdef class ntl_ZZ:
    cdef ZZ* x
    cdef set(self, void *y)
    cdef public int get_as_int(ntl_ZZ self)
    cdef public void set_from_int(ntl_ZZ self, int value)

cdef extern from "ntl_wrap.h":
    struct ZZX
    cdef void ZZX_setitem_from_int(ZZX* x, long i, int value)
    cdef int ZZX_getitem_as_int(ZZX* x, long i)
    cdef void ZZX_getitem_as_mpz(mpz_t* output, ZZX* x, long i)

cdef class ntl_ZZX:
    cdef ZZX* x
    cdef set(self, void *y)
    cdef void setitem_from_int(ntl_ZZX self, long i, int value)
    cdef int getitem_as_int(ntl_ZZX self, long i)

cdef extern from "ntl_wrap.h":
    struct ZZ_p
    cdef int ZZ_p_to_int(ZZ_p* x)
    cdef ZZ_p* int_to_ZZ_p(int value)
    cdef void ZZ_p_set_from_int(ZZ_p* x, int value)

cdef class ntl_ZZ_p:
    cdef ZZ_p* x
    cdef set(self, void *y)
    cdef public int get_as_int(ntl_ZZ_p self)
    cdef public void set_from_int(ntl_ZZ_p self, int value)

cdef extern from "ntl_wrap.h":
    struct ZZ_pX
    cdef void ZZ_pX_setitem_from_int(ZZ_pX* x, long i, int value)
    cdef int ZZ_pX_getitem_as_int(ZZ_pX* x, long i)

cdef class ntl_ZZ_pX:
    cdef ZZ_pX* x
    cdef set(self, void *y)
    cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value)
    cdef int getitem_as_int(ntl_ZZ_pX self, long i)

cdef extern from "ntl_wrap.h":
    struct mat_ZZ

cdef class ntl_mat_ZZ:
    cdef mat_ZZ* x
    cdef long __nrows, __ncols

cdef extern from "ntl_wrap.h":
    struct GF2X

cdef class ntl_GF2X:
    cdef GF2X* gf2x_x
    cdef set(self,void *y)

cdef extern from "ntl_wrap.h":
    struct GF2E

cdef class ntl_GF2E(ntl_GF2X):
    cdef GF2E* gf2e_x
    cdef set(self,void *y)

cdef extern from "ntl_wrap.h":
    struct mat_GF2E

cdef class ntl_mat_GF2E:
    cdef mat_GF2E* x
    cdef long __nrows, __ncols

cdef extern from "ntl_wrap.h":
    struct GF2EX

cdef class ntl_GF2EX:
    cdef GF2EX* x
    cdef set(self,void *y)
