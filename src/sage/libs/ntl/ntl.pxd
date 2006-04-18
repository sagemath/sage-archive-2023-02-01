cdef extern from "wrap.h":
    struct ZZ

cdef class ntl_ZZ:
    cdef ZZ* x
    cdef set(self, void *y)

cdef extern from "wrap.h":
    struct ZZX

cdef class ntl_ZZX:
    cdef ZZX* x
    cdef set(self, void *y)

cdef extern from "wrap.h":
    struct ZZ_p

cdef class ntl_ZZ_p:
    cdef ZZ_p* x
    cdef set(self, void *y)

cdef extern from "wrap.h":
    struct ZZ_pX

cdef class ntl_ZZ_pX:
    cdef ZZ_pX* x
    cdef set(self, void *y)

cdef extern from "wrap.h":
    struct mat_ZZ

cdef class ntl_mat_ZZ:
    cdef mat_ZZ* x
    cdef long __nrows, __ncols

cdef extern from "wrap.h":
    struct GF2X

cdef class ntl_GF2X:
    cdef GF2X* gf2x_x
    cdef set(self,void *y)

cdef extern from "wrap.h":
    struct GF2E

cdef class ntl_GF2E(ntl_GF2X):
    cdef GF2E* gf2e_x
    cdef set(self,void *y)

cdef extern from "wrap.h":
    struct mat_GF2E

cdef class ntl_mat_GF2E:
    cdef mat_GF2E* x
    cdef long __nrows, __ncols

cdef extern from "wrap.h":
    struct GF2EX

cdef class ntl_GF2EX:
    cdef GF2EX* x
    cdef set(self,void *y)
