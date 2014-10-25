from sage.libs.arb.arb cimport arb_t, arb_struct

cdef extern from "acb.h":
     ctypedef struct acb_struct:
         arb_struct real
         arb_struct imag
     ctypedef acb_struct[1] acb_t

     void acb_init(acb_t x)
     void acb_clear(acb_t x)
     void acb_set(acb_t z, const acb_t x)
     void acb_set_ui(acb_t z, unsigned long c)

     void acb_add(acb_t z, const acb_t x, const acb_t y, long prec)
     void acb_sub(acb_t z, const acb_t x, const acb_t y, long prec)
     void acb_mul(acb_t z, const acb_t x, const acb_t y, long prec)
     void acb_div(acb_t z, const acb_t x, const acb_t y, long prec)

     void acb_pow(acb_t z, const acb_t x, const acb_t y, long prec)
