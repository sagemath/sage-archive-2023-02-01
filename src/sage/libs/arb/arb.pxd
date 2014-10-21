cdef extern from "arb.h":
     ctypedef struct arb_struct:
         pass
     ctypedef arb_struct[1] arb_t

     void arb_init(arb_t x)
     void arb_clear(arb_t x)
     void arb_print(const arb_t x)
     void arb_set_ui(arb_t x, unsigned long y)
     void arb_zeta(arb_t z, const arb_t s, long prec)
