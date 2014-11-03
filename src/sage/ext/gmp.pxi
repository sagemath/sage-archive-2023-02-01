# gmp.inc -- misc. useful GMP functions that don't depend on
#            other things and can be included in other files
# USAGE
#
# Include the following at the top of client .pyx file.
#
#    include "gmp.pxi"
#
# to include this in a file.

include 'sage/ext/interrupt.pxi'

cimport libc.stdlib

############ The following is the "one global set of vars"
cdef extern from "gmp_globals.h":
    cdef mpz_t u, v, q, u0, u1, u2, v0, v1, v2, t0, t1, t2, x, y, ssqr, m2
    # changed sqr to ssqr due to a collision with ntl
    cdef mpq_t tmp

    cdef mpz_t a1, a2, mod1, sage_mod2, g, s, t, xx

    cdef mpz_t crtrr_a, crtrr_mod

    cdef mpz_t rand_val, rand_n, rand_n1

    cdef gmp_randstate_t rand_state

    void init_mpz_globals_c "init_mpz_globals"()
    void clear_mpz_globals_c "clear_mpz_globals"()

########################################################

cdef object mpz_to_str(mpz_t x):
    """
    Convert a GMP integer to a Python string.
    """
    cdef char *s
    sig_on()
    s = mpz_get_str(NULL, 10, x)
    t = str(s)
    # Emulate sage_free() to avoid needing to include stdsage.pxi
    sig_block()
    libc.stdlib.free(s)
    sig_unblock()
    sig_off()
    return t
