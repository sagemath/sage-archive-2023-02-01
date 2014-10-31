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

cdef int mpq_rational_reconstruction(mpq_t answer, mpz_t a, mpz_t m) except -1:
    """
    Set answer to the unique mpq is a modulo m such that ...

    We assume that answer has been mpq_init'd.
    If the rational reconstruction doesn't exist,
    raises a ValueError
    """
##     #debug
##     cdef char* s
##     s = mpz_get_str(NULL, 10, a)
##     t = str(s)
##     print 'a = ', t
##     s = mpz_get_str(NULL, 10, m)
##     t = str(s)
##     print 'm = ', t
##     #debug
    sig_on()
    mpz_mod(a, a, m)     # a = a % m
    if mpz_sgn(a) == 0 or mpz_sgn(m) == 0:    # a or m is zero
        mpq_set_si(answer, 0, 1)              # return 0
        sig_off()
        return 0
    if mpz_sgn(m) < 0:  # m negative
        mpz_neg(m, m)   # replace m by -m
    if mpz_sgn(a) < 0:  # a negative
        mpz_sub(a, m, a)  # replace a by m - a
    if mpz_cmp_si(a, 1) == 0:   # if a is 1
        mpq_set_si(answer, 1, 1)
        sig_off()
        return 0

    mpz_set(u, m)       # u = m
    mpz_set(v, a)       # v = a
    mpz_set_si(u0,1); mpz_set_si(u1,0); mpz_set(u2,u)
    mpz_set_si(v0,0); mpz_set_si(v1,1); mpz_set(v2,v)
    mpz_fdiv_q_ui(m2, m, 2)
    while 1:
        mpz_mul(ssqr, v2, v2)
        if mpz_cmp(ssqr, m2) <= 0:
            break
        mpz_fdiv_q(q, u2, v2)                  # q = floor of u2/v2
        mpz_mul(x, q, v0); mpz_sub(t0, u0, x)  # t0 = u0-q*v0
        mpz_mul(x, q, v1); mpz_sub(t1, u1, x)  # t0 = u1-q*v1
        mpz_mul(x, q, v2); mpz_sub(t2, u2, x)  # t0 = u2-q*v2
        mpz_set(u0,v0); mpz_set(u1,v1); mpz_set(u2,v2)   # permute
        mpz_set(v0,t0); mpz_set(v1,t1); mpz_set(v2,t2)

    mpz_abs(x, v1)
    mpz_set(y, v2)
    if mpz_sgn(v1)<0:
        mpz_neg(y, y)
    mpz_mul(ssqr, x, x)
    mpz_gcd(q, x, y)
    if mpz_cmp(ssqr, m2) <= 0 and mpz_cmp_si(q, 1) == 0:
        # return y/x
        mpq_set_z(answer, y)
        mpq_set_z(tmp, x)
        mpq_div(answer, answer, tmp)
        sig_off()
        return 0

    sig_off()
    raise ValueError, "Rational reconstruction of %s (mod %s) does not exist."%(mpz_to_str(a),mpz_to_str(m))
