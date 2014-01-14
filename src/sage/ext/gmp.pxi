# gmp.inc -- misc. useful GMP functions that don't depend on
#            other things and can be included in other files
# USAGE
#
# Include the following at the top of client .pyx file.
#
#    include "cdefs.pxi"
#    include "gmp.pxi"
#
# to include this in a file.

include 'sage/ext/interrupt.pxi'

cimport libc.stdlib
from cpython.mem cimport *

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

cdef int mpz_print(mpz_t x) except -1:
    print mpz_to_str(x)
    return 0

cdef int next_probab_prime_int(int p) except -1:
    """
    Use GMP to compute the next int prime after p.
    """
    cdef mpz_t p1, p2
    mpz_init_set_si(p1, p)
    mpz_init(p2)
    mpz_nextprime(p2, p1)
    cdef int p3
    p3 = mpz_get_si(p2)
    mpz_clear(p1)
    mpz_clear(p2)
    return p3

cdef int previous_probab_prime_int(int p) except -1:
    """
    Use GMP to compute the next int prime after p.
    """
    p = p - 1
    if p < 2:
        raise RuntimeError, "no prime < 2"
    if p == 2:
        return 2
    if p % 2 == 0:
        p = p - 1

    cdef mpz_t p1
    mpz_init_set_si(p1, p)
    while mpz_probab_prime_p(p1, 10) == 0:
        mpz_sub_ui(p1, p1, 2)
    cdef int p2
    p2 = mpz_get_si(p1)
    mpz_clear(p1)
    return p2

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


cdef int mpz_crt_vec(mpz_t **z, mpz_t *x, mpz_t *y, int r, mpz_t m, mpz_t n) except -1:
    """
    Allocates memory for z[0] and does exactly the same thing as mpz_crt to fill z in.
    Here r is the length of the arrays x and y.
    The caller of this function is responsible for using PyMem_Free to de-allocate
    the memory used by z.
    A MemoryError is raised if memory cannot be allocated.
    An ArithmeticError is raised if gcd(m,n) =/= 1.
    """
    cdef int i
    cdef mpz_t g, s, t, mn

    mpz_init(g); mpz_init(s); mpz_init(t); mpz_init(mn)
    mpz_mul(mn, m, n)
    mpz_gcdext(g, s, t, m, n)
    if mpz_cmp_si(g,1) != 0:
        mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(mn)
        raise ArithmeticError, "Moduli must be coprime (but gcd=%s)."%mpz_to_str(g)
    z[0] = <mpz_t*> PyMem_Malloc(sizeof(mpz_t)*r)
    if z[0] == <mpz_t*> 0:
        mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(mn)
        raise MemoryError
    # Now s*m + t*n = 1, so mod n, s*m = 1.
    # The answer is x + (y-x)*s*m, since mod m this is x, and mod n this is y
    for i from 0 <= i < r:
        mpz_init(z[0][i])
        mpz_sub(z[0][i], y[i], x[i])     # z = y-x
        mpz_mul(z[0][i], z[0][i], s)     # z = (y-x)*s
        mpz_mul(z[0][i], z[0][i], m)     # z = (y-x)*s*m
        mpz_add(z[0][i], z[0][i], x[i])  # z = z + (y-x)*s*m
        mpz_mod(z[0][i], z[0][i], mn)    # z = z (mod mn)

    mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(mn)
    return 0

cdef int mpz_crt_intvec(mpz_t* answer, mpz_t* modulus,
                        unsigned int *v, unsigned int *m, int n) except -1:
    cdef int i

    mpz_set_ui(a1, v[0])
    mpz_set_ui(mod1, m[0])
    for i from 1 <= i < n:
        mpz_set_ui(a2, v[i])
        mpz_set_ui(sage_mod2, m[i])
        mpz_gcdext(g, s, t, mod1, sage_mod2)
        if mpz_cmp_si(g,1) != 0:
            raise ArithmeticError, "Moduli must be coprime (gcd=%s)."%mpz_to_str(g)
        # Set a1 = a1 + (a2-a1) * s * mod1
        mpz_sub(xx, a2, a1)
        mpz_mul(xx, xx, s)
        mpz_mul(xx, xx, mod1)
        mpz_add(a1, a1, xx)
        # Set mod1 = mod1*sage_mod2
        mpz_mul(mod1, mod1, sage_mod2)

    # Now a1 is the CRT for the whole list v and mod1 is the modulus
    mpz_set(answer[0], a1)
    mpz_set(modulus[0], mod1)
    return 0   # everything OK

cdef int mpq_using_crt_and_rr(mpq_t answer,
                      unsigned int *v, unsigned int *m, int n) except -1:
    """
    INPUT:
        answer -- an mpq where the result goes; we assume it has been mpq_init'd.
        v -- an array of n integers
        m -- an array of n moduli
    OUTPUT:
        A single GMP rational the reduces to each v[i] modulo each n[i], if
        such an element exists, or raises a ValueError exception.
        You must mpz_clear this return value when you are done with it.
    """
    # Algorithm:
    #   1. Use the CRT to find a single mpz that reduces to each v[i], modulo each m[i].
    #   2. Apply rational reconstruction to that mpz modulo the product of the m[i].
    mpz_crt_intvec(&crtrr_a, &crtrr_mod, v, m, n)
    mpq_rational_reconstruction(answer, crtrr_a, crtrr_mod)
    return 0

cdef int mpz_randrange(mpz_t val, int n1, int n2) except -1:
    mpz_set_si(rand_n, n2-n1)
    mpz_set_si(rand_n1, n1)
    mpz_urandomm(val, rand_state, rand_n)
    mpz_add(val, val, rand_n1)



