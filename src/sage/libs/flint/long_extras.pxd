include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

from sage.libs.flint.flint cimport *


cdef extern from "FLINT/long_extras.h":

    cdef unsigned long z_randint(unsigned long limit)

    cdef unsigned long z_randbits(unsigned long bits)

    cdef unsigned long z_randprime(unsigned long bits)

    cdef double z_precompute_inverse(unsigned long n)

    cdef double z_precompute_inverse2(unsigned long n)

    cdef double z_ll_precompute_inverse(unsigned long n)

    cdef unsigned long z_addmod(unsigned long a, unsigned long b, unsigned long p)

    cdef unsigned long z_submod(unsigned long a, unsigned long b, unsigned long p)

    cdef unsigned long z_negmod(unsigned long a, unsigned long p)

    cdef unsigned long z_mod_precomp(unsigned long a, unsigned long n, double ninv)

    cdef unsigned long z_div2_precomp(unsigned long a, unsigned long n, double ninv)

    cdef unsigned long z_mod2_precomp(unsigned long a, unsigned long n, double ninv)

    cdef unsigned long z_ll_mod_precomp(unsigned long a_hi, unsigned long a_lo, unsigned long n, double ninv)

    cdef unsigned long z_mulmod_precomp(unsigned long a, unsigned long b, unsigned long n, double ninv)

    cdef unsigned long z_mulmod2_precomp(unsigned long a, unsigned long b, unsigned long n, double ninv)

    cdef unsigned long z_powmod(unsigned long a, long exp, unsigned long n)

    cdef unsigned long z_powmod2(unsigned long a, long exp, unsigned long n)

    cdef unsigned long z_powmod_precomp(unsigned long a, long exp, unsigned long n, double ninv)

    cdef unsigned long z_powmod2_precomp(unsigned long a, long exp, unsigned long n, double ninv)

    cdef int z_legendre_precomp(unsigned long a, unsigned long p, double pinv)

    cdef int z_jacobi(long x, unsigned long y)

    cdef int z_ispseudoprime_fermat(unsigned long n, unsigned long i)

    cdef int z_isprime(unsigned long n)

    cdef int z_isprime_precomp(unsigned long n, double ninv)

    cdef unsigned long z_nextprime(unsigned long n)

    cdef int z_isprime_pocklington(unsigned long n, unsigned long iterations)

    cdef int z_ispseudoprime_lucas_ab(unsigned long n, int a, int b)

    cdef int z_ispseudoprime_lucas(unsigned long n)

    cdef unsigned long z_pow(unsigned long a, unsigned long exp)

    cdef unsigned long z_sqrtmod(unsigned long a, unsigned long p)

    cdef unsigned long z_cuberootmod(unsigned long * cuberoot1, unsigned long a, unsigned long p)

    cdef unsigned long z_invert(unsigned long a, unsigned long p)

    cdef long z_gcd_invert(long* a, long x, long y)

    cdef long z_extgcd(long* a, long* b, long x, long y)

    cdef unsigned long z_gcd(long x, long y)

    cdef unsigned long z_intsqrt(unsigned long r)

    cdef int z_issquare(long x)

    cdef unsigned long z_CRT(unsigned long x1, unsigned long n1, unsigned long x2, unsigned long n2)

    cdef int z_issquarefree(unsigned long n)

    cdef int z_remove_precomp(unsigned long * n, unsigned long p, double pinv)

    cdef int z_remove(unsigned long * n, unsigned long p)

    cdef unsigned long z_factor_SQUFOF(unsigned long n)

    cdef unsigned long z_primitive_root(unsigned long p)

    cdef unsigned long z_primitive_root_precomp(unsigned long p, double p_inv)
