from sage.libs.gmp.all cimport mpz_t

ctypedef struct MPopts:
    long prec
    int rounding

cdef mpz_set_integer(mpz_t v, x)
cdef inline mpzi(mpz_t n)
cdef inline mpzl(mpz_t n)
cdef inline str rndmode_to_python(int rnd)
cdef inline rndmode_from_python(str rnd)

ctypedef struct MPF:
    mpz_t man
    mpz_t exp
    int special

cdef inline void MPF_init(MPF *x)
cdef inline void MPF_clear(MPF *x)
cdef inline void MPF_set(MPF *dest, MPF *src)
cdef inline void MPF_set_zero(MPF *x)
cdef inline void MPF_set_nan(MPF *x)
cdef inline void MPF_set_inf(MPF *x)
cdef inline void MPF_set_ninf(MPF *x)
cdef MPF_set_si(MPF *x, long n)
cdef MPF_set_int(MPF *x, n)
cdef MPF_set_man_exp(MPF *x, man, exp)
cdef MPF_set_tuple(MPF *x, tuple value)
cdef MPF_to_tuple(MPF *x)
cdef MPF_set_double(MPF *r, double x)
cdef double MPF_to_double(MPF *x, bint strict)
cdef MPF_to_fixed(mpz_t r, MPF *x, long prec, bint truncate)
cdef int MPF_sgn(MPF *x)
cdef void MPF_neg(MPF *r, MPF *s)
cdef void MPF_abs(MPF *r, MPF *s)
cdef MPF_normalize(MPF *x, MPopts opts)
cdef MPF_add(MPF *r, MPF *s, MPF *t, MPopts opts)
cdef MPF_sub(MPF *r, MPF *s, MPF *t, MPopts opts)
cdef bint MPF_eq(MPF *s, MPF *t)
cdef bint MPF_ne(MPF *s, MPF *t)
cdef int MPF_cmp(MPF *s, MPF *t)
cdef bint MPF_lt(MPF *s, MPF *t)
cdef bint MPF_le(MPF *s, MPF *t)
cdef bint MPF_gt(MPF *s, MPF *t)
cdef bint MPF_ge(MPF *s, MPF *t)
cdef MPF_mul(MPF *r, MPF *s, MPF *t, MPopts opts)
cdef MPF_div(MPF *r, MPF *s, MPF *t, MPopts opts)
cdef MPF_sqrt(MPF *r, MPF *s, MPopts opts)
cdef MPF_hypot(MPF *r, MPF *a, MPF *b, MPopts opts)
cdef MPF_pow_int(MPF *r, MPF *x, mpz_t n, MPopts opts)
cdef MPF_set_double(MPF *r, double x)
