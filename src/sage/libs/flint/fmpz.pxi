include "sage/libs/ntl/decl.pxi"

from sage.libs.gmp.mpz cimport mpz_t
from sage.libs.flint.nmod_poly cimport mp_srcptr, nmod_t

cdef extern from "flint/fmpz.h":

    ctypedef long fmpz
    ctypedef fmpz[1] fmpz_t

    ctypedef gmp_randstate_t fmpz_randstate_t

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        long n
        mp_bitcnt_t norm

    ctypedef fmpz_preinvn_struct[1] fmpz_preinvn_t

    void _fmpz_clear_mpz(fmpz f)

    void _fmpz_cleanup_mpz_content()

    void _fmpz_cleanup()

    void fmpz_init(fmpz_t f)

    void fmpz_init2(fmpz_t f, unsigned long limbs)

    void fmpz_init_set(fmpz_t f, fmpz_t g)

    void fmpz_init_set_ui(fmpz_t f, unsigned long g)

    void fmpz_clear(fmpz_t f)

    long fmpz_get_si(fmpz_t f)

    unsigned long fmpz_get_ui(fmpz_t f)

    void fmpz_set_si(fmpz_t f, long val)

    void fmpz_set_ui(fmpz_t f, unsigned long val)

    void fmpz_neg_ui(fmpz_t f, unsigned long val)

    void fmpz_set_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)

    void fmpz_neg_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)

    void fmpz_get_mpz(mpz_t x, fmpz_t f)

    void fmpz_set_mpz(fmpz_t f, mpz_t x)

    double fmpz_get_d(fmpz_t f)

    void fmpz_set_d(fmpz_t f, double c)

    int fmpz_set_str(fmpz_t f, char * str, int b)

    void flint_mpz_init_set_readonly(mpz_t z, fmpz_t f)

    void flint_mpz_clear_readonly(mpz_t z)

    void fmpz_init_set_readonly(fmpz_t f, mpz_t z)

    void fmpz_clear_readonly(fmpz_t f)

    int fmpz_abs_fits_ui(fmpz_t f)

    int fmpz_fits_si(fmpz_t f)

    void fmpz_zero(fmpz_t f)

    void fmpz_one(fmpz_t f)

    int fmpz_is_zero(fmpz_t f)

    int fmpz_is_one(fmpz_t f)
    int fmpz_sgn(fmpz_t f)


    void fmpz_neg(fmpz_t rop, fmpz_t op)
    void fmpz_get_mpz(mpz_t rop, const fmpz_t op)
    void fmpz_set_mpz(fmpz_t rop, const mpz_t op)

    void fmpz_add_ui(fmpz_t f, fmpz_t g, unsigned long c)

    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, const unsigned long g)

    void fmpz_abs(fmpz_t f1, const fmpz_t f2)
    void fmpz_set(fmpz_t f, const fmpz_t val)
    void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    int fmpz_sgn(const fmpz_t f)
    double fmpz_get_d(const fmpz_t f)
    double fmpz_get_d_2exp(long *exp, const fmpz_t f)
    void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)
    unsigned long fmpz_fdiv_ui(const fmpz_t g, unsigned long x)
    void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)

    size_t fmpz_sizeinbase(const fmpz_t f, int b)
    char * fmpz_get_str(char * str, int b, const fmpz_t f)
    int fmpz_set_str(fmpz_t f, const char * str, int b)
    int fmpz_cmp(const fmpz_t f, const fmpz_t g)
    int fmpz_cmp_si(const fmpz_t f, int g)
    int fmpz_cmp_ui(const fmpz_t f, unsigned int g)

    void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
