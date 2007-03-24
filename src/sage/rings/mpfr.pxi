#*****************************************************************************
# Headers.  When you past things in here from mpfr, be sure
# to remove const's, since those aren't allowed in pyrex.  Also, it can be
# challenging figuring out how to modify things from mpfr.h to be valid pyrex
# code.    Note that what is here is only used for generating the C code.
# The C compiler doesn't see any of this -- it only sees mpfr.h and stdlib.h
#*****************************************************************************

cdef extern from "mpfr.h":
    ctypedef void* mpq_t
    ctypedef void* mpz_t
    ctypedef struct __mpfr_struct:
        pass
    #ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct* mpfr_t
    ctypedef mpfr_t mpfr_ptr
    ctypedef mpfr_t mpfr_srcptr
    ctypedef enum mpfr_rnd_t:
        GMP_RNDN = 0
        GMP_RNDZ = 1
        GMP_RNDU = 2
        GMP_RNDD = 3
        GMP_RND_MAX = 4
        GMP_RNDNA = -1

    ctypedef mpfr_rnd_t mp_rnd_t
    ctypedef long int mp_exp_t
    ctypedef long mp_prec_t

    int MPFR_PREC_MIN, MPFR_PREC_MAX

    #mp_rnd_t GMP_RNDZ, GMP_RNDN, GMP_RNDU, GMP_RNDD

    void mpfr_init (mpfr_t x)
    void mpfr_init2 (mpfr_t x, mp_prec_t prec)
    void mpfr_clear (mpfr_t x)

    int mpfr_set (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_set_si (mpfr_t rop, long int op, mp_rnd_t rnd)
    int mpfr_set_ui (mpfr_t rop, long unsigned int op, mp_rnd_t rnd)
    int mpfr_set_str (mpfr_t rop, char *s, int base, mp_rnd_t rnd)
    void mpfr_set_inf (mpfr_t x, int sign)
    void mpfr_set_nan (mpfr_t x)
    int mpfr_set_z (mpfr_t rop, mpz_t op, mp_rnd_t rnd)
    int mpfr_set_q (mpfr_t rop, mpq_t op, mp_rnd_t rnd)
    int mpfr_set_d (mpfr_t rop, double op, mp_rnd_t rnd)

    char * mpfr_get_str (char *str, mp_exp_t *expptr, int base, size_t n, mpfr_t op, mp_rnd_t rnd)
    size_t mpfr_out_str (int *stream, int base, size_t n, mpfr_t op, mp_rnd_t rnd)
    void mpfr_free_str (char *str)

    void mpfr_get_z(mpz_t rop, mpfr_t op, mp_rnd_t rnd)
    mp_exp_t mpfr_get_z_exp(mpz_t rop, mpfr_t op)

    # Arithmetic
    int mpfr_add (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_sub (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_mul (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_div (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_mul_2exp (mpfr_t rop, mpfr_t op1, unsigned long int op2,mp_rnd_t rnd)
    int mpfr_div_2exp(mpfr_t rop, mpfr_t op1, unsigned long int op2,mp_rnd_t rnd)
    int mpfr_mul_2si (mpfr_t top, mpfr_t op1, long int opt2, mp_rnd_t rnd)

    # constants
    int mpfr_const_log2 (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_pi (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_euler (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_catalan (mpfr_t rop, mp_rnd_t rnd)

    # Special functions
    int mpfr_sqrt (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    #int mpfr_sqrt_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));
    int mpfr_cbrt (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_root (mpfr_t rop, mpfr_t op, unsigned long int k, mp_rnd_t rnd)

    int mpfr_log (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_log2 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_log10 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_exp (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_exp2 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_exp10 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_cos (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_sin (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_tan (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_sin_cos (mpfr_t rop, mpfr_t op, mpfr_t, mp_rnd_t rnd)

    int mpfr_acos (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_asin (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_atan (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_cosh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_sinh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_tanh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_atanh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_acosh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_asinh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_agm (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_gamma (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_zeta (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_erf (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_pow (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    #int mpfr_ui_pow _PROTO ((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t));
    #int mpfr_pow_si _PROTO ((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t));

    #int mpfr_min _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_max _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));

    int mpfr_fac_ui (mpfr_t rop, unsigned long int op, mp_rnd_t rnd)

    int mpfr_abs(mpfr_ptr rop, mpfr_srcptr op, mp_rnd_t rnd)
    int mpfr_sgn(mpfr_t op)
    #define mpfr_abs(a,b,r) mpfr_set4(a,b,r,1)
    #define mpfr_sgn(x) mpfr_cmp_ui(x,0)

    int mpfr_round (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_trunc (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_ceil (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_floor (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_frac (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    # Status functions
    int mpfr_nan_p (mpfr_t op)
    int mpfr_inf_p (mpfr_t op)
    int mpfr_number_p (mpfr_t op)
    int mpfr_zero_p (mpfr_t op)

    double mpfr_get_d (mpfr_t op, mp_rnd_t rnd)

    # Miscellaneous
    void mpfr_nexttoward (mpfr_t X, mpfr_t Y)
    void mpfr_nextabove (mpfr_t X)
    void mpfr_nextbelow (mpfr_t X)

    int mpfr_set_exp (mpfr_t op, mp_exp_t E)

    # Operators

    int mpfr_neg (mpfr_ptr rop, mpfr_srcptr op, mp_rnd_t rnd)
    # int mpfr_eq (mpfr_srcptr rop, mpfr_srcptr op, unsigned long i)
    int mpfr_less_p (mpfr_t op1, mpfr_t op2)
    int mpfr_lessequal_p (mpfr_t op1, mpfr_t op2)
    int mpfr_cmp (mpfr_t op1, mpfr_t op2)
