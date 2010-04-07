include 'sage/ext/cdefs.pxi'
include 'sage/ext/stdsage.pxi'
include 'sage/ext/interrupt.pxi'
include 'sage/ext/gmp.pxi'
include 'sage/libs/flint/fmpz_poly.pxi'

from sage.rings.rational_field import QQ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import primes, bernoulli
from sage.rings.integer cimport Integer
from sage.libs.flint.fmpz_poly cimport Fmpz_poly


cpdef eisenstein_series_poly(int k, int prec = 10) :
    r"""
    Return the q-expansion up to precision ``prec`` of the weight `k`
    Eisenstein series, as a FLINT :class:`~sage.libs.flint.fmpz_poly.Fmpz_poly`
    object, normalised so the coefficients are integers with no common factor.

    Used internally by the functions
    :func:`~sage.modular.modform.eis_series.eisenstein_series_qexp` and
    :func:`~sage.modular.modform.vm_basis.victor_miller_basis`; see the
    docstring of the former for further details.

    EXAMPLES :

        sage: from sage.modular.modform.eis_series_cython import eisenstein_series_poly
        sage: eisenstein_series_poly(12, prec=5)
        5  691 65520 134250480 11606736960 274945048560
    """
    cdef mpz_t *val = <mpz_t *>malloc(prec * sizeof(mpz_t))
    cdef mpz_t one, mult, term, last, term_m1, last_m1
    cdef unsigned long int expt
    cdef long ind, ppow, int_p
    cdef int i
    cdef Fmpz_poly res = PY_NEW(Fmpz_poly)

    if k%2 or k < 2:
        raise ValueError, "k (=%s) must be an even positive integer"%k
    if prec < 0:
        raise ValueError, "prec (=%s) must be an even nonnegative integer"%prec
    if (prec == 0):
        return PY_NEW(Fmpz_poly)

    mpz_init(one)
    mpz_init(term)
    mpz_init(last)
    mpz_init(mult)
    mpz_init(term_m1)
    mpz_init(last_m1)

    for i from 0 <= i < prec :
        mpz_init(val[i])
        mpz_set_si(val[i], 1)

    mpz_set_si(one, 1)

    expt = <unsigned long int>(k - 1)
    a0 = - bernoulli(k) / (2*k)

    for p in primes(1,prec) :
        int_p = int(p)
        ppow = <long int>int_p

        mpz_set_si(mult, int_p)
        mpz_pow_ui(mult, mult, expt)
        mpz_mul(term, mult, mult)
        mpz_set(last, mult)

        while (ppow < prec):
            ind = ppow
            mpz_sub(term_m1, term, one)
            mpz_sub(last_m1, last, one)
            while (ind < prec):
                mpz_mul(val[ind], val[ind], term_m1)
                mpz_fdiv_q(val[ind], val[ind], last_m1)
                ind += ppow
            ppow *= int_p
            mpz_set(last, term)
            mpz_mul(term, term, mult)

    mpz_clear(one)
    mpz_clear(term)
    mpz_clear(last)
    mpz_clear(mult)
    mpz_clear(term_m1)
    mpz_clear(last_m1)

    fmpz_poly_set_coeff_mpz(res.poly, prec-1, val[prec-1])
    for i from 1 <= i < prec - 1 :
        fmpz_poly_set_coeff_mpz(res.poly, i, val[i])

    fmpz_poly_scalar_mul_mpz(res.poly, res.poly, (<Integer>(a0.denominator())).value)
    fmpz_poly_set_coeff_mpz(res.poly, 0, (<Integer>(a0.numerator())).value)

    return res
