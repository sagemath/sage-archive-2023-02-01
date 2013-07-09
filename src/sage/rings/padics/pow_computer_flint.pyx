
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

cdef class PowComputer_flint(PowComputer_class):
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        self._initialized = 0
        fmpz_init(self.fprime)
        fmpz_init(self.ftmp)
        fmpz_init(self.ftmp2)
        mpz_init(self.temp_m)
        mpz_init(self.top_power)
        fmpz_set_mpz(self.fprime, prime.value)
        padic_ctx_init(self.ctx, self.fprime, prec_cap, PADIC_SERIES)
        mpz_pow_ui(self.top_power, prime.value, prec_cap)
        self._initialized = 1

    def __dealloc__(self):
        if self._initialized:
            fmpz_clear(self.fprime)
            fmpz_clear(self.ftmp)
            fmpz_clear(self.ftmp2)
            mpz_clear(self.temp_m)
            mpz_clear(self.top_power)
            padic_ctx_clear(self.ctx)

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n):
        cdef padic_ctx_struct ctx = (<padic_ctx_struct*>self.ctx)[0]
        if ctx.min <= n and n < ctx.max:
            self.ftmp[0] = (ctx.pow + (n - ctx.min))[0]
        else:
            fmpz_pow_ui(self.ftmp, self.fprime, n)
        return &self.ftmp

    cdef mpz_t* pow_mpz_t_tmp(self, unsigned long n):
        cdef padic_ctx_struct ctx = (<padic_ctx_struct*>self.ctx)[0]
        if ctx.min <= n and n < ctx.max:
            self.ftmp[0] = (ctx.pow + (n - ctx.min))[0]
        else:
            fmpz_pow_ui(self.ftmp, self.fprime, n)
        fmpz_get_mpz(self.temp_m, self.ftmp)
        return &(self.temp_m)

    cdef mpz_t* pow_mpz_t_top(self):
        return &self.top_power

    cdef unsigned long capdiv(self, unsigned long n):
        if self.e == 1: return n
        if n == 0: return 0
        return (n-1) / self.e + 1

cdef class PowComputer_flint_1step(PowComputer_flint):
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly, shift_seed=None):
        self._initialized = 0
        cdef Polynomial_integer_dense_flint poly = _poly
        cdef long length = fmpz_poly_length(poly.__poly)
        self.deg = length - 1
        fmpz_poly_init2(self.modulus, length)
        self._initialized = 1
        fmpz_poly_set(self.modulus, poly.__poly)
        fmpz_poly_init(self.tmp_poly)
        self._initialized = 2
        sig_on()
        self._moduli = <fmpz_poly_t*>sage_malloc(sizeof(fmpz_poly_t) * (cache_limit + 2))
        sig_off()
        if self._moduli == NULL:
            raise MemoryError
        self._initialized = 3

        cdef Py_ssize_t i
        cdef fmpz* coeffs = (<fmpz_poly_struct*>self._moduli[i])[0].coeffs
        fmpz_one(self.ftmp)
        for i from 1 <= i <= cache_limit:
            fmpz_mul(self.ftmp, self.ftmp, self.fprime)
            fmpz_poly_init2(self._moduli[i], length)
            self._initialized += 1
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[i])[0].coeffs, coeffs, length, self.ftmp)
            _fmpz_poly_set_length(self._moduli[i], length)
        # We use cache_limit + 1 as a temporary holder
        fmpz_poly_init2(self._moduli[cache_limit+1], length)
        self._initialized += 1
        _fmpz_poly_set_length(self._moduli[cache_limit+1], length)

    def __dealloc__(self):
        cdef int init = self._initialized
        if init > 0: fmpz_poly_clear(self.modulus)
        if init > 1: fmpz_poly_clear(self.tmp_poly)
        cdef Py_ssize_t i
        for i from 1 <= i <= self.cache_limit+1:
            if init >= 3+i: fmpz_poly_clear(self._moduli[i])
        if init > 2: sage_free(self._moduli)

    cdef fmpz_poly_t* get_modulus(self, unsigned long n):
        if n <= self.cache_limit:
            return &(self._moduli[n])
        else:
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[self.cache_limit+1])[0].coeffs,
                                      (<fmpz_poly_struct*>self.modulus)[0].coeffs,
                                      self.deg + 1,
                                      self.pow_fmpz_t_tmp(n)[0])
            return &(self._moduli[self.cache_limit+1])

    cdef fmpz_poly_t* get_modulus_capdiv(self, unsigned long n):
        return self.get_modulus(self.capdiv(n))

    def polynomial(self):
        cdef Polynomial_integer_dense_flint ans = PY_NEW(Polynomial_integer_dense_flint)
        fmpz_poly_set(ans.__poly, self.modulus)
        return ans

cdef class PowComputer_flint_unram(PowComputer_flint_1step):
    def __init__(self):
        self.e = 1
        self.f = fmpz_poly_degree(self.modulus)

cdef class PowComputer_flint_eis(PowComputer_flint_1step):
    def __init__(self):
        self.e = fmpz_poly_degree(self.modulus)
        self.f = 1
