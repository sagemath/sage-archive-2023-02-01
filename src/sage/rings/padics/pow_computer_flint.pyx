
from sage.rings.integer cimport Integer

cdef class PowComputer_flint_base(PowComputer_class):
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
