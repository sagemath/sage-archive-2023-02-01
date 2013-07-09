
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

cdef class PowComputer_flint(PowComputer_class):
    """
    A PowComputer for use in `p`-adics implemented via FLINT.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
        sage: A = PowComputer_flint_maker(5, 20, 20, 20, False); A
        FLINT PowComputer for 5
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False)
            sage: TestSuite(A).run()
        """
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
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False)
            sage: del A
        """
        if self._initialized:
            fmpz_clear(self.fprime)
            fmpz_clear(self.ftmp)
            fmpz_clear(self.ftmp2)
            mpz_clear(self.temp_m)
            mpz_clear(self.top_power)
            padic_ctx_clear(self.ctx)

    def _repr_(self):
        """
        String representation of this powcomputer.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False); A
            FLINT PowComputer for 5
        """
        return "FLINT PowComputer for %s"%(self.prime)

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n):
        """
        Returns a pointer to a FLINT ``fmpz_t`` holding `p^n`.

        Analogous to
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`
        but with FLINT ``fmpz_t`` rather than GMP ``mpz_t``.  The same
        important warnings apply.
        """
        cdef padic_ctx_struct ctx = (<padic_ctx_struct*>self.ctx)[0]
        if ctx.min <= n and n < ctx.max:
            self.ftmp[0] = (ctx.pow + (n - ctx.min))[0]
        else:
            fmpz_pow_ui(self.ftmp, self.fprime, n)
        return &self.ftmp

    cdef mpz_t* pow_mpz_t_tmp(self, unsigned long n):
        """
        Returns a pointer to an ``mpz_t`` holding `p^n`.

        See
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`
        for important warnings.
        """
        fmpz_get_mpz(self.temp_m, self.pow_fmpz_t_tmp(n)[0])
        return &(self.temp_m)

    cdef mpz_t* pow_mpz_t_top(self):
        """
        Returns a pointer to an ``mpz_t`` holding `p^N`, where `N` is
        the precision cap.
        """
        return &self.top_power

    cdef unsigned long capdiv(self, unsigned long n):
        """
        Returns ceil(n / e).
        """
        if self.e == 1: return n
        if n == 0: return 0
        return (n-1) / self.e + 1

cdef class PowComputer_flint_1step(PowComputer_flint):
    """
    A PowComputer for a `p`-adic extension defined by a single polynomial.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
        sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
        sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, _poly, shift_seed=None):
        """
        Memory initialization.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f)
            sage: TestSuite(A).run()
        """
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
        fmpz_init(self._an)
        self._initialized = 4
        fmpz_poly_get_coeff_fmpz(self._an, self.modulus, self.deg)
        if fmpz_is_one(self._an):
            self.is_monic = True
        else:
            sig_on()
            self._inv_an = <fmpz_t*>sage_malloc(sizeof(fmpz_t) * (cache_limit + 2))
            sig_off()
            if self._inv_an == NULL:
                raise MemoryError
        self._initialized = 5

        cdef Py_ssize_t i
        cdef fmpz* coeffs = (<fmpz_poly_struct*>self.modulus)[0].coeffs
        fmpz_one(self.ftmp)
        for i from 1 <= i <= cache_limit:
            fmpz_mul(self.ftmp, self.ftmp, self.fprime)
            fmpz_poly_init2(self._moduli[i], length)
            self._initialized += 1
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[i])[0].coeffs, coeffs, length, self.ftmp)
            _fmpz_poly_set_length(self._moduli[i], length)
            if not self.is_monic:
                fmpz_init(self._inv_an[i])
                # Could use Newton lifting here
                fmpz_invmod(self._inv_an[i], self._an, self.pow_fmpz_t_tmp(i)[0])
            self._initialized += 1
        # We use cache_limit + 1 as a temporary holder
        fmpz_poly_init2(self._moduli[cache_limit+1], length)
        _fmpz_poly_set_length(self._moduli[cache_limit+1], length)
        self._initialized += 1
        if not self.is_monic:
            fmpz_init(self._inv_an[cache_limit+1])
        self._initialized += 1

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f)
            sage: del A
        """
        cdef int init = self._initialized
        if init > 0: fmpz_poly_clear(self.modulus)
        if init > 1: fmpz_poly_clear(self.tmp_poly)
        cdef Py_ssize_t i
        for i from 1 <= i <= self.cache_limit+1:
            if init >= 3+2*i: fmpz_poly_clear(self._moduli[i])
            if init >= 4+2*i and not self.is_monic: fmpz_clear(self._inv_an[i])
        if init > 2: sage_free(self._moduli)
        if init > 3: fmpz_clear(self._an)
        if init > 4 and not self.is_monic: sage_free(self._inv_an)

    def _repr_(self):
        """
        String representation of this powcomputer.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f); A
            FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
        """
        return "FLINT PowComputer for %s with polynomial %s"%(self.prime, self.polynomial())

    def __cmp__(self, other):
        """
        Comparison.

        Lexicographic on class, prime, precision cap, cache_limit and polynomial.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2; g = x^3 - (8 + 5^22)*x - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f)
            sage: B = PowComputer_flint_maker(5, 20, 20, 20, False, g)
            sage: A == B
            False
        """
        c = PowComputer_flint.__cmp__(self, other)
        if c: return c
        cdef PowComputer_flint_1step o = other
        if fmpz_poly_equal(self.modulus, o.modulus): return 0
        #return cmp(self.polynomial(), other.polynomial())
        return 1

    cdef fmpz_poly_t* get_modulus(self, unsigned long k):
        """
        Returns the defining polynomial reduced modulo `p^k`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        cdef long c
        if k <= self.cache_limit:
            return &(self._moduli[k])
        else:
            c = self.cache_limit+1
            _fmpz_vec_scalar_mod_fmpz((<fmpz_poly_struct*>self._moduli[c])[0].coeffs,
                                      (<fmpz_poly_struct*>self.modulus)[0].coeffs,
                                      self.deg + 1,
                                      self.pow_fmpz_t_tmp(k)[0])
            return &(self._moduli[c])

    cdef fmpz_poly_t* get_modulus_capdiv(self, unsigned long k):
        """
        Returns the defining polynomial reduced modulo `p^a`, where
        `a` is the ceiling of `k/e`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        return self.get_modulus(self.capdiv(k))

    cdef fmpz_t* get_inv_an(self, unsigned long k):
        """
        Returns the inverse of the leading coefficient modulo `p^k`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        cdef long c
        if self.is_monic: return &(self._an)
        if k <= self.cache_limit:
            return &(self._inv_an[k])
        else:
            c = self.cache_limit+1
            fmpz_poly_get_coeff_fmpz(self._inv_an[c], self.modulus, self.deg)
            # Could use Newton lifting here
            fmpz_invmod(self._inv_an[c], self._inv_an[c], self.pow_fmpz_t_tmp(k)[0])
            return &(self._inv_an[c])

    cdef fmpz_t* get_inv_an_capdiv(self, unsigned long k):
        """
        Returns the inverse of the leading coefficient modulo `p^k`, where
        `k` is the ceiling of `n/e`.

        The same warnings apply as for
        :meth:`sage.rings.padics.pow_computer.PowComputer_class.pow_mpz_t_tmp`.
        """
        if self.is_monic: return &(self._an)
        return self.get_inv_an(self.capdiv(k))

    def polynomial(self, _n=None, var='x'):
        """
        Returns the polynomial attached to this ``PowComputer``, possibly reduced modulo a power of `p`.

        INPUT:

        - ``_n`` -- (default ``None``) an integer, the power of `p`
          modulo which to reduce.

        - ``var`` -- (default ``'x'``) the variable for the returned polynomial

        .. NOTE::

            From Cython you should use :meth:`get_modulus` instead for
            speed.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
            sage: R.<y> = ZZ[]; f = y^3 - 8*y - 2
            sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f)
            sage: A.polynomial()
            x^3 - 8*x - 2
            sage: A.polynomial(var='y')
            y^3 - 8*y - 2
            sage: A.polynomial(2)
            x^3 + 17*x + 23
        """
        from sage.rings.all import ZZ
        R = ZZ[var]
        x = R.gen()
        cdef Polynomial_integer_dense_flint ans = (<Polynomial_integer_dense_flint?>x)._new()
        if _n is None:
            fmpz_poly_set(ans.__poly, self.modulus)
        else:
            fmpz_poly_set(ans.__poly, self.get_modulus(_n)[0])
        return ans

cdef class PowComputer_flint_unram(PowComputer_flint_1step):
    """
    A PowComputer for a `p`-adic extension defined by an unramified polynomial.

    For a description of inputs see :func:`PowComputer_flint_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_flint import PowComputer_flint_maker
        sage: R.<x> = ZZ[]; f = x^3 - 8*x - 2
        sage: A = PowComputer_flint_maker(5, 20, 20, 20, False, f); A
        FLINT PowComputer for 5 with polynomial x^3 - 8*x - 2
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None):
        self.e = 1
        self.f = fmpz_poly_degree(self.modulus)

cdef class PowComputer_flint_eis(PowComputer_flint_1step):
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly=None, shift_seed=None):
        self.e = fmpz_poly_degree(self.modulus)
        self.f = 1

def PowComputer_flint_maker(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly):
    if poly is None:
        return PowComputer_flint(prime, cache_limit, prec_cap, ram_prec_cap, in_field)
    else:
        return PowComputer_flint_unram(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly)
