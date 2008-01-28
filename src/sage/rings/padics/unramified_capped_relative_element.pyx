

cdef class UnramifiedCappedRelativeElement(pAdicZZpXCRElement):


    cdef int _set_from_mpz_rel(self, mpz_t x, unsigned long relprec) except -1:
        if mpz_sgn(x) == 0:
            self._set_exact_zero()
            return
        cdef mpz_t tmp_m
        mpz_init(tmp_m)
        self.ordp = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        self.ordp *= self.prime_pow.e
        cdef ZZ_c tmp_z
        ZZ_construct(&tmp_z)
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        self.prime_pow.restore_context(relprec)
        ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
        ZZ_destruct(&tmp_z)
        self.relprec = relprec
        self._normalized = 1

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, unsigned long relprec) except -1:
        if mpz_sgn(x) == 0:
            self._set_inexact_zero(absprec)
        cdef long rprec
        cdef mpz_t tmp_m
        mpz_init(tmp_m)
        self.ordp = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        rprec = absprec - self.ordp
        if rprec > relprec:
            self.relprec = relprec
        elif rprec < 0:
            self._set_inexact_zero(absprec)
            return
        else:
            self.relprec = rprec
        cdef ZZ_c tmp_z
        ZZ_construct(&tmp_z)
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        self.prime_pow.restore_context(self.relprec)
        ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
        ZZ_destruct(&tmp_z)

    cdef int _set_from_mpq_rel(self, mpq_t x, unsigned long relprec) except -1:
        if mpq_sgn(x) == 0:
            self._set_exact_zero()
            return
        cdef mpz_t num_unit, den_unit
        cdef long num_ordp, den_ordp
        _sig_on
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(x), self.prime_pow.prime.value)
        _sig_off
        self.ordp = (num_ordp - den_ordp)
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError, "p divides the denominator"
        cdef ZZ_c num_zz, den_zz
        cdef ZZ_p_c tmp_zp
        ZZ_construct(&num_zz)
        ZZ_construct(&den_zz)
        mpz_to_ZZ(&num_zz, &num_unit)
        mpz_to_ZZ(&den_zz, &den_unit)
        self.prime_pow.restore_context(cap_div(relprec, self.prime_pow.e))
        ZZ_p_construct(&tmp_zp)
        ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
        ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
        ZZ_p_destruct(&tmp_zp)
        ZZ_destruct(&num_zz)
        ZZ_destruct(&den_zz)
        mpz_clear(num_unit)
        mpz_clear(den_unit)
        self.relprec = relprec
        self._normalized = 1

    cdef int _set_from_mpq_both(self, mpz_t x, long absprec, unsigned long relprec) except -1:
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return
        cdef mpz_t num_unit, den_unit
        cdef long num_ordp, den_ordp
        _sig_on
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(x), self.prime_pow.prime.value)
        _sig_off
        self.ordp = (num_ordp - den_ordp)
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError, "p divides the denominator"
        cdef ZZ_c num_zz, den_zz
        cdef ZZ_p_c tmp_zp
        ZZ_construct(&num_zz)
        ZZ_construct(&den_zz)
        mpz_to_ZZ(&num_zz, &num_unit)
        mpz_to_ZZ(&den_zz, &den_unit)
        self.relprec = absprec - self.ordp
        if relprec < self.relprec:
            self.relprec = relprec
        elif self.relprec < 0:
            self._set_inexact_zero(absprec)
        self.prime_pow.restore_context(cap_div(self.relprec, self.prime_pow.e))
        ZZ_p_construct(&tmp_zp)
        ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
        ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
        ZZ_p_destruct(&tmp_zp)
        ZZ_destruct(&num_zz)
        ZZ_destruct(&den_zz)
        mpz_clear(num_unit)
        mpz_clear(den_unit)
        self._normalized = 1
