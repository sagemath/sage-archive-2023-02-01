from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef class EisensteinCappedRelativeElement(pAdicZZpXCRElement):


    cdef int _set_from_mpz_rel(self, mpz_t x, unsigned long relprec) except -1:
        if mpz_sgn(x) == 0:
            self._set_exact_zero()
            return
        cdef mpz_t tmp_m
        cdef long val
        _sig_on
        mpz_init(tmp_m)
        self.ordp = mpz_remove(tmp_m, x, self.prime_pow.prime.value) * self.prime_pow.e
        mpz_clear(tmp_m)
        _sig_off
        self._set_absprec_rel(relprec)
        self._set_from_mpz_part2(x)

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, unsigned long relprec) except -1:
        if mpz_sgn(x) == 0:
            self._set_inexact_zero(absprec)
        cdef long aprec, val
        cdef mpz_t tmp_m
        mpz_init(tmp_m)
        val = mpz_remove(tmp_m, x, self.prime_pow.prime.value) * self.prime_pow.e

    cdef int _set_from_mpz_part2(self, mpz_t x) except -1:
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        ZZ_construct(&tmp_z)
        # I can't just put &x as the second argument in the mpz_to_ZZ call because it's a local variable.
        mpz_init(tmp_m)
        mpz_set(tmp_m, x)
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        # Note that we have set self.relprec to be self's absolute precision
        self.prime_pow.restore_context(cap_div(self.relprec, self.prime_pow.e))
        ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
        ZZ_destruct(&tmp_z)
        if self.ordp > 0:
            self._normalized = 0
            self._normalize()
        else:
            self._normalized = 1

    cdef int _set_from_mpq_rel(self, mpq_t x, unsigned long relprec) except -1:
        if mpq_sgn(x) == 0:
            self._set_exact_zero()
            return
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        self._set_prec_rel(relprec)
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_both(self, mpz_t x, long absprec, unsigned long relprec) except -1:
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        self._set_prec_both(absprec, relprec)
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_part1(self, mpz_t num_unit, mpz_t den_unit, mpq_t x) except -1:
        cdef long num_ordp, den_ordp
        _sig_on
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(x), self.prime_pow.prime.value)
        _sig_off
        if num_ordp > 0:
            # Because of the way _set_from_mpq_part2 works (ie by eis_shifting back)
            # we want to reset num_unit to include the power of p.
            mpz_set(num_unit, mpq_numref(x))
        self.ordp = (num_ordp - den_ordp) * self.prime_pow.e
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError, "p divides the denominator"

    cdef int _set_from_mpq_part2(self, mpz_t num_unit, mpz_t den_unit) except -1:
        cdef ZZ_c num_zz, den_zz
        cdef ZZ_p_c tmp_zp
        cdef ZZ_pX_c high_shifter
        cdef ZZ_pX_Modulus_c modulus
        cdef ntl_ZZ_pContext_class c
        cdef long i
        cdef mpz_t tmp_m
        mpz_init(tmp_m)
        cdef long val
        ZZ_construct(&num_zz)
        ZZ_construct(&den_zz)
        # I can't just put &num_unit as the second argument in the mpz_to_ZZ call since num_unit is passed in as an argument.
        mpz_set(tmp_m, num_unit)
        mpz_to_ZZ(&num_zz, &tmp_m)
        mpz_set(tmp_m, den_unit)
        mpz_to_ZZ(&den_zz, &tmp_m)
        if self.ordp >= 0:
            self.relprec += self.ordp
            self.prime_pow.restore_context(cap_div(self.relprec, self.prime_pow.e))
            # have to construct here, after the context has been restored
            ZZ_p_construct(&tmp_zp)
            # Note that num_zz actually includes the power of p because of our modification in _set_from_mpq_part1
            ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
            ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
            if self.ordp == 0:
                self._normalized == 1
            else:
                self.ordp = 0
                self._normalized = 0
                self._normalize()
        else:
            # In order to deal with a negative power of p, we use self.prime_pow.high_shifter[i], which stores (p/x^e)^(2^i)
            val = -self.ordp / self.prime_pow.e
            i = 0
            c = self.prime_pow.get_context(cap_div(self.relprec, self.prime_pow.e))
            c.restore_c()
            modulus = self.prime_pow.get_modulus(cap_div(self.relprec, self.prime_pow.e))[0]
            ZZ_pX_construct(&high_shifter)
            ZZ_p_construct(&tmp_zp)
            ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
            ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
            if val >= self.prime_pow.prec_cap:
                # high_shifter = p^(2^(high_length - 1))/x^(e*2^(high_length - 1))
                ZZ_pX_conv_modulus(high_shifter, self.prime_pow.high_shifter[self.prime_pow.high_length-1].val(), c.x)
                # if val = r + s * 2^(high_length - 1)
                # then high_shifter = p^(s*2^(high_length - 1))/x^(e*s*2^(high_length - 1))
                ZZ_pX_PowerMod_long_pre(high_shifter, high_shifter, (val / (1 << (self.prime_pow.high_length - 1))), modulus)
                ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus)
                # Now we only need to multiply self.unit by p^r/x^(e*r) where r < 2^(high_length - 1), which is tractible.
                val = val % (1 << (self.prime_pow.high_length - 1))
            while val > 0:
                if val & 1:
                    ZZ_pX_conv_modulus(high_shifter, self.prime_pow.high_shifter[i].val(), c.x)
                    ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus)
                val = val >> 1
                i += 1
            self._normalized = 1
        ZZ_p_destruct(&tmp_zp)
        ZZ_destruct(&num_zz)
        ZZ_destruct(&den_zz)
        mpz_clear(tmp_m)
        mpz_clear(num_unit)
        mpz_clear(den_unit)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, unsigned long relprec) except -1:
        """
        self.prime_pow must already be set.
        relprec should be in range 0 <= relprec <= self.prime_pow.ram_prec_cap
        """
        if ZZX_IsZero(poly):
            self._set_inexact_zero(absprec)
            return
        cdef long i = 0
        cdef long deg = ZZX_deg(poly)
        cdef long mini = -1
        cdef long minval
        cdef long curval
        cdef ZZ_c tmp_z
        cdef ZZ_c *ppow
        cdef ZZ_pX tmp_p
        ZZ_construct(&tmp_z)
        while mini == -1:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                minval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                mini = i
            i += 1
        while i <= deg:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                curval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                if curval < minval:
                    minval = curval
                    mini = i
            i += 1
        if self.prime_pow.e == 1: # unramified
            self.ordp = minval
            self.relprec = absprec - self.ordp
            if relprec < self.relprec:
                self.relprec = relprec
            elif self.relprec < 0:
                self.set_inexact_zero(absprec)
                return
            self.prime_pow.restore_context(relprec)
            ppow = self.prime_pow.pow_ZZ_tmp(minval)
            for i from 0 <= i <= deg:
                ZZ_div(tmp_z, ZZX_coeff(poly, i), ppow[0])
                ZZ_pX_SetCoeff(self.unit, i, ZZ_to_ZZ_p(tmp_z))
        else: # eisenstein
            self.ordp = minval * self.prime_pow.e + mini
            self.relprec = absprec - self.ordp
            if relprec < self.relprec:
                self.relprec = relprec
            elif self.relprec < 0:
                self.set_inexact_zero(absprec)
                return
            self.prime_pow.restore_context(minval + cap_div(relprec + mini, self.prime_pow.e))
            ZZ_pX_construct(&tmp_p)
            ZZX_to_ZZ_pX(tmp_p, poly)
            if PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                ZZ_pX_eis_shift(self.unit, tmp_p, self.ordp, (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).low_shifter, \
                                (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter, \
                                self.prime_pow.get_modulus(cap_div(relprec, self.prime_pow.e)), \
                                self.prime_pow.pow_ZZ_tmp(1)[0], \
                                self.prime_pow.get_context(cap_div(relprec, self.prime_pow.e)))
            elif PY_TYPE_CHECK(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                ZZ_pX_eis_shift(self.unit, tmp_p, self.ordp, (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).low_shifter, \
                                (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter, \
                                self.prime_pow.get_modulus(cap_div(relprec, self.prime_pow.e)), \
                                self.prime_pow.pow_ZZ_tmp(1)[0], \
                                self.prime_pow.get_context(cap_div(relprec, self.prime_pow.e)))
            else:
                raise RuntimeError, "prime_pow is of inconsistent type"
            ZZ_pX_destruct(&tmp_p)
        ZZ_destruct(&tmp_z)
        self._normalized = 1


    cdef int _set_prec_rel(self, unsigned long relprec) except -1:
        self.relprec = relprec

    cdef int _set_prec_both(self, long absprec, unsigned long relprec) except -1:
        self.relprec = absprec - self.ordp
        if relprec < self.relprec:
            self.relprec = relprec
        elif self.relprec <= 0:
            self._set_inexact_zero(absprec)

    cdef int _set_absprec_rel(self, unsigned long relprec) except -1:
        self.relprec = self.ordp + relprec

    cdef int _set_absprec_both(self, long absprec, unsigned long relprec) except -1:
        self.relprec = self.ordp + relprec
        if absprec < aprec:
            self.relprec = absprec
            if self.relprec < val:
                self._set_inexact_zero(absprec)

    cdef int _normalize(self) except -1:
        cdef long minval, mini, shift
        if self._normalized == 0:
            if self._is_exact_zero(): # this case should never happen
                #print "exact zero somehow"
                pass
            elif ZZ_pX_IsZero(self.unit):
                self.ordp += self.relprec
                self.relprec = 0
            else:
                ZZ_pX_min_val_coeff(minval, mini, self.unit, self.prime_pow.pow_ZZ_tmp(1)[0])
                shift = minval * self.prime_pow.e + mini
                self.relprec -= shift
                self.ordp += shift
                ZZ_pX_eis_shift(self.unit, self.unit, shift, self.prime_pow.low_shifter, self.prime_pow.high_shifter, self.prime_pow.get_modulus(cap_div(self.relprec, self.prime_pow.e)), self.prime_pow.pow_ZZ_tmp(1)[0], self.prime_pow.get_context(cap_div(self.relprec, self.prime_pow.e)))
            self._normalized = 1

cdef long cap_div(long a, long e):
    if e == 1:
        return a
    elif a % e == 0:
        return a / e
    else:
        return a / e + 1
