import sage.rings.integer_mod
import padic_base_generic_element
import sys

Mod = sage.rings.integer_mod.Mod
Integer = sage.rings.integer.Integer
#pAdicLazyGenericElement = sage.rings.padics.padic_lazy_generic_element.pAdicLazyGenericElement
infinity = sage.rings.infinity.infinity
HaltingError = sage.rings.padics.precision_error.HaltingError
PrecisionLimitError = sage.rings.padics.precision_error.PrecisionLimitError
PrecisionError = sage.rings.padics.precision_error.PrecisionError
pAdicBaseGenericElement = padic_base_generic_element.pAdicBaseGenericElement

class pAdicLazyElement(pAdicBaseGenericElement):
    def _add_(self, right):
        if isinstance(self, pAdicLazy_zero):
            return right
        if isinstance(right, pAdicLazy_zero):
            return self
        return pAdicLazy_add(self, right)

    def __cmp__(self, right, halt = None):
        if halt is None:
            halt = self.parent().halting_parameter()
        #print "self = %s, type = %s"%(self, type(self))
        #print "right = %s, type = %s"%(right, type(right))
        if self.valuation() < right.valuation():
            return -1
        elif self.valuation() > right.valuation():
            return 1
        elif self.valuation() is infinity and right.valuation() is infinity:
            return 0
        #comparing valuations can throw an exception if self and right are both congruent to zero modulo p^halt, but we want this to be passed upstream.  We now know that the valuations are equal and both less than halt, so both have type Integer now rather than Valuation.
        if self._cache.lift() % self.parent().prime_pow(halt - self.valuation()) < right._cache.lift() % self.parent().prime_pow(halt - right.valuation()):
            return -1
        elif self._cache.lift() % self.parent().prime_pow(halt - self.valuation()) > right._cache.lift() % self.parent().prime_pow(halt - right.valuation()):
            return 1
        if self._cache != right._cache:
            #According to spec, we only compare up to halt.  We now know that self and right are equal up to halt, but unequal above.  In order to maintain consistency with previous calls to __cmp__, we throw an exception here.
            raise PrecisionError, "I'm supposed to tell you these are equal, because they're congruent modulo p^halt, but "
        jump = 1
        start = min(self.precision_absolute(), right.precision_absolute())
        end = halt - start
        while jump <= end:
            self.set_precision_absolute(start + jump)
            right.set_precision_absolute(start + jump)
            if self._cache != right._cache:
                if self._cache.lift() % self.parent().prime_pow(start + jump) < right._cache.lift() % self.parent().prime_pow(start + jump):
                    return -1
                elif self._cache.lift() % self.parent().prime_pow(start + jump) < right._cache.lift() % self.parent().prime_pow(start + jump):
                    return 1
            jump = 2 * jump
        self.set_precision_absolute(halt)
        right.set_precision_absolute(halt)
        if self._cache != right._cache:
            if self._cache.lift() % self.parent().prime_pow(start + jump) < right._cache.lift() % self.parent().prime_pow(start + jump):
                return -1
            elif self._cache.lift() % self.parent().prime_pow(start + jump) < right._cache.lift() % self.parent().prime_pow(start + jump):
                return 1
        return 0

    def _div_(self, right):
        if isinstance(right, pAdicLazy_zero):
            raise ZeroDivisionError, "Cannot divide by zero"
        if isinstance(self, pAdicLazy_zero):
            return self
        return pAdicLazy_div(self, right)

    def __floordiv__(self, right):
        right = self.parent()(right)
        if self.parent().is_field():
            return self._div_(right)
        if isinstance(right, pAdicLazy_zero):
            raise ZeroDivisionError, "Cannot divide by zero"
        if isinstance(self, pAdicLazy_zero):
            return self
        try:
            right.set_precision_relative(1)
        except PrecisionError:
            raise ZeroDivisionError, "Cannot divide by a p-adic very close to zero.  Casting both into a lazy ring with higher halting parameter may help"
        return pAdicLazy_rshift(self / right.unit_part(), right.valuation())

    def __getitem__(self, n):
        if isinstance(n, slice):
            if n.start == 0 and self.valuation() < 0:
                raise ValueError, "Because of a limitation in Python, you need to use the function slice() or the syntax [::] instead of the syntax [:]"
            if n.start == 0:
                i = None
            else:
                i = n.start
            if n.start == sys.MAXINT:
                j = None
            else:
                j = n.start
            if n.step == None:
                k = 1
            else:
                k = n.step
            return pAdic_slice(self, i, j, k)
        self.set_precision_absolute(n + 1)
        ppow = self.parent().prime_pow(n - self._base_valuation)
        return (self._cache.lift() // ppow) % self.parent().prime()

    def __lshift__(self, shift):
        shift = Integer(shift)
        if shift < 0 and not self.parent().is_field():
            return pAdicLazy_rshift(self, -shift)
        return pAdicLazy_lshift(self, shift)

    def __rshift__(self, shift):
        shift = Integer(shift)
        if shift < 0 or self.parent().is_field():
            return pAdicLazy_lshift(self, -shift)
        return pAdicLazy_rshift(self, shift)

    def __mod__(self, right):
        right = self.parent()(right)
        if self.parent().is_field() or self.valuation() > right.valuation() or isinstance(self, pAdicLazy_zero):
            return pAdicLazy_zero(self.parent())
        #return self - right * (self // right) #can be made more efficient
        try:
            right.set_precision_relative(1)
        except PrecisionError:
            raise ZeroDivisionError, "Cannot divide by a p-adic very close to zero.  Casting both into a lazy ring with higher halting parameter may help"
        rval = right.valuation()
        self.set_precision_absolute(rval)
        right.set_precision_relative(rval - self.valuation())
        before = self._cache / right._cache
        intreturn = before.lift() % self.parent().prime_pow(right.valuation())
        if intreturn == 0:
            return pAdicLazy_zero(self.parent())
        return pAdicLazy_integer(self.parent(), intreturn)

    def _mul_(self, right):
        if isinstance(self, pAdicLazy_zero):
            return self
        if isinstance(right, pAdicLazy_zero):
            return right
        return pAdicLazy_mul(self, right)

    def _neg_(self):
        if isinstance(self, pAdicLazy_zero):
            return self
        return pAdicLazy_neg(self)

    def __pow__(self, right):
        #Add zero support
        return pAdicLazy_pow(self, right)

    def _repr_(self):
        return self._repr()

    def _repr(self, mode = None, do_latex = False):
        self._recompute()
        return pAdicBaseGenericElement._repr(self, mode, do_latex)

    def _sub_(self, right):
        if isinstance(right, pAdicLazy_zero):
            return self
        if isinstance(self, pAdicLazy_zero):
            return pAdicLazy_neg(right)
        return pAdicLazy_sub(self, right)

    def _set_cache(self, cache):
        self._cache = cache

    def _set_cache_prec(self, cache_prec):
        self._cache_prec = cache_prec

    def _set_base_valuation(self, base_valuation):
        self._base_valuation = base_valuation

    def exp(self):
        if (self.valuation() > 0 and self.parent().prime() != 2) or self.valuation() > 1:
            return pAdicLazy_exp(self)
        raise ValueError, "cannot take exp of a unit (or x == 2 (mod 4))"

    def exp_artin_hasse(self):
        #Do valuation checking
        return pAdicLazy_expah(self)

    def is_square(self):
        if isinstance(self, pAdicLazy_zero):
            return True
        if self.parent().prime() != 2:
            self.set_precision_relative(1)
            return (self.valuation() % 2 == 0) and Mod(self._cache, self.parent().prime()).is_square()
        self.set_precision_relative(3)
        return (self.valuation() % 2 == 0) and (Mod(self._cache, 8) == 1)

    def is_equal_to(self, right, prec):
        if self.precision_relative() > 0 and right.precision_relative() > 0:
            if self.valuation() != right.valuation():
                return False
            if self._cache != right._cache:
                return False
        elif self.precision_relative() > 0:
            if self.valuation() < right._min_valuation():
                return False
        elif right.precision_relative() > 0:
            if right.valuation() < self._min_valuation():
                return False
        n = min(self._min_valuation(), right._min_valuation())
        while n <= prec:
            self.set_precision_absolute(n)
            right.set_precision_absolute(n)
            if (self._cache != right._cache) or self._cache_prec != right._cache_prec:
                return False
        return True

    def is_zero(self, prec):
        return self.is_equal_to(pAdicLazy_zero(self.parent()), prec)

    def lift(self):
        return self.parent().prime_pow(self._base_valuation) * self._cache.lift()

    def lift_to_precision(self, absprec):
        self.set_precision_absolute(absprec)
        return copy.copy(self)

    def list(self):
        if isinstance(self, pAdicLazy_zero):
            return []
        def plist(n, p, prec):
            if prec == 0:
                return []
            else:
                return [n % p] + plist(n // p, p, prec - 1)
        if self.parent().is_field():
            return plist(self._cache, self.parent().prime(), self._cache_prec)
        return ([0] * self._base_valuation) + plist(self._cache, self.parent().prime())

    def log(self):
        if not self.is_unit():
            raise ValueError, "cannot take log of a non-unit"
        return pAdicLazy_log(self)

    def log_artin_hasse(self):
        #Do valuation checking
        return pAdic_logah(self)

    def padded_list(self, absprec = infinity, relprec = infinity):
        if absprec is infinity and relprec is infinity:
            raise ValueError, "must specify at least one of absprec and relprec"
        if self.parent().is_field():
            if self.valuation() is infinity:
                if relprec < infinity:
                    return [self.parent().residue_class_field()(0)]*relprec
                else:
                    return []
            relprec = min(relprec, absprec - self.valuation())
            return self.list()[:relprec] + [self.parent().residue_class_field()(0)]*(self.precision_relative() - relprec)
        if self.valuation() is infinity:
            if absprec < infinity:
                return [self.parent().residue_class_field()(0)]*absprec
            else:
                raise ValueError, "would return infinite list"
        absprec = min(absprec, relprec + self.valuation())
        return self.list()[:absprec] + [self.parent().residue_class_field()(0)]*(self.precision_absolute() - absprec)


    def precision_absolute(self):
        return self._cache_prec + self._base_valuation

    def precision_relative(self):
        return self._cache_prec

    ##Rational Reconstruction?

    def _recompute(self):
        pass

    def residue(self, n):
        self.set_precision_absolute(n)
        try:
            return Mod(self._cache.lift() * self.parent().prime_pow(self._base_valuation), self.parent().prime_pow(n))
        except ZeroDivisionError:
            raise ValueError, "element must have non-negative valuation in order to compute residue"

    def slice(self, i, j, k = 1): #despite the fact that slice is a keyword, this should still work
        return pAdicLazy_slice(self, i, j, k)

    def sqrt(self):
        if isinstance(self, pAdicLazy_zero):
            return self
        return pAdicLazy_sqrt(self)

    def square_root(self):
        if isinstance(self, pAdicLazy_zero):
            return self
        return pAdicLazy_sqrt(self)


    def _unit_part(self):
        if self.precision_relative() > 0:
            return self._cache
        return Mod(0,1)

    def unit_part(self):
        if isinstance(self, pAdicLazy_zero):
            raise ValueError, "Zero does not have a well defined unit part"
        return pAdicLazy_unitpart(self)

    def valuation(self, halt = None):
        if isinstance(self, pAdicLazy_zero):
            return infinity
        if self._cache != 0:
            return self._base_valuation
        if halt is None:
            halt = self.parent().halting_parameter()
        #The following is a stopgap measure until valuation.py is working.
        self.set_precision_relative(1)
        return self.valuation()

    def _min_valuation(self):
        if self._cache != 0:
            return self._base_valuation
        return self._base_valuation + self._cache_prec

    def _is_exact_zero(self):
        return False

# The following subclasses are used to create pAdicLazyElements from other data
# They need to implement the following methods:
# set_precision_absolute, set_precision_relative, __init__, (copy)

class pAdicLazy_integer(pAdicLazyElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity): #cannot call with input zero
        """
        sage: a = ZpL(5)(-3); loads(dumps(a)) == a
        True
        """
        pAdicLazyElement.__init__(self, parent)
        self._set_base_valuation(x.valuation(parent.prime()))
        self._x = x // parent.prime_pow(self._base_valuation)
        prec = min(parent.precision_cap(), relprec, absprec - self._base_valuation)
        self._set_cache(Mod(self._x, parent.prime_pow(prec)))
        self._set_cache_prec(prec)

    def set_precision_relative(self, n, halt = None):
        if n > self._cache_prec:
            self._set_cache_prec(n)
            self._set_cache(Mod(self._x, self.parent().prime_pow(n)))

    def set_precision_absolute(self, n, halt = None):
        if n > self._base_valuation + self._cache_prec:
            self._set_cache_prec(n - self._base_valuation)
            self._set_cache(Mod(self._x, self.parent().prime_pow(self._cache_prec)))


class pAdicLazy_rational(pAdicLazyElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        self._set_base_valuation(x.valuation(parent.prime()))
        self._x = x / parent.prime_pow(self._base_valuation)
        prec = min(parent.precision_cap(), relprec, absprec - self._base_valuation)
        self._set_cache(Mod(self._x, parent.prime_pow(prec)))
        self._set_cache_prec(prec)

    def set_precision_relative(self, n, halt = None):
        if n > self._cache_prec:
            self._set_cache_prec(n)
            self._set_cache(Mod(self._x, self.parent().prime_pow(n)))

    def set_precision_absolute(self, n, halt = None):
        if n > self._base_valuation + self._cache_prec:
            self._set_cache_prec(n - self._base_valuation)
            self._set_cache(Mod(self._x, self.parent().prime_pow(self._cache_prec)))


class pAdicLazy_otherpadic(pAdicLazyElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        self._set_base_valuation(x.valuation())
        self._x = x
        prec = min(parent.precision_cap(), relprec, absprec - self._base_valuation, x.precision_relative())
        self._set_cache_prec(prec)
        self._set_cache(Mod(x.unit_part().lift(), self.parent().prime_pow(prec)))

    def set_precision_relative(self, n, halt = None):
        if n > self._cache_prec:
            if n > x.precision_relative():
                raise PrecisionLimitError, "Cannot compute more p-adic digits"
            self._set_cache(Mod(x.unit_part().lift(), self.parent().prime_pow(n)))
            self._set_cache_prec(n)

    def set_precision_absolute(self, n, halt = None):
        if n > self._cache_prec + self._base_valuation:
            if n > self._x.precision_absolute():
                raise PrecisionLimitError, "Cannot compute more p-adic digits"
            self._set_cache_prec(n - self._base_valuation)
            self._set_cache(Mod(x.unit_part().lift(), self.parent().prime_pow(self._cache_prec)))

class pAdicLazy_mod(pAdicLazyElement):
    def __init__(self, parent, x, ppow):
        raise NotImplementedError

class pAdicLazy_zero(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._set_base_valuation(infinity)
        self._set_cache(Mod(0,1))
        self._set_cache_prec(0)

    def _repr(self, mode = None, do_latex = False):
        return "0"

    def set_precision_relative(self, n, halt = None):
        if n > 0:
            raise PrecisionLimitError, "Cannot set the relative precision of 0 to a positive value"

    def set_precision_absolute(self, n, halt = None):
        pass

    def _is_exact_zero(self):
        return True

class pAdicLazy_valpower(pAdicLazyElement):
    def __init__(self, parent, v):
        pAdicLazyElement.__init__(self, parent)
        raise NotImplementedError

class pAdicLazy_iterator(pAdicLazyElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        raise NotImplementedError

class pAdicLazy_teichmuller(pAdicLazyElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        p = parent.prime()
        prec = min(parent.precision_cap(), relprec, absprec)
        x = Mod(x,parent.prime_pow(prec))
        u = 1/Mod(1-p,parent.prime_pow(prec))
        delta = u*(1 - x**(p-1))
        xnew = x - x*delta*(1-p*delta)
        while x != xnew:
            x = xnew
            delta = u*(1 - x**(p-1))
            xnew = x - x*delta*(1-p*delta)
        self._set_cache(x)
        self._set_cache_prec(prec)
        self._set_base_valuation(Integer(0))

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            p = self.parent().prime()
            x = Mod(self._cache, self.parent().prime_pow(n))
            u = 1 / Mod(1 - p, self.parent().prime_pow(n))
            delta = u * (1 - x ** (p - 1))
            xnew = x - x*delta*(1 - p * delta)
            while x != xnew:
                x = xnew
                delta = u*(1-x**(p-1))
                xnew = x - x*delta*(1-p*delta)
            self._set_cache(xnew)
            self._set_cache_prec(n)

    def set_precision_absolute(self, n, halt = None):
        return self.set_precision_relative(n, halt)

class pAdicLazy_root(pAdicLazyElement):
    def __init__(self, parent, f, x, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        self._f = f
        self._x = x
        raise NotImplementedError

class pAdicLazy_integral(pAdicLazyElement):
    #I'm not sure exactly how this will work, but it's another way to get p-adics
    def __init__(self, parent, f, X, absprec=infinity, relprec=infinity):
        pAdicLazyElement.__init__(self, parent)
        raise NotImplementedError


# The following subclasses are used to create pAdicLazyElements from other pAdicLazyElements
#As above, they need to implement set_precision_relative , set_precision_absolute, __init__, (copy)

class pAdicLazy_bintype(pAdicLazyElement):
    def __init__(self, x, y, op):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._y = y
        self._op = op
        raise NotImplementedError


class pAdicLazy_addtype(pAdicLazyElement):
    def __init__(self, x, y): #neither x nor y should be exact zero
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._y = y
        self._recompute()

    def set_precision_absolute(self, n, halt = None):
        if n > self._cache_prec + self._base_valuation:
            self._x.set_precision_absolute(n, halt)
            self._y.set_precision_absolute(n, halt)
            self._recompute()

    def set_precision_relative(self, n, halt = None):
        if n > self._cache_prec:
            if self.precision_relative() > 0:
                self._x.set_precision_absolute(n + self._base_valuation, halt)
                self._y.set_precision_absolute(n + self._base_valuation, halt)
                self._recompute()
            else:
                if halt is None:
                    halt = self.parent().halting_parameter()
                while self._base_valuation < halt:
                    self._x.set_precision_absolute(self._base_valuation + n, halt)
                    self._y.set_precision_absolute(self._base_valuation + n, halt)
                    self._recompute()
                    if self.precision_relative() > 0:
                        self._x.set_precision_absolute(self._base_valuation + n, halt)
                        self._y.set_precision_absolute(self._base_valuation + n, halt)
                        self._recompute()
                        break
                if self.precision_relative() == 0:
                    raise HaltingError, "Stopped computing sum: set halting parameter higher if you want computation to continue"


class pAdicLazy_add(pAdicLazy_addtype):
    def _recompute(self):
        mv = min(self._x._base_valuation, self._y._base_valuation)
        value = self._x._cache.lift() * self.parent().prime_pow(self._x._base_valuation - mv) + self._y._cache.lift() * self.parent().prime_pow(self._y._base_valuation - mv)
        valuation = value.valuation(self.parent().prime())
        relprec = min(self._x._cache_prec + self._x._base_valuation, self._y._cache_prec + self._y._base_valuation) - valuation - mv
        if relprec <= 0:
            self._set_base_valuation(min(self._x._cache_prec - self._x._base_valuation, self._y._cache_prec - self._y._base_valuation))
            self._set_cache(Mod(0, 1))
            self._set_cache_prec(0)
        else:
            self._set_base_valuation(valuation + mv)
            self._set_cache_prec(relprec)
            self._set_cache(Mod(value // self.parent().prime_pow(valuation), self.parent().prime_pow(self._cache_prec)))


class pAdicLazy_sub(pAdicLazy_addtype):
    def _recompute(self):
        mv = min(self._x._base_valuation, self._y._base_valuation)
        value = self._x._cache.lift() * self.parent().prime_pow(self._x._base_valuation - mv) - self._y._cache.lift() * self.parent().prime_pow(self._y._base_valuation - mv)
        valuation = value.valuation(self.parent().prime())
        relprec = min(self._x._cache_prec + self._x._base_valuation, self._y._cache_prec + self._y._base_valuation) - valuation - mv
        if relprec <= 0:
            self._set_base_valuation(min(self._x._cache_prec - self._x._base_valuation, self._y._cache_prec - self._y._base_valuation))
            self._set_cache(Mod(0, 1))
            self._set_cache_prec(0)
        else:
            self._set_base_valuation(valuation + mv)
            self._set_cache_prec(relprec)
            self._set_cache(Mod(value // self.parent().prime_pow(valuation), self.parent().prime_pow(self._cache_prec)))

class pAdicLazy_multype(pAdicLazyElement):
    def __init__(self, x, y):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._y = y
        self._recompute()

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n, halt)
            self._y.set_precision_relative(n, halt)
            self._recompute()

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            if self.precision_relative() == 0: #We are currently indistinuishable from zero and thus must first try to make x and y separate from zero.
                if self._x.precision_relative() == 0 and self._y.precision_relative() == 0:
                    halfdiff = (n - self._base_valuation) // 2
                    try:
                        self._x.set_precision_absolute(self._x.precision_absolute() + halfdiff, halt)
                    except PrecisionError:
                        pass
                    try:
                        self._y.set_precision_absolute(n - self._x._min_valuation(), halt)
                    except PrecisionError:
                        self._x.set_precision_absolute(n - self._y._min_valuation(), halt)
                elif self._x.precision_relative() == 0:
                    self._x.precision_absolute(n - self._y._base_valuation, halt)
                elif self._y.precision_relative() == 0:
                    self._y.precision_absolute(n - self._x._base_valuation, halt)
                #Now that we have at least n valuation, either we're zero in which case we're done, or we have to
                if self._x.precision_relative() > 0 and self._y.precision_relative() > 0:
                    self._recompute()
                    self._x.set_precision_relative(n - self._base_valuation, halt)
                    self._y.set_precision_relative(n - self._base_valuation, halt)
                    self._recompute()
                else:
                    self._recompute()
            else:
                self.set_precision_relative(n - self._base_valuation, halt)


class pAdicLazy_mul(pAdicLazy_multype):
    def _recompute(self):
        self._set_base_valuation(self._x._base_valuation + self._y._base_valuation)
        self._set_cache_prec(min(self._x._cache_prec, self._y._cache_prec))
        self._set_cache(self._x._cache * self._y._cache)

class pAdicLazy_div(pAdicLazy_multype):
    def __init__(self, x, y, halt = None):
        pAdicLazyElement.__init__(self, x.parent().fraction_field())
        self._x = x
        self._y = y
        try:
            y.set_precision_relative(1, halt)
        except PrecisionError:
            raise ZeroDivisionError, "Cannot divide by a p-adic very close to zero.  Casting both into a lazy ring/field with higher halting parameter may help"
        self._recompute()

    def _recompute(self):
        self._set_base_valuation(self._x._base_valuation - self._y._base_valuation)
        self._set_cache_prec(min(self._x._cache_prec, self._y._cache_prec))
        self._set_cache(self._x._cache / self._y._cache)

class pAdicLazy_lshift(pAdicLazyElement):
    def __init__(self, x, shift): #self can be a ring or field element, shift positive or negative integer.  If a ring element, shift is positive
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._shift = shift
        self._recompute()

    def _recompute(self):
        self._set_base_valuation(self._x._base_valuation + self._shift)
        self._set_cache_prec(self._x._cache_prec)
        self._set_cache(self._x._cache)

    def set_precision_relative(self, n, halt = None):
        self._x.set_precision_relative(n, halt)
        self._recompute()

    def set_precision_absolute(self, n, halt = None):
        self._x.set_precision_absolute(n + self._shift, halt)

class pAdicLazy_rshift(pAdicLazyElement):
    def __init__(self, x, shift): #self is a ring element and shift is nonnegative
        pAdicLazyElement.__init__(self, x.parent())
        x.set_precision_relative(shift + 1) #I'm being a bit lazy here, so that the code in recompute is simpler
        self._x = x
        self._shift = shift
        self._recompute()

    def _recompute(self):
        val = self._x.valuation()
        if self._shift <= val:
            self._set_base_valuation(val - self._shift)
            self._set_cache(self._x._cache)
            self._set_cache_prec(self._x._cache_prec)
            self._suck = Integer(0)
        else:
            shift = self._shift
            unit = self._x._unit_part().lift() // self.parent().prime_pow(shift - val)
            val = unit.valuation(self.parent().prime())
            rprec = self._x.precision_absolute() - shift
            if val is infinity:
                self._suck = self._x.precision_relative() #this might not be right
                val = max(rprec, Integer(0))
                rprec = Integer(0)
                unit = Mod(0, 1)
            elif val > 0:
                self._suck = shift - self._x.valuation() + val
                rprec = rprec - val
                unit = Mod(unit // self.parent().prime_pow(val), self.parent().prime_pow(rprec))
            else:
                self._suck = shift - self._x.valuation()
                unit = Mod(unit, self.parent().prime_pow(rprec))
            self._set_base_valuation(val)
            self._set_cache(unit)
            self._set_cache_prec(rprec)

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n + self._suck, halt)
            self._recompute()

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            self._x.set_precision_absolute(n + self._shift, halt)

class pAdicLazy_pow(pAdicLazyElement):
    def __init__(self, x, y):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._y = y
        raise NotImplementedError

class pAdicLazy_uni(pAdicLazyElement):
    def __init__(self, x, op):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._op = op
        raise NotImplementedError

class pAdicLazy_neg(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._recompute()

    def _recompute(self):
        self._set_base_valuation(x._base_valuation)
        self._set_cache_prec(x._cache_prec)
        self._set_cache(-x._cache)

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            self._x.set_precision_absolute(n, halt)
            self._recompute()

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n, halt)
            self._recompute()

class pAdicLazy_invert(pAdicLazyElement):
    def __init__(self, x, halt = None):
        pAdicLazyElement.__init__(self, x.parent().fraction_field())
        self._x = x
        try:
            x.set_precision_relative(1, halt)
        except PrecisionError:
            raise ZeroDivisionError, "Cannot divide by a p-adic very close to zero.  Casting both into a lazy ring/field with higher halting parameter may help"
        self._recompute()

    def _recompute():
        self._set_base_valuation(-self._x._base_valuation)
        self._set_cache_prec(self._x._cache_prec)
        self._set_cache(self._x._cache.__invert__())

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n, halt)
            self._recompute()

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            self._x.set_precision_relative(n - self._base_valuation, halt)
            self._recompute()

class pAdicLazy_slice(pAdicLazyElement):
    def __init__(self, x, i, j, k):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        if not (i is None):
            self._i = i
        else:
            self._i = infinity
        if not (j is None):
            self._j = j
        else:
            self._j = x._min_valuation()
        self._recompute()

    def _recompute(self):
        pass

    #def set_precision_absolute(self, n, halt = None):
    #    if n > self.precision_absolute():
    #        if n > self._i:
    #            self._set_cache_prec(n)
    #            self._set_cache(Mod(self._cache, self.parent().prime_pow(self._cache_prec)))

class pAdicLazy_unitpart(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent().integer_ring())
        self._x = x
        self._set_base_valuation(Integer(0))
        self._recompute()

    def _recompute(self):
        self._set_cache(self._x._cache)
        self._set_cache_prec(self._x._cache_prec)

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n)
            self._recompute()

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_relative(n)
            self._recompute()

    def valuation(self):
        return Integer(0)


class pAdicLazy_log(pAdicLazyElement):
    def __init__(self, x): #x must be a unit
        pAdicLazyElement.__init__(self, x.parent())
        if x.residue(1) == 1:
            self._x = x
            self._n = Integer(1)
        else:
            self._n = x.residue(1).multiplicative_order()
            self._x = x.__pow__(self._n)
            self._n = x.parent().prime() - 1
        self._recompute()

    def _recompute(self):
        p = self.parent().prime()
        prec = self._x.precision_absolute()
        extra_prec = 0
        while extra_prec < Integer(prec + extra_prec).exact_log(p):
            extra_prec += 1
        from sage.rings.padics.zp import Zp
        working_ring = Zp(p, prec + extra_prec, type = 'capped-abs')
        x = working_ring(Integer(1) - self._x.residue(prec).lift())
        xpow = x
        ans = working_ring(0)
        for n in range(1, prec + extra_prec):
            ans -= xpow // working_ring(n)
            xpow *= x
        ans = ans // self._n
        ans = ans.residue(prec)
        val = ans.lift().valuation(self.parent().prime())
        ans = Mod(ans.lift() // self.parent().prime_pow(val), self.parent().prime_pow(prec - val))
        self._set_cache(ans)
        self._set_base_valuation(val)
        self._set_cache_prec(prec - self._base_valuation)

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            self._x.set_precision_absolute(n, halt)
            self._recompute()

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            v = self.valuation()
            self._x.set_precision_absolute(v + n, halt)
            self._recompute()

class pAdicLazy_exp(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        self._set_base_valuation(Integer(0))
        self._recompute()

    def _recompute(self):
        if self._x.precision_relative() == 0:
            if self.parent().prime() == 2:
                self._set_cache_prec(self._x.precision_absolute() - 1)
            else:
                self._set_cache_prec(self._x.precision_absolute())
            self._set_cache(Mod(1, self.parent().prime_pow(self._cache_prec)))
        else:
            val = self._x.valuation()
            p = self.parent().prime()
            prec = self._x.precision_absolute()
            max_term = ((p-1)*(prec-1)) // ((p-1) * val - 1) + 1
            extra_prec = max_term // (p - 1)
            from sage.rings.padics.zp import Zp
            working_ring = Zp(p, prec + extra_prec, type = 'capped-abs')
            x = working_ring(self._x.lift())
            term = ans = working_ring(1)
            for n in range(1, max_term):
                term *= x
                term //= working_ring(n)
                ans += term
            self._set_cache_prec(prec)
            self._set_cache(ans.residue(prec))

    def set_precision_absolute(self, n, halt = None):
        if n > self.precision_absolute():
            self._x.set_precision_absolute(n, halt) #need more when p == 2?
            self._recompute()

    def set_precision_relative(self, n, halt = None):
        if n > self.precision_relative():
            self._x.set_precision_absolute(n, halt) #need more when p == 2?
            self._recompute()

class pAdicLazy_logah(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        raise NotImplementedError

class pAdicLazy_expah(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        raise NotImplementedError

class pAdicLazy_sqrt(pAdicLazyElement):
    def __init__(self, x):
        pAdicLazyElement.__init__(self, x.parent())
        self._x = x
        raise NotImplementedError
