import sage.rings.commutative_ring_element
import sage.rings.extended_rational_field
import sage.rings.integer
import sage.rings.padics.precision_error
import sage.rings.infinity
import copy

Integer = sage.rings.integer.Integer
SignError = sage.rings.infinity.SignError
PrecisionError = sage.rings.padics.precision_error.PrecisionError
HaltingError = sage.rings.padics.precision_error.HaltingError
PrecisionLimitError = sage.rings.padics.precision_error.PrecisionLimitError
from sage.rings.infinity import infinity
infinity
minfinity = -infinity
CRE = sage.rings.commutative_ring_element.CommutativeRingElement

class Valuation(CRE):
    def _set_low(self, value):
        self._low = value

    def _set_high(self, value):
        self._high = value

    def _set_lowflex(self, value):
        self._lowflex = value

    def _set_highflex(self, value):
        self._highflex = value

    def _set_lowlimit(self, value):
        self._lowlimit = value

    def _set_highlimit(self, value):
        self._highlimit = value

    def _lowlimit(self):
        try:
            return self._lowlimit
        except AttributeError:
            return self._high

    def _highlimit(self):
        try:
            return self._highlimit
        except AttributeError:
            return self._low

    def _add_(self, right, halt = None):
        if isinstance(self, Valuation_infinity):
            if isinstance(right, Valuation_infinity) and right._low == -self._low:
                raise SignError, "cannot add infinity to minus infinity"
            if -self._low in right:
                right._fix(halt)
            return self
        elif isinstance(right, Valuation_infinity):
            if -right._low in self:
                self._fix(halt)
            return right
        elif self._low == minfinity and right._high == infinity or right._low == minfinity and self._high == infinity:
            _fix_one(self, right, halt)
        if self._low == self._high and right._low == right._high:
            return Valuation_point(self._low + right._low)
        return Valuation_add(self, right)

    def __cmp__(self):
        raise NotImplementedError

    def __contains__(self, x):
        return self._low <= x and x <= self._high

    def _div_(self, right):
        #Implementing division this way is an effort to reduce code complexity.
        return Valuation_mul(self, Valuation_invert(right))

    def __invert__(self):
        if 0 in self:
            if self._low == 0 or self._high == 0:
                if self._low == 0 and self._high == 0:
                    raise ZeroDivisionError, "cannot invert zero"
            else:
                if self._high > -self._low and self._lowflex:
                    try:
                        self._shrink_from_low_to(0)
                    except PrecisionError:
                        if self._low != 0 and self._highflex:
                            self._shrink_from_high_to(0)
                elif self._highflex:
                    try:
                        self._shrink_from_high_to(0)
                    except PrecisionError:
                        if self._high != 0 and self._lowflex:
                            self._shrink_from_low_to(0)
                elif self._lowflex:
                    try:
                        self._shrink_from_low_to(0)
                    except PrecisionError:
                        if self._low != 0 and self._highflex:
                            self._shrink_from_high_to(0)
                if 0 in self:
                    if self._low == 0 or self._high == 0:
                        if self._low == 0 and self._high == 0:
                            raise ZeroDivisionError, "cannot invert zero"
                    else:
                        #We tried being gentle, and it didn't work.  Now we hit self with the __cmp__ hammer.
                        if self == 0:
                            raise ZeroDivisionError, "cannot invert zero"
        if self._low == self._high:
            return Valuation_point(~self._low)
        return Valuation_invert(self)

    def _mul_(self, right, halt):
        if isinstance(self, Valuation_infinity):
            if right > 0:
                return self
            elif right < 0:
                return -self
            else:
                raise SignError, "cannot multiply infinity by zero"
        elif isinstance(right, Valuation_infinity):
            if self > 0:
                return right
            elif self < 0:
                return -right
            else:
                raise SignError, "cannot multiply infinity by zero"
        elif _is_unbounded(self) and 0 in right:
            if right == 0:
                self._fix(halt)
                return Valuation_point(Integer(0))
        elif _is_unbounded(right) and 0 in self:
            if self == 0:
                right._fix(halt)
                return Valuation_point(Integer(0))
        if self._low == self._high and right._low == right._high:
            return Valuation_point(self._low * right._low)
        return Valuation_mul(self, right)

    def _neg_(self):
        if isinstance(self, Valuation_infinity):
            return Valuation_infinity(-self._low)
        if self._low == self._high:
            return Valuation_point(self._low)
        return Valuation_neg(self)

    def __pow__(self, right):
        if not right:
            if isinstance(self, Valuation_infinity):
                raise ValueError, "infinity^0 not defined"
            elif not self:
                raise ArithmeticError, "0^0 is undefined."
            return Valuation_point(Integer(1))
        right = Integer(right)  #may later add support for Valuations, rationals and infinities in the exponent
        if isinstance(self, Valuation_infinity):
            if right > 0:
                return Valuation_infinity(self._low ** right)
            else: # right < 0
                return Valuation_point(self, Integer(0))
        if self._low == self._high:
            return Valuation_point(self, self._low ** right)
        return Valuation_pow(self, right)

    def _repr_(self, do_latex = False):
        self._recompute()
        if self._low == self._high:
            return "%s"%(self._low)
        return "[%s , %s]"%(self._low, self._high)

    def _sub_(self, right, halt):
        if isinstance(self, Valuation_infinity):
            if isinstance(right, Valuation_infinity) and right._low == self._low:
                raise SignError, "cannot add infinity to minus infinity"
            if self._low in right:
                right._fix(halt)
            return self
        elif isinstance(right, Valuation_infinity):
            if right._low in self:
                self._fix(halt)
            return -right
        elif self._low == minfinity and right._high == infinity or right._low == minfinity and self._high == infinity:
            _fix_one(self, right, halt)
        if self._low == self._high and right._low == right._high:
            return Valuation_point(self._low - right._low)
        return Valuation_sub(self, right)

    def _integer_(self, Z=None):
        self._fix()
        return Integer(self._low)

    def _rational_(self):
        self._fix()
        return Rational(self.low)

    def copy(self):
        return copy.copy(self)

    def lcm(self, other):
        if isinstance(self, Valuation_infinity):
            if other == 0:
                return Integer(0)
            else:
                return self.__abs__()
        elif isinstance(other, Valuation_infinity):
            if self == 0:
                return Integer(0)
            else:
                return other.__abs__()
        self._fix()
        other._fix()
        return lcm(Rational(self._low), Rational(other._low))

    def square_root(self):
        self._fix()
        return self._low.square_root()

    def sqrt(self):
        self._fix()
        return self._low.sqrt()

    def nth_root(self, n):
        self._fix()
        return self._low.nth_root(n)

    def __abs__(self):
        if isinstance(self, Valuation_infinity):
            return Valuation_infinity(infinity)
        if self._low == self._high:
            return Valuation_point(self._low)
        return Valuation_abs(self)

    def numerator(self):
        self._fix()
        return self._low.numerator()

    def numer(self):
        self._fix()
        return self._low.numer()

    def denominator(self):
        self._fix()
        return self._low.denominator()

    def denom(self):
        self._fix()
        return self._low.denom()

    def floor(self):
        if isinstance(self, Valuation_infinity):
            return self
        if self._low == self._high:
            return Valuation_point(self._low.floor())
        return Valuation_floor(self)

    def ceil(self):
        if isinstance(self, Valuation_infinity):
            return self
        if self._low == self._high:
            return Valuation_point(self._low.ceil())
        return Valuation_ceil(self)

class Valuation_infinity(Valuation):
    def __init__(self, value):
        self._set_low(value)
        self._set_high(value)
        self._set_lowflex(False)
        self._set_highflex(False)

    def _fix(self, halt = None):
        return

    def _recompute(self, chain = False):
        return

    def _recompute_low(self, chain = False):
        return

    def _recompute_high(self, chain = False):
        return

    def _flex(self):
        return False

class Valuation_point(Valuation):
    def __init__(self, parent, value):
        CRE.__init__(self, parent)
        self._set_low(value)
        self._set_high(value)
        self._set_lowflex(False)
        self._set_highflex(False)

    def _fix(self, halt = None):
        return

    def _recompute(self, chain = False):
        return

    def _recompute_low(self, chain = False):
        return

    def _recompute_high(self, chain = False):
        return

    def _flex(self):
        return False

class Valuation_padic(Valuation):
    def __init__(self, parent,  x):
        CRE.__init__(self, parent)
        self._x = x
        self._set_low(x._min_valuation())
        self._set_high(QQi.gen(1))
        self._set_lowflex(True)
        self._set_highflex(False)

    def _fix(self, halt = None):
        if halt is None:
            halt = 1 / (1 - self.parent().halting_parameter())
        self._x.set_precision_absolute(halt) #should switch to exponential backoff
        self._recompute()
        if self._low != self._high:
            raise HaltingError, "Cannot fix value; calling _fix with higher halt value may help"

    def _recompute(self, chain = False):
        self._set_low(self._x._min_valuation())
        if self._x._cache_prec > 0:
            self._set_high(self._x._min_valuation())
            self._set_lowflex(False)

    def _recompute_low(self, chain = False):
        self._set_low(self._x._min_valuation())
        if self._x._cache_prec > 0:
            self._set_lowflex(False)

    def _recompute_high(self, chain = False):
        if self._x._cache_prec > 0:
            self._set_high(self._x._min_valuation())

    def _flex(self):
        return self._lowflex

    def _budge(self, halt = None):
        error = False
        try:
            self._x._set_precision_absolute(self._x._min_valuation() + 1)
        except PrecisionLimitError:
            self._set_lowflex(False)
            error = True
            raise PrecisionLimitError
        finally:
            if not error:
                self._recompute()

    def _shrink_from_low_to(self, n):
        try:
            self._x._set_precision_absolute(n.floor() + 1) #should switch to exponential backoff
        except PrecisionLimitError:
            self._set_lowflex(False)
            raise PrecisionLimitError
        finally:
            self._recompute()

class Valuation_add(Valuation):
    def __init__(self, parent, x, y):
        CRE.__init__(self, parent)
        self._x = x
        self._y = y
        self._recompute()

    def _fix(self, halt = None):
        self._x.fix(halt)
        self._y.fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
            self._y._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = Valuation_point(self.parent(), self._x._low)
        if self._y._low == self._y._high:
            self._y = Valuation_point(self.parent(), self._y._low)

    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_low()
            self._y._recompute_low()
        self._set_low(self._x._low + self._y._low)
        self._set_lowflex(self._x._lowflex or self._y._lowflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_high()
            self._y._recompute_high()
        self._set_high(self._x._high + self._y._high)
        self._set_highflex(self._x._highflex or self._y._highflex)

    def _flex(self):
        return self._x._flex() or self._y._flex()

    def _budge(self, halt = None):
        """
        Preprocessing: self._flex() should return True.
        Postprocessing: ensure's self's state
        """
        if self._x._flex():
            try:
                self._x._budge(halt)
            except PrecisionLimitError:
                if self._y._flex():
                    self._y._budge(halt)
                else:
                    raise PrecisionLimitError
            finally:
                self._recompute(chain = False)
        else:
            try:
                self._y._budge(halt)
            finally:
                self._recompute(chain = False)


#This is the point I stopped writing.  There should be no compile issues, but this file is totally not ready.
    def _shrink_from_low_to(self, n, halt = None): #n must be finite, we assume lowflex and thus self._low finite, should call self._recompute_low() first
        """
        Preprocessing required: n must be finite, we assume that self._lowflex == True and thus self._low is finite.  Must call self._recompute_low(chain = True) before this function.  Must have self._low <= n < self._high
        Postprocessing: this function ensures that self's variables are updated correctly
        Errors: will raise a PrecisionLimitError if self cannot be found (ie cannot shrink one of the components from below) and self._low cannot be raised above n.  Note that this does not necessarily imply that self._lowflex = False.  self._high could have changed or there may be an internal limit point which would cause a halting error.
        """
        try:
            if self._x._lowflex and self._y._lowflex:
                if self._high is infinity:
                    if self._x._high is infinity and self._y._high is infinity:
                        try:
                            self._x._shrink_from_low_to(self._x._low + (n - self._low) / 2)
                        except PrecisionLimitError:
                            pass
                        try:
                            self._y._shrink_from_low_to(n - self._x._low)
                        except PrecisionLimitError:
                            if self._x._lowflex:
                                self._x._shrink_from_low_to(n - self._y._low)
                            else:
                                raise PrecisionLimitError
                    elif self._x._high is infinity:
                        try:
                            self._x._shrink_from_low_to(n - self._y._low)
                        except PrecisionLimitError:
                            pass
                        if n - self._x._low >= self._y._high:
                            raise PrecisionLimitError #this is the annoying case where we leave lowflex True but have to throw a PrecisionLimitError.  This will not be very well dealt with in the current version of __cmp__
                        else:
                            self._y._shrink_from_low_to(n - self._x._low)
        finally:
            self._recompute(chain = False)

    def _shrink_from_high_to(self, n):
        pass

class Valuation_sub(Valuation):
    def __init__(self, parent, x, y):
        CRE.__init__(self, parent)
        self._x = x
        self._y = y
        self._recompute()

    def _fix(self, halt = None):
        self._x.fix(halt)
        self._y.fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
            self._y._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)
        if self._y._low == self._y._high:
            self._y = self.parent()(self._y)


    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_low()
            self._y._recompute_high()
        self._set_low(self._x._low - self._y._high)
        self._set_lowflex(self._x._lowflex or self._y._highflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_high()
            self._y._recompute_low()
        self._set_high(self._x._high - self._y._low)
        self._set_highflex(self._x._highflex or self._y._lowflex)

    def _flex(self):
        return self._x._flex() or self._y._flex()

class Valuation_mul(Valuation):
    def __init__(self, parent, x, y):
        CRE.__init__(self, parent)
        self._x = x
        self._y = y
        self._recompute()

    def _fix(self, halt = None):
        self._x.fix(halt)
        self._y.fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
            self._y._recompute()
        self._recompute_signs()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)
        if self._y._low == self._y._high:
            self._y = self.parent()(self._y)


    def _recompute_signs(self):
        if self._x._low > 0:
            self._xsign = 1
        elif self._x._high < 0:
            self._xsign = -1
        else:
            self._xsign = 0
        if self._y._low > 0:
            self._ysign = 1
        elif self._y._high < 0:
            self._ysign = -1
        else:
            self._ysign = 0

    def _recompute_low(self, chain = True):
        if chain:
            xrecomputed = False
            if 0 in self.x:
                self._x._recompute()
                xrecomputed = True
            self._recompute_xsign()
            yrecomputed = False
            if 0 in self.y:
                self._y._recompute()
                yrecomputed = True
            self._recompute_ysign()
            if self._xsign == 1 and not yrecomputed:
                self._y._recompute_low()
            if self._ysign == 1 and not xrecomputed:
                self._x._recompute_low()
            if self._xsign == -1 and not yrecomputed:
                self._y._recompute_high()
            if self._ysign == -1 and not xrecomputed:
                self._x._recompute_high()
        if self._xsign == 1 and self._ysign == 1:
            self._set_low(self._x._low * self._y._low)
            self._set_lowflex(self._x._lowflex or self._y._lowflex)
        elif self._xsign == 1 and self._ysign == -1:
            self._set_low(self._x._high * self._y._low)
            self._set_lowflex(self._x._highflex or self._y._lowflex)
        elif self._xsign == -1 and self._ysign == 1:
            self._set_low(self._x._low * self._y._high)
            self._set_lowflex(self._x._lowflex or self._y._highflex)
        elif self._xsign == -1 and self._ysign == -1:
            self._set_low(self._x._high * self._y._high)
            self._set_lowflex(self._x._highflex or self._y._highflex)
        elif self._xsign == 0 and self._ysign == 1:
            self._set_low(self._x._low * self._y._high)
            self._set_lowflex(self._x._lowflex or self._y._highflex)
        elif self._xsign == 0 and self._ysign == -1:
            self._set_low(self._x._high * self._y._low)
            self._set_lowflex(self._x._highflex or self._y._lowflex)
        elif self._xsign == 1 and self._ysign == 0:
            self._set_low(self._x._high * self._y._low)
            self._set_lowflex(self._x._highflex or self._y._lowflex)
        elif self._xsign == -1 and self._ysign == 0:
            self._set_low(self._x._low * self._y._high)
            self._set_lowflex(self._x._lowflex or self._y._highflex)
        else:
            choice1 = self._x._low * self._y._high
            choice2 = self._x._high * self._y._low
            if choice1 < choice2:
                self._set_low(choice1)
                self._set_lowflex(self._x._lowflex or self._y._highflex)
            elif choice2 < choice1:
                self._set_low(choice2)
                self._set_lowflex(self._x._highflex or self._y._lowflex)
            else:
                self._set_low(choice1)
                self._set_lowflex(self._x._lowflex or self._x._highflex or self._y._lowflex or self._y._highflex)

    def _recompute_high(self, chain = True):
        if chain:
            xrecomputed = False
            if 0 in self.x:
                self._x._recompute()
                xrecomputed = True
            self._recompute_xsign()
            yrecomputed = False
            if 0 in self.y:
                self._y._recompute()
                yrecomputed = True
            self._recompute_ysign()
            if self._xsign == 1 and not yrecomputed:
                self._y._recompute_high()
            if self._ysign == 1 and not xrecomputed:
                self._x._recompute_high()
            if self._xsign == -1 and not yrecomputed:
                self._y._recompute_low()
            if self._ysign == -1 and not xrecomputed:
                self._x._recompute_low()
        if self._xsign == 1 and self._ysign == 1:
            self._set_high(self._x._high * self._y._high)
            self._set_highflex(self._x._highflex or self._y._highflex)
        elif self._xsign == 1 and self._ysign == -1:
            self._set_high(self._x._low * self._y._high)
            self._set_highflex(self._x._lowflex or self._y._highflex)
        elif self._xsign == -1 and self._ysign == 1:
            self._set_high(self._x._high * self._y._low)
            self._set_highflex(self._x._highflex or self._y._lowflex)
        elif self._xsign == -1 and self._ysign == -1:
            self._set_high(self._x._low * self._y.low)
            self._set_highflex(self._x._lowflex or self._y._lowflex)
        elif self._xsign == 0 and self._ysign == 1:
            self._set_high(self._x._high * self._y._high)
            self._set_highflex(self._x._highflex or self._y._highflex)
        elif self._xsign == 0 and self._ysign == -1:
            self._set_high(self._x._low * self._y._low)
            self._set_highflex(self._x._lowflex or self._y._lowflex)
        elif self._xsign == 1 and self._ysign == 0:
            self._set_high(self._x._high * self._y._high)
            self._set_highflex(self._x._highflex or self._y._highflex)
        elif self._xsign == -1 and self._ysign == 0:
            self._set_high(self._x._low * self._y._low)
            self._set_highflex(self._x._lowflex or self._y._lowflex)
        else:
            choice1 = self._x._low * self._y._low
            choice2 = self._x._high * self._y._high
            if choice1 > choice2:
                self._set_high(choice1)
                self._set_highflex(self._x._lowflex or self._y._lowflex)
            elif choice2 > choice1:
                self._set_high(choice2)
                self._set_highflex(self._x._highflex or self._y._highflex)
            else:
                self._set_high(choice1)
                self._set_highflex(self._x._lowflex or self._x._highflex or self._y._lowflex or self._y._highflex)

    def _flex(self):
        return self._x._flex() or self._y._flex()

class Valuation_pow(Valuation):
    def __init__(self, parent, x, n):
        CRE.__init__(self, parent)
        self._n = n
        if self._n > 0:
            self._x = x
        elif self._n < 0:
            self._x = x.__invert__()
            self._n = -self._n
        else:
            raise ValueError, "bug in valuation code"
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
        self._recompute_signs()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)

    def _recompute_signs(self):
        if self._x._low > 0:
            self._xsign = 1
        elif self._x._high < 0:
            self._xsign = -1
        else:
            self._xsign = 0

    def _recompute_low(self, chain = True):
        if chain:
            recomputed = False
            if 0 in self._x:
                self._x._recompute()
                recomputed = True
                self._recompute_sign()
            if (self._xsign == 1 or self._n % 2 == 1) and not recomputed:
                self._x._recompute_low()
            elif (self._xsign == -1) and not recomputed:
                self._x._recompute_high()
        if self._xsign == 1:
            self._set_low(self._x._low ** self._n)
            self._set_lowflex(self._x._lowflex)
        elif self._xsign == -1:
            if self._n % 2 == 0:
                self._set_low(self._x._high ** self._n)
                self._set_lowflex(self._x._highflex)
            else:
                self._set_low(self._x._low ** self._n)
                self._set_lowflex(self._x._lowflex)
        else:
            if self._n % 2 == 0:
                self._set_low(Integer(0))
                self._set_lowflex((self._x._low == 0 and self._x._lowflex) or (self._x._high == 0 and self._x._highflex))
            else:
                self._set_low(self._x._low ** self._n)
                self._set_lowflex(self._x._lowflex)

    def _recompute_high(self, chain = True):
        if chain:
            recomputed = False
            if 0 in self._x:
                self._x._recompute()
                recomputed = True
            self._recompute_sign()
            if (self._xsign == 1 or self._n % 2 == 1) and not recomputed:
                self._x._recompute_high()
            elif (self._xsign == -1) and not recomputed:
                self._x._recompute_low()
        if self._xsign == 1:
            self._set_high(self._x._high ** self._n)
            self._set_highflex(self._x._highflex)
        elif self._xsign == -1:
            if self._n % 2 == 0:
                self._set_high(self._x._low ** self._n)
                self._set_highflex(self._x._lowflex)
            else:
                self._set_high(self._x._high ** self._n)
                self._set_highflex(self._x._highflex)
        else:
            if self._n % 2 == 0:
                if self._x._high > -self._x._low:
                    self._set_high(self._x._high ** self._n)
                    self._set_highflex(self._x._highflex)
                elif -self._x._low > self._x._high:
                    self._set_high(self._x._low ** self._n)
                    self._set_highflex(self._x._lowflex)
                else:
                    self._set_high(self._x._high ** self._n)
                    self._set_highflex(self._x._lowflex or self._x._highflex)
            else:
                self._set_high(self._x._high ** self._n)
                self._set_highflex(self._x._highflex)

    def _flex(self):
        return self._x._flex()

class Valuation_neg(Valuation):
    def __init__(self, parent, x):
        CRE.__init__(self, parent)
        self._x = x
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)

    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_high()
        self._set_low(-self._x._high)
        self._set_lowflex(self._x._highflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_low()
        self._set_high(-self._x._low)
        self._set_highflex(self._x._lowflex)

    def _flex(self):
        return self._x._flex()

class Valuation_invert(Valuation):
    def __init__(self, parent, x):
        CRE.__init__(self, parent)
        self._x = x
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)


    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_high()
        if self._x._high == 0:
            self._set_low(minfinity)
        else:
            self._set_low(~self._x._high)
        self._set_lowflex(self._x._highflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_low()
        if self._x._low == 0:
            self._set_high(infinity)
        else:
            self._set_high(~self._x._low)
        self._set_highflex(self._x._lowflex)

    def _flex(self):
        return self._x._flex()

class Valuation_abs(Valuation):
    def __init__(self, parent, x):
        CRE.__init__(self, parent)
        self._x = x
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = false)

    def _recompute(self, chain = True):
        if chain:
            self._x.recompute()
        self._recompute_signs()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)

    def _recompute_signs(self):
        if self._x._low > 0:
            self._xsign = 1
        elif self._x._high < 0:
            self._xsign = -1
        else:
            self._xsign = 0

    def _recompute_low(self, chain = True):
        if chain:
            recomputed = False
            if 0 in self._x:
                self._x._recompute()
                recomputed = True
            self._recompute_signs()
            if self._xsign == 1 and not recomputed:
                self._x._recompute_low()
            elif self._xsign == -1 and not recomputed:
                self._x._recompute_high()
        if self._xsign == 1:
            self._set_low(self._x._low)
            self._set_lowflex(self._x._lowflex)
        elif self._xsign == -1:
            self._set_low(self._x._high)
            self._set_lowflex(self._x._highflex)
        else:
            self._set_low(Integer(0))
            self._set_lowflex((self._x._low == 0 and self._lowflex) or (self._x._high == 0 and self._highflex))

    def _recompute_high(self, chain = True):
        if chain:
            recomputed = False
            if 0 in self._x:
                self._x._recompute()
                recomputed = True
            self._recompute_signs()
            if self._xsign == 1 and not recomputed:
                self._x._recompute_high()
            elif self._xsign == -1 and not recomputed:
                self._x._recompute_low()
        if self._xsign == 1:
            self._set_high(self._x._high)
            self._set_highflex(self._x._highflex)
        elif self._xsign == -1:
            self._set_high(self._x._low)
            self._set_highflex(self._x._highflex)
        else:
            if self._x._high > -self._x._low:
                self._set_high(self._x._high)
                self._set_highflex(self._x._highflex)
            elif -self._x._low > self._x._high:
                self._set_high(-self._x._low)
                self._set_highflex(self._x._lowflex)
            else:
                self._set_high(self._x._high)
                self._set_highflex(self._x._lowflex or self._x._highflex)

    def _flex(self):
        return self._x._flex()

class Valuation_floor(Valuation):
    def __init__(self, parent, x):
        CRE.__init__(self, parent)
        self._x = x
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)

    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_low()
        self._set_low(self._x._low.floor())
        self._set_lowflex(self._x._lowflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_high()
        self._set_high(self._x._high.floor())
        self._set_highflex(self._x._highflex)

    def _flex(self):
        return self._x._flex()

class Valuation_ceil(Valuation):
    def __init__(self, parent, x):
        CRE.__init__(self, parent)
        self._x = x
        self._recompute()

    def _fix(self, halt = None):
        self._x._fix(halt)
        self._recompute(chain = False)

    def _recompute(self, chain = True):
        if chain:
            self._x._recompute()
        self._recompute_low(chain = False)
        self._recompute_high(chain = False)
        if self._x._low == self._x._high:
            self._x = self.parent()(self._x)

    def _recompute_low(self, chain = True):
        if chain:
            self._x._recompute_low()
        self._set_low(self._x._low.ceil())
        self._set_lowflex(self._x._lowflex)

    def _recompute_high(self, chain = True):
        if chain:
            self._x._recompute_high()
        self._set_high(self._x._high.ceil())
        self._set_highflex(self._x._highflex)

    def _flex(self):
        return self._x._flex()


def _fix_one(a, b, halt = None):
    """
    Fixes either a or b
    """
    try:
        a._fix(halt)
    except PrecisionError:
        b._fix(halt)

def _is_unbounded(a):
    return a._low == minfinity or a._high == infinity
