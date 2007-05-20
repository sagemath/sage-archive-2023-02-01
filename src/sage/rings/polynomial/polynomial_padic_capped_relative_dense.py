import sage.rings.polynomial.polynomial_element_generic
import sage.rings.polynomial.polynomial_element
import sage.rings.integer
import sage.rings.integer_ring
import sage.rings.padics.misc as misc
import sage.rings.padics.precision_error as precision_error
import sage.rings.fraction_field_element as fraction_field_element
import copy

from sage.libs.all import pari, pari_gen
from sage.libs.ntl.all import ZZX_class, ZZ_pX_class
from sage.structure.factorization import Factorization
from sage.rings.infinity import infinity

min = misc.min
ZZ = sage.rings.integer_ring.ZZ
PrecisionError = precision_error.PrecisionError
Integer = sage.rings.integer.Integer
Polynomial = sage.rings.polynomial.polynomial_element.Polynomial
is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial
Polynomial_generic_domain = sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_domain
Polynomial_integer_dense = sage.rings.polynomial.polynomial_element_generic.Polynomial_integer_dense

class Polynomial_padic_capped_relative_dense(Polynomial_generic_domain):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct = False, absprec = infinity, relprec = infinity):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if construct:
            (self._poly, self._valbase, self._relprecs, self._normalized, self._valaddeds, self._list) = x #the last two of these may be None
            return

        #First we list the types that are turned into Polynomials
        if isinstance(x, ZZX_class):
            from polynomial_ring_constructor import PolynomialRing
            x = Polynomial_integer_dense(PolynomialRing(ZZ, parent.variable_name()), x, construct = True)
        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
               x.denominator() == 1:
            #Currently we ignore precision information in the denominator.  This should be changed eventually
            x = x.numerator()

        #We now coerce various types into lists of coefficients.  There are fast pathways for some types of polynomials
        if isinstance(x, Polynomial):
            if x.parent() is self.parent():
                if not absprec is infinity or not relprec is infinity:
                    x._normalize()
                self._poly = x._poly
                self._valbase = x._valbase
                self._valaddeds = x._valaddeds
                self._relprecs = x._relprecs
                self._normalized = x._normalized
                self._list = x._list
                if not absprec is infinity or not relprec is infinity:
                    self._adjust_prec_info(absprec, relprec)
                return
            elif x.base_ring() is ZZ:
                self._poly = x
                self._valbase = Integer(0)
                self._comp_valaddeds()
                self._relprecs = [c + parent.base_ring().precision_cap() for c in self._valaddeds]
                self._normalized = len(self._valaddeds) == 0 or (min(self._valaddeds) == 0)
                self._list = None
                if not absprec is infinity or not relprec is infinity:
                    self._adjust_prec_info(absprec, relprec)
                return
            else:
                x = [parent.base_ring()(a) for a in x.list()]
                check = False
        elif isinstance(x, dict):
            zero = parent.base_ring()(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n + 1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v
        elif isinstance(x, pari_gen):
            x = [parent.base_ring()(w) for w in x.Vecrev()]
            check = False
        #The default behavior if we haven't already figured out what the type is is to assume it coerces into the base_ring as a constant polynomial
        elif not isinstance(x, list):
            x = [x] # constant polynomial

        if check:
            x = [parent.base_ring()(z) for z in x]

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        self._list = x
        self._valaddeds = [a.valuation() for a in x]
        self._valbase = sage.rings.padics.misc.min(self._valaddeds)
        if self._valbase is infinity:
            self._valaddeds = []
            self._relprecs = []
            self._poly = PolynomialRing(ZZ, parent.variable_name())()
            self._normalized = True
            if not absprec is infinity or not relprec is infinity:
                self._adjust_prec_info(absprec, relprec)
        else:
            self._valaddeds = [c - self._valbase for c in self._valaddeds]
            self._relprecs = [a.precision_absolute() - self._valbase for a in x]
            self._poly = PolynomialRing(ZZ, parent.variable_name())([a >> self._valbase for a in x])
            self._normalized = True
            if not absprec is infinity or not relprec is infinity:
                self._adjust_prec_info(absprec, relprec)

    def _normalize(self):
        # Currently slow: need to optimize
        if not self._normalized:
            if self._valaddeds is None:
                self._comp_valaddeds()
            val = sage.rings.padics.misc.min(self._valaddeds)
            prime_pow = self.base_ring().prime_pow
            selflist = self._poly.list()
            if val is infinity:
                pass
            elif val != 0:
                self._relprecs = [prec - val for prec in self._relprecs]
                self._poly.ntl_set_directly([Integer(0) if (e is infinity) else ((c // prime_pow(val)) % prime_pow(e)) for (c,e) in zip(selflist, self._relprecs)])
                self._valbase += val
                self._valaddeds = [c - val for c in self._valaddeds]
            else:
                self._poly.ntl_set_directly([Integer(0) if (e is infinity) else (c % prime_pow(e)) for (c,e) in zip(selflist, self._relprecs)])
            self._normalized = True

    def _reduce_poly(self):
        selflist = self._poly.list()
        prime_pow = self.base_ring().prime_pow
        self._poly.ntl_set_directly([Integer(0) if (e is infinity) else (c % prime_pow(e)) for (c, e) in zip(selflist, self._relprecs)])

    def _comp_list(self):
        polylist = self._poly.list()
        polylen = len(polylist)
        self._list = [self.base_ring()(polylist[i], absprec = self._relprecs[i]) << self._valbase for i in range(polylen)] \
                     + [self.base_ring()(0, absprec = self._relprecs[i] + self._valbase) for i in range(polylen, len(self._relprecs))]

    def _comp_valaddeds(self):
        self._valaddeds = []
        for i in range(self._poly.degree() + 1):
            tmp = self._poly.list()[i].valuation(self.parent().base_ring().prime())
            if tmp is infinity:
                self._valaddeds.append(self._relprecs[i])
            else:
                self._valaddeds.append(tmp)
        for i in range(self._poly.degree() + 1, len(self._relprecs)):
            self._valaddeds.append(self._relprecs[i])

    def _adjust_prec_info(self, absprec=infinity, relprec=infinity):
        r"""
        Assumes that self._poly, self._val and self._relprec are set initially and adjusts self._val and self._relprec to the termwise minimum of absprec and relprec.
        """
        return

#         min = sage.rings.padics.misc.min
#         slen = len(self._relprec)
#         if isinstance(absprec, list):
#             alen = len(absprec)
#         elif absprec is infinity:
#             alen = 0
#             absprec = []
#         else:
#             alen = 1
#         if isinstance(relprec, list):
#             rlen = len(relprec)
#         elif relprec is infinity:
#             rlen = 0
#             relprec = []
#         else:
#             rlen = 1
#         preclen = max(slen, rlen, alen)
#         if not isinstance(absprec, list):
#             absprec = [absprec] * preclen
#         if not isinstance(relprec, list):
#             relprec = [relprec] * preclen
#         vallist = [c.valuation(self.base_ring().prime()) + self._val for c in self._poly.list()] #######
#         vmin = min(vallist)
#         amin = min(absprec)
#         if amin < vmin:
#             vmin = amin
#         if vmin < self._val:
#             vadjust =

#         if not isinstance(absprec, list):
#             self._val = min(vallist + [absprec])
#             absprec = [absprec] * preclen
#         else:
#             self._val = padics.misc.min(vallist + absprec)
#             absprec = absprec + [infinity] * (preclen - len(absprec))
#         if self._val is infinity:
#             self._relprec = []
#             return
#         if not isinstance(relprec, list):
#             relprec = [relprec] * preclen
#         else:
#             relprec = relprec + [parent.base_ring().precision_cap()] * (preclen - len(relprec))
#         self._relprec = [min(a, v + r) - self._val for (a, r, v) in zip(absprec, relprec, vallist)]
#Remember to normalize at the end if self._normalized is true because you need to reduce mod p^n

    def _getprecpoly(self, n):
        one = Integer(1)
        return self._poly.parent()([(0 if (c is infinity) else (one << (n * c))) for c in self._relprecs])

    def _getvalpoly(self, n):
        one = Integer(1)
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self._poly.parent()([(0 if (c is infinity) else (one << (n * c))) for c in self._valaddeds] + \
                                   [(0 if (c is infinity) else (one << (n * c))) for c in self._relprecs[len(self._valaddeds):]])

    def list(self):
        if self._list is None:
            self._comp_list()
        return self._list

    def content(self):
        if self._normalized:
            return self._valbase
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self.base_ring()(self.base_ring().prime_pow(min(self._valaddeds) + self._valbase))

    def lift(self):
        return self._poly

    def __getitem__(self, n):
        if n >= len(self._relprecs):
            return self.base_ring()(0)
        if not self._list is None:
            return self._list[n]
        return self.base_ring()(self._poly[n], absprec = self._valbase + self._relprecs[n])

    def __getslice__(self, i, j):
        if j > len(self._relprecs):
            j = len(self._relprecs)
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly[i:j], self._valbase, [infinity]*i + self._relprecs[i:j], False, None if self._valaddeds is None else [infinity]*i + self._valaddeds[i:j], None if self._list is None else [self.base_ring()(0)] * i + self._list[i:j]), construct = True)

    def _add_(self, right):
        selfpoly = self._poly
        rightpoly = right._poly
        if self._valbase > right._valbase:
            selfpoly = selfpoly * self.base_ring().prime_pow(self._valbase - right._valbase)
            baseval = right._valbase
        elif self._valbase < right._valbase:
            rightpoly = rightpoly * self.base_ring().prime_pow(right._valbase - self._valbase)
            baseval = self._valbase
        else:
            baseval = self._valbase
        # Currently we don't reduce the coefficients of the answer modulo the appropriate power of p or normalize
        return Polynomial_padic_capped_relative_dense(self.parent(), \
                                                      (self._poly + right._poly, \
                                                       baseval, \
                                                       [min(a + self._valbase - baseval, b + right._valbase - baseval) for (a, b) in
                                                              zip(_extend_by_infinity(self._relprecs, max(len(self._relprecs), len(right._relprecs))), \
                                                                  _extend_by_infinity(right._relprecs, max(len(self._relprecs), len(right._relprecs))))], \
                                                       False, None, None), construct = True)

    def _sub_(self, right):
        selfpoly = self._poly
        rightpoly = right._poly
        if self._valbase > right._valbase:
            selfpoly = selfpoly * self.base_ring().prime_pow(self._valbase - right._valbase)
            baseval = right._valbase
        elif self._valbase < right._valbase:
            rightpoly = rightpoly * self.base_ring().prime_pow(right._valbase - self._valbase)
            baseval = self._valbase
        else:
            baseval = self._valbase
        # Currently we don't reduce the coefficients of the answer modulo the appropriate power of p or normalize
        return Polynomial_padic_capped_relative_dense(self.parent(), \
                                                      (self._poly - right._poly, \
                                                       baseval, \
                                                       [min(a + self._valbase - baseval, b + right._valbase - baseval) for (a, b) in
                                                              zip(_extend_by_infinity(self._relprecs, max(len(self._relprecs), len(right._relprecs))), \
                                                                  _extend_by_infinity(right._relprecs, max(len(self._relprecs), len(right._relprecs))))], \
                                                       False, None, None), construct = True)

    def _mul_(self, right):
        r"""
        Multiplies self and right.

        ALGORITHM: We use an algorithm thought up by Joe Wetherell to
        find the precisions of the product.  It works as follows:
        Suppose $f = \sum_i a_i x^i$ and $g = \sum_j b_j x^j$. Let $N
        = \max(\deg f, \deg g) + 1$ (in the actual implementation we
        use $N = 2^{\lfloor \log_2\max(\deg f, \deg g)\rfloor + 1}$).
        The valuations and absolute precisions of each coefficient
        contribute to the absolute precision of the kth coefficient of
        the product in the following way: for each $i + j = k$, you
        take the valuation of $a_i$ plus the absolute precision of
        $b_j$, and then take the valuation of $b_j$ plus the absolute
        precision of $a_i$, take the minimum of those two, and then
        take the minimum over all $i$, $j$ summing to $k$.

        You can simulate this as follows. Construct new polynomials of
        degree $N$:

        \begin{align*}
        A &= \sum_i N^{\mbox{valuation of $a_i$}} x^i \\
        B &= \sum_j N^{\mbox{absolute precision of $b_j$}} x^j \\
        C &= \sum_i N^{\mbox{absolute precision of $a_i$}} x^i \\
        D &= \sum_j N^{\mbox{valuation of $b_j$}} x^j \\
        \end{align*}

        Now you compute AB and CD. Because you're representing things
        'N-adically', you don't get any 'overflow', and you can just
        read off what the precisions of the product are. In fact it
        tells you more, it tells you exactly how many terms of each
        combination of valuation modulus contribute to each term of
        the product (though this feature is not currently exposed in
        our implementation.

        Since we're working 'N-adically' we can just consider
        $N^{\infty} = 0$.

        NOTE: When polynomials are normalized in arithmetic operations
        may very well change as we do more tests on the relative time
        requirements of these operations.
        """
        self._normalize()
        right._normalize()
        zzpoly = self._poly * right._poly
        n = Integer(len(self._relprecs) + len(right._relprecs) - 1).exact_log(2) + 1
        precpoly1 = self._getprecpoly(n) * right._getvalpoly(n)
        precpoly2 = self._getvalpoly(n) * right._getprecpoly(n)
        # These two will be the same length
        tn = Integer(1) << n
        preclist = [min(a, b).valuation(tn) for (a, b) in zip(precpoly1.list(), precpoly2.list())]
        answer = Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase + right._valbase, preclist, False, None, None), construct = True)
        answer._reduce_poly()
        return answer

    def _lmul_(self, right):
        print "lmul"
        val, unit = right._val_unit()
        self._comp_valaddeds()
        right_rprec = right.precision_relative()
        relprecs = [min(right_rprec + self._valaddeds[i], self._relprecs[i]) for i in range(len(self._relprecs))]
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly._lmul_(unit.lift()), self._valbase + val, relprecs, False, self._valaddeds, None), construct = True)

    def _rmul_(self, left):
        print "rmul"
        val, unit = left._val_unit()
        self._comp_valaddeds()
        left_rprec = left.precision_relative()
        relprecs = [min(left_rprec + self._valaddeds[i], self._relprecs[i]) for i in range(len(self._relprecs))]
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly._rmul_(unit.lift()), self._valbase + val, relprecs, False, self._valaddeds, None), construct = True)

    def _neg_(self):
        return Polynomial_padic_capped_relative_dense(self.parent(), (-self._poly, self._valbase, self._relprecs, False, self._valaddeds, None), construct = True)

    def lshift_coeffs(self, shift, no_list = False):
        if shift < 0:
            return self.rshift_coeffs(-shift, no_list)
        if no_list or self._list is None:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, None), construct = True)
        else:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, [c.__lshift__(shift) for c in self._list]), construct = True)

    def rshift_coeffs(self, shift, no_list = False):
        if shift < 0:
            return self.lshift_coeffs(-shift, no_list) # We can't just absorb this into the next if statement because we allow rshift to preserve _normalized
        if self.base_ring().is_field() or shift <= self._valbase:
            if no_list or self._list is None:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, None), construct = True)
            else:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, [c.__rshift__(shift) for c in self._list]), construct = True)
        else:
            fdiv = self.base_ring().prime_pow(shift - self._valbase)
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly // fdiv, 0, self._relprecs, False, None, None), construct = True)

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and right[0] in self.base_ring():
            d = self.base_ring()(right[0])
        elif (right in self.base_ring()):
            d = self.base_ring()(right)
        else:
            raise NotImplementedError
        return self._rmul_(self.base_ring()(~d.unit_part())).rshift_coeffs(d.valuation())

    def _unsafe_mutate(self, n, value):
        n = int(n)
        value = self.base_ring()(value)
        if self.is_gen():
            raise ValueError, "cannot modify generator"
        if n < 0:
            raise IndexError, "n must be >= 0"
        if self._valbase is infinity:
            if value._is_exact_zero():
                return
            self._valbase = value.valuation()
            if value != 0:
                self._poly._unsafe_mutate(self, n, value.unit_part().lift())
                self._relprecs = [infinity] * n + [value.precision_relative()]
            else:
                self._relprecs = [infinity] * n + [0]
            self._valaddeds = [infinity] * n + [0]
            zero = self.base_ring()(0)
            self._list = [zero] * n + [value]
            self._normalized = True
        elif value.valuation() >= self._valbase:
            # _valbase and _normalized stay the same
            if value != 0:
                self._poly._unsafe_mutate(self, n, (value.__rshift__(self._valbase)).lift())
            else:
                self._poly._unsafe_mutate(self, n, 0)
            if n < len(self._relprecs):
                self._relprecs[n] = value.precision_absolute() - self._valbase
                if not self._valaddeds is None:
                    self._valaddeds[n] = value.valuation() - self._valbase
                if not self._list is None:
                    self._list[n] = value
            else:
                self._relprecs.extend([infinity] * (n - len(self._relprecs)) + [value.precision_absolute() - self._valbase])
                if not self._valaddeds is None:
                    self._valaddeds.extend([infinity] * (n - len(self._relprecs)) + [value.valuation() - self._valbase])
                if not self._list is None:
                    zero = self.base_ring()(0)
                    self._list.extend([zero] * (n - len(self._relprecs)) + [value])
        else:
            basediff = self._valbase - value.valuation()
            self._valbase = value.valuation()
            if not self._valaddeds is None:
                self._valaddeds = [c + basediff for c in self._valaddeds]
            self._poly = self._poly * self.base_ring().prime_pow(basediff)
            if value != 0:
                self._poly._unsafe_mutate(self, n, value.unit_part().lift())
            else:
                self._poly._unsafe_mutate(self, n, 0)
            if n < len(self._relprecs):
                self._relprecs[n] = value.precision_relative()
            else:
                self._relprecs.extend([infinity] * (n - len(self._relprecs)) + [value.precision_relative()])
            self._normalized = False
            if not self._list is None:
                if n < len(self._list):
                    self._list[n] = value
                else:
                    zero = self._base_ring()(0)
                    self._list.extend([zero] * (n - len(self._list)) + [value])

    def _pari_(self, variable = None):
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.list()).Polrev(variable)

    def copy(self):
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.copy(), self._valbase, copy.copy(self._relprecs), self._normalized, copy.copy(self._valaddeds), copy.copy(self._list)), construct = True)

    def degree(self):
        return self._poly.degree()

    def prec_degree(self):
        return len(self._relprecs) - 1

    def precision_absolute(self, n = None):
        if n is None:
            return [c + self._valbase for c in self._relprecs]
        return self._relprecs[n] + self._valbase

    def precision_relative(self, n = None):
        if n is None:
            self._normalize()
            return copy.copy(self._relprecs)
        n = int(n)
        if n < 0 or n >= len(self._relprecs) or self._relprecs[n] is infinity:
            return Integer(0)
        if self._valaddeds is None:
            return self._relprecs[n] - self._poly[n].valuation(self.base_ring().prime())
        else:
            return self._relprecs[n] - self._valaddeds[n]

    def valuation_of_coefficient(self, n = None):
        if n is None:
            self._normalize()
            return [c + self._valbase for c in self._valadded]
        n = int(n)
        if n < 0 or n >= len(self._relprecs):
            return infinity
        if self._valaddeds is None:
            return self._valbase + self._poly[n].valuation(self.base_ring().prime())
        else:
            return self._valbase + self._valaddeds[n]

    def valuation(self, val_of_var = None):
        if val_of_var is None:
            return self._poly.valuation()
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self._valbase + min([self._valaddeds[i] + val_of_var * i for i in range(len(self._valaddeds))])

    def reverse(self, n = None):
        if n is None:
            n = self._poly.degree()
        zzlist = self._poly.list()[:(n+1)] + [0] * (n - self._poly.degree())
        zzlist.reverse()
        relprec = self._relprecs[:(n+1)] + [infinity] * (n - self.prec_degree())
        relprec.reverse()
        if self._valaddeds is None:
            valadded = None
        else:
            valadded = self._valaddeds[:(n+1)] + [infinity] * (n - self.prec_degree())
            valadded.reverse()
        if self._list is None:
            L = None
        else:
            L = self._list[:(n+1)] + [infinity] * (n - self.prec_degree())
            L.reverse()
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.parent()(zzlist), self._valbase, relprec, self._normalized, valadded, L), construct = True)

    def rescale(self, a):
        r"""
        Return f(a*X)
        """
        negval = False
        try:
            a = self.base_ring()(a)
        except ValueError, msg:
            if msg == "element has negative valuation.":
                negval = True
            else:
                raise ValueError, msg
        if negval:
            return parent.base_extend(self.base_ring().fraction_field())(self).rescale(a)
        if self.base_ring().is_field() and a.valuation() < 0:
            D = self.prec_degree()
            return a**D * self.reverse(D).rescale(~a).reverse(D)
        aval = a.valuation()
        arprec = a.precision_relative()
        if self._valaddeds is None:
            self._comp_valaddeds()
        valadded = [self._valaddeds[i] + aval * i for i in range(len(self._valaddeds))]
        relprec = [self._relprecs[0] if (i == 0) else (min(self._relprecs[i] - self._valaddeds[i], arprec) + aval * i) for i in range(len(self._relprecs))]
        if a == 0:
            zzpoly = self._poly.parent()(0)
        else:
            zzpoly = self._poly.rescale(Integer(a))
        return Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase, relprec, False, valadded, None), construct = True)

    def quo_rem(self, right):
        raise NotImplementedError

    def gcd(self, right):
        raise NotImplementedError

    def lcm(self, right):
        raise NotImplementedError

    def xgcd(self, right):
        raise NotImplementedError

    def discriminant(self):
        raise NotImplementedError

    def disc(self):
        return self.discriminant()

    def resultant(self):
        raise NotImplementedError

    def newton_slopes(self):
        polygon = self.newton_polygon()
        if polygon == []:
            return []
        answer = [infinity] * polygon[0][0]
        for m in range(1, len(polygon)):
            dx = polygon[m][0] - polygon[m - 1][0]
            dy = polygon[m][1] - polygon[m - 1][1]
            answer.extend([dy / dx] * dx)
        return answer

    def newton_polygon(self):
        r"""
        Returns a list of vertices of the Newton polygon of this polynomial.

        NOTES:
        The vertices are listed so that the first coordinates are strictly increasing, up to the polynomial's degree (not the limit of available precision information).  Also note that if some coefficients have very low precision an error is raised.
        """
        if self._poly == 0:
            return []
        for x in range(len(self._relprecs)):
            if not self._relprecs[x] is infinity:
                break
        self._comp_valaddeds()
        if self._poly[x] == 0:
            raise PrecisionError, "first term with non-infinite valuation must have determined valuation"
        yrel = self._valaddeds[x]
        answer = [(x, self._valbase + yrel)]
        xf = self._poly.degree()
        if xf == x:
            return answer
        yfrel = self._valaddeds[xf]
        curslope = (yfrel - yrel) / (xf - x)
        for i in range(x + 1, xf):
            yrel += curslope
            if self._valaddeds[i] < yrel:
                if self._relprecs[i] == self._valaddeds[i]:
                    raise PrecisionError, "not enough precision known in coefficient %s to compute newton polygon"%i
                yrel = self._valaddeds[i]
                answer.append((i, self._valbase + yrel))
                curslope = (yfrel - yrel) / (xf - i)
        answer.append((xf, self._valbase + self._valadded[xf]))
        return answer

    def hensel_lift(self, a):
        raise NotImplementedError

    def factor_mod(self):
        r"""
        Returns the factorization of self modulo p.
        """
        self._normalize()
        if self._valbase < 0:
            raise ValueError, "Polynomial does not have integral coefficients"
        elif self._valbase > 0:
            raise ValueError, "Factorization of the zero polynomial not defined"
        return self._poly.factor_mod(self.base_ring().prime())

    def factor(self):
        raise NotImplementedError


def _extend_by_infinity(L, n):
    return L + [infinity] * (n - len(L))
