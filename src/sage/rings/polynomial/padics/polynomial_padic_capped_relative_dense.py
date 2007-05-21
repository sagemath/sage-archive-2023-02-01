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
        """
        TESTS:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: R([K(13), K(1)])
        (1 + O(13^7))*t + (13 + O(13^8))
        sage: T.<t> = ZZ[]
        sage: R(t + 2)
        (1 + O(13^7))*t + (2 + O(13^7))
        """
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if construct:
            (self._poly, self._valbase, self._relprecs, self._normalized, self._valaddeds, self._list) = x #the last two of these may be None
            return
	elif is_gen:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            self._poly = PolynomialRing(ZZ, parent.variable_name()).gen()
            self._valbase = 0
            self._valaddeds = [infinity, 0]
            self._relprecs = [infinity, parent.base_ring().precision_cap()]
            self._normalized = True
            self._list = None
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
                p = parent.base_ring().prime()
                self._relprecs = [c.valuation(p) + parent.base_ring().precision_cap() for c in x.list()]
                self._comp_valaddeds()
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

        if len(x) == 1 and not x[0]:
            x = []
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
        if self.degree() == -1 and self._valbase == infinity:
            return []
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
        """
        Returns a list of coefficients of self.

        NOTE:
        The length of the list returned may be greater
        than expected since it includes any leading zeros
        that have finite absolute precision.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = 2*t^3 + 169*t - 1
        sage: a
        (2 + O(13^7))*t^3 + (13^2 + O(13^9))*t + (12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + O(13^7))
        sage: a.list()
        [12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + O(13^7),
         13^2 + O(13^9),
         0,
         2 + O(13^7)]
         """

        if self._list is None:
            self._comp_list()
        return list(self._list)

    def _repr(self, name=None):
        # TODO: what is new here (that doesn't come from parent class)?
        s = " "
        m = self.degree() + 1
        r = reversed(xrange(m))
        if name is None:
            name = self.parent().variable_name()
        coeffs = self.list()
        for n in r:
            x = coeffs[n]
            if x.valuation() < infinity:
                if n != m-1:
                    s += " + "
                x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        if s==" ":
            return "0"
        return s[1:]

    def content(self):
        """
        Returns the content of self.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = 13^7*t^3 + K(169,4)*t - 13^4
        sage: a.content()
        13^2 + O(13^9)
        """
        if self._normalized:
            return self._valbase
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self.base_ring()(self.base_ring().prime_pow(min(self._valaddeds) + self._valbase))

    def lift(self):
        """
        Returns an integer polynomial congruent to this one modulo the precision of each coefficient.

        NOTE: The lift that is returned will not necessarily be the same for polynomials with
              the same coefficients (ie same values and precisions): it will depend on how
              the polynomials are created.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = 13^7*t^3 + K(169,4)*t - 13^4
        sage: a.lift()
        62748517*t^3 + 169*t - 28561
        """
        return self.base_ring().prime_pow(self._valbase) * self._poly

    def __getitem__(self, n):
        """
        Returns the coefficient of x^n

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = 13^7*t^3 + K(169,4)*t - 13^4
        sage: a[1]
        13^2 + O(13^4)
        """
        if n >= len(self._relprecs):
            return self.base_ring()(0)
        if not self._list is None:
            return self._list[n]
        return self.base_ring()(self.base_ring().prime_pow(self._valbase) * self._poly[n], absprec = self._valbase + self._relprecs[n])

    def __getslice__(self, i, j):
        """
        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = 13^7*t^3 + K(169,4)*t - 13^4
        sage: a[1:2]
        (13^2 + O(13^4))*t
        """
        if j > len(self._relprecs):
            j = len(self._relprecs)
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly[i:j], self._valbase, [infinity]*i + self._relprecs[i:j], False, None if self._valaddeds is None else [infinity]*i + self._valaddeds[i:j], None if self._list is None else [self.base_ring()(0)] * i + self._list[i:j]), construct = True)

    def _add_(self, right):
        """
        Returns the sum of self and right.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = t^4 + 17*t^2 + 1
        sage: b = -t^4 + 9*t^2 + 13*t - 1
        sage: c = a + b; c
        (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + (O(13^7))
        sage: c.list()
        [O(13^7), 13 + O(13^8), 2*13 + O(13^7), 0, O(13^7)]
        """
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
                                                      (selfpoly + rightpoly, \
                                                       baseval, \
                                                       [min(a + self._valbase - baseval, b + right._valbase - baseval) for (a, b) in
                                                              zip(_extend_by_infinity(self._relprecs, max(len(self._relprecs), len(right._relprecs))), \
                                                                  _extend_by_infinity(right._relprecs, max(len(self._relprecs), len(right._relprecs))))], \
                                                       False, None, None), construct = True)

    def _sub_(self, right):
        """
        Returns the sum of self and right.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = t^4 + 17*t^2 + 1
        sage: b = t^4 - 9*t^2 - 13*t + 1
        sage: c = a - b; c
        (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + (O(13^7))
        sage: c.list()
        [O(13^7), 13 + O(13^8), 2*13 + O(13^7), 0, O(13^7)]
        """
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
                                                      (selfpoly - rightpoly, \
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

        NOTE: The timing of normalization in arithmetic operations
        may very well change as we do more tests on the relative time
        requirements of these operations.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = t^4 + 17*t^2 + 1
        sage: b = -t^4 + 9*t^2 + 13*t - 1
        sage: c = a + b; c
        (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + (O(13^7))
        sage: d = R([K(1,4), K(2, 6), K(1, 5)]); d
        (1 + O(13^5))*t^2 + (2 + O(13^6))*t + (1 + O(13^4))
        sage: e = c * d; e
        (2*13 + O(13^6))*t^4 + (5*13 + O(13^6))*t^3 + (4*13 + O(13^5))*t^2 + (13 + O(13^5))*t + (O(13^7))
        sage: e.list()
        [O(13^7),
         13 + O(13^5),
         4*13 + O(13^5),
         5*13 + O(13^6),
         2*13 + O(13^6),
         O(13^7),
         O(13^7)]
        """
        self._normalize()
        right._normalize()
        zzpoly = self._poly * right._poly
        if len(self._relprecs) == 0 or len(right._relprecs) == 0:
            return self.parent()(0)
        n = Integer(len(self._relprecs) + len(right._relprecs) - 1).exact_log(2) + 1
        precpoly1 = self._getprecpoly(n) * right._getvalpoly(n)
        precpoly2 = self._getvalpoly(n) * right._getprecpoly(n)
        # These two will be the same length
        tn = Integer(1) << n
        preclist = [min(a.valuation(tn), b.valuation(tn)) for (a, b) in zip(precpoly1.list(), precpoly2.list())]
        answer = Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase + right._valbase, preclist, False, None, None), construct = True)
        answer._reduce_poly()
        return answer

    def _lmul_(self, right):
        return self._rmul_(right)

    def _rmul_(self, left):
        """
        Returns self multiplied by a constant

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: a = t^4 + K(13,5)*t^2 + 13
        sage: K(13,7) * a
        (13 + O(13^7))*t^4 + (13^2 + O(13^6))*t^2 + (13^2 + O(13^8))
        """
        self._comp_valaddeds()
        if left != 0:
            val, unit = left._val_unit()
            left_rprec = left.precision_relative()
            relprecs = [min(left_rprec + self._valaddeds[i], self._relprecs[i]) for i in range(len(self._relprecs))]
        elif left._is_exact_zero():
            return Polynomial_padic_capped_relative_dense(self.parent(), [])
        else:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.parent()(0), self._valbase + left.valuation(), self._valaddeds, False, self._valaddeds, None), construct = True)
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly._rmul_(unit), self._valbase + val, relprecs, False, self._valaddeds, None), construct = True)

    def _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES:
        sage: K = Qp(13,2)
        sage: R.<t> = K[]
        sage: a = t^4 + 13*t^2 + 4
        sage: -a
        (12 + 12*13 + O(13^2))*t^4 + (12*13 + 12*13^2 + O(13^3))*t^2 + (9 + 12*13 + O(13^2))
        """
        return Polynomial_padic_capped_relative_dense(self.parent(), (-self._poly, self._valbase, self._relprecs, False, self._valaddeds, None), construct = True)

    def lshift_coeffs(self, shift, no_list = False):
        """
        Returns a new polynomials whose coefficients are multiplied by p^shift.

        EXAMPLES:
        sage: K = Qp(13, 4)
        sage: R.<t> = K[]
        sage: a = t + 52
        sage: a.lshift_coeffs(3)
        (13^3 + O(13^7))*t + (4*13^4 + O(13^8))
        """
        if shift < 0:
            return self.rshift_coeffs(-shift, no_list)
        if no_list or self._list is None:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, None), construct = True)
        else:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, [c.__lshift__(shift) for c in self._list]), construct = True)

    def rshift_coeffs(self, shift, no_list = False):
        """
        Returns a new polynomial whose coefficients are p-adiclly shifted to the right by shift.

        NOTES: Type Qp(5)(0).__rshift__? for more information.

        EXAMPLES:
        sage: K = Zp(13, 4)
        sage: R.<t> = K[]
        sage: a = t^2 + K(13,3)*t + 169; a
        (1 + O(13^4))*t^2 + (13 + O(13^3))*t + (13^2 + O(13^6))
        sage: b = a.rshift_coeffs(1); b
        (1 + O(13^2))*t + (13 + O(13^5))
        sage: b.list()
        [13 + O(13^5), 1 + O(13^2), O(13^3)]
        sage: b = a.rshift_coeffs(2); b
        (1 + O(13^4))
        sage: b.list()
        [1 + O(13^4), O(13^1), O(13^2)]
        """
        if shift < 0:
            return self.lshift_coeffs(-shift, no_list) # We can't just absorb this into the next if statement because we allow rshift to preserve _normalized
        if self.base_ring().is_field() or shift <= self._valbase:
            if no_list or self._list is None:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, None), construct = True)
            else:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, [c.__rshift__(shift) for c in self._list]), construct = True)
        else:
            shift = shift - self._valbase
            fdiv = self.base_ring().prime_pow(shift)
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly // fdiv, 0, [0 if a <= shift else a - shift for a in self._relprecs], False, None, None), construct = True)

    #def __floordiv__(self, right):
    #    if is_Polynomial(right) and right.is_constant() and right[0] in self.base_ring():
    #        d = self.base_ring()(right[0])
    #    elif (right in self.base_ring()):
    #        d = self.base_ring()(right)
    #    else:
    #        raise NotImplementedError
    #    return self._rmul_(self.base_ring()(~d.unit_part())).rshift_coeffs(d.valuation())

    def _unsafe_mutate(self, n, value):
        """
        It's a really bad idea to use this function for p-adic polynomials.  There are speed issues, and it may not be bug-free currently.
        """
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
        """
        Returns the degree of self, ie the largest n so that the coefficient of x^n is distinguishable from 0.
        """
	self._normalize()
        return self._poly.degree()

    def prec_degree(self):
        """
        Returns the largest n so that precision information is stored about the coefficient of x^n.

        Always greater than or equal to degree.
        """
        return len(self._relprecs) - 1

    def precision_absolute(self, n = None):
        """
        Returns absolute precision information about self.

        INPUT:
        self -- a p-adic polynomial
        n -- None or an integer (default None).

        OUTPUT:
        If n == None, returns a list of absolute precisions of coefficients.  Otherwise,
        returns the absolute precision of the coefficient of x^n.
        """
        if n is None:
            return [c + self._valbase for c in self._relprecs]
        return self._relprecs[n] + self._valbase

    def precision_relative(self, n = None):
        """
        Returns relative precision information about self.

        INPUT:
        self -- a p-adic polynomial
        n -- None or an integer (default None).

        OUTPUT:
        If n == None, returns a list of relative precisions of coefficients.  Otherwise,
        returns the relative precision of the coefficient of x^n.
        """
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
        """
        Returns valuation information about self's coefficients.

        INPUT:
        self -- a p-adic polynomial
        n -- None or an integer (default None).

        OUTPUT:
        If n == None, returns a list of valuations of coefficients.  Otherwise,
        returns the valuation of the coefficient of x^n.
        """
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
        """
        Returns the valuation of self

        INPUT:
        self -- a p-adic polynomial
        val_of_var -- None or a rational (default None).

        OUTPUT:
        If val_of_var == None, returns the largest power of the variable dividing self.  Otherwise,
        returns the valuation of self where the variable is assigned valuation val_of_var
        """
        if val_of_var is None:
            return self._poly.valuation()
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self._valbase + min([self._valaddeds[i] + val_of_var * i for i in range(len(self._valaddeds))])

    def reverse(self, n = None):
        """
        Returns a new polynomial whose coefficients are the reversed coefficients of self, where self is considered as a polynomial of degree n.

        If n is None, defaults to the degree of self.
        If n is smaller than the degree of self, some coefficients will be discarded.

        EXAMPLES:
        sage: K = Qp(13,7)
        sage: R.<t> = K[]
        sage: f = t^3 + 4*t; f
        (1 + O(13^7))*t^3 + (4 + O(13^7))*t
        sage: f.reverse()
        (4 + O(13^7))*t^2 + (1 + O(13^7))
        sage: f.reverse(3)
        (4 + O(13^7))*t^2 + (1 + O(13^7))
        sage: f.reverse(2)
        (4 + O(13^7))*t
        sage: f.reverse(4)
        (4 + O(13^7))*t^3 + (1 + O(13^7))*t
        sage: f.reverse(6)
        (4 + O(13^7))*t^5 + (1 + O(13^7))*t^3
        """
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
            L = self._list[:(n+1)] + [self.base_ring()(0)] * (n - self.prec_degree())
            L.reverse()
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.parent()(zzlist), self._valbase, relprec, self._normalized, valadded, L), construct = True)

    def rescale(self, a):
        r"""
        Return f(a*X)

        NOTE:  Need to write this function for integer polynomials before this works.

        EXAMPLES:
        sage.: K = Zp(13, 5)
        sage.: R.<t> = K[]
        sage.: f = t^3 + K(13, 3) * t
        sage.: f.rescale(2)
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
            return self.parent().base_extend(self.base_ring().fraction_field())(self).rescale(a)
        if self.base_ring().is_field() and a.valuation() < 0:
            D = self.prec_degree()
            return a**D * self.reverse(D).rescale(~a).reverse(D)
        aval = a.valuation()
        arprec = a.precision_relative()
        if self._valaddeds is None:
            self._comp_valaddeds()
        valadded = [self._valaddeds[i] + aval * i for i in range(len(self._valaddeds))]
        relprec = [infinity if (self._relprecs[i] is infinity) else (min(self._relprecs[i] - self._valaddeds[i], arprec) + aval * i + self._valaddeds[i]) for i in range(len(self._relprecs))]
        relprec[0] = self._relprecs[0]
        if a == 0:
            zzpoly = self._poly.parent()(0)
        else:
            zzpoly = self._poly.rescale(Integer(a))
        return Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase, relprec, False, valadded, None), construct = True)

    def quo_rem(self, right):
        return self._quo_rem_naive(right)

    def _quo_rem_naive(self, right):
        """
        An implementation of quo_rem that doesn't have good run-time or precision characteristics.
        """
        K = self.base_ring().fraction_field()
        f = self.base_extend(K)
        g = right.base_extend(K)
        if g == 0:
            raise ZeroDivisionError, "cannot divide by a polynomial indistinguishable from 0"
        x = f.parent().gen()
        quo = f.parent()(0)
        while f.degree() >= g.degree():
            a = f.leading_coefficient() / g.leading_coefficient()
            quo = quo + a * (x ** (f.degree() - g.degree()))
            f = f - a * (x ** (f.degree() - g.degree())) * g
        return (quo, f)

    #def gcd(self, right):
    #    raise NotImplementedError

    #def lcm(self, right):
    #    raise NotImplementedError

    def xgcd(self, right):
        return self._xgcd(right)

    #def discriminant(self):
    #    raise NotImplementedError

    def disc(self):
        return self.discriminant()

    #def resultant(self):
    #    raise NotImplementedError

    def newton_slopes(self):
        """
        Returns a list of the Newton slopes of this polynomial.  These are the valuations of the roots of this polynomial.

        EXAMPLES:
        sage: K = Qp(13)
        sage: R.<t> = K[]
        sage: f = t^4 + 13^5*t^2 + 4*13^2*t - 13^7
        sage: f.newton_polygon()
        [(0, 7), (1, 2), (4, 0)]
        sage: f.newton_slopes()
        [5, 2/3, 2/3, 2/3]
        """
        polygon = self.newton_polygon()
        if polygon == []:
            return []
        answer = [infinity] * polygon[0][0]
        for m in range(1, len(polygon)):
            dx = polygon[m][0] - polygon[m - 1][0]
            dy = polygon[m][1] - polygon[m - 1][1]
            answer.extend([-dy / dx] * dx)
        return answer

    def newton_polygon(self):
        r"""
        Returns a list of vertices of the Newton polygon of this polynomial.

        NOTES:
        The vertices are listed so that the first coordinates are strictly increasing, up to the polynomial's degree (not the limit of available precision information).  Also note that if some coefficients have very low precision an error is raised.

        EXAMPLES:
        sage: K = Qp(13)
        sage: R.<t> = K[]
        sage: f = t^4 + 13^5*t^2 + 4*13^2*t - 13^7
        sage: f.newton_polygon()
        [(0, 7), (1, 2), (4, 0)]
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
        answer.append((xf, self._valbase + self._valaddeds[xf]))
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
        elif min(self._relprecs) <= 0:
            raise PrecisionError, "Polynomial is not known to high enough precision"
        return self._poly.factor_mod(self.base_ring().prime())

    def factor(self):
        # This will eventually be improved.
        if self == 0:
            raise ValueError, "Factorization of the zero polynomial not defined"
	from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.padics.factory import ZpCA
        base = self.base_ring()
        #print self.list()
        m = min([x.precision_absolute() for x in self.list()])
        #print m
        R = ZpCA(base.prime(), prec = m)
        S = PolynomialRing(R, self.parent().variable_name())
        F = S(self).factor()
        return Factorization([(self.parent()(a), b) for (a, b) in F], base(F.unit()))

def _extend_by_infinity(L, n):
    return L + [infinity] * (n - len(L))
