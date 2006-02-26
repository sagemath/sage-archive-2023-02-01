"""
Univariate Polynomials
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

from sage.structure.element import Element, Element_cmp_
import sage.rings.rational_field
import sage.rings.integer_ring
import sage.rings.rational
import integer
import finite_field
import padic_field
import sage.rings.polynomial_ring
from sage.rings.coerce import bin_op, cmp as coerce_cmp
import arith
import sage.rings.ring_element as ring_element
import sage.rings.euclidean_domain_element as euclidean_domain_element
import sage.rings.integral_domain_element as integral_domain_element
import sage.rings.principal_ideal_domain_element as principal_ideal_domain_element
import integer_ring
import integer_mod_ring
import polynomial_pyx
import rational_field
import complex_field
import fraction_field_element
import fraction_field
from infinity import infinity
import sage.misc.misc as misc
from sage.libs.all import pari, gen
from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZX, ZZX_class, ZZ_p, ZZ_pX, ZZ_pX_class, set_modulus
import sage.misc.latex as latex
import sage.structure.factorization as factorization

from coerce import bin_op

QQ = rational_field.RationalField()

ZZ = integer_ring.IntegerRing()

def is_Polynomial(f):
    return isinstance(f, Polynomial)


class Polynomial(Element_cmp_, ring_element.RingElement):
    """
    Polynomial base class.
    """
    def __init__(self, parent, is_gen = False, construct=False):
        """
        The following examples illustrate creation of elements of
        polynomial rings, and some basic arithmetic.

        First we make a polynomial over the integers and do some arithmetic:
            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: f = x^5 + 2*x^2 + (-1); f
            x^5 + 2*x^2 - 1
            sage: f^2
            x^10 + 4*x^7 - 2*x^5 + 4*x^4 - 4*x^2 + 1

        Next we do arithmetic in a sparse polynomial ring over the integers:
            sage: R = PolynomialRing(IntegerRing(), "x"); x = R.gen(); R
            Univariate Polynomial Ring in x over Integer Ring
            sage: S = PolynomialRing(R, "Z"); Z = S.gen(); S
            Univariate Polynomial Ring in Z over Univariate Polynomial Ring in x over Integer Ring
            sage: f = Z^3 + (x^2-2*x+1)*Z - 3; f
            Z^3 + (x^2 - 2*x + 1)*Z + -3
            sage: f*f
            Z^6 + (2*x^2 - 4*x + 2)*Z^4 + (-6)*Z^3 + (x^4 - 4*x^3 + 6*x^2 - 4*x + 1)*Z^2 + (-6*x^2 + 12*x - 6)*Z + 9
            sage: f^3 == f*f*f
            True

        To have the element print as 'y', give 'y' as the
        second argument to the PolynomialRing constructor.
            sage: y = PolynomialRing(IntegerRing(), 'y').gen()
            sage: y^3 - 2*y
            y^3 - 2*y
        """
        ring_element.RingElement.__init__(self, parent)
        self._is_gen = is_gen

    def _add_(self, right):
        if self.degree() >= right.degree():
            x = list(self.list())
            y = right.list()
        else:
            x = list(right.list())
            y = self.list()

        for i in xrange(len(y)):
            x[i] += y[i]

        return self.polynomial(x)

    def __call__(self, *a):
        """
        Evaluate polynomial at x=a using Horner's rule

        INPUT:
            a -- ring element a; need not be in the coefficient
                 ring of the polynomial.

        OUTPUT:
            the value of f at a.

        EXAMPLES:
            sage: x = Q['x'].gen()
            sage: f = x/2 - 5
            sage: f(3)
            -7/2
            sage: x = Z['x'].gen()
            sage: f = (x-1)^5
            sage: f(2/3)
            -1/243

        AUTHORS:
            -- David Joyner, 2005-04-10
            -- William Stein, 2006-01-22; change so parent
               is determined by the arithmetic
        """
        a = a[0]
        if isinstance(a, tuple):
            a = a[0]
        d = self.degree()
        result = self[d]
        i = d - 1
        while i >= 0:
            result = result * a + self[i]
            i -= 1
        return result

    def _cmp_(self, other):
        """
        EXAMPLES:
            sage: x = QQ['x'].0
            sage: 3*x^3  + 5 > 10*x^2 + 19
            True
            sage: f = x^2 - 2*x + 1; g= x^2 - 1
            sage: f < g
            True
            sage: f > g
            False
            sage: g < f
            False
            sage: g > f
            True
        """
        #if not isinstance(other, Polynomial) or other.parent() != self.parent():
        #    return coerce_cmp(self, other)
        c = cmp(self.degree(), other.degree())
        if c: return c
        return cmp(list(reversed(self.list())), list(reversed(other.list())))

    def __getitem__(self, n):
        raise NotImplementedError

    def __hash__(self):
        return hash(tuple(self.list()))

    def __float__(self):
         if self.degree() > 0:
             raise TypeError, "cannot coerce nonconstant polynomial to float"
         return float(self[0])

    def __int__(self):
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial to int"
        return int(self[0])

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES:
            sage: R, x = PolynomialRing(ZZ).objgen()
            sage: H = Hom(R, QQ); H
            Set of Homomorphisms from Univariate Polynomial Ring in x over Integer Ring to Rational Field
            sage: f = H([5]); f
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn: x |--> 5
            sage: f(x)
            5
            sage: f(x^2 + 3)
            28
        """
        a = im_gens[0]
        P = a.parent()
        d = self.degree()
        result = P._coerce_(self[d])
        i = d - 1
        while i >= 0:
            result = result * a + P._coerce_(self[i])
            i -= 1
        return result

    def _integer_(self):
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial %s to int"%self
        return integer.Integer(self[0])

    def __invert__(self):
        return self.parent()(1)/self

    def inverse_of_unit(self):
        if self.degree() > 0:
            raise ValueError, "self is not a unit."
        return self.parent()(~(self[0]))

    def __long__(self):
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial to long"
        return long(self[0])

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: (x - 4)*(x^2 - 8*x + 16)
            x^3 - 12*x^2 + 48*x - 64
        """
        if right == 0 or self == 0:
            return self.polynomial(0)
        return self._mul_karatsuba(right)

    def __div__(self, right):
        """
        EXAMPLES:
            sage: x = QQ['x'].gen()
            sage: f = (x^3 + 5)/3; f
            1/3*x^3 + 5/3
            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field

        If we do the same over $\ZZ$ the result has to lie
        in the fraction field.

            sage: x  = ZZ['x'].gen()
            sage: f = (x^3 + 5)/3; f
            (x^3 + 5)/3
            sage: f.parent()
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent.
            sage: R, x = QQ['x'].objgen()
            sage: f = x^3 + 5
            sage: g = R(3)
            sage: h = f/g; h
            1/3*x^3 + 5/3
            sage: h.parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
        """
        try:
            if not isinstance(right, Element) or right.parent() != self.parent():
                R = self.parent().base_ring()
                x = R(right)
                return ~x * self
        except (TypeError, ValueError, ZeroDivisionError):
            pass
        return ring_element.RingElement.__div__(self, right)


    def _pow(self, right):
        if self.degree() <= 0:
            return self.parent()(self[0]**right)
        if right < 0:
            return (~self)**(-right)
        if self._is_gen:   # special case x**n should be faster!
            z = self.parent()(0)
            z[right] = 1
            return z
        return arith.generic_power(self, right)

    def _repr(self, name=None):
        s = " "
        m = self.degree() + 1
        r = reversed(xrange(m))
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x != 0:
                if n != m-1:
                    s += " + "
                x = str(x)
                if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if s==" ":
            return "0"
        return s[1:]

    def _repr_(self):
        return self._repr()

    def _latex_(self, name=None):
        s = " "
        m = self.degree() + 1
        r = reversed(xrange(m))
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x != 0:
                if n != m-1:
                    s += " + "
                x = latex.latex(x)
                if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "|%s^{%s}"%(name,n)
                elif n==1:
                    var = "|%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")
        if s==" ":
            return "0"
        return s[1:]


    def __setitem__(self, n, x):
        raise NotImplentedError


    def __floordiv__(self,right):
        """
        Quotient of division of self by other.  This is denoted //.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of self by other.
        EXAMPLES:
            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R

    def _is_atomic(self):
        return self.degree() == self.valuation()

    def _mul_generic(self, right):
        d1 = self.degree()
        d2 = right.degree()
        d = d1 + d2
        w = [sum([self[i]*right[k-i] for i in range(0,min(d1,k)+1) if \
                  i <= d1 and k-i <= d2 and self[i]!=0 and right[k-i]!=0]) \
                for k in range(d+1)]
        return Polynomial(self.parent(), w)

    def _mul_karatsuba(self, right):
        r"""
        Returns the product of two polynomials using the Karatsuba
        divide and conquer multiplication algorithm.  This is only
        used over a generic base ring.  (Special libraries like NTL
        are used, e.g., for the integers and rationals, which are much
        faster.)

        INPUT:
           self: Polynomial
           right: Polynomial (over same base ring as self)

        OUTPUT: Polynomial
           The product self*right.

        ALGORITHM:
           The basic idea is to use that
           $$
               (aX + b) (cX + d) = acX^2 + ((a+b)(c+d)-ac-bd)X + bd
           $$
           where ac=a*c and bd=b*d, which requires three
           multiplications instead of the naive four.  (In my examples,
           strangely just doing the above with four multiplications
           does tend to speed things up noticeably.)
           Given f and g of arbitrary degree bigger than one, let e
           be min(deg(f),deg(g))/2.  Write
           $$
                  f = a X^e + b   \text{ and }   g = c X^e + d
           $$
           and use the identity
           $$
                 (aX^e + b) (cX^e + d) = ac X^{2e} +((a+b)(c+d) - ac - bd)X^e + bd
           $$
           to recursively compute $fg$.

        TIMINGS:
        On a Pentium M 1.8Ghz laptop:
           f=R.random(1000,bound=100)
           g=R.random(1000,bound=100)
           time h=f._mul_karatsuba(g)
           Time: 0.42 seconds
           The naive multiplication algorithm takes 14.58 seconds.
           In contrast, MAGMA does this sort of product almost
           instantly, and can easily deal with degree 5000.  Basically
           MAGMA is 100 times faster at polynomial multiplication.

           Over Z using NTL, multiplying two polynomials constructed
           using R.random(10000,bound=100) takes 0.10 seconds.  Using
           MAGMA V2.11-10 the same takes 0.14 seconds.  So in this
           case NTL is somewhat faster than MAGMA.

           Over Q using PARI, multiplying two polynomials constructed
           using R.random(10000,bound=100) takes 1.23 seconds.  Not
           good!  TODO: use NTL polynomials over Z with a denominator
           instead of PARI.

        NOTES:
         * Karatsuba multiplication of polynomials is also implemented in PARI in
                src/basemath/polarit3.c
         * The MAGMA documentation appears to give no information about how
           polynomial multiplication is implemented.
        """

        def sum(v,w):
            if len(v)>=len(w):
                x = list(v)
                y = w
            else:
                x = list(w)
                y = v
            for i in range(len(y)):
                x[i] = x[i] + y[i]
            return x
        def dif(v,w):
            if len(v)>=len(w):
                x = list(v)
                y = w
            else:
                x = list(w)
                y = v
            for i in range(len(y)):
                x[i] -= y[i]
            return x
        def do_karatsuba(left, right):
            if len(left) == 0 or len(right) == 0:
                return []
            if len(left) == 1:
                return [left[0]*a for a in right]
            if len(right) == 1:
                return [right[0]*a for a in left]
            if len(left) == 2 and len(right) == 2:
                b = left[0]
                a = left[1]
                d = right[0]
                c = right[1]
                ac = a*c
                bd = b*d
                return [bd,(a+b)*(c+d)-ac-bd,ac]
            e = min(len(left), len(right))/2
            assert e>=1, "bug in karatsuba"
            a, b = left[e:], left[:e]
            c, d = right[e:], right[:e]
            ac = do_karatsuba(a,c)
            bd = do_karatsuba(b,d)
            zeros = [0 for _ in range(e)]
            t2 = zeros + zeros + ac
            t1 = zeros + dif(do_karatsuba(sum(a,b),sum(c,d)),sum(ac,bd))
            t0 = bd
            return sum(t0,sum(t1,t2))
        return self.parent()(do_karatsuba(self.list(), right.list()))

    def base_ring(self):
        """
        Return the base ring of the parent of self.

        EXAMPLES:
            sage: x = PolynomialRing(ZZ).gen()
            sage: x.base_ring()
            Integer Ring
            sage: (2*x+3).base_ring()
            Integer Ring
        """
        return self.parent().base_ring()

    def base_extend(self, R):
        """
        Return a copy of this polynomial but with coefficients in R.
        """
        S = sage.rings.polynomial_ring.PolynomialRing(R,
                                      name = self.parent().variable_name())
        return S(self)

    def copy(self):
        """
        Return a copy of self.

        EXAMPLES:
        We create the polynomial $f=x+3$, then set $g=f$, and change
        the coefficient of $x$ in $g$, which also changes the coefficient
        of $x$ in $f$.  If we instead copy $f$, then changing the
        coefficient of $x$ of $g$ does not change $f$.

            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: f = x+3
            sage: g = f
            sage: g[1]=3
            sage: f
            3*x + 3
            sage: g = f.copy()
            sage: g[1]=5
            sage: f
            3*x + 3
            sage: g
            5*x + 3
        """
        return self.polynomial(self)

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.

        EXAMPLES:
            sage: x = ZZ['x'].0
            sage: f = x^93 + 2*x + 1
            sage: f.degree()
            93
            sage: x = PolynomialRing(QQ, sparse=True).gen()
            sage: f = x^100000
            sage: f.degree()
            100000

            sage: x = QQ['x'].0
            sage: f = 2006*x^2006 - x^2 + 3
            sage: f.degree()
            2006
            sage: f = 0*x
            sage: f.degree()
            -1
            sage: f = x + 33
            sage: f.degree()
            1

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        raise NotImplementedError

    def denominator(self):
        """
        Return the least common multiple of the denominators of
        the entries of self, when this makes sense, i.e., when the
        coefficients have a denominator function.

        WARNING: This is not the denominator of the rational function
        defined by self, which would always be 1 since self is a polynomial.

        EXAMPLES:
        First we compute the denominator of a polynomial with integer
        coefficients, which is of course 1.
            sage: x = PolynomialRing(IntegerRing()).gen()
            sage: f = x^3 + 17*x + 1
            sage: f.denominator()
            1

        Next we compute the denominator of a polynomial with rational coefficients.
            sage: Q = RationalField()
            sage: x = PolynomialRing(Q).gen()
            sage: f = Q('1/17')*x^19 - Q('2/3')*x + Q('1/3'); f
            1/17*x^19 - 2/3*x + 1/3
            sage: f.denominator()
            51

        Finally, we try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements
        do not have a denominator method.
            sage: R = RealField()
            sage: x = PolynomialRing(R).gen()
            sage: f = x + R('0.3'); f
            1.0000000000000000*x + 0.29999999999999999
            sage: f.denominator()
            Traceback (most recent call last):
            ...
            AttributeError: 'mpfr.RealNumber' object has no attribute 'denominator'
        """
        if self.degree() == -1:
            return 1
        R = self.base_ring()
        x = self.list()
        d = x[0].denominator()
        for y in x:
            d = d.lcm(y.denominator())
        return d

    def derivative(self):
        return self.polynomial([self[n]*n for n in xrange(1,self.degree()+1)])

    def integral(self):
        try:
            return self.polynomial([0] + [self[n]/(n+1) for n in xrange(0,self.degree()+1)])
        except TypeError:
            raise ArithmeticError, "coefficients of integral of %s cannot be coerced into the base ring"%self


    def dict(self):
        X = {}
        Y = self.list()
        for i in xrange(len(Y)):
            X[i] = Y[i]
        return X

    def factor(self):
        r"""
        Return the factorization of self.

        INPUT:
            a polynomial

        OUTPUT:
            Factorization -- the factorization of self

        EXAMPLES:
        We factor some polynomials over $\Q$.
            sage: x = PolynomialRing(RationalField()).gen()
            sage: f = (x^3 - 1)^2
            sage: f.factor()
            (x - 1)^2 * (x^2 + x + 1)^2
            sage: f = 10*x^5 - 1
            sage: f.factor()
            (10*x^5 - 1)
            sage: f = 10*x^5 - 10
            sage: f.factor()
            (10) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)


        We factor a non-monic polynomial over the finite field $F_{25}$.
            sage: k, a = GF(25,'a').objgen()
            sage: R, x = PolynomialRing(k).objgen()
            sage: f = 2*x^10 + 2*x + 2*a
            sage: F = f.factor(); F
            (2) * (x + a + 2) * (x^2 + (a + 1)*x + a + 2) * (x^2 + 3*x + 4*a + 4) *
            (x^5 + (3*a + 4)*x^4 + (3*a + 3)*x^3 + 2*a*x^2 + (3*a + 1)*x + 3*a + 1)

        Notice that the unit factor is included when we multiply $F$ back out.
            sage: F.mul()
            2*x^10 + 2*x + 2*a

        Factorization also works even if the variable of the finite field is nefariously
        labeled "x".
            sage: R, x = PolynomialRing(GF(3^2, 'x')).objgen()
            sage: f = x^10 +7*x -13
            sage: f.factor()
            (x + 2*x + 1) * (x + x) * (x^4 + 2*x*x^3 + (x + 1)*x + 2) * (x^4 + (x + 2)*x^3 + (2*x + 2)*x + 2)
            sage: f.parent().base_ring().assign_names(['a'])
            sage: f.factor()
            (x + 2*a + 1) * (x + a) * (x^4 + 2*a*x^3 + (a + 1)*x + 2) * (x^4 + (a + 2)*x^3 + (2*a + 2)*x + 2)

            sage: k, a = GF(9,'x').objgen()
            sage: x = PolynomialRing(k,'x0').gen()
            sage: f = x^3 + x + 1
            sage: f.factor()
            (x0 + 2*x + 1) * (x0 + x) * (x0 + 2)

            sage: f = 0*x
            sage: f.factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

            sage: f = x^0
            sage: f.factor()
            1
        """

        # PERFORMANCE NOTE:
        #     In many tests with SMALL degree PARI is substantially
        #     better than NTL.  (And magma is better yet.)  And the
        #     timing difference has nothing to do with moving Python
        #     data to NTL and back.
        #     For large degree ( > 1500) in the one test I tried, NTL was
        #     *much* better than MAGMA, and far better than PARI.  So probably
        #     NTL's implementation is asymptotically better.  I could use
        #     PARI for smaller degree over other rings besides Z, and use
        #     NTL in general.

        R = self.parent().base_ring()
        if self.degree() < 0:
            raise ValueError, "factorization of 0 not defined"
        G = None

        from sage.rings.number_field.all import is_NumberField

        if integer_mod_ring.is_IntegerModRing(R) or finite_field.is_FiniteField(R) or \
               isinstance(R, (integer_ring.IntegerRing, rational_field.RationalField)):

            G = list(self._pari_('x').factor())

        elif is_NumberField(R):

            v = [x._pari_("a") for x in reversed(self.list())]
            f = pari(v).Pol()
            G = list(f.factor())


        if G is None:
            raise NotImplementedError
        return self._factor_pari_helper(G)

    def _factor_pari_helper(self, G, unit=None):
        pols = G[0]
        exps = G[1]
        F = []
        R = self.parent()
        c = R.base_ring()(1)
        for i in xrange(len(pols)):
            f = R(pols[i])
            e = int(exps[i])
            c *= f.leading_coefficient()
            F.append((f,e))
        if unit is None:
            unit = R.base_ring()(self.leading_coefficient()/c)
        if not unit.is_unit():
            F.append((R(unit), 1))
            unit = R.base_ring()(1)
        return factorization.Factorization(F, unit)

    def _lcm(self, other):
        """
        Let f and g be two polynomials.  Then this function
        returns the monic least common multiple of f and g.
        """
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q[q.degree()])*q  # make monic  (~ is inverse in python)

    def is_constant(self):
        return self.degree() <= 0

    def constant_coefficient(self):
        return self[0]

    def is_monic(self):
        """
        Returns True if this polynomial is monic.  The zero
        polynomial is by definition not monic.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = x + 33
            sage: f.is_monic()
            True
            sage: f = 0*x
            sage: f.is_monic()
            False
            sage: f = 3*x^3 + x^4 + x^2
            sage: f.is_monic()
            True
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.is_monic()
            False

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return not self.is_zero() and self[self.degree()] == 1

    def is_unit(self):
        if self.degree() > 0:
            return False
        return self[0].is_unit()

    def is_gen(self):
        return self._is_gen

    def is_zero(self):
        return self.degree() == -1

    def leading_coefficient(self):
        return self[self.degree()]

    def monic(self):
        """
        Return this polynomial divided by its leading coefficient.
        Does not change this polynomial.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.monic()
            x^5 + 1/56*x^3 + 1/28*x^2
            sage: f = (1/4)*x^2 + 3*x + 1
            sage: f.monic()
            x^2 + 12*x + 4

    The following happens because $f = 0$ cannot be made into a monic polynomial
            sage: f = 0*x
            sage: f.monic()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        Notice that the monic version of a polynomial over the
        integers is defined over the rationals.
            sage: x = ZZ['x'].0
            sage: f = 3*x^19 + x^2 - 37
            sage: g = f.monic(); g
            x^19 + 1/3*x^2 - 37/3
            sage: g.parent()
            Univariate Polynomial Ring in x over Rational Field


        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        if self.is_monic():
            return self
        a = ~self.leading_coefficient()
        R = self.parent()
        if a.parent() != R.base_ring():
            S = R.base_extend(a.parent())
            return a*S(self)
        else:
            return a*self


    def list(self):
        raise NotImplementedError

    def newton_raphson(self, n, x0):
        """
        Return a list of n iterative approximations to a root of this
        polynomial, computed using the Newton-Raphson method.

        The Newton-Raphson method is an iterative root-finding algorithm.
        For f(x) a polynomial, as is the case here, this is essentially
        the same as Horner's method.

        INPUT:
           n -- an integer (=the number of iterations),
           x0 -- an initial guess x0.

        OUTPUT:
           A list of numbers hopefully approximating a root of f(x)=0.

           ** If one of the iterates is a critical point of f then
              a ZeroDivisionError exception is raised.

        EXAMPLES:
            sage: x = PolynomialRing(RealField(), 'x').gen()
            sage: f = x^2 - 2
            sage: f.newton_raphson(4, 1)
            [1.5000000000000000, 1.4166666666666667, 1.4142156862745099, 1.4142135623746899]

        AUTHORS: David Joyner and William Stein (2005-11-28)
        """
        n = integer.Integer(n)
        df = self.derivative()
        def newton(z):
            return z -  self(z) / df(z)
        K = self.parent().base_ring()
        a = K(x0)
        L = []
        for i in range(n):
            a = newton(a)
            L.append(a)
        return L

    def polynomial(self, *args, **kwds):
        return self.parent()(*args, **kwds)

    def _pari_(self, variable=None):
        """
        Return polynomial as a PARI object.  Note that
        the variable will be "x" unless you explicitly specify
        otherwise, no matter what the polynomial indeterminate
        is.
        """
        try:
            return self.__pari
        except AttributeError:
            v = list(reversed(self.list()))
            try:
                v = [x._pari_() for x in v]
            except AttributeError:
                pass
            if variable is None:
                variable = self.parent().variable_name()
            self.__pari = pari(v).Pol(variable)
            return self.__pari

    def _pari_init_(self):
        return str(self._pari_())

    def resultant(self, other, flag=0):
        raise NotImplementedError

        ## This should be switched to use NTL, which can apparently compute
        ## resultants!
        ##        void XGCD(ZZ& r, ZZX& s, ZZX& t, const ZZX& a, const ZZX& b,
        ##          long deterministic=0);
        ##// r = resultant of a and b; if r != 0, then computes s and t such
        ##// that: a*s + b*t = r; otherwise s and t not affected.  if
        ##// !deterministic, then resultant computation may use a randomized
        ##// strategy that errs with probability no more than 2^{-80}.
        #m = magma.Magma()
        #cmd = "R<%s> := PolynomialRing(RationalField()); "%self.parent().variable_name() + \
        #      "Resultant(%s, %s);"%(self,other)
        #s = m.cmd(cmd)
        #i = s.find("\r")
        #return eval(s[:i])

    def reverse(self):
        v = list(self.list())
        v.reverse()
        return self.parent()(v)

    def roots(self):
        """
        Return all roots of this polynomial.

        EXAMPLES:
            sage: x = PolynomialRing(RationalField()).gen()
            sage: f = x^3 - 1
            sage: f.roots()
            [(1, 1)]
            sage: f = (x^3 - 1)^2
            sage: f.roots()
            [(1, 2)]

            sage: f = -19*x + 884736
            sage: f.roots()
            [(884736/19, 1)]
            sage: (f^20).roots()
            [(884736/19, 20)]
        """
        seq = []
        try:
            rts = self.factor()
        except NotImplementedError:
            raise NotImplementedError, "root finding for this polynomial not implemented"
        for fac in rts:
            g = fac[0]
            if g.degree() == 1:
                seq.append((-g[0]/g[1],fac[1]))
        return seq

    def valuation(self):
        r"""
        If $f = a_r x^r + a_{r+1}x^{r+1} + \cdots$, with $a_r$ nonzero,
        then the valuation of $f$ is $r$.  The valuation of the zero
        polynomial is $\infty$.
        """
        if self.is_zero():
            return infinity
        for i in xrange(self.degree()+1):
            if self[i] != 0:
                return i
        raise RuntimeError, "bug in computing valuation of polynomial"

    def name(self):
        return self.parent().variable_name()

    def _xgcd(self, other):
        r"""
        Extended gcd of self and polynomial other.

        Returns g, u, and v such that
              \code{g = u*self + v*other.}

        EXAMPLES:
            sage: P, x = PolynomialRing(QQ).objgen()
            sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (27*x^2 + 54, 1, -x^2 - 3*x - 9)
            sage: u*F + v*G
            27*x^2 + 54
            sage: x.xgcd(P(0))
            (1, 0, x)
            sage: f = P(0)
            sage: f.xgcd(x)
            (x, 0, 1)
        """
        if other.is_zero():
            R = self.parent()
            return R(1), R(0), self
        # Algorithm 3.2.2 of Cohen, GTM 138
        R = self.parent()
        A = self
        B = other
        U = R(1)
        G = A
        V1 = R(0)
        V3 = B
        while V3.is_nonzero():
            Q, R = G.quo_rem(V3)
            T = U - V1*Q
            U = V1
            G = V3
            V1 = T
            V3 = R
        V = (G-A*U)//B
        return G, U, V

    def is_irreducible(self):
        F = self.factor()
        if len(F) > 1 or F[0][1] > 1:
            return False
        return True


    def truncate(self, n):
        r"""
        Replace this polynomial by $\sum a_m x^m$ where the sum is
        over $m < n$.  The resulting polynomial is equivalent to self
        modulo $x^n$.
        """
        return self.parent()(self[:int(n)], check=False)


class Polynomial_generic_dense(Polynomial):
    """
    A generic dense polynomial.

    EXAMPLES:
        sage: R, x = PolynomialRing(PolynomialRing(Q)).objgen()
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_dense'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x == None:
            self.__coeffs = []
            return
        R = parent.base_ring()
#        if isinstance(x, Polynomial) and x.parent() == self.parent():
#            x = list(x.list())
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = list(x.list())
            elif x.parent() == R:
                x = [x]
            else:
                x = [R(a) for a in x.list()]
                check = False
                #raise TypeError, "Cannot coerce %s into %s."%(x, parent)
        elif isinstance(x, dict):
            zero = R(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v
        elif isinstance(x, gen):
            x = [R(w) for w in reversed(x.Vec())]
            check = True
        elif not isinstance(x, list):
            x = [x]   # constant polynomials
        if check:
            self.__coeffs = [R(z) for z in x]
        else:
            self.__coeffs = x
        if check:
            self.__normalize()

    def __normalize(self):
        x = self.__coeffs
        n = len(x)-1
        while n>=0 and x[n] == 0:
            del x[n]
            n -= 1

    def __getitem__(self,n):
        if n < 0 or n >= len(self.__coeffs):
            return self.base_ring()(0)
        return self.__coeffs[n]

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        v = self.__coeffs[i:j]
        zero = self.base_ring()(0)
        for k in xrange(len(self.__coeffs),j):
            v.append(zero)
        return v

    def __setitem__(self, n, value):
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        value = self.base_ring()(value)
        if n >= 0 and n < len(self.__coeffs):
            self.__coeffs[n] = value
            if n == len(self.__coeffs) and value == 0:
                self.__normalize()
        elif n < 0:
            raise IndexError, "polynomial coefficient index must be nonnegative"
        elif value != 0:
            zero = self.base_ring()(0)
            for _ in xrange(len(self.__coeffs), n):
                self.__coeffs.append(zero)
            self.__coeffs.append(value)

    def list(self):
        return self.__coeffs

    def degree(self):
        return len(self.__coeffs) - 1


class Polynomial_generic_sparse(Polynomial):
    """
    A generic sparse polynomial.

    EXAMPLES:
        sage: R, x = PolynomialRing(PolynomialRing(Q), sparse=True).objgen()
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_sparse'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x == None:
            self.__coeffs = {}
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = dict(x.dict())
            elif x.parent() == R:
                x = {0:x}
            else:
                w = {}
                for n, c in x.dict().iteritems():
                    w[n] = R(c)
                #raise TypeError, "Cannot coerce %s into %s."%(x, parent)
        elif isinstance(x, list):
            y = {}
            for i in xrange(len(x)):
                if x[i] != 0:
                    y[i] = x[i]
            x = y
        elif not isinstance(x, dict):
            x = {0:x}   # constant polynomials
        elif isinstance(x, gen):
            x = [R(w) for w in reversed(x.Vec())]
            check = True
        if check:
            self.__coeffs = {}
            for i, z in x.iteritems():
                self.__coeffs[i] = R(z)
        else:
            self.__coeffs = x
        if check:
            self.__normalize()


    def __normalize(self):
        x = self.__coeffs
        zero = self.base_ring()(0)
        D = [n for n, z in x.iteritems() if z == 0]
        for n in D:
            del x[n]

    def __getitem__(self,n):
        if not self.__coeffs.has_key(n):
            return self.base_ring()(0)
        return self.__coeffs[n]

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        zero = self.base_ring()(0)
        v = [zero for _ in xrange(i,j)]
        x = self.__coeffs
        for k in set(x.keys()).intersection(set(xrange(i,j))):
            v[k] = x[k]
        return v

    def __setitem__(self, n, value):
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        value = self.base_ring()(value)
        x = self.__coeffs
        if n < 0:
            raise IndexError, "polynomial coefficient index must be nonnegative"
        if value == 0:
            if x.has_key(n):
                del x[n]
        else:
            x[n] = value

    def list(self):
        zero = self.base_ring()(0)
        v = [zero for _ in xrange(self.degree()+1)]
        for n, x in self.__coeffs.iteritems():
            v[n] = x
        return v

    #def _pari_(self, variable=None):
    #    if variable is None:
    #        return self.__pari
    #    else:
    #        return self.__pari.subst('x',variable)

    def degree(self):
        v = self.__coeffs.keys()
        if len(v) == 0:
            return -1
        return max(v)


class Polynomial_generic_field(Polynomial, euclidean_domain_element.EuclideanDomainElement):
    def __init__(self, parent, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

    def quo_rem(self, other):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.

        EXAMPLES:
            sage: K = NumberField(x**ZZ(2)-ZZ(2),'t')
            sage: P, x = PolynomialRing(K).objgen()
            sage: x.quo_rem(K(1))
            (x, 0)
            sage: x.xgcd(K(1))
            (1, 0, 1)
        """
        other = self.parent()(other)
        if other.is_zero():
            raise ZeroDivisionError, "other (=%s) must be nonzero"%other

        # This is algorithm 3.1.1 in Cohen GTM 138
        A = self
        B = other
        R = A
        Q = self.polynomial(0)
        X = self.parent().gen()
        while R.degree() >= B.degree():
            S =  (R.leading_coefficient()/B.leading_coefficient()) * X**(R.degree()-B.degree())
            Q += S
            R -= S*B
        return (Q, R)

    def _gcd(self, other):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        g = euclidean_domain_element.EuclideanDomainElement._gcd(self, other)
        c = g.leading_coefficient()
        if c.is_unit():
            return (1/c)*g
        return g


class Polynomial_generic_sparse_field(Polynomial_generic_sparse, Polynomial_generic_field):
    """
    EXAMPLES:
        sage: R, x = PolynomialRing(RealField(), sparse=True).objgen()
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_sparse_field'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_sparse.__init__(self, parent, x, check, is_gen)


class Polynomial_generic_dense_field(Polynomial_generic_dense, Polynomial_generic_field):
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_dense.__init__(self, parent, x, check, is_gen)


class Polynomial_rational_dense(Polynomial_generic_field):
    """
    A dense polynomial over the rational numbers.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

        if construct:
            self.__poly = x
            return

        self.__poly = pari([]).Pol()

        if x is None:
            return         # leave initialized to 0 polynomial.


        if fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError, "denominator (=%s) must be 1"%x.denominator()
            else:
                x = x.numerator()

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = x.__poly.copy()
                return
            else:
                x = [QQ(a) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = QQ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, gen):
            f = x.Pol()
            self.__poly = f
            assert self.__poly.type() == "t_POL"
            return

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            x = [QQ(z) for z in x]

        self.__list = list(x)
        while len(self.__list) > 0 and self.__list[-1] == 0:
            del self.__list[-1]

        x.reverse()

        # NOTE: It is much faster to convert to string and let pari's parser at it,
        # which is why we pass str(x) in.
        self.__poly = pari(str(x)).Pol()
        assert self.__poly.type() == "t_POL"

    def _repr(self, name=None):
        if name is None:
            name = self.parent().variable_name()
        return str(self.__poly).replace("x", name)

    def _repr_(self):
        return self._repr()

    def __reduce__(self):
        return Polynomial_rational_dense, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        return QQ(self.__poly[n])

    def __getslice__(self, i, j):
        return [QQ(self.__poly[k]) for k in range(i,j)]

    def _pow(self, n):
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return Polynomial_rational_dense(self.parent(), self.__poly**n, construct=True)

    def _add_(self, right):
        return Polynomial_rational_dense(self.parent(),
                                         self.__poly + right.__poly, construct=True)

    def is_irreducible(self):
        try:
            return self.__poly.polisirreducible()
        except NotImplementedError:
            F = self.__poly.factor()
            if len(F) > 1 or F[0][1] > 1:
                return False
            return True

    def galois_group(self, pari_group=False, use_kash=False):
        r"""
        Return the Galois group of f as a permutation group.

        INPUT:
            self -- an irreducible polynomial

            pari_group -- bool (default: False); if True instead return
                          the Galois group as a PARI group.  This has
                          a useful label in it, and may be slightly faster
                          since it doesn't require looking up a group in
                          Gap.  To get a permutation group from a PARI
                          group P, type PermutationGroup(P).

            use_kash --   bool (default: False); if True use KASH's Galois
                          command instead of using the PARI C library.
                          An attempt is always made to use KASH if the
                          degree of the polynomial is >= 12.

        ALGORITHM: The Galois group is computed using PARI in C
        library mode, or possibly kash if available.

        \note{ The PARI documentation contains the following warning:
        The method used is that of resolvent polynomials and is
        sensitive to the current precision. The precision is updated
        internally but, in very rare cases, a wrong result may be
        returned if the initial precision was not sufficient.}

        EXAMPLES:
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(); G            # uses optional database_gap package
            Transitive group number 5 of degree 4
            sage: G.gens()                           # uses optional database_gap package
            ((1,2,3,4), (1,2))
            sage: G.order()                          # uses optional database_gap package
            24

        It is potentially useful to instead obtain the corresponding
        PARI group, which is little more than a $4$-tuple.  See the
        PARI manual for the exact details.  (Note that the third
        entry in the tuple is in the new standard ordering.)
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(pari_group=True); G
            PARI group [24, -1, 5, "S4"] of degree 4
            sage: PermutationGroup(G)                # uses optional database_gap package
            Transitive group number 5 of degree 4

        You can use KASH to compute Galois groups as well.  The
        avantage is that KASH can compute Galois groups of fields up
        to degree 23, whereas PARI only goes to degree 11.  (In my
        not-so-thorough experiments PARI is faster than KASH.)

            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(use_kash=true)      # requires optional KASH
            Transitive group number 5 of degree 4

        """
        from sage.groups.all import PariGroup, PermutationGroup, TransitiveGroup
        if not self.is_irreducible():
            raise ValueError, "polynomial (=%s) must be irreducible"%self
        if self.degree() > 11 or use_kash:
            # TODO -- maybe use KASH if available or print message that user should install KASH?
            try:
                from sage.interfaces.all import kash
                s = self._repr(name='X')
                G = kash('Galois(%s)'%s)
                d = int(kash.eval('%s.ext1'%G.name()))
                n = int(kash.eval('%s.ext2'%G.name()))
                return TransitiveGroup(d, n)
            except RuntimeError:
                raise NotImplementedError, "Sorry, computation of Galois groups of fields of degree bigger than 11 is not yet implemented.  Try installing the optional free (closed source) KASH package, which supports up to degree $23$."
        G = self.__poly.polgalois()
        H = PariGroup(G, self.degree())
        if pari_group:
            return H
        else:
            return PermutationGroup(H)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_rational_dense):
            right = self.parent()(right)
        if right.parent() != self.parent():
            raise TypeError
        v = self.__poly.divrem(right.__poly)
        return Polynomial_rational_dense(self.parent(), v[0], construct=True), \
               Polynomial_rational_dense(self.parent(), v[1], construct=True)


    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ).gen()
            sage: (x - QQ('2/3'))*(x^2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def __setitem__(self, n, value):
        try:
            del self.__list
        except AttributeError:
            pass
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        if n < 0:
            raise IndexError, "n (=%s) must be >= 0"%n
        if n <= self.__poly.poldegree():
            self.__poly[n] = QQ(value)
        else:
            self.__poly = self.__poly + pari('(%s)*x^%s'%(QQ(value),n))
        if hasattr(self, "__list"):
            del self.__list

    def complex_roots(self, flag=0):
        """
        Returns the complex roots of this polynomial.
        INPUT:
            flag -- optional, and can be
                    0: (default), uses Schonhage's method modified by Gourdon,
                    1: uses a modified Newton method.
        OUTPUT:
            list of complex roots of this polynomial, counted with multiplicities.

        NOTE: Calls the pari function polroots.

        EXAMPLE:
        We compute the roots of the characteristic polynomial of some Salem numbers:
            sage: R = PolynomialRing(QQ); x = R.gen()
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.complex_roots()[0]
            0.71363917353690087
        """
        R = self.__poly.polroots(flag)
        C = complex_field.CC
        return [C(a) for a in R]

    def copy(self):
        f = Polynomial_rational_dense(self.parent())
        f.__poly = self.__poly.copy()
        return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.poldegree(), -1)

    def discriminant(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ).gen()
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            -7911
        """
        return QQ(self.__poly.poldisc())

    def disc(self):
        """
        Same as discriminant().
        """
        return self.discriminant()

    def factor_mod(self, p):
        """
        Return the factorization of self modulo the prime p.

        INPUT:
            p -- prime

        OUTPUT:
            factorization of self reduced modulo p.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) must be prime"%p
        G = self._pari_().factormod(p)
        K = finite_field.FiniteField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, name=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, unit=K(self.leading_coefficient()))

    def factor_padic(self, p, prec=10):
        """
        Return p-adic factorization of self to given precision.

        INPUT:
            p -- prime
            prec -- integer; the precision

        OUTPUT:
            factorization of self reduced modulo p.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) must be prime"%p
        prec = integer.Integer(prec)
        if prec <= 0:
            raise ValueError, "prec (=%s) must be positive"%prec
        G = self._pari_().factorpadic(p, prec)
        K = padic_field.pAdicField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, name=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, K(self.leading_coefficient()))

    def list(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ).gen()
            sage: f = x^3 + 3*x - QQ('17/13')
            sage: f.list()
            [-17/13, 3, 0, 1]
        """
        try:
            return self.__list
        except AttributeError:
            self.__list = [QQ(x) for x in reversed(self.__poly.Vec())]
            return self.__list

##     def partial_fraction(self, g):
##         """
##         Return partial fraction decomposition of self/g, where g
##         has the same parent as self.
##         """
##         g = self.parent()(g)
##         from sage.interfaces.maxima import maxima
##         h = maxima(self)/maxima(g)
##         k = h.partfrac(self.parent().variable())

    def rescale(self, a):
        """
        Return f(a*X).
        """
        b = 1
        c = []
        for i in range(self.degree()+1):
            c.append(self[i]*b)
            b *= a
        return self.parent()(c)

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        NOTES:
            Implemented using pari's polresultant function.

        EXAMPLES:
            sage: x = PolynomialRing(RationalField()).gen()
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: f.resultant(g)
            -8
        """
        if not isinstance(other, Polynomial):
            other = self.polynomial(other)
        if self.parent() != other.parent():
            raise TypeError
        return QQ(str(self.__poly.polresultant(other.__poly, 0)))

    def hensel_lift(self, p, e):
        """
        Assuming that self factors modulo $p$ into distinct factors,
        computes the Hensel lifts of these factors modulo $p^e$.  We
        assume that $p$ has integer coefficients.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) must be prime"%p
        e = integer.Integer(e)
        if e < 1:
            raise ValueError, "e (=%s) must be at least 1"%e
        F = self.factor_mod(p)
        y = []
        for g, n in F:
            if n > 1:
                raise ArithmeticError, "The polynomial must be square free modulo p."
            y.append(g)
        H = self._pari_().polhensellift(y, p, e)
        R = integer_mod_ring.IntegerModRing(p**e)
        S = sage.rings.polynomial_ring.PolynomialRing(R, self.parent().variable_name())
        return [S(eval(str(m.Vec().Polrev().Vec()))) for m in H]

class Polynomial_integer_dense(Polynomial, integral_domain_element.IntegralDomainElement):
    """
    A dense polynomial over the integers.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if construct:
            if isinstance(x, ZZX_class):
                self.__poly = x
                return
            self.__poly = ZZX(x)
            return

        self.__poly = ZZX([])

        if x == None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = x.__poly.copy()
                return
            else:
                x = [ZZ(a) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = ZZ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, ZZX_class):
            self.__poly = x.copy()
            return

        elif isinstance(x, gen):
            x = list(reversed(x.list()))
            check = False

        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_integer_dense):
            if x.denominator() == 1:
                x = x.numerator().__poly
                check = False

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            x = [ZZ(z) for z in x]

        self.__poly = ZZX(x)

    def content(self):
        """
        Return the greatest common divisor of the coefficients of this
        polynomial.
        """
        return ZZ(self.__poly.content())

    def ntl_ZZX(self):
        """
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality may be available through
        this function.
        """
        return self.__poly

    def __reduce__(self):
        return Polynomial_integer_dense, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        return ZZ(self.__poly[n])

    def __getslice__(self, i, j):
        return [ZZ(self.__poly[k]) for k in range(i,j)]

    def _pow(self, n):
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    def _add_(self, right):
        return self.parent()(self.__poly + right.__poly, construct=True)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        v = self.__poly.quo_rem(right.__poly)
        return self.parent()(v[0], construct=True), \
               self.parent()(v[1], construct=True)

    def gcd(self, right):
        """
        Return the GCD of self and other.  The leading
        coefficient need not be 1.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.gcd(right.__poly)
        return self.parent()(g, construct=True)

    def lcm(self, right):
        """
        Return the LCM of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.lcm(right.__poly)
        return self.parent()(g, construct=True)

    def xgcd(self, right):
        """
        Return $g, u, v$ such that \code{g = u*self + v*right}.

        If self and right are coprime as polynomials over the
        rationals, then $g$ is guaranteed to be the resultant of self
        and right, as a constant polynomial.

        EXAMPLES:
            sage: P, x = PolynomialRing(ZZ).objgen()
            sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (27*x^2 + 54, 1, -x^2 - 3*x - 9)
            sage: u*F + v*G
            27*x^2 + 54
            sage: x.xgcd(P(0))
            (1, 0, x)
            sage: f = P(0)
            sage: f.xgcd(x)
            (x, 0, 1)
            sage: F = (x-3)^3; G = (x-15)^2
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (2985984, -432*x + 8208, 432*x^2 + 864*x + 14256)
            sage: u*F + v*G
            2985984
        """
        r, s, t = self.ntl_ZZX().xgcd(right.ntl_ZZX())
        K = self.base_ring()
        rr = K(str(r))   # optimize in future
        if rr == 0:
            QQ = sage.rings.rational_field.QQ
            f = self.base_extend(QQ)
            g, u, v = f.xgcd(right.base_extend(QQ))
            d = arith.lcm([g.denominator(), u.denominator(), v.denominator()])
            R = self.parent()
            return R(d*g), R(d*u), R(d*v)
        else:
            S = self.parent()
            return S(rr), S(s, construct=True), \
                   S(t, construct=True)


    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(ZZ).gen()
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 - 10*x^2 + 32*x - 32
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def __setitem__(self, n, value):
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        if n < 0:
            raise IndexError, "n (=%s) must be >= 0"%n
        self.__poly[n] = int(value)

    def complex_roots(self, flag=0):
        """
        Returns the complex roots of this polynomial.
        INPUT:
            flag -- optional, and can be
                    0: (default), uses Schonhage's method modified by Gourdon,
                    1: uses a modified Newton method.
        OUTPUT:
            list of complex roots of this polynomial, counted with multiplicities.

        NOTE: Calls the pari function polroots.

        EXAMPLE:
        We compute the roots of the characteristic polynomial of some Salem numbers:
            sage: R = PolynomialRing(ZZ); x = R.gen()
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.complex_roots()[0]    # todo: known bug in PARI 2.2.10 !!
            0.71363917353690087
        """
        QQ = sage.rings.rational_field.RationalField()
        R = sage.rings.polynomial_ring.PolynomialRing(QQ)
        return R(self.list()).complex_roots()

    def copy(self):
        f = Polynomial_integer_dense(self.parent())
        f.__poly = self.__poly.copy()
        return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def discriminant(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(ZZ).gen()
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            -7911
        """
        return ZZ(str(self.__poly.discriminant()))

    def _pari_(self, variable='x'):
        return pari(str(self.__poly).replace(' ',',')).Pol(variable).polrecip()

    def factor_mod(self, p):
        """
        Return the factorization of self modulo the prime p.

        INPUT:
            p -- prime

        OUTPUT:
            factorization of self reduced modulo p.

        EXAMPLES:
            sage: x = Z['x'].0
            sage: f = -3*x*(x-2)*(x-9) + x
            sage: f.factor_mod(3)
            0
            sage: f = -3*x*(x-2)*(x-9)
            sage: f.factor_mod(3)
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

            sage: f = 2*x*(x-2)*(x-9)
            sage: f.factor_mod(7)
            (2) * (x + 5)^2
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) must be prime"%p
        f = self._pari_()
        if f * pari('Mod(1,%s)'%p) == pari(0):
            raise ValueError, "factorization of 0 not defined"
        G = f.factormod(p)
        k = finite_field.FiniteField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(k, name=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, unit=k(self.leading_coefficient()))


    def factor_padic(self, p, prec=10):
        """
        Return p-adic factorization of self to given precision.

        INPUT:
            p -- prime
            prec -- integer; the precision

        OUTPUT:
            factorization of self reduced modulo p.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p (=%s) must be prime"%p
        prec = integer.Integer(prec)
        if prec <= 0:
            raise ValueError, "prec (=%s) must be positive"%prec
        G = self._pari_().factorpadic(p, prec)
        K = padic_field.pAdicField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, name=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, K(self.leading_coefficient()))

    def list(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(ZZ).gen()
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [-17, 3, 0, 1]
        """
        return [ZZ(str(self.__poly[i])) for i in xrange(self.degree()+1)]

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        NOTES:
            Implemented using NTL's polresultant function.

        EXAMPLES:
            sage: x = PolynomialRing(ZZ).gen()
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: f.resultant(g)
            -8
        """
        if not isinstance(other, Polynomial) or self.parent() != other.parent():
            other = self.polynomial(other)
        return ZZ(str(self.__poly.resultant(other.__poly, 0)))

    def ntl_set_directly(self, v):
        """
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers are stored internally using NTL's ZZX
        class.  Use this function to set the value of this polynomial using
        the NTL constructor, which is potentially quicker.   The input v
        is either a vector of ints or a string of the form '[ n1 n2 n3 ... ]'
        where the ni are integers and there are no commas between them.
        The optimal input format is the string format, since that's what NTL uses.

        EXAMPLES:
            sage: R = PolynomialRing(ZZ)
            sage: R([1,2,3])
            3*x^2 + 2*x + 1
            sage: f = R(0)
            sage: f.ntl_set_directly([1,2,3])
            sage: f
            3*x^2 + 2*x + 1
            sage: f.ntl_set_directly('[1 2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
        """
        if self._is_gen:
            raise TypeError, "Cannot change the value of the generator."
        self.__poly = ZZX(v)
        try:
            del self.__list
        except AttributeError:
            pass



class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite.

    EXAMPLES:
        sage: R, x = PolynomialRing(Integers(16)).objgen()
        sage: f = x^3 - x + 17
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

        if construct:
            if isinstance(x, ZZ_pX_class):
                self.__poly = x
                return
            parent._ntl_set_modulus()
            self.__poly = ZZ_pX(x)
            return

        self.__poly = ZZ_pX([])

        if x == None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                parent._ntl_set_modulus()
                self.__poly = x.__poly.copy()
                return
            else:
                x = [ZZ._coerce_(a) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = ZZ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, ZZX_class):
            self.__poly = x.copy()
            return

        elif isinstance(x, gen) and x.type() == 't_POL':
            x = list(reversed(eval(str(x.lift().list()))))
            check = False

        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator() == 1:
                x = x.numerator().__poly
                check = False

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            x = [ZZ(z) for z in x]

        parent._ntl_set_modulus()
        self.__poly = ZZ_pX(x)

    def __reduce__(self):
        return Polynomial_dense_mod_n, \
               (self.parent(), self.list(), False, self.is_gen())

    def int_list(self):
        return eval(str(self.__poly).replace(' ',','))

    def _pari_(self, variable='x'):
        return pari(str(list(reversed(self.int_list())))).Pol(variable) * \
               pari('Mod(1,%s)'%self.parent().base_ring().order())

    def ntl_ZZ_pX(self):
        """
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality is available through this
        function.
        """
        return self.__poly

    def __getitem__(self, n):
        return self.base_ring()(self.__poly[n])

    def __getslice__(self, i, j):
        R = self.base_ring()
        return [R(self.__poly[k]) for k in range(i,j)]

    def _pow(self, n):
        n = int(n)
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    def _add_(self, right):
        return self.parent()(self.__poly + right.__poly, construct=True)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_dense_mod_n):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        v = self.__poly.quo_rem(right.__poly)
        P = self.parent()
        return P(v[0], construct=True), P(v[1], construct=True)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(Integers(100)).gen()
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def __setitem__(self, n, value):
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        if n < 0:
            raise IndexError, "n (=%s) must be >= 0"%n
        self.parent()._ntl_set_modulus()
        self.__poly[n] = int(value)

    def copy(self):
        self.parent()._ntl_set_modulus()
        f = self.parent()()
        f.__poly = self.__poly.copy()
        return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def is_irreducible(self):
        return bool(self._pari_().polisirreducible())

    def list(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(Integers(100)).gen()
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [83, 3, 0, 1]
        """
        R = self.base_ring()
        return [R(x) for x in self.int_list()]

    def ntl_set_directly(self, v):
        r"""
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers modulo n are stored internally
        using NTL's ZZ_pX class.  Use this function to set the value
        of this polynomial using the NTL constructor, which is
        potentially \emph{very} fast.  The input v is either a vector
        of ints or a string of the form '[ n1 n2 n3 ... ]' where the
        ni are integers and there are no commas between them.  The
        optimal input format is the string format, since that's what
        NTL uses by default.

        EXAMPLES:
            sage: R = PolynomialRing(Integers(100))
            sage: R([1,-2,3])
            3*x^2 + 98*x + 1
            sage: f = R(0)
            sage: f.ntl_set_directly([1,-2,3])
            sage: f
            3*x^2 + 98*x + 1
            sage: f.ntl_set_directly('[1 -2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 98*x + 1
        """
        if self._is_gen:
            raise TypeError, "Cannot change the value of the generator."
        self.parent()._ntl_set_modulus()
        self.__poly = ZZ_pX(v)
        try:
            del self.__list
        except AttributeError:
            pass

class Polynomial_dense_mod_p(Polynomial_dense_mod_n, principal_ideal_domain_element.PrincipalIdealDomainElement):
    """
    A dense polynomial over the integers modulo p, where p is prime.
    """
    def __reduce__(self):
        return Polynomial_dense_mod_p, \
               (self.parent(), self.list(), False, self.is_gen())

    def _gcd(self, right):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_dense_mod_p):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.ntl_ZZ_pX().gcd(right.ntl_ZZ_pX())
        return self.parent()(g, construct=True)

    def _xgcd(self, right):
        """
        Return $g, u, v$ such that \code{g = u*self + v*right}.
        """
        r, s, t = self.ntl_ZZ_pX().xgcd(right.ntl_ZZ_pX())
        return self.parent()(r, construct=True), self.parent()(s, construct=True), \
               self.parent()(t, construct=True)

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        EXAMPLES:
            sage: x = PolynomialRing(GF(19)).gen()
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: f.resultant(g)
            11
        """
        if not isinstance(other, Polynomial) or self.parent() != other.parent():
            other = self.polynomial(other)
        self.parent()._ntl_set_modulus()
        return self.base_ring()(str(self.ntl_ZZ_pX().resultant(other.ntl_ZZ_pX())))

    def discriminant(self):
        """
        EXAMPLES:
            sage: x = PolynomialRing(GF(19)).gen()
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            12
        """
        self.parent()._ntl_set_modulus()
        return self.base_ring()(str(self.ntl_ZZ_pX().discriminant()))

    # PARI is way better than NTL for poly factor, and is called by default in the base class.
    #def factor(self, verbose=False):
    #    M = self.monic()
    #    self.parent()._ntl_set_modulus()
    #    F = [(self.parent()(f, construct=True), n) for f, n in M.ntl_ZZ_pX().factor(verbose)]
    #    return factorization.Factorization(F)


