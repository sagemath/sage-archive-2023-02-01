"""
Fraction Field Elements

AUTHOR: William Stein (input from David Joyner, David Kohel, and Joe Wetherell)
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import sage.rings.field_element as field_element
import integer_ring
from integer_ring import ZZ
from rational_field import QQ

import sage.misc.latex as latex
from sage.misc.misc import prod
from sage.misc.derivative import multi_derivative

def is_FractionFieldElement(x):
    """
    Returns whether or not x is of type FractionFieldElement

    EXAMPLES:
        sage: from sage.rings.fraction_field_element import is_FractionFieldElement
        sage: R.<x> = ZZ[]
        sage: is_FractionFieldElement(x/2)
        False
        sage: is_FractionFieldElement(2/x)
        True
        sage: is_FractionFieldElement(1/3)
        False
    """
    return isinstance(x, FractionFieldElement)

class FractionFieldElement(field_element.FieldElement):
    """
    EXAMPLES:
        sage: K, x = FractionField(PolynomialRing(QQ, 'x')).objgen()
        sage: K
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: loads(K.dumps()) == K
        True
        sage: f = (x^3 + x)/(17 - x^19); f
        (x^3 + x)/(-x^19 + 17)
        sage: loads(f.dumps()) == f
        True
    """

    def __init__(self, parent, numerator, denominator=1,
                 coerce=True, reduce=True):
        """
        EXAMPLES:
            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: K.<x> = Frac(ZZ['x'])
            sage: FractionFieldElement(K, x, 4)
            x/4
            sage: FractionFieldElement(K, x, x, reduce=False)
            x/x
            sage: f = FractionFieldElement(K, 'hi', 1, coerce=False, reduce=False)
            sage: f.numerator()
            'hi'
        """
        field_element.FieldElement.__init__(self, parent)
        if coerce:
            self.__numerator = parent.ring()(numerator)
            self.__denominator = parent.ring()(denominator)
        else:
            self.__numerator = numerator
            self.__denominator = denominator
        if reduce and parent.is_exact():
            try:
                self.reduce()
            except ArithmeticError:
                pass
        if self.__denominator.is_zero():
            raise ZeroDivisionError, "fraction field element division by zero"

    def reduce(self):
        """
        Divides out the gcd of the numerator and denominator.

        Automatically called for exact rings, but because it may be
        numerically unstable for inexact rings it must be called
        manually in that case.

        EXAMPLES:
            sage: R.<x> = RealField(10)[]
            sage: f = (x^2+2*x+1)/(x+1); f
            (1.0*x^2 + 2.0*x + 1.0)/(1.0*x + 1.0)
            sage: f.reduce(); f
            1.0*x + 1.0
        """
        try:
            g = self.__numerator.gcd(self.__denominator)
            if g != 1:
                numer, _ = self.__numerator.quo_rem(g)
                denom, _ = self.__denominator.quo_rem(g)
            else:
                numer = self.__numerator
                denom = self.__denominator
            if denom != 1 and denom.is_unit():
                try:
                    numer *= denom.inverse_of_unit()
                    denom = denom.parent()(1)
                except:
                    pass
            self.__numerator = numer; self.__denominator = denom
        except AttributeError, s:
            raise ArithmeticError, "unable to reduce because lack of gcd or quo_rem algorithm"
        except TypeError, s:
            raise ArithmeticError, "unable to reduce because gcd algorithm doesn't work on input"
        except NotImplementedError, s:
            raise ArithmeticError, "unable to reduce because gcd algorithm not implemented on input"

    def copy(self):
        """
        EXAMPLES:
            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.copy()
            (x + y)/y
        """
        return FractionFieldElement(self.parent(), self.__numerator,
                                    self.__denominator, coerce=False, reduce=False)

    def numerator(self):
        """
        EXAMPLES:
            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.numerator()
            x + y
        """
        return self.__numerator

    def denominator(self):
        """
        EXAMPLES:
            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.denominator()
            y
        """
        return self.__denominator

    def __hash__(self):
        """
        This function hashes in a special way to ensure that generators of a ring R
        and generators of a fraction field of R have the same hash.  This enables them
        to be used as keys interchangably in a dictionary (since \code{==} will claim them equal).
        This is particularly useful for methods like subs on \code{ParentWithGens} if you
        are passing a dictionary of substitutions.

        EXAMPLES:
            sage: R.<x>=ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: ((x+1)/(x^2+1)).subs({x:1})
            1
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x>=QQ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x,y,z>=ZZ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x,y,z>=QQ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: ((x+1)/(x^2+1)).subs({x:1})
            1
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: hash(R(1)/R(2))==hash(1/2)
            True
        """
        # This is same algorithm as used for members of QQ
        #cdef long n, d
        n = hash(self.__numerator)
        d = hash(self.__denominator)
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

    def partial_fraction_decomposition(self):
        """
        Decomposes fraction field element into a whole part and
        a list of fraction field elements over prime power denominators.

        The sum will be equal to the original fraction.

        AUTHOR:
             -- Robert Bradshaw (2007-05-31)

        EXAMPLES:
            sage: S.<t> = QQ[]
            sage: q = 1/(t+1) + 2/(t+2) + 3/(t-3); q
            (6*t^2 + 4*t - 6)/(t^3 - 7*t - 6)
            sage: whole, parts = q.partial_fraction_decomposition(); parts
            [3/(t - 3), 1/(t + 1), 2/(t + 2)]
            sage: sum(parts) == q
            True
            sage: q = 1/(t^3+1) + 2/(t^2+2) + 3/(t-3)^5
            sage: whole, parts = q.partial_fraction_decomposition(); parts
            [1/3/(t + 1), 3/(t^5 - 15*t^4 + 90*t^3 - 270*t^2 + 405*t - 243), (-1/3*t + 2/3)/(t^2 - t + 1), 2/(t^2 + 2)]
            sage: sum(parts) == q
            True

        We do the best we can over in-exact fields.
            sage: R.<x> = RealField(20)[]
            sage: q = 1/(x^2 + 2)^2 + 1/(x-1); q
            (1.0000*x^4 + 4.0000*x^2 + 1.0000*x + 3.0000)/(1.0000*x^5 - 1.0000*x^4 + 4.0000*x^3 - 4.0000*x^2 + 4.0000*x - 4.0000)
            sage: whole, parts = q.partial_fraction_decomposition(); parts
            [1.0000/(1.0000*x - 1.0000), 1.0000/(1.0000*x^4 + 4.0000*x^2 + 4.0000)]
            sage: sum(parts)
            (1.0000*x^4 + 4.0000*x^2 + 1.0000*x + 3.0000)/(1.0000*x^5 - 1.0000*x^4 + 4.0000*x^3 - 4.0000*x^2 + 4.0000*x - 4.0000)
        """
        denom = self.denominator()
        whole, numer = self.numerator().quo_rem(denom)
        factors = denom.factor()
        if factors.unit() != 1:
            numer *= ~factors.unit()
        if not self.parent().is_exact():
            # factors not grouped in this case
            # TODO: think about changing the factor code itself
            # (what side effects would this have this be bad?)
            all = {}
            for r in factors: all[r[0]] = 0
            for r in factors: all[r[0]] += r[1]
            factors = all.items()
            factors.sort() # for doctest consistency
        factors = [r**e for r,e in factors]
        parts = []
        for d in factors:
            n = numer * prod([r for r in factors if r != d]).inverse_mod(d) % d # we know the inverse exists as the two are relatively prime
            parts.append(n/d)
        return whole, parts

    def __call__(self, *x):
        """
        Evaluate the fraction at the given arguments.  This assumes
        that a call function is defined for the numerator and denominator.

        EXAMPLES:
            sage: x = PolynomialRing(RationalField(),'x',3).gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            -2*x1*x2 + x0 + x1
            sage: f(1,2,5)
            -17
            sage: h = f /(x[1] + x[2])
            sage: h
            (-2*x1*x2 + x0 + x1)/(x1 + x2)
            sage: h(1,2,5)
            -17/7
        """
        return self.__numerator(*x) / self.__denominator(*x)

    def _is_atomic(self):
        """
        EXAMPLES:
            sage: K.<x> = Frac(ZZ['x'])
            sage: x._is_atomic()
            True
            sage: f = 1/(x+1)
            sage: f._is_atomic()
            False
        """
        return self.__numerator._is_atomic() and self.__denominator._is_atomic()

    def _repr_(self):
        """
        EXAMPLES:
            sage: K.<x> = Frac(ZZ['x'])
            sage: repr(x+1)
            'x + 1'
            sage: repr((x+1)/(x-1))
            '(x + 1)/(x - 1)'
            sage: repr(1/(x-1))
            '1/(x - 1)'
            sage: repr(1/x)
            '1/x'
        """
        if self.is_zero():
            return "0"
        s = "%s"%self.__numerator
        if self.__denominator != 1:
            denom_string = str( self.__denominator )
            if self.__denominator._is_atomic() and not ('*' in denom_string or '/' in denom_string):
                s = "%s/%s"%(self.__numerator._coeff_repr(no_space=False),denom_string)
            else:
                s = "%s/(%s)"%(self.__numerator._coeff_repr(no_space=False),denom_string)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of this rational function.

        EXAMPLES:
            sage: R = PolynomialRing(QQ, 'x')
            sage: F = R.fraction_field()
            sage: x = F.gen()
            sage: a = x^2 / 1
            sage: latex(a)
            x^{2}
            sage: latex(x^2/(x^2+1))
            \frac{x^{2}}{x^{2} + 1}
            sage: a = 1/x
            sage: latex(a)
            \frac{1}{x}

        TESTS:
            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: z = FractionFieldElement(F, 0, R.gen(), coerce=False)
            sage: z.numerator() == 0
            True
            sage: z.denominator() == R.gen()
            True
            sage: latex(z)
            0
        """
        if self.is_zero():
            return "0"
        if self.__denominator == 1:
            return latex.latex(self.__numerator)
        return "\\frac{%s}{%s}"%(latex.latex(self.__numerator),
                                 latex.latex(self.__denominator))

    def _magma_init_(self, magma):
        """
        Return a string representation of self Magma can understand.

        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: magma((x^2 + x + 1)/(x + 1))          # optional - magma
            (x^2 + x + 1)/(x + 1)

            sage: R.<x,y> = QQ[]
            sage: magma((x+y)/x)                        # optional - magma
            (x + y)/x
        """
        pgens = magma(self.parent()).gens()

        s = self._repr_()
        for i, j in zip(self.parent().variable_names(), pgens):
            s = s.replace(i, j.name())

        return s

    def _add_(self, right):
        """
        EXAMPLES:
            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: x+y
            x + y
            sage: 1/x + 1/y
            (x + y)/(x*y)
            sage: 1/x + 1/(x*y)
            (y + 1)/(x*y)
            sage: Frac(CDF['x']).gen() + 3
            1.0*x + 3.0
        """
        if self.parent().is_exact():
            try:
                gcd_denom = self.__denominator.gcd(right.__denominator)
                if not gcd_denom.is_unit():
                    right_mul = self.__denominator // gcd_denom
                    self_mul = right.__denominator // gcd_denom
                    numer = self.__numerator * self_mul + right.__numerator * right_mul
                    denom = self.__denominator * self_mul
                    new_gcd = numer.gcd(denom)
                    if not new_gcd.is_unit():
                        numer = numer // new_gcd
                        denom = denom // new_gcd
                    return FractionFieldElement(self.parent(), numer, denom, coerce=False, reduce=False)
                # else: no reduction necessary
            except AttributeError: # missing gcd or quo_rem, don't reduce
                pass
            except NotImplementedError: # unimplemented gcd or quo_rem, don't reduce
                pass
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator + self.__denominator*right.__numerator,
           self.__denominator*right.__denominator,  coerce=False, reduce=False)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: K.<t> = Frac(GF(7)['t'])
            sage: t - 1/t
            (t^2 + 6)/t
        """
        if self.parent().is_exact():
            try:
                gcd_denom = self.__denominator.gcd(right.__denominator)
                if not gcd_denom.is_unit():
                    right_mul = self.__denominator // gcd_denom
                    self_mul = right.__denominator // gcd_denom
                    numer = self.__numerator * self_mul -  right.__numerator * right_mul
                    denom = self.__denominator * self_mul
                    new_gcd = numer.gcd(denom)
                    if not new_gcd.is_unit():
                        numer = numer // new_gcd
                        denom = denom // new_gcd
                    return FractionFieldElement(self.parent(), numer, denom, coerce=False, reduce=False)
                # else: no reduction necessary
            except AttributeError: # missing gcd or quo_rem, don't reduce
                pass
            except NotImplementedError: # unimplemented gcd or quo_rem, don't reduce
                pass
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator - self.__denominator*right.__numerator,
           self.__denominator*right.__denominator,  coerce=False, reduce=False)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: K.<t> = Frac(GF(7)['t'])
            sage: a = t/(1+t)
            sage: b = 3/t
            sage: a*b
            3/(t + 1)
        """
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__numerator,
           self.__denominator*right.__denominator, coerce=False, reduce=True)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: K.<x,y,z> = Frac(ZZ['x,y,z'])
            sage: a = (x+1)*(x+y)/(z-3)
            sage: b = (x+y)/(z-1)
            sage: a/b
            (x*z - x + z - 1)/(z - 3)
        """
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator,
           self.__denominator*right.__numerator, coerce=False, reduce=True)

    def derivative(self, *args):
        r"""
        The derivative of this rational function, with respect to variables
        supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        SEE ALSO:
            self._derivative()

        EXAMPLES:
            sage: F = FractionField(PolynomialRing(RationalField(),'x'))
            sage: x = F.gen()
            sage: (1/x).derivative()
            -1/x^2

            sage: (x+1/x).derivative(x, 2)
            2/x^3

            sage: F = FractionField(PolynomialRing(RationalField(),'x,y'))
            sage: x,y = F.gens()
            sage: (1/(x+y)).derivative(x,y)
            2/(x^3 + 3*x^2*y + 3*x*y^2 + y^3)

        """
        return multi_derivative(self, args)

    def _derivative(self, var=None):
        r"""
        Return the derivative of this rational function with respect to the
        variable var.

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: F = FractionField(PolynomialRing(RationalField(),'x'))
            sage: x = F.gen()
            sage: t = 1/x^2
            sage: t._derivative(x)
            -2/x^3
            sage: t.derivative()
            -2/x^3

            sage: F = FractionField(PolynomialRing(RationalField(),'x,y'))
            sage: x,y = F.gens()
            sage: t = (x*y/(x+y))
            sage: t._derivative(x)
            y^2/(x^2 + 2*x*y + y^2)
            sage: t._derivative(y)
            x^2/(x^2 + 2*x*y + y^2)
        """
        if var is None:
            bvar = None
        elif var in self.parent().gens():
            bvar = self.parent().ring()(var)
        else:
            bvar = var

        return (self.__numerator._derivative(bvar)*self.__denominator \
                    - self.__numerator*self.__denominator._derivative(bvar))/\
                    self.__denominator**2

    def __int__(self):
        """
        EXAMPLES:
            sage: K = Frac(ZZ['x'])
            sage: int(K(-3))
            -3
            sage: K.<x> = Frac(RR['x'])
            sage: x/x
            1.00000000000000*x/(1.00000000000000*x)
            sage: int(x/x)
            1
            sage: int(K(.5))
            0
        """
        if self.__denominator != 1:
            self.reduce()
        if self.__denominator == 1:
            return int(self.__numerator)
        else:
            raise TypeError, "denominator must equal 1"

    def _integer_(self, Z=ZZ):
        """
        EXAMPLES:
            sage: K = Frac(ZZ['x'])
            sage: K(5)._integer_()
            5
            sage: K.<x> = Frac(RR['x'])
            sage: ZZ(2*x/x)
            2
        """
        if self.__denominator != 1:
            self.reduce()
        if self.__denominator == 1:
            return Z(self.__numerator)
        raise TypeError, "no way to coerce to an integer."

    def _rational_(self, Q=QQ):
        """
        EXAMPLES:
            sage: K.<x> = Frac(QQ['x'])
            sage: K(1/2)._rational_()
            1/2
            sage: K(1/2 + x/x)._rational_()
            3/2
        """
        return Q(self.__numerator) / Q(self.__denominator)

    def __long__(self):
        """
        EXAMPLES:
            sage: K.<x> = Frac(QQ['x'])
            sage: long(K(3))
            3L
            sage: long(K(3/5))
            0L
        """
        return long(int(self))

    def __pow__(self, right):
        r"""
        Returns self raised to the $right^{th}$ power.

        Note that we need to check whether or not right is negative so we
        don't set __numerator or __denominator to an element of the
        fraction field instead of the underlying ring.

        EXAMPLES:
            sage: R = QQ['x','y']
            sage: FR = R.fraction_field()
            sage: x,y = FR.gens()
            sage: a = x^2; a
            x^2
            sage: type(a.numerator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: type(a.denominator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: a = x^(-2); a
            1/x^2
            sage: type(a.numerator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: type(a.denominator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: x^0
            1
            sage: ((x+y)/(x-y))^2
            (x^2 + 2*x*y + y^2)/(x^2 - 2*x*y + y^2)
            sage: ((x+y)/(x-y))^-2
            (x^2 - 2*x*y + y^2)/(x^2 + 2*x*y + y^2)
            sage: ((x+y)/(x-y))^0
            1
        """
        if right == 0:
            return FractionFieldElement(self.parent(), 1, 1, reduce=False)
        elif right > 0:
            return FractionFieldElement(self.parent(),
                                        self.__numerator**right,
                                        self.__denominator**right, coerce=False, reduce=False)
        else:
            right = -right
            return FractionFieldElement(self.parent(),
                                        self.__denominator**right,
                                        self.__numerator**right, coerce=False, reduce=False)

    def __neg__(self):
        """
        EXAMPLES:
            sage: K.<t> = Frac(GF(5)['t'])
            sage: f = (t^2+t)/(t+2); f
            (t^2 + t)/(t + 2)
            sage: -f
            (4*t^2 + 4*t)/(t + 2)
        """
        return FractionFieldElement(self.parent(), -self.__numerator, self.__denominator,
                                    coerce=False, reduce=False)

    def __abs__(self):
        """
        EXAMPLES:
            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: abs(FractionFieldElement(QQ, -2, 3, coerce=False))
            2/3
        """
        return abs(self.__numerator)/abs(self.__denominator)

    def __invert__(self):
        """
        EXAMPLES:
            sage: K.<t> = Frac(GF(7)['t'])
            sage: f = (t^2+5)/(t-1)
            sage: ~f
            (t + 6)/(t^2 + 5)
        """
        if self.is_zero():
            raise ZeroDivisionError, "Cannot invert 0"
        return FractionFieldElement(self.parent(),
           self.__denominator, self.__numerator, coerce=False, reduce=False)

    def __float__(self):
        """
        EXAMPLES:
            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: float(x/x + y/y)
            2.0
        """
        return float(self.__numerator) / float(self.__denominator)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: K.<t> = Frac(GF(7)['t'])
            sage: t/t == 1
            True
            sage: t+1/t == (t^2+1)/t
            True
            sage: t == t/5
            False
        """
        return cmp(self.__numerator * other.__denominator, self.__denominator*other.__numerator)

    def valuation(self):
        """
        Return the valuation of self, assuming that the numerator and
        denominator have valuation functions defined on them.

        EXAMPLES:
            sage: x = PolynomialRing(RationalField(),'x').gen()
            sage: f = (x**3 + x)/(x**2 - 2*x**3)
            sage: f
            (x^2 + 1)/(-2*x^2 + x)
            sage: f.valuation()
            -1
        """
        return self.__numerator.valuation() - self.__denominator.valuation()

    def is_zero(self):
        """
        Returns True if this element is equal to zero.

        EXAMPLES:
            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: t = F(0)/x
            sage: t.is_zero()
            True
            sage: u = 1/x - 1/x
            sage: u.is_zero()
            True
            sage: u.parent() is F
            True

        """
        return self.__numerator.is_zero()
