"""
Fraction Field Elements

AUTHORS:

- William Stein (input from David Joyner, David Kohel, and Joe
  Wetherell)
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

from sage.structure.element cimport FieldElement, ModuleElement, RingElement, \
        Element

import integer_ring
from integer_ring import ZZ
from rational_field import QQ

import sage.misc.latex as latex
from sage.misc.misc import prod
from sage.misc.derivative import multi_derivative

def is_FractionFieldElement(x):
    """
    Returns whether or not x is of type FractionFieldElement

    EXAMPLES::

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

cdef class FractionFieldElement(FieldElement):
    """
    EXAMPLES::

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
    cdef object __numerator
    cdef object __denominator

    def __init__(self, parent, numerator, denominator=1,
                 coerce=True, reduce=True):
        """
        EXAMPLES::

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
        FieldElement.__init__(self, parent)
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

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: K = GF(7)['a,b'].fraction_field()
            sage: a,b = K.gens()

        ::

            sage: phi = F.hom([a+b, a*b], K)
            sage: phi(x+y) # indirect doctest
            a*b + a + b

        ::

            sage: (x^2/y)._im_gens_(K, [a+b, a*b])
            (a^2 + 2*a*b + b^2)/(a*b)
            sage: (x^2/y)._im_gens_(K, [a, a*b])
            a/b
        """
        nnum = codomain.coerce(self.__numerator._im_gens_(codomain, im_gens))
        nden = codomain.coerce(self.__denominator._im_gens_(codomain, im_gens))
        return codomain.coerce(nnum/nden)

    def reduce(self):
        """
        Divides out the gcd of the numerator and denominator.

        Automatically called for exact rings, but because it may be
        numerically unstable for inexact rings it must be called manually
        in that case.

        EXAMPLES::

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
                    denom = denom.parent().one_element()
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
        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.copy()
            (x + y)/y
        """
        return self.__class__(self._parent, self.__numerator,
                self.__denominator, coerce=False, reduce=False)

    def numerator(self):
        """
        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.numerator()
            x + y
        """
        return self.__numerator

    def denominator(self):
        """
        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.denominator()
            y
        """
        return self.__denominator

    def __hash__(self):
        """
        This function hashes in a special way to ensure that generators of
        a ring R and generators of a fraction field of R have the same
        hash. This enables them to be used as keys interchangably in a
        dictionary (since ``==`` will claim them equal). This
        is particularly useful for methods like subs on
        ``ParentWithGens`` if you are passing a dictionary of
        substitutions.

        EXAMPLES::

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
        Decomposes fraction field element into a whole part and a list of
        fraction field elements over prime power denominators.

        The sum will be equal to the original fraction.

        AUTHORS:

        - Robert Bradshaw (2007-05-31)

        EXAMPLES::

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

        ::

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
        if not self._parent.is_exact():
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
        Evaluate the fraction at the given arguments. This assumes that a
        call function is defined for the numerator and denominator.

        EXAMPLES::

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
        EXAMPLES::

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
        EXAMPLES::

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

        EXAMPLES::

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

        TESTS::

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

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: magma((x^2 + x + 1)/(x + 1))          # optional - magma
            (x^2 + x + 1)/(x + 1)

        ::

            sage: R.<x,y> = QQ[]
            sage: magma((x+y)/x)                        # optional - magma
            (x + y)/x
        """
        pgens = magma(self._parent).gens()

        s = self._repr_()
        for i, j in zip(self._parent.variable_names(), pgens):
            s = s.replace(i, j.name())

        return s

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::

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
        r_numerator = (<FractionFieldElement>right).__numerator
        r_denominator = (<FractionFieldElement>right).__denominator
        if self._parent.is_exact():
            try:
                gcd_denom = self.__denominator.gcd(r_denominator)
                if not gcd_denom.is_unit():
                    right_mul = self.__denominator // gcd_denom
                    self_mul = r_denominator // gcd_denom
                    numer = self.__numerator * self_mul + \
                            r_numerator * right_mul
                    denom = self.__denominator * self_mul
                    new_gcd = numer.gcd(denom)
                    if not new_gcd.is_unit():
                        numer = numer // new_gcd
                        denom = denom // new_gcd
                    return self.__class__(self._parent, numer, denom,
                            coerce=False, reduce=False)
                # else: no reduction necessary
            except AttributeError: # missing gcd or quo_rem, don't reduce
                pass
            except NotImplementedError: # unimplemented gcd or quo_rem, don't reduce
                pass
        return self.__class__(self._parent,
           self.__numerator*r_denominator + self.__denominator*r_numerator,
           self.__denominator*r_denominator,  coerce=False, reduce=False)

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: t - 1/t
            (t^2 + 6)/t
        """
        r_numerator = (<FractionFieldElement>right).__numerator
        r_denominator = (<FractionFieldElement>right).__denominator
        if self._parent.is_exact():
            try:
                gcd_denom = self.__denominator.gcd(r_denominator)
                if not gcd_denom.is_unit():
                    right_mul = self.__denominator // gcd_denom
                    self_mul = r_denominator // gcd_denom
                    numer = self.__numerator * self_mul - \
                            r_numerator * right_mul
                    denom = self.__denominator * self_mul
                    new_gcd = numer.gcd(denom)
                    if not new_gcd.is_unit():
                        numer = numer // new_gcd
                        denom = denom // new_gcd
                    return self.__class__(self._parent, numer, denom,
                            coerce=False, reduce=False)
                # else: no reduction necessary
            except AttributeError: # missing gcd or quo_rem, don't reduce
                pass
            except NotImplementedError: # unimplemented gcd or quo_rem, don't reduce
                pass
        return self.__class__(self._parent,
           self.__numerator*r_denominator - self.__denominator*r_numerator,
           self.__denominator*r_denominator,  coerce=False, reduce=False)

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: a = t/(1+t)
            sage: b = 3/t
            sage: a*b
            3/(t + 1)
        """
        return self.__class__(self._parent,
           self.__numerator*(<FractionFieldElement>right).__numerator,
           self.__denominator*(<FractionFieldElement>right).__denominator,
           coerce=False, reduce=True)

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<x,y,z> = Frac(ZZ['x,y,z'])
            sage: a = (x+1)*(x+y)/(z-3)
            sage: b = (x+y)/(z-1)
            sage: a/b
            (x*z - x + z - 1)/(z - 3)
        """
        return self.__class__(self._parent,
           self.__numerator*(<FractionFieldElement>right).__denominator,
           self.__denominator*(<FractionFieldElement>right).__numerator,
           coerce=False, reduce=True)

    def derivative(self, *args):
        r"""
        The derivative of this rational function, with respect to variables
        supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. seealso::

           :meth:`_derivative`

        EXAMPLES::

            sage: F = FractionField(PolynomialRing(RationalField(),'x'))
            sage: x = F.gen()
            sage: (1/x).derivative()
            -1/x^2

        ::

            sage: (x+1/x).derivative(x, 2)
            2/x^3

        ::

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

        .. seealso::

           :meth:`derivative`

        EXAMPLES::

            sage: F = FractionField(PolynomialRing(RationalField(),'x'))
            sage: x = F.gen()
            sage: t = 1/x^2
            sage: t._derivative(x)
            -2/x^3
            sage: t.derivative()
            -2/x^3

        ::

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
        elif var in self._parent.gens():
            bvar = self._parent.ring()(var)
        else:
            bvar = var

        return (self.__numerator._derivative(bvar)*self.__denominator \
                    - self.__numerator*self.__denominator._derivative(bvar))/\
                    self.__denominator**2

    def __int__(self):
        """
        EXAMPLES::

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
        EXAMPLES::

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
        EXAMPLES::

            sage: K.<x> = Frac(QQ['x'])
            sage: K(1/2)._rational_()
            1/2
            sage: K(1/2 + x/x)._rational_()
            3/2
        """
        return Q(self.__numerator) / Q(self.__denominator)

    def __long__(self):
        """
        EXAMPLES::

            sage: K.<x> = Frac(QQ['x'])
            sage: long(K(3))
            3L
            sage: long(K(3/5))
            0L
        """
        return long(int(self))

    def __pow__(self, right, dummy):
        r"""
        Returns self raised to the `right^{th}` power.

        Note that we need to check whether or not right is negative so we
        don't set __numerator or __denominator to an element of the
        fraction field instead of the underlying ring.

        EXAMPLES::

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
        s_num = (<FractionFieldElement>self).__numerator
        s_den = (<FractionFieldElement>self).__denominator
        if right == 0:
            return self.__class__(self.parent(), 1, 1, reduce=False)
        elif right > 0:
            return self.__class__(self.parent(),
                                        s_num**right, s_den**right,
                                        coerce=False, reduce=False)
        else:
            right = -right
            return self.__class__(self.parent(),
                                        s_den**right, s_num**right,
                                        coerce=False, reduce=False)

    def __neg__(self):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: f = (t^2+t)/(t+2); f
            (t^2 + t)/(t + 2)
            sage: -f
            (4*t^2 + 4*t)/(t + 2)
        """
        return self.__class__(self._parent,
                -self.__numerator, self.__denominator,
                coerce=False, reduce=False)

    def __abs__(self):
        """
        EXAMPLES::

            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: abs(FractionFieldElement(QQ, -2, 3, coerce=False))
            2/3
        """
        return abs(self.__numerator)/abs(self.__denominator)

    def __invert__(self):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: f = (t^2+5)/(t-1)
            sage: ~f
            (t + 6)/(t^2 + 5)
        """
        if self.is_zero():
            raise ZeroDivisionError, "Cannot invert 0"
        return self.__class__(self._parent,
           self.__denominator, self.__numerator, coerce=False, reduce=False)

    def __float__(self):
        """
        EXAMPLES::

            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: float(x/x + y/y)
            2.0
        """
        return float(self.__numerator) / float(self.__denominator)

    def __richcmp__(left, right, int op):
        """
        EXAMPLES::

            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: x > y
            True
            sage: 1 > y
            False
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(self, Element other) except -2:
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: t/t == 1
            True
            sage: t+1/t == (t^2+1)/t
            True
            sage: t == t/5
            False
        """
        return cmp(self.__numerator * \
                (<FractionFieldElement>other).__denominator,
                self.__denominator*(<FractionFieldElement>other).__numerator)

    def valuation(self):
        """
        Return the valuation of self, assuming that the numerator and
        denominator have valuation functions defined on them.

        EXAMPLES::

            sage: x = PolynomialRing(RationalField(),'x').gen()
            sage: f = (x**3 + x)/(x**2 - 2*x**3)
            sage: f
            (x^2 + 1)/(-2*x^2 + x)
            sage: f.valuation()
            -1
        """
        return self.__numerator.valuation() - self.__denominator.valuation()

    def __nonzero__(self):
        """
        Returns True if this element is nonzero.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: t = F(0)/x
            sage: t.__nonzero__()
            False

        ::

            sage: (1/x).__nonzero__()
            True
        """
        return not self.__numerator.is_zero()

    def is_zero(self):
        """
        Returns True if this element is equal to zero.

        EXAMPLES::

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

    def is_one(self):
        """
        Returns True if this element is equal to one.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: (x/x).is_one()
            True
            sage: (x/y).is_one()
            False
        """
        return self.__numerator == self.__denominator

    def __reduce__(self):
        """
        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: f = F.random_element()
            sage: loads(f.dumps()) == f
            True
        """
        return (make_element,
                (self._parent, self.__numerator, self.__denominator))

def make_element(parent, numerator, denominator):
    """
    Used for unpickling FractionFieldElement objects (and subclasses).

    EXAMPLES::

        sage: from sage.rings.fraction_field_element import make_element
        sage: R = ZZ['x,y']
        sage: x,y = R.gens()
        sage: F = R.fraction_field()
        sage: make_element(F, 1+x, 1+y)
        (x + 1)/(y + 1)
    """

    return parent._element_class(parent, numerator, denominator)

def make_element_old(parent, cdict):
    """
    Used for unpickling old FractionFieldElement pickles.

    EXAMPLES::

        sage: from sage.rings.fraction_field_element import make_element_old
        sage: R.<x,y> = ZZ[]
        sage: F = R.fraction_field()
        sage: make_element_old(F, {'_FractionFieldElement__numerator':x+y,'_FractionFieldElement__denominator':x-y})
        (x + y)/(x - y)
    """
    return FractionFieldElement(parent,
            cdict['_FractionFieldElement__numerator'],
            cdict['_FractionFieldElement__denominator'],
            coerce=False, reduce=False)

