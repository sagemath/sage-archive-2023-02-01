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
import fraction_field
import integer_ring

import sage.misc.latex as latex

def is_FractionFieldElement(x):
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
        return FractionFieldElement(self.parent(), self.__numerator,
                                    self.__denominator, coerce=False, reduce=False)

    def numerator(self):
        return self.__numerator

    def denominator(self):
        return self.__denominator

    def __call__(self, *x):
        """
        Evaluate the fraction at the given arguments.  This assumes
        that a call function is defined for the numerator and denominator.

        EXAMPLES:
            sage: x = MPolynomialRing(RationalField(),'x',3).gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            x1 - 2*x1*x2 + x0
            sage: f(1,2,5)
            -17
            sage: h = f /(x[1] + x[2])
            sage: h
            (x1 - 2*x1*x2 + x0)/(x2 + x1)
            sage: h(1,2,5)
            -17/7
        """
        return self.__numerator(x) / self.__denominator(x)

    def _is_atomic(self):
        return self.__numerator._is_atomic() and self.__denominator._is_atomic()

    def _repr_(self):
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
        """
        Return a latex representation of this rational function.
        """
        if self.__denominator == 1:
            return latex.latex(self.__numerator)
        return "\\frac{%s}{%s}"%(latex.latex(self.__numerator),
                                 latex.latex(self.__denominator))

    def _add_(self, right):
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator + self.__denominator*right.__numerator,
           self.__denominator*right.__denominator, coerce=False)

    def _sub_(self, right):
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator - self.__denominator*right.__numerator,
           self.__denominator*right.__denominator,  coerce=False)

    def _mul_(self, right):
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__numerator,
           self.__denominator*right.__denominator, coerce=False, reduce=True)

    def _div_(self, right):
        return FractionFieldElement(self.parent(),
           self.__numerator*right.__denominator,
           self.__denominator*right.__numerator, coerce=False, reduce=True)

    def __int__(self):
        if self.__denominator == 1:
            return int(self.__numerator)
        else:
            raise TypeError, "denominator must equal 1"

    def _integer_(self):
        if self.__denominator == 1:
            try:
                return self.__numerator._integer_()
            except AttributeError:
                pass
        raise TypeError, "no way to coerce to an integer."

    def _rational_(self):
        Z = integer_ring.IntegerRing()
        try:
            return Z(self.__numerator) / Z(self.__denominator)
        except AttributeError:
            pass
        raise TypeError, "coercion to rational not defined"

    def __long__(self):
        if self.__denominator == 1:
            return long(self.__numerator)
        else:
            raise TypeError, "denominator must equal 1"

    def __pow__(self, right):
        return FractionFieldElement(self.parent(),
                                    self.__numerator**right,
                                    self.__denominator**right, coerce=False, reduce=False)

    def __neg__(self):
        return FractionFieldElement(self.parent(), -self.__numerator, self.__denominator,
                                    coerce=False, reduce=False)

    def __pos__(self):
        return self

    def __abs__(self):
        return abs(self.__numerator)/abs(self.__denominator)

    def __invert__(self):
        if self.is_zero():
            raise ZeroDivisionError, "Cannot invert 0"
        return FractionFieldElement(self.parent(),
           self.__denominator, self.__numerator)

    def __float__(self):
        return float(self.__numerator) / float(self.__denominator)

    def __cmp__(self, other):
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

