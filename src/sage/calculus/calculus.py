r"""
Symbolic Computation.

AUTHORS:
    Bobby Moretti and William Stein: 2006--2007

The \sage calculus module is loosely based on the \sage Enhancement Proposal
found at: \url{http://www.sagemath.org:9001/CalculusSEP.}

EXAMPLES:
    The basic units of the calculus package are symbolic expressions
    which are elements of the symbolic expression ring (SR). There are
    many subclasses of \class{SymbolicExpression}. The most basic of these is
    the formal indeterminate class, \class{SymbolicVariable}. To create a
    \class{SymbolicVariable} object in \sage, use the \code{var()} method, whose
    argument is the text of that variable.  Note that \sage is
    intelligent about {\LaTeX}ing variable names.

        sage: x1 = var('x1'); x1
        x1
        sage: latex(x1)
        x_{1}
        sage: theta = var('theta'); theta
        theta
        sage: latex(theta)
        \theta

    \sage predefines \code{x} to be a global indeterminate. Thus the following works:
        sage: x^2
        x^2
        sage: type(x)
        <class 'sage.calculus.calculus.SymbolicVariable'>

    More complicated expressions in \sage can be built up using
    ordinary arithmetic. The following are valid, and follow the rules
    of Python arithmetic: (The '=' operator represents assignment, and
    not equality)
        sage: var('x,y,z')
        (x, y, z)
        sage: f = x + y + z/(2*sin(y*z/55))
        sage: g = f^f; g
        (z/(2*sin(y*z/55)) + y + x)^(z/(2*sin(y*z/55)) + y + x)

    Differentiation and integration are available, but behind the
    scenes through maxima:

        sage: f = sin(x)/cos(2*y)
        sage: f.derivative(y)
        2*sin(x)*sin(2*y)/cos(2*y)^2
        sage: g = f.integral(x); g
        -cos(x)/cos(2*y)

    Note that these methods require an explicit variable name. If none
    is given, \sage will try to find one for you.
        sage: f = sin(x); f.derivative()
        cos(x)

    However when this is ambiguous, \sage will raise an exception:
        sage: f = sin(x+y); f.derivative()
        Traceback (most recent call last):
        ...
        ValueError: must supply an explicit variable for an expression containing more than one variable

    Substitution works similarly. We can substitute with a python dict:
        sage: f = sin(x*y - z)
        sage: f({x: var('t'), y: z})
        sin(t*z - z)

    Also we can substitute with keywords:
        sage: f = sin(x*y - z)
        sage: f(x = t, y = z)
        sin(t*z - z)

    If there is no ambiguity of variable names, we don't have to specify them:
        sage: f = sin(x)
        sage: f(y)
        sin(y)
        sage: f(pi)
        0

    However if there is ambiguity, we must explicitly state what
    variables we're substituting for:

        sage: f = sin(2*pi*x/y)
        sage: f(4)
        sin(8*pi/y)

    We can also make a \class{CallableSymbolicExpression}, which is a
    \class{SymbolicExpression} that is a function of specified
    variables in a fixed order. Each \class{SymbolicExpression} has a
    \code{function(...)} method that is used to create a
    \class{CallableSymbolicExpression}, as illustrated below:

        sage: u = log((2-x)/(y+5))
        sage: f = u.function(x, y); f
        (x, y) |--> log((2 - x)/(y + 5))

    There is an easier way of creating a \class{CallableSymbolicExpression}, which
    relies on the \sage preparser.

        sage: f(x,y) = log(x)*cos(y); f
        (x, y) |--> log(x)*cos(y)

    Then we have fixed an order of variables and there is no ambiguity
    substituting or evaluating:

        sage: f(x,y) = log((2-x)/(y+5))
        sage: f(7,t)
        log(-5/(t + 5))

    Some further examples:

        sage: f = 5*sin(x)
        sage: f
        5*sin(x)
        sage: f(2)
        5*sin(2)
        sage: f(pi)
        0
        sage: float(f(pi))             # random low order bits
        6.1232339957367663e-16

    Another example:
        sage: f = integrate(1/sqrt(9+x^2), x); f
        arcsinh(x/3)
        sage: f(3)
        arcsinh(1)
        sage: f.derivative(x)
        1/(3*sqrt(x^2/9 + 1))

    We compute the length of the parabola from 0 to 2:
        sage: x = var('x')
        sage: y = x^2
        sage: dy = derivative(y,x)
        sage: z = integral(sqrt(1 + dy^2), x, 0, 2)
        sage: print z
                             arcsinh(4) + 4 sqrt(17)
                             ---------------------
                                       4
        sage: n(z,200)
        4.6467837624329358733826155674904591885104869874232887508703
        sage: float(z)
        4.6467837624329356


    We test pickling:
        sage: x, y = var('x,y')
        sage: f = -sqrt(pi)*(x^3 + sin(x/cos(y)))
        sage: bool(loads(dumps(f)) == f)
        True

Coercion examples:

We coerce various symbolic expressions into the complex numbers:

    sage: CC(I)
    1.00000000000000*I
    sage: CC(2*I)
    2.00000000000000*I
    sage: ComplexField(200)(2*I)
    2.0000000000000000000000000000000000000000000000000000000000*I
    sage: ComplexField(200)(sin(I))
    1.1752011936438014568823818505956008151557179813340958702296*I
    sage: f = sin(I) + cos(I/2); f
    I*sinh(1) + cosh(1/2)
    sage: CC(f)
    1.12762596520638 + 1.17520119364380*I
    sage: ComplexField(200)(f)
    1.1276259652063807852262251614026720125478471180986674836290 + 1.1752011936438014568823818505956008151557179813340958702296*I
    sage: ComplexField(100)(f)
    1.1276259652063807852262251614 + 1.1752011936438014568823818506*I

We illustrate construction of an inverse sum where each denominator
has a new variable name:
    sage: f = sum(1/var('n%s'%i)^i for i in range(10))
    sage: print f
                 1     1     1     1     1     1     1     1    1
                --- + --- + --- + --- + --- + --- + --- + --- + -- + 1
                  9     8     7     6     5     4     3     2   n1
                n9    n8    n7    n6    n5    n4    n3    n2

Note that after calling var, the variables are immediately available for use:
    sage: (n1 + n2)^5
    (n2 + n1)^5

We can, of course, substitute:
    sage: print f(n9=9,n7=n6)
            1     1     1     1     1     1     1    1    387420490
           --- + --- + --- + --- + --- + --- + --- + -- + ---------
             8     6     7     5     4     3     2   n1   387420489
           n8    n6    n6    n5    n4    n3    n2

TESTS:

Substitution:
    sage: f = x
    sage: f(x=5)
    5

Simplifying expressions involving scientific notation:
    sage: k = var('k')
    sage: a0 = 2e-6; a1 = 12
    sage: c = a1 + a0*k; c
    2.000000000000000e-6*k + 12
    sage: sqrt(c)
    sqrt(2.000000000000000e-6*k + 12)
    sage: sqrt(c^3)
    sqrt((2.000000000000000e-6*k + 12)^3)

The symbolic Calculus package uses its own copy of maxima for
simplification, etc., which is separate from the default system-wide
version:
    sage: maxima.eval('[x,y]: [1,2]')
    '[1,2]'
    sage: maxima.eval('expand((x+y)^3)')
    '27'

If the copy of maxima used by the symbolic calculus package were
the same as the default one, then the following would return 27,
which would be very confusing indeed!
    sage: x, y = var('x,y')
    sage: expand((x+y)^3)
    y^3 + 3*x*y^2 + 3*x^2*y + x^3

Set x to be 5 in maxima:
    sage: maxima('x: 5')
    5
    sage: maxima('x + x + %pi')
    %pi+10

This simplification is done using maxima (behind the scenes):
    sage: x + x + pi
    2*x + pi

Note that \code{x} is still \var{x}, since the maxima used by the calculus package
is different than the one in the interactive interpreter.

Check to see that the problem with the variables method mentioned in Trac
ticket \#3779 is actually fixed:
    sage: f = function('F',x)
    sage: diff(f*SR(1),x)
    diff(F(x), x, 1)

"""

import weakref
import re

from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_MPolynomial, is_MPolynomialRing, is_FractionFieldElement,
                            is_RealNumber, is_ComplexNumber, RR,
                            Integer, Rational, CC, QQ, CDF,
                            QuadDoubleElement,
                            PolynomialRing, ComplexField,
                            algdep, Integer, RealNumber, RealIntervalField)

from sage.rings.real_mpfr import create_RealNumber

from sage.structure.element import RingElement, is_Element
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex, latex_variable_name
from sage.misc.misc import uniq as unique
from sage.structure.sage_object import SageObject

from sage.interfaces.maxima import MaximaElement, Maxima

import sage.numerical.optimize

# The calculus package uses its own copy of maxima, which is
# separate from the default system-wide version.
maxima = Maxima(init_code = ['display2d:false; domain: complex; keepfloat: true'])

from sage.misc.parser import Parser

from sage.calculus.equations import SymbolicEquation
from sage.rings.real_mpfr import RealNumber
from sage.rings.complex_number import ComplexNumber
from sage.rings.real_double import RealDoubleElement
from sage.rings.complex_double import ComplexDoubleElement
from sage.rings.real_mpfi import RealIntervalFieldElement
from sage.rings.infinity import InfinityElement

from sage.libs.pari.gen import pari, PariError, gen as PariGen

from sage.rings.complex_double import ComplexDoubleElement

import sage.functions.constants

from sage.misc.derivative import multi_derivative, derivative_parse

import math
import sage.functions.functions

import sage.ext.fast_eval as fast_float

# TODO: What the heck does this is_simplified thing do?
is_simplified = False

infixops = {operator.add: '+',
            operator.sub: '-',
            operator.mul: '*',
            operator.div: '/',
            operator.pow: '^'}

arc_functions =  ['asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh', 'acoth', 'asech', 'acsch', 'acot', 'acsc', 'asec']


def is_SymbolicExpression(x):
    """
    EXAMPLES:
        sage: is_SymbolicExpression(sin(x))
        True
        sage: is_SymbolicExpression(2/3)
        False
        sage: is_SymbolicExpression(sqrt(2))
        True
    """
    return isinstance(x, SymbolicExpression)

def is_SymbolicExpressionRing(x):
    """
    EXAMPLES:
        sage: is_SymbolicExpressionRing(QQ)
        False
        sage: is_SymbolicExpressionRing(SR)
        True
    """
    return isinstance(x, SymbolicExpressionRing_class)

cache = {}
class uniq(object):
    def __new__(cls):
        global cache
        if cache.has_key(cls):
            return cache[cls]
        O = CommutativeRing.__new__(cls)
        cache[cls] = O
        return O

class SymbolicExpressionRing_class(uniq, CommutativeRing):
    """
    The ring of all formal symbolic expressions.

    EXAMPLES:
        sage: SR
        Symbolic Ring
        sage: type(SR)
        <class 'sage.calculus.calculus.SymbolicExpressionRing_class'>

    TESTS:
    Test serializing:
        sage: loads(dumps(SR)) == SR
        True
    """
    def __init__(self, default_precision=53):
        """
        Create a symbolic expression ring.

        EXAMPLES:
            sage: from sage.calculus.calculus import SymbolicExpressionRing_class
            sage: SymbolicExpressionRing_class()
            Symbolic Ring
        """
        ParentWithBase.__init__(self, RR)

    def __cmp__(self, other):
        """
        Compare two symbolic expression rings.  They are equal if and
        only if they have the same type. Otherwise their types are
        compared.

        EXAMPLES:
            sage: cmp(SR, RR) #random
            1
            sage: cmp(RR, SymbolicExpressionRing()) #random
            -1
            sage: cmp(SR, SymbolicExpressionRing()) #random
            0
        """
        return cmp(type(self), type(other))

    def __contains__(self, x):
        r"""
        True if there is an element of the symbolic ring  that is equal
        to x under ==.

        EXAMPLES:
        The symbolic variable x is in the symbolic ring.
            sage: x.parent()
            Symbolic Ring
            sage: x in SR
            True

        2 is also in the symbolic ring since it is equal to something
        in SR, even though 2's parent is not SR.
            sage: 2 in SR
            True
            sage: parent(2)
            Integer Ring
            sage: 1/3 in SR
            True

        The finite field element 1 (in GF(3)) is not equal to
        anything in SR.
            sage: GF(3)(1) in SR
            False
        """
        try:
            x2 = self(x)
            return bool(x2 == x)
        except TypeError:
            return False

    def __call__(self, x):
        """
        Coerce $x$ into the symbolic expression ring SR.

        EXAMPLES:
            sage: a = SR(-3/4); a
            -3/4
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: a.parent()
            Symbolic Ring
            sage: type(SR(I))
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: is_SymbolicExpression(SR(I))
            True

        If $a$ is already in the symblic expression ring, coercing returns
        $a$ itself (not a copy):
            sage: SR(a) is a
            True

        A Python complex number:
            sage: SR(complex(2,-3))
            2.00000000000000 - 3.00000000000000*I
        """
        if is_SymbolicExpression(x):
            return x
        elif hasattr(x, '_symbolic_'):
            return x._symbolic_(self)
        elif isinstance(x, str):
            try:
                return symbolic_expression_from_string(x)
            except SyntaxError, err:
                msg, s, pos = err.args
                raise TypeError, "%s: %s !!! %s" % (msg, s[:pos], s[pos:])
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        """
        Used for implicit coercion.

        EXAMPLES:
            sage: x=var('x'); y0,y1=PolynomialRing(ZZ,2,'y').gens()
            sage: x+y0/y1
            y0/y1 + x
            sage: x.subs(x=y0/y1)
            y0/y1
        """
        if isinstance(x, CallableSymbolicExpression):
            return x._expr
        elif isinstance(x, SymbolicExpression):
            return x
        elif isinstance(x, MaximaElement):
            return symbolic_expression_from_maxima_element(x)
        # if "x" is a SymPy object, convert it to a SAGE object
        elif is_Polynomial(x) or is_MPolynomial(x):
            if x.base_ring() != self:  # would want coercion to go the other way
                return SymbolicPolynomial(x)
            else:
                raise TypeError, "Basering is Symbolic Ring, please coerce in the other direction."
        elif is_FractionFieldElement(x) and (is_Polynomial(x.numerator()) or is_MPolynomial(x.numerator())):
            if x.base_ring() != self:  # would want coercion to go the other way
                return SymbolicPolynomial(x.numerator()) / SymbolicPolynomial(x.denominator())
            else:
                raise TypeError, "Basering is Symbolic Ring, please coerce in the other direction."
        elif isinstance(x, (RealNumber,
                            RealDoubleElement,
                            RealIntervalFieldElement,
                            float,
                            sage.functions.constants.Constant,
                            Integer,
                            int,
                            Rational,
                            PariGen,
                            ComplexNumber,
                            ComplexDoubleElement,
                            QuadDoubleElement,
                            InfinityElement
                            )):
            return SymbolicConstant(x)
        elif isinstance(x, complex):
            return evaled_symbolic_expression_from_maxima_string('%s+%%i*%s'%(x.real,x.imag))

        from sympy.core.basic import Basic
        if isinstance(x, Basic):
            return self(x._sage_())
        else:
            raise TypeError, "cannot coerce type '%s' into a SymbolicExpression."%type(x)

    def _repr_(self):
        """
        Return string representation of this symbolic ring.

        EXAMPLES:
            sage: SR._repr_()
            'Symbolic Ring'
        """
        return 'Symbolic Ring'

    def _latex_(self):
        """
        Return latex representation of the symbolic ring.

        EXAMPLES:
            sage: latex(SR)
            \text{SR}
            sage: M = MatrixSpace(SR, 2); latex(M)
            \mbox{\rm Mat}_{2\times 2}(\text{SR})
        """
        return r'\text{SR}'

    def var(self, x):
        """
        Return the symbolic variable defined by x as an element of the
        symbolic ring.

        EXAMPLES:
            sage: zz = SR.var('zz'); zz
            zz
            sage: type(zz)
            <class 'sage.calculus.calculus.SymbolicVariable'>
            sage: t = SR.var('theta2'); t
            theta2
        """
        return var(x)

    def characteristic(self):
        """
        Return the characteristic of the symbolic ring, which is 0.

        OUTPUT:
            a Sage integer

        EXAMPLES:
            sage: c = SR.characteristic(); c
            0
            sage: type(c)
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(0)

    def _an_element_impl(self):
        """
        Return an element of the symbolic ring, which is used
        by the coercion model.

        EXAMPLES:
        Currently this function always returns 0.  That may change.
            sage: SR._an_element_impl()
            0
        """
        try:
            return self.__zero
        except AttributeError:
            self.__zero = SR(0)
        return self.__zero

    def is_field(self):
        """
        Returns True, since the symbolic expression ring is (for the
        most part) a field.

        EXAMPLES:
            sage: SR.is_field()
            True
        """
        return True

    def is_exact(self):
        """
        Return False, because there are approximate elements in
        the symbolic ring.

        EXAMPLES:
            sage: SR.is_exact()
            False

        Here is an inexact element.
            sage: SR(1.9393)
            1.93930000000000
        """
        return False

# Define the unique symbolic expression ring.
SR = SymbolicExpressionRing_class()

# The factory function that returns the unique SR.
def SymbolicExpressionRing():
    """
    Return the symbolic expression ring.  There is one
    globally defines symbolic expression ring in the
    calculus module.

    EXAMPLES:
        sage: SymbolicExpressionRing()
        Symbolic Ring
        sage: SymbolicExpressionRing() is SR
        True
    """
    return SR

class SymbolicExpression(RingElement):
    r"""
    A Symbolic Expression.

    EXAMPLES:
        Some types of \class{SymbolicExpression}s:

        sage: a = SR(2+2); a
        4
        sage: type(a)
        <class 'sage.calculus.calculus.SymbolicConstant'>

    """
    def __init__(self):
        """
        Create a symbolic expression.

        EXAMPLES:
        This example is mainly for testing purposes.

        We explicitly import the SymbolicExpression class.
            sage: from sage.calculus.calculus import SymbolicExpression

        Then we make an instance of it.  Note that it prints as a
        ``generic element'', since it doesn't even have a specific
        value!
            sage: a = SymbolicExpression(); a
            Generic element of a structure

        It's of the right type.
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicExpression'>

        And it has the right parent.
            sage: a.parent()
            Symbolic Ring
        """
        RingElement.__init__(self, SR)
        if is_simplified:
            self._simp = None

    def __hash__(self):
        """
        Returns the hash of the simplified string representation of
        this symbolic expression.

        EXAMPLES:
        We hash a symbolic polynomial:
            sage: hash(x^2 + 1) #random due to architecture dependence
            -832266011

        The default hashing strategy is to simply hash
        the string representation of the simplified form.
            sage: hash(repr(x^2+1)) #random due to architecture dependence
            -832266011

        In some cases a better hashing strategy is used.
            sage: hash(SR(3/1))
            3
            sage: hash(repr(SR(3/1))) #random due to architecture dependence
            -2061914958

        In this example hashing is important otherwise the answer is
        wrong:
            sage: uniq([x-x, -x+x])
            [0]
        """
        return hash(self._repr_(simplify=True))

    def __nonzero__(self):
        """
        Return True if this element is definitely not zero.

        EXAMPLES:
            sage: k = var('k')
            sage: pol = 1/(k-1) - 1/k -1/k/(k-1);
            sage: pol.is_zero()
            True

            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: f.is_zero()
            True
        """

        try:
            return self.__nonzero
        except AttributeError:
            ans = not bool(self == SR.zero_element())
            self.__nonzero = ans
        return ans

    def __str__(self):
        """
        Printing an object explicitly gives ASCII art:

        EXAMPLES:
            sage: var('x y')
            (x, y)
            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f
            y^2/(y + 1)^3 + x/(x - 1)^3
            sage: print f
                                              2
                                             y          x
                                          -------- + --------
                                                 3          3
                                          (y + 1)    (x - 1)

            sage: f = (exp(x)-1)/(exp(x/2)+1)
            sage: g = exp(x/2)-1
            sage: print f(10), g(10)
                                     10
                                    e   - 1
                                   --------
                                     5
                                    e  + 1
                                      5
                                     e  - 1
        """
        return '\n' + self.display2d(onscreen=False)

    def show(self):
        """
        Show this symbolic expression, i.e., typeset it nicely.

        EXAMPLES:
            sage: (x^2 + 1).show()
            {x}^{2}  + 1
        """
        from sage.misc.functional import _do_show
        return _do_show(self)

    def display2d(self, onscreen=True):
        r"""
        Display \code{self} using ASCII art.

        INPUT:
            onscreen -- string (optional, default True) If True,
                displays; if False, returns string.

        EXAMPLES:
        We display a fraction:
            sage: var('x,y')
            (x, y)
            sage: f = (x^3+y)/(x+3*y^2+1); f
            (y + x^3)/(3*y^2 + x + 1)
            sage: print f
                                                     3
                                                y + x
                                             ------------
                                                2
                                             3 y  + x + 1

        Use \code{onscreen=False} to get the 2d string:
             sage: f.display2d(onscreen=False)
             '                                         3\r\n                                    y + x\r\n                                 ------------\r\n                                    2\r\n                                 3 y  + x + 1'

        ASCII art really helps for the following integral:
            sage: f = integral(sin(x^2)); f
            sqrt(pi)*((sqrt(2)*I + sqrt(2))*erf((sqrt(2)*I + sqrt(2))*x/2) + (sqrt(2)*I - sqrt(2))*erf((sqrt(2)*I - sqrt(2))*x/2))/8
            sage: print f
                                                         (sqrt(2)  I + sqrt(2)) x
                   sqrt( pi) ((sqrt(2)  I + sqrt(2)) erf(------------------------)
                                                                    2
                                                               (sqrt(2)  I - sqrt(2)) x
                                  + (sqrt(2)  I - sqrt(2)) erf(------------------------))/8
                                                                          2
        """
        if not self._has_been_simplified():
            self = self.simplify()
        s = self._maxima_().display2d(onscreen=False)
        s = s.replace('%pi',' pi').replace('%i',' I').replace('%e', ' e')

        #Change asin, etc. to arcsin, etc
        for arc_function in arc_functions:
            s = s.replace("  "+arc_function, "arc"+arc_function[1:])
            s = s.replace(arc_function, "arc"+arc_function[1:])

        if onscreen:
            print s
        else:
            return s

    def _has_been_simplified(self):
        """
        Return True if this symbolic expression was constructed in
        such a way that it is known to already be simplified.

        WARNING: An expression that happens to be in a simplified form
        need not return True.

        EXAMPLES:
        This expression starts in simple form:
            sage: f = x^2 + 1; f
            x^2 + 1

        But it has not been simplified, i.e., it was constructed
        as a simplified expression.
            sage: f._has_been_simplified()
            False
            sage: g = f.simplify(); g
            x^2 + 1
            sage: g._has_been_simplified()
            True

        This function still does not return True, since f itself
        isn't constructed as a simplified expression.
            sage: f._has_been_simplified()
            False

        Here, though f looks simplified when printed, internally
        it is much more complicated and hence not simplified.
            sage: f = x + 1 - x; f
            1
            sage: f._has_been_simplified()
            False
            sage: type(f)
            <class 'sage.calculus.calculus.SymbolicArithmetic'>
            sage: f._operands
            [x + 1, x]
        """
        return hasattr(self, '_simp') and self._simp is None

    def _declare_simplified(self):
        """
        Call this function to 'convince' this symbolic expression that
        it is in fact simplified.  Basically this means it won't be
        simplified further before printing if you do this.

        This is mainly for internal use.

        EXAMPLES:
        We make an $x + 1 - x$ that prints as that expression with
        no simplification before printing:
            sage: f = x + 1 - x
            sage: f._declare_simplified()
            sage: f
            x + 1 - x
            sage: f._has_been_simplified()
            True

        Note that staying unsimplified as above does not persist
        when we do arithmetic with $f$.
            sage: f + 2
            3
        """
        self._simp = None

    def plot(self, *args, **kwds):
        """
        Plot a symbolic expression.

        All arguments are passed onto the standard plot command.

        EXAMPLES:
        This displays a straight line:
            sage: sin(2).plot((x,0,3))

        This draws a red oscillatory curve:
            sage: sin(x^2).plot((x,0,2*pi), rgbcolor=(1,0,0))

        Another plot using the variable theta:
            sage: var('theta')
            theta
            sage: (cos(theta) - erf(theta)).plot((theta,-2*pi,2*pi))

        A very thick green plot with a frame:
            sage: sin(x).plot((x,-4*pi, 4*pi), thickness=20, rgbcolor=(0,0.7,0)).show(frame=True)

        You can embed 2d plots in 3d space as follows:
            sage: plot(sin(x^2), (x,-pi, pi), thickness=2).plot3d(z = 1)

        A more complicated family:
            sage: G = sum([plot(sin(n*x), (x,-2*pi, 2*pi)).plot3d(z=n) for n in [0,0.1,..1]])
            sage: G.show(frame_aspect_ratio=[1,1,1/2])

        A plot involving the floor function:
            sage: plot(1.0 - x * floor(1/x), (x,0.00001,1.0))
        """
        from sage.plot.plot import plot

        # see if the user passed a variable in.
        if kwds.has_key('param'):
            param = kwds['param']
        else:
            param = None
            for i in range(len(args)):
                if isinstance(args[i], SymbolicVariable):
                    param = args[i]
                    args = args[:i] + args[i+1:]
                    break

        F = self.simplify()
        if isinstance(F, Symbolic_object):
            if hasattr(F._obj, '__call__'):
                f = lambda x: F.subs(x)
            else:
                y = float(F._obj)
                f = lambda x: y

        elif param is None:
            if isinstance(F, CallableSymbolicExpression):
                A = F.arguments()
                if len(A) == 0:
                    raise ValueError, "function has no input arguments"
                else:
                    param = A[0]
                try:
                    f = F._fast_float_(param)
                except NotImplementedError:
                    f = lambda x: F(x)
            else:
                A = F.variables()
                if len(A) == 0:
                    y = float(F)
                    f = lambda x: y
                else:
                    param = A[0]
                    try:
                        f = F._fast_float_(param)
                    except NotImplementedError:
                        return self.function(param)
        else:
            try:
                f = F._fast_float_(param)
            except NotImplementedError:
                return self.function(param)
        return plot(f, *args, **kwds)

    def __lt__(self, right):
        r"""
        Construct the symbolic inequality \code{self < right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:
            sage: x < x
            x < x
            sage: x < Mod(2,5) #random due to architecture dependence
            False
        """
        try:
            return SymbolicEquation(self, SR(right), operator.lt)
        except TypeError:
            return type(self) < type(right)

    def __le__(self, right):
        r"""
        Construct the symbolic inequality \code{self <= right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:
            sage: x <= x
            x <= x
            sage: x <= Mod(2,5) #random due to architecture dependence
            False
            sage: Mod(2,5) >= x #random due to architecture dependence
            False
            sage: Mod(2,5) <= x #random due to architecture dependence
            True
        """
        try:
            return SymbolicEquation(self, SR(right), operator.le)
        except TypeError:
            return type(self) <= type(right)

    def __eq__(self, right):
        r"""
        Construct the symbolic inequality \code{self == right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:

        """
        try:
            return SymbolicEquation(self, SR(right), operator.eq)
        except TypeError:
            return False

    def __ne__(self, right):
        r"""
        Construct the symbolic inequality \code{self != right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:
        """
        try:
            return SymbolicEquation(self, SR(right), operator.ne)
        except TypeError:
            return True

    def __ge__(self, right):
        r"""
        Construct the symbolic inequality \code{self >= right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:
        """
        try:
            return SymbolicEquation(self, SR(right), operator.ge)
        except TypeError:
            return type(self) >= type(right)

    def __gt__(self, right):
        r"""
        Construct the symbolic inequality \code{self > right}.

        NOTE: This only returns a Python bool if right does not coerce
        to the symbolic ring.  Otherwise it returns a symbolic equation.

        EXAMPLES:
        """
        try:
            return SymbolicEquation(self, SR(right), operator.gt)
        except TypeError:
            return type(self) > type(right)

    def __cmp__(self, right):
        """
        Compares self and right.

        This is by definition the comparison of the underlying Maxima
        objects, if right coerces to a symbolic (otherwise types are
        compared).  It is not used unless you explicitly call cmp,
        since all the other special comparison methods are overloaded.

        EXAMPLES:
        These two are equal:
            sage: cmp(e+e, e*2)
            0
            sage: cmp(SR(3), SR(5))
            -1
            sage: cmp(SR(5), SR(2))
            1

        Note that specifiec comparison operators do not call cmp.
            sage: SR(3) < SR(5)
            3 < 5
            sage: bool(SR(3) < SR(5))
            True

        We compare symbolic elements with non symbolic ones.
            sage: cmp(SR(3), 5)
            -1
            sage: cmp(3, SR(5))
            -1

        Here the underlying types are compared, since Mod(2,5)
        doesn't coerce to the symbolic ring.
            sage: cmp(SR(3), Mod(2,5)) #random due to architecture dependence
            1
            sage: cmp(type(SR(3)), type(Mod(2,5))) #random due to architecture dependence
            1
            sage: cmp(Mod(2,5), SR(3) ) #random due to architecture dependence
            -1

        Some comparisons are fairly arbitrary but consistent:
            sage: cmp(SR(3), x) #random due to architecture dependence
            -1
            sage: cmp(x, SR(3)) #random due to architecture dependence
            1
        """
        try:
            right = SR(right)
            return cmp(maxima(self), maxima(right))
        except TypeError:
            return cmp(type(self), type(right))

    def _richcmp_(left, right, op):
        """
        TESTS:
            sage: 3 < x
            3 < x
            sage: 3 <= x
            3 <= x
            sage: 3 == x
            3 == x
            sage: 3 >= x
            3 >= x
            sage: 3 > x
            3 > x
        """
        if op == 0:  #<
            return left < right
        elif op == 2: #==
            return left == right
        elif op == 4: #>
            return left > right
        elif op == 1: #<=
            return left <= right
        elif op == 3: #!=
            return left != right
        elif op == 5: #>=
            return left >= right


    def _neg_(self):
        r"""
        Return the formal negative of \code{self}.

        EXAMPLES:
            sage: var('a,x,y')
            (a, x, y)
            sage: -a
            -a
            sage: -(x+y)
            -y - x
        """
        return SymbolicArithmetic([self], operator.neg)

    ##################################################################
    # Coercions to interfaces
    ##################################################################
    # The maxima one is special:
    def _maxima_(self, session=None):
        r"""
        Method for coercing self as a Maxima \code{RingElement}.
        """
        if session is None:
            return RingElement._maxima_(self, maxima)
        else:
            return RingElement._maxima_(self, session)

    def _maxima_init_(self):
        return self._repr_(simplify=False)


    # The others all go through _sys_init_, which is defined below,
    # and does all interfaces in a unified manner.

    def _axiom_init_(self):
        return self._sys_init_('axiom')

    def _gp_init_(self):
        return self._sys_init_('pari')   # yes, gp goes through pari

    def _maple_init_(self):
        return self._sys_init_('maple')

    def _magma_init_(self):
        return '"%s"'%self.str()

    def _kash_init_(self):
        return self._sys_init_('kash')

    def _macaulay2_init_(self):
        return self._sys_init_('macaulay2')

    def _mathematica_init_(self):
        return self._sys_init_('mathematica')

    def _octave_init_(self):
        return self._sys_init_('octave')

    def _pari_init_(self):
        return self._sys_init_('pari')

    def _sys_init_(self, system):
        return repr(self)

    ##################################################################
    # These systems have no symbolic or numerical capabilities at all,
    # really, so we always just coerce to a string
    ##################################################################
    def _gap_init_(self):
        """
        Conversion of symbolic object to GAP always results in a GAP string.

        EXAMPLES:
            sage: gap(e+pi^2 + x^3)
            x^3 + pi^2 + e
        """
        return '"%s"'%repr(self)

    def _singular_init_(self):
        """
        Conversion of a symbolic object to Singular always results in a Singular string.

        EXAMPLES:
            sage: singular(e+pi^2 + x^3)
            x^3 + pi^2 + e
        """
        return '"%s"'%repr(self)

    ##################################################################
    # Non-canonical coercions to compiled built-in rings and fields
    ##################################################################
    def __int__(self):
        """
        EXAMPLES:
            sage: int(sin(2)*100)
            90
        """
        try:
            return int(repr(self))
        except (ValueError, TypeError):
            return int(float(self))

    def __long__(self):
        """
        EXAMPLES:
            sage: long(sin(2)*100)
            90L
        """
        return long(int(self))

    def numerical_approx(self, prec=None, digits=None):
        r"""
        Return a numerical approximation of \code{self} as either a real or
        complex number with at least the requested number of bits or
        digits of precision.

        NOTE: You can use \code{foo.n()} as a shortcut for
        \code{foo.numerical_approx()}.

        INPUT:
            prec -- an integer: the number of bits of precision
            digits -- an integer: digits of precision

        OUTPUT:
            A RealNumber or ComplexNumber approximation of self with
            prec bits of precision.

        EXAMPLES:
            sage: cos(3).numerical_approx()
            -0.989992496600445

        Use the \code{n()} shortcut:
            sage: cos(3).n()
            -0.989992496600445

        Higher precision:
            sage: cos(3).numerical_approx(200)
            -0.98999249660044545727157279473126130239367909661558832881409
            sage: numerical_approx(cos(3), digits=10)
            -0.9899924966
            sage: (i + 1).numerical_approx(32)
            1.00000000 + 1.00000000*I
            sage: (pi + e + sqrt(2)).numerical_approx(100)
            7.2740880444219335226246195788
        """
        if prec is None:
            if digits is None:
                prec = 53
            else:
                prec = int(digits * 3.4) + 2

        # make sure the field is of the right precision
        field = RealField(prec)

        try:
            approx = self._mpfr_(field)
        except TypeError:
            # try to return a complex result
            approx = self._complex_mpfr_field_(ComplexField(prec))

        return approx

    n = numerical_approx

    def minpoly(self, bits=None, degree=None, epsilon=0):
        r"""
        Return the minimal polynomial of self, if possible.

        INPUT:
            bits    -- the number of bits to use in numerical approx
            degree  -- the expected algebraic degree
            epsilon -- return without error as long as f(self) < epsilon,
                       in the case that the result cannot be proven.

            All of the above parameters are optional, with epsilon=0,
            bits and degree tested up to 1000 and 24 by default respectively.
            If these are known, it will be faster to give them explicitly.

        OUTPUT:
            The minimal polynomial of self. This is proved symbolically if
            epsilon=0 (default).

            If the minimal polynomial could not be found, two distinct kinds
            of errors are raised. If no reasonable candidate was found with
            the given bit/degree parameters, a \exception{ValueError} will be raised.
            If a reasonable candidate was found but (perhaps due to limits
            in the underlying symbolic package) was unable to be proved
            correct, a \exception{NotImplementedError} will be raised.

        ALGORITHM:
            Use the PARI algdep command on a numerical approximation of \code{self}
            to get a candidate minpoly $f$. Approximate $f(\code{self})$ to higher
            precision and if the result is still close enough to 0 then
            evaluate $f(\code{self})$ symbolically, attempting to prove vanishing.
            If this fails, and \var{epsilon} is non-zero, return $f$ as long as
            $f(\code{self}) < \var{epsilon}$. Otherwise raise an error.


        NOTE: Failure of this function does not prove self is
              not algebraic.


        EXAMPLES:

        First some simple examples:
            sage: sqrt(2).minpoly()
            x^2 - 2
            sage: a = 2^(1/3)
            sage: a.minpoly()
            x^3 - 2
            sage: (sqrt(2)-3^(1/3)).minpoly()
            x^6 - 6*x^4 + 6*x^3 + 12*x^2 + 36*x + 1

        Sometimes it fails.
            sage: sin(1).minpoly()
            Traceback (most recent call last):
            ...
            ValueError: Could not find minimal polynomial (1000 bits, degree 24).

        Note that simplification may be necessary.
            sage: a = sqrt(2)+sqrt(3)+sqrt(5)
            sage: f = a.minpoly(); f
            x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
            sage: f(a)
            (sqrt(5) + sqrt(3) + sqrt(2))^2*((sqrt(5) + sqrt(3) + sqrt(2))^2*((sqrt(5) + sqrt(3) + sqrt(2))^2*((sqrt(5) + sqrt(3) + sqrt(2))^2 - 40) + 352) - 960) + 576
            sage: f(a).simplify_radical()
            0

        Here we verify it gives the same result as the abstract number field.
            sage: (sqrt(2) + sqrt(3) + sqrt(6)).minpoly()
            x^4 - 22*x^2 - 48*x - 23
            sage: K.<a,b> = NumberField([x^2-2, x^2-3])
            sage: (a+b+a*b).absolute_minpoly()
            x^4 - 22*x^2 - 48*x - 23

        Works with trig functions too.
            sage: sin(pi/3).minpoly()
            x^2 - 3/4

        Here we show use of the \var{epsilon} parameter. That this result is
        actually exact can be shown using the addition formula for sin,
        but maxima is unable to see that.

            sage: a = sin(pi/5)
            sage: a.minpoly()
            Traceback (most recent call last):
            ...
            NotImplementedError: Could not prove minimal polynomial x^4 - 5/4*x^2 + 5/16 (epsilon 0.00000000000000e-1)
            sage: f = a.minpoly(epsilon=1e-100); f
            x^4 - 5/4*x^2 + 5/16
            sage: f(a).numerical_approx(100)
            0.00000000000000000000000000000

        The degree must be high enough (default tops out at 24).
            sage: a = sqrt(3) + sqrt(2)
            sage: a.minpoly(bits=100, degree=3)
            Traceback (most recent call last):
            ...
            ValueError: Could not find minimal polynomial (100 bits, degree 3).
            sage: a.minpoly(bits=100, degree=10)
            x^4 - 10*x^2 + 1

        Here we solve a cubic and then recover it from its complicated radical expansion.
            sage: f = x^3 - x + 1
            sage: a = f.solve(x)[0].rhs(); a
            (sqrt(3)*I/2 - 1/2)/(3*(sqrt(23)/(6*sqrt(3)) - 1/2)^(1/3)) + (sqrt(23)/(6*sqrt(3)) - 1/2)^(1/3)*(-sqrt(3)*I/2 - 1/2)
            sage: a.minpoly()
            x^3 - x + 1
        """
        bits_list = [bits] if bits else [100,200,500,1000]
        degree_list = [degree] if degree else [2,4,8,12,24]

        for bits in bits_list:
            a = self.numerical_approx(bits)
            check_bits = int(1.25 * bits + 80)
            aa = self.numerical_approx(check_bits)

            for degree in degree_list:

                f = algdep(a, degree) # TODO: use the known_bits parameter?
                # If indeed we have found a minimal polynomial,
                # it should be accurate to a much higher precision.
                error = abs(f(aa))
                dx = ~RR(Integer(1) << (check_bits - degree - 2))
                expected_error = abs(f.derivative()(CC(aa))) * dx

                if error < expected_error:
                    # Degree might have been an over-estimate, factor because we want (irreducible) minpoly.
                    ff = f.factor()
                    for g, e in ff:
                        lead = g.leading_coefficient()
                        if lead != 1:
                            g = g / lead
                        expected_error = abs(g.derivative()(CC(aa))) * dx
                        error = abs(g(aa))
                        if error < expected_error:
                            # See if we can prove equality exactly
                            if g(self).simplify_trig().simplify_radical() == 0:
                                return g
                            # Otherwise fall back to numerical guess
                            elif epsilon and error < epsilon:
                                return g
                            else:
                                raise NotImplementedError, "Could not prove minimal polynomial %s (epsilon %s)" % (g, RR(error).str(no_sci=False))

        raise ValueError, "Could not find minimal polynomial (%s bits, degree %s)." % (bits, degree)

    def _mpfr_(self, field):
        raise TypeError

    def _complex_mpfr_field_(self, field):
        raise TypeError

    def _complex_double_(self, C):
        raise TypeError

    def _real_double_(self, R):
        raise TypeError

    def _real_rqdf_(self, R):
        raise TypeError

    def _rational_(self):
        return Rational(repr(self))

    def __abs__(self):
        return abs_symbolic(self)

    def _integer_(self):
        """
        EXAMPLES:

        """
        return Integer(repr(self))

    def _add_(self, right):
        """
        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: x + y
            y + x
            sage: x._add_(y)
            y + x
        """
        return SymbolicArithmetic([self, right], operator.add)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: x - y
            x - y
        """
        return SymbolicArithmetic([self, right], operator.sub)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: x * y
            x*y
        """
        return SymbolicArithmetic([self, right], operator.mul)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: x / y
            x/y
        """
        return SymbolicArithmetic([self, right], operator.div)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: var('x,n')
            (x, n)
            sage: x^(n+1)
            x^(n + 1)
        """
        right = self.parent()(right)
        return SymbolicArithmetic([self, right], operator.pow)

    def variables(self):
        r"""
        Return sorted list of variables that occur in the simplified
        form of \code{self}.

        OUTPUT:
            a Python set

        EXAMPLES:
            sage: var('x,n')
            (x, n)
            sage: f = x^(n+1) + sin(pi/19); f
            x^(n + 1) + sin(pi/19)
            sage: f.variables()
            (n, x)

            sage: a = e^x
            sage: a.variables()
            (x,)
        """
        return tuple([])

    def arguments(self):
        r"""
        Return the arguments of self, if we view self as a callable
        function.

        This is the same as \code{self.variables()}.

        EXAMPLES:
            sage: x, theta, a = var('x, theta, a')
            sage: f = x^2 + theta^3 - a^x; f
            x^2 + theta^3 - a^x
            sage: f.arguments()
            (a, theta, x)
        """
        return self.variables()

    def number_of_arguments(self):
        """
        Returns the number of arguments the object can take.

        EXAMPLES:
            sage: a,b,c = var('a,b,c')
            sage: foo = function('foo', a,b,c)
            sage: foo.number_of_arguments()
            3
        """
        return len(self.variables())

    def _has_op(self, operator):
        r"""
        Recursively searches for the given operator in a
        \class{SymbolicExpression} object.

        INPUT:
            operator: the operator to search for

        OUTPUT:
            True or False

        EXAMPLES:
            sage: f = 4*(x^2 - 3)
            sage: f._has_op(operator.sub)
            True
            sage: f._has_op(operator.div)
            False
        """

        # if we *are* the operator, then return true right away
        try:
            if operator is self._operator:
                return True
        except AttributeError:
            pass

        # now try to look at this guy's operands
        try:
            ops = self._operands
        # if we don't have an operands, then we can return false
        except AttributeError:
            return False
        for oprnd in ops:
            if oprnd._has_op(operator): return True
            else: pass
        # if we get to this point, neither of the operands have the required
        # operator
        return False


    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
            sage: x,y=var('x,y')
            sage: f = x+y
            sage: f.arguments()
            (x, y)
            sage: f()
            y + x
            sage: f(3)
            y + 3
            sage: f(3,4)
            7
            sage: f(2,3,4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 2

            sage: f({x:3})
            y + 3
            sage: f({x:3,y:4})
            7
            sage: f(x=3)
            y + 3
            sage: f(x=3,y=4)
            7

            sage: a = (2^(8/9))
            sage: a(4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0


            sage: f = function('Gamma', var('z'), var('w')); f
            Gamma(z, w)
            sage: f(2)
            Gamma(2, w)
            sage: f(2,5)
            Gamma(2, 5)

        """
        if len(args) == 0:
            d = None
        elif len(args) == 1 and isinstance(args[0], dict):
            d = args[0]
        else:
            d = {}
            vars = self.arguments()
            for i in range(len(args)):
                try:
                    d[ vars[i] ] = args[i]
                except IndexError:
                    raise ValueError, "the number of arguments must be less than or equal to %s"%len(self.variables())

        return self.substitute(d, **kwds)

    def power_series(self, base_ring):
        """
        Return algebraic power series associated to this symbolic
        expression, which must be a polynomial in one variable, with
        coefficients coercible to the base ring.

        The power series is truncated one more than the degree.

        EXAMPLES:
            sage: theta = var('theta')
            sage: f = theta^3 + (1/3)*theta - 17/3
            sage: g = f.power_series(QQ); g
            -17/3 + 1/3*theta + theta^3 + O(theta^4)
            sage: g^3
            -4913/27 + 289/9*theta - 17/9*theta^2 + 2602/27*theta^3 + O(theta^4)
            sage: g.parent()
            Power Series Ring in theta over Rational Field
        """
        v = self.variables()
        if len(v) != 1:
            raise ValueError, "self must be a polynomial in one variable but it is in the variables %s"%tuple([v])
        f = self.polynomial(base_ring)
        from sage.rings.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, names=f.parent().variable_names())
        return R(f, f.degree()+1)

    def polynomial(self, base_ring):
        r"""
        Return \code{self} as an algebraic polynomial over the given
        base ring, if possible.

        The point of this function is that it converts purely symbolic
        polynomials into optimized algebraic polynomials over a given
        base ring.

        WARNING: This is different from \code{self.poly(x)} which is used
        to rewrite self as a polynomial in x.

        INPUT:
           base_ring -- a ring

        EXAMPLES:
            sage: f = x^2 -2/3*x + 1
            sage: f.polynomial(QQ)
            x^2 - 2/3*x + 1
            sage: f.polynomial(GF(19))
            x^2 + 12*x + 1

        Polynomials can be useful for getting the coefficients
        of an expression:
            sage: g = 6*x^2 - 5
            sage: g.coefficients()
            [[-5, 0], [6, 2]]
            sage: g.polynomial(QQ).list()
            [-5, 0, 6]
            sage: g.polynomial(QQ).dict()
            {0: -5, 2: 6}

            sage: f = x^2*e + x + pi/e
            sage: f.polynomial(RDF)
            2.71828182846*x^2 + 1.0*x + 1.15572734979
            sage: g = f.polynomial(RR); g
            2.71828182845905*x^2 + 1.00000000000000*x + 1.15572734979092
            sage: g.parent()
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: f.polynomial(RealField(100))
            2.7182818284590452353602874714*x^2 + 1.0000000000000000000000000000*x + 1.1557273497909217179100931833
            sage: f.polynomial(CDF)
            2.71828182846*x^2 + 1.0*x + 1.15572734979
            sage: f.polynomial(CC)
            2.71828182845905*x^2 + 1.00000000000000*x + 1.15572734979092

        We coerce a multivariate polynomial with complex symbolic coefficients:
            sage: x, y, n = var('x, y, n')
            sage: f = pi^3*x - y^2*e - I; f
            -1*e*y^2 + pi^3*x - I
            sage: f.polynomial(CDF)
            (-2.71828182846)*y^2 + 31.0062766803*x - 1.0*I
            sage: f.polynomial(CC)
            (-2.71828182845905)*y^2 + 31.0062766802998*x - 1.00000000000000*I
            sage: f.polynomial(ComplexField(70))
            (-2.7182818284590452354)*y^2 + 31.006276680299820175*x - 1.0000000000000000000*I

        Another polynomial:
            sage: f = sum((e*I)^n*x^n for n in range(5)); f
            e^4*x^4 - e^3*I*x^3 - e^2*x^2 + e*I*x + 1
            sage: f.polynomial(CDF)
            54.5981500331*x^4 + (-20.0855369232*I)*x^3 + (-7.38905609893)*x^2 + 2.71828182846*I*x + 1.0
            sage: f.polynomial(CC)
            54.5981500331442*x^4 + (-20.0855369231877*I)*x^3 + (-7.38905609893065)*x^2 + 2.71828182845905*I*x + 1.00000000000000

        A multivariate polynomial over a finite field:
            sage: f = (3*x^5 - 5*y^5)^7; f
            (3*x^5 - 5*y^5)^7
            sage: g = f.polynomial(GF(7)); g
            3*x^35 + 2*y^35
            sage: parent(g)
            Multivariate Polynomial Ring in x, y over Finite Field of size 7
        """
        vars = self.variables()
        if len(vars) == 0:
            vars = ['x']
        R = PolynomialRing(base_ring, names=vars)
        G = R.gens()
        V = R.variable_names()
        return self.substitute_over_ring(
             dict([(var(V[i]),G[i]) for i in range(len(G))]), ring=R)

    def _polynomial_(self, R):
        """
        Coerce this symbolic expression to a polynomial in $R$.

        EXAMPLES:
            sage: var('x,y,z,w')
            (x, y, z, w)

            sage: R = QQ[x,y,z]
            sage: R(x^2 + y)
            x^2 + y
            sage: R = QQ[w]
            sage: R(w^3 + w + 1)
            w^3 + w + 1
            sage: R = GF(7)[z]
            sage: R(z^3 + 10*z)
            z^3 + 3*z

        NOTE: If the base ring of the polynomial ring is the symbolic
        ring, then a constant polynomial is always returned.
            sage: R = SR[x]
            sage: a = R(sqrt(2) + x^3 + y)
            sage: a
            y + x^3 + sqrt(2)
            sage: type(a)
            <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
            sage: a.degree()
            0

        We coerce to a double precision complex polynomial ring:
            sage: f = e*x^3 + pi*y^3 + sqrt(2) + I; f
            pi*y^3 + e*x^3 + I + sqrt(2)
            sage: R = CDF[x,y]
            sage: R(f)
            2.71828182846*x^3 + 3.14159265359*y^3 + 1.41421356237 + 1.0*I

        We coerce to a higher-precision polynomial ring
            sage: R = ComplexField(100)[x,y]
            sage: R(f)
            2.7182818284590452353602874714*x^3 + 3.1415926535897932384626433833*y^3 + 1.4142135623730950488016887242 + 1.0000000000000000000000000000*I
        """
        vars = self.variables()
        B = R.base_ring()
        if B == SR:
            if is_MPolynomialRing(R):
                return R({tuple([0]*R.ngens()):self})
            else:
                return R([self])
        G = R.gens()
        sub = []
        for v in vars:
            r = repr(v)
            for g in G:
                if repr(g) == r:
                    sub.append((v,g))
        if len(sub) == 0:
            try:
                return R(B(self))
            except TypeError:
                if len(vars) == 1:
                    sub = [(vars[0], G[0])]
                else:
                    raise
        return self.substitute_over_ring(dict(sub), ring=R)

    def function(self, *args):
        r"""
        Return a \class{CallableSymbolicExpression}, fixing a variable order
        to be the order of args.

        EXAMPLES:
        We will use several symbolic variables in the examples below:
           sage: var('x, y, z, t, a, w, n')
           (x, y, z, t, a, w, n)

           sage: u = sin(x) + x*cos(y)
           sage: g = u.function(x,y)
           sage: g(x,y)
           x*cos(y) + sin(x)
           sage: g(t,z)
           t*cos(z) + sin(t)
           sage: g(x^2, x^y)
           x^2*cos(x^y) + sin(x^2)

            sage: f = (x^2 + sin(a*w)).function(a,x,w); f
            (a, x, w) |--> x^2 + sin(a*w)
            sage: f(1,2,3)
            sin(3) + 4

        Using the \method{function} method we can obtain the above function $f$,
        but viewed as a function of different variables:
            sage: h = f.function(w,a); h
            (w, a) |--> x^2 + sin(a*w)

        This notation also works:
            sage: h(w,a) = f
            sage: h
            (w, a) |--> x^2 + sin(a*w)

        You can even make a symbolic expression $f$ into a function by
        writing \code{f(x,y) = f}:
            sage: f = x^n + y^n; f
            y^n + x^n
            sage: f(x,y) = f
            sage: f
            (x, y) |--> y^n + x^n
            sage: f(2,3)
            3^n + 2^n
        """
        R = CallableSymbolicExpressionRing(args)
        return R(self)


    ###################################################################
    # derivative
    ###################################################################

    def _derivative(self, var=None):
        r"""
        Derivative of self with respect to var (a symbolic variable).

        If var is None, self must contain only one variable, and the result
        is the derivative with respect to that variable.

        If var does not appear in self, the result is zero.

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: x = var("x"); y = var("y"); z = var("z")
            sage: f = sin(x) * cos(y)
            sage: f._derivative(x)
            cos(x)*cos(y)
            sage: f._derivative(y)
            -sin(x)*sin(y)
            sage: f._derivative(z)
            0
            sage: f._derivative()
            Traceback (most recent call last):
            ...
            ValueError: must supply an explicit variable for an expression containing more than one variable

            sage: f = sin(x)
            sage: f._derivative()
            cos(x)

            sage: f._derivative(2)
            Traceback (most recent call last):
            ...
            TypeError: arguments must be SymbolicVariable objects

        """
        if var is None:
            # use default variable, if we can figure it out
            vars = self.variables()
            if len(vars) > 1:
                raise ValueError, "must supply an explicit variable for an " +\
                                  "expression containing more than one variable"
            if len(vars) == 0:
                # no variables in expression, derivative must be zero
                return self.parent()(0)
            var = vars[0]

        elif not isinstance(var, SymbolicVariable):
            raise TypeError, "arguments must be SymbolicVariable objects"

        t = maxima('diff(%s, %s)' % (self._maxima_().name(), repr(var)))
        f = self.parent()(t)
        return f


    def derivative(self, *args):
        r"""
        Derivative with respect to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        SEE ALSO:
            self._derivative()

        EXAMPLES:
            sage: h = sin(x)/cos(x)
            sage: derivative(h,x,x,x)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2
            sage: derivative(h,x,3)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2

            sage: var('x, y')
            (x, y)
            sage: u = (sin(x) + cos(y))*(cos(x) - sin(y))
            sage: derivative(u,x,y)
            sin(x)*sin(y) - cos(x)*cos(y)
            sage: f = ((x^2+1)/(x^2-1))^(1/4)
            sage: g = derivative(f, x); g # this is a complex expression
            x/(2*(x^2 - 1)^(1/4)*(x^2 + 1)^(3/4)) - x*(x^2 + 1)^(1/4)/(2*(x^2 - 1)^(5/4))
            sage: g.simplify_rational()
            -x/((x^2 - 1)^(5/4)*(x^2 + 1)^(3/4))

            sage: f = y^(sin(x))
            sage: derivative(f, x)
            cos(x)*y^sin(x)*log(y)

            sage: g(x) = sqrt(5-2*x)
            sage: g_3 = derivative(g, x, 3); g_3(2)
            -3

            sage: f = x*e^(-x)
            sage: derivative(f, 100)
            x*e^(-x) - 100*e^(-x)

            sage: g = 1/(sqrt((x^2-1)*(x+5)^6))
            sage: derivative(g, x)
            -3*(x + 5)^5/(((x + 5)^6)^(3/2)*sqrt(x^2 - 1)) - x/(sqrt((x + 5)^6)*(x^2 - 1)^(3/2))
        """
        # note: it would be simpler to use multi_derivative() here instead of all
        # the code below. The reason we do it this way is to reduce the number of
        # calls to maxima wherever possible.

        args = derivative_parse(args)

        if not args:
            # no differentation taking place
            return self

        # check all variables are really variables
        for arg in args:
            if arg is not None and not isinstance(arg, SymbolicVariable):
                raise TypeError, "arguments must be SymbolicVariable objects"

        vars = self.variables()

        if len(vars) == 0:
            # self has no variables, so result must be zero
            return self.parent()(0)

        if len(vars) == 1:
            # self has exactly one variable. If the argument list contains
            # a different variable somewhere, the result has to be zero.
            for arg in args:
                if arg is not None and arg is not vars[0]:
                    return self.parent()(0)

            # otherwise we're differentating with respect to that
            # variable n times (each None may be assumed to correspond
            # to that variable).
            t = maxima('diff(%s, %s, %d)' % (self._maxima_().name(), repr(vars[0]), len(args)))
            return self.parent()(t)

        # There's more than one variable in self.
        # If None appears anywhere in args, we'll just have to differentiate
        # one step at a time, since we can't tell in advance what the "default"
        # variable is going to be.
        if None in args:
            F = self
            for arg in args:
                F = F._derivative(arg)
            return F

        # The list of arguments is completely explicit, so we do it in
        # a single maxima call
        s = ""
        for arg in args:
            s = s + ", " + repr(arg) + ", 1"
        t = maxima('diff(%s%s)' % (self._maxima_().name(), s))
        return self.parent()(t)


    differentiate = derivative
    diff = derivative

    def gradient(self):
        r"""
        Compute the gradient of a symbolic function.
        This function returns a vector whose components are the
        derivatives of the original function.

        EXAMPLES:
            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.gradient()
            (2*x, 2*y)
        """

        from sage.modules.free_module_element import vector
        l=[self.derivative(x) for x in self.variables()]
        return vector(l)

    def hessian(self):
        r"""
        Compute the hessian of a function. This returns a matrix
        components are the 2nd partial derivatives of the original function.

        EXAMPLES:
            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.hessian()
            [2 0]
            [0 2]

        """

        from sage.matrix  import constructor
        grad=self.gradient()
        var_list=self.variables()
        l=[ [grad[i].derivative(x) for x in var_list] for i in xrange(len(grad))]
        return constructor.matrix(l)



    ###################################################################
    # Taylor series
    ###################################################################
    def taylor(self, v, a, n):
        r"""
        Expands \code{self} in a truncated Taylor or Laurent series in
        the variable $v$ around the point $a$, containing terms
        through $(x - a)^n$.

        INPUT:
            v -- variable
            a -- number
            n -- integer

        EXAMPLES:
            sage: var('a, x, z')
            (a, x, z)
            sage: taylor(a*log(z), z, 2, 3)
            log(2)*a + a*(z - 2)/2 - a*(z - 2)^2/8 + a*(z - 2)^3/24
            sage: taylor(sqrt (sin(x) + a*x + 1), x, 0, 3)
            1 + (a + 1)*x/2 - (a^2 + 2*a + 1)*x^2/8 + (3*a^3 + 9*a^2 + 9*a - 1)*x^3/48
            sage: taylor (sqrt (x + 1), x, 0, 5)
            1 + x/2 - x^2/8 + x^3/16 - 5*x^4/128 + 7*x^5/256
            sage: taylor (1/log (x + 1), x, 0, 3)
            1/x + 1/2 - x/12 + x^2/24 - 19*x^3/720
            sage: taylor (cos(x) - sec(x), x, 0, 5)
            -x^2 - x^4/6
            sage: taylor ((cos(x) - sec(x))^3, x, 0, 9)
            -x^6 - x^8/2
            sage: taylor (1/(cos(x) - sec(x))^3, x, 0, 5)
            -1/x^6 + 1/(2*x^4) + 11/(120*x^2) - 347/15120 - 6767*x^2/604800 - 15377*x^4/7983360
        """
        v = var(v)
        l = self._maxima_().taylor(v, SR(a), Integer(n))
        return self.parent()(l)

    ###################################################################
    # limits
    ###################################################################
    def limit(self, dir=None, taylor=False, **argv):
        r"""
        Return the limit as the variable $v$ approaches $a$ from the
        given direction.

        \begin{verbatim}
        expr.limit(x = a)
        expr.limit(x = a, dir='above')
        \end{verbatim}

        INPUT:
            dir -- (default: None); dir may have the value `plus' (or 'above')
                   for a limit from above, `minus' (or 'below') for a limit from
                   below, or may be omitted (implying a two-sided
                   limit is to be computed).
            taylor -- (default: False); if True, use Taylor series, which
                   allows more integrals to be computed (but may also crash
                   in some obscure cases due to bugs in Maxima).
            **argv -- 1 named parameter

        NOTE: The output may also use `und' (undefined), `ind'
        (indefinite but bounded), and `infinity' (complex infinity).

        EXAMPLES:
            sage: f = (1+1/x)^x
            sage: f.limit(x = oo)
            e
            sage: f.limit(x = 5)
            7776/3125
            sage: f.limit(x = 1.2)
            2.0696157546720...
            sage: f.limit(x = I, taylor=True)
            (1 - I)^I
            sage: f(1.2)
            2.0696157546720...
            sage: f(I)
            (1 - I)^I
            sage: CDF(f(I))
            2.06287223508 + 0.74500706218*I
            sage: CDF(f.limit(x = I))
            2.06287223508 + 0.74500706218*I

        More examples:
            sage: limit(x*log(x), x = 0, dir='above')
            0
            sage: lim((x+1)^(1/x),x = 0)
            e
            sage: lim(e^x/x, x = oo)
            +Infinity
            sage: lim(e^x/x, x = -oo)
            0
            sage: lim(-e^x/x, x = oo)
            -Infinity
            sage: lim((cos(x))/(x^2), x = 0)
            +Infinity
            sage: lim(sqrt(x^2+1) - x, x = oo)
            0
            sage: lim(x^2/(sec(x)-1), x=0)
            2
            sage: lim(cos(x)/(cos(x)-1), x=0)
            -Infinity
            sage: lim(x*sin(1/x), x=0)
            0

            sage: f = log(log(x))/log(x)
            sage: forget(); assume(x<-2); lim(f, x=0, taylor=True)
            und

        Here ind means "indefinite but bounded":
            sage: lim(sin(1/x), x = 0)
            ind
        """
        if len(argv) != 1:
            raise ValueError, "call the limit function like this, e.g. limit(expr, x=2)."
        else:
            k = argv.keys()[0]
            v = var(k, create=False)
            a = argv[k]
        if dir is None:
            if taylor:
                l = self._maxima_().tlimit(v, a)
            else:
                l = self._maxima_().limit(v, a)
        elif dir == 'plus' or dir == 'above':
            if taylor:
                l = self._maxima_().tlimit(v, a, 'plus')
            else:
                l = self._maxima_().limit(v, a, 'plus')
        elif dir == 'minus' or dir == 'below':
            if taylor:
                l = self._maxima_().tlimit(v, a, 'minus')
            else:
                l = self._maxima_().limit(v, a, 'minus')
        else:
            raise ValueError, "dir must be one of 'plus' or 'minus'"
        return self.parent()(l)

    ###################################################################
    # Laplace transform
    ###################################################################
    def laplace(self, t, s):
        r"""
        Attempts to compute and return the Laplace transform of
        \code{self} with respect to the variable $t$ and transform
        parameter $s$.  If this function cannot find a solution, a
        formal function is returned.

        The function that is returned may be be viewed as a function
        of $s$.

        DEFINITION:
        The Laplace transform of a function $f(t)$, defined for all
        real numbers $t \geq 0$, is the function $F(s)$ defined by
        $$
             F(s) = \int_{0}^{\infty} e^{-st} f(t) dt.
        $$

        EXAMPLES:
        We compute a few Laplace transforms:
            sage: var('x, s, z, t, t0')
            (x, s, z, t, t0)
            sage: sin(x).laplace(x, s)
            1/(s^2 + 1)
            sage: (z + exp(x)).laplace(x, s)
            z/s + 1/(s - 1)
            sage: log(t/t0).laplace(t, s)
            (-log(t0) - log(s) - euler_gamma)/s

        We do a formal calculation:
            sage: f = function('f', x)
            sage: g = f.diff(x); g
            diff(f(x), x, 1)
            sage: g.laplace(x, s)
            s*laplace(f(x), x, s) - f(0)

        EXAMPLE: A BATTLE BETWEEN the X-women and the Y-men (by David Joyner):
        Solve
        $$
          x' = -16y, x(0)=270,  y' = -x + 1, y(0) = 90.
        $$
        This models a fight between two sides, the "X-women"
        and the "Y-men", where the X-women have 270 initially and
        the Y-men have 90, but the Y-men are better at fighting,
        because of the higher factor of "-16" vs "-1", and also get
        an occasional reinforcement, because of the "+1" term.

            sage: var('t')
            t
            sage: t = var('t')
            sage: x = function('x', t)
            sage: y = function('y', t)
            sage: de1 = x.diff(t) + 16*y
            sage: de2 = y.diff(t) + x - 1
            sage: de1.laplace(t, s)
            16*laplace(y(t), t, s) + s*laplace(x(t), t, s) - x(0)
            sage: de2.laplace(t, s)
            s*laplace(y(t), t, s) + laplace(x(t), t, s) - 1/s - y(0)

        Next we form the augmented matrix of the above system:
            sage: A = matrix([[s, 16, 270],[1, s, 90+1/s]])
            sage: E = A.echelon_form()
            sage: xt = E[0,2].inverse_laplace(s,t)
            sage: yt = E[1,2].inverse_laplace(s,t)
            sage: print xt
				4 t	    - 4 t
			   91  e      629  e
        		 - -------- + ----------- + 1
			      2		   2
            sage: print yt
				 4 t	     - 4 t
			    91  e      629  e
        		    -------- + -----------
			       8	    8
            sage: p1 = plot(xt,0,1/2,rgbcolor=(1,0,0))
            sage: p2 = plot(yt,0,1/2,rgbcolor=(0,1,0))
            sage: (p1+p2).save()
        """
        return self.parent()(self._maxima_().laplace(var(t), var(s)))

    def inverse_laplace(self, t, s):
        r"""
        Attempts to compute the inverse Laplace transform of
        \code{self} with respect to the variable $t$ and transform
        parameter $s$.  If this function cannot find a solution, a
        formal function is returned.

        The function that is returned may be be viewed as a function
        of $s$.

        DEFINITION:
        The inverse Laplace transform of a function $F(s)$,
        is the function $f(t)$ defined by
        $$
             F(s) = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma + i\infty} e^{st} F(s) dt,
        $$
        where $\gamma$ is chosen so that the contour path of
        integration is in the region of convergence of $F(s)$.

        EXAMPLES:
            sage: var('w, m')
            (w, m)
            sage: f = (1/(w^2+10)).inverse_laplace(w, m); f
            sin(sqrt(10)*m)/sqrt(10)
            sage: laplace(f, m, w)
            1/(w^2 + 10)
        """
        return self.parent()(self._maxima_().ilt(var(t), var(s)))


    ###################################################################
    # a default variable, for convenience.
    ###################################################################
    def default_variable(self):
        """
        Return the default variable, which is by definition the first
        variable in self, or $x$ is there are no variables in self.
        The result is cached.

        EXAMPLES:
            sage: sqrt(2).default_variable()
            x
            sage: x, theta, a = var('x, theta, a')
            sage: f = x^2 + theta^3 - a^x
            sage: f.default_variable()
            a

        Note that this is the first \emph{variable}, not
        the first \emph{argument}:
            sage: f(theta, a, x) = a + theta^3
            sage: f.default_variable()
            a
            sage: f.variables()
            (a, theta)
            sage: f.arguments()
            (theta, a, x)
        """
        try:
            return self.__default_variable
        except AttributeError:
            pass
        v = self.variables()
        if len(v) == 0:
            ans = var('x')
        else:
            ans = v[0]
        self.__default_variable = ans
        return ans


    ###################################################################
    # integration
    ###################################################################
    def integral(self, v=None, a=None, b=None):
        r"""
        Returns the indefinite integral with respect to the variable
        $v$, ignoring the constant of integration. Or, if endpoints
        $a$ and $b$ are specified, returns the definite integral over
        the interval $[a, b]$.

        If \code{self} has only one variable, then it returns the
        integral with respect to that variable.

        INPUT:
            v -- (optional) a variable or variable name
            a -- (optional) lower endpoint of definite integral
            b -- (optional) upper endpoint of definite integral


        EXAMPLES:
            sage: h = sin(x)/(cos(x))^2
            sage: h.integral(x)
            1/cos(x)

            sage: f = x^2/(x+1)^3
            sage: f.integral()
            log(x + 1) + (4*x + 3)/(2*x^2 + 4*x + 2)

            sage: f = x*cos(x^2)
            sage: f.integral(x, 0, sqrt(pi))
            0
            sage: f.integral(a=-pi, b=pi)
            0

            sage: f(x) = sin(x)
            sage: f.integral(x, 0, pi/2)
            1

        Constraints are sometimes needed:
            sage: var('x, n')
            (x, n)
            sage: integral(x^n,x)
            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (use assume):
            Is  n+1  zero or nonzero?
            sage: assume(n > 0)
            sage: integral(x^n,x)
            x^(n + 1)/(n + 1)
            sage: forget()


        Note that an exception is raised when a definite integral is divergent.
            sage: integrate(1/x^3,x,0,1)
            Traceback (most recent call last):
            ...
            ValueError: Integral is divergent.
            sage: integrate(1/x^3,x,-1,3)
            Traceback (most recent call last):
            ...
            ValueError: Integral is divergent.

        NOTE: Above, putting assume(n == -1) does not yield the right behavior.
        Directly in maxima, doing

        The examples in the Maxima documentation:
            sage: var('x, y, z, b')
            (x, y, z, b)
            sage: integral(sin(x)^3)
            cos(x)^3/3 - cos(x)
            sage: integral(x/sqrt(b^2-x^2))
            x*log(2*sqrt(b^2 - x^2) + 2*b)
            sage: integral(x/sqrt(b^2-x^2), x)
            -sqrt(b^2 - x^2)
            sage: integral(cos(x)^2 * exp(x), x, 0, pi)
            3*e^pi/5 - 3/5
            sage: integral(x^2 * exp(-x^2), x, -oo, oo)
            sqrt(pi)/2

        We integrate the same function in both Mathematica and \sage (via Maxima):
            sage: f = sin(x^2) + y^z
            sage: g = mathematica(f)                           # optional  -- requires mathematica
            sage: print g                                      # optional
                      z        2
                     y  + Sin[x ]
            sage: print g.Integrate(x)                         # optional
                        z        Pi                2
                     x y  + Sqrt[--] FresnelS[Sqrt[--] x]
                                 2                 Pi
            sage: print f.integral(x)
                  z                                         (sqrt(2)  I + sqrt(2)) x
               x y  + sqrt( pi) ((sqrt(2)  I + sqrt(2)) erf(------------------------)
                                                                       2
                                                           (sqrt(2)  I - sqrt(2)) x
                              + (sqrt(2)  I - sqrt(2)) erf(------------------------))/8
                                                                      2

        We integrate the above function in maple now:
            sage: g = maple(f); g                             # optional -- requires maple
            sin(x^2)+y^z
            sage: g.integrate(x)                              # optional -- requires maple
            1/2*2^(1/2)*Pi^(1/2)*FresnelS(2^(1/2)/Pi^(1/2)*x)+y^z*x

        We next integrate a function with no closed form integral.  Notice that
        the answer comes back as an expression that contains an integral itself.
            sage: A = integral(1/ ((x-4) * (x^3+2*x+1)), x); A
            log(x - 4)/73 - integrate((x^2 + 4*x + 18)/(x^3 + 2*x + 1), x)/73
            sage: print A
                                     /  2
                                     [ x  + 4 x + 18
                                     I ------------- dx
                                     ]  3
                        log(x - 4)   / x  + 2 x + 1
                        ---------- - ------------------
                            73               73

        We now show that floats are not converted to rationals
        automatically since we by default have keepfloat: true in
        maxima.

            sage: integral(e^(-x^2),x, 0, 0.1)
            0.0562314580091424*sqrt(pi)


        ALIASES:
            integral() and integrate() are the same.


        EXAMPLES:
        Here is example where we have to use assume:
            sage: a,b = var('a,b')
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (use assume):
            Is  a  positive or negative?

        So we just assume that $a>0$ and the integral works:
            sage: assume(a>0)
            sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
            2*b^2*arctan((2*(b*x + a)^(1/3) + a^(1/3))/(sqrt(3)*a^(1/3)))/(3*sqrt(3)*a^(7/3)) - b^2*log((b*x + a)^(2/3) + a^(1/3)*(b*x + a)^(1/3) + a^(2/3))/(9*a^(7/3)) + 2*b^2*log((b*x + a)^(1/3) - a^(1/3))/(9*a^(7/3)) + (4*b^2*(b*x + a)^(5/3) - 7*a*b^2*(b*x + a)^(2/3))/(6*a^2*(b*x + a)^2 - 12*a^3*(b*x + a) + 6*a^4)
        """

        if v is None:
            v = self.default_variable()

        if not isinstance(v, SymbolicVariable):
            v = var(repr(v))
            #raise TypeError, 'must integrate with respect to a variable'
        if (a is None and (not b is None)) or (b is None and (not a is None)):
            raise TypeError, 'only one endpoint given'
        if a is None:
            return self.parent()(self._maxima_().integrate(v))
        else:
            try:
                return self.parent()(self._maxima_().integrate(v, a, b))
            except TypeError, error:
                s = str(error)
                if "divergent" in s or 'Principal Value' in s:
                    raise ValueError, "Integral is divergent."
                else:
                    raise TypeError, error


    integrate = integral

    def nintegral(self, x, a, b,
                  desired_relative_error='1e-8',
                  maximum_num_subintervals=200):
        r"""
        Return a floating point machine precision numerical
        approximation to the integral of \code{self} from $a$ to $b$, computed
        using floating point arithmetic via maxima.

        INPUT:
            x -- variable to integrate with respect to
            a -- lower endpoint of integration
            b -- upper endpoint of integration
            desired_relative_error -- (default: '1e-8') the desired
                 relative error
            maximum_num_subintervals -- (default: 200) maxima number
                 of subintervals

        OUTPUT:
            -- float: approximation to the integral
            -- float: estimated absolute error of the approximation
            -- the number of integrand evaluations
            -- an error code:
                  0 -- no problems were encountered
                  1 -- too many subintervals were done
                  2 -- excessive roundoff error
                  3 -- extremely bad integrand behavior
                  4 -- failed to converge
                  5 -- integral is probably divergent or slowly convergent
                  6 -- the input is invalid

        ALIAS:
            nintegrate is the same as nintegral

        REMARK:
            There is also a function \code{numerical_integral} that implements
            numerical integration using the GSL C library.  It is potentially
            much faster and applies to arbitrary user defined functions.

            Also, there are limits to the precision to which Maxima can compute
            the integral to due to limitations in quadpack.

            sage: f = x
            sage: f = f.nintegral(x,0,1,1e-14)
            Traceback (most recent call last):
            ...
            ValueError: Maxima (via quadpack) cannot compute the integral to that precision

        EXAMPLES:
            sage: f(x) = exp(-sqrt(x))
            sage: f.nintegral(x, 0, 1)
            (0.52848223531423055, 4.163...e-11, 231, 0)

        We can also use the \code{numerical_integral} function, which calls
        the GSL C library.
            sage: numerical_integral(f, 0, 1)       # random low-order bits
            (0.52848223225314706, 6.8392846084921134e-07)

        Note that in exotic cases where floating point evaluation of
        the expression leads to the wrong value, then the output
        can be completely wrong:
            sage: f = exp(pi*sqrt(163)) - 262537412640768744

        Despite appearance, $f$ is really very close to 0, but one
        gets a nonzero value since the definition of \code{float(f)} is
        that it makes all constants inside the expression floats, then
        evaluates each function and each arithmetic operation
        using float arithmetic:
            sage: float(f)
            -480.0

        Computing to higher precision we see the truth:
            sage: f.n(200)
            -7.4992740280181431112064614366622348652078895136533593355718e-13
            sage: f.n(300)
            -7.49927402801814311120646143662663009137292462589621789352095066181709095575681963967103004e-13

        Now numerically integrating, we see why the answer is wrong:
            sage: f.nintegrate(x,0,1)
            (-480.00000000000011, 5.3290705182007538e-12, 21, 0)

        It is just because every floating point evaluation of return
        -480.0 in floating point.

        Important note: using GP/PARI one can compute numerical
        integrals to high precision:
            sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
            '2.565728500561051482917356396 E-127'        # 32-bit
            '2.5657285005610514829173563961304785900 E-127'    # 64-bit
            sage: old_prec = gp.set_real_precision(50)
            sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
            '2.5657285005610514829173563961304785900147709554020 E-127'
            sage: gp.set_real_precision(old_prec)
            57

        Note that the input function above is a string in PARI
        syntax.
        """
        try:
            v = self._maxima_().quad_qags(var(x),
                                      a, b, desired_relative_error,
                                      maximum_num_subintervals)
        except TypeError, err:
            if "ERROR" in str(err):
                raise ValueError, "Maxima (via quadpack) cannot compute the integral to that precision"
            else:
                raise TypeError, err

        return float(v[0]), float(v[1]), Integer(v[2]), Integer(v[3])

    nintegrate = nintegral


    ###################################################################
    # Manipulating epxressions
    ###################################################################
    def coefficient(self, x, n=1):
        """
        Returns the coefficient of $x^n$ in self.

        INPUT:
            x -- variable, function, expression, etc.
            n -- integer, default 1.

        Sometimes it may be necessary to expand or factor first, since
        this is not done automatically.

        EXAMPLES:
            sage: var('a, x, y, z')
            (a, x, y, z)
            sage: f = (a*sqrt(2))*x^2 + sin(y)*x^(1/2) + z^z
            sage: f.coefficient(sin(y))
            sqrt(x)
            sage: f.coefficient(x^2)
            sqrt(2)*a
            sage: f.coefficient(x^(1/2))
            sin(y)
            sage: f.coefficient(1)
            0
            sage: f.coefficient(x, 0)
            z^z
        """
        return self.parent()(self._maxima_().coeff(x, n))
    coeff = coefficient

    def coefficients(self, x=None):
        r"""
        Coefficients of \code{self} as a polynomial in x.

        INPUT:
            x -- optional variable
        OUTPUT:
            list of pairs [expr, n], where expr is a symbolic
            expression and n is a power.

        EXAMPLES:
            sage: var('x, y, a')
            (x, y, a)
            sage: p = x^3 - (x-3)*(x^2+x) + 1
            sage: p.coefficients()
            [[1, 0], [3, 1], [2, 2]]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            x^2 - 2*sqrt(2)*a*x + x + 2*a^2 + 1
            sage: p.coefficients(a)
            [[x^2 + x + 1, 0], [-2*sqrt(2)*x, 1], [2, 2]]
            sage: p.coefficients(x)
            [[2*a^2 + 1, 0], [1 - 2*sqrt(2)*a, 1], [1, 2]]

        A polynomial with wacky exponents:
            sage: p = (17/3*a)*x^(3/2) + x*y + 1/x + x^x
            sage: p.coefficients(x)
            [[1, -1], [x^x, 0], [y, 1], [17*a/3, 3/2]]
        """
        f = self._maxima_()
        P = f.parent()
        P._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = var(x)
        G = f.coeffs(x)
        S = symbolic_expression_from_maxima_string(repr(G))
        return S[1:]

    coeffs = coefficients

    def poly(self, x=None):
        r"""
        Express \code{self} as a polynomial in $x$.  If \code{self} is not a polynomial
        in $x$, then some coefficients may be functions of $x$.

        WARNING: This is different from \code{self.polynomial()} which
        returns a \sage polynomial over a given base ring.

        EXAMPLES:
            sage: var('a, x')
            (a, x)
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            x^2 - 2*sqrt(2)*a*x + x + 2*a^2 + 1
            sage: p.poly(a)
            (x^2 + x + 1)*1 + -2*sqrt(2)*x*a + 2*a^2
            sage: bool(expand(p.poly(a)) == p)
            True
            sage: p.poly(x)
            (2*a^2 + 1)*1 + (1 - 2*sqrt(2)*a)*x + 1*x^2
        """
        f = self._maxima_()
        P = f.parent()
        P._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = var(x)
        G = f.coeffs(x)
        ans = None
        for i in range(1, len(G)):
            Z = G[i]
            coeff = SR(Z[0])
            n = SR(Z[1])
            if repr(coeff) != '0':
                if repr(n) == '0':
                    xpow = SR(1)
                elif repr(n) == '1':
                    xpow = x
                else:
                    xpow = x**n
                if ans is None:
                    ans = coeff*xpow
                else:
                    ans += coeff*xpow
        ans._declare_simplified()
        return ans

    def combine(self):
        r"""
        Simplifies \code{self} by combining all terms with the same
        denominator into a single term.

        EXAMPLES:
            sage: var('x, y, a, b, c')
            (x, y, a, b, c)
            sage: f = x*(x-1)/(x^2 - 7) + y^2/(x^2-7) + 1/(x+1) + b/a + c/a
            sage: print f
                                     2
                                    y      (x - 1) x     1     c   b
                                  ------ + --------- + ----- + - + -
                                   2         2         x + 1   a   a
                                  x  - 7    x  - 7
            sage: print f.combine()
                                     2
                                    y  + (x - 1) x     1     c + b
                                    -------------- + ----- + -----
                                         2           x + 1     a
                                        x  - 7
        """
        return self.parent()(self._maxima_().combine())

    def numerator(self):
        """
        EXAMPLES:
            sage: var('a,x,y')
            (a, x, y)
            sage: f = x*(x-a)/((x^2 - y)*(x-a))
            sage: print f
                                                  x
                                                ------
                                                 2
                                                x  - y
            sage: f.numerator()
            x
            sage: f.denominator()
            x^2 - y
        """
        return self.parent()(self._maxima_().num())

    def denominator(self):
        """
        Return the denominator of self.

        EXAMPLES:
            sage: var('x, y, z, theta')
            (x, y, z, theta)
            sage: f = (sqrt(x) + sqrt(y) + sqrt(z))/(x^10 - y^10 - sqrt(theta))
            sage: print f
                                      sqrt(z) + sqrt(y) + sqrt(x)
                                      ---------------------------
                                          10    10
                                       - y   + x   - sqrt(theta)
            sage: f.denominator()
            -y^10 + x^10 - sqrt(theta)
        """
        return self.parent()(self._maxima_().denom())

    def factor_list(self, dontfactor=[]):
        """
        Returns a list of the factors of self, as computed by the
        factor command.

        INPUT:
            self -- a symbolic expression
            dontfactor -- see docs for self.factor.

        REMARK: If you already have a factored expression and just
        want to get at the individual factors, use self._factor_list()
        instead.

        EXAMPLES:
            sage: var('x, y, z')
            (x, y, z)
            sage: f = x^3-y^3
            sage: f.factor()
            (x - y)*(y^2 + x*y + x^2)

        Notice that the -1 factor is separated out:
            sage: f.factor_list()
            [(x - y, 1), (y^2 + x*y + x^2, 1)]

        We factor a fairly straightforward expression:
            sage: factor(-8*y - 4*x + z^2*(2*y + x)).factor_list()
            [(z - 2, 1), (z + 2, 1), (2*y + x, 1)]

        This function also works for quotients:
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: g = f/(36*(1 + 2*y + y^2)); g
            (x^2*y^2 + 2*x*y^2 + y^2 - x^2 - 2*x - 1)/(36*(y^2 + 2*y + 1))
            sage: g.factor(dontfactor=[x])
            (x^2 + 2*x + 1)*(y - 1)/(36*(y + 1))
            sage: g.factor_list(dontfactor=[x])
            [(x^2 + 2*x + 1, 1), (y - 1, 1), (36, -1), (y + 1, -1)]

        An example, where one of the exponents is not an integer.
            sage: var('x, u, v')
            (x, u, v)
            sage: f = expand((2*u*v^2-v^2-4*u^3)^2 * (-u)^3 * (x-sin(x))^3)
            sage: f.factor()
            u^3*(2*u*v^2 - v^2 - 4*u^3)^2*(sin(x) - x)^3
            sage: g = f.factor_list(); g
            [(u, 3), (2*u*v^2 - v^2 - 4*u^3, 2), (sin(x) - x, 3)]

        This example also illustrates that the exponents do not have
        to be integers.
            sage: f = x^(2*sin(x)) * (x-1)^(sqrt(2)*x); f
            (x - 1)^(sqrt(2)*x)*x^(2*sin(x))
            sage: f.factor_list()
            [(x - 1, sqrt(2)*x), (x, 2*sin(x))]
        """
        return self.factor(dontfactor=dontfactor)._factor_list()

    def _factor_list(self):
        r"""
        Turn an expression already in factored form into a
        list of (prime, power) pairs.

        This is used, e.g., internally by the \code{factor_list}
        command.

        EXAMPLES:
            sage: g = factor(x^3 - 1); g
            (x - 1)*(x^2 + x + 1)
            sage: v = g._factor_list(); v
            [(x - 1, 1), (x^2 + x + 1, 1)]
            sage: type(v)
            <type 'list'>
        """
        if isinstance(self, SymbolicArithmetic):
            if self._operator == operator.mul:
                left, right = self._operands
                return left._factor_list() + right._factor_list()
            elif self._operator == operator.pow:
                left, right = self._operands
                return [(left, right)]
            elif self._operator == operator.div:
                left, right = self._operands
                return left._factor_list() + \
                       [(x,-y) for x, y in right._factor_list()]
            elif self._operator == operator.neg:
                expr = self._operands[0]
                v = expr._factor_list()
                return [(SR(-1),SR(1))] + v
        return [(self, 1)]


    ###################################################################
    # solve
    ###################################################################
    def roots(self, x=None, explicit_solutions=True):
        r"""
        Returns roots of \code{self} that can be found exactly, with
        multiplicities.  Not all root are guaranteed to be found.

        WARNING: This is \emph{not} a numerical solver -- use
        \code{find_root} to solve for self == 0 numerically on an
        interval.

        INPUT:
            x -- variable to view the function in terms of
                   (use default variable if not given)
            explicit_solutions -- bool (default True); require that
                roots be explicit rather than implicit
        OUTPUT:
            list of pairs (root, multiplicity)

        If there are infinitely many roots, e.g., a function
        like $\sin(x)$, only one is returned.

        EXAMPLES:
            sage: var('x, a')
            (x, a)

        A simple example:
            sage: ((x^2-1)^2).roots()
            [(-1, 2), (1, 2)]

        A complicated example.
            sage: f = expand((x^2 - 1)^3*(x^2 + 1)*(x-a)); f
            x^9 - a*x^8 - 2*x^7 + 2*a*x^6 + 2*x^3 - 2*a*x^2 - x + a

        The default variable is $a$, since it is the first in alphabetical order:
            sage: f.roots()
            [(x, 1)]

        As a polynomial in $a$, $x$ is indeed a root:
            sage: f.poly(a)
            (x^9 - 2*x^7 + 2*x^3 - x)*1 + (-x^8 + 2*x^6 - 2*x^2 + 1)*a
            sage: f(a=x)
            0

        The roots in terms of $x$ are what we expect:
            sage: f.roots(x)
            [(a, 1), (-1*I, 1), (I, 1), (1, 3), (-1, 3)]

        Only one root of $\sin(x) = 0$ is given:
            sage: f = sin(x)
            sage: f.roots(x)
            [(0, 1)]

        We derive the roots of a general quadratic polynomial:
            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: (a*x^2 + b*x + c).roots(x)
            [((-sqrt(b^2 - 4*a*c) - b)/(2*a), 1), ((sqrt(b^2 - 4*a*c) - b)/(2*a), 1)]


        By default, all the roots are required to be explicit rather than
        implicit.  To get implicit roots, pass \code{explicit_solutions=False}
        to \code{.roots()}
            sage: var('x')
            x
            sage: f = x^(1/9) + (2^(8/9) - 2^(1/9))*(x - 1) - x^(8/9)
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [((x^(8/9) - x^(1/9) + 2^(8/9) - 2^(1/9))/(2^(8/9) - 2^(1/9)), 1)]

        Another example, but involving a degree 5 poly whose roots
        don't get computed explicitly:
            sage: f = x^5 + x^3 + 17*x + 1
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [(x^5 + x^3 + 17*x + 1, 1)]
        """
        if x is None:
            x = self.default_variable()
        S, mul = self.solve(x, multiplicities=True, explicit_solutions=explicit_solutions)
        if len(mul) == 0 and explicit_solutions:
            raise RuntimeError, "no explicit roots found"
        else:
            return [(S[i].rhs(), mul[i]) for i in range(len(mul))]

    def solve(self, x, multiplicities=False, explicit_solutions=False):
        r"""
        Analytically solve the equation \code{self == 0} for the variable $x$.

        WARNING: This is not a numerical solver -- use
        \code{find_root} to solve for self == 0 numerically on an
        interval.

        INPUT:
            x -- variable to solve for
            multiplicities -- bool (default: False); if True, return corresponding multiplicities.
            explicit_solutions -- bool (default:False); if True, require that all solutions
                returned be explicit (rather than implicit)

        EXAMPLES:
            sage: z = var('z')
            sage: (z^5 - 1).solve(z)
            [z == e^(2*I*pi/5), z == e^(4*I*pi/5), z == e^(-(4*I*pi/5)), z == e^(-(2*I*pi/5)), z == 1]
        """
        x = var(x)
        return (self == 0).solve(x, multiplicities=multiplicities, explicit_solutions=explicit_solutions)

    def find_root(self, a, b, var=None, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
        """
        Numerically find a root of self on the closed interval [a,b]
        (or [b,a]) if possible, where self is a function in the one variable.

        INPUT:
            a, b -- endpoints of the interval
            var  -- optional variable
            xtol, rtol -- the routine converges when a root is known
                    to lie within xtol of the value return. Should be
                    >= 0.  The routine modifies this to take into
                    account the relative precision of doubles.
            maxiter -- integer; if convergence is not achieved in
                    maxiter iterations, an error is raised. Must be >= 0.
            full_output -- bool (default: False), if True, also return
                    object that contains information about convergence.

        EXAMPLES:
        Note that in this example both f(-2) and f(3) are positive, yet we still find a
        root in that interval.
            sage: f = x^2 - 1
            sage: f.find_root(-2, 3)
            1.0
            sage: z, result = f.find_root(-2, 3, full_output=True)
            sage: result.converged
            True
            sage: result.flag
            'converged'
            sage: result.function_calls
            11
            sage: result.iterations
            10
            sage: result.root
            1.0

        More examples:
            sage: (sin(x) + exp(x)).find_root(-10, 10)
            -0.588532743981862...

        An example with a square root:
            sage: f = 1 + x + sqrt(x+2); f.find_root(-2,10)
            -1.6180339887498949

       Some examples that Ted Kosan came up with:
            sage: t = var('t')
            sage: v = 0.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))
            sage: v.find_root(0, 0.002)
            0.001540327067911417...

             sage: a = .004*(8*e^(-(300*t)) - 8*e^(-(1200*t)))*(720000*e^(-(300*t)) - 11520000*e^(-(1200*t))) +.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))^2

        There is a 0 very close to the origin:
            sage: show(plot(a, 0, .002),xmin=0, xmax=.002)

        Using solve does not work to find it:
            sage: a.solve(t)
            []

        However \code{find_root} works beautifully:
            sage: a.find_root(0,0.002)
            0.0004110514049349341...

        We illustrate that root finding is only implemented
        in one dimension:
            sage: x, y = var('x,y')
            sage: (x-y).find_root(-2,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.
        """
        a = float(a); b = float(b)
        if var is None:
            w = self.variables()
            if len(w) > 1:
                raise NotImplementedError, "root finding currently only implemented in 1 dimension."
            if len(w) == 0:
                if bool(self == 0):
                    return a
                else:
                    raise RuntimeError, "no zero in the interval, since constant expression is not 0."
            var = repr(w[0])

        f = self._fast_float_(var)
        if a > b:
            a, b = b, a
        return sage.numerical.optimize.find_root(f, a=a,b=b, xtol=xtol, rtol=rtol,
                                                 maxiter=maxiter, full_output=full_output)

    def find_maximum_on_interval(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find the maximum of the expression \code{self} on
        the interval [a,b] (or [b,a]) along with the point at which
        the maximum is attained.

        See the documentation for \code{self.find_minimum_on_interval}
        for more details.

        EXAMPLES:
            sage: f = x*cos(x)
            sage: f.find_maximum_on_interval(0,5)
            (0.5610963381910451, 0.8603335890...)
            sage: f.find_maximum_on_interval(0,5, tol=0.1, maxfun=10)
            (0.561090323458081..., 0.857926501456)
        """
        minval, x = (-self).find_minimum_on_interval(a, b, var=var, tol=tol, maxfun=maxfun)
        return -minval, x

    def find_minimum_on_interval(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find the minimum of the expression \code{self} on
        the interval [a,b] (or [b,a]) and the point at which it
        attains that minimum.  Note that \code{self} must be a
        function of (at most) one variable.

        INPUT:
            var -- variable (default: first variable in self)
            a,b -- endpoints of interval on which to minimize self.
            tol -- the convergence tolerance
            maxfun -- maximum function evaluations

        OUTPUT:

            minval -- (float) the minimum value that self takes on in
                      the interval [a,b]

            x -- (float) the point at which self takes on the minimum value

        EXAMPLES:
            sage: f = x*cos(x)
            sage: f.find_minimum_on_interval(1, 5)
            (-3.2883713955908962, 3.42561846957)
            sage: f.find_minimum_on_interval(1, 5, tol=1e-3)
            (-3.288371361890984, 3.42575079030572)
            sage: f.find_minimum_on_interval(1, 5, tol=1e-2, maxfun=10)
            (-3.2883708459837844, 3.42508402203)
            sage: show(f.plot(0, 20))
            sage: f.find_minimum_on_interval(1, 15)
            (-9.4772942594797929, 9.52933441095)

        ALGORITHM: Uses \module{scipy.optimize.fminbound} which uses Brent's method.

        AUTHOR:
             -- William Stein (2007-12-07)
        """
        a = float(a); b = float(b)
        if var is None:
            var = first_var(self)
        f = lambda w: float(self.substitute({var:float(w)}))
        return sage.numerical.optimize.find_minimum_on_interval(f,
                                                 a=a,b=b,tol=tol, maxfun=maxfun )



    ###################################################################
    # simplify
    ###################################################################
    def simplify(self):
        r"""
        Return the simplified form of this symbolic expression.

        NOTE: Expressions always print simplified; a simplified
        expression is distinguished because the way it prints agrees
        with its underlyilng representation.

        OUTPUT:
            symbolic expression -- a simplified symbolic expression

        EXAMPLES:
        We create the symbolic expression $x - x + 1$, which is of
        course equal to $1$.
            sage: Z = x - x + 1

        It prints as $1$ but is really a symbolic arithmetic object,
        that has the information of the full expression:
            sage: Z
            1
            sage: type(Z)
            <class 'sage.calculus.calculus.SymbolicArithmetic'>

        Calling simplify returns a new object $W$ (it does not change $Z$),
        which actually simplified:
            sage: W = Z.simplify(); W
            1

        Thus $W$ is a single constant:
            sage: type(W)
            <class 'sage.calculus.calculus.SymbolicConstant'>

        The \code{_has_been_simplified} method tells whether an object was
        constructed by simplifying -- or at least is known to be already
        simplified:
            sage: Z._has_been_simplified()
            False
            sage: W._has_been_simplified()
            True

        Note that \code{__init__} automatically calls \code{_simp} when a
        symbolic expression is created:
            sage: f = x*x*x - (1+1+1)*x*x + 2 + 3 + 5 + 7 + 11; f
            x^3 - 3*x^2 + 28
            sage: x*x*x
            x^3

        \code{simplify}, however, can be called manually:
            sage: f
            x^3 - 3*x^2 + 28
            sage: f.simplify()
            x^3 - 3*x^2 + 28
        """
        try:
            if self._simp is None:
                return self
            return self._simp
        except AttributeError:
            S = evaled_symbolic_expression_from_maxima_string(self._maxima_init_())
            S._simp = None
            self._simp = S
            return S

    def simplify_full(self):
        """
        Applies simplify_trig, simplify_rational, and simplify_radical
        to self (in that order).

        ALIAS: simplfy_full and full_simplify are the same.

        EXAMPLES:
            sage: a = log(8)/log(2)
            sage: a.simplify_full()
            3

            sage: f = sin(x)^2 + cos(x)^2
            sage: f.simplify_full()
            1

            sage: f = sin(x/(x^2 + x))
            sage: f.simplify_full()
            sin(1/(x + 1))
        """
        x = self
        x = x.simplify_trig()
        x = x.simplify_rational()
        x = x.simplify_radical()
        return x


    full_simplify = simplify_full

    def simplify_trig(self):
        r"""
        First expands using trig_expand, then employs the identities
        $\sin(x)^2 + \cos(x)^2 = 1$ and $\cosh(x)^2 - \sin(x)^2 = 1$
        to simplify expressions containing tan, sec, etc., to sin,
        cos, sinh, cosh.

        ALIAS: trig_simplify and simplify_trig are the same

        EXAMPLES:
            sage: f = sin(x)^2 + cos(x)^2; f
            sin(x)^2 + cos(x)^2
            sage: f.simplify()
            sin(x)^2 + cos(x)^2
            sage: f.simplify_trig()
            1
        """
        # much better to expand first, since it often doesn't work
        # right otherwise!
        return self.parent()(self._maxima_().trigexpand().trigsimp())

    trig_simplify = simplify_trig

    def simplify_rational(self):
        """
        Simplify by expanding repeatedly rational expressions.

        ALIAS: rational_simplify and simplify_rational are the same

        EXAMPLES:
            sage: f = sin(x/(x^2 + x))
            sage: f
            sin(x/(x^2 + x))
            sage: f.simplify_rational()
            sin(1/(x + 1))

            sage: f = ((x - 1)^(3/2) - (x + 1)*sqrt(x - 1))/sqrt((x - 1)*(x + 1))
            sage: print f
                                          3/2
                                   (x - 1)    - sqrt(x - 1) (x + 1)
                                   --------------------------------
                                        sqrt((x - 1) (x + 1))
            sage: print f.simplify_rational()
                                              2 sqrt(x - 1)
                                            - -------------
                                                    2
                                              sqrt(x  - 1)
        """
        return self.parent()(self._maxima_().fullratsimp())

    rational_simplify = simplify_rational

    # TODO: come up with a way to intelligently wrap Maxima's way of
    # fine-tuning all simplificationsrational

    def simplify_radical(self):
        r"""
        Simplifies this symbolic expression, which can contain logs,
        exponentials, and radicals, by converting it into a form which
        is canonical over a large class of expressions and a given
        ordering of variables

        DETAILS: This uses the Maxima radcan() command. From the
        Maxima documentation: "All functionally equivalent forms are
        mapped into a unique form.  For a somewhat larger class of
        expressions, produces a regular form.  Two equivalent
        expressions in this class do not necessarily have the same
        appearance, but their difference can be simplified by radcan
        to zero.  For some expressions radcan is quite time
        consuming. This is the cost of exploring certain relationships
        among the components of the expression for simplifications
        based on factoring and partial fraction expansions of
        exponents."

        ALIAS: radical_simplify, simplify_radical, simplify_log,
        log_simplify, exp_simplify, simplify_exp are all the same

        EXAMPLES:
            sage: var('x,y,a')
            (x, y, a)

            sage: f = log(x*y)
            sage: f.simplify_radical()
            log(y) + log(x)

            sage: f = (log(x+x^2)-log(x))^a/log(1+x)^(a/2)
            sage: f.simplify_radical()
            log(x + 1)^(a/2)

            sage: f = (e^x-1)/(1+e^(x/2))
            sage: f.simplify_exp()
            e^(x/2) - 1
        """
        maxima.eval('domain: real$')
        res = self.parent()(self._maxima_().radcan())
        maxima.eval('domain: complex$')
        return res

    radical_simplify = simplify_log = log_simplify = simplify_radical
    simplify_exp = exp_simplify = simplify_radical

    ###################################################################
    # factor
    ###################################################################
    def factor(self, dontfactor=[]):
        """
        Factors self, containing any number of variables or functions,
        into factors irreducible over the integers.

        INPUT:
            self -- a symbolic expression
            dontfactor -- list (default: []), a list of variables with
                          respect to which factoring is not to occur.
                          Factoring also will not take place with
                          respect to any variables which are less
                          important (using the variable ordering
                          assumed for CRE form) than those on the
                          `dontfactor' list.

        EXAMPLES:
            sage: var('x, y, z')
            (x, y, z)

            sage: (x^3-y^3).factor()
            (x - y)*(y^2 + x*y + x^2)
            sage: factor(-8*y - 4*x + z^2*(2*y + x))
            (2*y + x)*(z - 2)*(z + 2)
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: F = factor(f/(36*(1 + 2*y + y^2)), dontfactor=[x])
            sage: print F
                                          2
                                        (x  + 2 x + 1) (y - 1)
                                        ----------------------
                                              36 (y + 1)

        If you are factoring a polynomial with rational coefficients
        (and dontfactor is empty) the factorization is done using
        Singular instead of Maxima, so the following is very fast instead
        of dreadfully slow:
            sage: var('x,y')
            (x, y)
            sage: (x^99 + y^99).factor()
            (y + x)*(y^2 - x*y + x^2)*(y^6 - x^3*y^3 + x^6)*...
        """
        if len(dontfactor) > 0:
            m = self._maxima_()
            name = m.name()
            cmd = 'block([dontfactor:%s],factor(%s))'%(dontfactor, name)
            return evaled_symbolic_expression_from_maxima_string(cmd)
        else:
            try:
                f = self.polynomial(QQ)
                w = repr(f.factor())
                return symbolic_expression_from_string(w, _vars)
            except TypeError:
                pass
            return self.parent()(self._maxima_().factor())

    ###################################################################
    # expand
    ###################################################################
    def expand(self):
        r"""
        Expand this symbolic expression.  Products of sums and
        exponentiated sums are multiplied out, numerators of rational
        expressions which are sums are split into their respective
        terms, and multiplications are distributed over addition at
        all levels.

        For polynomials one should usually use \code{expand_rational}
        which uses a more efficient algorithm.

        EXAMPLES:
        We expand the expression $(x-y)^5$ using both method and
        functional notation.
            sage: x,y = var('x,y')
            sage: a = (x-y)^5
            sage: a.expand()
            -y^5 + 5*x*y^4 - 10*x^2*y^3 + 10*x^3*y^2 - 5*x^4*y + x^5
            sage: expand(a)
            -y^5 + 5*x*y^4 - 10*x^2*y^3 + 10*x^3*y^2 - 5*x^4*y + x^5

        Note that \code{expand_rational} may be faster.
            sage: a.expand_rational()
            -y^5 + 5*x*y^4 - 10*x^2*y^3 + 10*x^3*y^2 - 5*x^4*y + x^5

        We expand some other expressions:
            sage: expand((x-1)^3/(y-1))
            x^3/(y - 1) - 3*x^2/(y - 1) + 3*x/(y - 1) - 1/(y - 1)
            sage: expand((x+sin((x+y)^2))^2)
            sin(y^2 + 2*x*y + x^2)^2 + 2*x*sin(y^2 + 2*x*y + x^2) + x^2
        """
        try:
            return self.__expand
        except AttributeError:
            e = self.parent()(self._maxima_().expand())
            self.__expand = e
        return e

    def expand_rational(self):
        """
        Expands self by multiplying out products of sums and
        exponentiated sums, combining fractions over a common
        denominator, cancelling the greatest common divisor of the
        numerator and denominator, then splitting the numerator (if a
        sum) into its respective terms divided by the denominator.

        EXAMPLES:
            sage: x,y = var('x,y')
            sage: a = (x-y)^5
            sage: a.expand_rational()
            -y^5 + 5*x*y^4 - 10*x^2*y^3 + 10*x^3*y^2 - 5*x^4*y + x^5
            sage: a.rational_expand()
            -y^5 + 5*x*y^4 - 10*x^2*y^3 + 10*x^3*y^2 - 5*x^4*y + x^5

        ALIAS: rational_expand and expand_rational are the same
        """
        return self.parent()(self._maxima_().ratexpand())

    rational_expand = expand_rational

    def expand_trig(self, full=False, half_angles=False, plus=True, times=True):
        """
        Expands trigonometric and hyperbolic functions of sums of angles
        and of multiple angles occurring in self.  For best results,
        self should already be expanded.

        INPUT:
            full -- (default: False) To enhance user control of
                    simplification, this function expands only one
                    level at a time by default, expanding sums of
                    angles or multiple angles.  To obtain full
                    expansion into sines and cosines immediately, set
                    the optional parameter full to True.
            half_angles -- (default: False) If True, causes
                    half-angles to be simplified away.
            plus -- (default: True) Controls the sum rule; expansion
                    of sums (e.g. `sin(x + y)') will take place only
                    if plus is True.
            times -- (default: True) Controls the product rule, expansion
                    of products (e.g. sin(2*x)) will take place only if
                    times is True.

        OUTPUT:
            -- a symbolic expression

        EXAMPLES:
            sage: sin(5*x).expand_trig()
            sin(x)^5 - 10*cos(x)^2*sin(x)^3 + 5*cos(x)^4*sin(x)

            sage: cos(2*x + var('y')).trig_expand()
            cos(2*x)*cos(y) - sin(2*x)*sin(y)

        We illustrate various options to this function:
            sage: f = sin(sin(3*cos(2*x))*x)
            sage: f.expand_trig()
            sin(x*(3*cos(cos(2*x))^2*sin(cos(2*x)) - sin(cos(2*x))^3))
            sage: f.expand_trig(full=True)
            sin(x*(3*(sin(cos(x)^2)*cos(sin(x)^2) - cos(cos(x)^2)*sin(sin(x)^2))*(sin(cos(x)^2)*sin(sin(x)^2) + cos(cos(x)^2)*cos(sin(x)^2))^2 - (sin(cos(x)^2)*cos(sin(x)^2) - cos(cos(x)^2)*sin(sin(x)^2))^3))
            sage: sin(2*x).expand_trig(times=False)
            sin(2*x)
            sage: sin(2*x).expand_trig(times=True)
            2*cos(x)*sin(x)
            sage: sin(2 + x).expand_trig(plus=False)
            sin(x + 2)
            sage: sin(2 + x).expand_trig(plus=True)
            cos(2)*sin(x) + sin(2)*cos(x)
            sage: sin(x/2).expand_trig(half_angles=False)
            sin(x/2)
            sage: sin(x/2).expand_trig(half_angles=True)
            sqrt(1 - cos(x))/sqrt(2)

        ALIAS: trig_expand and expand_trig are the same
        """
        M = self._maxima_()
        P = M.parent()
        opt = maxima_options(trigexpand=full, halfangles=half_angles,
                             trigexpandplus=plus, trigexpandtimes=times)
        cmd = 'trigexpand(%s), %s'%(M.name(), opt)
        ans = P(cmd)
        return self.parent()(ans)

    trig_expand = expand_trig


    ###################################################################
    # substitute
    ###################################################################
    def substitute(self, in_dict=None, **kwds):
        r"""
        Takes the symbolic variables given as dict keys or as keywords and
        replaces them with the symbolic expressions given as dict values or as
        keyword values.  Also run when you call a \class{SymbolicExpression}.

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        EXAMPLES:
            sage: x,y,t = var('x,y,t')
            sage: u = 3*y
            sage: u.substitute(y=t)
            3*t

            sage: u = x^3 - 3*y + 4*t
            sage: u.substitute(x=y, y=t)
            y^3 + t

            sage: f = sin(x)^2 + 32*x^(y/2)
            sage: f(x=2, y = 10)
            sin(2)^2 + 1024

            sage: f(x=pi, y=t)
            32*pi^(t/2)

            sage: f = x
            sage: f.variables()
            (x,)

            sage: f = 2*x
            sage: f.variables()
            (x,)

            sage: f = 2*x^2 - sin(x)
            sage: f.variables()
            (x,)
            sage: f(pi)
            2*pi^2
            sage: f(x=pi)
            2*pi^2

            sage: function('f',x)
            f(x)
            sage: (f(x)).substitute(f=log)
            log(x)
            sage: (f(x)).substitute({f:log})
            log(x)
            sage: (x^3 + 1).substitute(x=5)
            126
            sage: (x^3 + 1).substitute({x:5})
            126


        AUTHORS:
            -- Bobby Moretti: Initial version
        """
        X = self.simplify()
        kwds = self.__parse_in_dict(in_dict, kwds)
        kwds = self.__varify_kwds(kwds)
        return X._recursive_sub(kwds)

    def subs(self, *args, **kwds):
        return self.substitute(*args, **kwds)

    def _recursive_sub(self, kwds):
        raise NotImplementedError, "implement _recursive_sub for type '%s'!"%(type(self))

    def _recursive_sub_over_ring(self, kwds, ring):
        raise NotImplementedError, "implement _recursive_sub_over_ring for type '%s'!"%(type(self))

    def __parse_in_dict(self, in_dict, kwds):
        """
        EXAMPLES:
             sage: function('f',x)
             f(x)
             sage: f._SymbolicExpression__parse_in_dict({f:log},{})
             {f: <function log at 0x...>}
             sage: f._SymbolicExpression__parse_in_dict({},{'f':log})
             {'f': <function log at 0x...>}

        """
        if in_dict is None:
            return kwds

        if not isinstance(in_dict, dict):
            if len(kwds) > 0:
                raise ValueError, "you must not both give the variable and specify it explicitly when doing a substitution."
            in_dict = SR(in_dict)
            vars = self.variables()
            #if len(vars) > 1:
            #    raise ValueError, "you must specify the variable when doing a substitution, e.g., f(x=5)"
            if len(vars) == 0:
                return {}
            else:
                return {vars[0]: in_dict}

        # merge dictionaries
        kwds.update(in_dict)
        return kwds

    def __varify_kwds(self, kwds):
        """
        EXAMPLES:
            sage: function('f',x)
            f(x)
            sage: a = f._SymbolicExpression__parse_in_dict({f:log},{})
            sage: f._SymbolicExpression__varify_kwds(a)
            {f: <function log at 0x...>}
            sage: b = f._SymbolicExpression__parse_in_dict({},{'f':log})
            sage: f._SymbolicExpression__varify_kwds(b)
            {f: <function log at 0x...>}

        """
        return dict([(var(k) if isinstance(k,str) else k, w) for k,w in kwds.iteritems()])

    def substitute_over_ring(self, in_dict=None, ring=None, **kwds):
        X = self.simplify()
        kwds = self.__parse_in_dict(in_dict, kwds)
        kwds = self.__varify_kwds(kwds)
        if ring is None:
            return X._recursive_sub(kwds)
        else:
            return X._recursive_sub_over_ring(kwds, ring)


    ###################################################################
    # Expression substitution
    ###################################################################
    def subs_expr(self, *equations):
        """
        Given a dictionary of key:value pairs, substitute all occurences
        of key for value in self.

        WARNING: This is a formal pattern substitution, which may or
        may not have any mathematical meaning.  The exact rules used
        at present in Sage are determined by Maxima's subst command.
        Sometimes patterns are not replaced even though one would think
        they should be -- see examples below.

        EXAMPLES:
            sage: f = x^2 + 1
            sage: f.subs_expr(x^2 == x)
            x + 1

            sage: var('x,y,z'); f = x^3 + y^2 + z
            (x, y, z)
            sage: f.subs_expr(x^3 == y^2, z == 1)
            2*y^2 + 1

            sage: f = x^2 + x^4
            sage: f.subs_expr(x^2 == x)
            x^4 + x
            sage: f = cos(x^2) + sin(x^2)
            sage: f.subs_expr(x^2 == x)
            sin(x) + cos(x)

            sage: f(x,y,t) = cos(x) + sin(y) + x^2 + y^2 + t
            sage: f.subs_expr(y^2 == t)
            (x, y, t) |--> sin(y) + cos(x) + x^2 + 2*t

        The following seems really weird, but it *is* what maple does:
            sage: f.subs_expr(x^2 + y^2 == t)
            (x, y, t) |--> sin(y) + y^2 + cos(x) + x^2 + t
            sage: maple.eval('subs(x^2 + y^2 = t, cos(x) + sin(y) + x^2 + y^2 + t)')          # optional requires maple
            'cos(x)+sin(y)+x^2+y^2+t'
            sage: maxima.quit()
            sage: maxima.eval('cos(x) + sin(y) + x^2 + y^2 + t, x^2 + y^2 = t')
            'sin(y)+y^2+cos(x)+x^2+t'

        Actually Mathematica does something that makes more sense:
            sage: mathematica.eval('Cos[x] + Sin[y] + x^2 + y^2 + t /. x^2 + y^2 -> t')       # optional -- requires mathematica
            2 t + Cos[x] + Sin[y]
        """
        for x in equations:
            if not isinstance(x, SymbolicEquation):
                raise TypeError, "each expression must be an equation"
        R = self.parent()
        v = ','.join(['%s=%s'%(x.lhs()._maxima_init_(), x.rhs()._maxima_init_()) \
                      for x in equations])
        return R(self._maxima_().subst(v))

    ###################################################################
    # Real and imaginary parts
    ###################################################################
    def real(self):
        r"""
        Return the real part of \code{self}.

        EXAMPLES:
            sage: a = log(3+4*I)
            sage: print a
                                             log(4  I + 3)
            sage: print a.real()
                                                log(5)
            sage: print a.imag()
                                                     4
                                              arctan(-)
                                                     3

        Now make a and b symbolic and compute the general real part:
            sage: var('a,b')
            (a, b)
            sage: f = log(a + b*I)
            sage: f.real()
            log(b^2 + a^2)/2
        """
        return self.parent()(self._maxima_().real())

    def imag(self):
        r"""
        Return the imaginary part of \code{self}.

        EXAMPLES:
            sage: sqrt(-2).imag()
            sqrt(2)

        We simplify $\ln(\exp(z))$ to $z$ for $-\pi<{\rm Im}(z)<=\pi$:

            sage: z = var('z')
            sage: f = log(exp(z))
            sage: assume(-pi < imag(z))
            sage: assume(imag(z) <= pi)
            sage: print f
        			       z
            sage: forget()

        A more symbolic example:
            sage: var('a, b')
            (a, b)
            sage: f = log(a + b*I)
            sage: f.imag()
            arctan(b/a)
        """
        return self.parent()(self._maxima_().imag())

    def conjugate(self):
        r"""
        The complex conjugate of \code{self}.

        EXAMPLES:
            sage: a = 1 + 2*I
            sage: a.conjugate()
            1 - 2*I
            sage: a = sqrt(2) + 3^(1/3)*I; a
            3^(1/3)*I + sqrt(2)
            sage: a.conjugate()
            sqrt(2) - 3^(1/3)*I
        """
        return self.parent()(self._maxima_().conjugate())

    def norm(self):
        r"""
        The complex norm of \code{self}, i.e., \code{self} times its complex conjugate.

        EXAMPLES:
            sage: a = 1 + 2*I
            sage: a.norm()
            5
            sage: a = sqrt(2) + 3^(1/3)*I; a
            3^(1/3)*I + sqrt(2)
            sage: a.norm()
            3^(2/3) + 2
            sage: CDF(a).norm()
            4.08008382305
            sage: CDF(a.norm())
            4.08008382305
        """
        m = self._maxima_()
        return self.parent()((m*m.conjugate()).expand())

    ###################################################################
    # Partial fractions
    ###################################################################
    def partial_fraction(self, var=None):
        r"""
        Return the partial fraction expansion of \code{self} with respect to
        the given variable.

        INPUT:
            var -- variable name or string (default: first variable)

        OUTPUT:
            Symbolic expression

        EXAMPLES:
            sage: var('x')
            x
            sage: f = x^2/(x+1)^3
            sage: f.partial_fraction()
            1/(x + 1) - 2/(x + 1)^2 + 1/(x + 1)^3
            sage: print f.partial_fraction()
                                        1        2          1
                                      ----- - -------- + --------
                                      x + 1          2          3
                                              (x + 1)    (x + 1)

        Notice that the first variable in the expression is used by default:
            sage: var('y')
            y
            sage: f = y^2/(y+1)^3
            sage: f.partial_fraction()
            1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3

            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f.partial_fraction()
            y^2/(y^3 + 3*y^2 + 3*y + 1) + 1/(x - 1)^2 + 1/(x - 1)^3

        You can explicitly specify which variable is used.
            sage: f.partial_fraction(y)
            1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3 + x/(x^3 - 3*x^2 + 3*x - 1)
        """
        if var is None:
            var = self.default_variable()
        return self.parent()(self._maxima_().partfrac(var))

    ###################################################################
    # Fast Evaluation
    ###################################################################

    def _fast_float_(self, *vars):
        """
        EXAMPLES:
            sage: x,y,z = var('x,y,z')
            sage: f = 1 + sin(x)/x + sqrt(z^2+y^2)/cosh(x)
            sage: ff = f._fast_float_('x', 'y', 'z')
            sage: f(1.0,2.0,3.0)
            4.1780638977866...
            sage: ff(1.0,2.0,3.0)
            4.17806389778660...
        """
        try:
            return fast_float.fast_float_constant(float(self))
        except:
            raise NotImplementedError # return lambda x: float(self(x))


class Symbolic_object(SymbolicExpression):
    r"""
    A class representing a symbolic expression in terms of a \class{SageObject} (not
    necessarily a `constant').
    """
    def __init__(self, obj):
        SymbolicExpression.__init__(self)
        self._obj = obj

    def __hash__(self):
        return hash(self._obj)

    #def derivative(self, *args):
        # TODO: remove
    #    return self.parent().zero_element()

    def obj(self):
        """
        EXAMPLES:
        """
        return self._obj

    def _fast_float_(self, *vars):
        return fast_float.fast_float_constant(float(self))

    def __float__(self):
        """
        EXAMPLES:
        """
        return float(self._obj)

    def __complex__(self):
        """
        EXAMPLES:
        """
        return complex(self._obj)

    def _mpfr_(self, field):
        """
        EXAMPLES:
        """
        return field(self._obj)

    def _complex_mpfr_field_(self, C):
        """
        EXAMPLES:
        """
        return C(self._obj)

    def _complex_double_(self, C):
        """
        EXAMPLES:
        """
        return C(self._obj)

    def _real_double_(self, R):
        """
        EXAMPLES:
        """
        return R(self._obj)

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
        """
        return R(self._obj)

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
        """
        return repr(self._obj)

    def _latex_(self):
        """
        EXAMPLES:
        """
        return latex(self._obj)

    def str(self, bits=None):
        if bits is None:
            return str(self._obj)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self._obj))

    def _maxima_init_(self):
        return maxima_init(self._obj)

    def _sys_init_(self, system):
        return sys_init(self._obj, system)

    def number_of_arguments(self):
        """
        Returns the number of arguments this object can take.

        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(e)
            sage: a.number_of_arguments()
            0
        """
        return 0

def maxima_init(x):
    try:
        return x._maxima_init_()
    except AttributeError:
        return repr(x)

def sys_init(x, system):
    try:
        return x.__getattribute__('_%s_init_'%system)()
    except AttributeError:
        try:
            return x._system_init_(system)
        except AttributeError:
            return repr(x)

class SymbolicConstant(Symbolic_object):
    def __init__(self, x):
        from sage.rings.rational import Rational
        if isinstance(x, Rational):
            if x.is_integral() and x >= 0:
                self._precedence = 10**6
            else:
                self._precedence = 2000
        Symbolic_object.__init__(self, x)
##         if hasattr(x, 'parent'):
##             Symbolic_object.__init__(self, x)
##         elif isinstance(x, int):
##             Symbolic_object.__init__(self, Integer(x))
##         elif isinstance(x, float):
##             Symbolic_object.__init__(self, RealNumber(x))
##         else:
##             raise ValueError, "%s needs parent"%type(x)

    #def _is_atomic(self):
    #    try:
    #        return self._atomic
    #    except AttributeError:
    #        if isinstance(self, Rational):
    #            self._atomic = False
    #        else:
    #            self._atomic = True
    #        return self._atomic
    def _is_atomic(self):
        try:
            return self._atomic
        except AttributeError:
            try:
                return self._obj._is_atomic()
            except AttributeError:
                if isinstance(self._obj, int):
                    return True

    def _fast_float_(self, *vars):
        return fast_float.fast_float_constant(float(self))

    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
            sage: a = SR(5/6)
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: a(x=3)
            5/6
        """
        return self

    def _recursive_sub_over_ring(self, kwds, ring):
        return ring(self)

    def _algebraic_(self, field):
        """
        EXAMPLES:
            sage: a = SR(5/6)
            sage: AA(a)
            5/6
            sage: type(AA(a))
            <class 'sage.rings.qqbar.AlgebraicReal'>
            sage: QQbar(a)
            5/6
            sage: type(QQbar(a))
            <class 'sage.rings.qqbar.AlgebraicNumber'>
            sage: from sage.calculus.calculus import SymbolicConstant
            sage: i = SymbolicConstant(I)
            sage: AA(i)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
            sage: QQbar(i)
            1*I
            sage: phi = SymbolicConstant(golden_ratio)
            sage: AA(phi)
            1.618033988749895?
            sage: QQbar(phi)
            1.618033988749895?
        """

        # Out of the many kinds of things that can be in a SymbolicConstant,
        # we accept only rational numbers (or things that can be coerced
        # to rational) and a few instances of Constant.

        if isinstance(self._obj, sage.functions.constants.Constant) \
               and hasattr(self._obj, '_algebraic_'):
            return self._obj._algebraic_(field)
        return field(Rational(self._obj))

    def _add_(self, right):
        """
        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(2)
            sage: b = a+2; b
            4
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: b = sum([a for i in range(1000)]); b
            2000
            sage: type(_)
            <class 'sage.calculus.calculus.SymbolicConstant'>
        """
        if isinstance(right, SymbolicConstant):
            try:
                self_parent = self._obj.parent()
                right_parent = right._obj.parent()
                if self_parent != SR and right_parent != SR and ( self_parent.has_coerce_map_from(right_parent) or right_parent.has_coerce_map_from(self_parent) ):
                    return SymbolicConstant( operator.add(self._obj, right._obj) )

            except AttributeError:
                #Either self._obj or right._obj doesn't have a
                #parent method (like 'int' or 'float')
                constants = isinstance(self._obj, sage.functions.constants.Constant) or isinstance(right._obj, sage.functions.constants.Constant)
                if not constants:
                    return SymbolicConstant( operator.add(self._obj, right._obj) )



        return SymbolicArithmetic([self, right], operator.add)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(2)
            sage: b = a-2; b
            0
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: b = SR(2000)
            sage: for i in range(1000): b -= a;
            sage: b
            0
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>

        """
        if isinstance(right, SymbolicConstant):
            try:
                self_parent = self._obj.parent()
                right_parent = right._obj.parent()
                if self_parent != SR and right_parent != SR and ( self_parent.has_coerce_map_from(right_parent) or right_parent.has_coerce_map_from(self_parent) ):
                    return SymbolicConstant( operator.sub(self._obj, right._obj) )

            except AttributeError:
                #Either self._obj or right._obj doesn't have a
                #parent method (like 'int' or 'float')
                constants = isinstance(self._obj, sage.functions.constants.Constant) or isinstance(right._obj, sage.functions.constants.Constant)
                if not constants:
                    return SymbolicConstant( operator.sub(self._obj, right._obj) )



        return SymbolicArithmetic([self, right], operator.sub)



    def _mul_(self, right):
        """
        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(2)
            sage: b = a*2; b
            4
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: prod([a for i in range(1000)])
            10715086071862673209484250490600018105614048117055336074437503883703510511249361224931983788156958581275946729175531468251871452856923140435984577574698574803934567774824230985421074605062371141877954182153046474983581941267398767559165543946077062914571196477686542167660429831652624386837205668069376
            sage: type(_)
            <class 'sage.calculus.calculus.SymbolicConstant'>
        """
        if isinstance(right, SymbolicConstant):
            try:
                self_parent = self._obj.parent()
                right_parent = right._obj.parent()
                if self_parent != SR and right_parent != SR and ( self_parent.has_coerce_map_from(right_parent) or right_parent.has_coerce_map_from(self_parent) ):
                    return SymbolicConstant( operator.mul(self._obj, right._obj) )

            except AttributeError:
                #Either self._obj or right._obj doesn't have a
                #parent method (like 'int' or 'float')
                constants = isinstance(self._obj, sage.functions.constants.Constant) or isinstance(right._obj, sage.functions.constants.Constant)
                if not constants:
                    return SymbolicConstant( operator.mul(self._obj, right._obj) )

        return SymbolicArithmetic([self, right], operator.mul)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(2)
            sage: b = a/2; b
            1
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: b = SR(2^1000)
            sage: for i in range(1000): b /= a;
            sage: b
            1
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>

        """
        if isinstance(right, SymbolicConstant):
            try:
                self_parent = self._obj.parent()
                right_parent = right._obj.parent()
                if self_parent != SR and right_parent != SR and ( self_parent.has_coerce_map_from(right_parent) or right_parent.has_coerce_map_from(self_parent) ):
                    return SymbolicConstant( operator.div(self._obj, right._obj) )

            except AttributeError:
                #Either self._obj or right._obj doesn't have a
                #parent method (like 'int' or 'float')
                constants = isinstance(self._obj, sage.functions.constants.Constant) or isinstance(right._obj, sage.functions.constants.Constant)
                if not constants:
                    return SymbolicConstant( operator.div(self._obj, right._obj) )



        return SymbolicArithmetic([self, right], operator.div)

    def __neg__(self):
        return SymbolicConstant(-self._obj)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: SR = SymbolicExpressionRing()
            sage: a = SR(2)
            sage: b = a^2; b
            4
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicArithmetic'>
        """
        right = self.parent()(right)
        return SymbolicArithmetic([self, right], operator.pow)


class SymbolicPolynomial(Symbolic_object):
    """
    An element of a polynomial ring as a formal symbolic expression.

    EXAMPLES:
    A single variate polynomial:
        sage: R.<x> = QQ[]
        sage: f = SR(x^3 + x)
        sage: f(y=7)
        x^3 + x
        sage: f(x=5)
        130
        sage: f.integral(x)
        x^4/4 + x^2/2
        sage: f(x=var('y'))
        y^3 + y

    A multivariate polynomial:

        sage: R.<x,y,theta> = ZZ[]
        sage: f = SR(x^3 + x + y + theta^2); f
        x^3 + theta^2 + x + y
        sage: f(x=y, theta=y)
        y^3 + y^2 + 2*y
        sage: f(x=5)
        y + theta^2 + 130

    The polynomial must be over a field of characteristic 0.
        sage: R.<w> = GF(7)[]
        sage: f = SR(w^3 + 1)
        Traceback (most recent call last):
        ...
        TypeError: polynomial must be over a field of characteristic 0.
    """
    # for now we do nothing except pass the info on to the superconstructor. It's
    # not clear to me why we need anything else in this class -Bobby
    def __init__(self, p):
        if p.parent().base_ring().characteristic() != 0:
            raise TypeError, "polynomial must be over a field of characteristic 0."
        Symbolic_object.__init__(self, p)

    def _recursive_sub(self, kwds):
        return self._recursive_sub_over_ring(kwds, ring)

    def _recursive_sub_over_ring(self, kwds, ring):
        f = self._obj
        if is_Polynomial(f):
            # Single variable case
            v = f.parent().variable_name()
            if kwds.has_key(v):
                return ring(f(kwds[v]))
            else:
                if not ring is SR:
                    return ring(self)
        else:
            # Multivariable case
            t = []
            for g in f.parent().gens():
                s = repr(g)
                if kwds.has_key(s):
                    t.append(kwds[s])
                else:
                    t.append(g)
            return ring(f(*t))

    def variables(self):
        r"""
        Return sorted list of variables that occur in \code{self}.
        The ordering is alphabetic.

        EXAMPLES:
            sage: R.<x> = QQ[]; S.<y> = R[]
            sage: f = x+y*x+y^2
            sage: g = SR(f)
            sage: g.variables()
            (x, y)
        """
        P = self._obj.parent()
        variables = map(var, P.variable_names_recursive())
        variables.sort()
        return tuple(variables)

    def number_of_arguments(self):
        r"""
        Returns the number of arguments this object can take.  For
        \class{SymbolicPolynomial}s, this is just the number of variables
        of the polynomial.

        EXAMPLES:
            sage: R.<x> = QQ[]; S.<y> = R[]
            sage: f = x+y*x+y^2
            sage: g = SR(f)
            sage: g.number_of_arguments()
            2
        """
        return len(self.variables())

    def polynomial(self, base_ring):
        r"""
        Return \code{self} as a polynomial over the given base ring, if possible.

        INPUT:
           base_ring -- a ring

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = SR(x^2 -2/3*x + 1)
            sage: f.polynomial(QQ)
            x^2 - 2/3*x + 1
            sage: f.polynomial(GF(19))
            x^2 + 12*x + 1
        """
        f = self._obj
        if base_ring is f.base_ring():
            return f
        return f.change_ring(base_ring)

    def _fast_float_(self, *vars):
        # use Horners rule
        return self._obj._fast_float_(*vars)


##################################################################


zero_constant = SymbolicConstant(Integer(0))

class SymbolicOperation(SymbolicExpression):
    r"""
    A parent class representing any operation on \class{SymbolicExpression} objects.
    """
    def __init__(self, operands):
        SymbolicExpression.__init__(self)
        self._operands = operands   # don't even make a copy -- ok, since immutable.

    def variables(self):
        r"""
        Return sorted list of variables that occur in the simplified
        form of \code{self}.  The ordering is alphabetic.

        EXAMPLES:
            sage: var('x,y,z,w,a,b,c')
            (x, y, z, w, a, b, c)
            sage: f = (x - x) + y^2 - z/z + (w^2-1)/(w+1); f
            y^2 + (w^2 - 1)/(w + 1) - 1
            sage: f.variables()
            (w, y)

            sage: (x + y + z + a + b + c).variables()
            (a, b, c, x, y, z)

            sage: (x^2 + x).variables()
            (x,)
        """
        if not self._has_been_simplified():
            return self.simplify().variables()

        try:
            return self.__variables
        except AttributeError:
            pass
        vars = list(set(sum([list(op.variables()) for op in self._operands], [])))

        vars.sort(var_cmp)
        vars = tuple(vars)
        self.__variables = vars
        return vars

    def number_of_arguments(self):
        r"""
        Returns the number of arguments this object can take.

        EXAMPLES:
            sage: x,y,z = var('x,y,z')
            sage: (x+y).number_of_arguments()
            2
            sage: (x+1).number_of_arguments()
            1
            sage: (sin+1).number_of_arguments()
            1
            sage: (sin+x).number_of_arguments()
            1
            sage: (sin+x+y).number_of_arguments()
            2
            sage: (sin(z)+x+y).number_of_arguments()
            3
            sage: (sin+cos).number_of_arguments()
            1
            sage: (sin(x+y)).number_of_arguments()
            2

            sage: ( 2^(8/9) - 2^(1/9) )(x-1)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0

        Note that \code{self} is simplified first:
            sage: f = x + pi - x; f
            pi
            sage: f.number_of_arguments()
            0
        """
        try:
            return self.__number_of_args
        except AttributeError:
            pass
        variables = self.variables()
        if not self._has_been_simplified():
            n = self.simplify().number_of_arguments()
        else:
            # We need to do this maximum to correctly handle the case where
            # self is something like (sin+1)
            n = max( max(map(lambda i: i.number_of_arguments(), self._operands)+[0]), len(variables) )
        self.__number_of_args = n
        return n

def var_cmp(x,y):
    return cmp(repr(x), repr(y))

symbols = {operator.add:' + ', operator.sub:' - ', operator.mul:'*',
        operator.div:'/', operator.pow:'^', operator.neg:'-'}


class SymbolicArithmetic(SymbolicOperation):
    r"""
    Represents the result of an arithmetic operation on
    $f$ and $g$.
    """
    def __init__(self, operands, op):
        SymbolicOperation.__init__(self, operands)
        self._operator = op
        # assume a really low precedence by default
        self._precedence = -1
        # set up associativity and precedence rules
        if op is operator.neg:
            self._binary = False
            self._unary = True
            self._precedence = 2000
        else:
            self._binary = True
            self._unary = False
        if op is operator.pow:
            self._precedence = 3000
            self._l_assoc = False
            self._r_assoc = True
        elif op is operator.mul:
            self._precedence = 2000
            self._l_assoc = True
            self._r_assoc = True
        elif op is operator.div:
            self._precedence = 2000
            self._l_assoc = True
            self._r_assoc = False
        elif op is operator.sub:
            self._precedence = 1000
            self._l_assoc = True
            self._r_assoc = False
        elif op is operator.add:
            self._precedence = 1000
            self._l_assoc = True
            self._r_assoc = True

    def __call__(self, *args, **kwargs):
        """
        Method for handling a function call.

        EXAMPLES:
            sage: x,y,z=var('x,y,z')

            sage: h = sin + cos
            sage: h(1)
            sin(1) + cos(1)
            sage: h(x)
            sin(x) + cos(x)
            sage: h = 3*sin
            sage: h(1)
            3*sin(1)
            sage: h(x)
            3*sin(x)

            sage: (sin+cos)(1)
            sin(1) + cos(1)
            sage: (sin+1)(1)
            sin(1) + 1
            sage: (x+sin)(5)
            sin(5) + 5
            sage: (y+sin)(5)
            sin(5) + 5
            sage: (x+y+sin)(5)
            y + sin(5) + 5


            sage: f = x + 2*y + 3*z
            sage: f(1)
            3*z + 2*y + 1
            sage: f(0,1)
            3*z + 2
            sage: f(0,0,1)
            3

            sage: (sqrt(2) + 17)(x+2)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0
            sage: (I*17+3*5)(x+2)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0
        """
        if kwargs and args:
            raise ValueError, "args and kwargs cannot both be specified"

        if len(args) == 1 and isinstance(args[0], dict):
            kwargs = dict(map(lambda x: (repr(x[0]), x[1]), args[0].iteritems()))

        if kwargs:
            #Handle the case where kwargs are specified
            x = var('x')

            new_ops = []
            for op in self._operands:
                try:
                    new_ops.append( op(**kwargs) )
                except ValueError:
                    new_ops.append(op)



        else:
            #Handle the case where args are specified

            #Get all the variables
            variables = list( self.arguments() )

            if len(args) > self.number_of_arguments():
                raise ValueError, "the number of arguments must be less than or equal to %s"%self.number_of_arguments()

            new_ops = []
            for op in self._operands:
                try:
                    op_vars = op.variables()
                    if len(op_vars) == 0:
                        if len(args) != 0:
                            new_ops.append( op(args[0]) )
                        else:
                            new_ops.append( op )
                        continue
                    else:
                        indices = filter(lambda i: i < len(args), map(variables.index, op_vars))
                        if len(indices) == 0:
                            new_ops.append( op )
                        else:
                            new_ops.append( op(*[args[i] for i in indices]) )
                except ValueError:
                    new_ops.append( op )

        #Check to see if all of the new_ops are symbolic constants
        #If so, then we should return a symbolic constant.
        new_ops = map(SR, new_ops)
        is_constant = all(map(lambda x: isinstance(x, SymbolicConstant), new_ops))
        if is_constant:
            return SymbolicConstant( self._operator(*map(lambda x: x._obj, new_ops)) )
        else:
            return self._operator(*new_ops)


    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
            sage: var('x, y, z, w')
            (x, y, z, w)
            sage: f = (x - x) + y^2 - z/z + (w^2-1)/(w+1); f
            y^2 + (w^2 - 1)/(w + 1) - 1
            sage: f(y=10)
            (w^2 - 1)/(w + 1) + 99
            sage: f(w=1,y=10)
            99
            sage: f(y=w,w=y)
            (y^2 - 1)/(y + 1) + w^2 - 1

            sage: f = y^5 - sqrt(2)
            sage: f(10)
            100000 - sqrt(2)

            sage: a = x^2; b = a(2); b
            4
            sage: type(b)
            <class 'sage.calculus.calculus.SymbolicConstant'>
        """
        ops = self._operands
        new_ops = [SR(op._recursive_sub(kwds)) for op in ops]

        #Check to see if all of the new_ops are symbolic constants
        #If so, then we should return a symbolic constant.
        is_constant = all(map(lambda x: isinstance(x, SymbolicConstant), new_ops))
        if is_constant:
            return SymbolicConstant( self._operator(*map(lambda x: x._obj, new_ops)) )
        else:
            return self._operator(*new_ops)

    def _recursive_sub_over_ring(self, kwds, ring):
        ops = self._operands
        if self._operator == operator.pow:
            new_ops = [ops[0]._recursive_sub_over_ring(kwds, ring=ring), Integer(ops[1])]
        else:
            new_ops = [op._recursive_sub_over_ring(kwds, ring=ring) for op in ops]
        return ring(self._operator(*new_ops))

    def _fast_float_(self, *vars):
        """
        EXAMPLES:
            sage: x,y = var('x,y')
            sage: f = x*x-y
            sage: ff = f._fast_float_('x','y')
            sage: ff(2,3)
            1.0
        """
        fops = [op._fast_float_(*vars) for op in self._operands]
        return self._operator(*fops)

    def _convert(self, typ):
        """
        Convert self to the given type by converting each of the
        operands to that type and doing the arithmetic.

        EXAMPLES:
            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: f._convert(RDF)
            -1.40006081534
            sage: f._convert(float)
            -1.4000608153399503

        Converting to an int can have surprising consequences, since
        Python int is ``floor'' and one individual factor can floor to 0
        but the product doesn't:
            sage: int(f)
            -1
            sage: f._convert(int)
            0
            sage: int(sqrt(2))
            1
            sage: int(cos(3))
            0

        TESTS:
        This illustrates how the conversion works even when a type
        exception is raised, since here one operand is still x (in the
        unsimplified form):
            sage: f = sin(0)*x
            sage: f._convert(CDF)
            0
        """
        try:
            fops = [typ(op) for op in self._operands]
        except TypeError:
            g = self.simplify()
            if self is g:
                raise
            else:
                return typ(g)
        return self._operator(*fops)


    def __float__(self):
        """
        TESTS:
            sage: f=x*sin(0)
            sage: float(f(x=1))
            0.0
            sage: w = I - I
            sage: float(w)
            0.0
        """
        return self._convert(float)

    def __complex__(self):
        """
        EXAMPLES:
            sage: complex((-2)^(1/4))
            (0.840896415253714...+0.840896415253714...j)

        TESTS:
            sage: complex(I - I)
            0j
            sage: w = I-I; complex(w)
            0j
            sage: complex(w * x)
            0j
        """
        return self._convert(complex)

    def _mpfr_(self, field):
        """
        EXAMPLES:
            sage: RealField(200)(sqrt(2))
            1.4142135623730950488016887242096980785696718753769480731767

        TESTS:
            sage: w = I-I; RR(w)
            0.000000000000000
            sage: w = I-I; RealField(200)(w)
            0.00000000000000000000000000000000000000000000000000000000000
            sage: RealField(200)(x*sin(0))
            0.00000000000000000000000000000000000000000000000000000000000
        """
        return self._convert(field)

    def _complex_mpfr_field_(self, field):
        """
        EXAMPLES:
            sage: CC(sqrt(2))
            1.41421356237310
            sage: a = sqrt(-2); a
            sqrt(2)*I
            sage: CC(a)
            1.41421356237310*I
            sage: ComplexField(200)(a)
            1.4142135623730950488016887242096980785696718753769480731767*I
            sage: ComplexField(100)((-1)^(1/10))
            0.95105651629515357211643933338 + 0.30901699437494742410229341718*I

        TESTS:
            sage: CC(x*sin(0))
            0
        """
        return self._convert(field)

    def _complex_double_(self, field):
        """
        EXAMPLES:
            sage: CDF((-1)^(1/3))
            0.5 + 0.866025403784*I

        Watch out -- right now Maxima algebraically simplifies the above to -1:
            sage: (-1)^(1/3)
            (-1)^(1/3)

        So when doing this conversion it is the non-simplified form
        that is converted, for efficiency purposes, which can result
        in a different answer in some cases.  You can always force
        using the simplified form:

            sage: a = (-1)^(1/3)
            sage: CDF(a.simplify())
            0.5 + 0.866025403784*I

        TESTS:
            sage: CDF(x*sin(0))
            0
        """
        return self._convert(field)

    def _real_double_(self, field):
        """
        EXAMPLES:
            sage: RDF(sqrt(2))
            1.41421356237


        TESTS:
            sage: RDF(x*sin(0))
            0.0
        """
        return self._convert(field)

    def _real_rqdf_(self, field):
        """
        EXAMPLES:
            sage: RQDF(sqrt(2))
            1.414213562373095048801688724209698078569671875376948073176679738

        TESTS:
            sage: RQDF(x*sin(0))
            0.000000000000000000000000000000000000000000000000000000000000000
        """
        return self._convert(field)

    def _algebraic_(self, field):
        """
        Convert a symbolic expression to an algebraic number.

        EXAMPLES:
            sage: QQbar(sqrt(2) + sqrt(8))
            4.242640687119285?
            sage: AA(sqrt(2) ^ 4) == 4
            True
            sage: AA(-golden_ratio)
            -1.618033988749895?
            sage: QQbar((2*I)^(1/2))
            1.0000000000000000? + 1.0000000000000000?*I

        TESTS:
            sage: AA(x*sin(0))
            0
            sage: QQbar(x*sin(0))
            0
        """

        # We try to avoid simplifying, because maxima's simplify command
        # can change the value of a radical expression (by changing which
        # root is selected).
        try:
            if self._operator is operator.pow:
                base = field(self._operands[0])
                expt = Rational(self._operands[1])
                return field(field(base) ** expt)
            else:
                return self._convert(field)
        except TypeError:
            if self._has_been_simplified():
                raise
            else:
                return self.simplify()._algebraic_(field)

    def _is_atomic(self):
        try:
            return self._atomic
        except AttributeError:
            op = self._operator
            # todo: optimize with a dictionary
            if op is operator.add:
                self._atomic = False
            elif op is operator.sub:
                self._atomic = False
            elif op is operator.mul:
                self._atomic = True
            elif op is operator.div:
                self._atomic = False
            elif op is operator.pow:
                self._atomic = True
            elif op is operator.neg:
                self._atomic = False
            return self._atomic

    def _repr_(self, simplify=True):
        """
        TESTS:
            sage: var('r')
            r
            sage: a = (1-1/r)^(-1); a
            1/(1 - 1/r)
            sage: a.derivative(r)
            -1/((1 - 1/r)^2*r^2)

            sage: var('a,b')
            (a, b)
            sage: s = 0*(1/a) + -b*(1/a)*(1 + -1*0*(1/a))*(1/(a*b + -1*b*(1/a)))
            sage: s
            -b/(a*(a*b - b/a))
            sage: s(a=2,b=3)
            -1/3
            sage: -3/(2*(2*3-(3/2)))
            -1/3
            sage: (-1)^(1/4)
            (-1)^(1/4)

            sage: (-(x-1)/2)._latex_(simplify=False)
            '\\frac{-\\left( x - 1 \\right)}{2}'

        """
        if simplify:
            if hasattr(self, '_simp'):
                if self._simp is None:
                    return self._repr_(simplify=False)
                return self._simp._repr_(simplify=False)
            else:
                return self.simplify()._repr_(simplify=False)

        ops = self._operands
        op = self._operator
        s = [x._repr_(simplify=simplify) for x in ops]

        # if an operand is a rational number, trick SAGE into thinking it's an
        # operation
        li = []
        for o in ops:
            try:
                obj = o._obj
                # negative numbers are not handled correctly because _is_atomic has no sense of precedence
                if o is ops[0] and str(obj)[0] == '-':
                    temp = SymbolicConstant(obj)
                    temp._operator = operator.neg
                    temp._binary = False
                    temp._unary = True
                    temp._precedence = 2000
                    li.append(temp)
                elif isinstance(obj, Rational):
                    temp = SymbolicConstant(obj)
                    if not temp._obj.is_integral():
                        temp._operator = operator.div
                        temp._l_assoc = True
                        temp._r_assoc = False
                        temp._precedence = 2000
                        temp._binary = True
                        temp._unary = False
                    li.append(temp)
                else:
                    li.append(o)
            except AttributeError:
                li.append(o)

        ops = li

        rop = ops[0]
        if self._binary:
            lop = rop
            rop = ops[1]

        lparens = True
        rparens = True

        if self._binary:
            try:
                l_operator = lop._operator
            except AttributeError:
                # if it's not arithmetic on the left, see if it's atomic
                try:
                    prec = lop._precedence
                except AttributeError:
                    if lop._is_atomic():
                    # if it has no concept of precedence, leave the parens
                        lparens = False
                else:
                    # if it a higher precedence, don't draw parens
                    if self._precedence < lop._precedence:
                        lparens = False
            else:
                # if the left op is the same is this operator
                if op is l_operator:
                    # if it's left associative, get rid of the left parens
                    if self._l_assoc:
                        lparens = False
                # different operators, same precedence, get rid of the left parens
                elif self._precedence == lop._precedence:
                    if self._l_assoc:
                        lparens = False
                # if we have a lower precedence than the left, get rid of the parens
                elif self._precedence < lop._precedence:
                    lparens = False

        try:
            r_operator = rop._operator
        except AttributeError:
            try:
                prec = rop._precedence
            except AttributeError:
                if rop._is_atomic():
                    rparens = False
            else:
                if self._precedence < rop._precedence:
                    rparens = False
        else:
            if rop._binary:
                if op is r_operator:
                    try:
                        if self._r_assoc:
                            rparens = False
                    except AttributeError:
                        pass
                elif self._precedence == rop._precedence:
                    try:
                        if self._r_assoc:
                            rparens = False
                    except AttributeError:
                        pass
                # if the RHS has higher precedence, it comes first and parens are
                # redundant
                elif self._precedence < rop._precedence:
                    rparens = False
        if self._binary:
            if lparens:
                s[0] = '(%s)'% s[0]
            if rparens:
                s[1] = '(%s)'% s[1]

            return '%s%s%s' % (s[0], symbols[op], s[1])

        elif self._unary:
            if rparens:
                s[0] = '(%s)'%s[0]
            return '%s%s' % (symbols[op], s[0])

    def _latex_(self, simplify=True):
        """
        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: f=(x+y)*(x-y)*(x^2-2)*(y^2-3)
            sage: latex(f)
            {{{\left( {x}^{2}  - 2 \right) \left( x - y \right)} \left( y + x \right)} \left( {y}^{2}  - 3 \right)}
            sage: latex(cos*(x+1))
            {\left( x + 1 \right) \cos}
        """
        # if we are not simplified, return the latex of a simplified version
        if simplify and not self._has_been_simplified():
            return self.simplify()._latex_()
        op = self._operator
        ops = self._operands
        s = []
        for x in self._operands:
            try:
                s.append(x._latex_(simplify=simplify))
            except TypeError:
                s.append(x._latex_())
        ## what it used to be --> s = [x._latex_() for x in self._operands]

        if op is operator.add:
            return '%s + %s' % (s[0], s[1])
        elif op is operator.sub:
            return '%s - %s' % (s[0], s[1])
        elif op is operator.mul:
            for i in [0,1]:
                if isinstance(ops[i], SymbolicArithmetic) and ops[i]._operator in [operator.add, operator.sub]:
                    s[i] = r'\left( %s \right)' %s[i]
            return '{%s %s}' % (s[0], s[1])
            # used to be --> return '{%s \\cdot %s}' % (s[0], s[1])
        elif op is operator.div:
            return '\\frac{%s}{%s}' % (s[0], s[1])
        elif op is operator.pow:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub) \
               or ops[0]._has_op(operator.mul) or ops[0]._has_op(operator.div) \
               or ops[0]._has_op(operator.pow):
                s[0] = r'\left( %s \right)' % s[0]
            return '{%s}^{%s} ' % (s[0], s[1])
        elif op is operator.neg:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub):
                s[0] = r'\left( %s \right)'%s[0]
            return '-%s' % s[0]

    def _maxima_init_(self):
        ops = self._operands
        if self._operator is operator.neg:
            return '-(%s)' % ops[0]._maxima_init_()
        else:
            return '(%s) %s (%s)' % (ops[0]._maxima_init_(),
                             infixops[self._operator],
                             ops[1]._maxima_init_())

    def _sympy_(self):
        """Converts any expression to SymPy."""

        # Current implementation is a little fragile - it first converts the
        # expression to string, then preparses it, then gets rid of "Integer"
        # and then sympifies this string.

        # In order to make this more robust, one would have to implement
        # _sympy_ recursively in all expressions, but that would pollute Sage
        # quite a lot, and currently this oneliner does the job and it's easy
        # to understand.
        from sympy import sympify
        from sage.all import preparse
        return sympify(preparse(repr(self)).replace("Integer",""))

    def _sys_init_(self, system):
        ops = self._operands
        if self._operator is operator.neg:
            return '-(%s)' % sys_init(ops[0], system)
        else:
            return '(%s) %s (%s)' % (sys_init(ops[0], system),
                             infixops[self._operator],
                             sys_init(ops[1], system))

is_python_identifier = re.compile("[_a-zA-Z][_a-zA-Z0-9]*")

def is_SymbolicVariable(x):
    """
    Return True if $x$ is a symbolic variable.

    INPUT:
        x -- object
    OUTPUT:
        bool -- True precisely if x is a symbolic variable.

    EXAMPLES:
        sage: is_SymbolicVariable('x')
        False
        sage: is_SymbolicVariable(x)
        True
    """
    return isinstance(x, SymbolicVariable)

class SymbolicVariable(SymbolicExpression):
    """
    A symbolic variable, which is a calculus object.

    EXAMPLES:
        sage: z = var('z')
        sage: type(z)
        <class 'sage.calculus.calculus.SymbolicVariable'>
    """
    def __init__(self, name):
        SymbolicExpression.__init__(self)
        self._name = name
        if len(name) == 0:
            raise ValueError, "variable name must be nonempty"
        elif not is_python_identifier.match(name):
            raise ValueError, "variable name is not a valid Python identifier"

    def __hash__(self):
        """
        Return the hash of this symbolic variable, which is just the
        hash of the underlying name (as a string).

        EXAMPLES:
            sage: z = var('z')
            sage: hash(z) #random due to architecture dependence
            -1563822213
            sage: hash('z') #random due to architecture dependence
            -1563822213
        """
        return hash(self._name)

    def _fast_float_(self, *vars):
        r"""
        Returns a quickly-evaluating function with named parameters
        \code{vars}. Specifically, if \code{self} is the $n$-th parameter
        it returns a function extracting the $n$-th item out of a tuple.

        EXAMPLES:
            sage: f = x._fast_float_('x', 'y')
            sage: f(1,2)
            1.0
            sage: f = x._fast_float_('y', 'x')
            sage: f(1,2)
            2.0
            sage: sqrt(2)._fast_float_()(2)
            1.4142135623730951
        """
        if self._name in vars:
            return fast_float.fast_float_arg(list(vars).index(self._name))
        svars = [repr(x) for x in vars]
        if self._name in svars:
            return fast_float.fast_float_arg(list(svars).index(self._name))
        try:
            return fast_float.fast_float_constant(float(self))
        except TypeError:
            raise ValueError, "free variable: %s" % self._name

    def _recursive_sub(self, kwds):
        # do the replacement if needed
        if kwds.has_key(self):
            return kwds[self]
        else:
            return self

    def _recursive_sub_over_ring(self, kwds, ring):
        if kwds.has_key(self):
            return ring(kwds[self])
        else:
            return ring(self)

    def variables(self):
        r"""
        Return sorted list of variables that occur in the simplified
        form of \code{self}.
        """
        return (self, )

    def number_of_arguments(self):
        r"""
        Returns the number of arguments of \code{self}.

        EXAMPLES:
            sage: x = var('x')
            sage: x.number_of_arguments()
            1
        """
        return 1

    def __cmp__(self, right):
        if isinstance(right, SymbolicVariable):
            return cmp(self._repr_(), right._repr_())
        else:
            return SymbolicExpression.__cmp__(self, right)

    def _repr_(self, simplify=True):
        return self._name

    def __str__(self):
        return self._name

    def _latex_(self):
        try:
            return self.__latex
        except AttributeError:
            pass
        self.__latex = latex_variable_name(self._name)
        return self.__latex

    def _maxima_init_(self):
        return self._name

    def _sys_init_(self, system):
        return self._name


_vars = {}
def var(s, create=True):
    r"""
    Create a symbolic variable with the name \emph{s}.

    INPUTS:
        s -- a string, either a single variable name
             or a space or comma separated list of
             variable names

    NOTE: sage.calculus.calculus.var is better suited for using var in library
    code since it won't touch the global namespace. To create a new variable
    in the global namespace, use sage.calculus.var.var. That is, use
    sage.calculus.calculus.var when defining a symbolic variable for use in
    sage.calculus.calculus functions or library code.

    EXAMPLES:
    Note that sage.calculus.calculus.var defines a variable which is locally
    defined in the calculus namespace but is not globally defined:

        sage: alpha = 42; alpha
        42
        sage: sage.calculus.calculus.var('alpha')
        alpha
        sage: alpha
        42

    The variable is still of type SymbolicVariable and belongs to the
    symbolic expression ring:

        sage: type(alpha)
        <type 'sage.rings.integer.Integer'>
        sage: type(sage.calculus.calculus.var('alpha'))
        <class 'sage.calculus.calculus.SymbolicVariable'>
        sage: var('beta')
        beta
        sage: type(beta)
        <class 'sage.calculus.calculus.SymbolicVariable'>

    TESTS:
        sage: var('xx')
        xx
        sage: var('.foo')
        Traceback (most recent call last):
        ...
        ValueError: variable name is not a valid Python identifier
        sage: var('.foo/x')
        Traceback (most recent call last):
        ...
        ValueError: variable name is not a valid Python identifier
    """
    if isinstance(s, SymbolicVariable):
        return s
    s = str(s)
    if ',' in s:
        return tuple([var(x.strip()) for x in s.split(',')])
    elif ' ' in s:
        return tuple([var(x.strip()) for x in s.split()])
    try:
        v = _vars[s]
        _syms[s] = v
        return v
    except KeyError:
        if not create:
            raise ValueError, "the variable '%s' has not been defined"%var
        pass
    v = SymbolicVariable(s)
    _vars[s] = v
    _syms[s] = v
    return v


#########################################################################################
#  Callable functions
#########################################################################################
def is_CallableSymbolicExpressionRing(x):
    """
    Return True if x is a callable symbolic expression.

    INPUT:
        x -- object

    OUTPUT:
        bool

    EXAMPLES:
        sage: is_CallableSymbolicExpressionRing(QQ)
        False
        sage: var('x,y,z')
        (x, y, z)
        sage: is_CallableSymbolicExpressionRing(CallableSymbolicExpressionRing((x,y,z)))
        True
    """
    return isinstance(x, CallableSymbolicExpressionRing_class)

class CallableSymbolicExpressionRing_class(CommutativeRing):
    def __init__(self, args):
        self._args = args
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        """

        TESTS:
            sage: f(x) = x+1; g(y) = y+1
            sage: f.parent()(g)
            x |--> y + 1
            sage: g.parent()(f)
            y |--> x + 1
            sage: f(x) = x+2*y; g(y) = y+3*x
            sage: f.parent()(g)
            x |--> y + 3*x
            sage: g.parent()(f)
            y |--> 2*y + x
        """
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, SymbolicExpression):
            if isinstance(x, PrimitiveFunction):
                return CallableSymbolicExpression(self, x( *self._args ))
            if isinstance(x, CallableSymbolicExpression):
                x = x._expr
            return CallableSymbolicExpression(self, x)
        return self._coerce_try(x, [SR])

    def _repr_(self):
        """
        String representation of ring of callable symbolic expressions.

        EXAMPLES:
            sage: R = CallableSymbolicExpressionRing(var('x,y,theta'))
            sage: R._repr_()
            'Callable function ring with arguments (x, y, theta)'
        """
        return "Callable function ring with arguments %s"%(self._args,)

    def args(self):
        r"""
        Returns the arguments of \code{self}.  The order that the variables appear
        in \code{self.args()} is the order that is used in evaluating the elements
        of \code{self}.

        EXAMPLES:
            sage: x,y = var('x,y')
            sage: f(x,y) = 2*x+y
            sage: f.parent().args()
            (x, y)
            sage: f(y,x) = 2*x+y
            sage: f.parent().args()
            (y, x)
        """
        return self._args

    arguments = args

    def zero_element(self):
        """
        Return the zero element of the ring of callable symbolic expressions.

        EXAMPLES:
            sage: R = CallableSymbolicExpressionRing(var('x,y,theta'))
            sage: f = R.zero_element(); f
            (x, y, theta) |--> 0
            sage: f(2,3,4)
            0
        """
        try:
            return self.__zero_element
        except AttributeError:
            z = CallableSymbolicExpression(self, SR.zero_element())
            self.__zero_element = z
            return z

    def _an_element_impl(self):
        """
        Return an element of the ring of callabel symbolic expressions.
        This is used by the coercion model.

        EXAMPLES:
            sage: R = CallableSymbolicExpressionRing(var('x,y,theta'))
            sage: R._an_element_impl()
            (x, y, theta) |--> 0
        """
        return CallableSymbolicExpression(self, SR._an_element())


_cfr_cache = {}
def CallableSymbolicExpressionRing(args, check=True):
    if check:
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = args[0]
        for arg in args:
            if not isinstance(arg, SymbolicVariable):
                raise TypeError, "Must construct a function with a tuple (or list) of" \
                                +" SymbolicVariables."
        args = tuple(args)
    if _cfr_cache.has_key(args):
        R = _cfr_cache[args]()
        if not R is None:
            return R
    R = CallableSymbolicExpressionRing_class(args)
    _cfr_cache[args] = weakref.ref(R)
    return R

def is_CallableSymbolicExpression(x):
    r"""
    Returns true if \var{x} is a callable symbolic expression.

    EXAMPLES:
        sage: var('a x y z')
        (a, x, y, z)
        sage: f(x,y) = a + 2*x + 3*y + z
        sage: is_CallableSymbolicExpression(f)
        True
        sage: is_CallableSymbolicExpression(a+2*x)
        False
        sage: def foo(n): return n^2
        ...
        sage: is_CallableSymbolicExpression(foo)
        False
    """
    return isinstance(x, CallableSymbolicExpression)

class CallableSymbolicExpression(SymbolicExpression):
    r"""
    A callable symbolic expression that knows the ordered list of
    variables on which it depends.

    EXAMPLES:
        sage: var('a, x, y, z')
        (a, x,   y, z)
        sage: f(x,y) = a + 2*x + 3*y + z
        sage: f
        (x, y) |--> z + 3*y + 2*x + a
        sage: f(1,2)
        z + a + 8
    """
    def __init__(self, parent, expr):
        RingElement.__init__(self, parent)
        self._expr = expr

    def variables(self):
        """
        EXAMPLES:
            sage: a = var('a')
            sage: g(x) = sin(x) + a
            sage: g.variables()
            (a, x)
            sage: g.args()
            (x,)
            sage: g(y,x,z) = sin(x) + a - a
            sage: g
            (y, x, z) |--> sin(x)
            sage: g.args()
            (y, x, z)
        """
        return self._expr.variables()

    def integral(self, x=None, a=None, b=None):
        r"""
        Returns an integral of \code{self}.
        """
        if a is None:
            return SymbolicExpression.integral(self, x, None, None)
            # if l. endpoint is None, compute an indefinite integral
        else:
            if x is None:
                x = self.default_variable()
            if not isinstance(x, SymbolicVariable):
                x = var(repr(x))
                # if we supplied an endpoint, then we want to return a number.
            return SR(self._maxima_().integrate(x, a, b))

    integrate = integral

    def expression(self):
        """
        Return the underlying symbolic expression (i.e., forget the
        extra map structure).
        """
        return self._expr

    def args(self):
        return self.parent().args()

    def arguments(self):
        r"""
        Returns the arguments of \code{self}.  The order that the
        variables appear in \code{self.arguments()} is the order that
        is used in \code{self.__call__}.

        EXAMPLES:
            sage: x,y = var('x,y')
            sage: f(x,y) = 2*x+y
            sage: f.arguments()
            (x, y)
            sage: f(2)
            y + 4
            sage: f(2, 1)
            5

            sage: f(y,x) = 2*x+y
            sage: f.arguments()
            (y, x)
            sage: f(2)
            2*x + 2
            sage: f(2, 1)
            4
        """
        return self.args()

    def number_of_arguments(self):
        r"""
        Returns the number of arguments of \code{self}.

        EXAMPLES:
            sage: a = var('a')
            sage: g(x) = sin(x) + a
            sage: g.number_of_arguments()
            1
            sage: g(x,y,z) = sin(x) - a + a
            sage: g.number_of_arguments()
            3
        """
        return len(self.args())

    def _maxima_init_(self):
        return self._expr._maxima_init_()

    def _fast_float_(self, *vars):
        return self._expr._fast_float_(*vars)

    def __float__(self):
        return float(self._expr)

    def __complex__(self):
        return complex(self._expr)

    def _mpfr_(self, field):
        """
        Coerce to a multiprecision real number.

        EXAMPLES:
            sage: RealField(100)(SR(10))
            10.000000000000000000000000000
        """
        return (self._expr)._mpfr_(field)

    def _complex_mpfr_field_(self, field):
        return field(self._expr)

    def _complex_double_(self, C):
        return C(self._expr)

    def _real_double_(self, R):
        return R(self._expr)

    def _real_rqdf_(self, R):
        return R(self._expr)

    # TODO: should len(args) == len(vars)?
    def __call__(self, *args):
        vars = self.args()
        dct = dict( (vars[i], args[i]) for i in range(len(args)) )
        return self._expr.substitute(dct)

    def _repr_(self, simplify=True):
        args = self.args()
        if len(args) == 1:
            return "%s |--> %s" % (args[0], self._expr._repr_(simplify=simplify))
        else:
            args = ", ".join(map(str, args))
            return "(%s) |--> %s" % (args, self._expr._repr_(simplify=simplify))

    def _latex_(self):
        """
        Finds the LaTeX representation of this expression.

        EXAMPLES:
            sage: f(A, t, omega, psi) = A*cos(omega*t - psi)
            sage: f._latex_()
            '\\left(A, t, \\omega, \\psi \\right)\\ {\\mapsto}\\ {\\cos \\left( {\\omega t} - \\psi \\right) A}'

            sage: f(mu) =  mu^3
            sage: f._latex_()
            '\\mu \\ {\\mapsto}\\ {\\mu}^{3} '
        """
        args = self.args()
        args = [arg._latex_() for arg in args]
        if len(args) == 1:
            return "%s \\ {\mapsto}\\ %s" % (args[0],
                    self._expr._latex_())
        else:
            vars = ", ".join(args)
            # the weird TeX is to workaround an apparent JsMath bug
            return "\\left(%s \\right)\\ {\\mapsto}\\ %s" % (vars, self._expr._latex_())

    def _neg_(self):
        return CallableSymbolicExpression(self.parent(), -self._expr)

    def __add__(self, right):
        """
        EXAMPLES:
            sage: var('x y z n m')
            (x, y, z, n, m)
            sage: f(x,n,y) = x^n + y^m;  g(x,n,m,z) = x^n +z^m
            sage: f + g
            (x, n, m, y, z) |--> z^m + y^m + 2*x^n
            sage: g + f
            (x, n, m, y, z) |--> z^m + y^m + 2*x^n

            sage: f(x) = x^2
            sage: f+sin
            x |--> sin(x) + x^2
            sage: g(y) = y^2
            sage: g+sin
            y |--> sin(y) + y^2
            sage: h = g+sin
            sage: h(2)
            sin(2) + 4
        """
        if isinstance(right, CallableSymbolicExpression):
            if self.parent() is right.parent():
                return self._add_(right)
            else:
                args = self._unify_args(right)
                R = CallableSymbolicExpressionRing(args)
                return R(self) + R(right)
        else:
            return RingElement.__add__(self, right)

    def _add_(self, right):
        return CallableSymbolicExpression(self.parent(), self._expr + right._expr)

    def __sub__(self, right):
        """
        EXAMPLES:
            sage: var('x y z n m')
            (x, y, z, n, m)
            sage: f(x,n,y) = x^n + y^m;  g(x,n,m,z) = x^n +z^m
            sage: f - g
            (x, n, m, y, z) |--> y^m - z^m
            sage: g - f
            (x, n, m, y, z) |--> z^m - y^m
        """
        if isinstance(right, CallableSymbolicExpression):
            if self.parent() is right.parent():
                return self._sub_(right)
            else:
                args = self._unify_args(right)
                R = CallableSymbolicExpressionRing(args)
                return R(self) - R(right)
        else:
            return RingElement.__sub__(self, right)

    def _sub_(self, right):
        return CallableSymbolicExpression(self.parent(), self._expr - right._expr)

    def __mul__(self, right):
        """
        EXAMPLES:
            sage: var('x y z a b c n m')
            (x, y, z, a, b, c, n, m)

            sage: f(x) = x+2*y; g(y) = y+3*x
            sage: f*(2/3)
            x |--> 2*(2*y + x)/3
            sage: f*g
            (x, y) |--> (y + 3*x)*(2*y + x)
            sage: (2/3)*f
            x |--> 2*(2*y + x)/3

            sage: f(x,y,z,a,b) = x+y+z-a-b; f
            (x, y, z, a, b) |--> z + y + x - b - a
            sage: f * (b*c)
            (x, y, z, a, b) |--> b*c*(z + y + x - b - a)
            sage: g(x,y,w,t) = x*y*w*t
            sage: f*g
            (x, y, a, b, t, w, z) |--> t*w*x*y*(z + y + x - b - a)
            sage: (f*g)(2,3)
            6*t*w*(z - b - a + 5)

            sage: f(x,n,y) = x^n + y^m;  g(x,n,m,z) = x^n +z^m
            sage: f * g
            (x, n, m, y, z) |--> (y^m + x^n)*(z^m + x^n)
            sage: g * f
            (x, n, m, y, z) |--> (y^m + x^n)*(z^m + x^n)
        """
        if isinstance(right, CallableSymbolicExpression):
            if self.parent() is right.parent():
                return self._mul_(right)
            else:
                args = self._unify_args(right)
                R = CallableSymbolicExpressionRing(args)
                return R(self)*R(right)
        else:
            return RingElement.__mul__(self, right)

    def _mul_(self, right):
        return CallableSymbolicExpression(self.parent(), self._expr * right._expr)

    def __div__(self, right):
        """
        EXAMPLES:
            sage: var('x,y,z,m,n')
            (x, y, z, m, n)
            sage: f(x,n,y) = x^n + y^m;  g(x,n,m,z) = x^n +z^m
            sage: f / g
            (x, n, m, y, z) |--> (y^m + x^n)/(z^m + x^n)
            sage: g / f
            (x, n, m, y, z) |--> (z^m + x^n)/(y^m + x^n)
        """
        if isinstance(right, CallableSymbolicExpression):
            if self.parent() is right.parent():
                return self._div_(right)
            else:
                args = self._unify_args(right)
                R = CallableSymbolicExpressionRing(args)
                return R(self) / R(right)
        else:
            return RingElement.__div__(self, right)

    def _div_(self, right):
        return CallableSymbolicExpression(self.parent(), self._expr / right._expr)

    def __pow__(self, right):
        return CallableSymbolicExpression(self.parent(), self._expr ** right)

    def _unify_args(self, x):
        r"""
        Takes the variable list from another \class{CallableSymbolicExpression} object and
        compares it with the current \class{CallableSymbolicExpression} object's variable list,
        combining them according to the following rules:

        Let \code{a} be \code{self}'s variable list, let \code{b} be
        \code{y}'s variable list.

        \begin{enumerate}

        \item If \code{a == b}, then the variable lists are identical,
              so return that variable list.

        \item If \code{a} $\neq$ \code{b}, then check if the first $n$
              items in \code{a} are the first $n$ items in \code{b},
              or vice-versa. If so, return a list with these $n$
              items, followed by the remaining items in \code{a} and
              \code{b} sorted together in alphabetical order.

        \end{enumerate}

        Note: When used for arithmetic between \class{CallableSymbolicExpression}s,
        these rules ensure that the set of \class{CallableSymbolicExpression}s will have
        certain properties.  In particular, it ensures that the set is
        a \emph{commutative} ring, i.e., the order of the input
        variables is the same no matter in which order arithmetic is
        done.

        INPUT:
            x -- A CallableSymbolicExpression

        OUTPUT:
            A tuple of variables.

        EXAMPLES:
            sage: f(x, y, z) = sin(x+y+z)
            sage: f
            (x, y, z) |--> sin(z + y + x)
            sage: g(x, y) = y + 2*x
            sage: g
            (x, y) |--> y + 2*x
            sage: f._unify_args(g)
            (x, y, z)
            sage: g._unify_args(f)
            (x, y, z)

            sage: f(x, y, z) = sin(x+y+z)
            sage: g(w, t) = cos(w - t)
            sage: g
            (w, t) |--> cos(w - t)
            sage: f._unify_args(g)
            (t, w, x, y, z)

            sage: f(x, y, t) = y*(x^2-t)
            sage: f
            (x, y, t) |--> (x^2 - t)*y
            sage: g(x, y, w) = x + y - cos(w)
            sage: f._unify_args(g)
            (x, y, t, w)
            sage: g._unify_args(f)
            (x, y, t, w)
            sage: f*g
            (x, y, t, w) |--> (x^2 - t)*y*(y + x - cos(w))

            sage: f(x,y, t) = x+y
            sage: g(x, y, w) = w + t
            sage: f._unify_args(g)
            (x, y, t, w)
            sage: g._unify_args(f)
            (x, y, t, w)
            sage: f + g
            (x, y, t, w) |--> y + x + w + t

        AUTHORS:
            -- Bobby Moretti, thanks to William Stein for the rules

        """
        a = self.args()
        b = x.args()

        # Rule #1
        if [str(x) for x in a] == [str(x) for x in b]:
            return a

        # Rule #2
        new_list = []
        done = False
        i = 0
        while not done and i < min(len(a), len(b)):
            if var_cmp(a[i], b[i]) == 0:
                new_list.append(a[i])
                i += 1
            else:
                done = True

        temp = set([])
        # Sorting remaining variables.
        for j in range(i, len(a)):
            if not a[j] in temp:
                temp.add(a[j])

        for j in range(i, len(b)):
            if not b[j] in temp:
                temp.add(b[j])

        temp = list(temp)
        temp.sort(var_cmp)
        new_list.extend(temp)
        return tuple(new_list)


#########################################################################################
#  End callable functions
#########################################################################################

class SymbolicComposition(SymbolicOperation):
    r"""
    Represents the symbolic composition of $f \circ g$.
    """
    def __init__(self, f, g):
        """
        INPUT:
            f, g -- both must be in the symbolic expression ring.
        """
        SymbolicOperation.__init__(self, [f,g])

    def _polynomial_(self, R):
        """
        Symbolic compositions cannot be converted to polynomials unless
        they are constants.

        EXAMPLES:
            sage: sqrt(2).polynomial(RR)
            1.41421356237310

            sage: sqrt(2).polynomial(CC)
            1.41421356237310

            sage: cos(x).polynomial(QQ)
            Traceback (most recent call last):
            ....
            TypeError: cannot convert self (= cos(x)) to a polynomial

            sage: sqrt(x).polynomial(QQ)
            Traceback (most recent call last):
            ....
            TypeError: cannot convert self (= sqrt(x)) to a polynomial

            sage: K3.<a> = NumberField(sqrt(x))
            Traceback (most recent call last):
            ....
            TypeError: polynomial (=sqrt(x)) must be a polynomial.
        """
        if self.number_of_arguments() == 0:
            #Convert self into R's base ring and then into R since
            #self must be a constant.
            return R( R.base_ring()(self) )
        else:
            raise TypeError, "cannot convert self (= %s) to a polynomial"%str(self).strip()


    def number_of_arguments(self):
        r"""
        Returns the number of arguments that \code{self} can take.

        EXAMPLES:
            sage: sqrt(x).number_of_arguments()
            1
            sage: sqrt(2).number_of_arguments()
            0
        """
        try:
            return self.__number_of_args
        except AttributeError:
            pass
        variables = self.variables()
        if not self._has_been_simplified():
            n = self.simplify().number_of_arguments()
        else:
            # Note that we use self._operands[1:] so we don't include the
            # number of arguments that the function takes since it is
            # already being "called"
            n = max( max(map(lambda i: i.number_of_arguments(), self._operands[1:])+[0]), len(variables) )
        self.__number_of_args = n
        return n

    def _recursive_sub(self, kwds):
        ops = self._operands
        return ops[0](ops[1]._recursive_sub(kwds))

    def _recursive_sub_over_ring(self, kwds, ring):
        ops = self._operands
        return ring(ops[0](ops[1]._recursive_sub_over_ring(kwds, ring=ring)))

    def _is_atomic(self):
        return True

    def _repr_(self, simplify=True):
        if simplify:
            if hasattr(self, '_simp'):
                if self._simp is None:
                    return self._repr_(simplify=False)
                return self._simp._repr_(simplify=False)
            else:
                return self.simplify()._repr_(simplify=False)
        ops = self._operands
        try:
            return ops[0]._repr_evaled_(ops[1]._repr_(simplify=False))
        except AttributeError:
            return "%s(%s)"% (ops[0]._repr_(simplify=False), ops[1]._repr_(simplify=False))

    def _maxima_init_(self):
        ops = self._operands
        try:
            return ops[0]._maxima_init_evaled_(ops[1]._maxima_init_())
        except AttributeError:
            return '%s(%s)' % (ops[0]._maxima_init_(), ops[1]._maxima_init_())

    def _latex_(self):
        if not self._has_been_simplified():
            return self.simplify()._latex_()
        ops = self._operands

        #Check to see if the function has a _latex_composition method
        if hasattr(ops[0], '_latex_composition'):
            return ops[0]._latex_composition(ops[1])

        # certain functions (such as \sqrt) need braces in LaTeX
        if (ops[0]).tex_needs_braces():
            return r"%s{ %s }" % ( (ops[0])._latex_(), (ops[1])._latex_())
        # ... while others (such as \cos) don't
        return r"%s \left( %s \right)"%((ops[0])._latex_(),(ops[1])._latex_())

    def _sys_init_(self, system):
        ops = self._operands
        return '%s(%s)' % (sys_init(ops[0],system), sys_init(ops[1],system))

    def _mathematica_init_(self):
        system = 'mathematica'
        ops = self._operands
        return '%s[%s]' % (sys_init(ops[0],system).capitalize(), sys_init(ops[1],system))

    def __float__(self):
        f = self._operands[0]
        g = self._operands[1]
        return float(f._approx_(float(g)))

    def _fast_float_(self, *vars):
        f = self._operands[0]
        g = self._operands[1]._fast_float_(*vars)
        try:
            return f(g)
        except TypeError:
            if f is abs_symbolic:
                return abs(g) # special case
            else:
                return fast_float.fast_float_func(f, g)
        return lambda x: f(g(x))

    def __complex__(self):
        """
        Convert this symbolic composition to a Python complex number.

        EXAMPLES:
            sage: complex(cos(3))
            (-0.98999249660044542+0j)
            sage: complex(cos(3*I))
            (10.067661995777771+0j)
        """
        f = self._operands[0]
        g = self._operands[1]
        return complex(f._complex_approx_(complex(g)))

    def _mpfr_(self, field):
        """
        Coerce to a multiprecision real number.

        EXAMPLES:
            sage: RealField(100)(sin(2)+cos(2))
            0.49315059027853930839845163641

            sage: RR(sin(pi))
            1.22464679914735e-16

            sage: type(RR(sqrt(163)*pi))
            <type 'sage.rings.real_mpfr.RealNumber'>

            sage: RR(coth(pi))
            1.00374187319732
            sage: RealField(100)(coth(pi))
            1.0037418731973212882015526912
            sage: RealField(200)(arccos(1/10))
            1.4706289056333368228857985121870581235299087274579233690964
        """
        f = self._operands[0]
        g = self._operands[1]
        x = f(g._mpfr_(field))
        if isinstance(x, SymbolicExpression):
            if field.prec() <= 53:
                return field(float(x))
            else:
                raise TypeError, "precision loss"
        else:
            return x


    def _complex_mpfr_field_(self, field):
        """
        Coerce to a multiprecision complex number.

        EXAMPLES:
            sage: ComplexField(100)(sin(2)+cos(2)+I)
            0.49315059027853930839845163641 + 1.0000000000000000000000000000*I

        """
        f = self._operands[0]
        g = self._operands[1]
        x = f(g._complex_mpfr_field_(field))
        if isinstance(x, SymbolicExpression):
            if field.prec() <= 53:
                return field(complex(x))
            else:
                raise TypeError, "precision loss"
        else:
            return x

    def _complex_double_(self, field):
        """
        Coerce to a complex double.

        EXAMPLES:
            sage: CDF(sin(2)+cos(2)+I)
            0.493150590279 + 1.0*I
            sage: CDF(coth(pi))
            1.0037418732
        """
        f = self._operands[0]
        g = self._operands[1]
        z = f(g._complex_double_(field))
        if isinstance(z, SymbolicExpression):
            return field(complex(z))
        return z

    def _real_double_(self, field):
        """
        Coerce to a real double.

        EXAMPLES:
            sage: RDF(sin(2)+cos(2))
            0.493150590279
        """
        if not self._has_been_simplified():
            return self.simplify()._real_double_(field)
        f = self._operands[0]
        g = self._operands[1]
        z = f(g._real_double_(field))
        if isinstance(z, SymbolicExpression):
            return field(float(z))
        return z

    def _real_rqdf_(self, field):
        """
        Coerce to a real qdrf.

        EXAMPLES:

        """
        if not self._has_been_simplified():
            return self.simplify()._real_rqdf_(field)
        f = self._operands[0]
        g = self._operands[1]
        z = f(g._real_rqdf_(field))
        if isinstance(z, SymbolicExpression):
            raise TypeError, "precision loss"
        else:
            return z

    def _algebraic_(self, field):
        """
        Coerce to an algebraic number.

        EXAMPLES:
            sage: QQbar(sqrt(2))
            1.414213562373095?
            sage: AA(abs(1+I))
            1.414213562373095?
        """
        # We try to avoid simplifying, because maxima's simplify command
        # can change the value of a radical expression (by changing which
        # root is selected).
        f = self._operands[0]
        g = self._operands[1]
        try:
            return field(f(g._algebraic_(field)))
        except (TypeError, ValueError):
            if self._has_been_simplified():
                raise
            else:
                return self.simplify()._algebraic_(field)


class PrimitiveFunction(SymbolicExpression):
    def __init__(self, needs_braces=False):
        SymbolicExpression.__init__(self)
        self._tex_needs_braces = needs_braces

    def _recursive_sub(self, kwds):
        if kwds.has_key(self):
            return kwds[self]
        return self

    def _recursive_sub_over_ring(self, kwds, ring):
        if kwds.has_key(self):
            return kwds[self]
        return self

    def plot(self, *args, **kwds):
        f = self(var('x'))
        return SymbolicExpression.plot(f, *args, **kwds)

    def _is_atomic(self):
        return True

    def tex_needs_braces(self):
        return self._tex_needs_braces

    def __call__(self, x, *args):
        if isinstance(x, float):
            return self._approx_(x)
        try:
            return getattr(x, self._repr_())(*args)
        except AttributeError:
            return SymbolicComposition(self, SR(x))

    def _approx_(self, x):  # must *always* be called with a float x as input.
        s = '%s(%s), numer'%(self._repr_(), float(x))
        return float(maxima.eval(s))

    def _complex_approx_(self, x): # must be called with Python complex float as input
        """
        Given a Python complex $x$, evaluate self and return a complex value.

        EXAMPLES:
            sage: complex(cos(3*I))
            (10.067661995777771+0j)

        The following fails because we and Maxima haven't implemented
        erf yet for complex values:
            sage: complex(erf(3*I))
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to complex approximation
        """
        if x.imag == 0:
            return complex(self._approx_(x.real))
        s = '%s(%s+%s*%%i), numer'%(self._repr_(), x.real, x.imag)
        a = maxima.eval(s).replace('%i', '1j')
        if '(' in a:
            # unable to simplify to a complex -- still function calls there.
            raise TypeError, "unable to simplify to complex approximation"
        return complex(eval(a))

    def number_of_arguments(self):
        """
        Returns the number of arguments of self.

        EXAMPLES:
            sage: sin.variables()
            ()
            sage: sin.number_of_arguments()
            1
        """
        return 1

_syms = {}

class Function_erf(PrimitiveFunction):
    r"""
    The error function, defined as $\text{erf}(x) =
    \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt$.

    \sage currently only implements the error function (via a call to
    PARI) when the input is real.
    """

    def _repr_(self, simplify=True):
        return "erf"

    def _latex_(self):
        return "\\text{erf}"

    def _approx_(self, x):
        return float(1 - pari(float(x)).erfc())


erf = Function_erf()
_syms['erf'] = erf

class Function_abs(PrimitiveFunction):
    """
    The absolute value function.

    EXAMPLES:
        sage: var('x y')
        (x, y)
        sage: abs(x)
        abs(x)
        sage: abs(x^2 + y^2)
        y^2 + x^2
        sage: abs(-2)
        2
        sage: sqrt(x^2)
        sqrt(x^2)
        sage: abs(sqrt(x))
        sqrt(x)
    """
    def _repr_(self, simplify=True):
        return "abs"

    def _latex_(self):
        return "\\mathrm{abs}"

    def _latex_composition(self, x):
        """
        sage: f = sage.calculus.calculus.Function_abs()
        sage: latex(f)
        \mathrm{abs}
        sage: latex(abs(x))
        \left| x \right|
        """
        return "\\left| " + latex(x) + " \\right|"

    def _approx_(self, x):
        return float(x.__abs__())

    def _complex_approx_(self, x):
        """
        EXAMPLES:
            sage: complex(abs(3*I))
            (3+0j)
            sage: abs_symbolic._complex_approx_(complex(3*I))
            (3+0j)
        """
        return complex(x.__abs__())

    def __call__(self, x): # special case
        return SymbolicComposition(self, SR(x))

abs_symbolic = Function_abs()
_syms['abs'] = abs_symbolic



class Function_ceil(PrimitiveFunction):
    r"""
    The ceiling function.

    The ceiling of $x$ is computed in the following manner.
    \begin{enumerate}

    \item The \code{x.ceil()} method is called and returned if it is there.
         If it is not, then \sage checks if $x$ is one of Python's
         native numeric data types.  If so, then it calls
         and returns \code{Integer(int(math.ceil(x)))}.

    \item \sage tries to convert $x$ into a \class{RealIntervalField} with 53
          bits of precision. Next, the ceilings of the endpoints are computed.
          If they are the same, then that value is returned.  Otherwise, the
          precision of the \class{RealIntervalField} is increased until they
          do match up or it reaches \code{maximum_bits} of precision.

    \item If none of the above work, \sage returns a \class{SymbolicComposition}
         object.

    \end{enumerate}

    EXAMPLES:
        sage: a = ceil(2/5 + x)
        sage: a
        ceil(x + 2/5)
        sage: a(4)
        5
        sage: a(4.0)
        5
        sage: ZZ(a(3))
        4
        sage: a = ceil(x^3 + x + 5/2)
        sage: a
        ceil(x^3 + x + 1/2) + 2
        sage: a(x=2)
        13

        sage: ceil(log(8)/log(2))
        3

        sage: ceil(5.4)
        6
        sage: type(ceil(5.4))
        <type 'sage.rings.integer.Integer'>

        sage: ceil(factorial(50)/exp(1))
        11188719610782480504630258070757734324011354208865721592720336801
        sage: ceil(SR(10^50 + 10^(-50)))
        100000000000000000000000000000000000000000000000001
        sage: ceil(SR(10^50 - 10^(-50)))
        100000000000000000000000000000000000000000000000000
    """
    def _repr_(self, simplify=True):
        return "ceil"

    def _latex_(self):
        return "\\text{ceil}"

    def _maxima_init_(self):
        return "ceiling"

    _approx_ = math.ceil

    def __call__(self, x, maximum_bits=20000):
        try:
            return x.ceil()
        except AttributeError:
            if isinstance(x, (float, int, long, complex)):
                return Integer(int(math.ceil(x)))

        x_original = x
        if isinstance(x, SymbolicExpression):
            x = x.full_simplify()

        #If x can be coerced into a real interval, then we should
        #try increasing the number of bits of precision until
        #we get the ceiling at each of the endpoints is the same.
        #The precision will continue to be increased up to maximum_bits
        #of precision at which point it will raise a value error.
        bits = 53
        try:
            x_interval = RealIntervalField(bits)(x)
            upper_ceil = x_interval.upper().ceil()
            lower_ceil = x_interval.lower().ceil()

            while upper_ceil != lower_ceil and bits < maximum_bits:
                bits += 100
                x_interval = RealIntervalField(bits)(x)
                upper_ceil = x_interval.upper().ceil()
                lower_ceil = x_interval.lower().ceil()

            if bits < maximum_bits:
                return lower_ceil
            else:
                raise ValueError, "x (= %s) requires more than %s bits of precision to compute its ceiling"%(x, maximum_bits)

        except TypeError:
            #If x cannot be coerced into a RealField, then
            #it should be left as a symbolic expression.
            return SymbolicComposition(self, SR(x_original))

ceil = Function_ceil()
_syms['ceiling'] = ceil   # spelled ceiling in maxima


class Function_floor(PrimitiveFunction):
    r"""
    The floor function.

    The floor of $x$ is computed in the following manner.
    \begin{enumerate}
    \item The \code{x.floor()} method is called and returned if it is there.
         If it is not, then \sage checks if $x$ is one of Python's
         native numeric data types.  If so, then it calls
         and returns \code{Integer(int(math.floor(x)))}.

    \item \sage tries to convert $x$ into a \class{RealIntervalField} with
          53 bits of precision. Next, the floors of the endpoints are computed.
          If they are the same, then that value is returned.  Otherwise,
          the precision of the \class{RealIntervalField} is increased
          until they do match up or it reaches \code{maximum_bits} of precision.

    \item If none of the above work, \sage returns a \class{SymbolicComposition}
         object.

    \end{enumerate}

    EXAMPLES:
        sage: floor(5.4)
        5
        sage: type(floor(5.4))
        <type 'sage.rings.integer.Integer'>
        sage: var('x')
        x
        sage: a = floor(5.4 + x); a
        floor(x + 0.400000000000000) + 5
        sage: a(2)
        7

        sage: floor(log(8)/log(2))
        3

        sage: floor(factorial(50)/exp(1))
        11188719610782480504630258070757734324011354208865721592720336800
        sage: floor(SR(10^50 + 10^(-50)))
        100000000000000000000000000000000000000000000000000
        sage: floor(SR(10^50 - 10^(-50)))
        99999999999999999999999999999999999999999999999999
    """
    def _repr_(self, simplify=True):
        return "floor"

    def _latex_(self):
        return "\\text{floor}"

    def _maxima_init_(self):
        return "floor"

    _approx_ = math.floor

    def __call__(self, x, maximum_bits=20000):
        try:
            return x.floor()
        except AttributeError:
            if isinstance(x, (float, int, long, complex)):
                return Integer(int(math.floor(x)))

        x_original = x
        if isinstance(x, SymbolicExpression):
            x = x.full_simplify()

        #If x can be coerced into a real interval, then we should
        #try increasing the number of bits of precision until
        #we get the floor at each of the endpoints is the same.
        #The precision will continue to be increased up to maximum_bits
        #of precision at which point it will raise a value error.
        bits = 53
        try:
            x_interval = RealIntervalField(bits)(x)
            upper_floor = x_interval.upper().floor()
            lower_floor = x_interval.lower().floor()

            while upper_floor != lower_floor and bits < maximum_bits:
                bits += 100
                x_interval = RealIntervalField(bits)(x)
                upper_floor = x_interval.upper().floor()
                lower_floor = x_interval.lower().floor()

            if bits < maximum_bits:
                return lower_floor
            else:
                raise ValueError, "x (= %s) requires more than %s bits of precision to compute its floor"%(x, maximum_bits)

        except TypeError:
            #If x cannot be coerced into a RealField, then
            #it should be left as a symbolic expression.
            return SymbolicComposition(self, SR(x_original))


floor = Function_floor()
_syms['floor'] = floor   # spelled ceiling in maxima


class Function_sin(PrimitiveFunction):
    """
    The sine function
    """
    def _repr_(self, simplify=True):
        return "sin"

    def _latex_(self):
        return "\\sin"

    _approx_ = math.sin

    def __call__(self, x):
        try:
            return x.sin()
        except AttributeError:
            if isinstance(x, float):
                return math.sin(x)
        return SymbolicComposition(self, SR(x))

    def _fast_float_(self):
        return math.sin

sin = Function_sin()
_syms['sin'] = sin

class Function_cos(PrimitiveFunction):
    """
    The cosine function
    """
    def _repr_(self, simplify=True):
        return "cos"

    def _latex_(self):
        return "\\cos"

    _approx_ = math.cos

    def __call__(self, x):
        try:
            return x.cos()
        except AttributeError:
            if isinstance(x, float):
                return math.cos(x)
        return SymbolicComposition(self, SR(x))

    def _fast_float_(self):
        return math.cos


cos = Function_cos()
_syms['cos'] = cos

class Function_sec(PrimitiveFunction):
    """
    The secant function

    EXAMPLES:
        sage: sec(pi/4)
        sqrt(2)
        sage: RR(sec(pi/4))
        1.41421356237309
        sage: n(sec(pi/4),100)
        1.4142135623730950488016887242
        sage: sec(1/2)
        sec(1/2)
        sage: sec(0.5)
        1.13949392732455
    """
    def _repr_(self, simplify=True):
        return "sec"

    def _latex_(self):
        return "\\sec"

    def _approx_(self, x):
        return 1/math.cos(x)

sec = Function_sec()
_syms['sec'] = sec

class Function_csc(PrimitiveFunction):
    """
    The cosecant function.

    EXAMPLES:
        sage: csc(pi/4)
        sqrt(2)
        sage: RR(csc(pi/4))
        1.41421356237310
        sage: n(csc(pi/4),100)
        1.4142135623730950488016887242
        sage: csc(1/2)
        csc(1/2)
        sage: csc(0.5)
        2.08582964293349
    """
    def _repr_(self, simplify=True):
        return "csc"

    def _latex_(self):
        return "\\csc"

    def _approx_(self, x):
        return 1/math.sin(x)

csc = Function_csc()
_syms['csc'] = csc

class Function_cot(PrimitiveFunction):
    """
    The cotangent function.

    EXAMPLES:
        sage: cot(pi/4)
        1
        sage: RR(cot(pi/4))
        1.00000000000000
        sage: n(cot(pi/4),100)
        1.0000000000000000000000000000
        sage: cot(1/2)
        cot(1/2)
        sage: cot(0.5)
        1.83048772171245
    """
    def _repr_(self, simplify=True):
        return "cot"

    def _latex_(self):
        return "\\cot"

    def _approx_(self, x):
        return 1/math.tan(x)

cot = Function_cot()
_syms['cot'] = cot

class Function_tan(PrimitiveFunction):
    """
    The tangent function

    EXAMPLES:
        sage: tan(pi)
        0
        sage: tan(3.1415)
        -0.0000926535900581913
        sage: tan(3.1415/4)
        0.999953674278156
        sage: tan(pi/4)
        1
        sage: tan(1/2)
        tan(1/2)
        sage: RR(tan(1/2))
        0.546302489843790
    """
    def _repr_(self, simplify=True):
        return "tan"

    def _latex_(self):
        return "\\tan"

    def _approx_(self, x):
        return math.tan(x)

tan = Function_tan()
_syms['tan'] = tan

def arctan2(y, x):
    r"""
    Modified version of arctan function, since it is used by Maxima.

         \code{arctan2(y,x) = arctan(y/x)}

    This is mainly for internal use.

    TODO: entering 'atan2(1,2)' into Sage returns a  NameError that 'atan2'
    is not defined despite the two lines following this function definition.
    However, one can enter 'atan(1/2)' with no errors.

    EXAMPLES:
        sage: arctan2 = sage.calculus.calculus.arctan2
        sage: arctan2(1,2)
        arctan(1/2)
        sage: float(arctan2(1,2))
        0.46364760900080609
        sage: arctan2(2,3)
        arctan(2/3)
    """
    return arctan(y/x)

atan2 = arctan2
_syms['atan2'] = arctan2

class Function_arcsin(PrimitiveFunction):
    """
    The arcsine function

    EXAMPLES:
        sage: arcsin(0.5)
        0.523598775598299
        sage: arcsin(1/2)
        pi/6
        sage: arcsin(1 + I*1.0)
        1.061275061905036*I + 0.666239432492515
    """
    def _repr_(self, simplify=True):
        """
        Return string representation of arcsin.

        EXAMPLES:
            sage: arcsin._repr_()
            'arcsin'
        """
        return "arcsin"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arcsin._maxima_init_()
            'asin'
        """
        return "asin"

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES:
            sage: arcsin._latex_()
            '\\sin^{-1}'
        """
        return "\\sin^{-1}"

    def _approx_(self, x):
        """
        Return floating point approximation to the inverse of sine.

        EXAMPLES:
            sage: arcsin._approx_(0.5)
            0.52359877559829893
        """
        return math.asin(x)

arcsin = Function_arcsin()
asin = arcsin
_syms['asin'] = arcsin

class Function_arcsinh(PrimitiveFunction):
    """
    The inverse of the hyperbolic sine function.

    EXAMPLES:
        sage: arcsinh(0.5)
        0.481211825059603
        sage: arcsinh(1/2)
        arcsinh(1/2)
        sage: arcsinh(1 + I*1.0)
        0.666239432492515*I + 1.061275061905036
    """
    def _repr_(self, simplify=True):
        """
        Return string representation of arcsinh.

        EXAMPLES:
            sage: arcsinh._repr_()
            'arcsinh'
        """
        return "arcsinh"

    def _maxima_init_(self):
        """
        Return Maxima representation of this function.

        EXAMPLES:
            sage: arcsinh._maxima_init_()
            'asinh'
        """
        return "asinh"

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES:
            sage: arcsinh._latex_()
            '\\sinh^{-1}'
        """
        return "\\sinh^{-1}"

    def _approx_(self, x):
        """
        Return floating point numerical approximation to inverse hyperbolic sin at $x$.

        EXAMPLES:
            sage: arcsinh._approx_(0.5)
            0.48121182505960347
            sage: sinh(arcsinh._approx_(0.5))
            0.5
        """
        return float(pari(float(x)).asinh())

arcsinh = Function_arcsinh()
asinh = arcsinh
_syms['asinh'] = arcsinh

class Function_arccosh(PrimitiveFunction):
    """
    The inverse of the hyperbolic cosine function.

    EXAMPLES:
        sage: arccosh(1/2)
        arccosh(1/2)
        sage: arccosh(1 + I*1.0)
        0.904556894302381*I + 1.061275061905036

    Warning: If the input is real the output will be real or NaN:
        sage: arccosh(0.5)
        NaN

    But evaluation where the input is in the complex field yields a complex output:
        sage: arccosh(CC(0.5))
        1.04719755119660*I
    """
    def _repr_(self, simplify=True):
        """
        Return string representation of arccosh.

        EXAMPLES:
            sage: arccosh._repr_()
            'arccosh'
        """
        return "arccosh"

    def _maxima_init_(self):
        """
        Return Maxima representation of this function.

        EXAMPLES:
            sage: arccosh._maxima_init_()
            'acosh'
        """
        return "acosh"

    def _latex_(self):
        """
        Return latex representation of inverse cosine.

        EXAMPLES:
            sage: arccosh._latex_()
            '\\cosh^{-1}'
        """
        return "\\cosh^{-1}"

    def _approx_(self, x):
        """
        Return floating point approximation to arccosh.

        EXAMPLES:
            sage: float(arccosh(2))
            1.3169578969248168
            sage: cosh(float(arccosh(2)))
            2.0
        """
        return float(pari(float(x)).acosh())

arccosh = Function_arccosh()
acosh = arccosh
_syms['acosh'] = arccosh

class Function_arctanh(PrimitiveFunction):
    """
    The inverse of the hyperbolic tangent function.

    EXAMPLES:
        sage: arctanh(0.5)
        0.549306144334055
        sage: arctanh(1/2)
        arctanh(1/2)
        sage: arctanh(1 + I*1.0)
        1.017221967897851*I + 0.402359478108525
    """
    def _repr_(self, simplify=True):
        return "arctanh"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arctanh._maxima_init_()
            'atanh'
        """
        return "atanh"

    def _latex_(self):
        return "\\tanh^{-1}"

    def _approx_(self, x):
        return float(pari(float(x)).atanh())

arctanh = Function_arctanh()
atanh = arctanh
_syms['atanh'] = arctanh

class Function_arccoth(PrimitiveFunction):
    """
    The inverse of the hyperbolic cotangent function.

    EXAMPLES:
        sage: arccoth(2.)
        0.549306144334055
        sage: arccoth(2)
        arccoth(2)
        sage: arccoth(1 + I*1.0)
        0.402359478108525 - 0.553574358897045*I
    """
    def _repr_(self, simplify=True):
        return "arccoth"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arccoth._maxima_init_()
            'acoth'
        """
        return "acoth"

    def _latex_(self):
        return "\\coth^{-1}"

    def _approx_(self, x):
        return float(pari(float(1/x)).atanh())

arccoth = Function_arccoth()
acoth = arccoth
_syms['acoth'] = arccoth

class Function_arcsech(PrimitiveFunction):
    """
    The inverse of the hyperbolic secant function.

    EXAMPLES:
        sage: arcsech(.5)
        1.316957896924817
        sage: arcsech(1/2)
        arcsech(1/2)
        sage: arcsech(1 + I*1.0)
        0.530637530952518 - 1.118517879643706*I
    """
    def _repr_(self, simplify=True):
        return "arcsech"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arcsech._maxima_init_()
            'asech'
        """
        return "asech"

    def _latex_(self):
        return "\\sech^{-1}"

    def _approx_(self, x):
        return float(pari(float(1/x)).acosh())

arcsech = Function_arcsech()
asech = arcsech
_syms['asech'] = arcsech

class Function_arccsch(PrimitiveFunction):
    """
    The inverse of the hyperbolic cosecant function.

    EXAMPLES:
        sage: arccsch(2.)
        0.481211825059603
        sage: arccsch(2)
        arccsch(2)
        sage: arccsch(1 + I*1.0)
        0.530637530952518 - 0.452278447151191*I
    """
    def _repr_(self, simplify=True):
        return "arccsch"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arccsch._maxima_init_()
            'acsch'
        """
        return "acsch"

    def _latex_(self):
        return "\\csch^{-1}"

    def _approx_(self, x):
        return float(pari(float(1/x)).arcsinh())

arccsch = Function_arccsch()
acsch = arccsch
_syms['acsch'] = arccsch

class Function_arccos(PrimitiveFunction):
    """
    The arccosine function

    EXAMPLES:
        sage: arccos(0.5)
        1.04719755119660
        sage: arccos(1/2)
        pi/3
        sage: arccos(1 + I*1.0)
        0.904556894302381 - 1.061275061905036*I
    """
    def _repr_(self, simplify=True):
        return "arccos"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arccos._maxima_init_()
            'acos'
        """
        return "acos"

    def _latex_(self):
        return "\\cos^{-1}"

    def _approx_(self, x):
        return math.acos(x)

arccos = Function_arccos()
acos = arccos
_syms['acos'] = acos


class Function_arctan(PrimitiveFunction):
    """
    The arctangent function.

    EXAMPLES:
        sage: arctan(1/2)
        arctan(1/2)
        sage: RDF(arctan(1/2))
        0.463647609001
        sage: arctan(1 + I)
        arctan(I + 1)
    """
    def _repr_(self, simplify=True):
        return "arctan"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arctan._maxima_init_()
            'atan'
        """
        return "atan"

    def _latex_(self):
        return "\\tan^{-1}"

    def _approx_(self, x):
        return math.atan(x)

arctan = Function_arctan()
atan = arctan
_syms['atan'] = arctan

class Function_arccot(PrimitiveFunction):
    """
    The arccotangent function.

    EXAMPLES:
        sage: arccot(1/2)
        arccot(1/2)
        sage: RDF(arccot(1/2))
        1.10714871779
        sage: arccot(1 + I)
        arccot(I + 1)
    """
    def _repr_(self, simplify=True):
        return "arccot"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arccot._maxima_init_()
            'acot'
        """
        return "acot"

    def _latex_(self):
        return "\\cot^{-1}"

    def _approx_(self, x):
        return math.pi/2 - math.atan(x)

arccot = Function_arccot()
acot = arccot
_syms['acot'] = arccot

class Function_arccsc(PrimitiveFunction):
    """
    The arccosecant function.

    EXAMPLES:
        sage: arccsc(2)
        arccsc(2)
        sage: RDF(arccsc(2))
        0.523598775598
        sage: arccsc(1 + I)
        arccsc(I + 1)
    """
    def _repr_(self, simplify=True):
        return "arccsc"

    def _maxima_init_(self):
        """
        EXAMPLES:
            sage: arccsc._maxima_init_()
            'acsc'
        """
        return "acsc"

    def _latex_(self):
        return "\\csc^{-1}"

    def _approx_(self, x):
        return math.asin(1/x)

arccsc = Function_arccsc()
acsc = arccsc
_syms['acsc'] = arccsc

class Function_arcsec(PrimitiveFunction):
    """
    The arcsecant function.

    EXAMPLES:
        sage: arcsec(2)
        arcsec(2)
        sage: RDF(arcsec(2))
        1.0471975512
        sage: arcsec(1 + I)
        arcsec(I + 1)
    """
    def _repr_(self, simplify=True):
        return "arcsec"

    def _maxima_init_(self):
        """
        EXAMPLES:
            asec: arcsec._maxima_init_()
            'acsc'
        """
        return "asec"

    def _latex_(self):
        return "\\sec^{-1}"

    def _approx_(self, x):
        return math.acos(1/x)

arcsec = Function_arcsec()
asec = arcsec
_syms['asec'] = arcsec


#######
# Hyperbolic functions
#######

#tanh
class Function_tanh(PrimitiveFunction):
    """
    The hyperbolic tangent function.

    EXAMPLES:
        sage: tanh(pi)
        tanh(pi)
        sage: tanh(3.1415)
        0.996271386633702
        sage: float(tanh(pi))       # random low-order bits
        0.99627207622074987
        sage: tan(3.1415/4)
        0.999953674278156
        sage: tanh(pi/4)
        tanh(pi/4)
        sage: RR(tanh(1/2))
        0.462117157260010

        sage: CC(tanh(pi + I*e))
        0.997524731976164 - 0.00279068768100315*I
        sage: ComplexField(100)(tanh(pi + I*e))
        0.99752473197616361034204366446 - 0.0027906876810031453884245163923*I
        sage: CDF(tanh(pi + I*e))
        0.997524731976 - 0.002790687681*I
    """
    def _repr_(self, simplify=True):
        return "tanh"

    def _latex_(self):
        return "\\tanh"

    def _approx_(self, x):
        return math.tanh(x)

tanh = Function_tanh()
_syms['tanh'] = tanh

#sinh
class Function_sinh(PrimitiveFunction):
    """
    The hyperbolic sine function.

    EXAMPLES:
        sage: sinh(pi)
        sinh(pi)
        sage: sinh(3.1415)
        11.5476653707437
        sage: float(sinh(pi))              # random low-order bits
        11.548739357257748
        sage: RR(sinh(pi))
        11.5487393572577
    """
    def _repr_(self, simplify=True):
        return "sinh"

    def _latex_(self):
        return "\\sinh"

    def _approx_(self, x):
        return math.sinh(x)

sinh = Function_sinh()
_syms['sinh'] = sinh

#cosh
class Function_cosh(PrimitiveFunction):
    """
    The hyperbolic cosine function.

    EXAMPLES:
        sage: cosh(pi)
        cosh(pi)
        sage: cosh(3.1415)
        11.5908832931176
        sage: float(cosh(pi))       # random low order bits
        11.591953275521519
        sage: RR(cosh(1/2))
        1.12762596520638
    """
    def _repr_(self, simplify=True):
        return "cosh"

    def _latex_(self):
        return "\\cosh"

    def _approx_(self, x):
        return math.cosh(x)

cosh = Function_cosh()
_syms['cosh'] = cosh

#coth
class Function_coth(PrimitiveFunction):
    """
    The hyperbolic cotangent function.

    EXAMPLES:
        sage: coth(pi)
        coth(pi)
        sage: coth(3.1415)
        1.00374256795520
        sage: float(coth(pi))
        1.0037418731973213
        sage: RR(coth(pi))
        1.00374187319732
    """
    def _repr_(self, simplify=True):
        return "coth"

    def _latex_(self):
        return "\\coth"

    def _approx_(self, x):
        return 1/math.tanh(x)

coth = Function_coth()
_syms['coth'] = coth

#sech
class Function_sech(PrimitiveFunction):
    """
    The hyperbolic secant function.

    EXAMPLES:
        sage: sech(pi)
        sech(pi)
        sage: sech(3.1415)
        0.0862747018248192
        sage: float(sech(pi))    # random low order bits
        0.086266738334054432
        sage: RR(sech(pi))
        0.0862667383340544
    """
    def _repr_(self, simplify=True):
        return "sech"

    def _latex_(self):
        return "\\sech"

    def _approx_(self, x):
        return 1/math.cosh(x)

sech = Function_sech()
_syms['sech'] = sech


#csch
class Function_csch(PrimitiveFunction):
    """
    The hyperbolic cosecant function.

    EXAMPLES:
        sage: csch(pi)
        csch(pi)
        sage: csch(3.1415)
        0.0865975907592133
        sage: float(csch(pi))           # random low-order bits
        0.086589537530046945
        sage: RR(csch(pi))
        0.0865895375300470
    """
    def _repr_(self, simplify=True):
        return "csch"

    def _latex_(self):
        return "\\csch"

    def _approx_(self, x):
        return 1/math.sinh(x)

csch = Function_csch()
_syms['csch'] = csch

#############
# log
#############

class Function_log(PrimitiveFunction):
    """
    The log function.

    EXAMPLES:
        sage: log(e^2)
        2
        sage: log(1024, 2) # the following is ugly (for now)
        log(1024)/log(2)
        sage: log(10, 4)
        log(10)/log(4)

        sage: RDF(log(10,2))
        3.32192809489
        sage: RDF(log(8, 2))
        3.0
        sage: log(RDF(10))
        2.30258509299
        sage: log(2.718)
        0.999896315728952
    """
    def __init__(self):
        PrimitiveFunction.__init__(self)

    def _repr_(self, simplify=True):
        return "log"

    def _latex_(self):
        return "\\log"

    _approx_ = math.log

function_log = Function_log()

def ln(x):
    """
    The natural logarithm of x.

    INPUT:
        x -- positive real number

    OUTPUT:
        ln(x) -- real number

    EXAMPLES:
        sage: ln(e^2)
        2
        sage: ln(2)
        log(2)
        sage: ln(2.0)
        0.693147180559945
    """
    return function_log(x)

def log(x, base=None):
    """
    Return the logarithm of x to the given base.

    Calls the \code{log} method of the object x when computing the logarithm,
    thus allowing use of logarithm on any object containing a \code{log}
    method. In other words, log works on more than just real numbers.

    TODO: Add p-adic log example.

    EXAMPLES:
        sage: log(e^2)
        2
        sage: log(1024, 2); RDF(log(1024, 2))
        log(1024)/log(2)
        10.0
        sage: log(10, 4); RDF(log(10, 4))
        log(10)/log(4)
        1.66096404744

        sage: log(10, 2)
        log(10)/log(2)
        sage: n(log(10, 2))
        3.32192809488736
        sage: log(10, e)
        log(10)
        sage: n(log(10, e))
        2.30258509299405

    The log function also works in finite fields as long as the base is
    generator of the multiplicative group:
        sage: F = GF(13); g = F.multiplicative_generator(); g
        2
        sage: a = F(8)
        sage: log(a,g); g^log(a,g)
        3
        8
        sage: log(a,3)
        Traceback (most recent call last):
        ...
        ValueError: base (=3) for discrete log must generate multiplicative group
    """
    if base is None:
        try:
            return x.log()
        except AttributeError:
            return ln(x)
    else:
        try:
            return x.log(base)
        except AttributeError:
            return ln(x) / ln(base)

_syms['log'] = log
_syms['ln'] = log

#####################
# The polylogarithm
#####################


class Function_polylog(PrimitiveFunction):
    r"""
    The polylog function $\text{Li}_n(z) = \sum_{k=1}^{\infty} z^k / k^n$.

    INPUT:
        n -- object
        z -- object

    EXAMPLES:
        sage: f = polylog(1,x)._operands[0]; f
        polylog(1)
        sage: type(f)
        <class 'sage.calculus.calculus.Function_polylog'>
    """
    def __init__(self, n):
        """

        """
        PrimitiveFunction.__init__(self)
        self._n = n

    def _repr_(self, simplify=True):
        return "polylog(%s)"%self._n

    def _repr_evaled_(self, args):
        return 'polylog(%s, %s)'%(self._n, args)

    def _maxima_init_(self):
        """
        Return string representation of this polylog function in Maxima.

        EXAMPLES:
            sage: polylog(1,x)._operands[0]._maxima_init_()
            'li[1]'
            sage: polylog(2,x)._operands[0]._maxima_init_()
            'li[2]'
            sage: polylog(3,x)._operands[0]._maxima_init_()
            'li[3]'
            sage: polylog(4,x)._operands[0]._maxima_init_()
            'polylog(4)'
        """
        if self._n in [1,2,3]:
            return 'li[%s]'%self._n
        else:
            return 'polylog(%s)'%self._n

    def _maxima_init_evaled_(self, args):
        if self._n in [1,2,3]:
            return 'li[%s](%s)'%(self._n, args)
        else:
            return 'polylog(%s, %s)'%(self._n, args)

    def index(self):
        r"""
        Return the index of this polylogarithm, i.e., if this is $\text{Li}_n(z)$, then
        this function returns $n$.

        EXAMPLES:
            sage: a = polylog(5,x); a
            polylog(5, x)
            sage: a._operands
            [polylog(5), x]
            sage: a._operands[0].index()
            5
        """
        return self._n

    def _latex_(self):
        """
        Return Latex representation of this polylogarithm.

        EXAMPLES:
            sage: polylog(5,x)._operands[0]._latex_()
            '\\text{Li}_{5}'
        """
        return "\\text{Li}_{%s}"%(self._n)

    def _approx_(self, x):
        """
        Return real numerical approximation for this polylogarithm evaluated
        at $x$.

        EXAMPLES:
            sage: f = polylog(4,x)._operands[0]; f
            polylog(4)
            sage: f._approx_(1)
            1.0823232337111381
            sage: type(f._approx_(1))
            <type 'float'>
        """
        try:
            return float(pari(x).polylog(self._n))
        except PariError:
            raise TypeError, 'unable to coerce polylogarithm to float'

    def _complex_approx_(self, x):
        """
        Return real numerical approximation for this polylogarithm
        evaluated at $x$.

        EXAMPLES:
            sage: a = pari('1+I')
            sage: CDF(a)
            1.0 + 1.0*I
            sage: complex(polylog(4,2))
            (2.4278628067547032-0.17437130002545306j)
            sage: polylog(4,x)._operands[0]._complex_approx_(2)
            (2.4278628067547032-0.17437130002545306j)
        """
        try:
            # kind of lame using CDF here.
            return complex(CDF(pari(CDF(x)).polylog(self._n)))
        except PariError:
            raise TypeError, 'unable to coerce polylogarithm to float'

def polylog(n, z):
    r"""
    The polylogarithm function $\text{Li}_n(z) = \sum_{k=1}^{\infty} z^k / k^n$.

    EXAMPLES:
        sage: polylog(2,1)
        pi^2/6
        sage: polylog(2,x^2+1)
        polylog(2, x^2 + 1)
        sage: polylog(4,0.5)
        polylog(4, 0.500000000000000)
        sage: float(polylog(4,0.5))
        0.51747906167389934

        sage: var('z')
        z
        sage: polylog(2,z).taylor(z, 1/2, 3)
        -(6*log(2)^2 - pi^2)/12 + 2*log(2)*(z - 1/2) + (-2*log(2) + 2)*(z - 1/2)^2 + (8*log(2) - 4)*(z - 1/2)^3/3
    """
    return Function_polylog(n)(z)

def dilog(z):
    r"""
    The dilogarithm function $\text{Li}_2(z) = \sum_{k=1}^{\infty} z^k / k^2$.

    This is simply an alias for polylog(2, z).

    EXAMPLES:
        sage: dilog(1)
        pi^2/6
	sage: dilog(1/2)
	pi^2/12 - log(2)^2/2
	sage: dilog(x^2+1)
        polylog(2, x^2 + 1)
        sage: float(dilog(1))
	1.6449340668482264
	sage: var('z')
        z
	sage: dilog(z).diff(z, 2)
	1/((1 - z)*z) + log(1 - z)/z^2
	sage: dilog(z).taylor(z, 1/2, 3)
        -(6*log(2)^2 - pi^2)/12 + 2*log(2)*(z - 1/2) + (-2*log(2) + 2)*(z - 1/2)^2 + (8*log(2) - 4)*(z - 1/2)^3/3
    """
    return polylog(2, z)


_syms['polylog2'] = Function_polylog(2)
_syms['polylog3'] = Function_polylog(3)

##############################
# square root
##############################

class Function_sqrt(PrimitiveFunction):
    """
    The square root function. This is a symbolic square root.

    EXAMPLES:
        sage: sqrt(-1)
        I
        sage: sqrt(2)
        sqrt(2)
        sage: sqrt(x^2)
        sqrt(x^2)
    """
    def __init__(self):
        PrimitiveFunction.__init__(self, needs_braces=True)

    def _repr_(self, simplify=True):
        return "sqrt"

    def _latex_(self):
        return "\\sqrt"

    def _do_sqrt(self, x, prec=None, extend=True, all=False):
        if prec:
            return ComplexField(prec)(x).sqrt(all=all)
        z = SymbolicComposition(self, SR(x))
        if all:
            return [z, -z]
        return z

    def __call__(self, x, *args, **kwds):
        """
        INPUT:
            x -- a number
            prec -- integer (default: None): if None, returns an exact
                 square root; otherwise returns a numerical square
                 root if necessary, to the given bits of precision.
            extend -- bool (default: True); if True, return a square
                 root in an extension ring, if necessary. Otherwise,
                 raise a ValueError if the square is not in the base
                 ring.
            all -- bool (default: False); if True, return all square
                 roots of self, instead of just one.
        """
        if isinstance(x, float):
            return math.sqrt(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.sqrt(*args, **kwds)
            except AttributeError:
                pass
        return self._do_sqrt(x, *args, **kwds)

    def _approx_(self, x):
        return math.sqrt(x)

sqrt = Function_sqrt()
_syms['sqrt'] = sqrt

class Function_exp(PrimitiveFunction):
    r"""
    The exponential function, $\exp(x) = e^x$.

    EXAMPLES:
        sage: exp(-1)
        e^-1
        sage: exp(2)
        e^2
        sage: exp(x^2 + log(x))
        x*e^x^2
        sage: exp(2.5)
        12.1824939607035
        sage: exp(float(2.5))         # random low order bits
        12.182493960703473
        sage: exp(RDF('2.5'))
        12.1824939607
    """
    def __init__(self):
        PrimitiveFunction.__init__(self, needs_braces=True)

    def _repr_(self, simplify=True):
        return "exp"

    def _latex_(self):
        return "\\exp"

    def __call__(self, x):
        # if x is an integer or rational, never call the sqrt method
        if isinstance(x, float):
            return self._approx_(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.exp()
            except AttributeError:
                pass
        return SymbolicComposition(self, SR(x))


    def _approx_(self, x):
        return math.exp(x)

exp = Function_exp()
_syms['exp'] = exp


#######################################################
# Symbolic functions
#######################################################

class SymbolicFunction(PrimitiveFunction):
    """
    A formal symbolic function.

    EXAMPLES:
        sage: f = function('foo')
        sage: var('x,y,z')
        (x, y, z)
        sage: g = f(x,y,z)
        sage: g
        foo(x, y, z)
        sage: g(x=var('theta'))
        foo(theta, y, z)
    """
    def __init__(self, name):
        PrimitiveFunction.__init__(self, needs_braces=True)
        self._name = str(name)

    def __hash__(self):
        return hash(self._name)

    def _repr_(self, simplify=True):
        return self._name


    def _is_atomic(self):
        return True

    def _latex_(self):
        return "{\\rm %s}"%self._name

    def _maxima_init_(self):
        return "'%s"%self._name

    def _approx_(self, x):
        raise TypeError

    def __call__(self, *args):
        return SymbolicFunctionEvaluation(self, [SR(x) for x in args])



class SymbolicFunction_delayed(SymbolicFunction):
    def simplify(self):
        """
        Return the simplified form of this delayed function.  This
        always just returns this delayed function itself.

        OUTPUT:
            self

        EXAMPLES:
            sage: f = sage.calculus.calculus.symbolic_expression_from_maxima_string("?%jacobi_cd")
            sage: type(f)
            <class 'sage.calculus.calculus.SymbolicFunction_delayed'>
            sage: f.simplify()
            jacobi_cd
            sage: f.simplify() is f
            True
        """
        return self

    def _has_been_simplified(self):
        """
        Return True, since delayed symbolic functions are simplified
        by construction.

        OUTPUT:
            bool -- True

        EXAMPLES:
            sage: f = sage.calculus.calculus.symbolic_expression_from_maxima_string("?%jacobi_cd")
            sage: type(f)
            <class 'sage.calculus.calculus.SymbolicFunction_delayed'>
            sage: f._has_been_simplified()
            True
        """
        return True

    def _maxima_init_(self):
        """
        Return Maxima version of self.

        EXAMPLES:
            sage: f = sage.calculus.calculus.symbolic_expression_from_maxima_string("?%jacobi_cd")
            sage: f._maxima_init_()
            '?%jacobi_cd'
            sage: maxima(f)
            ?%jacobi_cd
        """
        return '?%%%s'%self._name

    def __call__(self, *args):
        """
        Call this delayed function evaluation at the given inputs.

        OUTPUT:
            a delayed function evaluation

        EXAMPLES:
            sage: f = sage.calculus.calculus.symbolic_expression_from_maxima_string("?%jacobi_cd")
            sage: f(2)
            jacobi_cd(2)
            sage: f(2,3)
            jacobi_cd(2, 3)
            sage: f(2,3,x)
            jacobi_cd(2, 3, x)
            sage: type(f(2,3,x))
            <class 'sage.calculus.calculus.SymbolicFunctionEvaluation_delayed'>
        """
        return SymbolicFunctionEvaluation_delayed(self, [SR(x) for x in args])



class SymbolicFunctionEvaluation(SymbolicExpression):
    """
    The result of evaluating a formal symbolic function.

    EXAMPLES:
        sage: h = function('gfun', x); h
        gfun(x)
        sage: k = h.integral(x); k
        integrate(gfun(x), x)
        sage: k(gfun=sin)
        -cos(x)
        sage: k(gfun=cos)
        sin(x)
        sage: k.diff(x)
        gfun(x)
    """
    def __init__(self, f, args=None):
        """
        INPUT:
            f -- symbolic function
            args -- a tuple or list of symbolic expressions, at which
                    f is formally evaluated.
        """
        SymbolicExpression.__init__(self)
        self._f = f
        if not args is None:
            if not isinstance(args, tuple):
                args = tuple(args)
        self._args = args

    def __float__(self):
        return float(maxima(self))

    def _is_atomic(self):
        return True

    def arguments(self):
        """
        Return arguments of self.

        EXAMPLES:
            sage: f = function('Gamma', var('z'), var('w'))
            sage: f.arguments()
            (z, w)
        """
        return tuple(self._args)

##     def keyword_arguments(self):
##         """
##         Return the keyword arguments to this formal function evaluation.

##         EXAMPLES:
##             sage: g = f(2,3,z=10); g
##             abc(2, 3, )
##             sage: type(g)
##             <class 'sage.calculus.calculus.SymbolicFunctionEvaluation'>
##             sage: g.keyword_arguments()
##             {'z': 10}
##         """
##         return self._kwds

    def _repr_(self, simplify=True):
        if simplify:
            return self.simplify()._repr_(simplify=False)
        else:
            args = ', '.join([x._repr_(simplify=simplify) for x in
                                                      self._args])
            #if not self._kwds is None:
            #    kwds = ', '.join(["%s=%s" %(x, y) for x,y in self._kwds.iteritems()])
            #    return '%s(%s, %s)' % (self._f._name, args, kwds)
            #else:
            return '%s(%s)' % (self._f._name, args)

    def _latex_(self):
        return "{\\rm %s}(%s)"%(self._f._name, ', '.join([x._latex_() for
                                                       x in self._args]))

    def _maxima_init_(self):
        r"""
        Return string that in Maxima evaluates to something
        equivalent to \code{self}.

        EXAMPLES:
            sage: f = function('Gamma', var('w'), var('theta')); f
            Gamma(w, theta)
            sage: f._maxima_init_()
            "'Gamma(w, theta)"
            sage: maxima(f(sqrt(2), theta+3))
            'Gamma(sqrt(2),theta+3)
        """
        try:
            return self.__maxima_init
        except AttributeError:
            n = self._f._name
            if not (n in ['integrate', 'diff']):
                n = "'" + n
            s = "%s(%s)"%(n, ', '.join([x._maxima_init_()
                                           for x in self._args]))
            self.__maxima_init = s
        return s

    def _mathematica_init_(self):
        r"""
        Return string that in Mathematica evaluates to something
        equivalent to \code{self}.

        EXAMPLES:
            sage: f = function('Gamma', var('z'))
            sage: mathematica(f) # optional
            Gamma[z]
            sage: f = function('Gamma', var('w'), var('z')); f
            Gamma(w, z)
            sage: f._mathematica_init_()
            'Gamma[w, z]'
            sage: mathematica(f(sqrt(2), z+1)) # optional
            Gamma[Sqrt[2], 1 + z]
        """
        n = self._f._name
        return "%s[%s]"%(n, ', '.join([x._mathematica_init_()
                                           for x in self._args]))

    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
            sage: y = var('y')
            sage: f = function('foo',x); f
            foo(x)
            sage: f(foo=sin)
            sin(x)
            sage: f(x+y)
            foo(y + x)
            sage: a = f(pi)
            sage: a.substitute(foo = sin)
            0
            sage: a = f(pi/2)
            sage: a.substitute(foo = sin)
            1

            sage: b = f(pi/3) + x + y
            sage: b
            y + x + foo(pi/3)
            sage: b(foo = sin)
            y + x + sqrt(3)/2
            sage: b(foo = cos)
            y + x + 1/2
            sage: b(foo = cos, x=y)
            2*y + 1/2
        """
        function_sub = False
        for x in kwds.keys():
            if repr(x) == self._f._name:
                g = kwds[x]
                function_sub = True
                break
        if function_sub:
            # Very important to make a copy, since we are mutating a dictionary
            # that will get used again by the calling function!
            kwds = dict(kwds)
            del kwds[x]

        arg = tuple([SR(x._recursive_sub(kwds)) for x in self._args])

        if function_sub:
            return g(*arg)
        else:
            return self.__class__(self._f, arg)

    def _recursive_sub_over_ring(self, kwds, ring):
        raise TypeError, "no way to coerce the formal function to ring."

    def variables(self):
        r"""
        Return the variables appearing in the simplified form of \code{self}.

        EXAMPLES:
            sage: foo = function('foo')
            sage: var('x,y,a,b,z,t')
            (x, y, a, b, z, t)
            sage: w = foo(x,y,a,b,z) + t
            sage: w
            foo(x, y, a, b, z) + t
            sage: w.variables()
            (a, b, t, x, y, z)
        """
        try:
            return self.__variables
        except AttributeError:
            pass
        vars = sum([list(op.variables()) for op in self._args], [])
        vars.sort(var_cmp)
        vars = tuple(vars)
        self.__variables = vars
        return vars


class SymbolicFunctionEvaluation_delayed(SymbolicFunctionEvaluation):
    """
    We create an example of such a thing.

    EXAMPLES:
        sage: f = sage.calculus.calculus.symbolic_expression_from_maxima_string("?%jacobi_cd")
        sage: type(f)
        <class 'sage.calculus.calculus.SymbolicFunction_delayed'>
        sage: g = f(10)
        sage: type(g)
        <class 'sage.calculus.calculus.SymbolicFunctionEvaluation_delayed'>
        sage: g
        jacobi_cd(10)
    """
    def simplify(self):
        return self

    def _has_been_simplified(self):
        return True

    def __float__(self):
        return float(self._maxima_())

    def __complex__(self):
        return complex(self._maxima_())

    def _real_double_(self, R):
        return R(float(self))

    def _real_rqdf_(self, R):
        raise TypeError

    def _complex_double_(self, C):
        return C(float(self))

    def _mpfr_(self, field):
        if field.prec() <= 53:
            return field(float(self))
        raise TypeError

    def _complex_mpfr_field_(self, field):
        if field.prec() <= 53:
            return field(float(self))
        raise TypeError

    def _maxima_init_(self):
        try:
            return self.__maxima_init
        except AttributeError:
            n = self._f._name
            s = "%s(%s)"%(n, ', '.join([x._maxima_init_()
                                           for x in self._args]))
            self.__maxima_init = s
        return s



_functions = {}
def function(s, *args):
    """
    Create a formal symbolic function with the name \emph{s}.

    EXAMPLES:
        sage: var('a, b')
        (a, b)
        sage: f = function('cr', a)
        sage: g = f.diff(a).integral(b)
        sage: g
        diff(cr(a), a, 1)*b
        sage: g(cr=cos)
        -sin(a)*b
        sage: g(cr=sin(x) + cos(x))
        (cos(a) - sin(a))*b

    Basic arithmetic:
        sage: x = var('x')
        sage: h = function('f',x)
        sage: 2*f
        2*f
        sage: 2*h
        2*f(x)
    """
    if len(args) > 0:
        return function(s)(*args)
    if isinstance(s, SymbolicFunction):
        return s
    s = str(s)
    if ',' in s:
        return tuple([function(x.strip()) for x in s.split(',')])
    elif ' ' in s:
        return tuple([function(x.strip()) for x in s.split()])
    try:
        f =  _functions[s]
        _syms[s] = f
        return f
    except KeyError:
        pass
    v = SymbolicFunction(s)
    _functions[s] = v
    _syms[s] = v
    return v

def dummy_limit(*args):
    """
    This function is called to create formal wrappers of limits that
    Maxima can't compute:

    EXAMPLES:
        sage: a = lim(exp(x^2)*(1-erf(x)), x=infinity); a
        limit(e^x^2 - e^x^2*erf(x), x, +Infinity)
        sage: a = sage.calculus.calculus.dummy_limit(sin(x)/x, x, 0);a
        limit(sin(x)/x, x, 0)
    """
    s = str(args[1])
    #return SymbolicFunctionEvaluation(function('limit'), args=(args[0],), kwds={s: args[2]})
    return SymbolicFunctionEvaluation(function('limit'), args=(args[0], var(s), SR(args[2])))

######################################i################




#######################################################

symtable = {'%pi':'pi', '%e': 'e', '%i':'I', '%gamma':'euler_gamma',
            'li[2]':'polylog2', 'li[3]':'polylog3'}

from sage.rings.infinity import infinity, minus_infinity

_syms['inf'] = infinity
_syms['minf'] = minus_infinity

from sage.misc.multireplace import multiple_replace

import re


maxima_tick = re.compile("'[a-z|A-Z|0-9|_]*")

maxima_qp = re.compile("\?\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

maxima_var = re.compile("\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

sci_not = re.compile("(-?(?:0|[1-9]\d*))(\.\d+)?([eE][-+]\d+)")

def symbolic_expression_from_maxima_string(x, equals_sub=False, maxima=maxima):
    """
    Given a string representation of a Maxima expression, parse it and
    return the corresponding Sage symbolic expression.

    INPUT:
        x -- a string
        equals_sub -- (default: False) if True, replace '=' by '==' in self
        maxima -- (default: the calculus package's Maxima) the Maxima
                  interpreter to use.

    EXAMPLES:
        sage: sage.calculus.calculus.symbolic_expression_from_maxima_string('x^%e + %e^%pi + %i + sin(0)')
        x^e + I + e^pi
    """
    syms = dict(_syms)

    if len(x) == 0:
        raise RuntimeError, "invalid symbolic expression -- ''"
    maxima.set('_tmp_',x)

    # This is inefficient since it so rarely is needed:
    #r = maxima._eval_line('listofvars(_tmp_);')[1:-1]

    s = maxima._eval_line('_tmp_;')

    formal_functions = maxima_tick.findall(s)
    if len(formal_functions) > 0:
        for X in formal_functions:
            syms[X[1:]] = function(X[1:])
        # You might think there is a potential very subtle bug if 'foo is in a string literal --
        # but string literals should *never* ever be part of a symbolic expression.
        s = s.replace("'","")

    delayed_functions = maxima_qp.findall(s)
    if len(delayed_functions) > 0:
        for X in delayed_functions:
            syms[X[2:]] = SymbolicFunction_delayed(X[2:])
        s = s.replace("?%","")

    s = multiple_replace(symtable, s)
    s = s.replace("%","")

    if equals_sub:
        s = s.replace('=','==')

    #replace all instances of Maxima's scientific notation
    #with regular notation
    search = sci_not.search(s)
    while not search is None:
        (start, end) = search.span()
        r = create_RealNumber(s[start:end]).str(no_sci=2, truncate=True)
        s = s.replace(s[start:end], r)
        search = sci_not.search(s)

    # have to do this here, otherwise maxima_tick catches it
    syms['limit'] = dummy_limit

    global is_simplified
    try:
        # use a global flag so all expressions obtained via
        # evaluation of maxima code are assumed pre-simplified
        is_simplified = True
        return symbolic_expression_from_string(s, syms, accept_sequence=True)
    except SyntaxError:
        raise TypeError, "unable to make sense of Maxima expression '%s' in SAGE"%s
    finally:
        is_simplified = False

def symbolic_expression_from_maxima_element(x):
    """
    Given an element of the calculus copy of the Maxima interface,
    create the corresponding Sage symbolic expression.

    EXAMPLES:
        sage: a = sage.calculus.calculus.maxima('x^(sqrt(y)+%pi) + sin(%e + %pi)')
        sage: sage.calculus.calculus.symbolic_expression_from_maxima_element(a)
        x^(sqrt(y) + pi) - sin(e)
    """
    return symbolic_expression_from_maxima_string(x.name())

def evaled_symbolic_expression_from_maxima_string(x):
    """
    Given a string expression that makes sense in Maxima, return the
    corresponding Sage symbolic expression.  This is used mainly
    internally by the Calculus module.

    EXAMPLES:
        sage: sage.calculus.calculus.evaled_symbolic_expression_from_maxima_string('2*x + x^3 + y*z*sin(sqrt(x)*erf(theta))')
        sin(erf(theta)*sqrt(x))*y*z + x^3 + 2*x
        sage: sage.calculus.calculus.evaled_symbolic_expression_from_maxima_string('x^%e + %e^%pi + %i')
        x^e + I + e^pi
    """
    return symbolic_expression_from_maxima_string(maxima.eval(x))

def first_var(expr):
    """
    Return the first variable in expr or `x' if there are no variables
    in expression.

    EXAMPLES:
        sage: var('a,x,y')
        (a, x, y)
        sage: sage.calculus.calculus.first_var(a + y^x)
        a
        sage: sage.calculus.calculus.first_var(y^x - x^3)
        x
    """
    v = expr.variables()
    if len(v) > 0:
        return v[0]
    else:
        return var('x')



# External access used by restore
syms_cur = _syms
syms_default = dict(syms_cur)

# Comma format options for Maxima
def mapped_opts(v):
    """
    Used internally when creating a string of options to pass to Maxima.

    INPUT:
        v -- an object
    OUTPUT:
        a string.

    The main use of this is to turn Python bools into lower case strings.

    EXAMPLES:
        sage: sage.calculus.calculus.mapped_opts(True)
        'true'
        sage: sage.calculus.calculus.mapped_opts(False)
        'false'
        sage: sage.calculus.calculus.mapped_opts('bar')
        'bar'
    """
    if isinstance(v, bool):
        return str(v).lower()
    return str(v)

def maxima_options(**kwds):
    """
    Used internally to create a string of options to pass to Maxima.

    EXAMPLES:
        sage: sage.calculus.calculus.maxima_options(an_option=True, another=False, foo='bar')
        'an_option=true,foo=bar,another=false'
    """
    return ','.join(['%s=%s'%(key,mapped_opts(val)) for key, val in kwds.iteritems()])


# Parser for symbolic ring elements

_augmented_syms = {}

def _find_var(name):
    try:
        return (_augmented_syms or _syms)[name]
    except KeyError:
        pass
    try:
        return SR(sage.all.__dict__[name])
    except (KeyError, TypeError):
        return var(name)

def _find_func(name):
    try:
        func = (_augmented_syms or _syms)[name]
        if not isinstance(func, (SymbolicConstant, SymbolicVariable)):
            return func
    except KeyError:
        pass
    try:
        func = SR(sage.all.__dict__[name])
        if not isinstance(func, (SymbolicConstant, SymbolicVariable)):
            return func
    except (KeyError, TypeError):
        return function(name)

SR_parser = Parser(make_int      = lambda x: SymbolicConstant(Integer(x)),
                   make_float    = lambda x: SymbolicConstant(create_RealNumber(x)),
                   make_var      = _find_var,
                   make_function = _find_func)

def symbolic_expression_from_string(s, syms=None, accept_sequence=False):
    parse_func = SR_parser.parse_sequence if accept_sequence else SR_parser.parse_expression
    if syms is None:
        return parse_func(s)
    else:
        try:
            global _augmented_syms
            _augmented_syms = syms
            return parse_func(s)
        finally:
            _augmented_syms = {}
