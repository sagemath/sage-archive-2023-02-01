r"""
Symbolic Computation.

AUTHORS:
    Bobby Moretti and William Stein: 2006--2007

The \sage calculus module is loosely based on the \sage Enhahcement Proposal
found at: http://www.sagemath.org:9001/CalculusSEP.

EXAMPLES:

    The basic units of the calculus package are symbolic expressions
    which are elements of the symbolic expression ring (SR). There are
    many subclasses of SymbolicExpression. The most basic of these is
    the formal indeterminate class, SymbolicVariable. To create a
    SymbolicVariable object in \sage, use the var() method, whose
    argument is the text of that variable.  Note that \sage is
    intelligent about {\latex}ing variable names.

        sage: x1 = var('x1'); x1
        x1
        sage: latex(x1)
        \mbox{x}_{1}
        sage: theta = var('theta'); theta
        theta
        sage: latex(theta)
        \theta

    \sage predefines upper and lowercase letters as global
    indeterminates. Thus the following works:
        sage: x^2
        x^2
        sage: type(x)
        <class 'sage.calculus.calculus.SymbolicVariable'>

    More complicated expressions in SAGE can be built up using
    ordinary arithmetic. The following are valid, and follow the rules
    of Python arithmetic: (The '=' operator represents assignment, and
    not equality)
        sage: f = x + y + z/(2*sin(y*z/55))
        sage: g = f^f; g
        (z/(2*sin(y*z/55)) + y + x)^(z/(2*sin(y*z/55)) + y + x)

    Differentiation and integration are available, but behind the
    scenes through maxima:

        sage: f = sin(x)/cos(2*y)
        sage: f.derivative(y)
        2*sin(x)*sin(2*y)/(cos(2*y)^2)
        sage: g = f.integral(x); g
        -cos(x)/(cos(2*y))

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
        sage: f({x: t, y: z})
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

    We can also make CallableSymbolicExpressions, which is a SymbolicExpression
    that are functions of variables in a fixed order. Each
    SymbolicExpression has a function() method used to create a
    CallableSymbolicExpression.

        sage: u = log((2-x)/(y+5))
        sage: f = u.function(x, y); f
        (x, y) |--> log((2 - x)/(y + 5))

    There is an easier way of creating a CallableSymbolicExpression, which
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

COERCION EXAMPLES:

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
We test pickling:
    sage: f = -sqrt(pi)*(x^3 + sin(x/cos(y)))
    sage: bool(loads(dumps(f)) == f)
    True

Substitution:
    sage: f = x
    sage: f(x=5)
    5

"""

import weakref

from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_MPolynomial, is_MPolynomialRing,
                            is_RealNumber, is_ComplexNumber, RR,
                            Integer, Rational, CC,
                            PolynomialRing)

from sage.structure.element import RingElement, is_Element
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex
from sage.structure.sage_object import SageObject

from sage.interfaces.maxima import MaximaElement, Maxima
from sage.interfaces.all import maxima

from sage.misc.sage_eval import sage_eval

from sage.calculus.equations import SymbolicEquation
from sage.rings.real_mpfr import RealNumber
from sage.rings.complex_number import ComplexNumber
from sage.rings.real_double import RealDoubleElement
from sage.rings.complex_double import ComplexDoubleElement
from sage.rings.real_mpfi import RealIntervalFieldElement
from sage.rings.infinity import InfinityElement

from sage.libs.pari.gen import pari, gen as PariGen

from sage.rings.complex_double import ComplexDoubleElement

import sage.functions.constants

import math

is_simplified = False

infixops = {operator.add: '+',
            operator.sub: '-',
            operator.mul: '*',
            operator.div: '/',
            operator.pow: '^'}


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


class SymbolicExpressionRing_class(CommutativeRing):
    """
    The ring of all formal symbolic expressions.

    EXAMPLES:
        sage: SR
        Symbolic Ring
        sage: type(SR)
        <class 'sage.calculus.calculus.SymbolicExpressionRing_class'>
    """
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        """
        Coerce x into the symbolic expression ring SR.

        EXAMPLES:
            sage: a = SR(-3/4); a
            -3/4
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: a.parent()
            Symbolic Ring

        If a is already in the symblic expression ring, coercing returns
        a itself (not a copy):
            sage: SR(a) is a
            True

        A Python complex number:
            sage: SR(complex(2,-3))
            2.00000000000000 - 3.00000000000000*I
        """
        if is_Element(x) and x.parent() is self:
            return x
        elif hasattr(x, '_symbolic_'):
            return x._symbolic_(self)
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, CallableSymbolicExpression):
            return x._expr
        elif isinstance(x, SymbolicExpression):
            return x
        elif isinstance(x, MaximaElement):
            return symbolic_expression_from_maxima_element(x)
        elif is_Polynomial(x) or is_MPolynomial(x):
            return SymbolicPolynomial(x)
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
                            InfinityElement
                            )):
            return SymbolicConstant(x)
        elif isinstance(x, complex):
            return evaled_symbolic_expression_from_maxima_string('%s+%%i*%s'%(x.real,x.imag))
        else:
            raise TypeError, 'cannot coerce %s into a SymbolicExpression.'%x

    def _repr_(self):
        return 'Symbolic Ring'

    def _latex_(self):
        return 'SymbolicExpressionRing'

    def var(self, x):
        return var(x)

    def characteristic(self):
        return Integer(0)

    def _an_element_impl(self):
        return zero_constant

    def is_field(self):
        return True

    def is_exact(self):
        return False

# Define the unique symbolic expression ring.
SR = SymbolicExpressionRing_class()

# The factory function that returns the unique SR.
def SymbolicExpressionRing():
    """
    Return the symbolic expression ring.

    EXAMPLES:
        sage: SymbolicExpressionRing()
        Symbolic Ring
        sage: SymbolicExpressionRing() is SR
        True
    """
    return SR

class SymbolicExpression(RingElement):
    """
    A Symbolic Expression.

    EXAMPLES:

    """
    def __init__(self):
        RingElement.__init__(self, SR)
        if is_simplified:
            self._simp = self

    def __nonzero__(self):
        # Best to error on side of being nonzero in most cases.
        return not bool(self == SR.zero_element())

    def __str__(self):
        """
        Printing an object explicitly gives ASCII art:

        EXAMPLES:
            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f
            y^2/(y + 1)^3 + x/(x - 1)^3
            sage: print f
                                              2
                                             y          x
                                          -------- + --------
                                                 3          3
                                          (y + 1)    (x - 1)

        """
        return self.display2d(onscreen=False)

    def display2d(self, onscreen=True):
        """
        Display self using ASCII art.

        INPUT:
            onscreen -- string (optional, default True) If True,
                displays; if False, returns string.

        EXAMPLES:
        We display a fraction:
            sage: f = (x^3+y)/(x+3*y^2+1); f
            (y + x^3)/(3*y^2 + x + 1)
            sage: print f
                                                     3
                                                y + x
                                             ------------
                                                2
                                             3 y  + x + 1

        Use onscreen=False to get the 2d string:
             sage: f.display2d(onscreen=False)
             '\t\t\t\t\t 3\r\n\t\t\t\t    y + x\r\n         \t\t\t ------------\r\n\t\t\t\t    2\r\n\t\t\t\t 3 y  + x + 1'

        ASCII art is really helps for the following integral:
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
        s = self._maxima_().display2d(onscreen=False)
        s = s.replace('%pi',' pi').replace('%i',' I').replace('%e', ' e')
        if onscreen:
            print s
        else:
            return s

    def is_simplified(self):
        return hasattr(self, '_simp') and self._simp is self

    def _declare_simplified(self):
        self._simp = self

    def hash(self):
        return hash(self._repr_(simplify=False))

    def plot(self, *args, **kwds):
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
                f = lambda x: F._obj(x)
            else:
                y = float(F._obj)
                f = lambda x: y

        elif param is None:
            if isinstance(self, CallableSymbolicExpression):
                A = self.arguments()
                if len(A) == 0:
                    raise ValueError, "function has no input arguments"
                else:
                    param = A[0]
                f = lambda x: self(x)
            else:
                A = self.variables()
                if len(A) == 0:
                    y = float(self)
                    f = lambda x: y
                else:
                    param = A[0]
                f = self.function(param)
        else:
            f = self.function(param)
        return plot(f, *args, **kwds)

    def __lt__(self, right):
        return SymbolicEquation(self, SR(right), operator.lt)

    def __le__(self, right):
        return SymbolicEquation(self, SR(right), operator.le)

    def __eq__(self, right):
        return SymbolicEquation(self, SR(right), operator.eq)

    def __ne__(self, right):
        return SymbolicEquation(self, SR(right), operator.ne)

    def __ge__(self, right):
        return SymbolicEquation(self, SR(right), operator.ge)

    def __gt__(self, right):
        return SymbolicEquation(self, SR(right), operator.gt)

    def __cmp__(self, right):
        """
        Compares self and right.

        This is by definition the comparison of the underlying Maxima
        objects.

        EXAMPLES:
        These two are equal:
            sage: cmp(e+e, e*2)
            0
        """
        return cmp(maxima(self), maxima(right))

    def _neg_(self):
        """
        Return the formal negative of self.

        EXAMPLES:
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
            "x^3 + pi^2 + e"
        """
        return '"%s"'%repr(self)

    def _singular_init_(self):
        """
        Conversion of symbolic object to Singular always results in a Singular string.

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


    def _mpfr_(self, field):
        raise TypeError

    def _complex_mpfr_field_(self, field):
        raise TypeError

    def _complex_double_(self, C):
        raise TypeError

    def _real_double_(self, R):
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
            sage: x + y
            y + x
            sage: x._add_(y)
            y + x
        """
        return SymbolicArithmetic([self, right], operator.add)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: x - y
            x - y
        """
        return SymbolicArithmetic([self, right], operator.sub)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x * y
            x*y
        """
        return SymbolicArithmetic([self, right], operator.mul)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: x / y
            x/y
        """
        return SymbolicArithmetic([self, right], operator.div)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: x^(n+1)
            x^(n + 1)
        """
        right = self.parent()(right)
        return SymbolicArithmetic([self, right], operator.pow)

    def variables(self, vars=tuple([])):
        """
        Return sorted list of variables that occur in the simplified
        form of self.

        OUTPUT:
            a Python set

        EXAMPLES:
            sage: f = x^(n+1) + sin(pi/19); f
            x^(n + 1) + sin(pi/19)
            sage: f.variables()
            (n, x)

            sage: a = e^x
            sage: a.variables()
            (x,)
        """
        return vars

    def _first_variable(self):
        try:
            return self.__first_variable
        except AttributeError:
            pass
        v = self.variables()
        if len(v) == 0:
            ans = var('x')
        else:
            ans = v[0]
        self.__first_variable = ans
        return ans

    def _has_op(self, operator):
        """
        Recursively searches for the given operator in a SymbolicExpression
        object.

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


    def __call__(self, dict=None, **kwds):
        return self.substitute(dict, **kwds)

    def power_series(self, base_ring):
        """
        Return algebraic power series associated to this symbolic
        expression, which must be a polynomial in one variable, with
        coefficients coercible to the base ring.

        The power series is truncated one more than the degree.

        EXAMPLES:
            sage: var('theta')
            theta
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
        Return self as an algebraic polynomial over the given base ring, if
        possible.

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
            sage: f = pi^3*x - y^2*e - I; f
            -1*e*y^2 + pi^3*x - I
            sage: f.polynomial(CDF)
            -1.0*I + (-2.71828182846)*y^2 + 31.0062766803*x
            sage: f.polynomial(CC)
            -1.00000000000000*I + (-2.71828182845905)*y^2 + 31.0062766802998*x
            sage: f.polynomial(ComplexField(70))
            -1.0000000000000000000*I + (-2.7182818284590452354)*y^2 + 31.006276680299820175*x

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
            2*y^35 + 3*x^35
            sage: parent(g)
            Polynomial Ring in x, y over Finite Field of size 7
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
        Coerce this symbolic expression to a polynomial in R.

        EXAMPLES:
            sage: R = QQ[x,y,z]
            sage: R(x^2 + y)
            y + x^2
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
            <class 'sage.rings.polynomial_element_generic.Polynomial_generic_dense'>
            sage: a.degree()
            0

        We coerce to a double precision complex polynomial ring:
            sage: f = e*x^3 + pi*y^3 + sqrt(2) + I; f
            pi*y^3 + e*x^3 + I + sqrt(2)
            sage: R = CDF[x,y]
            sage: R(f)
            1.41421356237 + 1.0*I + 3.14159265359*y^3 + 2.71828182846*x^3

        We coerce to a higher-precision polynomial ring
            sage: R = ComplexField(100)[x,y]
            sage: R(f)
            1.4142135623730950066967437806 + 1.0000000000000000000000000000*I + 3.1415926535897932384626433833*y^3 + 2.7182818284590452353602874714*x^3

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
            return R(B(self))
        return self.substitute_over_ring(dict(sub), ring=R)

    def function(self, *args):
        """
        Return a CallableSymbolicExpression, fixing a variable order
        to be the order of args.

        EXAMPLES:
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

        Using the function method we can obtain the above function f,
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
    def derivative(self, *args):
        """
        Returns the derivative of itself. If self has exactly one variable, then
        it differentiates with respect to that variable. If there is more than one
        variable in the expression, then you must explicitly supply a variable.
        If you supply a variable $x$ followed by a number $n$, then it will
        differentiate with respect to $n$ times with respect to $n$.

        You may supply more than one variable. Each variable may optionally be
        followed by a positive integer. Then SAGE will differentiate with
        respect to the first variable $n$ times, where $n$ is the number
        immediately following the variable in the parameter list. If the
        variable is not followed by an integer, then SAGE will differentiate
        once. Then SAGE will differentiate by the second variables, and if that
        is followed by a number $m$, it will differentiate $m$ times, and so on.

        EXAMPLES:
            sage: h = sin(x)/cos(x)
            sage: diff(h,x,x,x)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2
            sage: diff(h,x,3)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2

            sage: u = (sin(x) + cos(y))*(cos(x) - sin(y))
            sage: diff(u,x,y)
            sin(x)*sin(y) - cos(x)*cos(y)
            sage: f = ((x^2+1)/(x^2-1))^(1/4)
            sage: g = diff(f, x); g # this is a complex expression
            x/(2*(x^2 - 1)^(1/4)*(x^2 + 1)^(3/4)) - (x*(x^2 + 1)^(1/4)/(2*(x^2 - 1)^(5/4)))
            sage: g.simplify_rational()
            -x/((x^2 - 1)^(5/4)*(x^2 + 1)^(3/4))

            sage: f = y^(sin(x))
            sage: diff(f, x)
            cos(x)*y^sin(x)*log(y)

            sage: g(x) = sqrt(5-2*x)
            sage: g_3 = diff(g, x, 3); g_3(2)
            -3

            sage: f = x*e^(-x)
            sage: diff(f, 100)
            x*e^(-x) - 100*e^(-x)

            sage: g = 1/(sqrt((x^2-1)*(x+5)^6))
            sage: diff(g, x)
            -3/((x + 5)^3*sqrt(x^2 - 1)*abs(x + 5)) - (x/((x^2 - 1)^(3/2)*abs(x + 5)^3))
        """
        # check each time
        s = ""
        # see if we can implicitly supply a variable name
        try:
            a = args[0]
        except IndexError:
            # if there were NO arguments, try assuming
            a = 1
        if a is None or isinstance(a, (int, long, Integer)):
            vars = self.variables()
            if len(vars) == 1:
                s = "%s, %s" % (vars[0], a)
            else:
                raise ValueError, "must supply an explicit variable for an " +\
                                "expression containing more than one variable"
        for i in range(len(args)):
            if isinstance(args[i], SymbolicVariable):
                s = s + '%s, ' %repr(args[i])
                # check to see if this is followed by an integer
                try:
                    if isinstance(args[i+1], (int, long, Integer)):
                        s = s + '%s, ' %repr(args[i+1])
                    else:
                        s = s + '1, '
                except IndexError:
                    s = s + '1'
            elif isinstance(args[i], (int, long, Integer)):
                if args[i] == 0:
                    return self
                if args[i] < 0:
                    raise ValueError, "cannot take negative derivative"
            else:
                raise TypeError, "arguments must be integers or " +\
                                 "SymbolicVariable objects"

        try:
            if s[-2] == ',':
                s = s[:-2]
        except IndexError:
            pass
        t = maxima('diff(%s, %s)'%(self._maxima_().name(), s))
        f = self.parent()(t)
        return f

    differentiate = derivative
    diff = derivative


    ###################################################################
    # Taylor series
    ###################################################################
    def taylor(self, v, a, n):
        """
        Expands self in a truncated Taylor or Laurent series in the
        variable v around the point a, containing terms through $(x - a)^n$.

        INPUT:
            v -- variable
            a -- number
            n -- integer

        EXAMPLES:
            sage: taylor(a*log(z), z, 2, 3)
            log(2)*a + a*(z - 2)/2 - (a*(z - 2)^2/8) + a*(z - 2)^3/24
            sage: taylor(sqrt (sin(x) + a*x + 1), x, 0, 3)
            1 + (a + 1)*x/2 - ((a^2 + 2*a + 1)*x^2/8) + (3*a^3 + 9*a^2 + 9*a - 1)*x^3/48
            sage: taylor (sqrt (x + 1), x, 0, 5)
            1 + x/2 - (x^2/8) + x^3/16 - (5*x^4/128) + 7*x^5/256
            sage: taylor (1/log (x + 1), x, 0, 3)
            1/x + 1/2 - (x/12) + x^2/24 - (19*x^3/720)
            sage: taylor (cos(x) - sec(x), x, 0, 5)
            -x^2 - (x^4/6)
            sage: taylor ((cos(x) - sec(x))^3, x, 0, 9)
            -x^6 - (x^8/2)
            sage: taylor (1/(cos(x) - sec(x))^3, x, 0, 5)
            -1/x^6 + 1/(2*x^4) + 11/(120*x^2) - 347/15120 - (6767*x^2/604800) - (15377*x^4/7983360)
        """
        v = var(v)
        l = self._maxima_().taylor(v, SR(a), Integer(n))
        return self.parent()(l)

    ###################################################################
    # limits
    ###################################################################
    def limit(self, dir=None, **argv):
        """
        Return the limit as the variable v approaches a from the
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
            **argv -- 1 named parameter

        NOTE: Output it may also use `und' (undefined), `ind'
        (indefinite but bounded), and `infinity' (complex infinity).

        EXAMPLES:
            sage: f = (1+1/x)^x
            sage: f.limit(x = oo)
            e
            sage: f.limit(x = 5)
            7776/3125
            sage: f.limit(x = 1.2)
            2.069615754672029
            sage: f.limit(x = I)
            e^(I*log(1 - I))
            sage: f(1.2)
            2.069615754672029
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

            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (use assume):
            Is  x  positive or negative?

            sage: f = log(log(x))/log(x)
            sage: forget(); assume(x<-2); lim(f, x=0)
            limit(log(log(x))/log(x), x=0)

        The following means "indefinite but bounded":
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
            l = self._maxima_().limit(v, a)
        elif dir == 'plus' or dir == 'above':
            l = self._maxima_().limit(v, a, 'plus')
        elif dir == 'minus' or dir == 'below':
            l = self._maxima_().limit(v, a, 'minus')
        else:
            raise ValueError, "dir must be one of 'plus' or 'minus'"
        return self.parent()(l)

    ###################################################################
    # Laplace transform
    ###################################################################
    def laplace(self, t, s):
        r"""
        Attempts to compute and return the Laplace transform of self
        with respect to the variable t and transform parameter s.  If
        Laplace cannot find a solution, a formal function is returned.

        The function that is returned maybe be viewed as a function of s.

        DEFINITION:
        The Laplace transform of a function $f(t)$, defined for all
        real numbers $t \geq 0$, is the function $F(s)$ defined by
        $$
             F(s) = \int_{0}^{\infty} e^{-st} f(t) dt.
        $$

        EXAMPLES:
        We compute a few Laplace transforms:
            sage: sin(x).laplace(x, s)
            1/(s^2 + 1)
            sage: (z + exp(x)).laplace(x, s)
            z/s + 1/(s - 1)

            sage: var('t0')
            t0
            w
            sage: log(t/t0).laplace(t, s)
            (-log(t0) - log(s) - euler_gamma)/s

        We do a formal calculation:
            sage: f = function('f', x)
            sage: g = f.diff(x); g
            diff(f(x), x, 1)
            sage: g.laplace(x, s)
            s*laplace(f(x), x, s) - f(0)

        """
        return self.parent()(self._maxima_().laplace(var(t), var(s)))

    def inverse_laplace(self, t, s):
        r"""
        Attempts to compute the inverse Laplace transform of self with
        respect to the variable t and transform parameter s.  If
        Laplace cannot find a solution, a formal function is returned.

        The function that is returned maybe be viewed as a function of s.


        DEFINITION:
        The inverse Laplace transform of a function $F(s)$,
        is the function $f(t)$ defined by
        $$
             F(s) = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma + i\infty} e^{st} F(s) dt,
        $$
        where $\gamma$ is chosen so that the contour path of
        integration is in the region of convergence of $F(s)$.

        EXAMPLES:
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
        vars = self.variables()
        if len(vars) < 1:
            return var('x')
        else:
            return vars[0]

    ###################################################################
    # integration
    ###################################################################
    def integral(self, v=None, a=None, b=None):
        """
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

        Constraints are sometimes needed:
            sage: integral(x^n,x)
            Traceback (most recent call last):
            ...
            TypeError: Computation failed since Maxima requested additional constraints (use assume):
            Is  n+1  zero or nonzero?
            sage: assume(n > 0)
            sage: integral(x^n,x)
            x^(n + 1)/(n + 1)
            sage: forget()

        NOTE: Above, putting assume(n == -1) does not yield the right behavior.
        Directly in maxima, doing

        The examples in the Maxima documentation:
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

        We integrate the same function in both Mathematica and SAGE (via Maxima):
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
            log(x - 4)/73 - (integrate((x^2 + 4*x + 18)/(x^3 + 2*x + 1), x)/73)
            sage: print A
                                     /  2
                                     [ x  + 4 x + 18
                                     I ------------- dx
                                     ]  3
                        log(x - 4)   / x  + 2 x + 1
                        ---------- - ------------------
                            73               73
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
            return self.parent()(self._maxima_().integrate(v, a, b))

    integrate = integral

    def nintegral(self, x, a, b,
                  desired_relative_error='1e-8',
                  maximum_num_subintervals=200):
        r"""
        Return a numerical approximation to the integral of self from
        a to b.

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

        EXAMPLES:
            sage: exp(-sqrt(x)).nintegral(x, 0, 1)
            (0.52848223531423055, 4.1633141378838452e-11, 231, 0)

        We can also use the \code{numerical_integral} function, which calls
        the GSL C library.
            sage: numerical_integral(exp(-sqrt(x)), 0, 1)             # random low-order bits
            (0.52848223225314706, 6.8392846084921134e-07)
        """
        v = self._maxima_().quad_qags(var(x),
                                      a, b, desired_relative_error,
                                      maximum_num_subintervals)
        return float(v[0]), float(v[1]), Integer(v[2]), Integer(v[3])

    nintegrate = nintegral


    ###################################################################
    # Manipulating epxressions
    ###################################################################
    def coeff(self, x, n=1):
        """
        Returns the coefficient of $x^n$ in self.

        INPUT:
            x -- variable, function, expression, etc.
            n -- integer, default 1.

        Sometimes it may be necessary to expand or factor first, since
        this is not done automatically.

        EXAMPLES:
            sage: f = (a*sqrt(2))*x^2 + sin(y)*x^(1/2) + z^z
            sage: f.coeff(sin(y))
            sqrt(x)
            sage: f.coeff(x^2)
            sqrt(2)*a
            sage: f.coeff(x^(1/2))
            sin(y)
            sage: f.coeff(1)
            0
            sage: f.coeff(x, 0)
            z^z
        """
        return self.parent()(self._maxima_().coeff(x, n))

    def coeffs(self, x=None):
        """
        Coefficients of self as a polynomial in x.

        INPUT:
            x -- optional variable
        OUTPUT:
            list of pairs [expr, n], where expr is a symbolic expression and n is a power.

        EXAMPLES:
            sage: p = x^3 - (x-3)*(x^2+x) + 1
            sage: p.coeffs()
            [[1, 0], [3, 1], [2, 2]]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            x^2 - 2*sqrt(2)*a*x + x + 2*a^2 + 1
            sage: p.coeffs(a)
            [[x^2 + x + 1, 0], [-2*sqrt(2)*x, 1], [2, 2]]
            sage: p.coeffs(x)
            [[2*a^2 + 1, 0], [1 - 2*sqrt(2)*a, 1], [1, 2]]

        A polynomial with wacky exponents:
            sage: p = (17/3*a)*x^(3/2) + x*y + 1/x + x^x
            sage: p.coeffs(x)
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

    def poly(self, x=None):
        r"""
        Express self as a polynomial in x.  Is self is not a polynomial
        in x, then some coefficients may be functions of x.

        WARNING: This is different from \code{self.polynomial()} which
        returns a SAGE polynomial over a given base ring.

        EXAMPLES:
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
        """
        Simplifies self by combining all terms with the same
        denominator into a single term.

        EXAMPLES:
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
        EXAMPLES:
            sage: f = (sqrt(x) + sqrt(y) + sqrt(z))/(x^10 - y^10 - sqrt(var('theta')))
            sage: print f
                                      sqrt(z) + sqrt(y) + sqrt(x)
                                      ---------------------------
                                          10    10
                                       - y   + x   - sqrt(theta)
            sage: f.denominator()
            -y^10 + x^10 - sqrt(theta)
        """
        return self.parent()(self._maxima_().denom())

    ###################################################################
    # solve
    ###################################################################
    def roots(self, x=None):
        r"""
        Returns the roots of self, with multiplicities.

        INPUT:
            x -- variable to view the function in terms of
                   (use default variable if not given)
        OUTPUT:
            list of pairs (root, multiplicity)

        If there are infinitely many roots, e.g., a function
        like sin(x), only one is returned.

        EXAMPLES:
        A simple example:
            sage: ((x^2-1)^2).roots()
            [(-1, 2), (1, 2)]

        A complicated example.
            sage: f = expand((x^2 - 1)^3*(x^2 + 1)*(x-a)); f
            x^9 - a*x^8 - 2*x^7 + 2*a*x^6 + 2*x^3 - 2*a*x^2 - x + a

        The default variable is a, since it is the first in alphabetical order:
            sage: f.roots()
            [(x, 1)]

        As a polynomial in a, x is indeed a root:
            sage: f.poly(a)
            (x^9 - 2*x^7 + 2*x^3 - x)*1 + (-x^8 + 2*x^6 - 2*x^2 + 1)*a
            sage: f(a=x)
            0

        The roots in terms of x are what we expect:
            sage: f.roots(x)
            [(a, 1), (-1*I, 1), (I, 1), (1, 3), (-1, 3)]

        Only one root of $\sin(x) = 0$ is given:
            sage: f = sin(x)
            sage: f.roots(x)
            [(0, 1)]

        We derive the roots of a general quadratic polynomial:
            sage: (a*x^2 + b*x + c).roots(x)
            [((-sqrt(b^2 - 4*a*c) - b)/(2*a), 1), ((sqrt(b^2 - 4*a*c) - b)/(2*a), 1)]
        """
        if x is None:
            x = self.default_variable()
        S, mul = self.solve(x, multiplicities=True)
        return [(S[i].rhs(), mul[i]) for i in range(len(mul))]

    def solve(self, x, multiplicities=False):
        r"""
        Solve the equation \code{self == 0} for the variable x.

        INPUT:
            x -- variable to solve for
            multiplicities -- bool (default: False); if True, return corresponding multiplicities.

        EXAMPLES:
            sage: (z^5 - 1).solve(z)
            [z == e^(2*I*pi/5), z == e^(4*I*pi/5), z == e^(-(4*I*pi/5)), z == e^(-(2*I*pi/5)), z == 1]
        """
        x = var(x)
        return (self == 0).solve(x, multiplicities=multiplicities)

    ###################################################################
    # simplify
    ###################################################################
    def simplify(self):
        try:
            return self._simp
        except AttributeError:
            S = evaled_symbolic_expression_from_maxima_string(self._maxima_init_())
            S._simp = S
            self._simp = S
            return S

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
        Wraps the Maxima radcan() command. From the Maxima documentation:

            Simplifies this expression, which can contain logs, exponentials,
            and radicals, by converting it into a form which is canonical over a
            large class of expressions and a given ordering of variables; that
            is, all functionally equivalent forms are mapped into a unique form.
            For a somewhat larger class of expressions, produces a regular form.
            Two equivalent expressions in this class do not necessarily have the
            same appearance, but their difference can be simplified by radcan to
            zero.

            For some expressions radcan is quite time consuming. This is the
            cost of exploring certain relationships among the components of the
            expression for simplifications based on factoring and partial
            fraction expansions of exponents.

        ALIAS: radical_simplify, simplify_log, log_simplify, exp_simplify,
        simplify_exp are all the same

        EXAMPLES:

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
        return self.parent()(self._maxima_().radcan())

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
            sage: (x^3-y^3).factor()
            (-(y - x))*(y^2 + x*y + x^2)
            sage: factor(-8*y - 4*x + z^2*(2*y + x))
            (2*y + x)*(z - 2)*(z + 2)
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: F = factor(f/(36*(1 + 2*y + y^2)), dontfactor=[x])
            sage: print F
                                          2
                                        (x  + 2 x + 1) (y - 1)
                                        ----------------------
                                              36 (y + 1)
        """
        if len(dontfactor) > 0:
            m = self._maxima_()
            name = m.name()
            cmd = 'block([dontfactor:%s],factor(%s))'%(dontfactor, name)
            return evaled_symbolic_expression_from_maxima_string(cmd)
        else:
            return self.parent()(self._maxima_().factor())

    ###################################################################
    # expand
    ###################################################################
    def expand(self):
        """
        """
        return self.parent()(self._maxima_().expand())

    def expand_trig(self):
        """
        EXAMPLES:
            sage: sin(5*x).expand_trig()
            sin(x)^5 - 10*cos(x)^2*sin(x)^3 + 5*cos(x)^4*sin(x)

            sage: cos(2*x + y).trig_expand()
            cos(2*x)*cos(y) - sin(2*x)*sin(y)

        ALIAS: trig_expand and expand_trig are the same
        """
        return self.parent()(self._maxima_().trigexpand())

    trig_expand = expand_trig


    ###################################################################
    # substitute
    ###################################################################
    def substitute(self, in_dict=None, **kwds):
        """
        Takes the symbolic variables given as dict keys or as keywords and
        replaces them with the symbolic expressions given as dict values or as
        keyword values.  Also run when you call a SymbolicExpression.

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        EXAMPLES:
            sage: u = (x^3 - 3*y + 4*t)
            sage: u.substitute(x=y, y=t)
            y^3 + t

            sage: f = sin(x)^2 + 32*x^(y/2)
            sage: f(x=2, y = 10)
            sin(2)^2 + 1024

            sage: f(x=pi, y=t)
            32*pi^(t/2)

            sage: f = 2*x^2 - sin(x)
            sage: f(pi)
            2*pi^2

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
        return dict([(var(k),w) for k,w in kwds.iteritems()])

    def substitute_over_ring(self, in_dict=None, ring=None, **kwds):
        X = self.simplify()
        kwds = self.__parse_in_dict(in_dict, kwds)
        kwds = self.__varify_kwds(kwds)
        if ring is None:
            return X._recursive_sub(kwds)
        else:
            return X._recursive_sub_over_ring(kwds, ring)


    ###################################################################
    # Real and imaginary parts
    ###################################################################
    def real(self):
        """
        Return the real part of self.

        EXAMPLES:
            sage: a = log(3+4*I)
            sage: print a
                                             log(4  I + 3)
            sage: print a.real()
                                                log(5)
            sage: print a.imag()
                                                     4
                                                atan(-)
                                                     3

        Now make a and b symbolic and compute the general real part:
            sage: restore('a,b')
            sage: f = log(a + b*I)
            sage: f.real()
            log(b^2 + a^2)/2
        """
        return self.parent()(self._maxima_().real())

    def imag(self):
        """
        Return the imaginary part of self.

        EXAMPLES:
            sage: sqrt(-2).imag()
            sqrt(2)

        We simplify Ln(Exp(z)) to z for -Pi<Im(z)<=Pi:

            sage: f = log(exp(z))
            sage: assume(-pi < imag(z))
            sage: assume(imag(z) <= pi)
            sage: print f
        			       z
            sage: forget()

        A more symbolic example:
            sage: f = log(a + b*I)
            sage: f.imag()
            atan(b/a)
        """
        return self.parent()(self._maxima_().imag())

    def conjugate(self):
        """
        The complex conjugate of self.

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
        """
        The complex norm of self, i.e., self times its complex conjugate.

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
        """
        Return the partial fraction expansion of self with respect to
        the given variable.

        INPUT:
            var -- variable name or string (default: first variable)

        OUTPUT:
            Symbolic expression

        EXAMPLES:
            sage: f = x^2/(x+1)^3
            sage: f.partial_fraction()
            1/(x + 1) - (2/(x + 1)^2) + 1/(x + 1)^3
            sage: print f.partial_fraction()
                                        1        2          1
                                      ----- - -------- + --------
                                      x + 1          2          3
                                              (x + 1)    (x + 1)

        Notice that the first variable in the expression is used by default:
            sage: f = y^2/(y+1)^3
            sage: f.partial_fraction()
            1/(y + 1) - (2/(y + 1)^2) + 1/(y + 1)^3

            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f.partial_fraction()
            y^2/(y^3 + 3*y^2 + 3*y + 1) + 1/(x - 1)^2 + 1/(x - 1)^3

        You can explicitly specify which variable is used.
            sage: f.partial_fraction(y)
            1/(y + 1) - (2/(y + 1)^2) + 1/(y + 1)^3 + x/(x^3 - 3*x^2 + 3*x - 1)
        """
        if var is None:
            var = self._first_variable()
        return self.parent()(self._maxima_().partfrac(var))




class Symbolic_object(SymbolicExpression):
    r"""
    A class representing a symbolic expression in terms of a SageObject (not
    necessarily a 'constant').
    """
    def __init__(self, obj):
        SymbolicExpression.__init__(self)
        self._obj = obj

    def obj(self):
        """
        EXAMPLES:
        """
        return self._obj

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
        Symbolic_object.__init__(self, x)

    def _is_atomic(self):
        try:
            return self._atomic
        except AttributeError:
            if isinstance(self, Rational):
                self._atomic = False
            else:
                self._atomic = True
            return self._atomic

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
        sage: f(x=y)
        y^3 + y

    A multivariate polynomial:

        sage: R.<x,y,theta> = ZZ[]
        sage: f = SR(x^3 + x + y + theta^2); f
        theta^2 + y + x + x^3
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

    def polynomial(self, base_ring):
        """
        Return self as a polynomial over the given base ring, if possible.

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

##################################################################


zero_constant = SymbolicConstant(Integer(0))

class SymbolicOperation(SymbolicExpression):
    r"""
    A parent class representing any operation on SymbolicExpression objects.
    """
    def __init__(self, operands):
        SymbolicExpression.__init__(self)
        self._operands = operands   # don't even make a copy -- ok, since immutable.

    def variables(self, vars=tuple([])):
        """
        Return sorted list of variables that occur in the simplified
        form of self.  The ordering is alphabetic.

        EXAMPLES:
            sage: f = (x - x) + y^2 - z/z + (w^2-1)/(w+1); f
            y^2 + (w^2 - 1)/(w + 1) - 1
            sage: f.variables()
            (w, y)

            sage: (x + y + z + a + b + c).variables()
            (a, b, c, x, y, z)

            sage: (x^2 + x).variables()
            (x,)
        """
        if not self.is_simplified():
            return self.simplify().variables(vars)

        try:
            return self.__variables
        except AttributeError:
            pass
        if vars is None:
            vars = []
        else:
            vars = list(vars)
        vars = list(set(sum([list(op.variables()) for op in self._operands], [])))
        vars.sort(var_cmp)
        vars = tuple(vars)
        self.__variables = vars
        return vars

def var_cmp(x,y):
    return cmp(str(x), str(y))

symbols = {operator.add:' + ', operator.sub:' - ', operator.mul:'*',
            operator.div:'/', operator.pow:'^'}


class SymbolicArithmetic(SymbolicOperation):
    r"""
    Represents the result of an arithemtic operation on
    $f$ and $g$.
    """
    def __init__(self, operands, op):
        SymbolicOperation.__init__(self, operands)
        self._operator = op

    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
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
        """
        ops = self._operands
        new_ops = [SR(op._recursive_sub(kwds)) for op in ops]
        return self._operator(*new_ops)

    def _recursive_sub_over_ring(self, kwds, ring):
        ops = self._operands
        if self._operator == operator.pow:
            new_ops = [ops[0]._recursive_sub_over_ring(kwds, ring=ring), Integer(ops[1])]
        else:
            new_ops = [op._recursive_sub_over_ring(kwds, ring=ring) for op in ops]
        return ring(self._operator(*new_ops))

    def __float__(self):
        fops = [float(op) for op in self._operands]
        return self._operator(*fops)

    def __complex__(self):
        fops = [complex(op) for op in self._operands]
        return self._operator(*fops)

    def _mpfr_(self, field):
        if not self.is_simplified():
            return self.simplify()._mpfr_(field)
        rops = [op._mpfr_(field) for op in self._operands]
        return self._operator(*rops)

    def _complex_mpfr_field_(self, field):
        if not self.is_simplified():
            return self.simplify()._complex_mpfr_field_(field)
        rops = [op._complex_mpfr_field_(field) for op in self._operands]
        return self._operator(*rops)

    def _complex_double_(self, field):
        if not self.is_simplified():
            return self.simplify()._complex_double_(field)
        rops = [op._complex_double_(field) for op in self._operands]
        return self._operator(*rops)

    def _real_double_(self, field):
        if not self.is_simplified():
            return self.simplify()._real_double_(field)
        rops = [op._real_double_(field) for op in self._operands]
        return self._operator(*rops)

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
            sage: a = (1-1/r)^(-1); a
            1/(1 - (1/r))
            sage: a.derivative(r)
            -1/((1 - (1/r))^2*r^2)

            sage: reset('a,b')
            sage: s = 0*(1/a) + -b*(1/a)*(1 + -1*0*(1/a))*(1/(a*b + -1*b*(1/a)))
            sage: s
            -b/(a*(a*b - (b/a)))
            sage: s(a=2,b=3)
            -1/3
            sage: -3/(2*(2*3-(3/2)))
            -1/3
        """
        if simplify:
            if hasattr(self, '_simp'):
                return self._simp._repr_(simplify=False)
            else:
                return self.simplify()._repr_(simplify=False)

        ops = self._operands
        op = self._operator

        ###############
        # some bugs here in parenthesis -- exposed by above doctest
        ###############

        s = [o._repr_(simplify=False) for o in ops]

        # for the left operand, we need to surround it in parens when the
        # operator is mul/div/pow, and when the left operand contains an
        # operation of lower precedence
        if op in [operator.mul, operator.div]:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub):
                if not ops[0]._is_atomic():
                    s[0] = '(%s)' % s[0]
            else:
                try:
                    if isinstance(ops[0]._obj, Rational):
                        s[0] = '(%s)' % s[0]
                except AttributeError:
                    pass
                try:
                    if isinstance(ops[1]._obj, Rational):
                        s[1] = '(%s)' % s[1]
                except AttributeError:
                    pass

        # for the right operand, we need to surround it in parens when
        # the operation is mul/div/sub, and when the right operand
        # contains a + or -.
        if op in [operator.mul, operator.sub]:
                # avoid drawing parens if s1 an atomic operation
                if not ops[1]._is_atomic():
                    s[1] = '(%s)' % s[1]

        elif op is operator.div:
            if not ops[1]._is_atomic() or ops[1]._has_op(operator.mul):
                s[1] = '(%s)' % s[1]

        elif op is operator.pow:
            if not ops[0]._is_atomic():
                s[0] = '(%s)'% s[0]
            if not ops[1]._is_atomic() or ('/' in s[1] or '*' in s[1]):
                s[1] = '(%s)'% s[1]

        if op is operator.neg:
            if ops[0]._is_atomic():
                return '-%s' % s[0]
            else:
                return '-(%s)'%s[0]
        else:
            return '%s%s%s' % (s[0], symbols[op], s[1])

    def _latex_(self):
        # if we are not simplified, return the latex of a simplified version
        if not self.is_simplified():
            return self.simplify()._latex_()
        op = self._operator
        ops = self._operands
        s = [x._latex_() for x in self._operands]

        if op is operator.add:
            return '%s + %s' % (s[0], s[1])
        elif op is operator.sub:
            return '%s - %s' % (s[0], s[1])
        elif op is operator.mul:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub):
                s[0] = r'\left( %s \right)' %s[0]
            return '{%s \\cdot %s}' % (s[0], s[1])
        elif op is operator.div:
            return '\\frac{%s}{%s}' % (s[0], s[1])
        elif op is operator.pow:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub) \
               or ops[0]._has_op(operator.mul) or ops[0]._has_op(operator.div) \
               or ops[0]._has_op(operator.pow):
                s[0] = r'\left( %s \right)' % s[0]
            return '{%s}^{%s} ' % (s[0], s[1])
        elif op is operator.neg:
            return '-%s' % s[0]

    def _maxima_init_(self):
        ops = self._operands
        if self._operator is operator.neg:
            return '-(%s)' % ops[0]._maxima_init_()
        else:
            return '(%s) %s (%s)' % (ops[0]._maxima_init_(),
                             infixops[self._operator],
                             ops[1]._maxima_init_())

    def _sys_init_(self, system):
        ops = self._operands
        if self._operator is operator.neg:
            return '-(%s)' % sys_init(ops[0], system)
        else:
            return '(%s) %s (%s)' % (sys_init(ops[0], system),
                             infixops[self._operator],
                             sys_init(ops[1], system))

import re

class SymbolicVariable(SymbolicExpression):
    def __init__(self, name):
        SymbolicExpression.__init__(self)
        self._name = name
        if len(name) == 0:
            raise ValueError, "variable name must be nonempty"

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

    def variables(self, vars=tuple([])):
        """
        Return sorted list of variables that occur in the simplified
        form of self.
        """
        return (self, )

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
        a = self._name
        if len(a) > 1:
            m = re.search('(\d|[.,])+$',a)
            if m is None:
                a = tex_varify(a)
            else:
                b = a[:m.start()]
                a = '%s_{%s}'%(tex_varify(b), a[m.start():])

        self.__latex = a
        return a

    def _maxima_init_(self):
        return self._name

    def _sys_init_(self, system):
        return self._name

common_varnames = ['alpha',
                   'beta',
                   'gamma',
                   'Gamma',
                   'delta',
                   'Delta',
                   'epsilon',
                   'zeta',
                   'eta',
                   'theta',
                   'Theta'
                   'iota',
                   'kappa',
                   'lambda',
                   'Lambda',
                   'mu',
                   'nu',
                   'xi',
                   'Xi'
                   'pi',
                   'Pi'
                   'rho',
                   'sigma',
                   'Sigma'
                   'tau',
                   'upsilon',
                   'varphi',
                   'chi',
                   'psi',
                   'Psi'
                   'omega',
                   'Omega']


def tex_varify(a):
    if a in common_varnames:
        return "\\" + a
    else:
        return '\\mbox{%s}'%a

_vars = {}
def var(s, create=True):
    r"""
    Create a symbolic variable with the name \emph{s}.

    EXAMPLES:
        sage: var('xx')
        xx
    """
    if isinstance(s, SymbolicVariable):
        return s
    s = str(s)
    if ',' in s:
        return tuple([var(x.strip()) for x in s.split(',')])
    elif ' ' in s:
        return tuple([var(x.strip()) for x in s.split()])
    try:
        X = _vars[s]()
        if not X is None:
            return X
    except KeyError:
        if not create:
            raise ValueError, "the variable '%s' has not been defined"%var
        pass
    v = SymbolicVariable(s)
    _vars[s] = weakref.ref(v)
    _syms[s] = v
    return v


#########################################################################################
#  Callable functions
#########################################################################################
def is_CallableSymbolicExpressionRing(x):
    """
    EXAMPLES:
        sage: is_CallableSymbolicExpressionRing(QQ)
        False
        sage: is_CallableSymbolicExpressionRing(CallableSymbolicExpressionRing((x,y,z)))
        True
    """
    return isinstance(x, CallableSymbolicExpressionRing_class)

class CallableSymbolicExpressionRing_class(CommutativeRing):
    def __init__(self, args):
        self._default_precision = 53 # default precision bits
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
            if isinstance(x, CallableSymbolicExpression):
                x = x._expr
            return CallableSymbolicExpression(self, x)
        return self._coerce_try(x, [SR])

    def _repr_(self):
        return "Callable function ring with arguments %s"%(self._args,)

    def args(self):
        return self._args

    def arguments(self):
        return self.args()

    def zero_element(self):
        try:
            return self.__zero_element
        except AttributeError:
            z = CallableSymbolicExpression(SR.zero_element(), self._args)
            self.__zero_element = z
            return z

    def _an_element_impl(self):
        return self.zero_element()


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
    """
    Returns true if x is a callable symbolic expression.

    EXAMPLES:
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

    def expression(self):
        """
        Return the underlying symbolic expression (i.e., forget the
        extra map structure).
        """
        return self._expr

    def args(self):
        return self.parent().args()

    def _maxima_init_(self):
        return self._expr._maxima_init_()

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
        args = self.args()
        if len(args) == 1:
            return "%s \\ {\mapsto}\\ %s" % (args[0],
                    self._expr._latex_())
        else:
            vars = ", ".join(map(latex, args))
            # the weird TeX is to workaround an apparent JsMath bug
            return "\\left(%s \\right)\\ {\\mapsto}\\ %s" % (args, self._expr._latex_())

    def _neg_(self):
        return CallableSymbolicExpression(self.parent(), -self._expr)

    def __add__(self, right):
        """
        EXAMPLES:
            sage: f(x,n,y) = x^n + y^m;  g(x,n,m,z) = x^n +z^m
            sage: f + g
            (x, n, m, y, z) |--> z^m + y^m + 2*x^n
            sage: g + f
            (x, n, m, y, z) |--> z^m + y^m + 2*x^n
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
        Takes the variable list from another CallableSymbolicExpression object and
        compares it with the current CallableSymbolicExpression object's variable list,
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

        Note: When used for arithmetic between CallableSymbolicExpressions,
        these rules ensure that the set of CallableSymbolicExpressions will have
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
                return self._simp._repr_(simplify=False)
            else:
                return self.simplify()._repr_(simplify=False)
        ops = self._operands
        return "%s(%s)"% (ops[0]._repr_(simplify=False), ops[1]._repr_(simplify=False))

    def _latex_(self):
        if not self.is_simplified():
            return self.simplify()._latex_()
        ops = self._operands
        # certain functions (such as \sqrt) need braces in LaTeX
        if (ops[0]).tex_needs_braces():
            return r"%s{ %s }" % ( (ops[0])._latex_(), (ops[1])._latex_())
        # ... while others (such as \cos) don't
        return r"%s \left( %s \right)"%((ops[0])._latex_(),(ops[1])._latex_())

    def _maxima_init_(self):
        ops = self._operands
        return '%s(%s)' % (ops[0]._maxima_init_(), ops[1]._maxima_init_())

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

    def __complex__(self):
        f = self._operands[0]
        g = self._operands[1]
        return complex(f._approx_(float(g)))

    def _mpfr_(self, field):
        """
        Coerce to a multiprecision real number.

        EXAMPLES:
            sage: RealField(100)(sin(2)+cos(2))
            0.49315059027853930839845163641

            sage: RR(sin(pi))
            0.000000000000000

            sage: type(RR(sqrt(163)*pi))
            <type 'sage.rings.real_mpfr.RealNumber'>

            sage: RR(coth(pi))
            1.00374187319732
            sage: RealField(100)(coth(pi))
            1.0037418731973212882015526912
            sage: RealField(200)(acos(1/10))
            1.4706289056333368228857985121870581235299087274579233690964
        """
        if not self.is_simplified():
            return self.simplify()._mpfr_(field)
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
        if not self.is_simplified():
            return self.simplify()._complex_mpfr_field_(field)
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
        if not self.is_simplified():
            return self.simplify()._complex_double_(field)
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
        if not self.is_simplified():
            return self.simplify()._real_double_(field)
        f = self._operands[0]
        g = self._operands[1]
        z = f(g._real_double_(field))
        if isinstance(z, SymbolicExpression):
            return field(float(z))
        return z

class PrimitiveFunction(SymbolicExpression):
    def __init__(self, needs_braces=False):
        self._tex_needs_braces = needs_braces

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

_syms = {}

class Function_erf(PrimitiveFunction):
    r"""
    The error function, defined as $\text{erf}(x) =
    \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2} dt$.
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
        sage: abs(x)
        abs(x)
        sage: abs(x^2 + y^2)
        y^2 + x^2
        sage: abs(-2)
        2
        sage: sqrt(x^2)
        abs(x)
        sage: abs(sqrt(x))
        sqrt(x)
    """
    def _repr_(self, simplify=True):
        return "abs"

    def _latex_(self):
        return "\\abs"

    def _approx_(self, x):
        return float(x.__abs__())

    def __call__(self, x): # special case
        return SymbolicComposition(self, SR(x))

abs_symbolic = Function_abs()
_syms['abs'] = abs_symbolic



class Function_ceil(PrimitiveFunction):
    """
    The ceiling function.

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

        sage: ceil(5.4)
        6
        sage: type(ceil(5.4))
        <type 'sage.rings.integer.Integer'>
    """
    def _repr_(self, simplify=True):
        return "ceil"

    def _latex_(self):
        return "\\text{ceil}"

    def _maxima_init_(self):
        return "ceiling"

    _approx_ = math.ceil

    def __call__(self, x):
        try:
            return x.ceil()
        except AttributeError:
            if isinstance(x, float):
                return math.ceil(x)
        return SymbolicComposition(self, SR(x))

ceil = Function_ceil()
_syms['ceiling'] = ceil   # spelled ceiling in maxima


class Function_floor(PrimitiveFunction):
    """
    The floor function.

    EXAMPLES:
        sage: floor(5.4)
        5
        sage: type(floor(5.4))
        <type 'sage.rings.integer.Integer'>
        sage: a = floor(5.4 + x); a
        floor(x + 0.4000000000000004) + 5
        sage: a(2)
        7
    """
    def _repr_(self, simplify=True):
        return "floor"

    def _latex_(self):
        return "\\text{floor}"

    def _maxima_init_(self):
        return "floor"

    _approx_ = math.floor

    def __call__(self, x):
        try:
            return x.floor()
        except AttributeError:
            if isinstance(x, float):
                return math.floor(x)
        return SymbolicComposition(self, SR(x))

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


cos = Function_cos()
_syms['cos'] = cos

class Function_sec(PrimitiveFunction):
    """
    The secant function

    EXAMPLES:
        sage: sec(pi/4)
        sqrt(2)
        sage: RR(sec(pi/4))
        1.41421356237310
        sage: sec(1/2)
        sec(1/2)
        sage: sec(0.5)
        1.139493927324549
    """
    def _repr_(self, simplify=True):
        return "sec"

    def _latex_(self):
        return "\\sec"

    def _approx_(self, x):
        return 1/math.cos(x)

sec = Function_sec()
_syms['sec'] = sec

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

def atan2(y, x):
    return atan(y/x)

_syms['atan2'] = atan2

class Function_asin(PrimitiveFunction):
    """
    The arcsine function

    EXAMPLES:
        sage: asin(0.5)
        0.523598775598299
        sage: asin(1/2)
        pi/6
        sage: asin(1 + I*1.0)
        1.061275061905036*I + 0.6662394324925153
    """
    def _repr_(self, simplify=True):
        return "asin"

    def _latex_(self):
        return "\\sin^{-1}"

    def _approx_(self, x):
        return math.asin(x)

asin = Function_asin()
_syms['asin'] = asin

class Function_acos(PrimitiveFunction):
    """
    The arccosine function

    EXAMPLES:
        sage: acos(0.5)
        1.04719755119660
        sage: acos(1/2)
        pi/3
        sage: acos(1 + I*1.0)
        0.9045568943023813 - 1.061275061905036*I
    """
    def _repr_(self, simplify=True):
        return "acos"

    def _latex_(self):
        return "\\cos^{-1}"

    def _approx_(self, x):
        return math.acos(x)

acos = Function_acos()
_syms['acos'] = acos


class Function_atan(PrimitiveFunction):
    """
    The arctangent function.

    EXAMPLES:
        sage: atan(1/2)
        atan(1/2)
        sage: RDF(atan(1/2))
        0.463647609001
        sage: atan(1 + I)
        atan(I + 1)
    """
    def _repr_(self, simplify=True):
        return "atan"

    def _latex_(self):
        return "\\tan^{-1}"

    def _approx_(self, x):
        return math.atan(x)

atan = Function_atan()
_syms['atan'] = atan


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
        sage: float(cosh(pi))
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
        sage: float(sech(pi))
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
# log and exp
#############

class Function_log(PrimitiveFunction):
    """
    The log funtion. This is a symbolic logarithm.

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
ln = function_log

def log(x, base=None):
    if base is None:
        try:
            return x.log()
        except AttributeError:
            return function_log(x)
    else:
        try:
            return x.log(base)
        except AttributeError:
            return function_log(x) / function_log(base)

class Function_sqrt(PrimitiveFunction):
    """
    The square root function. This is a symbolic square root.

    EXAMPLES:
        sage: sqrt(-1)
        I
        sage: sqrt(2)
        sqrt(2)
        sage: sqrt(x^2)
        abs(x)
    """
    def __init__(self):
        PrimitiveFunction.__init__(self, needs_braces=True)

    def _repr_(self, simplify=True):
        return "sqrt"

    def _latex_(self):
        return "\\sqrt"

    def __call__(self, x):
        # if x is an integer or rational, never call the sqrt method
        if isinstance(x, float):
            return self._approx_(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.sqrt()
            except AttributeError:
                pass
        return SymbolicComposition(self, SR(x))


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
        sage: exp(float(2.5))
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
        sage: g = f(x,y,z)
        sage: g
        foo(x, y, z)
        sage: g(x=var('theta'))
        foo(theta, y, z)
    """
    def __init__(self, name):
        PrimitiveFunction.__init__(self, needs_braces=True)
        self._name = str(name)

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

    def __call__(self, *args, **kwds):
        return SymbolicFunctionEvaluation(self, [SR(x) for x in args])

class SymbolicFunction_delayed(SymbolicFunction):
    def simplify(self):
        return self

    def is_simplified(self):
        return True

    def _maxima_init_(self):
        return "%s"%self._name

    def __call__(self, *args):
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
    def __init__(self, f, args=None, kwds=None):
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
        self._kwds = kwds

    def _is_atomic(self):
        return True

    def arguments(self):
        return tuple(self._args)

    def keyword_arguments(self):
        return self._kwds

    def _repr_(self, simplify=True):
        if simplify:
            return self.simplify()._repr_(simplify=False)
        else:
            args = ', '.join([x._repr_(simplify=simplify) for x in
                                                      self._args])
            if not self._kwds is None:
                kwds = ', '.join(["%s=%s" %(x, y) for x,y in self._kwds.iteritems()])
                return '%s(%s, %s)' % (self._f._name, args, kwds)
            else:
                return '%s(%s)' % (self._f._name, args)

    def _latex_(self):
        return "{\\rm %s}(%s)"%(self._f._name, ', '.join([x._latex_() for
                                                       x in self._args]))

    def _maxima_init_(self):
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

    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
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
            del kwds[x]

        arg = tuple([SR(x._recursive_sub(kwds)) for x in self._args])

        if function_sub:
            return g(*arg)
        else:
            return self.__class__(self._f, arg)

    def _recursive_sub_over_ring(self, kwds, ring):
        raise TypeError, "no way to coerce the formal function to ring."

    def variables(self):
        """
        Return the variables appearing in the simplified form of self.

        EXAMPLES:
            sage: foo = function('foo')
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
    def simplify(self):
        return self

    def is_simplified(self):
        return True

    def __float__(self):
        return float(self._maxima_())

    def __complex__(self):
        return complex(self._maxima_())

    def _real_double_(self, R):
        return R(float(self))

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
        sage: f = function('cr', a)
        sage: g = f.diff(a).integral(b)
        sage: g
        diff(cr(a), a, 1)*b
        sage: g(cr=cos)
        -sin(a)*b
        sage: g(cr=sin(x) + cos(x))
        (cos(a) - sin(a))*b
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
        X = _functions[s]()
        if not X is None:
            return X
    except KeyError:
        pass
    v = SymbolicFunction(s)
    _functions[s] = weakref.ref(v)
    return v

#######################################################
def dummy_limit(*args):
    s = str(args[1])
    return SymbolicFunctionEvaluation(function('limit'), args=(args[0],),kwds={s: args[2]})
######################################i################




#######################################################

symtable = {'%pi':'pi', '%e': 'e', '%i':'I', '%gamma':'euler_gamma'}

from sage.rings.infinity import infinity, minus_infinity

_syms['inf'] = infinity
_syms['minf'] = minus_infinity

from sage.misc.multireplace import multiple_replace

maxima_tick = re.compile("'[a-z|A-Z|0-9|_]*")

maxima_qp = re.compile("\?\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

maxima_var = re.compile("\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

def symbolic_expression_from_maxima_string(x, equals_sub=False, maxima=maxima):
    global _syms

    if len(x) == 0:
        raise RuntimeError, "invalid symbolic expression -- ''"
    maxima.set('_tmp_',x)

    # This is inefficient since it so rarely is needed:
    #r = maxima._eval_line('listofvars(_tmp_);')[1:-1]

    s = maxima._eval_line('_tmp_;')

    formal_functions = maxima_tick.findall(s)
    if len(formal_functions) > 0:
        for X in formal_functions:
            _syms[X[1:]] = function(X[1:])
        # You might think there is a potential very subtle bug if 'foo is in a string literal --
        # but string literals should *never* ever be part of a symbolic expression.
        s = s.replace("'","")

    delayed_functions = maxima_qp.findall(s)
    if len(delayed_functions) > 0:
        for X in delayed_functions:
            _syms[X[2:]] = SymbolicFunction_delayed(X[2:])
        s = s.replace("?%","")

    s = multiple_replace(symtable, s)
    s = s.replace("%","")

    if equals_sub:
        s = s.replace('=','==')

    # have to do this here, otherwise maxima_tick catches it
    _syms['limit'] = dummy_limit

    global is_simplified
    try:
        # use a global flag so all expressions obtained via
        # evaluation of maxima code are assumed pre-simplified
        is_simplified = True
        last_msg = ''
        while True:
            try:
                w = sage_eval(s, _syms)
            except NameError, msg:
                if msg == last_msg:
                    raise NameError, msg
                msg = str(msg)
                last_msg = msg
                i = msg.find("'")
                j = msg.rfind("'")
                nm = msg[i+1:j]
                _syms[nm] = var(nm)
            else:
                break
        if isinstance(w, (list, tuple)):
            return w
        else:
            x = SR(w)
        return x
    except SyntaxError:
        raise TypeError, "unable to make sense of Maxima expression '%s' in SAGE"%s
    finally:
        is_simplified = False

def symbolic_expression_from_maxima_element(x):
    return symbolic_expression_from_maxima_string(x.name())

def evaled_symbolic_expression_from_maxima_string(x):
    return symbolic_expression_from_maxima_string(maxima.eval(x))
