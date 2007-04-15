r"""
The \sage calculus and elementary functions module.

AUTHORS:
    Bobby Moretti and William Stein: 2006--2007

The \sage calculus module is loosely based on the \sage Enhahcement Proposal
found at: http://www.sagemath.org:9001/CalculusSEP.

EXAMPLES:
    The basic units of the calculus suite are
    SymbolicExpressions. There are many subclasses of
    SymbolicExpression. The most basic of these is the formal
    indeterminate class, SymbolicVariable. To create a
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
        (z/2*sin(y*z/55) + y + x)^(z/2*sin(y*z/55) + y + x)

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
        Traceback (most recent call last):
        ...
        TypeError: must give explicit variable names to substitute for

    We can also make CallableFunctions, which is a SymbolicExpression
    that are functions of variables in a fixed order. Each
    SymbolicExpression has a function() method used to create a
    CallableFunction.

        sage: u = log((2-x)/(y+5))
        sage: f = u.function(x, y); f
        (x, y) |--> log((2 - x)/(y + 5))

    There is an easier way of creating a CallableFunction, which
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
    I*((exp(1) - exp(-1))/2) + (exp(1/2) + exp(-1/2))/2
    sage: CC(f)
    1.12762596520638 + 1.17520119364380*I
    sage: ComplexField(200)(f)
    1.1276259652063807852262251614026720125478471180986674836290 + 1.1752011936438014568823818505956008151557179813340958702296*I
    sage: ComplexField(100)(f)
    1.1276259652063807852262251614 + 1.1752011936438014568823818506*I

"""



from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_MPolynomial,
                            is_RealNumber, is_ComplexNumber, RR,
                            Integer, Rational, CC)

from sage.structure.element import RingElement
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex
from sage.structure.sage_object import SageObject

from sage.interfaces.maxima import MaximaElement, Maxima
from sage.interfaces.all import maxima

from sage.misc.sage_eval import sage_eval

from sage.functions.constants import Constant
import sage.functions.constants as c

from sage.calculus.equations import SymbolicEquation
from sage.rings.real_mpfr import RealNumber
from sage.rings.complex_number import ComplexNumber
from sage.rings.real_double import RealDoubleElement
from sage.rings.real_mpfi import RealIntervalFieldElement

from sage.libs.pari.gen import pari, gen as PariGen

from sage.rings.complex_double import ComplexDoubleElement
def is_DoublePrecisionNumber(x):
    if isinstance(x, (float, RealDoubleElement, ComplexDoubleElement)):
        return True
    if isinstance(x, RealNumber) and x.prec() == 53:
        return True
    if isinstance(x, ComplexNumber) and x.prec() == 53:
        return True
    return False

import math

# There will only ever be one instance of this class
class SymbolicExpressionRing_class(CommutativeRing):
    """
    The ring of all formal symbolic expressions.

    EXAMPLES:
        sage: SER
        Ring of Symbolic Expressions
        sage: type(SER)
        <class 'sage.calculus.calculus.SymbolicExpressionRing_class'>
    """
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        """
        Coerce x into the symbolic expression ring SER.

        EXAMPLES:
            sage: a = SER(-3/4); a
            -3/4
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: a.parent()
            Ring of Symbolic Expressions

        If a is already in the symblic expression ring, coercing returns
        a itself (not a copy):
            sage: SER(a) is a
            True
        """
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, CallableFunction):
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
                            Constant,
                            Integer,
                            int,
                            Rational,
                            PariGen,
                            ComplexNumber)):
            return SymbolicConstant(x)
        elif isinstance(x, Constant):
            return SymbolicConstant(x)
        else:
            raise TypeError, 'cannot coerce %s into a SymbolicExpression.'%x

    def _repr_(self):
        return  'Ring of Symbolic Expressions'

    def _latex_(self):
        return 'SymbolicExpressionRing'

    def characteristic(self):
        return Integer(0)

    def _an_element_impl(self):
        return zero_constant

    def is_field(self):
        return True

# Define the unique symbolic expression ring.
SER = SymbolicExpressionRing_class()

# The factory function that returns the unique SER.
def SymbolicExpressionRing():
    """
    Return the symbolic expression ring.

    EXAMPLES:
        sage: SymbolicExpressionRing()
        Ring of Symbolic Expressions
        sage: SymbolicExpressionRing() is SER
        True
    """
    return SER

class SymbolicExpression(RingElement):
    """
    A Symbolic Expression.

    EXAMPLES:

    """
    def __init__(self):
        RingElement.__init__(self, SER)

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
            sage: f.display2d()
                                                     3
                                                y + x
                                             ------------
                                                2
                                             3 y  + x + 1

        Use onscreen=False to get the 2d string:
             sage: print f.display2d(onscreen=False)
                                         3
                                    y + x
                                 ------------
                                    2
                                 3 y  + x + 1


        ASCII are really helps for the following integral:
            sage: f = integral(sin(x^2)); f
            sqrt(pi)*((sqrt(2)*I + sqrt(2))*erf((sqrt(2)*I + sqrt(2))*x/2) + (sqrt(2)*I - sqrt(2))*erf((sqrt(2)*I - sqrt(2))*x/2))/8
            sage: f.display2d()
                                                         (sqrt(2)  I + sqrt(2)) x
                   sqrt( Pi) ((sqrt(2)  I + sqrt(2)) erf(------------------------)
                                                                    2
                                                               (sqrt(2)  I - sqrt(2)) x
                                  + (sqrt(2)  I - sqrt(2)) erf(------------------------))/8
                                                                          2
        """
        s = self._maxima_().display2d(onscreen=False)
        s = s.replace('%pi',' Pi').replace('%i',' I').replace('%e', ' E')
        if onscreen:
            print s
        else:
            return s

    def is_simplified(self):
        return hasattr(self, '_simp')

    def hash(self):
        return hash(self._repr_(simplify=False))

    def plot(self, **kwds):
        from sage.plot.plot import plot
        # see if the user passed a param
        try:
            param = kwds['param']
        except KeyError:
            if isinstance(self.simplify(), SymbolicConstant):
                return plot(self.simplify()._obj)

            vars = list(self._get_vars())
            if len(vars) == 1:
                return plot(self.function(vars[0]), **kwds)

            raise TypeError, "must give an explicit parameter to plot this"\
                            + " expression."

        del kwds['param']
        return plot(self.function(param), **kwds)


    def __eq__(self, right):
        return SymbolicEquation(self, right)

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

    def _maxima_init_(self):
        return self._repr_(simplify=False)

    def _mpfr_(self, field):
        raise TypeError

    def _complex_mpfr_field_(self, field):
        raise TypeError

    def _complex_double_(self, C):
        raise TypeError

    def _real_double_(self, R):
        raise TypeError

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

    def _get_vars(self, vars=None):
        """
        Returns the set of variables appearing in self.

        OUTPUT:
            a Python set

        EXAMPLES:
            sage: f = x^(n+1) + sin(pi/19); f
            x^(n + 1) + sin(pi/19)
            sage: f._get_vars()
            set([x, n])
        """
        if vars is None:
            vars = set([])
        return vars

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

    def function(self, *args):
        """
        Return a CallableFunction, fixing a variable order to be the order of
        args.

        EXAMPLES:
           sage: u = sin(x) + x*cos(y)
           sage: g = u.function(x,y)
           sage: g(x,y)
           x*cos(y) + sin(x)
           sage: g(t,z)
           t*cos(z) + sin(t)
           sage: g(x^2, x^y)
           x^2*cos(x^y) + sin(x^2)
        """
        return CallableFunction(self, args)


    ###################################################################
    # derivative
    ###################################################################
    def derivative(self, *args):
        """
        EXAMPLES:
            sage: h = sin(x)/cos(x)
            sage: diff(h,x,x,x)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2
            sage: diff(h,x,3)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2

            sage: u = (sin(x) + cos(y))*(cos(x) - sin(y))
            sage: diff(u,x,y)
            sin(x)*sin(y) - cos(x)*cos(y)
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
            vars = list(self._get_vars())
            if len(vars) == 1:
                s = "%s, %s" % (vars[0], a)
            else:
                raise ValueError, "must supply an explicit variable for an " +\
                                "expression containing more than one variable"
        for i in range(len(args)):
            if isinstance(args[i], SymbolicVariable):
                s = s + '%s, ' %str(args[i])
                # check to see if this is followed by an integer
                try:
                    if isinstance(args[i+1], (int, long, Integer)):
                        s = s + '%s, ' %str(args[i+1])
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


    ###################################################################
    # integral
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

            sage: f = x*cos(x^2)
            sage: f.integral(x, 0, sqrt(pi))
            0
            sage: f.integral(a=-pi, b=pi)
            0

        We integrate a spherical Bessel function:
            sage: f = spherical_bessel_Y(2,x); f
            (-3*sin(x)/x - (1 - (3/x^2))*cos(x))/x
            sage: f.integral(x)
            ???
        """
        if v is None:
            vars = self._get_vars()
            if len(vars) == 1:
                v = list(vars)[0]
        if not isinstance(v, SymbolicVariable):
            v = var(str(v))
            #raise TypeError, 'must integrate with respect to a variable'
        if (a is None and (not b is None)) or (b is None and (not a is None)):
            raise TypeError, 'only one endpoint given'
        if a is None:
            return self.parent()(self._maxima_().integrate(v))
        else:
            return self.parent()(self._maxima_().integrate(v, a, b))


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
        Employs the identities $\sin(x)^2 + \cos(x)^2 = 1$ and
        $\cosh(x)^2 - \sin(x)^2 = 1$ to simplify expressions
        containing tan, sec, etc., to sin, cos, sinh, cosh.

        ALIAS: trig_simplify and simplify_trig are the same

        EXAMPLES:
            sage: f = sin(x)^2 + cos(x)^2; f
            sin(x)^2 + cos(x)^2
            sage: f.simplify()
            sin(x)^2 + cos(x)^2
            sage: f.simplify_trig()
            1
        """
        return self.parent()(self._maxima_().trigsimp())

    trig_simplify = simplify_trig

    def simplify_rational(self):
        """

        ALIAS: rational_simplify and simplify_rational are the same
        """
        return self.parent()(self._maxima_().fullratsimp())

    rational_simplify = simplify_rational

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
            (-y - x)*(y^2 + x*y + x^2)
            sage: factor(-8*y - 4*x + z^2*(2*y + x))
            (2*y + x)*(z - 2)*(z + 2)
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: F = factor(f/(36*(1 + 2*y + y^2)), dontfactor=[x])
            sage: F.display2d()
                                          2
                                        (x  + 2 x + 1) (y - 1)
                                        ----------------------
                                              36 (y + 1)
        """
        if len(dontfactor) > 0:
            m = self._maxima_()
            name = m.name()
            cmd = 'block ([dontfactor: %s], factor (%s))'%(dontfactor,
                                               name)
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
        try:
            in_dict = SER(in_dict)
        except TypeError:
            pass

        if isinstance(in_dict, SymbolicExpression):
            vars = list(self._get_vars())
            if len(vars) == 1:
                in_dict = {vars[0]: in_dict}

            elif not ((isinstance(in_dict, dict) or in_dict is None)):
               raise TypeError, "must give explicit variable names to substitute for"

        if in_dict is not None:
            for k, v in in_dict.iteritems():
               kwds[str(k)] = v
        # find the keys from the keywords
        return self._recursive_sub(kwds)

    def _recursive_sub(self, kwds):
        raise NotImplementedError, "implement _recursive_sub for type '%s'!"%(type(self))

class CallableFunction(RingElement):
    r"""
    A callable symbolic function that knows the variables on which it depends.
    """
    def __init__(self, expr, args):
        RingElement.__init__(self, CallableFunctionRing)
        if args == [] or args == () or args is None:
            raise ValueError, "A CallableFunction must know at least one of" \
                               +" its variables."
        for arg in args:
            if not isinstance(arg, SymbolicExpression):
                raise TypeError, "Must construct a function with a list of" \
                                +" SymbolicVariables."

            self._varlist = args
            self._expr = expr

        def __float__(self):
            return float(self.obj)

        def _mpfr_(self, field):
            """
            Coerce to a multiprecision real number.

            EXAMPLES:
                sage: RealField(100)(SER(10))
                10.000000000000000000000000000
            """
            return (self.obj)._mpfr_(field)

    # TODO: should len(args) == len(vars)?
    def __call__(self, *args):
        vars = self._varlist

        dct = {}
        for i in range(len(args)):
            dct[vars[i]] = args[i]

        return self._expr.substitute(dct)

    def _repr_(self, simplify=True):
        if len(self._varlist) == 1:
            return "%s |--> %s" % (self._varlist[0], self._expr._repr_(simplify=simplify))
        else:
            vars = ", ".join(map(str, self._varlist))
            return "(%s) |--> %s" % (vars, self._expr._repr_(simplify=simplify))

    def _latex_(self):
        if len(self._varlist) == 1:
            return "%s \\ {\mapsto}\\ %s" % (self._varlist[0],
                    self._expr._latex_())
        else:
            vars = ", ".join(map(latex, self._varlist))
            # the weird TeX is to workaround an apparent JsMath bug
            return "\\left(%s \\right)\\ {\\mapsto}\\ %s" % (vars, self._expr._latex_())

    def derivative(self, *args):
        """
        Returns the derivative of this symbolic function with respect to the
        given variables. If the function has only one variable, it returns the
        derivative with respect to that variable.
        """
        # if there are no args and only one var, differentiate wrt that var
        if args == () and len(self._varlist) == 1:
            return CallableFunction(self._expr.derivative(self._varlist[0]), self._varlist)
        # if there's only one var and the arg is an integer n, diff. n times wrt
        # that var
        elif isinstance(args[0], Integer) and len(self._varlist) == 1:
            f = self._expr.derivative(self._varlist[0], args[0])
            return CallableFunction(f, self._varlist)
        # else just take the derivative
        else:
            return CallableFunction(self._expr.derivative(*args), self._varlist)

    def integral(self, dx):
        return CallableFunction(self._expr.integral(dx), self._varlist)

    def substitute(self, dict=None, **kwds):
        return CallableFunction(self._expr.substitute(dict, kwds), self._varlist)

    def simplify(self):
        return CallableFunction(self._expr.simplify(), self._varlist)

    def trig_simplify(self):
        return CallableFunction(self._expr.trig_simplify(), self._varlist)

    simplify_trig = trig_simplify

    def rational_simplify(self):
        return CallableFunction(self._expr.rational_simplify(), self._varlist)

    simplify_rational = rational_simplify

    def expand(self):
        return CallableFunction(self._expr.expand(), self._varlist)

    def trig_expand(self):
        return CallableFunction(self._expr.trig_expand(), self._varlist)

    def _neg_(self):
        result = SymbolicArithmetic([self._expr], operator.neg)
        return CallableFunction(result, self._unify_varlists(right))

    def _add_(self, right):
        result = SymbolicArithmetic([self._expr, right._expr], operator.add)
        return CallableFunction(result, self._unify_varlists(right))

    def _sub_(self, right):
        result = SymbolicArithmetic([self._expr, right._expr], operator.sub)
        return CallableFunction(result, self._unify_varlists(right))

    def _mul_(self, right):
        result = SymbolicArithmetic([self._expr, right._expr], operator.mul)
        return CallableFunction(result, self._unify_varlists(right))

    def _div_(self, right):
        result = SymbolicArithmetic([self._expr, right._expr], operator.div)
        return CallableFunction(result, self._unify_varlists(right))

    def __pow__(self, right):
        try:
            result = SymbolicArithmetic([self._expr, right._expr], operator.pow)
        except AttributeError:
            result = SymbolicArithmetic([self._expr, right], operator.pow)
        return CallableFunction(result, self._varlist)


    def _unify_varlists(self, x):
        r"""
        Takes the variable list from another CallableFunction object and
        compares it with the current CallableFunction object's variable list,
        combining them according to the following rules:

        Let \code{a} be \code{self}'s variable list, let \code{b} be
        \code{y}'s variable list.

        \begin{enumerate}
        \item If \code{a == b}, then the variable lists are identical, so return
        that variable list.
        \item If \code{a} $\neq$ \code{b}, then check if the first $n$ items in
        \code{a} are the first $n$ items in \code{b}, or vice-versa. If so,
        return a list with these $n$ items, followed by the remaining items in
        $\code{b}$ and $\code{a}$.
        \item If none of the previous conditions are met, then just return
        \code{a} $\cup$ \code{b}, sorted alphabetically.

        Note: When used for arithmetic between CallableFunctions, these rules
        ensure that the set of CallableFunctions will have certain properties.
        In particular, it ensures that the set is a commutative ring.

        INPUT:
            x -- A CallableFunction

        OUTPUT:
            A combination of \code{self}'s variables and \code{x}'s variables
            sorted in a nice way.

        EXAMPLES:
            sage: f(x, y, z) = sin(x+y+z)
            sage: f
            (x, y, z) |--> sin(z + y + x)
            sage: g(x, y) = y + 2*x
            sage: g
            (x, y) |--> y + 2*x
            sage: f._unify_varlists(g)
            (x, y, z)
            sage: g._unify_varlists(f)
            (x, y, z)

            sage: f(x, y, z) = sin(x+y+z)
            sage: g(w, t) = cos(w - t)
            sage: g
            (w, t) |--> cos(w - t)
            sage: f._unify_varlists(g)
            (t, w, x, y, z)

            sage: f(x, y, t) = y*(x^2-t)
            sage: f
            (x, y, t) |--> (x^2 - t)*y
            sage: g(x, y, w) = x + y - cos(w)
            sage: f._unify_varlists(g)
            (x, y, t, w)
            sage: g._unify_varlists(f)
            (x, y, t, w)
            sage: f*g
            (x, y, t, w) |--> (x^2 - t)*y*(y + x - cos(w))

            sage: f(x,y, t) = x+y
            sage: g(x, y, w) = w + t
            sage: f._unify_varlists(g)
            (x, y, t, w)
            sage: g._unify_varlists(f)
            (x, y, t, w)
            sage: f + g
            (x, y, t, w) |--> y + x + w + t

        AUTHORS:
            -- Bobby Moretti, thanks to William Stein for the rules

        """
        a = self._varlist
        b = x._varlist

        # Rule #1
        if a == b:
            return self._varlist

        # Rule #2
        new_list = []
        done = False
        i = 0
        while not done and i < min(len(a), len(b)):
            if a[i] == b[i]:
                new_list.append(a[i])
                i += 1
            else:
                done = True

        temp = []
        # Rule #3, plus the rest of Rule #4
        for j in range(i, len(a)):
            if not a[j] in temp:
                temp.append(a[j])

        for j in range(i, len(b)):
            if not b[j] in temp:
                temp.append(b[j])

        temp.sort()
        new_list.extend(temp)
        return tuple(new_list)


class Symbolic_object(SymbolicExpression):
    r"""
    A class representing a symbolic expression in terms of a SageObject (not
    necessarily a 'constant')
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
        return str(self._obj)

    def _latex_(self):
        """
        EXAMPLES:
        """
        return latex(self._obj)

    def __pow__(self, right):
        """
        EXAMPLES:
        """
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj ** right._obj)
            except TypeError:
                pass
        return SymbolicExpression.__pow__(self, right)

    def str(self, bits=None):
        if bits is None:
            return str(self._obj)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self._obj))

    def _maxima_init_(self):
        return maxima_init(self._obj)

def maxima_init(x):
    try:
        return x._maxima_init_()
    except AttributeError:
        return str(x)

class SymbolicConstant(Symbolic_object):
    def __init__(self, x):
        Symbolic_object.__init__(self, x)

    def _is_atomic(self):
        try:
            return self._atomic
        except AttributeError:
            if isinstance(x, Rational):
                self._atomic = False
            else:
                self._atomic = True
            return self._atomic

    def _recursive_sub(self, kwds):
        """
        EXAMPLES:
            sage: a = SER(5/6)
            sage: type(a)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: a(x=3)
            5/6
        """
        return self



class SymbolicPolynomial(Symbolic_object):
    """
    An element of a polynomial ring as a formal symbolic expression.

    EXAMPLES:
    A single variate polynomial:
        sage: R.<x> = QQ[]
        sage: f = SER(x^3 + x)
        sage: f(y=7)
        x^3 + x
        sage: f(x=5)
        130
        sage: f.integral(x)
        x^4/4 + x^2/2
        sage: f(x=y)
        y*(y^2 + 1)

    A multivariate polynomial:

        sage: R.<x,y,theta> = ZZ[]
        sage: f = SER(x^3 + x + y + theta^2); f
        theta^2 + y + x + x^3
        sage: f(x=y, theta=y)
        2*y + y^2 + y^3
        sage: f(x=5)
        130 + theta^2 + y

    The polynomial must be over a field of characteristic 0.
        sage: R.<w> = GF(7)[]
        sage: f = SER(w^3 + 1)
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
        f = self._obj
        if is_Polynomial(f):
            # Single variable case
            v = f.parent().variable_name()
            if kwds.has_key(v):
                return SER(f(kwds[v]))
            else:
                return self
        else:
            # Multivariable case
            t = []
            for g in f.parent().gens():
                s = str(g)
                if kwds.has_key(s):
                    t.append(kwds[s])
                else:
                    t.append(g)
            return SER(f(*t))

##################################################################


zero_constant = SymbolicConstant(Integer(0))

class SymbolicOperation(SymbolicExpression):
    r"""
    A parent class representing any operation on SymbolicExpression objects.
    """
    def __init__(self, operands):
        SymbolicExpression.__init__(self)
        self._operands = operands   # don't even make a copy -- ok, since immutable.

    def _get_vars(self, vars=None):
        if vars is None:
            vars = set([])
        for op in self._operands:
            opvars = op._get_vars(vars)
            for v in opvars:
                vars.add(v)
        return vars

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
            sage: ???
        """
        ops = self._operands
        new_ops = [op._recursive_sub(kwds) for op in ops]
        arith = SymbolicArithmetic(new_ops, self._operator)
        return  arith._operator(*(arith._operands))

    def __float__(self):
        fops = [float(op) for op in self._operands]
        return self._operator(*fops)

    def _mpfr_(self, field):
        rops = [op._mpfr_(field) for op in self._operands]
        return self._operator(*rops)

    def _complex_mpfr_field_(self, C):
        rops = [op._complex_mpfr_field_(C) for op in self._operands]
        return self._operator(*rops)

    def _complex_double_(self, C):
        rops = [op._complex_double_(C) for op in self._operands]
        return self._operator(*rops)

    def _real_double_(self, R):
        rops = [op._real_double_(R) for op in self._operands]
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
        if op in [operator.mul, operator.div, operator.sub]:
                # avoid drawing parens if s1 an atomic operation
                if not ops[1]._is_atomic():
                    s[1] = '(%s)' % s[1]

        elif op is operator.pow:
            if not ops[0]._is_atomic():
                s[0] = '(%s)'% s[0]
            if not ops[1]._is_atomic():
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
        infixops = {operator.add: '+',
                    operator.sub: '-',
                    operator.mul: '*',
                    operator.div: '/',
                    operator.pow: '^'}

        if self._operator is operator.neg:
            return '-(%s)' % ops[0]._maxima_init_()
        else:
            return '(%s) %s (%s)' % (ops[0]._maxima_init_(),
                                 infixops[self._operator],
                                 ops[1]._maxima_init_())

import re

class SymbolicVariable(SymbolicExpression):
    def __init__(self, name):
        SymbolicExpression.__init__(self)
        self._name = name
        if len(name) == 0:
            raise ValueError, "variable name must be nonempty"

    def _recursive_sub(self, kwds):
        s = str(self)
        # do the replacement if needed
        if s in kwds:
            return kwds[s]
        else:
            return self

    def _get_vars(self, vars=None):
        if vars is None:
            vars = set([])
        vars.add(self)
        return vars

    def __cmp__(self, right):
        return cmp(self._repr_(), right._repr_())

    def __call__(self, x):
        raise AttributeError, "A symbolic variable is not callable."

    def _repr_(self, simplify=True):
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
def var(s):
    r"""
    Create a symbolic variable with the name \emph{s}.
    """
    if ',' in s:
        return tuple([var(x.strip()) for x in s.split(',')])
    elif ' ' in s:
        return tuple([var(x.strip()) for x in s.split()])
    try:
        return _vars[s]
    except KeyError:
        v = SymbolicVariable(s)
        _vars[s] = v
        return v

class CallableFunctionRing_class(CommutativeRing):
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, (Integer, Rational, Constant, int)):
            return CallableFunction(x, SER(x))
        else:
            raise NotImplementedError, "cannot coerce this (yet)"

CallableFunctionRing = CallableFunctionRing_class()
CFR = CallableFunctionRing

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
        """
        EXAMPLES:
            sage: ??
        """
        ops = self._operands
        return ops[0](ops[1]._recursive_sub(kwds))

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
        if isinstance(ops[0], Function_log):
            base = ops[0]._base
            if not base is None:
                return 'log(%s) / log(%s)' % (ops[1]._maxima_init_(),str(base))
        return '%s(%s)' % (ops[0]._maxima_init_(), ops[1]._maxima_init_())

    def __float__(self):
        f = self._operands[0]
        g = self._operands[1]
        #return float(f(float(g)))
        return f.__call__(float(g))

    def _mpfr_(self, field):
        """
        Coerce to a multiprecision real number.

        EXAMPLES:
            sage: RealField(100)(sin(2)+cos(2))
            0.49315059027853930839845163641
        """
        f = self._operands[0]
        g = self._operands[1]
        return f(g._mpfr_(field))

    def _complex_mpfr_field_(self, field):
        """
        Coerce to a multiprecision complex number.

        EXAMPLES:
            sage: ComplexField(100)(sin(2)+cos(2)+I)
            0.49315059027853930839845163641 + 1.0000000000000000000000000000*I
        """
        f = self._operands[0]
        g = self._operands[1]
        return f(g._complex_mpfr_field_(field))

    def _complex_double_(self, field):
        """
        Coerce to a complex double.

        EXAMPLES:
            sage: CDF(sin(2)+cos(2)+I)
            0.493150590279 + 1.0*I
        """
        f = self._operands[0]
        g = self._operands[1]
        return f(g._complex_double_(field))

    def _real_double_(self, field):
        """
        Coerce to a real double.

        EXAMPLES:
            sage: RDF(sin(2)+cos(2))
            0.493150590279
        """
        f = self._operands[0]
        g = self._operands[1]
        return f(g._real_double_(field))


class PrimitiveFunction(SymbolicExpression):
    def __init__(self, needs_braces=False):
        self._tex_needs_braces = needs_braces
        pass

    def tex_needs_braces(self):
        return self._tex_needs_braces

    def __call__(self, x, *args):
        if isinstance(x, float):
            return self._approx_(x)
        try:
            return getattr(x, self._repr_())(*args)
        except AttributeError:
            return SymbolicComposition(self, SER(x))

    def _approx_(self, x):
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

    def _is_atomic(self):
        return True

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

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        return x.__abs__()

    def __call__(self, x):
        try:
            return x.__abs__()
        except (AttributeError, TypeError):
            return SymbolicComposition(self, SER(x))

abs = Function_abs()
_syms['abs'] = abs

class Function_sin(PrimitiveFunction):
    """
    The sine function

    EXAMPLES:
        sage: ???

    """
    def _repr_(self, simplify=True):
        return "sin"

    def _latex_(self):
        return "\\sin"

    def _is_atomic(self):
        return True

    _approx_ = math.sin

    def __call__(self, x):
        try:
            return x.sin()
        except AttributeError:
            if isinstance(x, float):
                return math.sin(x)
        return SymbolicComposition(self, SER(x))

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

    def _is_atomic(self):
        return True

    _approx_ = math.cos

    def __call__(self, x):
        try:
            return x.cos()
        except AttributeError:
            if isinstance(x, float):
                return math.cos(x)
        return SymbolicComposition(self, SER(x))


cos = Function_cos()
_syms['cos'] = cos

class Function_sec(PrimitiveFunction):
    """
    The secant function

    EXAMPLES:
        sage: sec(pi/4)
        sqrt(2)
        sage: RR(sec(pi/4))
        1.414213562373094
        sage: sec(1/2)
        sec(1/2)
        sage: sec(0.5)
        1.139493927324549
    """
    def _repr_(self, simplify=True):
        return "sec"

    def _latex_(self):
        return "\\sec"

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        try:
            return x.cos()
        except AttributeError:
            if isinstance(x, float):
                return float(1)/float(cos(x))
            else:
                return SymbolicComposition(self, SER(x))

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

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        try:
            return x.tan()
        except AttributeError:
            if isinstance(x, float):
                return math.tan(x)
            else:
                return SymbolicComposition(self, SER(x))

tan = Function_tan()
_syms['tan'] = tan

class Function_asin(PrimitiveFunction):
    """
    The arcsine function

    EXAMPLES:
        sage: asin(0.5)
        0.523598775598299
        sage: asin(1/2)
        (pi/6)
        sage: asin(1 + I*1.0)
        0.666239432492515 + 1.06127506190504*I

    """
    def _repr_(self, simplify=True):
        return "asin"

    def _latex_(self):
        return "\\sin^{-1}"

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        try:
            return x.asin()
        except AttributeError:
            if isinstance(x, float):
                return math.asin(x)
            else:
                return SymbolicComposition(self, SER(x))

asin = Function_asin()
_syms['asin'] = asin

class Function_acos(PrimitiveFunction):
    """
    The arccosine function

    EXAMPLES:
        sage: acos(0.5)
        1.04719755119660
        sage: acos(1/2)
        (pi/3)
        sage: acos(1 + I*1.0)
        0.904556894302381 - 1.06127506190504*I

    """
    def _repr_(self, simplify=True):
        return "acos"

    def _latex_(self):
        return "\\cos^{-1}"

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        try:
            return x.acos()
        except AttributeError:
            if isinstance(x, float):
                return math.acos(x)
            else:
                return SymbolicComposition(self, SER(x))

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
        1.01722196789785 + 0.402359478108525*I
    """
    def _repr_(self, simplify=True):
        return "atan"

    def _latex_(self):
        return "\\tan^{-1}"

    def _is_atomic(self):
        return True

    def _approx_(self, x):
        try:
            return x.atan()
        except AttributeError:
            if isinstance(x, float):
                return math.atan(x)
            else:
                return SymbolicComposition(self, SER(x))

atan = Function_atan()
_syms['atan'] = atan


class Function_log(PrimitiveFunction):
    """
    The log funtion.

    EXAMPLES:
        sage: log(e^2)
        2
        sage: log(1024, 2) # the following is ugly (for now)
        log(1024)/log(2)
        sage: log(10, 4)
        log(10)/log(4)
        sage: log(1024.0, 2.0)
        10.0000000000000

        sage: RDF(log(10,2))
        3.32192809489
        sage: RDF(log(8, 2))
        3.0
        sage: log(RDF(10))
        2.30258509299
        sage: log(2.718)
        0.999896315728952
    """
    def __init__(self, base=None):
        self._base = base
        PrimitiveFunction.__init__(self)

    def _repr_(self, simplify=True):
        if self._base is None:
            return "log"
        else:
            if self._base._is_atomic():
                return "log_%s" % self._base
            else:
                return "log_(%s)" % self._base

    def _latex_(self):
        if self._base is None:
            return "\\log"
        else:
            return "\\log_{%s}" % self._base

    def _is_atomic(self):
        return True

    def __call__(self, x, *args):
        # see if we have a base that's None
        base = self._base
        if base is None:
            # first try getting the base from what was passed
            try:
                base = args[0]
            except IndexError:
            # if that fails, get the base from the object
                base = self._base

        try:
            if base is None:
                return x.log()
            else:
                return x.log(base)
        except AttributeError:
            if isinstance(x, float):
                return self._approx_(x, self._base)

        # if the base is None, we behave as before
        if base is None:
            return SymbolicComposition(self, SER(x))
        # else construct a new log function with the correct base
        else:
            return SymbolicComposition(Function_log(base), SER(x))

    _approx_ = math.log

log = Function_log()
_syms['log'] = log

class Function_sqrt(PrimitiveFunction):
    """
    The square root function.

    EXAMPLES:
        sage: sqrt(-1)
        1.00000000000000*I
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

    def _is_atomic(self):
        return True

    def __call__(self, x):
        # if x is an integer or rational, never call the sqrt method
        if isinstance(x, float):
            return self._approx_(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.sqrt()
            except AttributeError:
                pass
        return SymbolicComposition(self, SER(x))


    def _approx_(self, x):
        try:
            return x._obj.sqrt()
        except AttributeError:
            try:
                return x.sqrt()
            except AttributeError:
                if isinstance(x, float):
                    return math.sqrt(x)
                else:
                    return SymbolicComposition(self, SER(x))

sqrt = Function_sqrt()
_syms['sqrt'] = sqrt

class Function_exp(PrimitiveFunction):
    """
    The square root function.

    EXAMPLES:
        sage: exp(-1)
        ?
        sage: exp(2)
        ?
        sage: exp(x^2 + log(x))
        ?
        sage: exp(2.5)
        sage: exp(float(2.5))
        sage: exp(RDF('2.5'))

    """
    def __init__(self):
        PrimitiveFunction.__init__(self, needs_braces=True)

    def _repr_(self, simplify=True):
        return "exp"

    def _latex_(self):
        return "\\exp"

    def _is_atomic(self):
        return True

    def __call__(self, x):
        # if x is an integer or rational, never call the sqrt method
        if isinstance(x, float):
            return self._approx_(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.exp()
            except AttributeError:
                pass
        return SymbolicComposition(self, SER(x))


    def _approx_(self, x):
        try:
            return x._obj.exp()
        except AttributeError:
            try:
                return x.exp()
            except AttributeError:
                if isinstance(x, float):
                    return math.exp(x)
                else:
                    return SymbolicComposition(self, SER(x))

exp = Function_exp()
_syms['exp'] = exp

#######################################################
symtable = {'%pi':'_Pi_', '%e': '_E_', '%i':'_I_'}
import sage.functions.constants as c
_syms['_Pi_'] = SER(c.pi)
_syms['_E_'] = SER(c.e)
_syms['_I_'] = SER(c.I)

def symbolic_expression_from_maxima_string(x):
    global _syms
    maxima.eval('listdummyvars: false')
    maxima.eval('_tmp_: %s'%x)
    r = maxima.eval('listofvars(_tmp_)')[1:-1]
    if len(r) > 0:
        # Now r is a list of all the indeterminate variables that
        # appear in the expression x.
        v = r.split(',')
        for a in v:
            _syms[a] = var(a)
    s = maxima.eval('_tmp_')
    for z, w in symtable.iteritems():
        s = s.replace(z, w)
    try:
        return SER(sage_eval(s, _syms))
    except SyntaxError:
        raise TypeError, "unable to make sense of Maxima expression '%s' in SAGE"%s

def symbolic_expression_from_maxima_element(x):
    return symbolic_expression_from_maxima_string(x.name())

def evaled_symbolic_expression_from_maxima_string(x):
    return symbolic_expression_from_maxima_string(maxima.eval(x))

a = var('a')
b = var('b')
c = var('c')
d = var('d')
e = var('e')
f = var('f')
g = var('g')
h = var('h')
i = var('i')
j = var('j')
k = var('k')
l = var('l')
m = var('m')
n = var('n')
o = var('o')
p = var('p')
q = var('q')
r = var('r')
s = var('s')
t = var('t')
u = var('u')
v = var('v')
w = var('w')
x = var('x')
y = var('y')
z = var('z')
A = var('A')
B = var('B')
C = var('C')
D = var('D')
E = var('E')
F = var('F')
G = var('G')
H = var('H')
J = var('J')
K = var('K')
L = var('L')
M = var('M')
N = var('N')
O = var('O')
P = var('P')
Q = var('Q')
R = var('R')
S = var('S')
T = var('T')
U = var('U')
V = var('V')
W = var('W')
X = var('X')
Y = var('Y')
Z = var('Z')

