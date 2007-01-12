from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_RealNumber, is_ComplexNumber, RR,
                            Integer, Rational)
#import sage.rings.rational
from sage.structure.element import RingElement
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex
from sage.structure.sage_object import SageObject

from sage.interfaces.maxima import MaximaElement, Maxima
from sage.interfaces.all import maxima

from sage.misc.sage_eval import sage_eval

from sage.functions.constants import Constant

import pdb

class PrimitiveFunction(SageObject):
    def __init__(self):
        # nothing so far
        pass

    def __call__(self, x):
        # if we're calling with a symbolic expression, do function composition
        if isinstance(x, SymbolicExpression) and not isinstance(x, Constant):
            return SymbolicComposition(self, x)

        # if x is a polynomial object, first turn it into a function and
        # compose self with x
        elif is_Polynomial(x):
            return Function_composition(self, SymbolicPolynomial(x))

        # if we can't figure out what x is, return a symbolic version of x
        elif isinstance(x, (Integer, Rational,
                            int, long, float, complex)):
            return self(Constant_object(x))

        else:
            return self(Symbolic_object(x))

# There will only ever be one instance of this class
class SymbolicExpressionRing_class(CommutativeRing):
    '''
    A Ring of formal expressions.
    '''
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, SymbolicExpression):
            return x
        elif isinstance(x, MaximaElement):
            return symbolic_expression_from_maxima_element(x)
        elif is_Polynomial(x):
            return SymbolicPolynomial(x)
        elif isinstance(x, Integer):
            return Constant_object(x)
        else:
            return Symbolic_object(x)

    def _repr_(self):
        return  "Ring of Symbolic Expressions"

    def _latex_(self):
        return "SymbolicExpressionRing"

    def characteristic(self):
        return Integer(0)

    def _an_element_impl(self):
        return zero_constant

# ... and here it is:
SymbolicExpressionRing = SymbolicExpressionRing_class()
SER = SymbolicExpressionRing

class SymbolicExpression(RingElement):
    r"""
    A Symbolic Expression in acoordance with SEP #1.

    """
    # conversions is the dict of the form system:command
    def __init__(self, conversions={}):
        RingElement.__init__(self, SymbolicExpressionRing)

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            m = maxima(self._maxima_init_())
            self.__maxima = m
            return m

    def _neg_(self):
        return -1*self

    def __cmp__(self, right):
        return cmp(maxima(self), maxima(right))

    def _add_(self, right):
        return SymbolicArithmetic(self, right, operator.add)

    def _sub_(self, right):
        return SymbolicArithmetic(self, right, operator.sub)

    def _mul_(self, right):
        return SymbolicArithmetic(self, right, operator.mul)

    def _div_(self, right):
        return SymbolicArithmetic(self, right, operator.div)

    def __pow__(self, right):
        right = self.parent()(right)
        return SymbolicArithmetic(self, right, operator.pow)

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
                    raise ValueError, "can not take negative derivative"
            else:
                raise TypeError, "arguments must be integers or \
                SymbolicVariable objects"

        try:
            if s[-2] == ',':
                s = s[:-2]
        except IndexError:
            pass
        t = maxima('diff(%s, %s)'%(self._maxima_().name(), s))
        f = self.parent()(t)
        return f



    ###################################################################
    # simplify
    ###################################################################
    def simplify(self):
        return self.parent()(self._maxima_())

    def simplify_trig(self):
        r"""
        Employs the identities $\sin(x)^2 + \cos(x)^2 = 1$ and
        $\cosh(x)^2 - \sin(x)^2 = 1$ to simplify expressions
        containing tan, sec, etc., to sin, cos, sinh, cosh.

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

        ALIAS: trig_expand
        """
        return self.parent()(self._maxima_().trigexpand())

    trig_expand = expand_trig




class Symbolic_object(SymbolicExpression):
    r'''
    A class representing a symbolic expression in terms of a SageObject.
    '''

    def __init__(self, obj):
        SymbolicExpression.__init__(self)
        self._obj = obj

    def obj(self):
        return self._obj

    def _repr_(self):
        return str(self._obj)

    def _latex_(self):
        return latex(self._obj)

    # TODO: do _pow_?
    def _add_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj + right._obj)
            except TypeError:
                pass
        return SymbolicArithmetic(self, right, operator.add)

    def _sub_(self, right):
        if isinstnace(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj - right._obj)
            except TypeError:
                pass
        return SymbolicArithmetic(self, right, operator.sub)

    def _mul_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj * right._obj)
            except TypeError:
                pass
        return SymbolicArithmetic(self, right, operator.mul)

    def _div_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj / right._obj)
            except TypeError:
                pass
        return SymbolicArithmetic(self, right, operator.div)

    def str(self, bits=None):
        if bits is None:
            return str(self._obj)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self._obj))

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            m = maxima(self._obj)
            self.__maxima = m
            return m

class SymbolicPolynomial(Symbolic_object):
    "An element of a polynomial ring as a formal symbolic expression."

    # for now we do nothing except pass the info on to the supercontructor. It's
    # not clear to me why we need anything else in this class -Bobby
    def __init__(self, p):
       Symbolic_object.__init__(self, p)

class SymbolicConstant(SymbolicExpression):
    def __call__(self, x):
        return self

class Constant_object(SymbolicConstant, Symbolic_object):
    pass

zero_constant = Constant_object(Integer(0))

class SymbolicOperation(SymbolicExpression):
    r"""
    A parent class representing any operation on SymbolicExpression objects.
    """
    def __init__(self, operands):
        SymbolicExpression.__init__(self)
        self._operands = operands

class SymbolicComposition(SymbolicOperation):
    r'''
    Represents the symbolic composition of $f \circ g$.
    '''
    def __init__(self, f, g):
        SymbolicOperation.__init__(self, (f,g))

    def _repr_(self):
        ops = self._operands
        return "%s(%s)"% (ops[0]._repr_(), ops[1]._repr_())

    def _latex_(self):
        ops = self._operands
        return r"%s \left( %s \right)"% (latex(ops[0]), latex(ops[1]))

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            ops = self._operands
            m = ops[0]._maxima_(maxima)(ops[1]._maxima_(maxima))
            self.__maxima = m
            return m


symbols = {operator.add:' + ', operator.sub:' - ', operator.mul:'*',
            operator.div:'/', operator.pow:'^'}


class SymbolicArithmetic(SymbolicOperation):
    r'''
    Represents the result of an arithemtic operation on
    $f$ and $g$.
    '''

    def __init__(self, f, g, op):
        if not isinstance(f, SymbolicExpression) or not isinstance(g,
                SymbolicExpression):
            raise TypeError, "Symbolic Arithmetic is only defined on"+\
            " SymbolicExpression objects."

        self._operator = op
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

        SymbolicOperation.__init__(self, (f, g))

    def _repr_(self):
        ops = self._operands
        s0 = str(ops[0])
        s1 = str(ops[1])
        return self.parenthesize(s0, s1, self._operator)

    def _is_atomic(self):
        return self._atomic

    def parenthesize(self, s0, s1, op, op_symbols=symbols, paren_symbols=('(',')')):
        '''
        Properly add parentheses with correct operator precedence in a string
        of ASCII math or latex.
        '''
        if op in [operator.mul, operator.div, operator.pow]:
            if '+' in s0 or '-' in s0:
                # avoid drawing parens if s0 is atomic, since they are redundant
                if not self._operands[0]._is_atomic():
                    s0 = '%s%s%s' % (paren_symbols[0],s0,paren_symbols[1])
        if op in [operator.mul, operator.div, operator.sub]:
            if '+' in s1 or '-' in s1:
                # avoid drawing parens if s1 is atomic
                if not self._operands[1]._is_atomic():
                    s1 = '%s%s%s' % (paren_symbols[0],s1,paren_symbols[1])
        return "%s%s%s" % (s0, op_symbols[op], s1)

    def _latex_(self):
        parens = ('\\left(', '\\right)')
        ops = [x._latex_() for x in self._operands]
        op_symbols = {operator.mul:' \\cdot ', operator.div:''}

        if self._operator is operator.add:
            return self.parenthesize(ops[0], ops[1], self._operator,
                    paren_symbols=parens)
        elif self._operator is operator.sub:
            return self.parenthesize(ops[0], ops[1], self._operator,
                    paren_symbols=parens)
        elif self._operator == operator.mul:
            return self.parenthesize('{%s}'%ops[0],'{%s}'%ops[1], self._operator,
                    op_symbols, paren_symbols=parens)
        elif self._operator == operator.div:
            return self.parenthesize('\\frac{%s}'%ops[0],'{%s}'%ops[1],
                    self._operator, op_symbols=op_symbols, paren_symbols=parens)
        elif self._operator == operator.pow:
            return self.parenthesize(ops[0],'{%s}'%ops[1], self._operator,
                    paren_symbols=parens)
        else: raise NotImplementedError, 'Operator %s unkown' % self._operator


    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            ops = self._operands
            m = self._operator(ops[0]._maxima_(maxima),
                                  ops[1]._maxima_(maxima))
            self.__maxima = m
            return m


import re

class SymbolicVariable(SymbolicExpression):
    def __init__(self, name):
        SymbolicExpression.__init__(self)
        self._name = name
        if len(name) == 0:
            raise ValueError, "variable name must be nonempty"

    def __call__(self, x):
        raise AttributeError, "A symbolic variable is not callable."

    def _repr_(self):
        return self._name

    def _latex_(self):
        try:
            return self.__latex
        except AttributeError:
            pass
        a = self._name
        if len(a) > 1:
            m = re.search('\d+$',a)
            if m is None:
                a = tex_varify(a)
            else:
                b = a[:m.start()]
                a = '%s_{%s}'%(tex_varify(b), a[m.start():])

        self.__latex = a
        return a

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            m = maxima(self._name)
            self.__maxima = m
            return m

common_varnames = ['alpha',
                   'beta',
                   'gamma',
                   'delta',
                   'epsilon',
                   'zeta',
                   'eta',
                   'theta',
                   'iota',
                   'kappa',
                   'lambda',
                   'mu',
                   'nu',
                   'xi',
                   'pi',
                   'rho',
                   'sigma',
                   'tau',
                   'upsilon',
                   'phi',
                   'chi',
                   'psi',
                   'omega']
common_varnames.append([i.capitalize() for i in common_varnames])


def tex_varify(a):
    # TODO: add more
    if a in ['theta', 'eta', 'alpha']:
        return "\\" + a
    else:
        return '\\mbox{%s}'%a

_vars = {}
def var(s):
    try:
        return _vars[s]
    except KeyError:
        v = SymbolicVariable(s)
        _vars[s] = v
        return v

t = var('t')
x = var('x')
y = var('y')
z = var('z')
w = var('w')

class Function_sin(PrimitiveFunction):
    '''
    The sine function
    '''
    def _repr_(self):
        return "sin"

    def _latex_(self):
        return "\\sin"

sin = Function_sin()

class Function_cos(PrimitiveFunction):
    '''
    The cosine function
    '''
    def _repr_(self):
        return "cos"

    def _latex_(self):
        return "\\cos"

cos = Function_cos()

class Function_sec(PrimitiveFunction):
    '''
    The secant function
    '''
    def _repr_(self):
        return "sec"

    def _latex_(self):
        return "\\sec"

sec = Function_sec()


#######################################################

def symbolic_expression_from_maxima_string(x):
    maxima.eval('listdummyvars: false')
    maxima.eval('_tmp_: %s'%x)
    r = maxima.eval('listofvars(_tmp_)')[1:-1]
    if len(r) == 0:
        d = {}
    else:
        # Now r is a list of all the indeterminate variables that
        # appear in the expression x.
        v = r.split(',')
        d = dict([(a, var(a)) for a in v])
    s = maxima.eval('_tmp_')
    return SymbolicExpressionRing(sage_eval(s, d))

def symbolic_expression_from_maxima_element(x):
    return symbolic_expression_from_maxima_string(x.name())

