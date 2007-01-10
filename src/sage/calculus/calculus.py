from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_RealNumber, is_ComplexNumber, RR)
#import sage.rings.rational
from sage.rings.integer import Integer
from sage.structure.element import RingElement
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex
from sage.structure.sage_object import SageObject

# DEBUGGER!!
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
        elif isinstance(x, (sage.rings.integer.Integer,
                            sage.rings.rationalRational,
                            int, long, float, complex)):
            return Constant_object(x)

        else:
            return Symbolic_object(x)

# There will only ever be one instance of this class
class SymbolicExpressionRing_class(CommutativeRing):
    '''
    A Ring of formal expressions.
    '''
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def _coerce_impl(self, x):
        if is_Polynomial(x):
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
        return 0

# ... and here it is:
SymbolicExpressionRing = SymbolicExpressionRing_class()


class SymbolicExpression(RingElement):
    r"""
    A Symbolic Expression in acoordance with SEP #1.

    """
    # conversions is the dict of the form system:command
    def __init__(self, conversions={}):
        #db.set_trace()
        RingElement.__init__(self, SymbolicExpressionRing)

    def _add_(self, right):
        # let's try to nip some simplification in the bud
        if right is self:
            return SymbolicArithmetic(Constant_object(2), self, operator.mul)
        return SymbolicArithmetic(self, right, operator.add)


    def _sub_(self, right):
        return SymbolicArithmetic(self, right, operator.sub)

    def _mul_(self, right):
        return SymbolicArithmetic(self, right, operator.mul)

    def _div_(self, right):
        return SymbolicArithmetic(self, right, operator.div)

    def __pow__(self, right):
       return SymbolicArithmetic(self, right, operator.pow)

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

    def str(self, bits=None):
        if bits is None:
            return str(self._obj)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self._obj))

    def _maxima_(self, maxima):
        return maxima(self._obj)

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

    def _maxima_(self, mxma_inst):
        ops = self._operands
        return ops[0]._maxima_(mxma_inst)(ops[1]._maxima_(mxma_inst))


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
        SymbolicOperation.__init__(self, (f, g))

    def _repr_(self):
        ops = self._operands
        return "(%s%s%s)"% (ops[0], symbols[self._operator],ops[1])

    def _latex_(self):
        ops = self._operands
        if self._operator == operator.add:
            return '%s + %s:' % (ops[0], ops[1])
        elif self._operator == operator.sub:
            return '%s - %s:' % (ops[0], ops[1])
        elif self._operator == operator.mul:
            return '%s \\cdot %s:' % (ops[0], ops[1])
        elif self._operator == operator.div:
            return '\\frac{%s}{%s}' % (ops[0], ops[1])
        elif self._operator == operator.pow:
            return '%s^{%s}' % (ops[0], ops[1])
        else: raise NotImplementedError, 'Operator %s unkown' % self._operator

    def _maxima_(self, mxma_inst):
        ops = self._operands
        return self._operator(ops[0]._maxima_(mxma_inst),
            ops[1]._maxima_(mxma_inst))


class SymbolicVariable(SymbolicExpression):
    def __init__(self, name):
        SymbolicExpression.__init__(self)
        self._name = name

    def __call__(self, x):
        raise AttributeError, "A symbolic variable is not callable."

    def _repr_(self):
        return self._name

    def _latex_(self):
        return self._name

t = SymbolicVariable('t')
x = SymbolicVariable('x')

class Function_sin(PrimitiveFunction):
    '''
    The sine function
    '''
    def __init__(self):
        PrimitiveFunction.__init__(self)

    def _repr_(self):
        return "sin"

    def _latex_(self):
        return "\\sin"

    def __call__(self, x):
        return SymbolicComposition(self, x)

class Function_cos(PrimitiveFunction):
    '''
    The cosine function
    '''
    def __init__(self):
        PrimitiveFunction.__init__(self)

    def _repr_(self):
        return "cos"

    def _latex_(self):
        return "\\cos"

    def __call__(self, x):
        return SymbolicComposition(self, x)

sin = Function_sin()
cos = Function_cos()
