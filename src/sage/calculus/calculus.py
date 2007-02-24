from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_RealNumber, is_ComplexNumber, RR,
                            Integer, Rational, CC)
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
import sage.functions.constants as c

from sage.calculus.equations import SymbolicEquation
from sage.rings.real_mpfr import RealNumber
from sage.rings.complex_number import ComplexNumber

import pdb

# There will only ever be one instance of this class
class SymbolicExpressionRing_class(CommutativeRing):
    """
    A Ring of formal expressions.
    """
    def __init__(self):
        self._default_precision = 53 # default precision bits
        ParentWithBase.__init__(self, RR)

    def __call__(self, x):
        return self._coerce_impl(x)

    def _coerce_impl(self, x):
        if isinstance(x, CallableFunction):
            return x._expr
        elif isinstance(x, SymbolicExpression):
            return x
        elif isinstance(x, MaximaElement):
            return symbolic_expression_from_maxima_element(x)
        elif is_Polynomial(x):
            return SymbolicPolynomial(x)
        elif isinstance(x, (RealNumber, float, Constant, Integer, int, Rational, ComplexNumber)):
            return Constant_object(x)
        else:
            # we really want to raise a type error, here
            #raise TypeError, 'cannot coerce %s into a SymbolicExpression.'%x
            pass

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

# ... and here it is:
SymbolicExpressionRing = SymbolicExpressionRing_class()
SER = SymbolicExpressionRing
# conversions is the dict of the form system:command

class SymbolicExpression(RingElement):
    """
    A Symbolic Expression.
    """
    def __init__(self, conversions={}):
        RingElement.__init__(self, SymbolicExpressionRing)
        self._simp = None

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            m = maxima(self._maxima_init_())
            self.__maxima = m
            return m

    def hash(self):
        return hash(maxima(self))

    def __float__(self):
        return float(self._obj)

    def _plot_(self, x=None, *args, **kwds):
        from sage.plot.plot import plot
        # see if the user passed a param
        try:
            param = kwds['param']
        except KeyError:
            return SageObject.plot(self, *args, **kwds)
        del kwds['param']
        return plot(self.function(param), *args, **kwds)


    def __eq__(self, right):
        return SymbolicEquation(self, right)

    def __cmp__(self, right):
        return cmp(maxima(self), maxima(right))

    def _neg_(self):
        return SymbolicArithmetic([self], operator.neg)

    def _add_(self, right):
        # if we are adding a negation, instead subtract the operand of negation
        #if isinstance(right, SymbolicArithmetic):
        #    if right._operator is operator.neg:
        #        return SymbolicArithmetic([self, right._operands[0]], operator.sub)
        #elif isinstance(right, Symbolic_object) and right < 0:
        #    return SymbolicArithmetic([self, SER(abs(right._obj))], operator.sub)
        #else:
        return SymbolicArithmetic([self, right], operator.add)

    def _sub_(self, right):
        return SymbolicArithmetic([self, right], operator.sub)

    def _mul_(self, right):
        # do some simplification... pull out negatives from the operands and put it
        # in front of this multiplication
        #if isinstance(self, SymbolicArithmetic) and isinstance(right, SymbolicArithmetic):
        #    if self._operator is operator.neg and right._operator is operator.neg:
        #        s_unneg = self._operands[0]
        #        r_unneg = right._operands[0]
        #        return SymbolicArithmetic([s_unneg, r_unneg], operator.mul)
        #if isinstance(self, SymbolicArithmetic):
        #    if self._operator is operator.neg and (not isinstance(right, SymbolicArithmetic) \
        #    or not (right._operator is operator.neg)):
        #        s_unneg = self._operands[0]
        #        return -SymbolicArithmetic([s_unneg, right], operator.mul)
        #if isinstance(right, SymbolicArithmetic):
        #    if not isinstance(self, SymbolicArithmetic) or (not (self._operator is operator.neg)) \
        #    and right._operator is operator.neg:
        #        r_unneg = right._operands[0]
        #        return -SymbolicArithmetic([self, r_unneg], operator.mul)

        return SymbolicArithmetic([self, right], operator.mul)

    def _div_(self, right):
        # do some simplification... pull out negatives from the operands and put it
        # in front of this division
        #if isinstance(self, SymbolicArithmetic) and isinstance(right, SymbolicArithmetic):
        #    if self._operator is operator.neg and right._operator is operator.neg:
        #        s_unneg = self._operands[0]
        #        r_unneg = right._operands[0]
        #        return SymbolicArithmetic([s_unneg, r_unneg], operator.div)
        #elif isinstance(self, SymbolicArithmetic):
        #    if self._operator is operator.neg and (not isinstance(right, SymbolicArithmetic) \
        #    or not (right._operator is operator.neg)):
        #        s_unneg = self._operands[0]
        #        return -SymbolicArithmetic([s_unneg, right], operator.div)
        #elif isinstance(right, SymbolicArithmetic):
        #    if not isinstance(self, SymbolicArithmetic) or (not (self._operator is operator.neg)) \
        #    and right._operator is operator.neg:
        #        r_unneg = right._operands[0]
        #        return -SymbolicArithmetic([self, r_unneg], operator.div)

        return SymbolicArithmetic([self, right], operator.div)

    def __pow__(self, right):
        right = self.parent()(right)
        return SymbolicArithmetic([self, right], operator.pow)

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
           sin(x) + x*cos(y)
           sage: g(t,z)
           sin(t) + t*cos(z)
           sage: g(x^2, x^y)
           sin(x^2) + x^2*cos(x^y)
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
        for i in range(len(args)):
            #pdb.set_trace()
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
    def integral(self, v):
        """
        EXAMPLES:
            sage: h = sin(x)/(cos(x))^2
            sage: h.integral(x)
            1/cos(x)
        """
        if not isinstance(v, SymbolicVariable):
            raise TypeError, "second argument to diff must be a variable"
        return self.parent()(self._maxima_().integrate(v))


    ###################################################################
    # simplify
    ###################################################################
    def simplify(self):
        if self._simp is None:
            self._simp = self.parent()(self._maxima_())
        return self._simp

    def simplify_trig(self):
        r"""
        Employs the identities $\sin(x)^2 + \cos(x)^2 = 1$ and
        $\cosh(x)^2 - \sin(x)^2 = 1$ to simplify expressions
        containing tan, sec, etc., to sin, cos, sinh, cosh.

        trig_simplify() and simplify_trig() are the same method.

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
        return self.parent()(self._maxima_().fullratsimp())

    rational_simplify = simplify_rational

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



    ###################################################################
    # substitute
    ###################################################################
    def substitute(self, dict=None, **kwds):
        """
        Takes the symbolic variables given as dict keys or as keywords and
        replaces them with the symbolic expressions given as dict values or as
        keyword values. Also run when you call a SymbolicExpression.

        EXAMPLES:
            sage: u = (x^3 - 3*y + 4*t)
            sage: u.substitute(x=y, y=t)
            y^3 - 3*t + 4*t

            sage: f = sin(x)^2 + 32*x^(y/2)
            sage: f.substitute(x=2, y = 10)
            sin(2)^2 + 1024

            sage: f(x=pi, y=t)
            sin(pi)^2 + 32*pi^(t/2)

        AUTHORS:
            -- Bobby Moretti: Initial version
        """
        if dict is not None:
            for k, v in dict.iteritems():
               kwds[str(k)] = v
        # find the keys from the keywords
        return self._recursive_sub(kwds)

    def _recursive_sub(self, kwds):
        # if we have operands, call ourself on each operand
        try:
            ops = self._operands
        except AttributeError:
            pass
        # if we are a symbolic variable, we're at a leaf node
        if isinstance(self, SymbolicVariable):
            s = str(self)
            # do the replacement if needed
            if s in kwds:
                return kwds[s]
            else:
                return self
        elif isinstance(self, SymbolicConstant):
            return self
        elif isinstance(self, SymbolicArithmetic):
            new_ops = [op._recursive_sub(kwds) for op in ops]
            arith = SymbolicArithmetic(new_ops, self._operator)
            return  arith._operator(*(arith._operands))
        elif isinstance(self, SymbolicComposition):
            return SymbolicComposition(ops[0], ops[1]._recursive_sub(kwds))

class PrimitiveFunction(SymbolicExpression):
    def __init__(self):
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
            raise TypeError, 'cannot coerce %s into a SymbolicExpression.'%x

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

class CallableFunction(RingElement):
    r'''
    A callable, symbolic function that knows the variables on which it depends.
    '''
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
            if len(self._varlist) != 1:
                raise TypeError, "Cannot coerce %s to a float" % self

            return float(self.obj)

    # TODO: should len(args) == len(vars)?
    def __call__(self, *args):
        vars = self._varlist

        dct = {}
        for i in range(len(args)):
            dct[vars[i]] = args[i]

        return self._expr.substitute(dct)

    def _repr_(self):
        if len(self._varlist) == 1:
            return "%s |--> %s" % (self._varlist[0], self._expr._repr_())
        else:
            vars = ", ".join(map(str, self._varlist))
            return "(%s) |--> %s" % (vars, self._expr._repr_())

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
        return CallableFunction(result, self._unify_varlit(right))

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
            (x, y, z) |--> sin(x + y + z)
            sage: g(x, y) = 2*x + y
            sage: g
            (x, y) |--> 2*x + y
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
            (x, y, t) |--> y*(x^2 - t)
            sage: g(x, y, w) = x + y - cos(w)
            sage: f._unify_varlists(g)
            (x, y, t, w)
            sage: g._unify_varlists(f)
            (x, y, t, w)

            sage: f(x,y, t) = x+y
            sage: g(x, y, w) = w + t
            sage: f._unify_varlists(g)
            (x, y, t, w)
            sage: g._unify_varlists(f)
            (x, y, t, w)
            sage: f + g
            (x, y, t, w) |--> x + y + w + t

        AUTHORS:
            -- Bobby Moretti, thanks to William Stein for the rules

        """
        a = self._varlist
        b = x._varlist

        # Rule #1
        if a == b:
            return self._varlist

        # Rule #dd
        #elif set(a).issubset(set(b)):
        #    return b
        #elif set(b).issubset(set(a)):
        #    return a

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

    def __pow__(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj ** right._obj)
            except TypeError:
                pass
        return SymbolicExpression.__pow__(self, right)

    def _add_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj + right._obj)
            except TypeError:
                pass
        return SymbolicExpression._add_(self, right)

    def _sub_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj - right._obj)
            except TypeError:
                pass
        return SymbolicExpression._sub_(self, right)

    def _mul_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj * right._obj)
            except TypeError:
                pass
        return SymbolicExpression._mul_(self, right)

    def _div_(self, right):
        if isinstance(right, Symbolic_object):
            try:
                return Symbolic_object(self._obj / right._obj)
            except TypeError:
                pass
        return SymbolicExpression._div_(self, right)

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

    # for now we do nothing except pass the info on to the superconstructor. It's
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
        self._operands = [SER(op) for op in operands]

class SymbolicComposition(SymbolicOperation):
    r'''
    Represents the symbolic composition of $f \circ g$.
    '''
    def __init__(self, f, g):
        SymbolicOperation.__init__(self, [f,g])

    def _is_atomic(self):
        return True

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

    def __init__(self, operands, op):
        SymbolicOperation.__init__(self, operands)

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

        self._operator = op

    def _is_atomic(self):
        return self._atomic

    def _repr_(self):

        ops = self._operands
        op = self._operator

        s = [str(o) for o in ops]
        # for the left operand, we need to surround it in parens when the
        # operator is mul/div/pow, and when the left operand contains an
        # operation of lower precedence
        if op in [operator.mul, operator.div, operator.pow]:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub):
                if not ops[0]._is_atomic():
                    s[0] = '(%s)' % s[0]
        # for the right operand, we need to surround it in parens when the
        # operation is mul/div/sub and the, and when the right operand contains
        # a + or -.
        if op in [operator.mul, operator.div, operator.sub]:
            if ops[1]._has_op(operator.add) or ops[1]._has_op(operator.sub):
                # avoid drawing parens if s1 an atomic operation
                if not ops[1]._is_atomic():
                    s[1] = '(%s)' % s[1]
        # if we have a compound expression on the right, then we need parens on
        # the right... but in TeX, this is expressed by font size and position
        elif op is operator.pow:
            if ops[1]._has_op(operator.add) or  \
            ops[1]._has_op(operator.sub) or  \
            ops[1]._has_op(operator.mul) or  \
            ops[1]._has_op(operator.div) or  \
            ops[1]._has_op(operator.pow):
                s[1] = '(%s)' % s[1]
        if op is operator.neg:
            return '-%s' % s[0]
        else:
            return '%s%s%s' % (s[0], symbols[op], s[1])

    def _latex_(self):
        op = self._operator
        ops = self._operands
        ops_tex = [x._latex_() for x in self._operands]

        s = [o._latex_() for o in ops]

        # follow a very similar logic to _repr_, but we don't parenthesize
        # exponents
        if op in [operator.mul, operator.div, operator.pow]:
            if ops[0]._has_op(operator.add) or ops[0]._has_op(operator.sub):
                if not ops[0]._is_atomic():
                    s[0] = '\\left( %s \\right)' % s[0]
        if op in [operator.mul, operator.div, operator.sub]:
            if ops[1]._has_op(operator.add) or ops[1]._has_op(operator.sub):
                if not ops[1]._is_atomic():
                    s[1] = '\\left( %s \\right)' % s[1]

        if op is operator.add:
            return '%s + %s' % (s[0], s[1])
        elif op is operator.sub:
            return '%s - %s' % (s[0], s[1])
        elif op is operator.mul:
            return '{%s \\cdot %s}' % (s[0], s[1])
        elif op is operator.div:
            return '\\frac{%s}{%s}' % (s[0], s[1])
        elif op is operator.pow:
            return '{%s}^{%s} ' % (s[0], s[1])
        elif op is operator.neg:
            return '-%s' % s[0]

    def _maxima_(self, maxima=maxima):
        try:
            return self.__maxima
        except AttributeError:
            ops = self._operands
            if self._operator is operator.neg:
                m = self._operator(ops[0]._maxima_(maxima))
            else:
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

    def __cmp__(self, right):
        return cmp(self._repr_(), right._repr_())

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
            m = re.search('(\d|[.,])+$',a)
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
    r''' Create a symbolic variable with the name \emph{s}'''
    try:
        return _vars[s]
    except KeyError:
        v = SymbolicVariable(s)
        _vars[s] = v
        return v

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
I = var('I')
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

_syms = {}

class Function_erf(PrimitiveFunction):
    r'''
    The error function, defined as $\text{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x
    e^{-t^2} dt$.
    '''

    def _repr_(self):
        return "erf"

    def _latex_(self):
        return "\\text{erf}"

    def _is_atomic(self):
        return True

erf = Function_erf()
_syms['erf'] = erf

class Function_sin(PrimitiveFunction):
    '''
    The sine function
    '''
    def _repr_(self):
        return "sin"

    def _latex_(self):
        return "\\sin"

    def _is_atomic(self):
        return True

sin = Function_sin()
_syms['sin'] = sin

class Function_cos(PrimitiveFunction):
    '''
    The cosine function
    '''
    def _repr_(self):
        return "cos"

    def _latex_(self):
        return "\\cos"

    def _is_atomic(self):
        return True

cos = Function_cos()
_syms['cos'] = cos

class Function_sec(PrimitiveFunction):
    '''
    The secant function
    '''
    def _repr_(self):
        return "sec"

    def _latex_(self):
        return "\\sec"

    def _is_atomic(self):
        return True

sec = Function_sec()
_syms['sec'] = sec

class Function_log(PrimitiveFunction):
    '''
    The log function
    '''
    def _repr_(self):
        return "log"

    def _latex_(self):
        return "\\log"

    def _is_atomic(self):
        return True

log = Function_log()
_syms['log'] = log

class Function_sqrt(PrimitiveFunction):
    '''
    The (positive) square root function.
    '''
    def _repr_(self):
        return "sqrt"

    def _latex_(self):
        return "\\sqrt"

    def _is_atomic(self):
        return True

sqrt = Function_sqrt()
_syms['sqrt'] = sqrt

#######################################################
symtable = {'%pi':'_Pi_', '%e': '_E_', '%i':'_I_'}
import sage.functions.constants as c
_syms['_Pi_'] = SER(c.Pi)
_syms['_E_'] = SER(c.E)
_syms['_I_'] = SER(CC.gen(0))  # change when we create a symbolic I.

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
    for x, y in symtable.iteritems():
        s = s.replace(x, y)
    return SymbolicExpressionRing(sage_eval(s, _syms))

def symbolic_expression_from_maxima_element(x):
    return symbolic_expression_from_maxima_string(x.name())
