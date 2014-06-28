"""
Conversion of symbolic expressions to other types

This module provides routines for converting new symbolic expressions
to other types.  Primarily, it provides a class :class:`Converter`
which will walk the expression tree and make calls to methods
overridden by subclasses.
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

import operator as _operator
from sage.symbolic.ring import SR,var
from sage.symbolic.pynac import I
from sage.functions.all import exp
from sage.symbolic.operators import arithmetic_operators, relation_operators, FDerivativeOperator
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic
from functools import reduce
GaussianField = I.pyobject().parent()

class FakeExpression(object):
    r"""
    Pynac represents `x/y` as `xy^{-1}`.  Often, tree-walkers would prefer
    to see divisions instead of multiplications and negative exponents.
    To allow for this (since Pynac internally doesn't have division at all),
    there is a possibility to pass use_fake_div=True; this will rewrite
    an Expression into a mixture of Expression and FakeExpression nodes,
    where the FakeExpression nodes are used to represent divisions.
    These nodes are intended to act sufficiently like Expression nodes
    that tree-walkers won't care about the difference.
    """

    def __init__(self, operands, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: FakeExpression([x, y], operator.div)
            FakeExpression([x, y], <built-in function div>)
        """
        self._operands = operands
        self._operator = operator

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: FakeExpression([x, y], operator.div)
            FakeExpression([x, y], <built-in function div>)
        """
        return "FakeExpression(%r, %r)"%(self._operands, self._operator)

    def pyobject(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: f.pyobject()
            Traceback (most recent call last):
            ...
            TypeError: self must be a numeric expression
        """
        raise TypeError('self must be a numeric expression')

    def operands(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: f.operands()
            [x, y]
        """
        return self._operands

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: f[0]
            x
        """
        return self._operands[i]

    def operator(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: f.operator()
            <built-in function div>
        """
        return self._operator

    def _fast_callable_(self, etb):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: fast_callable(f, vars=['x','y']).op_list()
            [('load_arg', 0), ('load_arg', 1), 'div', 'return']
        """
        return fast_callable(self, etb)

    def _fast_float_(self, *vars):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.div)
            sage: fast_float(f, 'x', 'y').op_list()
            [('load_arg', 0), ('load_arg', 1), 'div', 'return']
        """
        return fast_float(self, *vars)

class Converter(object):
    def __init__(self, use_fake_div=False):
        """
        If use_fake_div is set to True, then the converter will try to
        replace expressions whose operator is operator.mul with the
        corresponding expression whose operator is operator.div.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c.use_fake_div
            True
        """
        self.use_fake_div = use_fake_div

    def __call__(self, ex=None):
        """
        .. note::

           If this object does not have an attribute *ex*, then an argument
           must be passed into :meth`__call__`::

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c(SR(2))
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
            sage: c(x+2)
            Traceback (most recent call last):
            ...
            NotImplementedError: arithmetic
            sage: c(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: symbol
            sage: c(x==2)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
            sage: c(sin(x))
            Traceback (most recent call last):
            ...
            NotImplementedError: composition
            sage: c(function('f', x).diff(x))
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative

        We can set a default value for the argument by setting
        the ``ex`` attribute::

            sage: c.ex = SR(2)
            sage: c()
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
        """
        if ex is None:
            ex = self.ex

        try:
            obj = ex.pyobject()
            return self.pyobject(ex, obj)
        except TypeError as err:
            if 'self must be a numeric expression' not in err:
                raise err

        operator = ex.operator()
        if operator is None:
            return self.symbol(ex)

        if operator in arithmetic_operators:
            if getattr(self, 'use_fake_div', False) and operator is _operator.mul:
                div = self.get_fake_div(ex)
                return self.arithmetic(div, div.operator())
            return self.arithmetic(ex, operator)
        elif operator in relation_operators:
            return self.relation(ex, operator)
        elif isinstance(operator, FDerivativeOperator):
            return self.derivative(ex, operator)
        else:
            return self.composition(ex, operator)

    def get_fake_div(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c.get_fake_div(sin(x)/x)
            FakeExpression([sin(x), x], <built-in function div>)
            sage: c.get_fake_div(-1*sin(x))
            FakeExpression([sin(x)], <built-in function neg>)
            sage: c.get_fake_div(-x)
            FakeExpression([x], <built-in function neg>)
            sage: c.get_fake_div((2*x^3+2*x-1)/((x-2)*(x+1)))
            FakeExpression([2*x^3 + 2*x - 1, FakeExpression([x + 1, x - 2], <built-in function mul>)], <built-in function div>)

        Check if #8056 is fixed, i.e., if numerator is 1.::

            sage: c.get_fake_div(1/pi/x)
            FakeExpression([1, FakeExpression([pi, x], <built-in function mul>)], <built-in function div>)
        """
        d = []
        n = []
        for arg in ex.operands():
            ops = arg.operands()
            try:
                if arg.operator() is _operator.pow and repr(ops[1]) == '-1':
                    d.append(ops[0])
                else:
                    n.append(arg)
            except TypeError:
                n.append(arg)

        len_d = len(d)
        if len_d == 0:
            repr_n = map(repr, n)
            if len(n) == 2 and "-1" in repr_n:
                a = n[0] if repr_n[1] == "-1" else n[1]
                return FakeExpression([a], _operator.neg)
            else:
                return ex
        elif len_d == 1:
            d = d[0]
        else:
            d = FakeExpression(d, _operator.mul)

        if len(n) == 0:
            return FakeExpression([SR.one_element(), d], _operator.div)
        elif len(n) == 1:
            n = n[0]
        else:
            n = FakeExpression(n, _operator.mul)

        return FakeExpression([n,d], _operator.div)

    def pyobject(self, ex, obj):
        """
        The input to this method is the result of calling
        :meth:`pyobject` on a symbolic expression.

        .. note::

           Note that if a constant such as ``pi`` is encountered in
           the expression tree, its corresponding pyobject which is an
           instance of :class:`sage.symbolic.constants.Pi` will be
           passed into this method.  One cannot do arithmetic using
           such an object.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = SR(1)
            sage: Converter().pyobject(f, f.pyobject())
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
        """
        raise NotImplementedError("pyobject")

    def symbol(self, ex):
        """
        The input to this method is a symbolic expression which
        corresponds to a single variable.  For example, this method
        could be used to return a generator for a polynomial ring.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: Converter().symbol(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: symbol
        """
        raise NotImplementedError("symbol")

    def relation(self, ex, operator):
        """
        The input to this method is a symbolic expression which
        corresponds to a relation.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: import operator
            sage: Converter().relation(x==3, operator.eq)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
            sage: Converter().relation(x==3, operator.lt)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
        """
        raise NotImplementedError("relation")

    def derivative(self, ex, operator):
        """
        The input to this method is a symbolic expression which
        corresponds to a relation.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: a = function('f', x).diff(x); a
            D[0](f)(x)
            sage: Converter().derivative(a, a.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative
        """
        raise NotImplementedError("derivative")

    def arithmetic(self, ex, operator):
        """
        The input to this method is a symbolic expression and the
        infix operator corresponding to that expression. Typically,
        one will convert all of the arguments and then perform the
        operation afterward.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = x + 2
            sage: Converter().arithmetic(f, f.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: arithmetic
        """
        raise NotImplementedError("arithmetic")

    def composition(self, ex, operator):
        """
        The input to this method is a symbolic expression and its
        operator.  This method will get called when you have a symbolic
        function application.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = sin(2)
            sage: Converter().composition(f, f.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: composition
        """
        raise NotImplementedError("composition")

class InterfaceInit(Converter):
    def __init__(self, interface):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: a = pi + 2
            sage: m(a)
            '(%pi)+(2)'
            sage: m(sin(a))
            'sin((%pi)+(2))'
            sage: m(exp(x^2) + pi + 2)
            '(%pi)+(exp((x)^(2)))+(2)'

        """
        self.name_init = "_%s_init_"%interface.name()
        self.interface = interface
        self.relation_symbols = interface._relation_symbols()

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.symbol(x)
            'x'
            sage: f(x) = x
            sage: m.symbol(f)
            'x'
        """
        return repr(SR(ex))

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: ii = InterfaceInit(gp)
            sage: f = 2+I
            sage: ii.pyobject(f, f.pyobject())
            'I + 2'

            sage: ii.pyobject(SR(2), 2)
            '2'

            sage: ii.pyobject(pi, pi.pyobject())
            'Pi'
        """
        if (self.interface.name() in ['pari','gp'] and
            isinstance(obj, NumberFieldElement_quadratic) and
            obj.parent() == GaussianField):
            return repr(obj)
        try:
            return getattr(obj, self.name_init)()
        except AttributeError:
            return repr(obj)

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.relation(x==3, operator.eq)
            'x = 3'
            sage: m.relation(x==3, operator.lt)
            'x < 3'
        """
        return "%s %s %s"%(self(ex.lhs()), self.relation_symbols[operator],
                           self(ex.rhs()))

    def derivative(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: f = function('f')
            sage: a = f(x).diff(x); a
            D[0](f)(x)
            sage: print m.derivative(a, a.operator())
            diff('f(x), x, 1)
            sage: b = f(x).diff(x, x)
            sage: print m.derivative(b, b.operator())
            diff('f(x), x, 2)

        We can also convert expressions where the argument is not just a
        variable, but the result is an "at" expression using temporary
        variables::

            sage: y = var('y')
            sage: t = (f(x*y).diff(x))/y
            sage: t
            D[0](f)(x*y)
            sage: m.derivative(t, t.operator())
            "at(diff('f(t0), t0, 1), [t0 = x*y])"
        """
        #This code should probably be moved into the interface
        #object in a nice way.
        from sage.symbolic.ring import is_SymbolicVariable
        if self.name_init != "_maxima_init_":
            raise NotImplementedError
        args = ex.operands()
        if (not all(is_SymbolicVariable(v) for v in args) or
            len(args) != len(set(args))):
            # An evaluated derivative of the form f'(1) is not a
            # symbolic variable, yet we would like to treat it like
            # one. So, we replace the argument `1` with a temporary
            # variable e.g. `t0` and then evaluate the derivative
            # f'(t0) symbolically at t0=1. See trac #12796.
            temp_args=[var("t%s"%i) for i in range(len(args))]
            f = operator.function()
            params = operator.parameter_set()
            params = ["%s, %s"%(temp_args[i], params.count(i)) for i in set(params)]
            subs = ["%s = %s"%(t,a) for t,a in zip(temp_args,args)]
            return "at(diff('%s(%s), %s), [%s])"%(f.name(),
                ", ".join(map(repr,temp_args)),
                ", ".join(params),
                ", ".join(subs))

        f = operator.function()
        params = operator.parameter_set()
        params = ["%s, %s"%(args[i], params.count(i)) for i in set(params)]

        return "diff('%s(%s), %s)"%(f.name(),
                                    ", ".join(map(repr, args)),
                                    ", ".join(params))

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.arithmetic(x+2, operator.add)
            '(x)+(2)'
        """
        args = ["(%s)"%self(op) for op in ex.operands()]
        return arithmetic_operators[operator].join(args)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.composition(sin(x), sin)
            'sin(x)'
            sage: m.composition(ceil(x), ceil)
            'ceiling(x)'

            sage: m = InterfaceInit(mathematica)
            sage: m.composition(sin(x), sin)
            'Sin[x]'
        """
        ops = ex.operands()
        #FIXME: consider stripping pyobjects() in ops
        if hasattr(operator, self.name_init + "evaled_"):
            return getattr(operator, self.name_init + "evaled_")(*ops)
        else:
            ops = map(self, ops)
        try:
            op = getattr(operator, self.name_init)()
        except (TypeError, AttributeError):
            op = repr(operator)

        return self.interface._function_call_string(op,ops,[])

#########
# Sympy #
#########
class SympyConverter(Converter):
    """
    Converts any expression to SymPy.

    EXAMPLE::

        sage: import sympy
        sage: var('x,y')
        (x, y)
        sage: f = exp(x^2) - arcsin(pi+x)/y
        sage: f._sympy_()
        exp(x**2) - asin(x + pi)/y
        sage: _._sage_()
        -arcsin(pi + x)/y + e^(x^2)

        sage: sympy.sympify(x) # indirect doctest
        x

    TESTS:

    Make sure we can convert I (trac #6424)::

        sage: bool(I._sympy_() == I)
        True
        sage: (x+I)._sympy_()
        x + I

    """
    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter
            sage: s = SympyConverter()
            sage: f = SR(2)
            sage: s.pyobject(f, f.pyobject())
            2
            sage: type(_)
            <class 'sympy.core.numbers.Integer'>
        """
        try:
            return obj._sympy_()
        except AttributeError:
            return obj

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter
            sage: s = SympyConverter()
            sage: f = x + 2
            sage: s.arithmetic(f, f.operator())
            x + 2
        """
        import sympy
        operator = arithmetic_operators[operator]
        ops = [sympy.sympify(self(a), evaluate=False) for a in ex.operands()]
        if operator == "+":
            return sympy.Add(*ops)
        elif operator == "*":
            return sympy.Mul(*ops)
        elif operator == "-":
            return sympy.Sub(*ops)
        elif operator == "/":
            return sympy.Div(*ops)
        elif operator == "^":
            return sympy.Pow(*ops)
        else:
            raise NotImplementedError

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter
            sage: s = SympyConverter()
            sage: s.symbol(x)
            x
            sage: type(_)
            <class 'sympy.core.symbol.Symbol'>
        """
        import sympy
        return sympy.symbols(repr(ex))

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter
            sage: s = SympyConverter()
            sage: f = sin(2)
            sage: s.composition(f, f.operator())
            sin(2)
            sage: type(_)
            sin
            sage: f = arcsin(2)
            sage: s.composition(f, f.operator())
            asin(2)
        """
        f = operator._sympy_init_()
        g = ex.operands()
        import sympy

        f_sympy = getattr(sympy, f, None)
        if f_sympy:
            return f_sympy(*sympy.sympify(g, evaluate=False))
        else:
            raise NotImplementedError("SymPy function '%s' doesn't exist" % f)

sympy = SympyConverter()

#############
# Algebraic #
#############
class AlgebraicConverter(Converter):
    def __init__(self, field):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: a.field
            Algebraic Field
            sage: a.reciprocal_trig_functions['cot']
            tan
        """
        self.field = field

        from sage.functions.all import reciprocal_trig_functions
        self.reciprocal_trig_functions = reciprocal_trig_functions

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: f = SR(2)
            sage: a.pyobject(f, f.pyobject())
            2
            sage: _.parent()
            Algebraic Field
        """
        return self.field(obj)

    def arithmetic(self, ex, operator):
        """
        Convert a symbolic expression to an algebraic number.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: f = 2^(1/2)
            sage: a = AlgebraicConverter(QQbar)
            sage: a.arithmetic(f, f.operator())
            1.414213562373095?

        TESTS::

            sage: f = pi^6
            sage: a = AlgebraicConverter(QQbar)
            sage: a.arithmetic(f, f.operator())
            Traceback (most recent call last):
            ...
            TypeError: unable to convert pi^6 to Algebraic Field
        """
        # We try to avoid simplifying, because maxima's simplify command
        # can change the value of a radical expression (by changing which
        # root is selected).
        try:
            if operator is _operator.pow:
                from sage.rings.all import Rational
                base, expt = ex.operands()
                base = self.field(base)
                expt = Rational(expt)
                return self.field(base**expt)
            else:
                return reduce(operator, map(self, ex.operands()))
        except TypeError:
            pass

        if operator is _operator.pow:
            from sage.symbolic.constants import e, pi, I
            base, expt = ex.operands()
            if base == e and expt / (pi*I) in QQ:
                return exp(expt)._algebraic_(self.field)

        raise TypeError("unable to convert %s to %s"%(ex, self.field))

    def composition(self, ex, operator):
        """
        Coerce to an algebraic number.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: a.composition(exp(I*pi/3), exp)
            0.500000000000000? + 0.866025403784439?*I
            sage: a.composition(sin(pi/5), sin)
            0.5877852522924731? + 0.?e-18*I

        TESTS::

            sage: QQbar(zeta(7))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert zeta(7) to Algebraic Field
        """
        func = operator
        operand, = ex.operands()

        QQbar = self.field.algebraic_closure()
        # Note that comparing functions themselves goes via maxima, and is SLOW
        func_name = repr(func)
        if func_name == 'exp':
            rat_arg = (operand.imag()/(2*ex.parent().pi()))._rational_()
            if rat_arg == 0:
                # here we will either try and simplify, or return
                raise ValueError("Unable to represent as an algebraic number.")
            real = operand.real()
            if real:
                mag = exp(operand.real())._algebraic_(QQbar)
            else:
                mag = 1
            res = mag * QQbar.zeta(rat_arg.denom())**rat_arg.numer()
        elif func_name in ['sin', 'cos', 'tan']:
            exp_ia = exp(SR(-1).sqrt()*operand)._algebraic_(QQbar)
            if func_name == 'sin':
                res = (exp_ia - ~exp_ia)/(2*QQbar.zeta(4))
            elif func_name == 'cos':
                res = (exp_ia + ~exp_ia)/2
            else:
                res = -QQbar.zeta(4)*(exp_ia - ~exp_ia)/(exp_ia + ~exp_ia)
        elif func_name in ['sinh', 'cosh', 'tanh']:
            exp_a = exp(operand)._algebraic_(QQbar)
            if func_name == 'sinh':
                res = (exp_a - ~exp_a)/2
            elif func_name == 'cosh':
                res = (exp_a + ~exp_a)/2
            else:
                res = (exp_a - ~exp_a) / (exp_a + ~exp_a)
        elif func_name in self.reciprocal_trig_functions:
            res = ~self.reciprocal_trig_functions[func_name](operand)._algebraic_(QQbar)
        else:
            res = func(operand._algebraic_(self.field))
            #We have to handle the case where we get the same symbolic
            #expression back.  For example, QQbar(zeta(7)).  See
            #ticket #12665.
            if cmp(res, ex) == 0:
                raise TypeError("unable to convert %s to %s"%(ex, self.field))
        return self.field(res)

def algebraic(ex, field):
    """
    Returns the symbolic expression *ex* as a element of the algebraic
    field *field*.

    EXAMPLES::

        sage: a = SR(5/6)
        sage: AA(a)
        5/6
        sage: type(AA(a))
        <class 'sage.rings.qqbar.AlgebraicReal'>
        sage: QQbar(a)
        5/6
        sage: type(QQbar(a))
        <class 'sage.rings.qqbar.AlgebraicNumber'>
        sage: QQbar(i)
        1*I
        sage: AA(golden_ratio)
        1.618033988749895?
        sage: QQbar(golden_ratio)
        1.618033988749895?
        sage: QQbar(sin(pi/3))
        0.866025403784439?

        sage: QQbar(sqrt(2) + sqrt(8))
        4.242640687119285?
        sage: AA(sqrt(2) ^ 4) == 4
        True
        sage: AA(-golden_ratio)
        -1.618033988749895?
        sage: QQbar((2*I)^(1/2))
        1 + 1*I
        sage: QQbar(e^(pi*I/3))
        0.500000000000000? + 0.866025403784439?*I

        sage: AA(x*sin(0))
        0
        sage: QQbar(x*sin(0))
        0
    """
    return AlgebraicConverter(field)(ex)

##############
# Polynomial #
##############
class PolynomialConverter(Converter):
    def __init__(self, ex, base_ring=None, ring=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: x, y = var('x,y')
            sage: p = PolynomialConverter(x+y, base_ring=QQ)
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Multivariate Polynomial Ring in x, y over Rational Field

            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Univariate Polynomial Ring in x over Rational Field

            sage: p = PolynomialConverter(x, ring=QQ['x,y'])
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Multivariate Polynomial Ring in x, y over Rational Field

            sage: p = PolynomialConverter(x+y, ring=QQ['x'])
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Univariate Polynomial Ring in x over Rational Field


        """
        if not (ring is None or base_ring is None):
            raise TypeError("either base_ring or ring must be specified, but not both")
        self.ex = ex

        if ring is not None:
            base_ring = ring.base_ring()
            self.varnames = ring.variable_names_recursive()
            for v in ex.variables():
                if repr(v) not in self.varnames and v not in base_ring:
                    raise TypeError("%s is not a variable of %s" %(v, ring))
            self.ring = ring
            self.base_ring = base_ring
        elif base_ring is not None:
            self.base_ring = base_ring
            vars = self.ex.variables()
            if len(vars) == 0:
                vars = ['x']
            from sage.rings.all import PolynomialRing
            self.ring = PolynomialRing(self.base_ring, names=vars)
            self.varnames = self.ring.variable_names()
        else:
            raise TypeError("either a ring or base ring must be specified")

    def symbol(self, ex):
        """
        Returns a variable in the polynomial ring.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: p.symbol(x)
            x
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.ring(repr(ex))

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: f = SR(2)
            sage: p.pyobject(f, f.pyobject())
            2
            sage: _.parent()
            Rational Field
        """
        return self.base_ring(obj)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: a = sin(2)
            sage: p = PolynomialConverter(a*x, base_ring=RR)
            sage: p.composition(a, a.operator())
            0.909297426825682
        """
        return self.base_ring(ex)

    def relation(self, ex, op):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import PolynomialConverter

            sage: x, y = var('x, y')
            sage: p = PolynomialConverter(x, base_ring=RR)

            sage: p.relation(x==3, operator.eq)
            x - 3.00000000000000
            sage: p.relation(x==3, operator.lt)
            Traceback (most recent call last):
            ...
            ValueError: Unable to represent as a polynomial

            sage: p = PolynomialConverter(x - y, base_ring=QQ)
            sage: p.relation(x^2 - y^3 + 1 == x^3, operator.eq)
            -x^3 - y^3 + x^2 + 1
        """
        import operator
        if op == operator.eq:
            return self(ex.lhs()) - self(ex.rhs())
        else:
            raise ValueError("Unable to represent as a polynomial")

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import PolynomialConverter

            sage: x, y = var('x, y')
            sage: p = PolynomialConverter(x, base_ring=RR)
            sage: p.arithmetic(pi+e, operator.add)
            5.85987448204884
            sage: p.arithmetic(x^2, operator.pow)
            x^2

            sage: p = PolynomialConverter(x+y, base_ring=RR)
            sage: p.arithmetic(x*y+y^2, operator.add)
            x*y + y^2

            sage: p = PolynomialConverter(y^(3/2), ring=SR['x'])
            sage: p.arithmetic(y^(3/2), operator.pow)
            y^(3/2)
            sage: _.parent()
            Symbolic Ring
        """
        if not any(repr(v) in self.varnames for v in ex.variables()):
            return self.base_ring(ex)
        elif operator == _operator.pow:
            from sage.rings.all import Integer
            base, exp = ex.operands()
            return self(base)**Integer(exp)
        else:
            ops = [self(a) for a in ex.operands()]
            return reduce(operator, ops)

def polynomial(ex, base_ring=None, ring=None):
    """
    Returns a polynomial from the symbolic expression *ex*.  Either a
    base ring *base_ring* or a polynomial ring *ring* can be specified
    for the parent of result.  If just a base ring is given, then the variables
    of the base ring will be the variables of the expression *ex*.

    EXAMPLES::

         sage: from sage.symbolic.expression_conversions import polynomial
         sage: f = x^2 + 2
         sage: polynomial(f, base_ring=QQ)
         x^2 + 2
         sage: _.parent()
         Univariate Polynomial Ring in x over Rational Field

         sage: polynomial(f, ring=QQ['x,y'])
         x^2 + 2
         sage: _.parent()
         Multivariate Polynomial Ring in x, y over Rational Field

         sage: x, y = var('x, y')
         sage: polynomial(x + y^2, ring=QQ['x,y'])
         y^2 + x
         sage: _.parent()
         Multivariate Polynomial Ring in x, y over Rational Field

         sage: s,t=var('s,t')
         sage: expr=t^2-2*s*t+1
         sage: expr.polynomial(None,ring=SR['t'])
         t^2 - 2*s*t + 1
         sage: _.parent()
         Univariate Polynomial Ring in t over Symbolic Ring

         sage: polynomial(y - sqrt(x), ring=SR[y])
         y - sqrt(x)
         sage: _.list()
         [-sqrt(x), 1]

    The polynomials can have arbitrary (constant) coefficients so long as
    they coerce into the base ring::

         sage: polynomial(2^sin(2)*x^2 + exp(3), base_ring=RR)
         1.87813065119873*x^2 + 20.0855369231877
    """
    converter = PolynomialConverter(ex, base_ring=base_ring, ring=ring)
    res = converter()
    return converter.ring(res)

##############
# Fast Float #
##############

class FastFloatConverter(Converter):
    def __init__(self, ex, *vars):
        """
        Returns an object which provides fast floating point
        evaluation of the symbolic expression *ex*.  This is an class
        used internally and is not meant to be used directly.

        See :mod:`sage.ext.fast_eval` for more information.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: f = 1 + sin(x)/x + sqrt(z^2+y^2)/cosh(x)
            sage: ff = f._fast_float_('x', 'y', 'z')
            sage: f(x=1.0,y=2.0,z=3.0).n()
            4.1780638977...
            sage: ff(1.0,2.0,3.0)
            4.1780638977...

        Using _fast_float_ without specifying the variable names is
        deprecated::

            sage: f = x._fast_float_()
            doctest:...: DeprecationWarning: Substitution using
            function-call syntax and unnamed arguments is deprecated
            and will be removed from a future release of Sage; you
            can use named arguments instead, like EXPR(x=..., y=...)
            See http://trac.sagemath.org/5930 for details.
            sage: f(1.2)
            1.2

        Using _fast_float_ on a function which is the identity is
        now supported (see Trac 10246)::

            sage: f = symbolic_expression(x).function(x)
            sage: f._fast_float_(x)
            <sage.ext.fast_eval.FastDoubleFunc object at ...>
            sage: f(22)
            22
        """
        self.ex = ex

        if vars == ():
            try:
                vars = ex.arguments()
            except AttributeError:
                vars = ex.variables()

            if vars:
                from sage.misc.superseded import deprecation
                deprecation(5930, "Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)")


        self.vars = vars

        import sage.ext.fast_eval as fast_float
        self.ff = fast_float

        Converter.__init__(self, use_fake_div=True)

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: ff = fast_float(x == 2, 'x')
            sage: ff(2)
            0.0
            sage: ff(4)
            2.0
            sage: ff = fast_float(x < 2, 'x')
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if operator is not _operator.eq:
            raise NotImplementedError
        return self(ex.lhs() - ex.rhs())

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: f = SR(2)._fast_float_()
            sage: f(3)
            2.0
        """
        try:
            return obj._fast_float_(*self.vars)
        except AttributeError:
            return self.ff.fast_float_constant(float(obj))

    def symbol(self, ex):
        r"""
        EXAMPLES::

            sage: f = x._fast_float_('x', 'y')
            sage: f(1,2)
            1.0
            sage: f = x._fast_float_('y', 'x')
            sage: f(1,2)
            2.0
        """
        if self.vars == ():
            return self.ff.fast_float_arg(0)

        vars = list(self.vars)
        name = repr(ex)
        if name in vars:
            return self.ff.fast_float_arg(vars.index(name))
        svars = [repr(x) for x in vars]
        if name in svars:
            return self.ff.fast_float_arg(svars.index(name))

        if ex.is_symbol(): # case of callable function which is the variable, like f(x)=x
            name = repr(SR(ex)) # this gets back just the 'output' of the function
            if name in svars:
                return self.ff.fast_float_arg(svars.index(name))

        try:
            return self.ff.fast_float_constant(float(ex))
        except TypeError:
            raise ValueError("free variable: %s" % repr(ex))

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = x*x-y
            sage: ff = f._fast_float_('x','y')
            sage: ff(2,3)
            1.0

            sage: a = x + 2*y
            sage: f = a._fast_float_('x', 'y')
            sage: f(1,0)
            1.0
            sage: f(0,1)
            2.0

            sage: f = sqrt(x)._fast_float_('x'); f.op_list()
            ['load 0', 'call sqrt(1)']

            sage: f = (1/2*x)._fast_float_('x'); f.op_list()
            ['load 0', 'push 0.5', 'mul']
        """
        operands = ex.operands()
        if operator is _operator.neg:
            return operator(self(operands[0]))

        from sage.rings.all import Rational
        if operator is _operator.pow and operands[1] == Rational(((1,2))):
            from sage.functions.all import sqrt
            return sqrt(self(operands[0]))
        fops = map(self, operands)
        return reduce(operator, fops)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: f = sqrt(x)._fast_float_('x')
            sage: f(2)
            1.41421356237309...
            sage: y = var('y')
            sage: f = sqrt(x+y)._fast_float_('x', 'y')
            sage: f(1,1)
            1.41421356237309...

        ::

            sage: f = sqrt(x+2*y)._fast_float_('x', 'y')
            sage: f(2,0)
            1.41421356237309...
            sage: f(0,1)
            1.41421356237309...
        """
        f = operator
        g = map(self, ex.operands())
        try:
            return f(*g)
        except TypeError:
            from sage.functions.other import abs_symbolic
            if f is abs_symbolic:
                return abs(*g) # special case
            else:
                return self.ff.fast_float_func(f, *g)

def fast_float(ex, *vars):
    """
    Returns an object which provides fast floating point evaluation of
    the symbolic expression *ex*.

    See :mod:`sage.ext.fast_eval` for more information.

    EXAMPLES::

        sage: from sage.symbolic.expression_conversions import fast_float
        sage: f = sqrt(x+1)
        sage: ff = fast_float(f, 'x')
        sage: ff(1.0)
        1.4142135623730951
    """
    return FastFloatConverter(ex, *vars)()

#################
# Fast Callable #
#################

class FastCallableConverter(Converter):
    def __init__(self, ex, etb):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FastCallableConverter
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x'])
            sage: f = FastCallableConverter(x+2, etb)
            sage: f.ex
            x + 2
            sage: f.etb
            <sage.ext.fast_callable.ExpressionTreeBuilder object at 0x...>
            sage: f.use_fake_div
            True
        """
        self.ex = ex
        self.etb = etb
        Converter.__init__(self, use_fake_div=True)

    def pyobject(self, ex, obj):
        r"""
        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x'])
            sage: pi._fast_callable_(etb)
            pi
            sage: etb = ExpressionTreeBuilder(vars=['x'], domain=RDF)
            sage: pi._fast_callable_(etb)
            3.14159265359
        """
        from sage.symbolic.constants import Constant
        if isinstance(obj, Constant):
            obj = obj.expression()
        return self.etb.constant(obj)

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: ff = fast_callable(x == 2, vars=['x'])
            sage: ff(2)
            0
            sage: ff(4)
            2
            sage: ff = fast_callable(x < 2, vars=['x'])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if operator is not _operator.eq:
            raise NotImplementedError
        return self(ex.lhs() - ex.rhs())

    def arithmetic(self, ex, operator):
        r"""
        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: var('x,y')
            (x, y)
            sage: (x+y)._fast_callable_(etb)
            add(v_0, v_1)
            sage: (-x)._fast_callable_(etb)
            neg(v_0)
            sage: (x+y+x^2)._fast_callable_(etb)
            add(add(ipow(v_0, 2), v_0), v_1)

        TESTS:

        Check if rational functions with numerator 1 can be converted. #8056::

            sage: (1/pi/x)._fast_callable_(etb)
            div(1, mul(pi, v_0))

            sage: etb = ExpressionTreeBuilder(vars=['x'], domain=RDF)
            sage: (x^7)._fast_callable_(etb)
            ipow(v_0, 7)
            sage: f(x)=1/pi/x; plot(f,2,3)
        """
        # This used to convert the operands first.  Doing it this way
        # instead gives a chance to notice powers with an integer
        # exponent before the exponent gets (potentially) converted
        # to another type.
        operands = ex.operands()
        if operator is _operator.pow:
            exponent = operands[1]
            if exponent == -1:
                return self.etb.call(_operator.div, 1, operands[0])
            elif exponent == 0.5:
                from sage.functions.all import sqrt
                return self.etb.call(sqrt, operands[0])
            elif exponent == -0.5:
                from sage.functions.all import sqrt
                return self.etb.call(_operator.div, 1, self.etb.call(sqrt, operands[0]))
        elif operator is _operator.neg:
            return self.etb.call(operator, operands[0])
        return reduce(lambda x,y: self.etb.call(operator, x,y), operands)

    def symbol(self, ex):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x, y, z = var('x,y,z')
            sage: x._fast_callable_(etb)
            v_0
            sage: y._fast_callable_(etb)
            v_1
            sage: z._fast_callable_(etb)
            Traceback (most recent call last):
            ...
            ValueError: Variable 'z' not found
        """
        return self.etb.var(SR(ex))

    def composition(self, ex, function):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x,y = var('x,y')
            sage: sin(sqrt(x+y))._fast_callable_(etb)
            sin(sqrt(add(v_0, v_1)))
            sage: arctan2(x,y)._fast_callable_(etb)
            {arctan2}(v_0, v_1)
        """
        return self.etb.call(function, *ex.operands())


def fast_callable(ex, etb):
    """
    Given an ExpressionTreeBuilder *etb*, return an Expression representing
    the symbolic expression *ex*.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=['x','y'])
        sage: x,y = var('x,y')
        sage: f = y+2*x^2
        sage: f._fast_callable_(etb)
        add(mul(ipow(v_0, 2), 2), v_1)

        sage: f = (2*x^3+2*x-1)/((x-2)*(x+1))
        sage: f._fast_callable_(etb)
        div(add(add(mul(ipow(v_0, 3), 2), mul(v_0, 2)), -1), mul(add(v_0, 1), add(v_0, -2)))

    """
    return FastCallableConverter(ex, etb)()

class RingConverter(Converter):
    def __init__(self, R, subs_dict=None):
        """
        A class to convert expressions to other rings.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF, subs_dict={x:2})
            sage: R.ring
            Real Interval Field with 53 bits of precision
            sage: R.subs_dict
            {x: 2}
            sage: R(pi+e)
            5.85987448204884?
            sage: loads(dumps(R))
            <sage.symbolic.expression_conversions.RingConverter object at 0x...>
        """
        self.subs_dict = {} if subs_dict is None else subs_dict
        self.ring = R

    def symbol(self, ex):
        """
        All symbols appearing in the expression must appear in *subs_dict*
        in order for the conversion to be successful.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF, subs_dict={x:2})
            sage: R(x+pi)
            5.141592653589794?

            sage: R = RingConverter(RIF)
            sage: R(x+pi)
            Traceback (most recent call last):
            ...
            TypeError
        """
        try:
            return self.ring(self.subs_dict[ex])
        except KeyError:
            raise TypeError

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF)
            sage: R(SR(5/2))
            2.5000000000000000?
        """
        return self.ring(obj)

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: P.<z> = ZZ[]
            sage: R = RingConverter(P, subs_dict={x:z})
            sage: a = 2*x^2 + x + 3
            sage: R(a)
            2*z^2 + z + 3
        """
        if operator not in [_operator.add, _operator.mul, _operator.pow]:
            raise TypeError

        operands = ex.operands()
        if operator is _operator.pow:
            from sage.all import Integer, Rational
            base, expt = operands

            if expt == Rational(((1,2))):
                from sage.functions.all import sqrt
                return sqrt(self(base))
            try:
                expt = Integer(expt)
            except TypeError:
                pass

            base = self(base)
            return base ** expt
        else:
            return reduce(operator, map(self, operands))

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF)
            sage: R(cos(2))
            -0.4161468365471424?
        """
        res =  operator(*map(self, ex.operands()))
        if res.parent() is not self.ring:
            raise TypeError
        else:
            return res

class SubstituteFunction(Converter):
    def __init__(self, ex, original, new):
        """
        A class that walks the tree and replaces occurrences of a
        function with another.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: s(1/foo(foo(x)) + foo(2))
            1/bar(bar(x)) + bar(2)
        """
        self.original = original
        self.new = new
        self.ex = ex

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: s.symbol(x)
            x
        """
        return ex

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: f = SR(2)
            sage: s.pyobject(f, f.pyobject())
            2
            sage: _.parent()
            Symbolic Ring
        """
        return ex

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: eq = foo(x) == x
            sage: s.relation(eq, eq.operator())
            bar(x) == x
        """
        return operator(self(ex.lhs()), self(ex.rhs()))

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: f = x*foo(x) + pi/foo(x)
            sage: s.arithmetic(f, f.operator())
            x*bar(x) + pi/bar(x)
        """
        return reduce(operator, map(self, ex.operands()))

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: f = foo(x)
            sage: s.composition(f, f.operator())
            bar(x)
            sage: f = foo(foo(x))
            sage: s.composition(f, f.operator())
            bar(bar(x))
            sage: f = sin(foo(x))
            sage: s.composition(f, f.operator())
            sin(bar(x))
            sage: f = foo(sin(x))
            sage: s.composition(f, f.operator())
            bar(sin(x))
        """
        if operator == self.original:
            return self.new(*map(self, ex.operands()))
        else:
            return operator(*map(self, ex.operands()))

    def derivative(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: f = foo(x).diff(x)
            sage: s.derivative(f, f.operator())
            D[0](bar)(x)

        TESTS:

        We can substitute functions under a derivative operator,
        :trac:`12801`::

            sage: f = function('f')
            sage: g = function('g')
            sage: f(g(x)).diff(x).substitute_function(g, sin)
            cos(x)*D[0](f)(sin(x))

        """
        if operator.function() == self.original:
            return operator.change_function(self.new)(*map(self,ex.operands()))
        else:
            return operator(*map(self, ex.operands()))
