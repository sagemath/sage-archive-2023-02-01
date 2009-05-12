r"""
Symbolic Equations and Inequalities.

Sage can solve symbolic equations and express inequalities. For
example, we derive the quadratic formula as follows::

    sage: a,b,c = var('a,b,c')
    sage: qe = (a*x^2 + b*x + c == 0)
    sage: print qe
                                     2
                                  a x  + b x + c == 0
    sage: print solve(qe, x)
    [
                                          2
                                  - sqrt(b  - 4 a c) - b
                              x == ----------------------
                                           2 a,
                                         2
                                   sqrt(b  - 4 a c) - b
                               x == --------------------
                                           2 a
    ]

AUTHORS:

- Bobby Moretti: initial version (based on a trick that Robert
  Bradshaw suggested).

- William Stein: second version

- William Stein (2007-07-16): added arithmetic with symbolic equations

EXAMPLES::

    sage: x,y,a = var('x,y,a')
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(1 - y^2), x == sqrt(1 - y^2)]

::

    sage: f = x^5 + a
    sage: solve(f==0,x)
    [x == e^(2*I*pi/5)*(-a)^(1/5), x == e^(4*I*pi/5)*(-a)^(1/5), x == e^(-(4*I*pi/5))*(-a)^(1/5), x == e^(-(2*I*pi/5))*(-a)^(1/5), x == (-a)^(1/5)]

You can also do arithmetic with inequalities, as illustrated
below::

    sage: var('x y')
    (x, y)
    sage: f = x + 3 == y - 2
    sage: f
    x + 3 == y - 2
    sage: g = f - 3; g
    x == y - 5
    sage: h =  x^3 + sqrt(2) == x*y*sin(x)
    sage: h
    x^3 + sqrt(2) == x*sin(x)*y
    sage: h - sqrt(2)
    x^3 == x*sin(x)*y - sqrt(2)
    sage: h + f
    x^3 + x + sqrt(2) + 3 == x*sin(x)*y + y - 2
    sage: f = x + 3 < y - 2
    sage: g = 2 < x+10
    sage: f - g
    x + 1 < y - x - 12
    sage: f + g
    x + 5 < y + x + 8
    sage: f*(-1)
    -x - 3 > 2 - y

TESTS: We test serializing symbolic equations::

    sage: eqn = x^3 + 2/3 >= x
    sage: loads(dumps(eqn))
    x^3 + 2/3 >= x
    sage: loads(dumps(eqn)) == eqn
    True
"""

_assumptions = []

from sage.structure.sage_object import SageObject
from sage.structure.sequence    import Sequence

from calculus                   import maxima

import operator

symbols = {operator.lt:' < ', operator.le:' <= ', operator.eq:' == ',
           operator.ne:' != ',
            operator.ge:' >= ', operator.gt:' > '}

maxima_symbols = dict(symbols)
maxima_symbols[operator.eq] = '='
maxima_symbols[operator.ne] = '#'


latex_symbols = {operator.lt:' < ', operator.le:' \\leq ', operator.eq:' = ',
                 operator.ne:' \\neq ',
            operator.ge:' \\geq ', operator.gt:' > '}

comparisons = {operator.lt:set([-1]), operator.le:set([-1,0]), operator.eq:set([0]),
               operator.ne:set([-1,1]),  operator.ge:set([0,1]), operator.gt:set([1])}

opposite_op = {operator.lt:operator.gt, operator.le:operator.ge,
               operator.ne:operator.ne, operator.eq:operator.eq,
               operator.gt: operator.lt, operator.ge:operator.le}

def var_cmp(x,y):
    """
    Return comparison of the two variables x and y, which is just the
    comparison of the underlying string representations of the
    variables. This is used internally by the Calculus package.

    INPUT:


    -  ``x, y`` - symbolic variables


    OUTPUT: Python integer; either -1, 0, or 1.

    EXAMPLES::

        sage: sage.calculus.equations.var_cmp(x,x)
        0
        sage: sage.calculus.equations.var_cmp(x,var('z'))
        -1
        sage: sage.calculus.equations.var_cmp(x,var('a'))
        1
    """
    return cmp(str(x), str(y))

def is_SymbolicEquation(x):
    r"""
    Return True if x is a symbolic equation.

    EXAMPLES: The following two examples are symbolic equations::

        sage: from sage.calculus.equations import is_SymbolicEquation
        sage: is_SymbolicEquation(sin(x) == x)
        True
        sage: is_SymbolicEquation(sin(x) < x)
        True

    This is not, since ``2==3`` evaluates to the boolean
    ``False``::

        sage: is_SymbolicEquation(2 == 3)
        False

    However here since both 2 and 3 are coerced to be symbolic, we
    obtain a symbolic equation::

        sage: is_SymbolicEquation(SR(2) == SR(3))
        True
    """
    return isinstance(x, SymbolicEquation)

class SymbolicEquation(SageObject):
    """
    A symbolic equation, which consists of a left hand side, an
    operator and a right hand side.

    EXAMPLES:
    """
    def __init__(self, left, right, op):
        r"""
        Create a symbolic expression.

        Internally a symbolic expression is simply a left side
        (``self._left``), operator
        (``self._op``), and a right hand side
        (``self._right``), where the left and right hand sides
        are symbolic expressions and the operator is a Python equation or
        inequality operator, e.g., ``operator.le``.

        EXAMPLES:

        One should not call the SymbolicEquation constructor directly,
        since it does no type checking. However, we illustrate how to do so
        below.

        Bad illustrative usage::

            sage: eqn = sage.calculus.equations.SymbolicEquation(-x, x^2 + 1, operator.gt); eqn
            -x > x^2 + 1

        Really bad usage!

        ::

            sage: eqn.__init__(x, 2*x+pi, operator.lt)
            sage: eqn                 # cripes!
            x < 2*x + pi
        """
        self._left = left
        self._right = right
        self._op = op

    def __call__(self, *args, **argv):
        """
        Substitute both sides of this equation

        This is very slow currently since we piggy-back off of the symbolic
        matrix functionality.

        EXAMPLES::

            sage: var('theta')
            theta
            sage: eqn =   (x^3 + theta < sin(x*theta))
            sage: eqn(x = 5)
            theta + 125 < sin(5*theta)
            sage: eqn(theta=x, x=0)
            x < 0
            sage: var('y')
            y
            sage: eqn = x^3 < sin(y)
            sage: eqn(2)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            8 < sin(y)
            sage: eqn(x=2)
            8 < sin(y)
            sage: eqn(2,3)
            8 < sin(3)
            sage: eqn(x=2,y=3)
            8 < sin(3)
            sage: eqn = x^3 < 2
            sage: eqn(2)
            8 < 2
        """
        from sage.matrix.all import matrix
        from sage.calculus.all import SR
        m = matrix(SR, 1, 2, [self._left, self._right])
        left,right = m(*args, **argv)[0]
        return self._op(left, right)


    def _maple_(self, maple=None):
        """
        Returns a Maple version of self.

        EXAMPLES::

            sage: eq = x == 2
            sage: maple(eq)   #optional
            x = 2
        """
        if maple is None:
            from sage.interfaces.maple import maple
        lhs = maple(self.lhs())
        rhs = maple(self.rhs())
        return maple("%s = %s"%(lhs.name(), rhs.name()))

    def __getitem__(self, i):
        """
        Return the ith part of this equation:

        OUTPUT:


        -  ``self[0]`` - left hand side

        -  ``self[1]`` - operator

        -  ``self[2]`` - right hand side


        EXAMPLES::

            sage: eqn = x^2 + sin(x) < cos(x^2)
            sage: eqn[0]
            sin(x) + x^2
            sage: eqn[1]
            <built-in function lt>
            sage: eqn[2]
            cos(x^2)
            sage: eqn[-1]
            cos(x^2)
        """
        return [self._left, self._op, self._right][i]

    def _scalar(self, scalar, op, checksign=True):
        """
        INPUT:


        -  ``scalar`` - number

        -  ``op`` - operation to perform

        -  ``checksign`` - (default: True) boolean; if True and
           op is multiply or divides, switch direction of inequality of x is
           negative; otherwise direction will not switch.


        EXAMPLES::

            sage: var('x y')
            (x, y)
            sage: f = x + 3 < y - 2
            sage: f*-1
            -x - 3 > 2 - y
            sage: f._scalar(-1, operator.mul, checksign=True)
            -x - 3 > 2 - y
            sage: f._scalar(-1, operator.mul, checksign=False)
            -x - 3 < 2 - y
            sage: f * 5
            5*(x + 3) < 5*(y - 2)
            sage: f - 3
            x < y - 5
            sage: f + 2
            x + 5 < y
        """
        SR = self._left.parent()
        x = SR(scalar)

        # There are no subtleties if both sides are equal or
        # we are adding or subtracting from both sides.
        if not checksign or (self._op == operator.eq or op in [operator.add, operator.sub]):
            return SymbolicEquation(op(self._left, x),
                                    op(self._right, x),
                                    self._op)

        # Now we are multiplying both sides by a scalar and we have an inequality, i.e.,
        # one of <, <=, >=, > or !=.
        if x == 0:
            return SymbolicEquation(SR(0), SR(0), operator.eq)
        elif x > 0:
            return SymbolicEquation(op(self._left, x), op(self._right, x), self._op)
        elif x < 0:
            # Now we are multiplying or dividing both sides by a negative number.
            op2 = opposite_op[self._op]
            return SymbolicEquation(op(self._left, x), op(self._right, x), op2)
        else:
            raise ValueError, "unable to multiply or divide both sides of an inequality by a number whose sign can't be determined."

    def multiply_both_sides(self, x):
        """
        Multiply both sides of this inequality by `x`.

        EXAMPLES::

            sage: var('x,y'); f = x + 3 < y - 2
            (x, y)
            sage: f.multiply_both_sides(7)
            7*(x + 3) < 7*(y - 2)
            sage: f.multiply_both_sides(-1/2)
            (-x - 3)/2 > (2 - y)/2
            sage: f*(-2/3)
            -2*(x + 3)/3 > -2*(y - 2)/3
            sage: f*(-pi)
            -1*pi*(x + 3) > -1*pi*(y - 2)
            sage: f*(1+I)
            Traceback (most recent call last):
            ...
            ValueError: unable to multiply or divide both sides of an inequality by a number whose sign can't be determined.

        Multiplying by complex numbers works only if it's an equality::

            sage: f = sqrt(2) + x == y^3
            sage: f.multiply_both_sides(I)
            I*(x + sqrt(2)) == I*y^3
            sage: f.multiply_both_sides(-1)
            -x - sqrt(2) == -y^3

        Some further examples::

            sage: (x^3 + 1 > 2*sqrt(3)) * (-1)
            -x^3 - 1 < -2*sqrt(3)
            sage: (x^3 + 1 >= 2*sqrt(3)) * (-1)
            -x^3 - 1 <= -2*sqrt(3)
            sage: (x^3 + 1 <= 2*sqrt(3)) * (-1)
            -x^3 - 1 >= -2*sqrt(3)
        """
        return self._scalar(x, operator.mul)

    def divide_both_sides(self, x, checksign=True):
        """
        Divide both sides of the inequality by `x`.

        INPUT:


        -  ``x`` - number

        -  ``checksign`` - (default: True) boolean; if True,
           switch direction of inequality of x is negative; otherwise
           direction will not switch.


        EXAMPLES::

            sage: var('theta')
            theta
            sage: eqn =   (x^3 + theta < sin(x*theta))
            sage: eqn.divide_both_sides(theta, checksign=False)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: assume(theta > 0)
            sage: eqn.divide_both_sides(theta)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: eqn/theta
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: forget(theta > 0)
            sage: eqn.divide_both_sides(theta)
            Traceback (most recent call last):
            ...
            ValueError: unable to multiply or divide both sides of an inequality by a number whose sign can't be determined.

        As a shorthand you can just use the divides notation::

            sage: (x^3 + 1 > x^2 - 1) / (-1)
            -x^3 - 1 < 1 - x^2

        The quantity `x^2 - 1` could be either negative or positive
        depending on `x`, so dividing by it is not defined.

        ::

            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            Traceback (most recent call last):
            ...
            ValueError: unable to multiply or divide both sides of an inequality by a number whose sign can't be determined.

        If we specify that `x^2 - 1> 0`, then dividing is defined.

        ::

            sage: assume(x^2 - 1 > 0)
            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            (x^3 + 1)/(x^2 - 1) > 1
            sage: forget()

        We can also specify that `x^2 - 1 < 0`. Note that now the
        inequality direction changes.

        ::

            sage: assume(x^2 - 1 < 0)
            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            (x^3 + 1)/(x^2 - 1) < 1
            sage: forget()
        """
        return self._scalar(x, operator.div, checksign=checksign)

    def add_to_both_sides(self, x):
        """
        Add `x` to both sides of this symbolic equation.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: eqn = x^2 + y^2 + z^2 <= 1
            sage: eqn.add_to_both_sides(-z^2)
            y^2 + x^2 <= 1 - z^2
            sage: eqn.add_to_both_sides(I)
            z^2 + y^2 + x^2 + I <= I + 1
        """
        return self._scalar(x, operator.add)

    def subtract_from_both_sides(self, x):
        """
        Subtract `x` from both sides of this symbolic equation.

        EXAMPLES::

            sage: eqn = x*sin(x)*sqrt(3) + sqrt(2) > cos(sin(x))
            sage: eqn.subtract_from_both_sides(sqrt(2))
            sqrt(3)*x*sin(x) > cos(sin(x)) - sqrt(2)
            sage: eqn.subtract_from_both_sides(cos(sin(x)))
            -cos(sin(x)) + sqrt(3)*x*sin(x) + sqrt(2) > 0
        """
        return self._scalar(x, operator.sub)

    def _arith(self, right, op):
        """
        This function is called internally to implement arithmetic
        operations on symbolic expressions.

        INPUT:


        -  ``self`` - a symbolic equation

        -  ``right`` - a symbolic equation

        -  ``op`` - an operation, e.g., operator.add


        EXAMPLES: We create two symbolic equations and add them::

            sage: e1 = x^3 + x < sin(2*x)
            sage: e2 = x^2 - x < cos(x)
            sage: e1._arith(e2, operator.add)
            x^3 + x^2 < sin(2*x) + cos(x)

        We try to multiply them, which doesn't really make sense::

            sage: e1._arith(e2, operator.mul)
            Traceback (most recent call last):
            ...
            ValueError: cannot multiply or divide inequalities.

        We can multiply equalities though::

            sage: e1 = x^3 + x == sin(2*x)
            sage: e2 = x^2 - x == cos(x)
            sage: f = e1._arith(e2, operator.mul); f
            (x^2 - x)*(x^3 + x) == cos(x)*sin(2*x)

        By the way, we can expand the above product by calling the
        ``expand`` method::

            sage: f.expand()
            x^5 - x^4 + x^3 - x^2 == cos(x)*sin(2*x)
        """
        if not isinstance(right, SymbolicEquation):
            return self._scalar(right, op)

        if self._op != right._op:
            raise ValueError, "can only do arithmetic with symbolic equations with the same equality or inequality"

        if not (op in [operator.add, operator.sub]):
            if not (self._op in [operator.eq, operator.ne]):
                raise ValueError, "cannot multiply or divide inequalities."

        return SymbolicEquation(op(self._left, right._left),
                                op(self._right, right._right),
                                self._op)

    def __add__(self, right):
        """
        Add two symbolic equations.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: m = 144 == -10 * a + b
            sage: n = 136 == 10 * a + b
            sage: m + n
            280 == 2*b
        """
        return self._arith(right, operator.add)

    def __radd__(self, left):
        """
        Add two symbolic equations.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: m = 144 == -10 * a + b
            sage: n = 136 == 10 * a + b
            sage: int(-144) + m
            0 == b - 10*a - 144
        """
        # ok to use this, since everything is commutative in the symbolic ring.
        return self._arith(left, operator.add)

    def __sub__(self, right):
        """
        Subtract two symbolic equations.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: m = 144 == 20 * a + b
            sage: n = 136 == 10 * a + b
            sage: m - n
            8 == 10*a
        """
        return self._arith(right, operator.sub)

    def __rsub__(self, left):
        """
        Subtract two symbolic equations.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: m = 144 == -10 * a + b
            sage: n = 136 == 10 * a + b
            sage: int(144) - m
            0 == b - 10*a - 144
        """
        # ok to use this, since everything is commutative in the symbolic ring.
        return self._arith(left, operator.sub)

    def __mul__(self, right):
        """
        Multiply two symbolic equations.

        EXAMPLES::

            sage: m = x == 5*x + 1
            sage: n = sin(x) == sin(x+2*pi)
            sage: m * n
            x*sin(x) == (5*x + 1)*sin(x)
        """
        return self._arith(right, operator.mul)

    def __rmul__(self, left):
        """
        Multiply two symbolic equations.

        ::

            sage: m = 2*x == 3*x^2 - 5
            sage: int(-1) * m
            -2*x == 5 - 3*x^2
        """
        return self._arith(left, operator.mul)

    def __div__(self, right):
        """
        Divide two symbolic equations.

        EXAMPLES::

            sage: m = x == 5*x + 1
            sage: n = sin(x) == sin(x+2*pi)
            sage: m / n
            x/sin(x) == (5*x + 1)/sin(x)
            sage: m = x != 5*x + 1
            sage: n = sin(x) != sin(x+2*pi)
            sage: m / n
            x/sin(x) != (5*x + 1)/sin(x)
        """
        return self._arith(right, operator.div)


    # The maxima one is special:
    def _maxima_(self, session=None):
        """
        Return version of this symbolic expression but in the given Maxima
        session.

        EXAMPLES::

            sage: e1 = x^3 + x == sin(2*x)
            sage: z = e1._maxima_()
            sage: z.parent() is sage.calculus.calculus.maxima
            True
            sage: z = e1._maxima_(maxima)
            sage: z.parent() is maxima
            True
            sage: z = maxima(e1)
            sage: z.parent() is maxima
            True
        """
        if session is None:
            return SageObject._maxima_(self, sage.calculus.calculus.maxima)
        else:
            return SageObject._maxima_(self, session)

    def substitute(self, *args, **kwds):
        """
        Do the given symbolic substitution to both sides of the equation.
        The notation is the same for substitute on a symbolic expression.

        EXAMPLES::

            sage: var('a')
            a
            sage: e = (x^3 + a == sin(x/a)); e
            x^3 + a == sin(x/a)
            sage: e.substitute(x=5*x)
            125*x^3 + a == sin(5*x/a)
            sage: e.substitute(a=1)
            x^3 + 1 == sin(x)
            sage: e.substitute(a=x)
            x^3 + x == sin(1)
            sage: e.substitute(a=x, x=1)
            x + 1 == sin(1/x)
            sage: e.substitute({a:x, x:1})
            x + 1 == sin(1/x)
        """
        return self.__call__(*args, **kwds)

    subs = substitute

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: (x>0) == (x>0)
            True
            sage: (x>0) == (x>1)
            False
            sage: (x>0) != (x>1)
            True
        """
        if not isinstance(right, SymbolicEquation):
            return cmp(type(self), type(right))
        c = cmp(self._op, right._op)
        if c: return c
        i = bool(self._left == right._left)
        if not i: return -1
        i = bool(self._right == right._right)
        if not i: return 1
        return 0

    def _repr_(self):
        r"""
        Return non-ASCII art string representation of this symbolic
        equation. This is called implicitly when displaying an equation
        (without using print).

        EXAMPLES: We create an inequality `f` and called the
        ``_repr_`` method on it, and note that this produces
        the same string as just displaying `f`::

            sage: f = x^3 + 1/3*x - sqrt(2) <= sin(x)
            sage: f._repr_()
            'x^3 + x/3 - sqrt(2) <= sin(x)'
            sage: f
            x^3 + x/3 - sqrt(2) <= sin(x)

        When using print the ``__str__`` method is called
        instead, which results in ASCII art::

            sage: print f
                               3   x
                              x  + - - sqrt(2) <= sin(x)
                                   3
        """
        return "%r%s%r" %(self._left, symbols[self._op], self._right)

    def variables(self):
        """
        Return the variables appearing in this symbolic equation.

        OUTPUT:


        -  ``tuple`` - tuple of the variables in this equation
           (the result of calling this is cached).


        EXAMPLES::

            sage: var('x,y,z,w')
            (x, y, z, w)
            sage: f =  (x+y+w) == (x^2 - y^2 - z^3);   f
            y + x + w == -z^3 - y^2 + x^2
            sage: f.variables()
            (w, x, y, z)
        """
        try:
            return self.__variables
        except AttributeError:
            v = list(set(list(self._left.variables()) + list(self._right.variables())))
            v.sort(var_cmp)
            v = tuple(v)
            self.__variables = v
            return v

    def operator(self):
        """
        Return the operator in this equation.

        EXAMPLES::

            sage: eqn = x^3 + 2/3 >= x - pi
            sage: eqn.operator()
            <built-in function ge>
            sage: (x^3 + 2/3 < x - pi).operator()
            <built-in function lt>
            sage: (x^3 + 2/3 == x - pi).operator()
            <built-in function eq>
        """
        return self._op

    def left(self):
        r"""
        Return the left hand side of this equation.

        EXAMPLES::

            sage: eqn = x^3 + 2/3 >= x - pi
            sage: eqn.lhs()
            x^3 + 2/3
            sage: eqn.left()
            x^3 + 2/3
            sage: eqn.left_hand_side()
            x^3 + 2/3

        SYNONYMS: ``lhs``, ``left_hand_side``
        """
        return self._left
    lhs = left
    left_hand_side = left

    def right(self):
        r"""
        Return the right hand side of this equation.

        EXAMPLES::

            sage: (x + sqrt(2) >= sqrt(3) + 5/2).right()
            sqrt(3) + 5/2
            sage: (x + sqrt(2) >= sqrt(3) + 5/2).rhs()
            sqrt(3) + 5/2
            sage: (x + sqrt(2) >= sqrt(3) + 5/2).right_hand_side()
            sqrt(3) + 5/2

        SYNONYMS: ``rhs``, ``right_hand_side``
        """
        return self._right
    rhs = right
    right_hand_side = right

    def __str__(self):
        r"""
        Return the string representation of this equation, in 2-d ASCII
        art.

        OUTPUT: string

        EXAMPLES::

            sage: f =  (x^2 - x == 0)
            sage: f
            x^2 - x == 0
            sage: print f
                                               2
                                              x  - x == 0

        Here we call ``__str__`` explicitly::

            sage: (x > 2/3).__str__()
            '                                         2\r\n                                     x > -\r\n                                         3'
        """
        s = self._maxima_().display2d(onscreen=False)
        s = s.replace('%pi','pi').replace('%i',' I').replace('%e', ' e').replace(' = ',' == ')
        return s

    def _latex_(self):
        r"""
        Return latex representation of this symbolic equation.

        This is obtained by calling the ``_latex_`` method on
        both the left and right hand sides, and typesetting the operator
        symbol correctly.

        OUTPUT:


        -  ``string`` - a string


        EXAMPLES: The output is a strig with backslashes, so prints funny::

            sage: (x^(3/5) >= pi)._latex_()
            '{x}^{\\frac{3}{5}}   \\geq  \\pi'

        Call the latex method to get an object that prints more nicely::

            sage: latex(x^(3/5) >= pi)
            {x}^{\frac{3}{5}}   \geq  \pi
        """
        return "%s %s %s" %(self._left._latex_(), latex_symbols[self._op],
                            self._right._latex_())

    def __nonzero__(self):
        """
        Return True if this (in)equality is definitely true. Return False
        if it is false or the algorithm for testing (in)equality is
        inconclusive.

        EXAMPLES::

            sage: k = var('k')
            sage: pol = 1/(k-1) - 1/k -1/k/(k-1);
            sage: bool(pol == 0)
            True
            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: bool(f == 0)
            True
            sage: bool( x == x )
            True
            sage: bool( x != x )
            False
            sage: bool( x > x )
            False
            sage: bool( x^2 > x )
            False
            sage: bool( x + 2 > x )
            True
            sage: bool( x - 2 > x )
            False
        """
        m = self._maxima_()

        #Handle some basic cases first
        if repr(m) in ['0=0']:
            return True
        elif repr(m) in ['0#0', '1#1']:
            return False

        try:
            s = m.parent()._eval_line('is (%s)'%m.name())
        except TypeError, msg:
            raise ValueError, "unable to evaluate the predicate '%s'"%repr(self)

        if s == 'true':
            return True
        elif s == 'unknown':
            return False

        if self.operator() != operator.eq:
            return False

        difference = self._left - self._right
        if repr(difference) == '0':
            return True

        #Try to apply some simplifications to see if left - right == 0
        simp_list = [difference.simplify, difference.simplify_log, difference.simplify_rational, difference.simplify_exp,difference.simplify_radical,difference.simplify_trig]
        for f in simp_list:
            try:
                if repr( f() ).strip() == "0":
                    return True
                    break
            except:
                pass
        return False

    def _maxima_init_(self, maxima=maxima, assume=False):
        """
        Return string representation for this symbolic equation in a form
        suitable for evaluation in Maxima.

        EXAMPLES::

            sage: (x^(3/5) >= pi^2 + e^i)._maxima_init_()
            '((x) ^ (3/5)) >= (((%pi) ^ (2)) + ((%e) ^ (%i)))'
            sage: (x == 0)._maxima_init_(assume=True)
            'equal(x, 0)'
            sage: (x != 0)._maxima_init_(assume=True)
            'notequal(x, 0)'
        """
        l = self._left._maxima_init_()
        r = self._right._maxima_init_()
        if assume:
            if  self._op == operator.eq:
                return 'equal(%s, %s)'%(l, r)
            elif self._op == operator.ne:
                return 'notequal(%s, %s)'%(l, r)
        return '(%s)%s(%s)' % (l, maxima_symbols[self._op], r)

    def assume(self):
        r"""
        Assume that this equation holds. This is relevant for symbolic
        integration, among other things.

        EXAMPLES: We call the assume method to assume that `x>2`::

            sage: (x > 2).assume()

        Bool returns True below if the inequality is *definitely* known to
        be True.

        ::

            sage: bool(x > 0)
            True
            sage: bool(x < 0)
            False

        This may or may not be True, so bool returns False::

            sage: bool(x > 3)
            False

        TESTS::

            sage: v,c = var('v,c')
            sage: assume(c != 0)
            sage: integral((1+v^2/c^2)^3/(1-v^2/c^2)^(3/2),v)
            -75*sqrt(c^2)*arcsin(sqrt(c^2)*v/c^2)/8 - v^5/(4*c^4*sqrt(1 - v^2/c^2)) - 17*v^3/(8*c^2*sqrt(1 - v^2/c^2)) + 83*v/(8*sqrt(1 - v^2/c^2))
        """
        if not self in _assumptions:
            m = self._maxima_init_(assume=True)
            maxima.assume(m)
            _assumptions.append(self)

    def find_root(self, *args, **kwds):
        r"""
        If this is a symbolic equality with an equals sign
        ``==`` find numerically a single root of this equation
        in a given interval. Otherwise raise a ``ValueError``.
        See the documentation for the global ``find_root``
        method for more about the options to this function.

        Note that this symbolic expression must involve at most one
        variable.

        EXAMPLES::

            sage: (x == sin(x)).find_root(-2,2)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            0.0
            sage: (x^5 + 3*x + 2 == 0).find_root(-2,2,x)
            -0.63283452024215225
            sage: (cos(x) == sin(x)).find_root(10,20)
            19.634954084936208

        We illustrate some valid error conditions::

            sage: (cos(x) != sin(x)).find_root(10,20)
            Traceback (most recent call last):
            ...
            ValueError: Symbolic equation must be an equality.
            sage: (SR(3)==SR(2)).find_root(-1,1)
            Traceback (most recent call last):
            ...
            RuntimeError: no zero in the interval, since constant expression is not 0.

        There must be at most one variable::

            sage: x, y = var('x,y')
            sage: (x == y).find_root(-2,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.
        """
        if self._op != operator.eq:
            raise ValueError, "Symbolic equation must be an equality."
        return (self._left - self._right).find_root(*args, **kwds)

    def forget(self):
        """
        Forget the given constraint.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: forget()
            sage: assume(x>0, y < 2)
            sage: assumptions()
            [x > 0, y < 2]
            sage: forget(y < 2)
            sage: assumptions()
            [x > 0]
        """
        m = self._maxima_()
        m.parent().forget(m)
        try:
            _assumptions.remove(self)
        except ValueError:
            pass

    def solve(self, x=None, multiplicities=False, solution_dict=False, explicit_solutions=False):
        """
        Symbolically solve for the given variable.

        .. warning::

           In many cases, only one solution is computed.

        INPUT:


        -  ``x`` - a SymbolicVariable object (if not given, the
           first in the expression is used)

        -  ``multiplicities`` - (default: False) if True, also
           returns the multiplicities of each solution, in order.

        -  ``solution_dict`` - (default: False) if True,
           return the solution as a dictionary rather than an equation.

        -  ``explicit_solution`` - (default: False); if True,
           require that all solutions returned be explicit (rather than
           implicit)


        OUTPUT: A list of SymbolicEquations with the variable to solve for
        on the left hand side.

        EXAMPLES::

            sage: S = solve(x^3 - 1 == 0, x)
            sage: S
            [x == (sqrt(3)*I - 1)/2, x == (-sqrt(3)*I - 1)/2, x == 1]
            sage: S[0]
            x == (sqrt(3)*I - 1)/2
            sage: S[0].right()
            (sqrt(3)*I - 1)/2
            sage: S = solve(x^3 - 1 == 0, x, solution_dict=True)
            sage: S
            [{x: (sqrt(3)*I - 1)/2}, {x: (-sqrt(3)*I - 1)/2}, {x: 1}]
            sage: z = 5
            sage: solve(z^2 == sqrt(3),z)
            Traceback (most recent call last):
            ...
            TypeError: 5 is not a valid variable.

        We illustrate finding multiplicities of solutions::

            sage: f = (x-1)^5*(x^2+1)
            sage: solve(f == 0, x)
            [x == -1*I, x == I, x == 1]
            sage: solve(f == 0, x, multiplicities=True)
            ([x == -1*I, x == I, x == 1], [1, 1, 5])
        """
        if not self._op is operator.eq:
            raise NotImplementedError, "solving only implemented for equalities"
        if x is None:
            v = self.variables()
            if len(v) == 0:
                if multiplicities:
                    return [], []
                else:
                    return []
            x = v[0]

        from sage.calculus.calculus import SymbolicVariable, SymbolicFunction, SymbolicFunctionEvaluation
        if not isinstance(x,(SymbolicVariable,SymbolicFunction,SymbolicFunctionEvaluation)):
            raise TypeError, "%s is not a valid variable."%x

        m = self._maxima_()
        P = m.parent()
        if explicit_solutions:
            P.eval('solveexplicit: true')
        s = m.solve(x).str()
        if explicit_solutions:
            P.eval('solveexplicit: false')

        X = string_to_list_of_solutions(s)
        if solution_dict==True:
            X=[dict([[sol.left(),sol.right()]]) for sol in X]

        if multiplicities:
            if len(X) == 0:
                return X, []
            else:
                return X, [int(e) for e in str(P.get('multiplicities'))[1:-1].split(',')]
        else:
            return X


    def expand(self, side=None):
        r"""
        Expands one or both sides of the equation.

        If side is not specified, then both sides of the equation are
        expanded by calling ``expand()`` on the corresponding
        ``SymbolicExpression``.

        If side is 'left' (or 'right'), then only the left (or right) side
        of the equation is expanded.

        EXAMPLES::

            sage: a = (16*x-13)/6 == (3*x+5)/2 - (4-x)/3
            sage: a.expand()
            8*x/3 - 13/6 == 11*x/6 + 7/6
            sage: a.expand('left')
            8*x/3 - 13/6 == (3*x + 5)/2 - (4 - x)/3
            sage: a.expand('right')
            (16*x - 13)/6 == 11*x/6 + 7/6
        """
        if side is None:
            return SymbolicEquation(self._left.expand(), self._right.expand(), self._op)
        elif side == 'left':
            return SymbolicEquation(self._left.expand(), self._right, self._op)
        elif side == 'right':
            return SymbolicEquation(self._left, self._right.expand(), self._op)
        else:
            raise ValueError, "side must be 'left', 'right', or None"

class GenericDeclaration(SageObject):

    def __init__(self, var, assumption):
        """
        This class represents generic assumptions, such as a variable being
        an integer or a function being increasing. It passes such
        information to maxima's declare (wrapped in a context so it is able
        to forget).

        INPUT:


        -  ``var`` - the variable about which assumptions are
           being made

        -  ``assumption`` - a maxima feature, either user
           defined or in the list given by maxima('features')


        EXAMPLES::

            sage: from sage.calculus.equations import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'integer')
            sage: decl.assume()
            sage: sin(x*pi)
            0
            sage: decl.forget()
            sage: sin(x*pi)
            sin(pi*x)

        Here is the list of acceptable features::

            sage: maxima('features')
            [integer,noninteger,even,odd,rational,irrational,real,imaginary,complex,analytic,increasing,decreasing,oddfun,evenfun,posfun,commutative,lassociative,rassociative,symmetric,antisymmetric,integervalued]
        """
        self._var = var
        self._assumption = assumption
        self._context = None

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.calculus.equations import GenericDeclaration
            sage: GenericDeclaration(x, 'foo')
            x is foo
        """
        return "%s is %s" % (self._var, self._assumption)

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.calculus.equations import GenericDeclaration as GDecl
            sage: var('y')
            y
            sage: GDecl(x, 'integer') == GDecl(x, 'integer')
            True
            sage: GDecl(x, 'integer') == GDecl(x, 'rational')
            False
            sage: GDecl(x, 'integer') == GDecl(y, 'integer')
            False
        """
        if isinstance(self, GenericDeclaration) and isinstance(other, GenericDeclaration):
            return cmp( (self._var, self._assumption),
                        (other._var, other._assumption) )
        else:
            return cmp(type(self), type(other))

    def assume(self):
        """
        TEST::

            sage: from sage.calculus.equations import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'even')
            sage: decl.assume()
            sage: cos(x*pi)
            1
            sage: decl.forget()
        """
        if self._context is None:
            # We get the list here because features may be added with time.
            valid_features = list(maxima("features"))
            if self._assumption not in [repr(x).strip() for x in list(valid_features)]:
                raise ValueError, "%s not a valid assumption, must be one of %s" % (self._assumption, valid_features)
            cur = maxima.get("context")
            self._context = maxima.newcontext(maxima._next_var_name())
            maxima.eval("declare(%s, %s)" % (self._var._name, self._assumption))
            maxima.set("context", cur)

        if not self in _assumptions:
            maxima.activate(self._context)
            _assumptions.append(self)

    def forget(self):
        """
        TEST::

            sage: from sage.calculus.equations import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'odd')
            sage: decl.assume()
            sage: cos(x*pi)
            -1
            sage: decl.forget()
            sage: cos(x*pi)
            cos(pi*x)
        """
        try:
            self = _assumptions.pop(_assumptions.index(self))
        except ValueError:
            pass
        if self._context is not None:
            maxima.deactivate(self._context)

def preprocess_assumptions(args):
    """
    Turns a list of the form (var1, var2, ..., 'property') into a
    sequence of declarations (var1 is property), (var2 is property),
    ...

    EXAMPLES::

        sage: from sage.calculus.equations import preprocess_assumptions
        sage: preprocess_assumptions([x, 'integer', x > 4])
        [x is integer, x > 4]
        sage: var('x,y')
        (x, y)
        sage: preprocess_assumptions([x, y, 'integer', x > 4, y, 'even'])
        [x is integer, y is integer, x > 4, y is even]
    """
    args = list(args)
    last = None
    for i, x in reversed(list(enumerate(args))):
        if isinstance(x, str):
            del args[i]
            last = x
        elif not hasattr(x, 'assume') and last is not None:
            args[i] = GenericDeclaration(x, last)
        else:
            last = None
    return args

def assume(*args):
    """
    Make the given assumptions.

    INPUT:


    -  ``*args`` - assumptions


    EXAMPLES::

        sage: assume(x > 0)
        sage: bool(sqrt(x^2) == x)
        True
        sage: forget()
        sage: bool(sqrt(x^2) == x)
        False

    An integer constraint::

        sage: var('n, P, r, r2')
        (n, P, r, r2)
        sage: assume(n, 'integer')
        sage: c = P*e^(r*n)
        sage: d = P*(1+r2)^n
        sage: solve(c==d,r2)
        [r2 == e^r - 1]

    ::

        sage: sin(n*pi)
        0
        sage: forget()
        sage: sin(n*pi)
        sin(pi*n)
    """
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            assume(*x)
        else:
            try:
                x.assume()
            except KeyError:
                raise TypeError, "assume not defined for objects of type '%s'"%type(x)

def forget(*args):
    """
    Forget the given assumption, or call with no arguments to forget
    all assumptions.

    Here an assumption is some sort of symbolic constraint.

    INPUT:


    -  ``*args`` - assumptions (default: forget all
       assumptions)


    EXAMPLES: We define and forget multiple assumptions::

        sage: var('x,y,z')
        (x, y, z)
        sage: assume(x>0, y>0, z == 1, y>0)
        sage: assumptions()
        [x > 0, y > 0, z == 1]
        sage: forget(x>0, z==1)
        sage: assumptions()
        [y > 0]
        sage: assume(y, 'even')
        sage: assumptions()
        [y > 0, y is even]
        sage: cos(y*pi)
        1
        sage: forget()
        sage: cos(y*pi)
        cos(pi*y)
        sage: assumptions()
        []
    """
    if len(args) == 0:
        _forget_all()
        return
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            assume(*x)
        else:
            try:
                x.forget()
            except KeyError:
                raise TypeError, "forget not defined for objects of type '%s'"%type(x)

def assumptions():
    """
    List all current symbolic assumptions.

    EXAMPLES::

        sage: var('x,y,z, w')
        (x, y, z, w)
        sage: forget()
        sage: assume(x^2+y^2 > 0)
        sage: assumptions()
        [y^2 + x^2 > 0]
        sage: forget(x^2+y^2 > 0)
        sage: assumptions()
        []
        sage: assume(x > y)
        sage: assume(z > w)
        sage: assumptions()
        [x > y, z > w]
        sage: forget()
        sage: assumptions()
        []
    """
    return list(_assumptions)

def _forget_all():
    """
    Forget all symbolic assumptions.

    This is called by ``forget()``.

    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: assume(x > 0, y < 0)
        sage: bool(x*y < 0)      # means definitely true
        True
        sage: bool(x*y > 0)      # might not be true
        False
        sage: forget()    # implicitly calls _forget_all
        sage: bool(x*y < 0)      # might not be true
        False
        sage: bool(x*y > 0)      # might not be true
        False
    """
    global _assumptions
    if len(_assumptions) == 0:
        return
    try:
        maxima._eval_line('forget(facts());')
    except TypeError:
        pass
    #maxima._eval_line('forget([%s]);'%(','.join([x._maxima_init_() for x in _assumptions])))
    for x in _assumptions:
        if isinstance(x, GenericDeclaration):
            # these don't show up in facts()
            x.forget()
    _assumptions = []

def solve(f, *args, **kwds):
    r"""
    Algebraically solve an equation of system of equations for given
    variables.

    INPUT:


    -  ``f`` - equation or system of equations (given by a
       list or tuple)

    -  ``*args`` - variables to solve for.

    -  ``solution_dict = True`` - return a list of
       dictionaries containing the solutions.


    EXAMPLES::

        sage: x, y = var('x, y')
        sage: solve([x+y==6, x-y==4], x, y)
        [[x == 5, y == 1]]
        sage: solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y)
        [[x == (-sqrt(3)*I - 1)/2, y == -sqrt(3 - sqrt(3)*I)/sqrt(2)],
         [x == (-sqrt(3)*I - 1)/2, y == sqrt(3 - sqrt(3)*I)/sqrt(2)],
         [x == (sqrt(3)*I - 1)/2, y == -sqrt(sqrt(3)*I + 3)/sqrt(2)],
         [x == (sqrt(3)*I - 1)/2, y == sqrt(sqrt(3)*I + 3)/sqrt(2)],
         [x == 0, y == -1],
         [x == 0, y == 1]]
        sage: solutions=solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y, solution_dict=True); solutions
        [{y: -sqrt(3 - sqrt(3)*I)/sqrt(2), x: (-sqrt(3)*I - 1)/2},
         {y: sqrt(3 - sqrt(3)*I)/sqrt(2), x: (-sqrt(3)*I - 1)/2},
         {y: -sqrt(sqrt(3)*I + 3)/sqrt(2), x: (sqrt(3)*I - 1)/2},
         {y: sqrt(sqrt(3)*I + 3)/sqrt(2), x: (sqrt(3)*I - 1)/2},
         {y: -1, x: 0},
         {y: 1, x: 0}]
        sage: for solution in solutions: print solution[x].n(digits=3), ",", solution[y].n(digits=3)
        -0.500 - 0.866*I , -1.27 + 0.341*I
        -0.500 - 0.866*I , 1.27 - 0.341*I
        -0.500 + 0.866*I , -1.27 - 0.341*I
        -0.500 + 0.866*I , 1.27 + 0.341*I
        0.000 , -1.00
        0.000 , 1.00
        sage: z = 5
        sage: solve([8*z + y == 3, -z +7*y == 0],y,z)
        Traceback (most recent call last):
        ...
        TypeError: 5 is not a valid variable.

    If ``True`` appears in the list of equations it is
    ignored, and if ``False`` appears in the list then no
    solutions are returned. E.g., note that the first
    ``3==3`` evaluates to ``True``, not to a
    symbolic equation.

    ::

        sage: solve([3==3, 1.00000000000000*x^3 == 0], x)
        [x == 0]
        sage: solve([1.00000000000000*x^3 == 0], x)
        [x == 0]

    Here, the first equation evaluates to ``False``, so
    there are no solutions::

        sage: solve([1==3, 1.00000000000000*x^3 == 0], x)
        []

    ::

        sage: var('s,i,b,m,g')
        (s, i, b, m, g)
        sage: sys = [ m*(1-s) - b*s*i, b*s*i-g*i ];
        sage: solve(sys,s,i)
        [[s == 1, i == 0], [s == g/b, i == (b - g)*m/(b*g)]]
        sage: solve(sys,(s,i))
        [[s == 1, i == 0], [s == g/b, i == (b - g)*m/(b*g)]]
        sage: solve(sys,[s,i])
        [[s == 1, i == 0], [s == g/b, i == (b - g)*m/(b*g)]]
    """
    try:
        return f.solve(*args,**kwds)
    except AttributeError:

        try:
            variables = tuple(args[0])
        except TypeError:
            variables = args

        from sage.calculus.calculus import SymbolicVariable, SymbolicFunction, SymbolicFunctionEvaluation
        for v in variables:
            if not isinstance(v,(SymbolicVariable,SymbolicFunction,SymbolicFunctionEvaluation)):
                raise TypeError, "%s is not a valid variable."%v

        try:
            f = [s for s in f if s is not True]
        except TypeError:
            raise ValueError, "Unable to solve %s for %s"%(f, args)

        if any(s is False for s in f):
            return []

        m = maxima(f)

        try:
            s = m.solve(variables)
        except:
            raise ValueError, "Unable to solve %s for %s"%(f, args)
        a = repr(s)
        sol_list = string_to_list_of_solutions(a)
        if 'solution_dict' in kwds and kwds['solution_dict']==True:
            sol_dict=[dict([[eq.left(),eq.right()] for eq in solution]) for solution in sol_list]
            return sol_dict
        else:
            return sol_list

import sage.categories.all
objs = sage.categories.all.Objects()

def string_to_list_of_solutions(s):
    r"""
    Used internally by the symbolic solve command to convert the output
    of Maxima's solve command to a list of solutions in Sage's symbolic
    package.

    EXAMPLES: We derive the (monic) quadratic formula::

        sage: var('x,a,b')
        (x, a, b)
        sage: solve(x^2 + a*x + b == 0, x)
        [x == (-sqrt(a^2 - 4*b) - a)/2, x == (sqrt(a^2 - 4*b) - a)/2]

    Behind the scenes when the above is evaluated the function
    ``string_to_list_of_solutions`` is called with
    input the string `s` below::

        sage: s = '[x=-(sqrt(a^2-4*b)+a)/2,x=(sqrt(a^2-4*b)-a)/2]'
        sage: sage.calculus.equations.string_to_list_of_solutions(s)
        [x == (-sqrt(a^2 - 4*b) - a)/2, x == (sqrt(a^2 - 4*b) - a)/2]
    """
    from sage.calculus.calculus import symbolic_expression_from_maxima_string
    v = symbolic_expression_from_maxima_string(s, equals_sub=True)
    return Sequence(v, universe=objs, cr_str=True)


############################################################
# Solving modulo N

def solve_mod(eqns, modulus, solution_dict = False):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    By default the solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``modulus`` - an integer

    -  ``solution_dict`` - (default: False) if True,  return a list of
       dictionaries containing the solutions.

    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We can solve with respect to a bigger modulus if it consists only of small prime factors::

        sage: solve_mod([5*x + y == 3, 2*x - 3*y == 9], 3*5*7*11*19*23*29, solution_dict = True)
        [{y: 8610183, x: 12915279}]

    We solve an simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]

    .. warning::

       The current implementation splits the modulus into prime powers,
       then naively enumerates all possible solutions and finally combines
       the solution using the Chinese Remainder Theorem.
       The interface is good, but the algorithm is horrible if the modulus
       has some larger prime factors! Sage {does} have the ability to do
       something much faster in certain cases at least by using Groebner
       basis, linear algebra techniques, etc. But for a lot of toy problems
       this function as is might be useful. At least it establishes an interface.
    """
    from sage.rings.all import Integer, Integers, PolynomialRing, factor, crt_basis
    from sage.misc.all import cartesian_product_iterator
    from sage.modules.all import vector
    from sage.matrix.all import matrix

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    modulus = Integer(modulus)
    if modulus < 1:
         raise ValueError, "the modulus must be a positive integer"
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(cmp = lambda x,y: cmp(repr(x), repr(y)))
    n = len(vars)

    factors = [p**i for p,i in factor(modulus)]
    crt_basis = vector(Integers(modulus), crt_basis(factors))
    solutions = [solve_mod_enumerate(eqns, p) for p in factors]

    ans = []
    for solution in cartesian_product_iterator(solutions):
        solution_mat = matrix(Integers(modulus), solution)
        ans.append(tuple(c.dot_product(crt_basis) for c in solution_mat.columns()))

    if solution_dict == True:
        sol_dict = [dict(zip(vars, solution)) for solution in ans]
        return sol_dict
    else:
        return ans

def solve_mod_enumerate(eqns, modulus):
    r"""
    Return all solutions to an equation or list of equations modulo the
    given integer modulus. Each equation must involve only polynomials
    in 1 or many variables.

    The solutions are returned as `n`-tuples, where `n`
    is the number of variables appearing anywhere in the given
    equations. The variables are in alphabetical order.

    INPUT:


    -  ``eqns`` - equation or list of equations

    -  ``modulus`` - an integer


    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(4, 2), (4, 6), (4, 9), (4, 13)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5::

        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    We solve an simple equation modulo 2::

        sage: x,y = var('x,y')
        sage: solve_mod([x == y], 2)
        [(0, 0), (1, 1)]


    .. warning::

       Currently this naively enumerates all possible solutions.  The
       interface is good, but the algorithm is horrible if the modulus
       is at all large! Sage {does} have the ability to do something
       much faster in certain cases at least by using the Chinese
       Remainder Theorem, Groebner basis, linear algebra techniques,
       etc. But for a lot of toy problems this function as is might be
       useful. At least it establishes an interface.
    """
    from sage.rings.all import Integer, Integers, PolynomialRing
    from calculus import is_SymbolicExpression
    from sage.misc.all import cartesian_product_iterator

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    modulus = Integer(modulus)
    if modulus < 1:
         raise ValueError, "the modulus must be a positive integer"
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort(cmp = lambda x,y: cmp(repr(x), repr(y)))
    n = len(vars)
    R = Integers(modulus)
    S = PolynomialRing(R, len(vars), vars)
    eqns_mod = [S(eq) if is_SymbolicExpression(eq) else
                S(eq.lhs() - eq.rhs()) for eq in eqns]
    ans = []
    for t in cartesian_product_iterator([R]*len(vars)):
        is_soln = True
        for e in eqns_mod:
            if e(t) != 0:
                is_soln = False
                break
        if is_soln:
            ans.append(t)

    return ans
