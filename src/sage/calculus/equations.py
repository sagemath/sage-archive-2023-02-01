r"""
Symbolic Equations and Inequalities.

\sage can solve symbolic equations and express inequalities.
For example, we derive the quadratic formula as follows:

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
    -- Bobby Moretti: initial version
    -- William Stein: second version
    -- William Stein (2007-07-16): added arithmetic with symbolic equations

EXAMPLES:
    sage: x,y,a = var('x,y,a')
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(1 - y^2), x == sqrt(1 - y^2)]

    sage: f = x^5 + a
    sage: solve(f==0,x)
    [x == e^(2*I*pi/5)*(-a)^(1/5), x == e^(4*I*pi/5)*(-a)^(1/5), x == e^(-(4*I*pi/5))*(-a)^(1/5), x == e^(-(2*I*pi/5))*(-a)^(1/5), x == (-a)^(1/5)]



You can also do arithmetic with inequalities, as illustrated below:
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
"""

_assumptions = []

from sage.structure.sage_object import SageObject
from sage.structure.sequence    import Sequence
from sage.misc.sage_eval        import sage_eval

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
    return cmp(str(x), str(y))

def paren(x):
    r = repr(x)
    if x._is_atomic() or not ('+' in r or '-' in r):
        return r
    else:
        return '(%s)'%r

def is_SymbolicEquation(x):
    r"""
    Return True if x is a symbolic equation.

    EXAMPLES:
    The following two examples are symbolic equations:
        sage: is_SymbolicEquation(sin(x) == x)
        True
        sage: is_SymbolicEquation(sin(x) < x)
        True

    This is not, since \code{2==3} evaluates to the boolean \code{False}:
        sage: is_SymbolicEquation(2 == 3)
        False

    However here since both 2 and 3 are coerced to be symbolic, we obtain
    a symbolic equation:
        sage: is_SymbolicEquation(SR(2) == SR(3))
        True
    """
    return isinstance(x, SymbolicEquation)

class SymbolicEquation(SageObject):
    def __init__(self, left, right, op):
        self._left = left
        self._right = right
        self._op = op

    def __call__(self, *args, **argv):
        return self._op(self._left(*args, **argv), self._right(*args,**argv))

    def __getitem__(self, i):
        return [self._left, self._op, self._right][i]

    def _scalar(self, scalar, op):
        """
        TESTS:
            sage: var('x y')
            (x, y)
            sage: f = x + 3 < y - 2
            sage: f*-1
            -x - 3 > 2 - y
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
        if self._op == operator.eq or op in [operator.add, operator.sub]:
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
        Multiply both sides of this inequality by $x$.

        EXAMPLES:
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

        Multiplying by complex numbers works only if it's an equality:
            sage: f = sqrt(2) + x == y^3
            sage: f.multiply_both_sides(I)
            I*(x + sqrt(2)) == I*y^3
            sage: f.multiply_both_sides(-1)
            -x - sqrt(2) == -y^3

        Some further examples:
            sage: (x^3 + 1 > 2*sqrt(3)) * (-1)
            -x^3 - 1 < -2*sqrt(3)
            sage: (x^3 + 1 >= 2*sqrt(3)) * (-1)
            -x^3 - 1 <= -2*sqrt(3)
            sage: (x^3 + 1 <= 2*sqrt(3)) * (-1)
            -x^3 - 1 >= -2*sqrt(3)
        """
        return self._scalar(x, operator.mul)

    def divide_both_sides(self, x):
        """
        Divide both sides of the inequality by $x$.

        EXAMPLES:
            sage: (x^3 + 1 > x^2 - 1) / (-1)
            -x^3 - 1 < 1 - x^2

        The quantity $x^2 - 1$ could be either negative or positive depending on $x$, so
        dividing by it is not defined.
            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            Traceback (most recent call last):
            ...
            ValueError: unable to multiply or divide both sides of an inequality by a number whose sign can't be determined.

        If we specify that $x^2 - 1> 0$, then dividing is defined.
            sage: assume(x^2 - 1 > 0)
            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            (x^3 + 1)/(x^2 - 1) > 1
            sage: forget()

        We can also specify that $x^2 - 1 < 0$.  Note that now the inequality direction changes.
            sage: assume(x^2 - 1 < 0)
            sage: (x^3 + 1 > x^2 - 1) / (x^2 - 1)
            (x^3 + 1)/(x^2 - 1) < 1
            sage: forget()
        """
        return self._scalar(x, operator.div)

    def add_to_both_sides(self, x):
        """
        Add $x$ to both sides of this symbolic equation.

        EXAMPLES:
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
        Subtract $x$ from both sides of this symbolic equation.

        EXAMPLES:
            sage: eqn = x*sin(x)*sqrt(3) + sqrt(2) > cos(sin(x))
            sage: eqn.subtract_from_both_sides(sqrt(2))
            sqrt(3)*x*sin(x) > cos(sin(x)) - sqrt(2)
            sage: eqn.subtract_from_both_sides(cos(sin(x)))
            -cos(sin(x)) + sqrt(3)*x*sin(x) + sqrt(2) > 0
        """
        return self._scalar(x, operator.sub)

    def _arith(self, right, op):
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

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
            sage: m = x == 5*x + 1
            sage: n = sin(x) == sin(x+2*pi)
            sage: m * n
            x*sin(x) == (5*x + 1)*sin(x)
        """
        return self._arith(right, operator.mul)

    def __rmul__(self, left):
        """
        Multiply two symbolic equations.
            sage: m = 2*x == 3*x^2 - 5
            sage: int(-1) * m
            -2*x == 5 - 3*x^2
        """
        return self._arith(left, operator.mul)

    def __div__(self, right):
        """
        Divide two symbolic equations.

        EXAMPLES:
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
        if session is None:
            return SageObject._maxima_(self, sage.calculus.calculus.maxima)
        else:
            return SageObject._maxima_(self, session)

    def substitute(self, *args, **kwds):
        return self.__call__(*args, **kwds)

    subs = substitute

    def __cmp__(self, right):
        """
        EXAMPLES:
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
        return "%r%s%r" %(self._left, symbols[self._op], self._right)

    def variables(self):
        """
        EXAMPLES:
            sage: var('x,y,z,w')
            (x, y, z, w)
            sage: f =  (x+y+w) == (x^2 - y^2 - z^3);   f
            y + x + w == -z^3 - y^2 + x^2
            sage: f.variables()
            [w, x, y, z]
        """
        try:
            return self.__variables
        except AttributeError:
            v = list(set(list(self._left.variables()) + list(self._right.variables())))
            v.sort(var_cmp)
            self.__variables = v
            return v

    def operator(self):
        return self._op

    def left(self):
        return self._left
    lhs = left
    left_hand_side = left

    def right(self):
        return self._right
    rhs = right
    right_hand_side = right

    def __str__(self):
        """
        EXAMPLES:
            sage: f =  (x^2 - x == 0)
            sage: f
            x^2 - x == 0
            sage: print f
                                               2
                                              x  - x == 0
        """
        s = self._maxima_().display2d(onscreen=False)
        s = s.replace('%pi','pi').replace('%i',' I').replace('%e', ' e').replace(' = ',' == ')
        return s

    def _latex_(self):
        return "%s %s %s" %(self._left._latex_(), latex_symbols[self._op],
                            self._right._latex_())

    # this is an excellent idea by Robert Bradshaw
    #def __nonzero__(self):
    #    result = self._left.__cmp__(self._right)
    #    return result in comparisons[self._op]
    def __nonzero__(self):
        """
        Return True if this (in)equality is definitely true.  Return False
        if it is false or the algorithm for testing (in)equality is
        inconclusive.

        EXAMPLES:
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
        l = self._left._maxima_init_()
        r = self._right._maxima_init_()
        if assume and self._op == operator.eq:
            return 'equal(%s, %s)'%(l, r)
        return '(%s)%s(%s)' % (l, maxima_symbols[self._op], r)

    def assume(self):
        if not self in _assumptions:
            m = self._maxima_init_(assume=True)
            maxima.assume(m)
            _assumptions.append(self)

    def find_root(self, *args, **kwds):
        r"""
        If this is a symbolic equality with an equals sign \code{==}
        find numerically a single root of this equation in a given
        interval.  Otherwise raise a \code{ValueError}.  See the
        documentation for the global \code{find_root} method for more
        about the options to this function.

        Note that this symbolic expression must involve at most one
        variable.

        EXAMPLES:
            sage: (x == sin(x)).find_root(-2,2)
            0.0
            sage: (x^5 + 3*x + 2 == 0).find_root(-2,2)
            -0.63283452024215225
            sage: (cos(x) == sin(x)).find_root(10,20)
            19.634954084936208

        We illustrate some valid error conditions:
            sage: (cos(x) != sin(x)).find_root(10,20)
            Traceback (most recent call last):
            ...
            ValueError: Symbolic equation must be an equality.
            sage: (SR(3)==SR(2)).find_root(-1,1)
            Traceback (most recent call last):
            ...
            RuntimeError: no zero in the interval, since constant expression is not 0.

        There must be at most one variable:
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

        WARNING: In many cases, only one solution is computed.

        INPUT:
            x -- a SymbolicVariable object (if not given, the first in
                 the expression is used)

            multiplicities -- (default: False) if True, also returns
                          the multiplicities of each solution, in order.

            solution_dict -- (default: False) if True, return the
                        solution as a dictionary rather than an equation.

            explicit_solution -- (default: False); if True, require
                that all solutions returned be explicit (rather than
                implicit) OUTPUT: A list of SymbolicEquations with the
                variable to solve for on the left hand side.

        EXAMPLES:
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

        We illustrate finding multiplicities of solutions:
            sage: f = (x-1)^5*(x^2+1)
            sage: solve(f == 0, x)
            [x == -1*I, x == I, x == 1]
            sage: solve(f == 0, x, multiplicities=True)
            ([x == -1*I, x == I, x == 1], [1, 1, 5])
        """
        if not self._op is operator.eq:
            raise ValueError, "solving only implemented for equalities"
        if x is None:
            v = self.variables()
            if len(v) == 0:
                if multiplicities:
                    return [], []
                else:
                    return []
            x = v[0]

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
                return X, sage_eval(P.get('multiplicities'))
        else:
            return X


    def expand(self, side=None):
        r"""
        Expands one or both sides of the equation.

        If side is not specified, then both sides of the equation
        are expanded by calling \code{expand()} on the corresponding
        \class{SymbolicExpression}.

        If side is `left' (or `right'), then only the left (or right)
        side of the equation is expanded.

        EXAMPLES:
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


def assume(*args):
    """
    Make the given assumptions.

    INPUT:
        *args -- assumptions

    EXAMPLES:
        sage: assume(x > 0)
        sage: bool(sqrt(x^2) == x)
        True
        sage: forget()
        sage: bool(sqrt(x^2) == x)
        False

    An integer constraint (todo: this needs to be made possible with
    just the assume command!):
        sage: from sage.calculus.calculus import maxima as calcmaxima
        sage: calcmaxima.eval('declare(n,integer)')
        'done'
        sage: var('n, P, r, r2')
        (n, P, r, r2)
        sage: c = P*e^(r*n)
        sage: d = P*(1+r2)^n
        sage: solve(c==d,r2)
        [r2 == e^r - 1]
    """
    for x in args:
        if isinstance(x, (tuple, list)):
            for y in x:
                assume(y)
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
        *args -- assumptions (default: forget all assumptions)

    EXAMPLES:
    We define and forget multiple assumptions:
        sage: var('x,y,z')
        (x, y, z)
        sage: assume(x>0, y>0, z == 1, y>0)
        sage: assumptions()
        [x > 0, y > 0, z == 1]
        sage: forget(x>0, z==1)
        sage: assumptions()
        [y > 0]
    """
    if len(args) == 0:
        forget_all()
        return
    for x in args:
        if isinstance(x, (tuple, list)):
            for y in x:
                assume(y)
        else:
            try:
                x.forget()
            except KeyError:
                raise TypeError, "forget not defined for objects of type '%s'"%type(x)

def assumptions():
    """
    List all current symbolic assumptions.

    EXAMPLES:
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

def forget_all():
    global _assumptions
    if len(_assumptions) == 0:
        return
    try:
        maxima._eval_line('forget(facts());')
    except TypeError:
        pass
    #maxima._eval_line('forget([%s]);'%(','.join([x._maxima_init_() for x in _assumptions])))
    _assumptions = []

def solve(f, *args, **kwds):
    r"""
    Algebraically solve an equation of system of equations for given variables.

    INPUT:
        f -- equation or system of equations (given by a list or tuple)
        *args -- variables to solve for.
	solution_dict = True -- return a list of dictionaries containing the solutions.

    EXAMPLES:
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

    If \code{True} appears in the list of equations it is ignored, and if
    \code{False} appears in the list then no solutions are returned.  E.g.,
    note that the first \code{3==3} evaluates to \code{True}, not to a symbolic
    equation.

        sage: solve([3==3, 1.00000000000000*x^3 == 0], x)
        [x == 0]
        sage: solve([1.00000000000000*x^3 == 0], x)
        [x == 0]

    Here, the first equation evaluates to \code{False}, so there are no solutions:
        sage: solve([1==3, 1.00000000000000*x^3 == 0], x)
        []

    """
    if isinstance(f, (list, tuple)):
        f = [s for s in f if s is not True]
        for s in f:
            if s is False:
                return []

        m = maxima(list(f))
        try:
            s = m.solve(args)
        except:
            raise ValueError, "Unable to solve %s for %s"%(f, args)
        a = repr(s)
	sol_list = string_to_list_of_solutions(a)
	if 'solution_dict' in kwds and kwds['solution_dict']==True:
            sol_dict=[dict([[eq.left(),eq.right()] for eq in solution]) for solution in sol_list]
            return sol_dict
        else:
            return sol_list
    else:
        return f.solve(*args, **kwds)

import sage.categories.all
objs = sage.categories.all.Objects()

def string_to_list_of_solutions(s):
    from sage.calculus.calculus import symbolic_expression_from_maxima_string
    v = symbolic_expression_from_maxima_string(s, equals_sub=True)
    return Sequence(v, universe=objs, cr_str=True)


############################################################
# Solving modulo N

def solve_mod(eqns, modulus):
    r"""
    Return all solutions to an equation or list of equations modulo
    the given integer modulus.  Each equation must involve only
    polynomials in 1 or many variables.

    The solutions are returned as $n$-tuples, where $n$ is the
    number of variables appearing anywhere in the given equations.
    The variables are in alphabetical order.


    INPUT:
        eqns -- equation or list of equations
        modulus -- an integer

    EXAMPLES:
        sage: var('x,y')
        (x, y)
        sage: solve_mod([x^2 + 2 == x, x^2 + y == y^2], 14)
        [(2, 4), (6, 4), (9, 4), (13, 4)]
        sage: solve_mod([x^2 == 1, 4*x  == 11], 15)
        [(14,)]

    Fermat's equation modulo 3 with exponent 5:
        sage: var('x,y,z')
        (x, y, z)
        sage: solve_mod([x^5 + y^5 == z^5], 3)
        [(0, 0, 0), (0, 1, 1), (0, 2, 2), (1, 0, 1), (1, 1, 2), (1, 2, 0), (2, 0, 2), (2, 1, 0), (2, 2, 1)]

    WARNING:
        Currently this naively enumerates all possible solutions.
        The interface is good, but the algorithm is horrible if the
        modulus is at all large!   \sage \strong{does} have the ability to do
        something much faster in certain cases at least by using
        the Chinese Remainder Theorem, Gr\"obner basis, linear algebra
        techniques, etc.  But for a lot of toy problems this function
        as is might be useful.  At least it establishes an interface.
    """
    from sage.rings.all import Integer, Integers, MPolynomialRing
    from calculus import is_SymbolicExpression
    from sage.misc.all import cartesian_product_iterator

    if not isinstance(eqns, (list, tuple)):
        eqns = [eqns]
    modulus = Integer(modulus)
    if modulus < 1:
         raise ValueError, "the modulus must be a positive integer"
    vars = list(set(sum([list(e.variables()) for e in eqns], [])))
    vars.sort()
    n = len(vars)
    R = Integers(modulus)
    S = MPolynomialRing(R, len(vars), vars)
    eqns_mod = [S(eq) if is_SymbolicExpression(eq) else \
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
