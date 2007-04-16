"""
Preliminary support for equations and solutions in \sage.

AUTHOR:
    -- Bobby Moretti initial version
    -- William Stein: rewrite; add inequalities

EXAMPLES:
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(1 - y^2), x == sqrt(1 - y^2)]

    sage: f = x^5 + a
    sage: solve(f==0,x)
    [x == -e^(2*I*pi/5)*a^(1/5), x == -e^(4*I*pi/5)*a^(1/5), x == -e^(-(4*I*pi/5))*a^(1/5), x == -e^(-(2*I*pi/5))*a^(1/5), x == (-a^(1/5))]
"""

_assumptions = []

from sage.structure.all import SageObject
from sage.interfaces.maxima import maxima


from sage.misc.sage_eval import sage_eval

import operator

symbols = {operator.lt:' < ', operator.le:' <= ', operator.eq:' == ',
           operator.ne:' != ',
            operator.ge:' >= ', operator.gt:' > '}

maxima_symbols = dict(symbols)
maxima_symbols[operator.eq] = '='


latex_symbols = {operator.lt:' < ', operator.le:' \\leq ', operator.eq:' = ',
                 operator.ne:' \\neq ',
            operator.ge:' \\geq ', operator.gt:' > '}

comparisons = {operator.lt:set([-1]), operator.le:set([-1,0]), operator.eq:set([0]),
               operator.ne:set([-1,1]),  operator.ge:set([0,1]), operator.gt:set([1])}


def paren(x):
    if x._is_atomic():
        return repr(x)
    else:
        return '(%s)'%repr(x)

class SymbolicEquation(SageObject):
    def __init__(self, left, right, op):
        self._left = left
        self._right = right
        self._op = op

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
        return "%s%s%s" %(paren(self._left), symbols[self._op], paren(self._right))

    def variables(self):
        """
        EXAMPLES:
            sage: f =  (x+y+w) == (x^2 - y^2 - z^3);   f
            (y + x + w) == (-z^3 - y^2 + x^2)
            sage: f.variables()
            [w, x, y, z]
        """
        try:
            return self.__variables
        except AttributeError:
            v = list(set(list(self._left.variables()) + list(self._right.variables())))
            v.sort()
            self.__variables = v
            return v

    def operator(self):
        return self._op

    def left(self):
        return self._left

    def right(self):
        return self._right

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
        s = s.replace('%pi',' Pi').replace('%i',' I').replace('%e', ' E').replace(' = ',' == ')
        return s

    def _latex_(self):
        return "%s %s %s" %(self._left._latex_(), latex_symbols[self._op],
                            self._right._latex_())

    # this is an excellent idea by Robert Bradshaw
    #def __nonzero__(self):
    #    result = self._left.__cmp__(self._right)
    #    return result in comparisons[self._op]
    def __nonzero__(self):
        m = self._maxima_()
        try:
            s = m.parent()._eval_line('is (%s)'%m.name())
        except TypeError, msg:
            #raise ValueError, "unable to evaluate the predicate '%s'"%repr(self)
            return cmp(self._left._maxima_() , self._right._maxima_()) == 0
        return s == 'true'

    def _maxima_init_(self, maxima=maxima):
        l = self._left._maxima_init_()
        r = self._right._maxima_init_()
        return '(%s)%s(%s)' % (l, maxima_symbols[self._op], r)

    def assume(self):
        if not self in _assumptions:
            m = self._maxima_()
            m.parent().assume(m)
            _assumptions.append(self)

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

    def solve(self, x=None, multiplicities=False):
        """
        Symbolically solve for the given variable.

        WARNING: In many cases, only one solution is computed.

        INPUT:
            x -- a SymbolicVariable object (if not given, the first in the expression is used)
            multiplicities -- (default: False) if True, also returns the multiplicities
                          of each solution, in order.

        OUTPUT:
            A list of SymbolicEquations with the variable to solve for on the
            left hand side.

        EXAMPLES:
            sage: S = solve(x^3 - 1 == 0, x)
            sage: S
            [x == ((sqrt(3)*I - 1)/2), x == ((-sqrt(3)*I - 1)/2), x == 1]
            sage: S[0]
            x == ((sqrt(3)*I - 1)/2)
            sage: S[0].right()
            (sqrt(3)*I - 1)/2

        We illustrate finding multiplicities of solutions:
            sage: f = (x-1)^5*(x^2+1)
            sage: solve(f == 0, x)
            [x == -1*I, x == I, x == 1]
            sage: solve(f == 0, x, multiplicities=True)
            ([x == -1*I, x == I, x == 1], [1, 1, 5])
            sage: solve(g == 0, x)
            []
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
        s = m.solve(x).str()

        X = string_to_list_of_solutions(s)
        if multiplicities:
            if len(X) == 0:
                return X, []
            else:
                return X, sage_eval(P.get('multiplicities'))
        else:
            return X



def assume(*args):
    """
    Make the given assumptions.

    INPUT:
        *args -- assumptions

    EXAMPLES:
        sage:
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
        sage: assume(x > 0)
        sage: sqrt(x^2)
        x
        sage: forget()
        Forgetting all assumptions.
        sage: sqrt(x^2)
        abs(x)

    We define and forget multiple assumptions:
        sage: assume(x>0, y>0, z == 1, y>0)
        sage: assumptions()
        [x > 0, y > 0, z == 1]
        sage: forget(x>0, z==1)
        sage: assumptions()
        [y > 0]
    """
    if len(args) == 0:
        print "Forgetting all assumptions."
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
    EXAMPLES:
        sage: assume(x^2+y^2 > 0)
        sage: assumptions()
        [(y^2 + x^2) > 0]
        sage: forget(x^2+y^2 > 0)
        sage: assumptions()
        []
        sage: assume(x > y)
        sage: assume(z > w)
        sage: assumptions()
        [x > y, z > w]
        sage: forget()
        Forgetting all assumptions.
        sage: assumptions()
        []
    """
    return list(_assumptions)

def forget_all():
    global _assumptions
    maxima._eval_line('forget([%s]);'%(','.join([x._maxima_init_() for x in _assumptions])))
    _assumptions = []

def solve(f, *args, **kwds):
    """
    Algebraically solve an equation of system of equations for given variables.

    INPUT:
        f -- equation or system of equations (given by a list or tuple)
        *args -- variables to solve for.

    EXAMPLES:
        sage: solve([x+y==6, x-y==4], x, y)
        [[x == 5, y == 1]]
        sage: solve([x^2+y^2 == 1, y^2 == x^3 + x + 1], x, y)
        [[x == ((-sqrt(3)*I - 1)/2), y == ((-sqrt(3 - sqrt(3)*I))/sqrt(2))],
         [x == ((-sqrt(3)*I - 1)/2), y == (sqrt(3 - sqrt(3)*I)/sqrt(2))],
         [x == ((sqrt(3)*I - 1)/2), y == ((-sqrt(sqrt(3)*I + 3))/sqrt(2))],
         [x == ((sqrt(3)*I - 1)/2), y == (sqrt(sqrt(3)*I + 3)/sqrt(2))],
         [x == 0, y == -1],
         [x == 0, y == 1]]
    """
    if isinstance(f, (list, tuple)):
        m = maxima(list(f))
        try:
            s = m.solve(args)
        except:
            raise ValueError, "Unable to solve %s for %s"%(f, args)
        a = repr(s)
        return string_to_list_of_solutions(a)
    else:
        return f.solve(*args, **kwds)

def string_to_list_of_solutions(s):
    from sage.calculus.calculus import symbolic_expression_from_maxima_string
    return symbolic_expression_from_maxima_string(s, equals_sub=True)
