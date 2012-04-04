import operator
from sage.symbolic.ring import is_SymbolicVariable, var

arithmetic_operators = {operator.add: '+',
                        operator.sub: '-',
                        operator.mul: '*',
                        operator.div: '/',
                        operator.pow: '^'}

relation_operators = {operator.eq:'==',
                      operator.lt:'<',
                      operator.gt:'>',
                      operator.ne:'!=',
                      operator.le:'<=',
                      operator.ge:'>='}

class FDerivativeOperator(object):
    def __init__(self, function, parameter_set):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: loads(dumps(op))
            D[0, 1](foo)
        """
        self._f = function
        self._parameter_set = map(int, parameter_set)

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: x,y = var('x,y')
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op(x,y)
            D[0, 1](foo)(x, y)
            sage: op(x,x^2)
            D[0, 1](foo)(x, x^2)
        """
        if (not all(is_SymbolicVariable(x) for x in args) or
                len(args) != len(set(args))):
            temp_args=[var("t%s"%i) for i in range(len(args))]
            vars=[temp_args[i] for i in self._parameter_set]
            return self._f(*temp_args).diff(*vars).function(*temp_args)(*args)
        vars = [args[i] for i in self._parameter_set]
        return self._f(*args).diff(*vars)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1]); op
            D[0, 1](foo)
        """
        return "D[%s](%s)"%(", ".join(map(repr, self._parameter_set)), self._f)

    def function(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.function()
            foo
        """
        return self._f

    def change_function(self, new):
        """
        Returns a new FDerivativeOperator with the same parameter set
        for a new function.

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: b = function('bar')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.change_function(bar)
            D[0, 1](bar)
        """
        return FDerivativeOperator(new, self._parameter_set)

    def parameter_set(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.parameter_set()
            [0, 1]
        """
        return self._parameter_set
