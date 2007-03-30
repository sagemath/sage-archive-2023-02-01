"""
Preliminary support for equations and solutions in \sage.

AUTHOR:
    -- Bobby Moretti initial version

EXAMPLES:
    sage: f = x^2 + y^2 == 1
    sage: f.solve(x)
    [x == -sqrt(1 - y^2), x == sqrt(1 - y^2)]

"""

from sage.structure.all import SageObject
from sage.interfaces.maxima import maxima
 #                                  symbolic_expression_from_maxima_element

class SymbolicEquation(SageObject):
    def __init__(self, left, right):
        self._left = left
        self._right = right

    def _repr_(self):
        return "%s == %s" %(self._left, self._right)

    def _latex_(self):
        return "%s = %s" %(self._left._latex_(), self._right._latex_())

    # this is an excellent idea by Robert Bradshaw
    def __nonzero__(self):
        result = self._left.__cmp__(self._right)

        if result == 0:
            return True
        else:
            return False

    def _maxima_(self, maxima=maxima):
        l = str(self._left._maxima_())
        r = str(self._right._maxima_())
        return maxima('%s = %s' % (l, r))

    def solve(self, x):
        """
        Uses Maxima to symbolically solve for the given variable.

        INPUT:
            x -- a SymbolicVariable object

        OUTPUT:
            A list of SymbolicEquations with the variable to solve for on the
            left hand side.
        """

        # get the string result from solving wiht maxima
        s =  str(self._maxima_().solve(x._maxima_()))

        # get rid of the []'s, split the string around each comma
        sols = (s[1:-1]).split(',')

        # maxima returns a list of comma-separated solutions, we return
        sols = [sol.split('=') for sol in sols]

        # now our solutions are a list of list of strings, each string being a
        # symbolic expression
        from sage.calculus.calculus import symbolic_expression_from_maxima_string

        # apply the symexpr_from_maxima method to every element of every
        # solution
        symexprs=[map(symbolic_expression_from_maxima_string, sol) for sol in sols]

        # and finally construct the new SymbolicEquation
        return [SymbolicEquation(expr[0], expr[1]) for expr in symexprs]
