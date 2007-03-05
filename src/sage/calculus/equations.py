"""nodoctest"""

from sage.structure.all import SageObject
from sage.interfaces.maxima import maxima
 #                                  symbolic_expression_from_maxima_element

class SymbolicEquation(SageObject):
    def __init__(self, left, right):
        self._left = left
        self._right = right

    def _repr_(self):
        return "%s == %s" %(self._left, self._right)

    def __nonzero__(self):
        result = self._left.__cmp__(self._right)

        if result == 0:
            return True
        else:
            return False

    def _maxima_(self):
        l = self._left._maxima_()._name
        r = self._right._maxima_()._name
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
