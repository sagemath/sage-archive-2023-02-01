"""nodoctest"""

from sage.structure.all import SageObject
from sage.interfaces.maxima import maxima
#from sage.calculus.calculus import symbolic_expression_from_maxima_string,
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
        """
        return self._maxima_().solve(x._maxima_())
