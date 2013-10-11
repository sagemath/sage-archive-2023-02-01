"""
Base class for Linear Expressions

A linear expression is just a linear polynomial in some (fixed)
variables. This class only implements linear expressions for others to
use.


EXAMPLES::

    sage: from sage.geometry.linear_expression import LinearExpressionModule
    sage: L.<x,y,z> = LinearExpressionModule(QQ)
    

"""

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method



class LinearExpression(ModuleElement):
    
    def __init__(self, parent, coefficients, constant):
        """
        A linear expression is a linear polynomial.
        """
        super(LinearExpression, self).__init__(parent)
        self._coeffs = coefficients
        self._const = constant
        


class LinearExpressionModule(Parent, UniqueRepresentation):

    Element = LinearExpression

    def __init__(self, base_ring, names=tuple()):
        """
        The parent of linear expressions.

        This is the module of linear polynomials.
        """
        from sage.categories.modules import Modules
        super(LinearExpressionModule, self).__init__(self, category=Modules(base_ring))
        self._names = names
        
    @cached_method
    def ngens(self):
        return len(self._names)

    @cached_method
    def gens(self):
        from sage.matrix.constructor import identity_matrix
        identity = identity_matrix(self.base_ring(), self.ngens())
        return tuple(self(e, 0) for e in identity.rows())
                     
            
