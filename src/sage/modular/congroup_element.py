"""
Element of a congruence subgroup.
"""

from sage.ext.element import MultiplicativeGroupElement
from sage.rings.all import IntegerRing
from sage.matrix.matrix_space import MatrixSpace

M2Z = MatrixSpace(IntegerRing(), 2)

class CongruenceSubgroupElement(MultiplicativeGroupElement):
    def __init__(self, parent, x, check=True):
        MultiplicativeGroupElement.__init__(self, parent)
        if check:
            x = M2Z(x)
            if x.determinant() != 1:
                raise TypeError, "matrix must have determinant 1"
        x.set_immutable()
        self.__x = x

    def __repr__(self):
        return "%s"%self.__x

    def _mul_(self, right):
        return self.parent()(self.__x * right.__x, check=False)

    def __invert__(self):
        I = M2Z([self.__x[1,1], -self.__x[0,1], -self.__x[1,0], self.__x[0,0]])
        return self.parent()(I, check=False)

    def matrix(self):
        return self.__x

    def determinant(self):
        return IntegerRing()(1)

    det = determinant

    def a(self):
        return self.__x[0,0]

    def b(self):
        return self.__x[0,1]

    def c(self):
        return self.__x[1,0]

    def d(self):
        return self.__x[1,1]
