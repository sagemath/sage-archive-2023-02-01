from sage.structure.sage_object import SageObject

import sage.rings.rational_field as rational_field

import number_field_element

QQ = rational_field.RationalField()

class MapVectorSpaceToNumberField(SageObject):
    def __init__(self, V, K):
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__V, self.__K)

    def __call__(self, v):
        f = self.__R(self.__V(v).list())
        return number_field_element.NumberFieldElement(self.__K, f)

    def domain(self):
        return self.__V

    def codomain(self):
        return self.__K

class MapNumberFieldToVectorSpace(SageObject):
    def __init__(self, K, V):
        self.__V = V
        self.__K = K
        self.__zero = QQ(0)
        self.__n = K.degree()

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__K, self.__V)

    def __call__(self, x):
        x = self.__K(x)
        v = x._coefficients()
        w = v + [self.__zero]*(self.__n - len(v))
        return self.__V(w)

    def domain(self):
        return self.__K

    def codomain(self):
        return self.__V
