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
        return number_field_element.NumberFieldElement_absolute(self.__K, f)

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


class IdentityMap(SageObject):
    def __init__(self, K):
        self.__K = K
    def _repr_(self):
        return "Identity map on %s"%self.__K
    def __call__(self, x):
        return self.__K(x)
    def domain(self):
        return self.__K
    def codomain(self):
        return self.__K

class MapRelativeToAbsoluteNumberField(SageObject):
    def __init__(self, R, A):
        self.__R = R          # relative field
        self.__A = A          # absolute field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__R, self.__A)

    def __call__(self, x):
        f = self.__R(x).polynomial()
        return number_field_element.NumberFieldElement_absolute(self.__A, f)

    def domain(self):
        return self.__R

    def codomain(self):
        return self.__A

class MapAbsoluteToRelativeNumberField(SageObject):
    def __init__(self, A, R):
        self.__A = A          # absolute field
        self.__R = R          # relative field
        self.__poly_ring = self.__A.polynomial_ring()
        self.__zero = QQ(0)
        self.__n = A.degree()

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__A, self.__R)

    def __call__(self, x):
        f = self.__A(x).polynomial()
        return number_field_element.NumberFieldElement_relative(self.__R, f)

    def domain(self):
        return self.__A

    def codomain(self):
        return self.__R
