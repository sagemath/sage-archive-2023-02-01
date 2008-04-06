"""
TESTS:
   sage: L.<cuberoot2, zeta3> = CyclotomicField(3).extension(x^3 - 2); K = L.absolute_field('a'); from_K,to_K = K.structure()

"""


from sage.structure.sage_object import SageObject

import sage.rings.rational_field as rational_field

from sage.libs.pari.all import pari

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
        return self.__K._element_class(self.__K, f)

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
        y = self.__K(x)
        v = y._coefficients()
        w = v + [self.__zero]*(self.__n - len(v))
        return self.__V(w)

    def domain(self):
        return self.__K

    def codomain(self):
        return self.__V

class MapRelativeVectorSpaceToRelativeNumberField(SageObject):
    def __init__(self, V, K):
        self.__V = V
        self.__K = K
        self.__R = K.polynomial_ring()
        self.__rnf = K.pari_rnf()
        self.__B = K.base_field().absolute_field('a')

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__V, self.__K)

    def __call__(self, v):
        # Given a relative vector space element, we have to
        # compute the corresponding number field element, in terms
        # of an absolute generator.
        w = self.__V(v).list()

        # First, construct from w a PARI polynomial in x with coefficients
        # that are polynomials in y:
        B = self.__B
        _, to_B = B.structure()
        # Apply to_B, so now each coefficient is in an absolute field,
        # and is expressed in terms of a polynomial in y, then make
        # the PARI poly in x.
        w = [pari(to_B(a).polynomial('y')) for a in w]
        h = pari(w).Polrev()

        # Next we write the poly in x over a poly in y in terms
        # of an absolute polynomial for the rnf.
        g = self.__R(self.__rnf.rnfeltreltoabs(h))
        return self.__K._element_class(self.__K, g)

    def domain(self):
        return self.__V

    def codomain(self):
        return self.__K

class MapRelativeNumberFieldToRelativeVectorSpace(SageObject):
    def __init__(self, K, V):
        self.__V = V
        self.__K = K
        self.__rnf = K.pari_rnf()
        self.__zero = QQ(0)
        self.__n = K.degree()
        self.__x = pari('x')
        self.__y = pari('y')
        self.__B = K.absolute_base_field()

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__K, self.__V)

    def __call__(self, alpha):
        # An element of a relative number field is represented
        # internally by an absolute polynomial over QQ.
        alpha = self.__K(alpha)
        f = alpha.polynomial('x')
        # f is the absolute polynomial that defines this number field element
        if self.__K.degree() == 1:
            # Pari doesn't return a valid relative polynomial from an
            # absolute one in the case of relative degree one. (Bug
            # submitted to Pari, fixed in svn unstable as of 1/22/08
            # according to Karim Belabas. Trac #1891 is a reminder to
            # remove this workaround once we update Pari in Sage.)
            # However, in this case, the polynomial we want for g is
            # just the polynomial we have for f, but as a polynomial
            # in the generator of the base field of this relative
            # extension. That generator is just -self.__rnf[0][0], so
            # we can plug it in and it works.
            g = f(-1*self.__rnf[0][0])
        else:
            g = self.__rnf.rnfeltabstorel(pari(f))
        # Now g is a relative polynomial that defines this element.
        # This g is a polynomial in a pari variable x with
        # coefficients polynomials in a variable y.
        # These coefficients define the coordinates of the
        # vector we are constructing.

        # The list v below has the coefficients that are the
        # components of the vector we are constructing, but each is
        # converted into polynomials in a variable x, which we will
        # use to define elements of the base field.
        (x, y) = (self.__x, self.__y)
        v = [g.polcoeff(i).subst(x,y) for i in range(self.__n)]
        B,from_B, _ = self.__B
        w = [from_B(B(z)) for z in v]

        # Now w gives the coefficients.
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

class NameChangeMap(SageObject):
    def __init__(self, K, L):
        self.__K = K
        self.__L = L

    def _repr_(self):
        return "Number field isomorphism from %s to %s given by variable name change"%(self.__K, self.__L)

    def __call__(self, x):
        y = self.__K(x)
        z = y.__copy__()
        z._set_parent(self.__L)
        return z

    def domain(self):
        return self.__K

    def codomain(self):
        return self.__L

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
        return self.__A._element_class(self.__A, f)

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
        return self.__R._element_class(self.__R, f)

    def domain(self):
        return self.__A

    def codomain(self):
        return self.__R


class MapVectorSpaceToRelativeNumberField(SageObject):
    def __init__(self, V, L, from_V, from_K):
        self.__V = V
        self.__L = L
        self.__from_V = from_V
        self.__from_K = from_K

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__V, self.__L)

    def __call__(self, x):
        return self.__from_K(self.__from_V(x))

    def domain(self):
        return self.__V

    def codomain(self):
        return self.__L

class MapRelativeNumberFieldToVectorSpace(SageObject):
    def __init__(self, L, V, to_K, to_V):
        self.__L = L
        self.__V = V
        self.__to_K = to_K
        self.__to_V = to_V

    def _repr_(self):
        return "Isomorphism from %s to %s"%(self.__L, self.__V)

    def __call__(self, x):
        return self.__to_V(self.__to_K(x))

    def domain(self):
        return self.__V

    def codomain(self):
        return self.__L

