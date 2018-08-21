from sage.rings.morphism import RingMap, RingHomomorphism
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
from sage.categories.map import Map
from sage.categories.all import Rings
from sage.categories.homset import End 


class RingDerivation(RingMap):
    """
    Derivation of rings
    """

    def __init__(self, parent, theta):
        self._ring = parent.domain()
        self._parent = parent
        if self._ring is not parent.codomain():
            raise NotImplementedError
        
        
        RingMap.__init__(self, parent)
        if (not (isinstance(theta, Map) and theta.category_for().is_subcategory(Rings()))
            or theta.domain() is not self._ring or theta.codomain() is not self._ring):
            raise TypeError("theta must be an homomorphism of %s" % self._ring)
        self._theta = theta



class RingDerivation_im_gens(RingDerivation):
    def __init__(self, parent, theta, images):
        RingDerivation.__init__(self, parent, theta)
        self._images = [ self._ring(x) for x in images]
        if len(self._images) != len(self._ring.gens()):
            raise ValueError("Number of images is incorrect")


class RingDerivation_polynomial(RingDerivation_im_gens):
    def __init__(self, parent, theta, images):
        if not theta.is_identity():
            raise NotImplementedError
        RingDerivation_im_gens.__init__(self, parent, theta, images)
        if not isinstance(self._ring, (PolynomialRing_general, MPolynomialRing_generic)):
            raise TypeError("The ring is not a polynomial ring")
        
    def _coerce_map_from_(self, R):
        if self.base_ring().has_coerce_map_from(R):
            return True
        
    def _call_(self, P):
        res = self._ring(0)
        for i in range(len(self._images)):
            res += P.derivative(self._ring.gen(i))*self._images[i]
        return res
            
    def _repr_(self):
        s = ""
        for i in range(len(self._images)):
            c = self._images[i]
            if c == 0:
                continue
            ddx = "d/d%s" %self._ring.gen(i)
            if c == 1:
                s += " + " + ddx
            elif c == -1:
                s += " - " + ddx
            elif c._is_atomic():
                s += " + %s*%s" %(c, ddx)
            elif (-c)._is_atomic():
                s += " - %s*%s" %(-c, ddx)
            else:
                s += " + (%s)*%s" %(c, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        return s



    
    def _add_(self, other):
        if self._ring != other._ring:
            raise TypeError("Rings are not the same")
        im = [self._images[i] + other._images[i] for i in range(len(self._ring.gens()))]
        return RingDerivation_polynomial(self._parent, self._theta, im)

    def scal_mult(self, factor):
        if self._ring != factor.parent():
            raise TypeError("The ring of the factor is not the right ring")
        res = [factor * x  for x in self._images]
        return RingDerivation_polynomial(self._parent, self._theta, res)


class RingDerivation_fraction():
    def __init__(self, parent, image):
        self._parent = parent
        self._images = image

    def __call__(self, frac):
        res = frac.derivative()*self._images
        return res
    def _repr_(self):
        s = ""
        for i in range(len(self._images)):
            c = self._images[i]
            if c == 0:
                continue
            ddx = "d/d%s" %self._ring.gen(i)
            if c == 1:
                s += " + " + ddx
            elif c == -1:
                s += " - " + ddx
            elif c._is_atomic():
                s += " + %s*%s" %(c, ddx)
            elif (-c)._is_atomic():
                s += " - %s*%s" %(-c, ddx)
            else:
                s += " + (%s)*%s" %(c, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        return s
