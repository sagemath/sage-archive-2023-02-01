from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.module import Module
from sage.structure.element import ModuleElement

from sage.categories.map import Map
from sage.categories.all import Rings, Algebras
from sage.rings.morphism import RingMap, RingHomomorphism



class RingDerivationModule(Module, UniqueRepresentation):
    def __init__(self, domain, codomain, twist=None, element_class=None):
        if not domain in Rings().Commutative():
            raise TypeError("The domain must be a commutative ring")
        if not (codomain in Rings().Commutative() and codomain.has_coerce_map_from(domain)):
            raise TypeError("The codomain must be an algebra over the domain")
        if twist is None:
            twist = codomain.coerce_map_from(domain)
        elif (not (isinstance(twist, Map) and twist.category_for().is_subcategory(Rings()))
            or twist.domain() is not domain or twist.codomain() is not codomain):
            raise TypeError("The twisting morphism must be an homomorphism of rings")
        self._domain = domain
        self._codomain = codomain
        self._twist = twist
        if element_class is None:
            self.element_class = RingDerivation_function
        else:
            self.element_class = element_class
        Module.__init__(self, codomain)

    def _repr_(self):
        if self._twist.is_identity():
            s = "Module of derivations"
            t = ""
        else:
            s = "Module of twisted derivations"
            t = "\nTwisting morphism: %s" % self._twist
        if self._domain is self._codomain:
            s += " over %s" % self._domain
        else:
            s += " from %s to %s" % (self._domain, self._codomain)
        return s + t

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def twist(self):
        return self._twist

    def ngens(self):
        return self.domain().ngens()

    def gens(self):
        domain = self.domain()
        return tuple([ self(x) for x in domain.gens() ])

    def gen(self, n=0):
        return self(self.domain().gen(n))


class RingDerivation(ModuleElement):
    """
    Derivation of rings
    """
    def __call__(self, x):
        arg = self.parent().domain()(x)
        return self._call_(arg)


class RingDerivation_im_gens(RingDerivation):
    def __init__(self, parent, images):
        RingDerivation.__init__(self, parent)
        codomain = parent.codomain()
        self._images = [ codomain(x) for x in images ]
        if len(self._images) != parent.domain().ngens():
            raise ValueError("Number of images is incorrect")

    def _add_(self, other):
        im = [ self._images[i] + other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.__class__(self.parent(), im)

    def _sub_(self, other):
        im = [ self._images[i] - other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.__class__(self.parent(), im)

    def _rmul_(self, factor):
        factor = self.parent().codomain()(factor)
        im = [ factor*x  for x in self._images ]
        return self.__class__(self.parent(), im)

    def _lmul_(self, factor):
        return self._rmul_(factor)


class RingDerivation_function(RingDerivation_im_gens):
    def __init__(self, parent, arg):
        twist = parent.twist()
        if not twist.is_identity():
            raise NotImplementedError
        domain = parent.domain()
        ngens = domain.ngens()
        if isinstance(arg, (list, tuple)):
            images = arg
        else:
            for i in range(ngens):
                if arg == domain.gen(i):
                    images = ngens * [0]
                    images[i] = 1
                    break
            else:
                raise ValueError("You must give a generator or a list of scalars")
        RingDerivation_im_gens.__init__(self, parent, images)
        
    def _coerce_map_from_(self, R):
        if self.base_ring().has_coerce_map_from(R):
            return True
        
    def _call_(self, P):
        res = self.parent().codomain()(0)
        domain = self.parent().domain()
        for i in range(len(self._images)):
            res += P.derivative(domain.gen(i)) * self._images[i]
        return res
            
    def _repr_(self):
        s = ""
        domain = self.parent().domain()
        for i in range(len(self._images)):
            c = self._images[i]
            if c == 0:
                continue
            ddx = "d/d%s" % domain.gen(i)
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
