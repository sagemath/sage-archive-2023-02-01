"""
Add doctest here
"""

#############################################################################
#    Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

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
        if twist is not None:
            if not (isinstance(twist, Map) and twist.category_for().is_subcategory(Rings())):
                raise TypeError("The twisting homorphism must be an homomorphism of rings")
            if twist.domain() is not domain:
                map = domain.coerce_map_from(twist.domain())
                if map is None:
                    raise TypeError("The domain of the derivation must coerce to the domain of the twisting homomorphism")
                twist = twist * map
            if twist.codomain() is not codomain:
                map = twist.codomain().coerce_map_from(codomain)
                if map is None:
                    raise TypeError("The codomain of the twisting homomorphism must coerce to the codomain of the derivation")
                twist = map * twist
        self._domain = domain
        self._codomain = codomain
        self._twist = twist
        if element_class is None:
            if twist is None:
                self.element_class = RingDerivationWithoutTwist_im_gens
            else:
                self.element_class = RingDerivationWithTwist_generic
        else:
            self.element_class = element_class
        Module.__init__(self, codomain)

    def _repr_(self):
        if self._twist is None:
            s = "Module of derivations"
            t = ""
        else:
            s = "Module of twisted derivations"
            t = " (twisting morphism: %s)" % self._twist._repr_short()
        if self._domain is self._codomain:
            s += " over %s" % self._domain
        else:
            s += " from %s to %s" % (self._domain, self._codomain)
        return s + t

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def twisting_homomorphism(self):
        if self._twist is None:
            return self._codomain.coerce_map_from(self._domain)
        else:
            return self._twist

    def ngens(self):
        if self._twist is not None:
            raise NotImplementedError("Generators are not implemented for twisted derivations")
        return self.domain().ngens()

    def gens(self):
        if self._twist is not None:
            raise NotImplementedError("Generators are not implemented for twisted derivations")
        domain = self.domain()
        return tuple([ self(x) for x in domain.gens() ])

    def gen(self, n=0):
        if self._twist is not None:
            raise NotImplementedError("Generators are not implemented for twisted derivations")
        return self(self.domain().gen(n))


class RingDerivation(ModuleElement):
    """
    Derivation of rings
    """
    def __call__(self, x):
        arg = self.parent().domain()(x)
        return self._call_(arg)


class RingDerivationWithoutTwist_im_gens(RingDerivation):
    def __init__(self, parent, arg):
        twist = parent.twisting_homomorphism()
        if not twist.is_identity():
            raise TypeError("This class only supports untwisted derivations")
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
        codomain = parent.codomain()
        self._images = [ codomain(x) for x in images ]
        if len(self._images) != parent.domain().ngens():
            raise ValueError("Number of images is incorrect")
        RingDerivation.__init__(self, parent)

    def _coerce_map_from_(self, R):
        if self.base_ring().has_coerce_map_from(R):
            return True
        
    def _repr_(self):
        s = ""
        domain = self.parent().domain()
        for i in range(len(self._images)):
            c = self._images[i]
            sc = str(c)
            if sc == "0":
                continue
            ddx = "d/d%s" % domain.gen(i)
            if sc == "1":
                s += " + " + ddx
            elif sc == "-1":
                s += " - " + ddx
            elif c._is_atomic():
                s += " + %s*%s" % (sc, ddx)
            elif (-c)._is_atomic():
                s += " - %s*%s" % (-c, ddx)
            else:
                s += " + (%s)*%s" % (sc, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        else:
            return s

    def _add_(self, other):
        im = [ self._images[i] + other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.parent()(im)

    def _sub_(self, other):
        im = [ self._images[i] - other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.parent()(im)

    def _rmul_(self, factor):
        factor = self.parent().codomain()(factor)
        im = [ factor*x  for x in self._images ]
        return self.parent()(im)

    def _lmul_(self, factor):
        return self._rmul_(factor)

    def _call_(self, P):
        res = self.parent().codomain()(0)
        domain = self.parent().domain()
        for i in range(len(self._images)):
            res += P.derivative(domain.gen(i)) * self._images[i]
        return res


class RingDerivationWithTwist_generic(RingDerivation):
    def __init__(self, parent, scalar=0):
        codomain = parent.codomain()
        self._scalar = codomain(scalar)
        RingDerivation.__init__(self, parent)

    def _repr_(self):
        scalar = self._scalar
        sc = str(scalar)
        if sc == "0":
            return "0"
        if sc == "1":
            s = ""
        elif sc == "-1":
            s = "-"
        elif scalar._is_atomic():
            s = "%s*" % sc
        elif (-scalar)._is_atomic():
            s = "-%s*" % (-c)
        else:
            s = "(%s)*" % sc
        return s + "([%s] - id)" % self.parent().twisting_homomorphism()._repr_short();

    def _add_(self, other):
        return self.parent()(self._scalar + other._scalar)

    def _sub_(self, other):
        return self.parent()(self._scalar - other._scalar)

    def _rmul_(self, factor):
        return self.parent()(factor * self._scalar)

    def _lmul_(self, factor):
        return self._rmul_(factor)

    def _call_(self, x):
        return self._scalar * (self.parent().twisting_homomorphism()(x) - x)
