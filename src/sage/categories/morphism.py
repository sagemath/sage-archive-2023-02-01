"""
Morphisms

AUTHORS:
    -- William Stein: initial version
    -- David Joyner (12-17-2005): added examples
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

import homset

import sage.rings.coerce
import sage.rings.arith as arith

from sage.structure.all import Element, Element_cmp_

def is_Morphism(x):
    return isinstance(x, Morphism)

class Morphism(Element, Element_cmp_):
    def __init__(self, parent):
        if not isinstance(parent, homset.Homset):
            raise TypeError, "parent (=%s) must be a Homspace"%parent
        Element.__init__(self, parent)

    def _repr_type(self):
        return "Generic"

    def _repr_defn(self):
        return ""

    def _repr_(self):
        if self.is_endomorphism():
            s = "%s endomorphism of %s"%(self._repr_type(), self.domain())
        else:
            s = "%s morphism:"%self._repr_type()
            s += "\n  From: %s"%self.domain()
            s += "\n  To:   %s"%self.codomain()
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    def domain(self):
        return self.parent().domain()

    def codomain(self):
        return self.parent().codomain()

    def category(self):
        return self.parent().category()

    def is_endomorphism(self):
        return self.parent().is_endomorphism_set()

    def __invert__(self):  # notation in python is (~f) for the inverse of f.
        raise NotImplementedError

    def __call__(self, x):
        try:
            y = self.domain()(x)
        except TypeError:
            raise TypeError, "%s must be coercible into %s"%(x,self.domain())
        return self._call_(y)

    def _call_(self, x):
        raise NotImplementedError

    def __mul__(self, right):
        r"""
        The multiplication * operator is operator composition.

        INPUT:
            self -- Morphism
            right -- Morphism

        OUTPUT:
            The morphism $x \mapsto self(right(x))$.
        """
        if not isinstance(right, Morphism):
            raise TypeError, "right must be a morphism"
        if right.codomain() != self.domain():
            raise TypeError, "self domain must equal right codomain"
        H = homset.Hom(right.domain(), self.codomain(), self.parent().category())
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        return FormalCompositeMorphism(homset, right, self)

    def __pow__(self, n):
        if not self.is_endomorphism():
            raise TypeError, "self must be an endomorphism."
        return arith.generic_power(self, n)

class FormalCoercionMorphism(Morphism):
    def __init__(self, parent):
        Morphism.__init__(self, parent)
        try:
            self.codomain()._coerce_(self.domain().gen(0))
        except TypeError:
            raise TypeError, "Natural coercion morphism from %s to %s not defined."%(self.domain(), self.codomain())

    def _repr_type(self):
        return "Coercion"

    def _call_(self, x):
        return self.codomain()._coerce_(self.domain()._coerce_(x))

class FormalCompositeMorphism(Morphism):
    def __init__(self, parent, first, second):
        Morphism.__init__(self, parent)
        self.__first = first
        self.__second = second

    def _call_(self, x):
        return self.__second(self.__first(x))

    def _repr_type(self):
        return "Composite"

    def _repr_defn(self):
        return "  %s\nthen\n  %s"%(self.__first, self.__second)

    def first(self):
        """
        The first morphism in the formal composition, where the
        composition is x|--> second(first(x)).

        """
        return self.__first

    def second(self):
        """
        The second morphism in the formal composition, where the
        composition is x|--> second(first(x)).
        """
        return self.__second

