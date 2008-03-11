"""
Morphisms

AUTHORS:
    -- William Stein: initial version
    -- David Joyner (12-17-2005): added examples
    -- Robert Bradshaw (2007-06-25) Pyrexification
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

include "../ext/stdsage.pxi"
from sage.structure.element cimport Element
from sage.structure.element import generic_power

def make_morphism(_class, parent, _dict, _slots):
    # from element.pyx
    cdef Morphism mor = _class.__new__(_class)
    mor._set_parent(parent)
    mor._update_slots(_slots)
    if HAS_DICTIONARY(mor):
        mor.__dict__ = _dict
    return mor

def is_Morphism(x):
    return isinstance(x, Morphism)

cdef class Morphism(Element):
    def __init__(Morphism self, parent):
        if not isinstance(parent, homset.Homset):
            raise TypeError, "parent (=%s) must be a Homspace"%parent
        Element.__init__(self, parent)
        self._domain = parent.domain()
        self._codomain = parent.codomain()

    cdef _update_slots(self, _slots):
        self._domain = _slots['_domain']
        self._codomain = _slots['_codomain']

    def _test_update_slots(self, _slots):
        self._update_slots(_slots)

    cdef _extra_slots(self, _slots):
        _slots['_domain'] = self._domain
        _slots['_codomain'] = self._codomain
        return _slots

    def _test_extra_slots(self, _slots):
        return self._extra_slots(_slots)

    def __reduce__(self):
        if HAS_DICTIONARY(self):
            _dict = self.__dict__
        else:
            _dict = {}
        return make_morphism, (self.__class__, self._parent, _dict, self._extra_slots({}))

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
        return self._domain

    def codomain(self):
        return self.parent().codomain()

    def category(self):
        return self.parent().category()

    def is_endomorphism(self):
        return self.parent().is_endomorphism_set()

    def __invert__(self):  # notation in python is (~f) for the inverse of f.
        raise NotImplementedError

    def __call__(self, x):
        """
        Apply this morphism to x.

        INPUT:
            x -- an element coercible to self; also objects like
                 ideals are supported in some cases

        OUTPUT:
            an element (or ideal, etc.)

        EXAMPLES:
            sage: R.<x,y> = QQ[]; phi=R.hom([y,x])
            sage: phi(y)
            x

        We take the image of an ideal:
            sage: I = ideal(x,y); I
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: phi(I)
            Ideal (y, x) of Multivariate Polynomial Ring in x, y over Rational Field
        """
        if not PY_TYPE_CHECK(x, Element):
            try:
                return self._call_c(x)
            except TypeError:
                raise TypeError, "%s must be coercible into %s (and is not an element)"%(x, self._domain)
        elif (<Element>x)._parent is not self._domain:
            try:
                x = self._domain(x)
            except TypeError:
                try:
                    return self.pushforward(x)
                except TypeError:
                    raise TypeError, "%s must be coercible into %s"%(x, self._domain)
        return self._call_c(x)

    def _call_(self, x):
        return self._call_c_impl(x)

    cdef Element _call_c(self, x):
        if HAS_DICTIONARY(self):
            return self._call_(x)
        else:
            return self._call_c_impl(x)

    cdef Element _call_c_impl(self, Element x):
        raise NotImplementedError

    def pushforward(self, I):
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
            raise TypeError, "right (=%s) must be a morphism to multiply it by %s"%(right, self)
        if right.codomain() != self.domain():
            raise TypeError, "self (=%s) domain must equal right (=%s) codomain"%(self, right)
        H = homset.Hom(right.domain(), self.codomain(), self.parent().category())
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        return FormalCompositeMorphism(homset, right, self)

    def __pow__(self, n, dummy):
        if not self.is_endomorphism():
            raise TypeError, "self must be an endomorphism."
        # todo -- what about the case n=0 -- need to specify the identity map somehow.
        return generic_power(self, n)

cdef class Section(Morphism):
    def __init__(self, morphism):
        from sage.categories.homset import Hom
        from sage.categories.category_types import SetsWithPartialMaps
        Morphism.__init__(self, Hom(morphism.codomain(), morphism.domain(), SetsWithPartialMaps()))
        self._morphism = morphism

    def _repr_type(self):
        return "Section"

cdef class FormalCoercionMorphism(Morphism):
    def __init__(self, parent):
        Morphism.__init__(self, parent)
        if not self.codomain().has_coerce_map_from(self.domain()):
            raise TypeError, "Natural coercion morphism from %s to %s not defined."%(self.domain(), self.codomain())

    def _repr_type(self):
        return "Coercion"

    # We need to override _call_c in this special case so that FormalCoercionMorphisms can operate on things that are not elements.
    cdef Element _call_c(self, x):
        return self._codomain._coerce_(x)

cdef class CallMorphism(Morphism):

    def _repr_type(self):
        return "Call"

    # We need to override _call_c in this special case so that CallMorphisms can operate on things that are not elements.
    cdef Element _call_c(self, x):
        return self._codomain(x)

cdef class IdentityMorphism(Morphism):

    def __init__(self, parent):
        if not isinstance(parent, homset.Homset):
            parent = homset.Hom(parent, parent)
        Morphism.__init__(self, parent)

    def _repr_type(self):
        return "Identity"

    cdef Element _call_c_impl(self, Element x):
        return x

    def __mul__(left, right):
        if not isinstance(right, Morphism):
            raise TypeError, "right (=%s) must be a morphism to multiply it by %s"%(right, left)
        if not isinstance(left, Morphism):
            raise TypeError, "left (=%s) must be a morphism to multiply it by %s"%(left, right)
        if right.codomain() != left.domain():
            raise TypeError, "self (=%s) domain must equal right (=%s) codomain"%(left, right)
        if isinstance(left, IdentityMorphism):
            return right
        else:
            return left

    def __pow__(self, n, dummy):
        return self

    def __invert__(self):
        return self


cdef class FormalCompositeMorphism(Morphism):
    def __init__(self, parent, first, second):
        Morphism.__init__(self, parent)
        self.__first = first
        self.__second = second

    cdef _update_slots(self, _slots):
        self.__first = _slots['__first']
        self.__second = _slots['__second']
        Morphism._update_slots(self, _slots)

    cdef _extra_slots(self, _slots):
        _slots['__first'] = self.__first
        _slots['__second'] = self.__second
        return Morphism._extra_slots(self, _slots)

    cdef Element _call_c_impl(self, Element x):
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

