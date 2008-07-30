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

cdef class Morphism(Map):

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

    def category(self):
        return self.parent().category()

    def is_endomorphism(self):
        return self.parent().is_endomorphism_set()

    def __invert__(self):  # notation in python is (~f) for the inverse of f.
        raise NotImplementedError

    def pushforward(self, I):
        raise NotImplementedError

cdef class FormalCoercionMorphism(Morphism):
    def __init__(self, parent):
        Morphism.__init__(self, parent)
        if not self.codomain().has_coerce_map_from(self.domain()):
            raise TypeError, "Natural coercion morphism from %s to %s not defined."%(self.domain(), self.codomain())

    def _repr_type(self):
        return "Coercion"

    cpdef Element _call_(self, x):
        return self._codomain.coerce(x)

cdef class CallMorphism(Morphism):

    def _repr_type(self):
        return "Call"

    cpdef Element _call_(self, x):
        return self._codomain(x)

cdef class IdentityMorphism(Morphism):

    def __init__(self, parent):
        if not isinstance(parent, homset.Homset):
            parent = homset.Hom(parent, parent)
        Morphism.__init__(self, parent)

    def _repr_type(self):
        return "Identity"

    cpdef Element _call_(self, x):
        return x

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        if len(args) == 0 and len(kwds) == 0:
            return x
        elif self._codomain._element_init_pass_parent:
            return self._codomain._element_class(self._codomain, x, *args, **kwds)
        else:
            return self._codomain._element_class(x, *args, **kwds)

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

cdef class SetMorphism(Morphism):
    def __init__(self, parent, function):
        """
        INPUT:
        parent -- a Homset
        function -- a Python function that takes elements of the domain as input
                    and returns elements of the domain.
        """
        self._function = function

    cpdef Element _call_(self, x):
        return self._function(x)

