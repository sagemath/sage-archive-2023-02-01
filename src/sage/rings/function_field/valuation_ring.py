r"""
Valuation rings

This module provides valuation rings at places of function fields.

EXAMPLES::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: p = L.places_finite()[0]
    sage: R = p.valuation_ring()
    sage: R
    Valuation ring at Place (x, x*y)

AUTHORS:

- Kwankyu Lee (2016): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.structure.parent import Parent

from sage.categories.homset import Hom
from sage.categories.rings import Rings

from sage.modules.free_module_element import vector
from sage.rings.morphism import RingHomomorphism

lazy_import('sage.matrix.constructor', 'matrix')

class FunctionFieldValuationRing(Parent):
    """
    """
    def __init__(self, field, place):
        Parent.__init__(self, category=Rings())

        self._field = field
        self._place = place

    def _element_constructor_(self, x):
        """
        Construct an element of the function field belonging to this
        valuation ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: y in R
            False
            sage: 1/y in R
            True
            sage: x + y in R
            False
            sage: 1/(x + y) in R
            True
        """
        x = self._field(x)
        if x.valuation(self._place) >= 0:
            return x
        else:
            raise TypeError

    def _repr_(self):
        return 'Valuation ring at {}'.format(self._place)

    def place(self):
        """
        Return the place corresponding to the valuation ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: v = p.valuation_ring()
            sage: p == v.place()
            True
        """
        return self._place

class FunctionFieldValuationRing_global(FunctionFieldValuationRing):
    """
    """
    def residue_field(self):
        """
        Return the residue field of the valuation ring along with
        the maps from and to it.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: k, fr_k, to_k = R.residue_field()
            sage: k
            Finite Field of size 2
            sage: fr_k
            Ring morphism:
              From: Finite Field of size 2
              To:   Valuation ring at Place (x, x*y)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x, x*y)
              To:   Finite Field of size 2
            sage: to_k(y)
            Traceback (most recent call last):
            ...
            TypeError: ...
            sage: to_k(1/y)
            0
            sage: to_k(y/(1+y))
            1
        """
        k, from_k, to_k = self._place._residue_field()
        mor_from_k = Morphism(Hom(k,self), from_k)
        mor_to_k = Morphism(Hom(self,k), to_k)
        return k, mor_from_k, mor_to_k

class Morphism(RingHomomorphism):
    def __init__(self, parent, from_def):
        RingHomomorphism.__init__(self, parent)

        self._map = from_def

    def _call_(self, x):
        return self._map(x)

