"""
Group, ring, etc. actions on objects.

The terminology and notation used is suggestive of groups
acting on sets, but this framework can be used for modules,
algebras, etc.

A group action $G \times S \rightarrow S$ is a functor from $G$ to Sets.

AUTHORS:
    -- Robert Bradshaw: initial version
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functor cimport Functor
from morphism cimport Morphism

import homset

include "../ext/stdsage.pxi"

cdef class Action(Functor):

    def __init__(self, G, S, bint is_left = 1):
        from category_types import Groupoid
        Functor.__init__(self, Groupoid(G), S.category())
        self._G = G
        self._S = S
        self._is_left = is_left

    def _apply_functor(self, x):
        return self(x)

    def __call__(self, *args):
        if len(args) == 1:
            g = args[0]
            if g in self._G:
                return ActionEndomorphism(self, self._G(g))
            elif g == self._G:
                return self._S
            else:
                raise TypeError, "%s not an element of %s"%(g, self._G)
        elif len(args) == 2:
            if self._is_left:
                return self._call_c(self._G(args[0]), self._S(args[1]))
            else:
                return self._call_c(self._S(args[0]), self._G(args[1]))

    def _call_(self, a, b):
        return self._call_c_impl(a, b)

    cdef Element _call_c(self, a, b):
        if HAS_DICTIONARY(self):
            return self._call_(a, b)
        else:
            return self._call_c_impl(a, b)

    cdef Element _call_c_impl(self, Element a, Element b):
        raise NotImplementedError, "Action not implemented."

    def act(self, g, a):
        """
        This is a consistant interface for acting on a by g,
        irregardless of whether its a left or right action.
        """
        if self._is_left:
            return self._call_c(g, a)
        else:
            return self._call_c(a, g)

    def is_left(self):
        return self._is_left

    def __repr__(self):
        side = "Left" if self._is_left else "Right"
        return "%s action of %r on %r"%(side, self._G, self._S)

cdef class ActionEndomorphism(Morphism):

    def __init__(self, Action action, g):
        Morphism.__init__(self, homset.Hom(action._S, action._S))
        self._action = action
        self._g = g

    cdef Element _call_c_impl(self, Element x):
        if self._action._is_left:
            return self._action._call_c(self._g, x)
        else:
            return self._action._call_c(x, self._g)

    def _repr_(self):
        return "Action of %s on %s under %s."%(self._g, self._action._S, self._action)

    def __mul__(left, right):
        cdef ActionEndomorphism left_c, right_c
        if PY_TYPE_CHECK(left, ActionEndomorphism) and PY_TYPE_CHECK(right, ActionEndomorphism):
            left_c = left
            right_c = right
            if left_c._action is right_c._action:
                if left_c._action._is_left:
                    return ActionEndomorphism(left_c._action, left_c._g * right_c._g)
                else:
                    return ActionEndomorphism(left_c._action, right_c._g * left_c._g)
        return Morphism.__mul__(left, right)

    def __invert__(self):
            return ActionEndomorphism(self._action, ~self._g)




