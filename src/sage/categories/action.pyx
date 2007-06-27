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
import sage.structure.element

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

    def __invert__(self):
        return InverseAction(self)

    def is_left(self):
        return self._is_left

    def __repr__(self):
        side = "Left" if self._is_left else "Right"
        return "%s action of %r on %r"%(side, self._G, self._S)

    def actor(self):
        return self._G

    def codomain(self):
        return self._S

    def left_domain(self):
        if self._is_left:
            return self._G
        else:
            return self._S

    def right_domain(self):
        if self._is_left:
            return self._S
        else:
            return self._G


cdef class InverseAction(Action):
    """
    An action whose acts as the inverse of the given action.
    """
    def __init__(self, Action action):
        Action.__init__(action._G, action._S, action._is_left)
        self._action = action

    cdef Element _call_c(self, a, b):
        if self._action._is_left:
            return self._action._call_c(~a, b)
        else:
            return self._action._call_c(a, ~b)

    def __invert__(self):
        return self._action


cdef class PrecomposedAction(Action):

    def __init__(self, Action action, Morphism left_precomposition, Morphism right_precomposition):
        left = action.left_domain()
        right = action.right_domain()
        if left_precomposition is not None:
            if left_precomposition._codomain is not left:
                left_precomposition = homset.Hom(left_precomposition._codomain, left).natural_map() * left_precomposition
            left = left_precomposition._domain
        if right_precomposition is not None:
            if right_precomposition._codomain is not right:
              right_precomposition = homset.Hom(right_precomposition._codomain, right).natural_map() * right_precomposition
            right = right_precomposition._domain
        if action._is_left:
            G = left
            S = right
        else:
            G = right
            S = left
        Action.__init__(G, S, action._is_left)
        self._action = action
        self._left_precomposition = left_precomposition
        self._right_precomposition = right_precomposition

    cdef Element _call_c(self, a, b):
        if self._left_precomposition is not None:
            a = self._left_precomposition._call_c(a)
        if self._right_precomposition is not None:
            b = self._right_precomposition._call_c(b)
        return self._action._call_c(a, b)


cdef class ActionEndomorphism(Morphism):

    def __init__(self, Action action, g):
        Morphism.__init__(self, homset.Hom(action._S, action._S))
        self._action = action
        self._g = g

    cdef Element _call_c(self, x):
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
            inv_g = ~self._g
            if sage.structure.element.parent(inv_g) is sage.structure.element.parent(self._g):
                return ActionEndomorphism(self._action, inv_g)
            else:
                return (~self._action)(self._g)


