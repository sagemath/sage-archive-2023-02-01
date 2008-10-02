r"""
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
from map cimport Map

import homset
import sage.structure.element

include "../ext/stdsage.pxi"


cdef class Action(Functor):

    def __init__(self, G, S, bint is_left = 1, op=None):
        from category_types import Groupoid
        Functor.__init__(self, Groupoid(G), S.category())
        self.G = G
        self.S = S
        self._is_left = is_left
        self.op = op

    def _apply_functor(self, x):
        return self(x)

    def __call__(self, *args):
        if len(args) == 1:
            g = args[0]
            if g in self.G:
                return ActionEndomorphism(self, self.G(g))
            elif g == self.G:
                return self.S
            else:
                raise TypeError, "%s not an element of %s"%(g, self.G)
        elif len(args) == 2:
            if self._is_left:
                return self._call_(self.G(args[0]), self.S(args[1]))
            else:
                return self._call_(self.S(args[0]), self.G(args[1]))

    cpdef Element _call_(self, a, b):
        raise NotImplementedError, "Action not implemented."

    def act(self, g, a):
        """
        This is a consistent interface for acting on a by g,
        regardless of whether it's a left or right action.
        """
        if self._is_left:
            return self._call_(g, a)
        else:
            return self._call_(a, g)

    def __invert__(self):
        return InverseAction(self)

    def is_left(self):
        return self._is_left

    def __repr__(self):
        side = "Left" if self._is_left else "Right"
        return "%s %s by %r on %r"%(side, self._repr_name_(), self.G, self.S)

    def _repr_name_(self):
        return "action"

    def actor(self):
        return self.G

    def codomain(self):
        return self.S

    def domain(self):
        return self.S

    def left_domain(self):
        if self._is_left:
            return self.G
        else:
            return self.domain()

    def right_domain(self):
        if self._is_left:
            return self.domain()
        else:
            return self.G

    def operation(self):
        return self.op


cdef class InverseAction(Action):
    """
    An action that acts as the inverse of the given action.

    TESTS:
    This illustrates a shortcoming in the current coercion model.
    See the comments in _call_ below.

        sage: x = polygen(QQ,'x')
        sage: a = 2*x^2+2; a
        2*x^2 + 2
        sage: a / 2
        x^2 + 1
        sage: a /= 2
        sage: a
        x^2 + 1
    """
    def __init__(self, Action action):
        G = action.G
        try:
            from sage.groups.group import Group
            # We must be in the case that parent(~a) == parent(a)
            # so we can invert in call_c code below.
            if (PY_TYPE_CHECK(G, Group) and G.is_multiplicative()) or G.is_field():
                Action.__init__(self, G, action.S, action._is_left)
                self._action = action
                return
            else:
                K = G._pseudo_fraction_field()
                Action.__init__(self, K, action.S, action._is_left)
                self._action = action
                return
        except (AttributeError, NotImplementedError):
            pass
        raise TypeError, "No inverse defined for %r." % action

    cpdef Element _call_(self, a, b):
        if self._action._is_left:
            if self.S_precomposition is not None:
                b = self.S_precomposition(b)
            return self._action._call_(~a, b)
        else:
            if self.S_precomposition is not None:
                a = self.S_precomposition(a)
            return self._action._call_(a, ~b)

    def codomain(self):
        return self._action.codomain()

    def __invert__(self):
        return self._action

    def _repr_name_(self):
        return "inverse action"

cdef class PrecomposedAction(Action):

    def __init__(self, Action action, Map left_precomposition, Map right_precomposition):
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
            Action.__init__(self, left, action.S, 1)
        else:
            Action.__init__(self, right, action.S, 0)
        self._action = action
        self.left_precomposition = left_precomposition
        self.right_precomposition = right_precomposition

    cpdef Element _call_(self, a, b):
        if self.left_precomposition is not None:
            a = self.left_precomposition._call_(a)
        if self.right_precomposition is not None:
            b = self.right_precomposition._call_(b)
        return self._action._call_(a, b)

    def domain(self):
        if self._is_left and self.right_precomposition is not None:
            return self.right_precomposition.domain()
        elif not self._is_left and self.left_precomposition is not None:
            return self.left_precomposition.domain()
        else:
            return self._action.domain()

    def codomain(self):
        return self._action.codomain()

    def __invert__(self):
        return PrecomposedAction(~self._action, self.left_precomposition, self.right_precomposition)

    def __repr__(self):
        s = repr(self._action)
        if self.left_precomposition is not None:
            s += "\nwith precomposition on left by %r" % self.left_precomposition
        if self.right_precomposition is not None:
            s += "\nwith precomposition on right by %r" % self.right_precomposition
        return s


cdef class ActionEndomorphism(Morphism):

    def __init__(self, Action action, g):
        Morphism.__init__(self, homset.Hom(action.S, action.S))
        self._action = action
        self._g = g

    cpdef Element _call_(self, x):
        if self._action._is_left:
            return self._action._call_(self._g, x)
        else:
            return self._action._call_(x, self._g)

    def _repr_(self):
        return "Action of %s on %s under %s."%(self._g, self._action.S, self._action)

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


