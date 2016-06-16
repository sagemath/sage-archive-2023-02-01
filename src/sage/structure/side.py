"""
This module implements a class ``Side`` with only two instances
``Left`` and ``Right``.

EXAMPLES::

Variables ``Left`` and ``Right`` are predefined::

    sage: Left
    Left
    sage: Right
    Right
    sage: type(Right)
    <class 'sage.structure.side.Side'>

If we loose them, we can recover them as follows::

    sage: Left = 1
    sage: Left
    1
    sage: Left = sage.structure.side.Side("left")
    sage: Left
    Left
    sage: type(Left)
    <class 'sage.structure.side.Side'>

``Left`` and ``Right`` are unique: if we try to create two times,
we get the same object::

    sage: x = sage.structure.side.Side("right")
    sage: y = sage.structure.side.Side("right")
    sage: x is y
    True
    sage: id(x) == id(y)
    True

A method ``opposite`` is implemented. It returns the opposite side::

    sage: Left.opposite()
    Right

AUTHOR:

- Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage_object import SageObject
from sage.misc.classcall_metaclass import ClasscallMetaclass

class Side(SageObject):
    """
    This class has only two instances: ``Left`` and ``Right``.
    """
    __metaclass__ = ClasscallMetaclass
    _left = None
    _right = None

    @staticmethod
    def __classcall__(cls,side):
        """
        We override the ``__classcall__`` method in order to ensure
        unique representation.
        """
        if isinstance(side,str):
            side = side.lower()
            if side == "left":
                if cls._left is None:
                    cls._left = cls.__new__(Side)
                    cls._left._repr = "Left"
                return cls._left
            elif side == "right":
                if cls._right is None:
                    cls._right = cls.__new__(Side)
                    cls._right._repr = "Right"
                return cls._right
        raise ValueError("side must be 'left' or 'right'")

    def _repr_(self):
        """
        Return a string representation of this side

        EXAMPLES::

            sage: s = Left
            sage: s._repr_()
            'Left'
        """
        return self._repr

    def opposite(self):
        """
        Return the opposite side

        EXAMPLES::

            sage: Left.opposite()
            Right
            sage: Right.opposite()
            Left
        """
        if self is self._left:
            return Side("right")
        elif self is self._right:
            return Side("left")

    def __neg__(self):
        """
        Return the opposite side

        EXAMPLES::

            sage: -Left
            Right
            sage: -Right
            Left
        """
        return self.opposite()

    def __invert__(self):
        """
        Return the opposite side

        EXAMPLES::

            sage: ~Left
            Right
            sage: ~Right
            Left
        """
        return self.opposite()


Left = Side("left")
Right = Side("right")
