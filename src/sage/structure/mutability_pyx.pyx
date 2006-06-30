"""
Mutability Pyrex Implementation
"""

##########################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################

cdef class Mutability:

    def __init__(self, is_immutable=False):
        self._is_immutable = is_immutable

    def _require_mutable(self):
        if self._is_immutable:
            raise ValueError, "object is immutable; please change a copy instead."%self

    def set_immutable(self):
        """
        Make this object immutable, so it can never again be changed.

        EXAMPLES:
            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.set_immutable()
            sage: v[3] = 7
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        self._is_immutable = True

    def is_immutable(self):
        """
        Return True if this object is immutable (can not be changed)
        and False if it is not.

        To make this object immutable use self.set_immutable().

        EXAMPLE:
            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.is_immutable()
            False
            sage: v.set_immutable()
            sage: v.is_immutable()
            True
        """
        try:
            return self._is_immutable
        except AttributeError:
            return False

    def is_mutable(self):
        try:
            return not self._is_immutable
        except AttributeError:
            return True
