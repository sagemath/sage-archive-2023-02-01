r"""
Convex Sets
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.abstract_method import abstract_method

class ConvexSet_base(SageObject):

    """
    Abstract base class for convex sets.
    """

    def is_empty(self):
        r"""
        Test whether ``self`` is the empty set

        OUTPUT:

        Boolean.
        """
        return self.dim() < 0

    def is_universe(self):
        r"""
        Test whether ``self`` is the whole ambient space

        OUTPUT:

        Boolean.
        """
        if not self.is_full_dimensional():
            return False
        raise NotImplementedError

    @abstract_method
    def dim(self):
        r"""
        Return the dimension of ``self``.
        """

    @abstract_method
    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.
        """

    def is_full_dimensional(self):
        r"""
        Return whether ``self`` is full dimensional.

        OUTPUT:

        Boolean. Whether the polyhedron is not contained in any strict
        affine subspace.
        """
        return self.dim() == self.ambient_dim()

    @abstract_method
    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        if self.is_open():
            return True
        raise NotImplementedError

    @abstract_method
    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.

        """

    def is_compact(self):
        r"""
        Return whether ``self`` is compact.

        OUTPUT:

        Boolean.

        """
        if not self.is_closed():
            return False
        if self.dimension() < 1:
            return True
        raise NotImplementedError

    def closure(self):
        r"""
        Return the topological closure of ``self``.
        """
        if self.is_closed():
            return self
        raise NotImplementedError

    def interior(self):
        r"""
        Return the topological interior of ``self``.
        """
        if self.is_open():
            return self
        raise NotImplementedError

    @abstract_method(optional=True)
    def affine_hull(self):
        r"""
        Return the affine hull of ``self``.
        """


class ConvexSet_closed(ConvexSet_base):

    r"""
    Abstract base class for closed convex sets.
    """

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.
        """
        return True

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        return self.is_empty() or self.is_universe()


class ConvexSet_compact(ConvexSet_closed):

    r"""
    Abstract base class for compact convex sets.
    """

    def is_universe(self):
        r"""
        Return whether ``self`` is the whole ambient space

        OUTPUT:

        Boolean.
        """
        return self.ambient_dim() == 0 and not self.is_empty()

    def is_compact(self):
        r"""
        Return whether ``self`` is compact.

        OUTPUT:

        Boolean.

        """
        return True


class ConvexSet_relatively_open(ConvexSet_base):

    r"""
    Abstract base class for relatively open sets.
    """

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT:

        Boolean.

        """
        return True

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        return self.is_full_dimensional()
