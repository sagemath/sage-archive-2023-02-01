r"""
Homsets with libgap.

Homsets for libgap based groups

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Simon Brandhorst (2018-02-08): initial version

"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.group_homset import GroupHomset_generic
from sage.groups.libgap_morphism import LibGAPGroupMorphism
from sage.categories.groups import Groups

class GroupHomset(GroupHomset_generic):

    def __init__(self, G, H, category=Groups()):
        r"""
        Return the homset of two matrix groups.

        INPUT:

        - ``G`` -- a matrix group

        - ``H`` -- a matrix group

        OUTPUT:

        The homset of two matrix groups.

        EXAMPLES::

        """
        from sage.categories.homset import Homset
        Homset.__init__(self, G, H, category)

    def __call__(self, im_gens, check=True):
        """
        Return the homomorphism defined by images of generators.

        INPUT:

        - ``im_gens`` -- iterable, the list of images of the
          generators of the domain

        - ``check`` -- bool (optional, default: ``True``), whether to
          check if images define a valid homomorphism

        OUTPUT:

        Group homomorphism.

        EXAMPLES::


        """
        from sage.groups.libgap_morphism import LibGAPGroupMorphism
        try:
            return LibGAPGroupMorphism(self, im_gens, check=check)
        except (NotImplementedError, ValueError) as err:
            raise TypeError('images do not define a group homomorphism')
