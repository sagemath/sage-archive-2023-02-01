"""
Constructor for permutations

This module contains the generic constructor to build element of the symmetric
groups (or more general permutation groups) called
:class:`~sage.groups.perm_gps.permgroup_element.PermutationGroupElement`. These
objects have a more group theoretic flavor than the more combinatorial
:class:`~sage.combinat.permutation.Permutation`.
"""
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Joyner
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from . import permgroup_element

def PermutationGroupElement(g, parent=None, check=True):
    r"""
    Builds a permutation from ``g``.

    INPUT:

    - ``g`` -- either

      - a list of images

      - a tuple describing a single cycle

      - a list of tuples describing the cycle decomposition

      - a string describing the cycle decomposition

    - ``parent`` -- (optional) an ambient permutation group for the result;
      it is mandatory if you want a permutation on a domain different
      from `\{1, \ldots, n\}`

    - ``check`` -- (default: ``True``) whether additional check are performed;
      setting it to ``False`` is likely to result in faster code

    EXAMPLES:

    Initialization as a list of images::

        sage: p = PermutationGroupElement([1,4,2,3])
        sage: p
        (2,4,3)
        sage: p.parent()
        Symmetric group of order 4! as a permutation group

    Initialization as a list of cycles::

        sage: p = PermutationGroupElement([(3,5),(4,6,9)])
        sage: p
        (3,5)(4,6,9)
        sage: p.parent()
        Symmetric group of order 9! as a permutation group

    Initialization as a string representing a cycle decomposition::

        sage: p = PermutationGroupElement('(2,4)(3,5)')
        sage: p
        (2,4)(3,5)
        sage: p.parent()
        Symmetric group of order 5! as a permutation group

    By default the constructor assumes that the domain is `\{1, \dots, n\}`
    but it can be set to anything via its second ``parent`` argument::

        sage: S = SymmetricGroup(['a', 'b', 'c', 'd', 'e'])
        sage: PermutationGroupElement(['e', 'c', 'b', 'a', 'd'], S)
        ('a','e','d')('b','c')
        sage: PermutationGroupElement(('a', 'b', 'c'), S)
        ('a','b','c')
        sage: PermutationGroupElement([('a', 'c'), ('b', 'e')], S)
        ('a','c')('b','e')
        sage: PermutationGroupElement("('a','b','e')('c','d')", S)
        ('a','b','e')('c','d')

    But in this situation, you might want to use the more direct::

        sage: S(['e', 'c', 'b', 'a', 'd'])
        ('a','e','d')('b','c')
        sage: S(('a', 'b', 'c'))
        ('a','b','c')
        sage: S([('a', 'c'), ('b', 'e')])
        ('a','c')('b','e')
        sage: S("('a','b','e')('c','d')")
        ('a','b','e')('c','d')
    """
    if isinstance(g, permgroup_element.PermutationGroupElement):
        if parent is None or g.parent() is parent:
            return g

    if parent is None:
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.combinat.permutation import from_cycles

        try:
            v = permgroup_element.standardize_generator(g, None)
        except KeyError:
            raise ValueError("invalid permutation vector: %s" % g)
        degree = max([1] + [max(cycle+(1,)) for cycle in v])
        v = from_cycles(degree, v)
        parent = SymmetricGroup(len(v))
        # We have constructed the parent from the element and already checked
        #   that it is a valid permutation
        check = False

    return parent.element_class(g, parent, check)

