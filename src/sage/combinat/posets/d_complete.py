# ****************************************************************************
#       Copyright (C) 2020 Stefan Grosser <stefan.grosser1@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import Poset, FinitePoset
from sage.misc.lazy_attribute import lazy_attribute

class DCompletePoset(FinitePoset):
    r"""
    A (finite) `n`-element d-complete poset constructed from a directed acyclic graph.

    INPUT:

    - ``hasse_diagram`` -- an instance of
      :class:`~sage.combinat.posets.posets.FinitePoset`, or a
      :class:`DiGraph` that is transitively-reduced, acyclic,
      loop-free, and multiedge-free.

    - ``elements`` -- an optional list of elements, with ``element[i]``
      corresponding to vertex ``i``. If ``elements`` is ``None``, then it is
      set to be the vertex set of the digraph. Note that if this option is set,
      then ``elements`` is considered as a specified linear extension of the poset
      and the `linear_extension` attribute is set.

    - ``category`` -- :class:`FinitePosets`, or a subcategory thereof.

    - ``facade`` -- a boolean or ``None`` (default); whether the
      :class:`~sage.combinat.posets.posets.DCompletePoset`'s elements should be
      wrapped to make them aware of the Poset they belong to.

      * If ``facade = True``, the
        :class:`~sage.combinat.posets.posets.DCompletePoset`'s elements are exactly
        those given as input.

      * If ``facade = False``, the
        :class:`~sage.combinat.posets.posets.DCompletePoset`'s elements will become
        :class:`~sage.combinat.posets.posets.PosetElement` objects.

      * If ``facade = None`` (default) the expected behaviour is the behaviour
        of ``facade = True``, unless the opposite can be deduced from the
        context (i.e. for instance if a
        :class:`~sage.combinat.posets.posets.DCompletePoset` is built from another
        :class:`~sage.combinat.posets.posets.DCompletePoset`, itself built with
        ``facade = False``)

    - ``key`` -- any hashable value (default: ``None``).

    """
    @lazy_attribute
    def _hooks(self):
        """
        Calculates the hook lengths of the elements of self
        """
        hooks = {}

        min_diamond = {} # Maps max of double-tailed diamond to min of double-tailed diamond
        max_diamond = {} # Maps min of double-tailed diamond to max of double-tailed diamond

        diamonds = [] # Tuples of four elements that are diamonds
        
        diamond_index = {} # Map max elmt of double tailed diamond to index of diamond
        # Find all diamonds
        for e in self._elements:
          covers = self.upper_covers(e)
          for p1, p2 in [(i,j) for i in covers for j in covers if i < j]:
            top = self.common_upper_covers([p1, p2])
            if len(top) != 1:
              continue

            else:
              diamonds.append((e, p1, p2, top[0]))
              min_diamond[top[0]] = e
              max_diamond[e] = top[0]
              diamond_index[top[0]] = len(diamonds) - 1
              break
        # Find all the double-tailed diamonds and map the mins and maxes. 
        for index, d in enumerate(diamonds):
          min_elmt = d[0]
          max_elmt = d[3]

          while True:
            potential_min = self.lower_covers(min_elmt)
            
            potential_max = self.upper_covers(max_elmt)

            # Check if any of these make a longer double tailed diamond
            found_diamond = False
            for (mn, mx) in [(i,j) for i in potential_min for j in potential_max]:
              if len(self.lower_covers(mx)) != 1:
                continue
              if len(self._hasse_diagram.all_paths(self._element_to_vertex(mn), self._element_to_vertex(mx))) == 2:
                # Success
                min_elmt = mn
                max_elmt = mx

                min_diamond[mx] = mn
                max_diamond[mn] = mx
                diamond_index[mx] = index
                found_diamond = True
                break
            if not found_diamond:
              break

        # Compute the hooks
        queue = list(self.minimal_elements())
        enqueued = set()
        while queue:
            elmt = queue.pop(0)
            
            if elmt not in diamond_index:
              hooks[elmt] = len(self.order_ideal([elmt]))
            else:
              diamond = diamonds[diamond_index[elmt]]
              side1 = diamond[1]
              side2 = diamond[2]
              hooks[elmt] = hooks[side1] + hooks[side2] - hooks[min_diamond[elmt]]
            enqueued.add(elmt)

            for c in self.upper_covers(elmt):
              if c not in enqueued:
                queue.append(c)
                enqueued.add(c)

        return hooks
                
    def __init__(self, hasse_diagram, elements, category, facade, key):
        FinitePoset.__init__(self, hasse_diagram=hasse_diagram, elements=elements, category=category, facade=facade, key=key)

    def get_hook(self, elmt):
        """
        Get the hook length of a specific element
        """
        return self._hooks[elmt]

    def get_hooks(self):
        """
        Get all the hook lengths returned in a dictionary
        """
        return dict(self._hooks)

    @staticmethod
    def is_d_complete(cls, poset):
        """
        Check if a poset is d-complete
        """
        pass

    def linear_extensions(self, facade=False):
        from .linear_extensions import LinearExtensionsOfPosetWithHooks
        return LinearExtensionsOfPosetWithHooks(self, facade=facade)