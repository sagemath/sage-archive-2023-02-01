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
        queue = list(self.minimal_elements())
        
        while queue:
            elmt = queue.pop(0)
            print(elmt)
            children = self.lower_covers(elmt)
            
            if len(children) != 2:
                hooks[elmt] = sum(hooks[child] for child in children) + 1
            
            else: 
                child_meets = set(self.lower_covers(children[0])).intersection(set(self.lower_covers(children[1])))
                if len(child_meets) == 1:
                    hooks[elmt] = hooks[children[0]] + hooks[children[1]] - hooks[child_meets.pop()]

                else:
                    hooks[elmt] = sum(hooks[child] for child in children) + 1
                
            queue.extend(list(self.upper_covers(elmt)))
            print(queue)
            
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

    def linear_extensions(self, facade):
        from .linear_extensions import LinearExtensionsOfPosetWithHooks
        return LinearExtensionsOfPosetWithHooks(self, facade=facade)