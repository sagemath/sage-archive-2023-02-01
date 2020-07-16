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
from .linear_extensions import LinearExtensionsOfPosetWithHooks
from .lattices import FiniteJoinSemilattice
from collections import deque

class DCompletePoset(FiniteJoinSemilattice):
    r"""
    D-complete posets are a class of posets introduced by Proctor in [Proc1999].
    It includes common families such as shapes, shifted shapes, and rooted forests. Proctor showed in [PDynk1999]
    that d-complete posets have decompositions in ``irreducible`` posets, and showed in [Proc2014] that
    d-complete posets admit a hook-length formula (see [1]). A complete proof of the hook-length formula
    can be found in [KY2019].
  
    EXAMPLES::

        sage: from sage.combinat.posets.poset_examples import Posets
        sage: P = Posets.DoubleTailedDiamond(2)
        sage: type(P)
        <class 'sage.combinat.posets.d_complete.DCompletePoset_with_category'>

    The additional internal data structure consists of:

    - the hook lengths of the elements of the poset

        sage: P._hooks
        {1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 5}

    TESTS::

        sage: TestSuite(P).run()

    See also the other tests in the class documentation.
    
    REFERENCES::
    
        .. [1] :wikipedia:`Hook_length_formula`, accessed 27th
                 June 2020.
    """
    
    _lin_ext_type = LinearExtensionsOfPosetWithHooks
    _desc = "Finite d-complete poset"

    @lazy_attribute
    def _hooks(self):
        r"""
        The hook lengths of the elements of the d-complete poset. For the definition 
        of hook lengths for d-complete posets, see [KY2019].

        TESTS::

            sage: from sage.combinat.posets.d_complete import DCompletePoset
            sage: P = DCompletePoset(DiGraph({0: [1, 2], 1: [3], 2: [3], 3: []}))
            sage: P._hooks
            {0: 1, 1: 2, 2: 2, 3: 3}
            sage: from sage.combinat.posets.poset_examples import Posets
            sage: P = DCompletePoset(Posets.YoungDiagramPoset(Partition([3,2,1]))._hasse_diagram.reverse()) 
            sage: P._hooks
            {0: 5, 1: 3, 2: 1, 3: 3, 4: 1, 5: 1}
        """
        hooks = {}

        min_diamond = {} # Maps max of double-tailed diamond to min of double-tailed diamond
        max_diamond = {} # Maps min of double-tailed diamond to max of double-tailed diamond
        
        H = self._hasse_diagram

        diamonds, _ = H.diamonds() # Tuples of four elements that are diamonds
        
        diamond_index = {} # Map max elmt of double tailed diamond to index of diamond

        # Find all the double-tailed diamonds and map the mins and maxes. 
        for index, d in enumerate(diamonds):
            min_diamond[d[3]] = d[0]
            max_diamond[d[0]] = d[3]
            diamond_index[d[3]] = index

            min_elmt = d[0]
            max_elmt = d[3]

            while True:
                potential_min = H.neighbors_in(min_elmt)
                
                potential_max = H.neighbors_out(max_elmt)

                # Check if any of these make a longer double tailed diamond
                found_diamond = False
                for (mn, mx) in [(i,j) for i in potential_min for j in potential_max]:
                    if len(H.neighbors_in(mx)) != 1:
                        continue
                    if len(H.all_paths(mn, mx)) == 2:
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
        queue = deque(H.sources())
        enqueued = set()
        while queue:
            elmt = queue.popleft()
            if elmt not in diamond_index:
                hooks[elmt] = H.order_ideal_cardinality([elmt])
            else:
                diamond = diamonds[diamond_index[elmt]]
                side1 = diamond[1]
                side2 = diamond[2]
                hooks[elmt] = hooks[side1] + hooks[side2] - hooks[min_diamond[elmt]]
            enqueued.add(elmt)

            for c in H.neighbors_out(elmt):
                if c not in enqueued:
                    queue.append(c)
                    enqueued.add(c)

        poset_hooks = {self._vertex_to_element(key): value for (key, value) in hooks.items()}
        return poset_hooks
                
    def get_hook(self, elmt):
        r"""
        Get the hook length of a specific element.

        TESTS::

            sage: from sage.combinat.posets.d_complete import DCompletePoset
            sage: P = DCompletePoset(DiGraph({0: [1], 1: [2]}))
            sage: P.get_hook(1)
            2
        """
        return self._hooks[elmt]

    def get_hooks(self):
        """
        Get all the hook lengths returned in a dictionary
        
        TESTS::

            sage: from sage.combinat.posets.d_complete import DCompletePoset
            sage: P = DCompletePoset(DiGraph({0: [1, 2], 1: [3], 2: [3], 3: []}))
            sage: P.get_hooks()
            {0: 1, 1: 2, 2: 2, 3: 3}
            sage: from sage.combinat.posets.poset_examples import Posets
            sage: P = DCompletePoset(Posets.YoungDiagramPoset(Partition([3,2,1]))._hasse_diagram.reverse()) 
            sage: P.get_hooks()
            {0: 5, 1: 3, 2: 1, 3: 3, 4: 1, 5: 1}

        """
        return dict(self._hooks)
    