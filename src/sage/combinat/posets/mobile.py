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
from sage.combinat.posets.d_complete import DCompletePoset
from sage.misc.lazy_attribute import lazy_attribute

#from .linear_extensions import LinearExtensionsOfForest

class MobilePoset(FinitePoset):
    r"""
    Mobile posets are an extension of d-complete posets which permit a determinant
    formula for counting linear extensions. 
    """
    
    #_lin_ext_type = LinearExtensionsOfForest
    _desc = 'Finite mobile poset'
    
    
    def __init__(self, hasse_diagram, elements, category, facade, key, ribbon=None):
        FinitePoset.__init__(self, hasse_diagram=hasse_diagram, elements=elements, category=category, facade=facade, key=key)
        if ribbon and self._test_valid_ribbon(ribbon):
            self._ribbon = ribbon
    
    def _test_valid_ribbon(self, ribbon):
        r"""
        Returns True if a ribbon has at most one anchor and that no vertex has two or more anchors
        
        INPUT:

            - ``ribbon`` -- a list of elements that form a ribbon in your poset   
        """
        ribbon = list(map(lambda x: self._element_to_vertex(x), ribbon))
        H = self._hasse_diagram
        R = H.subgraph(ribbon)
        num_anchors = 0
        for r in ribbon:
            anchor_neighbors = set(H.neighbors_out(r)).difference(set(R.neighbors_out(r)))
            if len(anchor_neighbors) == 1:
                num_anchors += 1
            elif len(anchor_neighbors) > 1:
                return False
        
        return num_anchors <= 1
        
    
    @lazy_attribute
    def _anchor(self):
        r"""
        The anchor of the mobile poset.
        
        TESTS::
        
            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M = MobilePoset(DiGraph([[0,1,2,3,4,5,6,7,8], [(1,0),(3,0),(2,1),(2,3),(4,3), (5,4),(5,6),(7,4),(7,8)]]))
            sage: M._anchor
            (4, 3)
            sage: M2 = MobilePoset(Poset([[0,1,2,3,4,5,6,7,8], [(1,0),(3,0),(2,1),(2,3),(4,3),(5,4),(7,4),(7,8)]]))   
            sage: M2._anchor
            (4, 3)
            sage: M3 = MobilePoset(Posets.RibbonPoset(5, [1,2]))
            sage: M3._anchor is None
            True
        """
        ribbon = list(map(lambda x: self._element_to_vertex(x), self._ribbon))
        H = self._hasse_diagram
        R = H.subgraph(ribbon)
        
        anchor = None
        
        # Find the anchor vertex, if it exists, and return the edge
        for r in ribbon:
            anchor_neighbors = set(H.neighbors_out(r)).difference(set(R.neighbors_out(r)))
            if len(anchor_neighbors) == 1:
                anchor = (r, anchor_neighbors.pop())
                break
        return (self._vertex_to_element(anchor[0]), self._vertex_to_element(anchor[1])) if not anchor is None else None
        
        
        
    @lazy_attribute
    def _ribbon(self):
        r"""
        The ribbon of the mobile poset.
        
        TESTS::
        
            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M = MobilePoset(DiGraph([[0,1,2,3,4,5,6,7,8], [(1,0),(3,0),(2,1),(2,3),(4,3), (5,4),(5,6),(7,4),(7,8)]]))
            sage: sorted(M._ribbon)
            [4, 5, 6, 7, 8]
            sage: M._test_valid_ribbon(M._ribbon)
            True
            sage: M2 = MobilePoset(Poset([[0,1,2,3,4,5,6,7,8], [(1,0),(3,0),(2,1),(2,3),(4,3),(5,4),(7,4),(7,8)]]))   
            sage: sorted(M2._ribbon)
            [4, 7, 8]
            sage: M2._test_valid_ribbon(M2._ribbon)
            True
            sage: M3 = MobilePoset(Posets.RibbonPoset(5, [1,2]))
            sage: sorted(M3._ribbon)
            [1, 2, 3, 4]
            sage: M3._test_valid_ribbon(M3._ribbon)
            True
        """
        return list(map(lambda x: self._vertex_to_element(x), self._compute_ribbon()))
    
    def _compute_ribbon(self):
        r"""
        The helper function of _ribbon for computing the ribbon.
        """
        H = self._hasse_diagram
        H_un = H.to_undirected()
        max_elmts = H.sinks()
        
        # Compute anchor, ribbon
        ribbon = [] # In order list of elements on zigzag
        
        if len(max_elmts) == 1:
            return [max_elmts[0]]
        # Compute max element tree by merging shortest paths
        start = max_elmts[0]
        
        zigzag_elmts = set()
        for m in max_elmts[1:]:
            sp = H_un.shortest_path(start, m)
            zigzag_elmts.update(sp)
        
        max_elmt_graph = H.subgraph(zigzag_elmts)
        G = max_elmt_graph.to_undirected()
        
        if G.is_path():
            # Check if there is a anchor by seeing if there is more than one acyclic path to the next max
            ends = max_elmt_graph.vertices_with_degree(1)
            # Form ribbon
            ribbon = G.shortest_path(ends[0], ends[1])
            for end_count, end in enumerate(ends):
                if not (H_un.is_cut_vertex(end) or H_un.degree(end) == 1):
                    traverse_ribbon = ribbon if end_count == 0 else ribbon[::-1]
                    for ind, p in enumerate(traverse_ribbon):
                        if H_un.is_cut_edge(p, traverse_ribbon[ind+1]):
                            return G.shortest_path(ends[(end_count + 1) % 2], traverse_ribbon[ind+1]) 
            return ribbon
                            
            
        else:
            # First check path counts between ends and deg3 vertex
            # Then check if more than one max elmt on way to degree 3 vertex. 
            # Arbitrarily choose between ones with just 1
            
            ends = max_elmt_graph.vertices_with_degree(1)
            deg3 = max_elmt_graph.vertices_with_degree(3)[0]
            
            anchoredEnd = None
            for end in ends:
                if not (H_un.is_cut_vertex(end) or H_un.degree(end) == 1):
                    anchoredEnd = end
                    break
            
            if not anchoredEnd is None:
                path = H.shortest_path(deg3, anchoredEnd)
                ends.remove(anchoredEnd)
                
                return G.shortest_path(ends[0], ends[1])
                
            else:
                possible_anchors = ends[:]
                for end in ends:
                    path = G.shortest_path(end, deg3)
                    if not sum(map(lambda z: z in max_elmts, path)) == 1:
                        possible_anchors.remove(end)
                
                anchoredEnd = possible_anchors[0]
                ends.remove(anchoredEnd)
                return  G.shortest_path(ends[0], ends[1])
        
    def get_ribbon(self):
        r"""
        Returns the ribbon of the mobile poset
        """
        return self._ribbon
    
    def get_anchor(self):
        r"""
        Returns the anchor of the mobile poset
        """
        return self._anchor
                
                
            
                

                
            
            
        
        
            
        
        
            
            
            
            
            
            
            
            

        
        
            
            
    
    
    
