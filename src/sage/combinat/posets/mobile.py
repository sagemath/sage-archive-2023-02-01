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
#from .linear_extensions import LinearExtensionsOfForest

class MobilePoset(FinitePoset):
    r"""
    Mobile posets are an extension of d-complete posets which permit a determinant
    formula for counting linear extensions. 

    """
    
    #_lin_ext_type = LinearExtensionsOfForest
    _desc = 'Finite mobile poset'
    
    def _mobile_structure(self):
        H = self._hasse_diagram
        H_un = H.to_undirected()
        max_elmts = H.sinks()
        
        # Compute plant, ribbon
        ribbon = [] # In order list of elements on zigzag
        
        plant = (0, 0) # The cut edge separating the plant from the zigzag

        # Compute max element tree by merging shortest paths
        
        start = max_elmts[0]
        
        zigzag_elmts = set()
        for m in max_elmts[1:]:
            sp = G.shortest_path(start, m)
            zigzag_elmts.add(sp)
        
        max_elmt_graph = H.subgraph(zigzag_elmts)
        G = max_elmt_graph.to_undirected()
        
        if G.is_path():
            # Check if there is a plant by seeing if there is more than one acyclic path to the next max
            ends = max_elmt_graph.vertices_with_degree(1)
            
            # Form ribbon
            ribbon = G.shortest_path(ends[0], ends[1])
            
            for end in ends:
                path = []
                nextElmt = end
                # Get next maximal element in zigzag
                while True:
                    path.append(nextElmt)
                    nbrs = G.neighbors(nextElmt)
                    nextElmt = nbrs[0] if nbrs[0] != nextElmt else nbrs[1]
                    if nextElmt in max_elmts:
                        break
                
                # Find spot where we enter zigzag
                foundPlant = False
                
                for i, v in enumerate(path):
                    if len(H_un.all_paths(v, nextElmt)) > 1:
                        foundPlant = True
                        continue
                    elif i == 0:
                        break
                    if foundPlant:
                        plant = (v, path[i-1])
                        # Shorten zigzag
                        
                        endIndex = ribbon.index(end)
                        plantIndex = ribbon.index(v)
                        
                        if plantIndex < endIndex:
                            ribbon = ribbon[:plantIndex]
                        else:
                            ribbon = ribbon[plantIndex:]
                            
                        break
                    
                if foundPlant:
                    break
            
        else:
            # First check path counts between ends and deg3 vertex
            # Then check if more than one max elmt on way to degree 3 vertex. 
            # Arbitrarily choose between ones with just 1
            
            ends = max_elmt_graph.vertices_with_degree(1)
            deg3 = max_elmt_graph.vertices_with_degree(3)[0]
            
            plantedEnd = None
            for end in ends:
                if H_un.all_paths(end, deg3) > 1:
                    plantedEnd = end
                    break
            
            if not plantedEnd is None:
                path = H.shortest_path(deg3, plantedEnd)
                plant = (path[0], path[1])
                
                ends.remove(plantedEnd)
                ribbon = max_elmt_graph.shortest_path(end[0], end[1])
                
            else:
                possible_plants = ends[:]
                for end in ends:
                    path = G.shortest_path(end, deg3)
                    if not reduce(lambda x,y: x^y, map(lambda z: z in max_elmts, path)):
                        possible_plants.remove(end)
                
                plantedEnd = possible_plants[0]
                ends.remove(plantedEnd)
                ribbon = G.shortest_path(ends[0], ends[1])
                
                plant = (deg3, G.shortest_path(deg3, plantedEnd)[1])
            
        self._ribbon = ribbon
        self._plant = plant
                
                
            
                

                
            
            
        
        
            
        
        
            
            
            
            
            
            
            
            

        
        
            
            
    
    
    
