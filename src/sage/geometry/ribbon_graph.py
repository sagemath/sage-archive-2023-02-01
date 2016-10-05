r"""
Ribbon Graphs

This file implements the object 'ribbon graph'.These are graphs 
together with a cyclic ordering of the half edges adjacent to each 
vertex. This data allows us to unambiguosly "thicken" the ribbon 
graph to an orientable surface with boundary. 

AUTHORS:

- Pablo Portilla (2016): initial version

"""

#*****************************************************************************
#       Copyright (C) 2016 Pablo Portilla  <p.portilla89@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.groups.perm_gps.permgroup_element import *
import copy

#Auxiliary functions that will be used in the classes.

def find(l, k):
    r"""
    Return the two coordinates of the element k in the list of lists l.

    INPUT:

    - ''l'' a list of lists.

    - ``k`` a candidate to be in a list in l.

    OUTPUT:

    A list with two integers describing the position of the first
    instance of k in l.

    EXAMPLES::

        sage: from sage.combinat.ribbon_graph import find
        sage: A = [[2,3,4],[4,5,2],[8,7]]
        sage: find(A,2)
        [0, 0]
        sage: find(A,7)
        [2, 1]
        sage: find(A,5)
        [1, 1]
    """
    pos=[]
    found = False
    i=0
    while i<len(l) and found!=True:
        j=0
        while j<len(l[i]) and found!=True:
            if (k == l[i][j]):
                pos.append(i)
                pos.append(j)
                found = True
            else:
                j+=1
        i+=1
    return pos



class RibbonGraph(SageObject):
    r"""
    A ribbon graph as two elements of a symmetric group.

    INPUT:

    - ''sigma'' a permutation in G which is a product of disjoint

    - ''rho'' a permutation in G which is a product of disjoint
     2-cycles and singletons (1-cycles)

    -''sing'' a boolean variable that tells if we consider vertices of
     valency 1 or not in our ribbon graph

    EXAMPLES::

        sage: s0 = PermutationGroupElement('(1)(2)')
        sage: r0 = PermutationGroupElement('(1,2)')
        sage: R0 = RibbonGraph(s0,r0); R0
        Sigma: [[1], [2]] 
        Rho: [[1, 2]]

        sage: s1 = PermutationGroupElement('(1,3,5,8)(2,4,6)')
        sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)(8,15)')
        sage: R1 = RibbonGraph(s1,r1); R1
        Sigma: [[1, 3, 5, 8], [2, 4, 6], [15]] 
        Rho: [[1, 2], [3, 4], [5, 6], [8, 15]]

        sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
        sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)')
        sage: R2 = RibbonGraph(s2,r2); R2
        Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18]] 
        Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12]]
    """

    def __init__(self, sigma, rho):
        """
        Initialize ``self``.

        """
        #with the 'trick' (sigma*rho)*rho**(-1) we get that if rho was
        #iniated with PermutationGroupElement and was naturally in a 
        #bigger permutation group than sigma, then automatically 
        #self.sigma lies in this bigger group. This avoids some
        #patologies. Note that the reverse thing cannot happen since
        #in the permutation rho all darts appear always.
        self.sigma=(sigma*rho)*rho**(-1)
        self.rho=rho
        #have to add lines to control error such as not admissible sigma, rho
        #things that are not ribbon graphs, etc.

    def _repr_(self):
        r'''
        Return string representation of the two permutations that define
        the ribbon graph.

        '''

        #On the first lines, we compute the vertices of valency 1 to add
        #them to the list repr_sigma.
        
        repr_sigma = [list(x) for 
                     x in self.sigma.cycle_tuples()]
        repr_rho = [list(x) for 
                   x in self.rho.cycle_tuples()]
        darts_rho = [j for i in range(len(repr_rho)) 
                        for j in repr_rho[i]]
        darts_sigma = [j for i in range(len(repr_sigma)) 
                        for j in repr_sigma[i]]
        val_one = [x for x in darts_rho if x not in darts_sigma]
        for i in range(len(val_one)):
            repr_sigma += [[val_one[i]]]
            
        return ('Sigma: %r \n' \
                'Rho: %r'% (repr_sigma, 
                                repr_rho
                            )
                )


    def n_boundary(self):
        r'''
        Return number of boundary components of the thickening of the
        ribbon graph.

        EXAMPLES:

        The first example is the ribbon graph corresponding to the torus
        with one hole::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1)
            sage: R1.n_boundary()
            1

        This example is constructed by taking the bipartite graph of 
        type '(3,3)'::

            sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
            sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)')
            sage: R2 = RibbonGraph(s2,r2)
            sage: R2.n_boundary()
            3
        '''
        return len((self.rho*self.sigma).cycle_tuples())


    def contract_edge(self, k):
        r"""
        Return the ribbon graph resulting from the contraction of
        the edge k in the ribbon graph given by sigma, rho. 
        That is, contracts the edge corresponding to the k-th 
        transposition of rho.

        INPUT:

        - ''k''-- non-negative integer. The position in rho of the
        transposition that is going to be contracted.

        OUTPUT:

        -A ribbon graph resulting from the contraction of that edge.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1);R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: S1=R1.contract_edge(1); S1
            Sigma: [[1, 6, 2, 5]] 
            Rho: [[1, 2], [5, 6]]

        However, this ribbon graphs is formed only by loops and hence
        it cannot be longer reduced, we get an error if we try to
        contract a loop::

            sage: S1.contract_edge(1)
            Traceback (most recent call last):
            ...
            ValueError: The edge is a loop and cannot be contracted

        In this example, we consider a graph that has one edge '(19,20)'
        such that one of its ends is a vertex of valency '1'. This is 
        the vertex '(20)' that is not specified when defining sigma.
        We contract precisely this edge and get a ribbon graph with no
        vertices of valency 1::

            sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18,19)')
            sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)(19,20)')
            sage: R2 = RibbonGraph(s2,r2);R2
            Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18, 19], [20]] 
            Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12], [19, 20]]
            sage: R2.contract_edge(9)
            Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18]] 
            Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12]]
            """

        #the following two lines convert the list of tuples to list of lists
        aux_sigma = [list(x) for 
                     x in self.sigma.cycle_tuples(singletons = 1)]
        aux_rho = [list(x) for 
                   x in self.rho.cycle_tuples()]
        #The following if rules out the cases when we would be 
        #contracting a loop (which is not admissible since we would 
        #lose the topological type of the graph).
        if (find(aux_sigma, aux_rho[k][0])[0] == 
                find(aux_sigma, aux_rho[k][1])[0]):
            raise ValueError("The edge is a loop and " \
                              "cannot be contracted")
        #We store in auxiliary variables the positions of the vertices
        #that are the ends of the edge to be contracted and we delete
        #from them the darts corresponding to the edge that is going
        #to be contracted. We also delete the contracted edge 
        #from aux_rho
        pos1 = find(aux_sigma, aux_rho[k][0])
        pos2 = find(aux_sigma, aux_rho[k][1])
        del aux_sigma[pos1[0]][pos1[1]]
        del aux_sigma[pos2[0]][pos2[1]]
        del aux_rho[k]

        #Now we insert in one of the two vertices, the darts of the other
        #vertex that appears after. We make sure that we don't
        #change the topological type of the thickening of the graph by
        #preserving the cyclic ordering.
        n = len(aux_sigma[pos2[0]])
        for i in range(n):
            aux_sigma[pos1[0]].insert(
                pos1[1]+i,
                aux_sigma[pos2[0]][(pos2[1]+i) % n]
            )
        #Finally we delete the vertex from which we copied all the darts.
        del aux_sigma[pos2[0]]

        #Now we convert this data that is on the form of lists of lists
        #to actual permutations that form a ribbon graph.
        return RibbonGraph(
                           PermutationGroupElement([tuple(x) for x in aux_sigma]), 
                           PermutationGroupElement([tuple(x) for x in aux_rho])
                           )

        

    def genus(self):
        r'''
        Return the genus of the thickening of the ribbon graph.
        
        '''
        
        vertices = len(self.sigma.cycle_tuples(singletons = 1))
        edges = len(self.rho.cycle_tuples())
        
        #formula for the genus using that the thickening is homotopically 
        #equivalent to the graph
        g = (vertices - edges + self.n_boundary()-2)/(-2)
        
        return g
    
    def mu(self):
        r'''
        Return the rank of the first homology group of the thickening
        of the ribbon graph.

        '''
        return 2*self.genus() + self.n_boundary() - 1


    def boundary(self):
        r'''
        Return the labeled boundaries of the ribbon graph.

        '''

        #initialize and empty list to hold the labels of the boundaries
        bound = []
        
        #since lists of tuples are not modifiable, we change the data to a
        #list of lists 
        aux_perm = (self.rho*self.sigma).cycle_tuples()
        
        #the cycles of the permutation rho*sigma are in 1:1 correspondence with 
        #the boundary components of the thickening (see function n_boundary())
        #but they are not the labeled boundary components.
        #With the next for, we convert the cycles of rho*sigma to actually 
        #the labelling of the edges. Each edge, therefore, should appear twice
        
        for i in range(len(aux_perm)):
            bound = bound + [[]]
            for j in range(len(aux_perm[i])):
                bound[i].append(aux_perm[i][j])
                bound[i].append(self.rho(aux_perm[i][j]))

        #finally the function returns a List of lists. Each list contains
        #a sequence of  numbers and each number corresponds to a half-edge.
        return bound

    #The following function takes the ribbon graph and returns another ribbon 
    #graph where the permutation sigma has just one cycle, (so the graph has 
    #1 vertex) and the permutation rho is therefore formed by mu 2-cycles,
    #where mu is the first betti number of the surface


    def reduced(self):
        r'''
        Return a ribbon graph with 1 vertex and mu edges (where $\mu$ is
        the first betti number of the graph).

        OUTPUT:
        
        -A ribbon graph whose sigma permutation has only 1 non-trivial
        cycle and whose rho permutation is a product of $\mu$ disjoint
        2-cycles.
        '''
        
        #the following two lines convert the list of tuples to list of lists
        #we have to contract exactly n edges
        aux_ribbon = copy.deepcopy(self)
        rho_self = [list(x) for 
                    x in aux_ribbon.rho.cycle_tuples()]
        n = len(rho_self)- self.mu()
        for i in range(n):
            aux_sigma = [list(x) for 
                     x in aux_ribbon.sigma.cycle_tuples(singletons = 1)]
            aux_rho = [list(x) for 
                    x in aux_ribbon.rho.cycle_tuples()]
            for j in range(len(aux_rho)):
                if (find(aux_sigma, aux_rho[j][0])[0] != 
                        find(aux_sigma, aux_rho[j][1])[0]):
                    aux_ribbon = aux_ribbon.contract_edge(j)
                    break
                else:
                    continue
        #finally we change the data to a list of tuples and return the
        #information as a ribbon graph. 
        return aux_ribbon
    
    #the next function computes a basis of homology, it uses
    #the previous function.
    
    def homology_basis(self):
        r'''
        Return an oriented basis of the firs homology group of the graph.
        
        OUTPUT:
        
        -A LIST of Lists of lists. Each List corresponds to an element 
        of the basis and each list in a List is just a 2-tuple which 
        corresponds to an 'ordered' edge of rho.
        '''
        
        aux_sigma = [list(x) for 
                     x in self.sigma.cycle_tuples(singletons = 1)]
        
        basis = [[list(x)] for 
                 x in self.reduced().rho.cycle_tuples()]
        
        #Now we define center as the set of edges that were contracted 
        #in reduced() this set is contractible and can be define as the 
        #complement of reduced_rho in rho
        
        center = [list(x) for x 
                  in self.rho.cycle_tuples(singletons = 1) 
                  if (x not in self.reduced().rho.cycle_tuples())]
        
        #We define an auxiliary list 'vertices' that will contain the
        #vertices (cycles of sigma) corresponding to each half edge. 
        
        vertices=[]
        
        for i in range(len(basis)):
            
            vertices = vertices + [[]]
            basis[i].extend(copy.deepcopy(center))
            
            for j in range (len(basis[i])):
                
                vertices[i].append(find(aux_sigma,basis[i][j][0])[0])
                vertices[i].append(find(aux_sigma,basis[i][j][1])[0])
            k = 0
            
            while k < (len(vertices[i])):
                
                if (vertices[i].count(vertices[i][k]) == 1):
                    
                    m = int(floor(k/2))
                    del basis[i][m]
                    del vertices[i][2*m:2*m+2]
                    k = 0
                    
                else:
                    k+=1

        for i in range (len(basis)):
            
            for j in range (1,len(basis[i])):
                
                n=[t for t, n in list(enumerate(vertices[i])) 
                   if n == vertices[i][2*j-1]][1]
                
                ind = int(floor(n/2))
                
                if (j != ind):
                    
                    basis[i][j], basis[i][ind] = basis[i][ind], basis[i][j]
                    
                    vertices[i][2*j], vertices[i][2*ind] = \
                    vertices[i][2*ind], vertices[i][2*j]
                
                    vertices[i][2*j+1], vertices[i][2*ind+1] = \
                    vertices[i][2*ind+1], vertices[i][2*j+1]

                if (vertices[i][2*j-1] != vertices[i][2*j]):
                    
                    vertices[i][2*j], vertices[i][2*j+1] = \
                    vertices[i][2*j+1], vertices[i][2*j]
                
                    basis[i][j][0], basis[i][j][1] = \
                    basis[i][j][1], basis[i][j][0]
                    
        #the variable basis is a LIST of Lists of lists. Each List 
        #corresponds to an element of the basis and each list in a List
        #is just a 2-tuple which corresponds to an 'ordered' edge of rho.
        
        return basis


