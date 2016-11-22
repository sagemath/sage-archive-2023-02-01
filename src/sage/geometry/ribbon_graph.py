r"""
Ribbon Graphs

This file implements the object ''ribbon graph''.These are graphs 
together with a cyclic ordering of the darts adjacent to each 
vertex. This data allows us to unambiguosly ''thicken'' the ribbon 
graph to an orientable surface with boundary. Also, every orientable
surface with non-empty boundary is the thickening of a ribbon graph.

A comprenhensive introduction on the topic can be found in the beginning
of [GIR]_ Chapter 4. More concretely, we will use a variation of what
is called in the reference ''The permutation representation pair of a
dessin''. Note that in that book, ribbon graphs are called ''dessins
d'enfant''. For the sake on completeness we reproduce an adapted version
of that introduction here.

**Brief introduction**

Let $\Sigma$ be an orientable surface with non-empty boundary and let 
$\Gamma$ be the topological realization of a graph that is embedded in 
$\Sigma$ in such a way that the graph is a strong deformation retract of
the surface. 

Let $v(\Gamma)$ be the set of vertices of $\Gamma$, suppose that these
are ''white'' vertices. Now we mark a ''black'' vertices in an interior 
point of each edge. In this way we get a bipartite graph where all the
black vertices have valency 2 and there is no restriction on the valency
of the white vertices. We call the edges of this new graph ''darts''
(sometimes they are also called ''half eldges''  of the original graph).
Observe that each edge of the original graph is formed by two darts.

Given a white vertex $v \in v(\Gamma)$, let $d(v)$ be the set of darts adjacent
to $v$. Let $D(\Gamma)$ be the set of all the darts of
$\Gamma$ and suppose that we enumerate the set $D(\Gamma)$ and that it
has $n$ elements.

With the orientation of the surface and the embedding of the graph in 
the surface we can produce two permutations:

    - A permutation that we will call ''sigma'' and denote by `\sigma`. 
      This permutation is a product of as many cycles as white vertices 
      (that is vertices in `\Gamma`). For each vertex consider a small 
      topological circle around it in `\Sigma`. This circle intersects 
      each adjacent dart once. The circle has an orientation induced by
      the orientation on `\Sigma` and so defines a cycle that sends the
      number associated to one dart to the number associated to the next
      dart in the positive orientation of the circle.

    - A permutation that we will call ''rho'' and denote by `\rho`. 
      This permutation is a product of as many $2$-cycles as edges has 
      `\Gamma`. It just tells which two darts belong to the same edge.

One can also define a ribbon graph abstractly:

    Consider a graph $\Gamma$ (not a priori embedded in any surface). 
    Now we can again consider one vertex in the interior of each edge 
    splitting each edge in two darts. We label the darts with numbers.

    We say that a ribbon structure on $\Gamma$ is a set of two 
    permutations `(\sigma, \rho)`. Where `\sigma` is formed by as many
    disjoint cycles as vertices had `\Gamma`. And each cycle is a 
    cyclic ordering of the darts adjacent to a vertex. The permutation
    `\rho` just tell us which two darts belong to the same edge.

    For any two such permutations there is a way of ''thickening'' the
    graph to a surface with boundary in such a way that the surface
    retracts (by a strong deformation retract) to the graph and hence
    the graph is embedded in the surface in a such a way that we could
    recover `\sigma` and `\rho`.


EXAMPLES:

Consider a graph that has $2$ vertices of valency $3$ (and hence $3$
edges). That is represented by  the  following two permutations::

    sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
    sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
    sage: R1 = RibbonGraph(s1,r1);R1
    Sigma: [[1, 3, 5], [2, 4, 6]] 
    Rho: [[1, 2], [3, 4], [5, 6]]

By drawing the picture in a piece of paper as explained in the 
introduction, one can see that its thickening has only $1$ boundary 
component. Since the the thickening is homotopically equivalent to the
graph and the graph has euler characteristic $-1$, we find that the
thickening has genus $1$::

    sage: R1.n_boundary()
    1
    sage: R1.genus()
    1

The next example is the complete bipartite graph of type (3,3) where we
have added an edge that ends at a vertex of valency 1::

    sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18,19)')
    sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)(19,20)')
    sage: R2 = RibbonGraph(s2,r2);R2
    Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18, 19], [20]] 
    Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12], [19, 20]]
    sage: R2.contract_edge(9)
    Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18]] 
    Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12]]
    sage: S2=R2.reduced();S2
    Sigma: [[5, 6, 8, 9, 14, 15, 11, 12]] 
    Rho: [[5, 14], [6, 11], [8, 15], [9, 12]]
    sage: R2.genus(); S2.genus()
    1
    1
    sage: R2.boundary()
    [[1, 16, 17, 4, 5, 14, 15, 8, 9, 12, 10, 3],
    [2, 13, 14, 5, 6, 11, 12, 9, 7, 18, 19, 20, 20, 19, 16, 1],
    [3, 10, 11, 6, 4, 17, 18, 7, 8, 15, 13, 2]]
    sage: S2.boundary()
    [[5, 14, 15, 8, 9, 12], [6, 11, 12, 9, 14, 5], [8, 15, 11, 6]]
    sage: R2.homology_basis()
    [[[5, 14], [13, 2], [1, 16], [17, 4]],
    [[6, 11], [10, 3], [1, 16], [17, 4]],
    [[8, 15], [13, 2], [1, 16], [18, 7]],
    [[9, 12], [10, 3], [1, 16], [18, 7]]]
    sage: S2.homology_basis()
    [[[5, 14]], [[6, 11]], [[8, 15]], [[9, 12]]]

AUTHORS:

- Pablo Portilla (2016)

REFERENCES:

.. [GIR] E. Girondo, G. Gonzalez-Diez, Introduction to Compact 
   Riemann surfaces and Dessins d'enfant, London Mathematical Society,
   SStudent Text 79.

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
from sage.functions.other import floor
import copy

#Auxiliary functions that will be used in the classes.

def _find(l, k):
    r"""
    Return the two coordinates of the element k in the list of lists l.

    INPUT:

    - ``l`` a list of lists.

    - ``k`` a candidate to be in a list in l.

    OUTPUT:

    A list with two integers describing the position of the first
    instance of k in l.

    TESTS::

        sage: from sage.geometry.ribbon_graph import _find
        sage: A = [[2,3,4],[4,5,2],[8,7]]
        sage: _find(A,2)
        [0, 0]
        sage: _find(A,7)
        [2, 1]
        sage: _find(A,5)
        [1, 1]

    """
    for i,lst in enumerate(l):
        if k in lst:
            return [i, lst.index(k)]
    raise ValueError("element {} not found".format(k))

def _clean(l):
    r"""
    Return a list where empty sublists of l have been removed.

    INPUT:

    - ``l`` -- A List of lists.

    OUTPUT:

    - A list which is a copy of l with all empty sublists removed.

    EXAMPLES::

        sage: from sage.geometry.ribbon_graph import _clean
        sage: A = [[1,2],[], [2,1,7],[],[],[1]]
        sage: _clean(A)
        [[1, 2], [2, 1, 7], [1]]

    """
    return [copy(elt) for elt in l if len(elt)>0]

class RibbonGraph(SageObject):
    r"""
    A ribbon graph codified as two elements of a certaing permutation 
    group.

    INPUT:

    - ``sigma`` -- a permutation a product of disjoint cycles of any
      length. Singletons (vertices of valency 1) need not be specified.

    - ``rho`` -- a permutation which is a product of disjoint
      2-cycles.


    EXAMPLES:

    Consider the ribbon graph consisting of just $1$ edge and $2$
    vertices of valency $1$::

        sage: s0 = PermutationGroupElement('(1)(2)')
        sage: r0 = PermutationGroupElement('(1,2)')
        sage: R0 = RibbonGraph(s0,r0); R0
        Sigma: [[1], [2]] 
        Rho: [[1, 2]]

    The following example corresponds to the  complete bipartite graph
    of type $(2,3)$, where we have added one more edge $(8,15)$ that
    ends at a vertex of valency $1$. Observe that it is not necessary
    to specify the  vertex $(15)$ of valency $1$ when we define sigma::

        sage: s1 = PermutationGroupElement('(1,3,5,8)(2,4,6)')
        sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)(8,15)')
        sage: R1 = RibbonGraph(s1,r1); R1
        Sigma: [[1, 3, 5, 8], [2, 4, 6], [15]] 
        Rho: [[1, 2], [3, 4], [5, 6], [8, 15]]

    This example is constructed by taking the bipartite graph of 
    type '(3,3)'::

        sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
        sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)')
        sage: R2 = RibbonGraph(s2,r2); R2
        Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18]] 
        Rho: [[1, 16], [2, 13], [3, 10], [4, 17], [5, 14], [6, 11], [7, 18], [8, 15], [9, 12]]

    The labeling of the darts can omit some numbers::

        sage: s3 = PermutationGroupElement('(3,5,10,12)')
        sage: r3 = PermutationGroupElement('(3,10)(5,12)')
        sage: R3= RibbonGraph(s3,r3);R3
        Sigma: [[3, 5, 10, 12]] 
        Rho: [[3, 10], [5, 12]]

    """

    def __init__(self, sigma, rho):
        r"""
        Initialize ``self``.

        """
        #with the 'trick' (sigma*rho)*rho**(-1) we get that if rho was
        #iniated with PermutationGroupElement and was naturally in a 
        #bigger permutation group than sigma, then automatically 
        #self.sigma lies in this bigger group. This avoids some
        #patologies. Note that the reverse thing cannot happen since
        #in the permutation rho all darts of the graph appear.
        self.sigma=(sigma*rho)*rho**(-1)
        self.rho=rho


    def _repr_(self):
        r"""
        Return string representation of the two permutations that define
        the ribbon graph.

        """

        #On the first lines, we compute the vertices of valency 1 to add
        #them to the list repr_sigma.This is was implemented in order
        #to avoid having the user to input more data. For example, if 
        #the order of the permutation group where sigma lies is bigger
        #than the number of darts, it could be interpreted that sigma
        #has more vertices of valency 1 that it really has. These lines
        #compute which of the singletons of sigma are actually vertices
        #of valency 1.
        repr_sigma = [list(x) for x in self.sigma.cycle_tuples()]
        repr_rho = [list(x) for x in self.rho.cycle_tuples()]
        darts_rho = flatten(repr_rho)
        darts_sigma = flatten(repr_sigma)
        val_one = [x for x in darts_rho if x not in darts_sigma]
        for i in range(len(val_one)):
            repr_sigma += [[val_one[i]]]

        return ('Sigma: %r \n' \
                'Rho: %r'% (repr_sigma, 
                            repr_rho
                            )
                )


    def n_boundary(self):
        r"""
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

        """
        #it might seem a bit overkill to call boundary() here but it is
        #necessary to either call it or do similar computations here.
        #The function boundary() avoids some patologies with boundaries
        #formed by just one loop.
        return len(self.boundary())


    def contract_edge(self, k):
        r"""
        Return the ribbon graph resulting from the contraction of
        the edge k in the ribbon graph given by sigma, rho. 
        That is, contracts the edge corresponding to the k-th 
        transposition of rho.

        INPUT:

        - ``k`` -- non-negative integer. The position in rho of the 
          transposition that is going to be contracted.

        OUTPUT:

        - A ribbon graph resulting from the contraction of that edge.

        EXAMPLES:

        We start again with the  one-holed torus ribbon graph::

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
        #The following ''if'' rules out the cases when we would be 
        #contracting a loop (which is not admissible since we would 
        #lose the topological type of the graph).
        if (_find(aux_sigma, aux_rho[k][0])[0] == 
                _find(aux_sigma, aux_rho[k][1])[0]):
            raise ValueError("The edge is a loop and " \
                              "cannot be contracted")
        #We store in auxiliary variables the positions of the vertices
        #that are the ends of the edge to be contracted and we delete
        #from them the darts corresponding to the edge that is going
        #to be contracted. We also delete the contracted edge 
        #from aux_rho
        pos1 = _find(aux_sigma, aux_rho[k][0])
        pos2 = _find(aux_sigma, aux_rho[k][1])
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

    def extrude_edge(self, vertex, dart1, dart2):
        r"""
        Return a ribbon graph resulting from extruding an edge from a
        vertex, pulling from it, all darts from dart1 to dart2 including
        both.

        INPUT:

        - ``vertex`` -- the position of the vertex in the permutation
          sigma when represented by the method _repr_. It has to be
          a vertex of valency at least 2.

        - ``dart1`` -- the position of the first in the
          cycle corresponding to ``vertex``.

        - ``dart2`` -- the position of the second dart inn the cycle
          corresponding to ``vertex``.

        OUTPUT:

        - A ribbon graph resulting from extruding a new edge that 
          ''pulls'' from ``vertex`` a new vertex that is, now, adjacent
          to all the darts from ``dart1``to ``dart2`` (not including
          dart2) in the cyclic ordering given by the cycle corresponding
          to ``vertex``. Note that ``dart1`` may be equal to ``dart2``
          allowing thus to extrude a contractible edge from a vertex.

        EXAMPLES:

        We try several possibilities in the same graph::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1);R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: R1.extrude_edge(1,1,2)
            Sigma: [[1, 3, 5], [2, 8, 6], [4, 7]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8]]
            sage: R1.extrude_edge(1,1,3)
            Sigma: [[1, 3, 5], [2, 8], [4, 6, 7]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8]]

        We can also extrude a contractible edge from a vertex. This
        new edge will end at a vertex of valency 1::

            sage: R1.extrude_edge(0,0,0)
            Sigma: [[1, 3, 5, 8], [2, 4, 6], [7]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8]]

        In the following example we first extrude one edge from a vertex
        of valency 3 generating a new vertex of valency 2. Then we 
        extrude a new edge from this vertex of valency 2::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1); R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: E1 = R1.extrude_edge(0,0,1); E1
            Sigma: [[1, 7], [2, 4, 6], [3, 5, 8]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8]]
            sage: F1 = E1.extrude_edge(0,0,1); F1 
            Sigma: [[1, 9], [2, 4, 6], [3, 5, 8], [7, 10]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]

        """

        #We first compute the vertices of valency 1 as in _repr_
        repr_sigma = [list(x) for x in self.sigma.cycle_tuples()]
        repr_rho = [list(x) for x in self.rho.cycle_tuples()]
        darts_rho = flatten(repr_rho)
        darts_sigma = flatten(repr_sigma)
        val_one = [x for x in darts_rho if x not in darts_sigma]
        for i in range(len(val_one)):
            repr_sigma += [[val_one[i]]]

        #we find which is the highes value a dart has, in order to 
        #add new darts that do not conflict with previous ones.
        k = max(darts_rho)

        #We create the new vertex and append it to sigma.
        new_vertex = [repr_sigma[vertex][j] for j in range(dart1,dart2)]
        new_vertex.insert(0,k+1)
        repr_sigma.append(new_vertex)

        #We add the  new dart at the vertex from which we are extruding
        #an edge. Also we delete the darts that have been extruded.
        repr_sigma[vertex].insert(dart1,k+2)
        del repr_sigma[vertex][dart1+1:dart2+1]

        #We update rho
        repr_rho.append([k+1,k+2])

        return RibbonGraph(
                           PermutationGroupElement([tuple(x) for x in repr_sigma]), 
                           PermutationGroupElement([tuple(x) for x in repr_rho])
                           )


    def genus(self):
        r"""
        Return the genus of the thickening of the ribbon graph.

        OUTPUT:

        - ``g`` -- non-negative integer representing the genus of the
          thickening of the ribbon graph.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1);R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: R1.genus()
            1

            sage: s3=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15,16)(17,18,19,20)(21,22,23,24)')
            sage: r3=PermutationGroupElement('(1,21)(2,17)(3,13)(4,22)(7,23)(5,18)(6,14)(8,19)(9,15)(10,24)(11,20)(12,16)')
            sage: R3 = RibbonGraph(s3,r3);R3.genus()
            3

        """
        #We now use the same procedure as in _repr_ to get the vertices
        #of valency 1 and distinguish them from the extra singletons of
        #the permutation sigma. 
        repr_sigma = [list(x) for x in self.sigma.cycle_tuples()]
        repr_rho = [list(x) for x in self.rho.cycle_tuples()]
        darts_rho = flatten(repr_rho)
        darts_sigma = flatten(repr_sigma)
        val_one = [x for x in darts_rho if x not in darts_sigma]

        #the total number of vertices of sigma is its number of cycles
        #of length >1 plus the number of singletons that are actually
        #vertices of valency 1
        
        vertices = len(self.sigma.cycle_tuples()) + len(val_one)
        edges = len(self.rho.cycle_tuples())
        #formula for the genus using that the thickening is homotopically 
        #equivalent to the graph
        g = (vertices - edges + self.n_boundary()-2)/(-2)

        return g

    def mu(self):
        r"""
        Return the rank of the first homology group of the thickening
        of the ribbon graph.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1);R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: R1.mu()
            2

        """
        return 2*self.genus() + self.n_boundary() - 1


    def boundary(self):
        r"""
        Return the labeled boundaries of the ribbon graph.

        If you cut the thickening of the graph along the graph. you
        get a collection of cylinders (recall that the graph was a
        strong deformation retract of the thickening). In each cylinder
        one of the boundary components has a labelling of its edges
        induced by the labelling of the darts.

        OUTPUT:

        - A List of lists -- the length of the List is the number of 
          boundary components of the surface. And each list in the List
          consists of an ordered tuple of numbers, each number comes
          from the number assigned to the corresponding dart before
          cutting.

        EXAMPLES:

        We start with a ribbon graph whose thickening has one boundary
        component. We compute its labeled boundary, then
        reduce it and compute the labeled boundary of the reduced 
        ribbon graph::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1);R1
            Sigma: [[1, 3, 5], [2, 4, 6]] 
            Rho: [[1, 2], [3, 4], [5, 6]]
            sage: R1.boundary()
            [[1, 2, 4, 3, 5, 6, 2, 1, 3, 4, 6, 5]]
            sage: H1 = R1.reduced(); H1
            Sigma: [[3, 5, 4, 6]] 
            Rho: [[3, 4], [5, 6]]
            sage: H1.boundary()
            [[3, 4, 6, 5, 4, 3, 5, 6]]

        We now consider a ribbon graph whose thickening has 3 boundary
        components. Also observe that in one of the labeled boundary
        components, a numbers appears twice in a row. That is because
        the ribbon graph has a vertex of valency 1::

            sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18,19)')
            sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)(19,20)')
            sage: R2 = RibbonGraph(s2,r2)
            sage: R2.n_boundary()
            3
            sage: R2.boundary()
            [[1, 16, 17, 4, 5, 14, 15, 8, 9, 12, 10, 3],
            [2, 13, 14, 5, 6, 11, 12, 9, 7, 18, 19, 20, 20, 19, 16, 1],
            [3, 10, 11, 6, 4, 17, 18, 7, 8, 15, 13, 2]]

        """

        #initialize and empty list to hold the labels of the boundaries
        bound = []

        #since lists of tuples are not modifiable, we change the data to a
        #list of lists 
        aux_perm = (self.rho*self.sigma).cycle_tuples(singletons = 1)

        #the cycles of the permutation rho*sigma are in 1:1 correspondence with 
        #the boundary components of the thickening (see function n_boundary())
        #but they are not the labeled boundary components.
        #With the next for, we convert the cycles of rho*sigma to actually 
        #the labelling of the edges. Each edge, therefore, should appear twice

        for i in range(len(aux_perm)):
            bound = bound + [[]]
            for j in range(len(aux_perm[i])):
                if self.rho(aux_perm[i][j]) != aux_perm[i][j]:
                    bound[i].append(aux_perm[i][j])
                    bound[i].append(self.rho(aux_perm[i][j]))
                else:
                    continue

        #finally the function returns a List of lists. Each list contains
        #a sequence of  numbers and each number corresponds to a half-edge.
        return _clean(bound)


    def reduced(self):
        r"""
        Return a ribbon graph with 1 vertex and $\mu$ edges (where $\mu$
        is the first betti number of the graph).

        OUTPUT:

        - A ribbon graph whose sigma permutation has only 1 non-singleton
          cycle and whose rho permutation is a product of $\mu$ disjoint
          2-cycles.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R1 = RibbonGraph(s1,r1)
            sage: R1.reduced()
            Sigma: [[3, 5, 4, 6]]
            Rho: [[3, 4], [5, 6]]

            sage: s2=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18,19)')
            sage: r2=PermutationGroupElement('(1,16)(2,13)(3,10)(4,17)(5,14)(6,11)(7,18)(8,15)(9,12)(19,20)')
            sage: R2 = RibbonGraph(s2,r2)
            sage: R2.reduced()
            Sigma: [[5, 6, 8, 9, 14, 15, 11, 12]] 
            Rho: [[5, 14], [6, 11], [8, 15], [9, 12]]

            sage: s3=PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15,16)(17,18,19,20)(21,22,23,24)')
            sage: r3=PermutationGroupElement('(1,21)(2,17)(3,13)(4,22)(7,23)(5,18)(6,14)(8,19)(9,15)(10,24)(11,20)(12,16)')
            sage: R3 = RibbonGraph(s3,r3);R3
            Sigma: [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15, 16], [17, 18, 19, 20], [21, 22, 23, 24]] 
            Rho: [[1, 21], [2, 17], [3, 13], [4, 22], [5, 18], [6, 14], [7, 23], [8, 19], [9, 15], [10, 24], [11, 20], [12, 16]]
            sage: R3.reduced()
            Sigma: [[5, 6, 8, 9, 11, 12, 18, 19, 20, 14, 15, 16]] 
            Rho: [[5, 18], [6, 14], [8, 19], [9, 15], [11, 20], [12, 16]]

        """

        #the following two lines convert the list of tuples to list of lists
        #we have to contract exactly n edges
        aux_ribbon = copy.deepcopy(self)
        aux_rho = [list(x) for 
                    x in aux_ribbon.rho.cycle_tuples()]

        #Observe that in the end we will have $\mu$ edges, so we
        #know exactly how many steps we will iterate
        while len(aux_rho) > self.mu():
            aux_sigma = [list(x) for 
                     x in aux_ribbon.sigma.cycle_tuples(singletons = 1)]
            aux_rho = [list(x) for 
                    x in aux_ribbon.rho.cycle_tuples()]
            for j in range(len(aux_rho)):
                if (_find(aux_sigma, aux_rho[j][0])[0] != 
                        _find(aux_sigma, aux_rho[j][1])[0]):
                    aux_ribbon = aux_ribbon.contract_edge(j)
                    aux_rho = [list(x) for 
                    x in aux_ribbon.rho.cycle_tuples()]
                    break
                else:
                    continue
        #finally we change the data to a list of tuples and return the
        #information as a ribbon graph. 
        return aux_ribbon
    
    #the next function computes a basis of homology, it uses
    #the previous function.

    def make_generic(self):
        r"""
        Return a ribbon graph equivalent to ''self'' but where every
        vertex has valency 3.

        OUTPUT:

        - A ribbon graph that is equivalent to ''self'' but is generic
          in the sense that all vertices have valency 3.

        EXAMPLES::

            sage: R = make_ribbon(1,3); R; R.genus(); R.n_boundary()
            Sigma: [[1, 2, 3, 9, 7], [4, 8, 10, 5, 6]] 
            Rho: [[1, 4], [2, 5], [3, 6], [7, 8], [9, 10]]
            1
            3
            sage: G = R.make_generic(); G
            Sigma: [[2, 3, 11], [5, 6, 13], [7, 8, 15], [9, 16, 17], [10, 14, 19], [12, 18, 21], [20, 22]] 
            Rho: [[2, 5], [3, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16], [17, 18], [19, 20], [21, 22]]
            sage: R.genus() == G.genus(); R.n_boundary() == G.n_boundary()
            True
            True

            sage: R = make_ribbon(5,4); R; R.genus(); R.n_boundary()
            Sigma: [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 27, 25, 23], [12, 24, 26, 28, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]] 
            Rho: [[1, 12], [2, 13], [3, 14], [4, 15], [5, 16], [6, 17], [7, 18], [8, 19], [9, 20], [10, 21], [11, 22], [23, 24], [25, 26], [27, 28]]
            5
            4
            sage: G = R.reduced(); G
            Sigma: [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 27, 25, 23, 24, 26, 28, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]] 
            Rho: [[2, 13], [3, 14], [4, 15], [5, 16], [6, 17], [7, 18], [8, 19], [9, 20], [10, 21], [11, 22], [23, 24], [25, 26], [27, 28]]
            sage: G.genus()==R.genus(); G.n_boundary() == R.n_boundary()
            True
            True

            sage: R = make_ribbon(0,6); R; R.genus(); R.n_boundary()
            Sigma: [[1, 11, 9, 7, 5, 3], [2, 4, 6, 8, 10, 12]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]]
            0
            6
            sage: G = R.reduced(); G
            Sigma: [[3, 4, 6, 8, 10, 12, 11, 9, 7, 5]] 
            Rho: [[3, 4], [5, 6], [7, 8], [9, 10], [11, 12]]
            sage: G.genus()==R.genus(); G.n_boundary() == R.n_boundary()
            True
            True

        """
        aux_ribbon = self.reduced()

        for i in range(2*aux_ribbon.mu() - 2):
            aux_ribbon = aux_ribbon.extrude_edge(i,0,2)

        return aux_ribbon



    def homology_basis(self):
        r"""
        Return an oriented basis of the first homology group of the 
        graph.

        OUTPUT:

        - A LIST of Lists of lists. The length of the LIST is `\mu` .
          Each List in the LIST corresponds to an element 
          of the basis. Each list in a List is just a 2-tuple which 
          corresponds to an ordered edge of the graph, that is a 
          2-cycle as an oriented tuple. The union of all the lists in
          a List is a circle contained in the graph.

        EXAMPLES::

            sage: R = make_ribbon(0,6); R; R.mu()
            Sigma: [[1, 11, 9, 7, 5, 3], [2, 4, 6, 8, 10, 12]] 
            Rho: [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]]
            5
            sage: R.homology_basis()
            [[[3, 4], [2, 1]],
            [[5, 6], [2, 1]],
            [[7, 8], [2, 1]],
            [[9, 10], [2, 1]],
            [[11, 12], [2, 1]]]

            sage: R = make_ribbon(1,1); R; R.mu()
            Sigma: [[1, 2, 3], [4, 5, 6]] 
            Rho: [[1, 4], [2, 5], [3, 6]]
            2
            sage: R.homology_basis()
            [[[2, 5], [4, 1]], [[3, 6], [4, 1]]]
            sage: H = R.reduced(); H
            Sigma: [[2, 3, 5, 6]] 
            Rho: [[2, 5], [3, 6]]
            sage: H.homology_basis()
            [[[2, 5]], [[3, 6]]]

            sage: s3 = PermutationGroupElement('(1,2,3,4,5,6,7,8,9,10,11,27,25,23)(12,24,26,28,13,14,15,16,17,18,19,20,21,22)')
            sage: r3 = PermutationGroupElement('(1,12)(2,13)(3,14)(4,15)(5,16)(6,17)(7,18)(8,19)(9,20)(10,21)(11,22)(23,24)(25,26)(27,28)')
            sage: R3 = RibbonGraph(s3,r3); R3; R3.mu()
            Sigma: [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 27, 25, 23], [12, 24, 26, 28, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]] 
            Rho: [[1, 12], [2, 13], [3, 14], [4, 15], [5, 16], [6, 17], [7, 18], [8, 19], [9, 20], [10, 21], [11, 22], [23, 24], [25, 26], [27, 28]]
            13
            sage: R3.homology_basis()
            [[[2, 13], [12, 1]],
            [[3, 14], [12, 1]],
            [[4, 15], [12, 1]],
            [[5, 16], [12, 1]],
            [[6, 17], [12, 1]],
            [[7, 18], [12, 1]],
            [[8, 19], [12, 1]],
            [[9, 20], [12, 1]],
            [[10, 21], [12, 1]],
            [[11, 22], [12, 1]],
            [[23, 24], [12, 1]],
            [[25, 26], [12, 1]],
            [[27, 28], [12, 1]]]
            sage: H3 = R3.reduced();H3
            Sigma: [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 27, 25, 23, 24, 26, 28, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]] 
            Rho: [[2, 13], [3, 14], [4, 15], [5, 16], [6, 17], [7, 18], [8, 19], [9, 20], [10, 21], [11, 22], [23, 24], [25, 26], [27, 28]]
            sage: H3.homology_basis()
            [[[2, 13]],
            [[3, 14]],
            [[4, 15]],
            [[5, 16]],
            [[6, 17]],
            [[7, 18]],
            [[8, 19]],
            [[9, 20]],
            [[10, 21]],
            [[11, 22]],
            [[23, 24]],
            [[25, 26]],
            [[27, 28]]]

        """

        aux_sigma = [list(x) for 
                     x in self.sigma.cycle_tuples(singletons = 1)]

        basis = [[list(x)] for 
                 x in self.reduced().rho.cycle_tuples()]

        #Now we define center as the set of edges that were contracted 
        #in reduced() this set is contractible and can be define as the 
        #complement of reduced_rho in rho

        center = [list(x) for x 
                  in self.rho.cycle_tuples() 
                  if (x not in self.reduced().rho.cycle_tuples())]

        #We define an auxiliary list 'vertices' that will contain the
        #vertices (cycles of sigma) corresponding to each half edge. 

        vertices=[]

        for i in range(len(basis)):
            vertices = vertices + [[]]
            basis[i].extend(copy.deepcopy(center))

            for j in range (len(basis[i])):
                vertices[i].append(_find(aux_sigma,basis[i][j][0])[0])
                vertices[i].append(_find(aux_sigma,basis[i][j][1])[0])
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
        #is just a 2-tuple which corresponds to an ''ordered'' edge of rho.

        return basis


def make_ribbon(g, r):
    r"""
    Return a ribbon graph whose thickening has genus g and r boundary
    components

    INPUT:

        - ``g`` -- non-negative integer representing the genus of the
          thickening.

        - ``r`` -- positive integer representing the number of boundary
          components of the thickening.

    OUTPUT:

        - A ribbon graph that has 2 vertices (two non-trivial cycles 
          in its sigma permutation) of valency `2g + r` and it has 
          `2*g + r` edges (and hence `4*g + 2*r` darts).

    EXAMPLES::

        sage: R = make_ribbon(0,1); R; R.genus(); R.n_boundary();
        Sigma: [[1], [2]] 
        Rho: [[1, 2]]
        0
        1

        sage: R = make_ribbon(0,5); R; R.genus(); R.n_boundary();
        Sigma: [[1, 9, 7, 5, 3], [2, 4, 6, 8, 10]] 
        Rho: [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]
        0
        5

        sage: R = make_ribbon(1,1); R; R.genus(); R.n_boundary();
        Sigma: [[1, 2, 3], [4, 5, 6]] 
        Rho: [[1, 4], [2, 5], [3, 6]]
        1
        1

        sage: R = make_ribbon(7,3); R; R.genus(); R.n_boundary();
        Sigma: [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 33, 31], [16, 32, 34, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]] 
        Rho: [[1, 16], [2, 17], [3, 18], [4, 19], [5, 20], [6, 21], [7, 22], [8, 23], [9, 24], [10, 25], [11, 26], [12, 27], [13, 28], [14, 29], [15, 30], [31, 32], [33, 34]]
        7
        3

    """
    #Initialize the two vertices of sigma and the edge joining them
    repr_sigma= [[1],[2*g+2]]
    repr_rho = [[1,2*g+2]]

    #We first generate the surface of genus g and boundary component.
    #This is done by considering the usual planar representation of
    #a surface as a poligon of 4*g+2 edges with identifications. (see
    #any topology  book on the classification of surfaces)
    for i in range(2*g):
        repr_sigma[0].append(i+2)
        repr_sigma[1].append(i+(2*g+2)+1)
        repr_rho +=[[i+2,i+(2*g+2)+1]]

    #finally we add an edge for each aditional boundary component. 
    max_dart = 4*g+2
    for j in range(r-1):
        repr_sigma[0].insert(0,max_dart+2*(j+1)-1)
        repr_sigma[1].insert(j+1,max_dart+2*(j+1))
        repr_rho +=[[max_dart+2*(j+1)-1,max_dart+2*(j+1)]]

    return RibbonGraph(
                       PermutationGroupElement([tuple(x) for x in repr_sigma]), 
                       PermutationGroupElement([tuple(x) for x in repr_rho])
                      )