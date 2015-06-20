# distutils: libraries = gmp
r"""
Cutwidth

This module implements several algorithms to compute the cutwidth of a graph and
the corresponding ordering of the vertices. It also implements tests functions
for evaluation the width of a linear ordering (or layout).

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n-1} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |\{(u,w)\in E(G)\mid u\in S\text{ and }w\in V(G)\backslash S\}|

The *cutwidth* of a graph `G` is equal to the minimum cost of an ordering of its
vertices.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`cutwidth` | Returns the cutwidth of the given graph and the corresponding linear ordering of the vertices
    :meth:`cutwidth_dyn` | Computes the cutwidth of `G` using an exponential time and space algorithm based on dynamic programming
    :meth:`width_of_cut_decomposition` | Returns the width of the cut decomposition induced by the linear ordering `L` of the vertices of `G`


Exponential algorithm for cutwidth
----------------------------------

This algorithm differs from
:meth:`sage.graphs.graph_decompositions.vertex_separation.vertex_separation_exp`
only in the cost function. See the :mod:`vertex separation module's
documentation <sage.graphs.graph_decompositions.vertex_separation>` for more
details on this algorithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on strictly less than 32 vertices. This can be changed to 64 if necessary,
    but 32 vertices already require 4GB of memory.


Authors
-------

- David Coudert (2015-06): Initial version


Methods
-------
"""
#*****************************************************************************
#          Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#      as published by the Free Software Foundation; either version 2 of
#              the License, or (at your option) any later version.
#                        http://www.gnu.org/licenses/
#*****************************************************************************




################################################################################
# Function for testing the validity of a linear vertex ordering
#
# A linear ordering `L` of the vertices of a graph `G` is valid if all vertices
# of `G` are in `L`, and if `L` contains no other vertex and no duplicated
# vertices.
################################################################################
from sage.graphs.graph_decompositions.vertex_separation import is_valid_ordering


################################################################################
# Measurement function of the width of some layout
################################################################################

def width_of_cut_decomposition(G, L):
    r"""
    Returns the width of the cut decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `cw_L(G)` of the cut decomposition induced by the linear
    ordering `L` of the vertices of `G`.

    .. MATH::

        cw_L(G) =  \max_{0\leq i< |V|-1} |\{(u,w)\in E(G)\mid u\in L[:i]\text{ and }w\in V(G)\setminus L[:i]\}|

    INPUT:

    - ``G`` -- a Graph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Cut decomposition of a Cycle graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices()
        sage: cutwidth.width_of_cut_decomposition(G, L)
        2

    Cut decomposition of a Path graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: P = graphs.PathGraph(6)
        sage: cutwidth.width_of_cut_decomposition(P, [0, 1, 2, 3, 4, 5])
        1
        sage: cutwidth.width_of_cut_decomposition(P, [5, 0, 1, 2, 3, 4])
        2
        sage: cutwidth.width_of_cut_decomposition(P, [0, 2, 4, 1, 3, 5])
        5

    TESTS:

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.width_of_cut_decomposition(Graph(), ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: The input linear vertex ordering L is not valid for G.
    """
    if not is_valid_ordering(G, L):
        raise ValueError("The input linear vertex ordering L is not valid for G.")

    position = dict( (u,i) for i,u in enumerate(L) )

    cpt = [0]*len(L)
    for u,v in G.edge_iterator(labels=None):
        x,y = position[u],position[v]
        if x>y:
            x,y = y,x
        for i in range(x,y):
            cpt[i] += 1

    return max(cpt)


