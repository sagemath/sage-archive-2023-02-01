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

