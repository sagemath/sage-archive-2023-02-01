r"""
Families of graphs
==================

This file gathers generators for some families of graphs.
- :meth:`RingedTree <GraphGenerators.RingedTree>`

AUTHORS:

- David Coudert (2012) Ringed Trees
"""

################################################################################
#           Copyright (C) 2012 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################

# import from Sage library
from sage.graphs.graph import Graph

def RingedTree(self, k, vertex_labels = True):
    r"""
    Return the ringed tree on k-levels.

    A ringed tree of level `k` is a binary tree with `k` levels (counting
    the root as a level), in which all vertices at the same level are connected
    by a ring.

    More precisely, in each layer of the binary tree (i.e. a layer is the set of
    vertices `[2^i...2^{i+1}-1]`) two vertices `u,v` are adjacent if `u=v+1` or
    if `u=2^i` and `v=`2^{i+1}-1`.

    Ringed trees are defined in [CFHM12]_.

    INPUT:

    - ``k`` -- the number of levels of the ringed tree.

    - ``vertex_labels`` (boolean) -- whether to label vertices as binary words
      (default) or as integers.

    EXAMPLE::

        sage: G = graphs.RingedTree(5)
        sage: P = G.plot(vertex_labels=False, vertex_size=10)
        sage: P.show() # long time
        sage: G.vertices()
        ['', '0', '00', '000', '0000', '0001', '001', '0010', '0011', '01',
         '010', '0100', '0101', '011', '0110', '0111', '1', '10', '100',
         '1000', '1001', '101', '1010', '1011', '11', '110', '1100', '1101',
         '111', '1110', '1111']

    TEST::

        sage: G = graphs.RingedTree(-1)
        Traceback (most recent call last):
        ...
        ValueError: The number of levels must be >= 1.
        sage: G = graphs.RingedTree(5, vertex_labels = False)
        sage: G.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

    REFERENCES:

    .. [CFHM12] On the Hyperbolicity of Small-World and Tree-Like Random Graphs
      Wei Chen, Wenjie Fang, Guangda Hu, Michael W. Mahoney
      http://arxiv.org/abs/1201.1717
    """
    if k<1:
        raise ValueError('The number of levels must be >= 1.')

    from sage.graphs.graph_generators import _circle_embedding
    from sage.graphs.graph_generators import GraphGenerators

    # Creating the Balanced tree, which contains most edges already
    g = GraphGenerators().BalancedTree(2,k-1)
    g.name('Ringed Tree on '+str(k)+' levels')

    # We consider edges layer by layer
    for i in range(1,k):
        vertices = range(2**(i)-1,2**(i+1)-1)

        # Add the missing edges
        g.add_cycle(vertices)

        # And set the vertices' positions
        radius = i if i <= 1 else 1.5**i
        shift = -2**(i-2)+.5 if i > 1 else 0
        _circle_embedding(g, vertices, radius = radius, shift = shift)

    # Specific position for the central vertex
    g.get_pos()[0] = (0,0.2)

    # Relabel vertices as binary words
    if not vertex_labels:
        return g

    vertices = ['']
    for i in range(k-1):
        for j in range(2**(i)-1,2**(i+1)-1):
            v = vertices[j]
            vertices.append(v+'0')
            vertices.append(v+'1')

    g.relabel(vertices)

    return g
