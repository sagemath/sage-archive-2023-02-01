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

# import from Python standard library
from math import sin, cos, pi

# import from Sage library
from sage.graphs.graph import Graph


def RingedTree(self, k):
    r"""
    Return the ringed tree on k-levels.

    A ringed tree of level `k` is a fully binary tree with `k` levels (counting
    the root as a level), in which all vertices at the same level are connected
    by a ring. More precisely, we can use a binary string to represent each
    vertex in the tree, such that the root (at level 0) is represented by an
    empty string, and the left child and the right child of a vertex with string
    `\sigma` are represented as `\sigma0` and \sigma1`, respectively. Then, at
    each level `i = 1, 2, \cdots , k-1`, we connect two vertices `u` and `v`
    represented by binary strings `\sigma_u` and `\sigma_v` if `(\sigma_u + 1)
    \pmod{2^i}\equiv \sigma_v`, where the addition treats the binary strings as
    the integers they represent.

    INPUT:

    - ``k`` -- the number of levels of the ringed tree.

    EXAMPLE:

        sage: G = graphs.RingedTree(5)
        sage: P = G.plot(vertex_labels=False, vertex_size=10)
        sage: P.show() # long time

    TEST:

        sage: G = graphs.RingedTree(-1)
        Traceback (most recent call last):
        ...
        ValueError: The number of levels must be >= 1.
    """
    if k<1:
        raise ValueError('The number of levels must be >= 1.')

    E = []
    W = {}
    # We first add the tree edges
    W[0] = ['']
    for l in xrange(k-1):
        W[l+1] = []
        for r in W[l]:
            E.append( (r,r+'0') )
            E.append( (r,r+'1') )
            W[l+1].append(r+'0')
            W[l+1].append(r+'1')

    # We then add the ring edges and fix positions of vertices for nice drawing
    pos_dict = {'':(0,0.2)}
    for l in xrange(1,k):
        twol = 1<<l
        r = 1 if l==1 else 1.5**l
        theta0 = pi/twol
        pos_dict[ W[l][twol-1] ] = ( r*sin(-theta0), r*cos(-theta0) )
        if l>1: # to avoid adding twice edge ('0','1')
            E.append( (W[l][0],W[l][twol-1]) )
        for i in xrange(twol-1):
            theta = (2*i+1)*theta0
            pos_dict[ W[l][i] ] = ( r*sin(theta), r*cos(theta) )
            E.append( (W[l][i],W[l][i+1]) )

    return Graph(E, pos=pos_dict, name='Ringed Tree on '+str(k)+' levels' )

