r"""
mutation_class

This file contains helper functions for compute the mutation class of a cluster algebra or quiver.

For the compendium on the cluster algebra and quiver package see

        http://arxiv.org/abs/1102.4844

AUTHORS:

- Gregg Musiker
- Christian Stump
"""

#*****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import time
from sage.groups.perm_gps.partn_ref.refinement_graphs import *
from sage.graphs.generic_graph import graph_isom_equivalent_non_edge_labeled_graph
from copy import copy
from sage.rings.all import ZZ, infinity
from sage.graphs.all import Graph, DiGraph
from sage.matrix.all import matrix
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _edge_list_to_matrix

#################################################
#                                               #
#   Methods to construct the mutation class     #
#      of a quiver and its dig6 data            #
#                                               #
#################################################

def _digraph_mutate( dg, k, n, m ):
    """
    Returns a digraph obtained from dg with n+m vertices by mutating at vertex k.

    INPUT:

    - ``dg`` -- a digraph with integral edge labels with ``n+m`` vertices
    - ``k`` -- the vertex at which ``dg`` is mutated

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_mutate
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: dg.edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]
        sage: _digraph_mutate(dg,2,4,0).edges()
        [(0, 1, (1, -1)), (1, 2, (1, -1)), (3, 2, (1, -1))]
    """
    edges = dict( ((v1,v2),label) for v1,v2,label in dg._backend.iterator_in_edges(dg,True) )
    in_edges = [ (v1,v2,edges[(v1,v2)]) for (v1,v2) in edges if v2 == k ]
    out_edges = [ (v1,v2,edges[(v1,v2)]) for (v1,v2) in edges if v1 == k ]
    in_edges_new = [ (v2,v1,(-label[1],-label[0])) for (v1,v2,label) in in_edges ]
    out_edges_new = [ (v2,v1,(-label[1],-label[0])) for (v1,v2,label) in out_edges ]
    diag_edges_new = []
    diag_edges_del = []

    for (v1,v2,label1) in in_edges:
        for (w1,w2,label2) in out_edges:
            l11,l12 = label1
            l21,l22 = label2
            if (v1,w2) in edges:
                diag_edges_del.append( (v1,w2,edges[(v1,w2)]) )
                a,b = edges[(v1,w2)]
                a,b = a+l11*l21, b-l12*l22
                diag_edges_new.append( (v1,w2,(a,b)) )
            elif (w2,v1) in edges:
                diag_edges_del.append( (w2,v1,edges[(w2,v1)]) )
                a,b = edges[(w2,v1)]
                a,b = b+l11*l21, a-l12*l22
                if a<0:
                    diag_edges_new.append( (w2,v1,(b,a)) )
                elif a>0:
                    diag_edges_new.append( (v1,w2,(a,b)) )
            else:
                a,b = l11*l21,-l12*l22
                diag_edges_new.append( (v1,w2,(a,b)) )

    del_edges = in_edges + out_edges + diag_edges_del
    new_edges = in_edges_new + out_edges_new + diag_edges_new
    new_edges += [ (v1,v2,edges[(v1,v2)]) for (v1,v2) in edges if not (v1,v2,edges[(v1,v2)]) in del_edges ]

    dg_new = DiGraph()
    for v1,v2,label in new_edges:
        dg_new._backend.add_edge(v1,v2,label,True)
    if dg_new.order() < n+m:
        dg_new_vertices = [ v for v in dg_new ]
        for i in [ v for v in dg if v not in dg_new_vertices ]:
            dg_new.add_vertex(i)

    return dg_new

def _dg_canonical_form( dg, n, m ):
    """
    Turns the digraph ``dg`` into its canonical form, and returns the corresponding isomorphism, and the vertex orbits of the automorphism group.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_canonical_form
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _dg_canonical_form(dg,4,0); dg.edges()
        ({0: 0, 1: 3, 2: 1, 3: 2}, [[0], [3], [1], [2]])
        [(0, 3, (1, -1)), (1, 2, (1, -2)), (1, 3, (1, -1))]
    """
    vertices = [ v for v in dg ]
    if m > 0:
        partition = [ vertices[:n], vertices[n:] ]
    else:
        partition = [ vertices ]
    partition_add, edges = _graph_without_edge_labels(dg,vertices)
    partition += partition_add
    automorphism_group, obsolete, iso = search_tree(dg, partition=partition, lab=True, dig=True, certify=True)
    orbits = get_orbits( automorphism_group, n+m )
    orbits = [ [ iso[i] for i in orbit] for orbit in orbits ]
    for v in iso.keys():
        if v >= n+m:
            del iso[v]
            v1,v2,label1 = dg._backend.iterator_in_edges([v],True).next()
            w1,w2,label2 = dg._backend.iterator_out_edges([v],True).next()
            dg._backend.del_edge(v1,v2,label1,True)
            dg._backend.del_edge(w1,w2,label2,True)
            dg._backend.del_vertex(v)
            add_index = True
            index = 0
            while add_index:
                l = partition_add[index]
                if v in l:
                    dg._backend.add_edge(v1,w2,edges[index],True)
                    add_index = False
                index += 1
    dg._backend.relabel( iso, True )
    return iso, orbits

def _graph_without_edge_labels(dg,vertices):
    """
    Replaces edge labels in dg other than ``(1,-1)`` by this edge label, and returns the corresponding partition of the edges.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _graph_without_edge_labels
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _graph_without_edge_labels(dg,range(3)); dg.edges()
        ([[5]], [(1, -2)])
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 5, (1, -1)), (5, 3, (1, -1))]
    """
    edges = [ edge for edge in dg._backend.iterator_in_edges(dg,True) ]
    edge_labels = sorted([ label for v1,v2,label in edges if not label == (1,-1)])
    i = 1
    while i < len(edge_labels):
        if edge_labels[i] == edge_labels[i-1]:
            edge_labels.pop(i)
        else:
            i += 1
    edge_partition = [[] for _ in xrange(len(edge_labels))]
    i = 0
    new_vertices = []
    for u,v,l in edges:
        while i in vertices or i in new_vertices:
            i += 1
        new_vertices.append(i)
        if not l == (1,-1):
            index = edge_labels.index(l)
            edge_partition[index].append(i)
            dg._backend.add_edge(u,i,(1,-1),True)
            dg._backend.add_edge(i,v,(1,-1),True)
            dg._backend.del_edge(u,v,l,True)
    return [a for a in edge_partition if a != []], edge_labels

#################################################
#                                               #
#   Old methods to translate between matrices,  #
#        digraphs, dig6 data, and quivers       #
#                                               #
#################################################

def _matrix_to_digraph( M ):
    """
    Returns the digraph obtained from the matrix M.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _matrix_to_digraph

        sage: mat = matrix(3,[0,1,0,-1,0,-1,0,1,0]); mat
        [ 0  1  0]
        [-1  0 -1]
        [ 0  1  0]

        sage: _matrix_to_digraph(mat)
        Digraph on 3 vertices
    """
    n = M.ncols()

    dg = DiGraph(sparse=True)
    for i,j in M.nonzero_positions():
        if i >= n: a,b = M[i,j],-M[i,j]
        else: a,b = M[i,j],M[j,i]
        if a > 0:
            dg._backend.add_edge(i,j,(a,b),True)
        elif i >= n:
            dg._backend.add_edge(j,i,(-a,-b),True)
    if dg.order() < M.nrows():
        for i in [ index for index in xrange(M.nrows()) if index not in dg ]:
            dg.add_vertex(i)
    return dg

def _principal_part( mat ):
    """
    Returns the principal part of a matrix.

    INPUT:

    - ``mat`` - a matrix with at least as many rows as columns

    OUTPUT:

    The top square part of the matrix ``mat``.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part
        sage: M = Matrix([[1,2],[3,4],[5,6]]); M
        [1 2]
        [3 4]
        [5 6]
        sage: _principal_part(M)
        [1 2]
        [3 4]
    """
    n, m = mat.ncols(), mat.nrows()-mat.ncols()
    if m < 0:
        raise ValueError, 'The input matrix has more columns than rows.'
    elif m == 0:
        return mat
    else:
        return mat.submatrix(0,0,n,n)

