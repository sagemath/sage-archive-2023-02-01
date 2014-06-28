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
        raise ValueError('The input matrix has more columns than rows.')
    elif m == 0:
        return mat
    else:
        return mat.submatrix(0,0,n,n)

def _digraph_mutate( dg, k, n, m ):
    """
    Returns a digraph obtained from dg with n+m vertices by mutating at vertex k.

    INPUT:

    - ``dg`` -- a digraph with integral edge labels with ``n+m`` vertices
    - ``k`` -- the vertex at which ``dg`` is mutated

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_mutate
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
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

def _matrix_to_digraph( M ):
    """
    Returns the digraph obtained from the matrix ``M``.
    In order to generate a quiver, we assume that ``M`` is skew-symmetrizable.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _matrix_to_digraph
        sage: _matrix_to_digraph(matrix(3,[0,1,0,-1,0,-1,0,1,0]))
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

def _dg_canonical_form( dg, n, m ):
    """
    Turns the digraph ``dg`` into its canonical form, and returns the corresponding isomorphism, and the vertex orbits of the automorphism group.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_canonical_form
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
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

def _mutation_class_iter( dg, n, m, depth=infinity, return_dig6=False, show_depth=False, up_to_equivalence=True, sink_source=False ):
    """
    Returns an iterator for mutation class of dg with respect to several parameters.

    INPUT:

    - ``dg`` -- a digraph with n+m vertices
    - ``depth`` -- a positive integer or infinity specifying (roughly) how many steps away from the initial seed to mutate
    - ``return_dig6`` -- indicates whether to convert digraph data to dig6 string data
    - ``show_depth`` -- if True, indicates that a running count of the depth is to be displayed
    - ``up_to_equivalence``  -- if True, only one digraph for each graph-isomorphism class is recorded
    - ``sink_source`` -- if True, only mutations at sinks or sources are applied

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _mutation_class_iter
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',[1,2],1]).digraph()
        sage: itt = _mutation_class_iter(dg, 3,0)
        sage: itt.next()[0].edges()
        [(0, 1, (1, -1)), (0, 2, (1, -1)), (1, 2, (1, -1))]
        sage: itt.next()[0].edges()
        [(0, 2, (1, -1)), (1, 0, (2, -2)), (2, 1, (1, -1))]
    """
    timer = time.time()
    depth_counter = 0
    if up_to_equivalence:
        iso, orbits = _dg_canonical_form( dg, n, m )
        iso_inv = dict( (iso[a],a) for a in iso )

    dig6 = _digraph_to_dig6( dg, hashable=True )
    dig6s = {}
    if up_to_equivalence:
        orbits = [ orbit[0] for orbit in orbits ]
        dig6s[dig6] = [ orbits, [], iso_inv ]
    else:
        dig6s[dig6] = [ range(n), [] ]
    if return_dig6:
        yield (dig6, [])
    else:
        yield (dg, [])

    gets_bigger = True
    if show_depth:
        timer2 = time.time()
        dc = str(depth_counter)
        dc += ' ' * (5-len(dc))
        nr = str(len(dig6s))
        nr += ' ' * (10-len(nr))
        print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)

    while gets_bigger and depth_counter < depth:
        gets_bigger = False
        keys = dig6s.keys()
        for key in keys:
            mutation_indices = [ i for i in dig6s[key][0] if i < n ]
            if mutation_indices:
                dg = _dig6_to_digraph( key )
            while mutation_indices:
                i = mutation_indices.pop()
                if not sink_source or _dg_is_sink_source( dg, i ):
                    dg_new = _digraph_mutate( dg, i, n, m )
                    if up_to_equivalence:
                        iso, orbits = _dg_canonical_form( dg_new, n, m )
                        i_new = iso[i]
                        iso_inv = dict( (iso[a],a) for a in iso )
                    else:
                        i_new = i
                    dig6_new = _digraph_to_dig6( dg_new, hashable=True )
                    if dig6_new in dig6s:
                        if i_new in dig6s[dig6_new][0]:
                            dig6s[dig6_new][0].remove(i_new)
                    else:
                        gets_bigger = True
                        if up_to_equivalence:
                            orbits = [ orbit[0] for orbit in orbits if i_new not in orbit ]
                            iso_history = dict( (a, dig6s[key][2][iso_inv[a]]) for a in iso )
                            i_history = iso_history[i_new]
                            history = dig6s[key][1] + [i_history]
                            dig6s[dig6_new] = [orbits,history,iso_history]
                        else:
                            orbits = range(n)
                            del orbits[i_new]
                            history = dig6s[key][1] + [i_new]
                            dig6s[dig6_new] = [orbits,history]
                        if return_dig6:
                            yield (dig6_new,history)
                        else:
                            yield (dg_new,history)
        depth_counter += 1
        if show_depth and gets_bigger:
            timer2 = time.time()
            dc = str(depth_counter)
            dc += ' ' * (5-len(dc))
            nr = str(len(dig6s))
            nr += ' ' * (10-len(nr))
            print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)

def _digraph_to_dig6( dg, hashable=False ):
    """
    Returns the dig6 and edge data of the digraph dg.

    INPUT:

    - ``dg`` -- a digraph
    - ``hashable`` -- (Boolean; optional; default:False) if ``True``, the edge labels are turned into a dict.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: _digraph_to_dig6(dg)
        ('COD?', {})
    """
    dig6 = dg.dig6_string()
    D = {}
    for E in dg._backend.iterator_in_edges(dg,True):
        if E[2] != (1,-1):
            D[ (E[0],E[1]) ] = E[2]
    if hashable:
        D = tuple( sorted( D.items() ) )
    return (dig6,D)

def _dig6_to_digraph( dig6 ):
    """
    Returns the digraph obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dig6_to_digraph
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: data = _digraph_to_dig6(dg)
        sage: _dig6_to_digraph(data)
        Digraph on 4 vertices
        sage: _dig6_to_digraph(data).edges()
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]
    """
    dig6, edges = dig6
    dg = DiGraph( dig6 )
    if not isinstance(edges, dict):
        edges = dict( edges )
    for edge in dg._backend.iterator_in_edges(dg,False):
        if edge in edges:
            dg.set_edge_label( edge[0],edge[1],edges[edge] )
        else:
            dg.set_edge_label( edge[0],edge[1], (1,-1) )
    return dg

def _dig6_to_matrix( dig6 ):
    """
    Returns the matrix obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6, _dig6_to_matrix
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: data = _digraph_to_dig6(dg)
        sage: _dig6_to_matrix(data)
        [ 0  1  0  0]
        [-1  0 -1  0]
        [ 0  1  0  1]
        [ 0  0 -1  0]
    """
    dg = _dig6_to_digraph( dig6 )
    return _edge_list_to_matrix( dg.edges(), dg.order(), 0 )

def _dg_is_sink_source( dg, v ):
    """
    Returns True iff the digraph dg has a sink or a source at vertex v.

    INPUT:

    - ``dg`` -- a digraph
    - ``v`` -- a vertex of dg

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_is_sink_source
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',[1,2],1]).digraph()
        sage: _dg_is_sink_source(dg, 0 )
        True
        sage: _dg_is_sink_source(dg, 1 )
        True
        sage: _dg_is_sink_source(dg, 2 )
        False
    """
    in_edges = [ edge for edge in dg._backend.iterator_in_edges([v],True) ]
    out_edges = [ edge for edge in dg._backend.iterator_out_edges([v],True) ]
    return not ( in_edges and out_edges )

def _graph_without_edge_labels(dg,vertices):
    """
    Replaces edge labels in dg other than ``(1,-1)`` by this edge label, and returns the corresponding partition of the edges.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _graph_without_edge_labels
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
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

def _has_two_cycles( dg ):
    """
    Returns True if the input digraph has a 2-cycle and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _has_two_cycles
        sage: _has_two_cycles( DiGraph([[0,1],[1,0]]))
        True
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: _has_two_cycles( ClusterQuiver(['A',3]).digraph() )
        False
    """
    edge_set = dg.edges(labels=False)
    for (v,w) in edge_set:
        if (w,v) in edge_set:
            return True
    return False

def _is_valid_digraph_edge_set( edges, frozen=0 ):
    """
    Returns True if the input data is the edge set of a digraph for a quiver (no loops, no 2-cycles, edge-labels of the specified format), and returns False otherwise.

    INPUT:

    - ``frozen`` -- (integer; default:0) The number of frozen vertices.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _is_valid_digraph_edge_set
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)]] )
        The given digraph has edge labels which are not integral or integral 2-tuples.
        False
        sage: _is_valid_digraph_edge_set( [[0,1,None],[2,3,(1,-1)]] )
        True
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)],[3,2,(1,-1)]] )
        The given digraph or edge list contains oriented 2-cycles.
        False
    """
    try:
        dg = DiGraph()
        dg.allow_multiple_edges(True)
        dg.add_edges( edges )

        # checks if the digraph contains loops
        if dg.has_loops():
            print "The given digraph or edge list contains loops."
            return False

        # checks if the digraph contains oriented 2-cycles
        if _has_two_cycles( dg ):
            print "The given digraph or edge list contains oriented 2-cycles."
            return False

        # checks if all edge labels are 'None', positive integers or tuples of positive integers
        if not all( i is None or ( i in ZZ and i > 0 ) or ( isinstance(i, tuple) and len(i) == 2 and i[0] in ZZ and i[1] in ZZ ) for i in dg.edge_labels() ):
            print "The given digraph has edge labels which are not integral or integral 2-tuples."
            return False

        # checks if all edge labels for multiple edges are 'None' or positive integers
        if dg.has_multiple_edges():
            for e in set( dg.multiple_edges(labels=False) ):
                if not all( i is None or ( i in ZZ and i > 0 ) for i in dg.edge_label( e[0], e[1] ) ):
                    print "The given digraph or edge list contains multiple edges with non-integral labels."
                    return False

        n = dg.order() - frozen
        if n < 0:
            print "The number of frozen variables is larger than the number of vertices."
            return False

        if [ e for e in dg.edges(labels=False) if e[0] >= n] != []:
            print "The given digraph or edge list contains edges within the frozen vertices."
            return False

        return True
    except Exception:
        print "Could not even build a digraph from the input data."
        return False
