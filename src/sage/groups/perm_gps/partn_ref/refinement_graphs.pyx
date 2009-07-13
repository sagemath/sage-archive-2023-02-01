"""
Graph-theoretic partition backtrack functions

DOCTEST:
    sage: import sage.groups.perm_gps.partn_ref.refinement_graphs

REFERENCE:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

"""

#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'data_structures_pyx.pxi' # includes bitsets

def isomorphic(G1, G2, partn, ordering2, dig, use_indicator_function, sparse=False):
    """
    Tests whether two graphs are isomorphic.

    sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic

    sage: G = Graph(2)
    sage: H = Graph(2)
    sage: isomorphic(G, H, [[0,1]], [0,1], 0, 1)
    {0: 0, 1: 1}
    sage: isomorphic(G, H, [[0,1]], [0,1], 0, 1)
    {0: 0, 1: 1}
    sage: isomorphic(G, H, [[0],[1]], [0,1], 0, 1)
    {0: 0, 1: 1}
    sage: isomorphic(G, H, [[0],[1]], [1,0], 0, 1)
    {0: 1, 1: 0}

    sage: G = Graph(3)
    sage: H = Graph(3)
    sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
    {0: 0, 1: 1, 2: 2}
    sage: G.add_edge(0,1)
    sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
    False
    sage: H.add_edge(1,2)
    sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
    {0: 1, 1: 2, 2: 0}

    """
    cdef int **part
    cdef int *output, *ordering
    cdef CGraph G
    cdef GraphStruct GS1 = GraphStruct()
    cdef GraphStruct GS2 = GraphStruct()
    cdef GraphStruct GS
    cdef int i, j, n = -1

    from sage.graphs.graph import GenericGraph, Graph, DiGraph
    which_G = 1
    for G_in in [G1, G2]:
        if which_G == 1:
            GS = GS1
            first=True
        else:
            GS = GS2
            first=False
        if isinstance(G_in, GenericGraph):
            if n == -1:
                n = G_in.num_verts()
            elif n != G_in.num_verts():
                return False
            G_in = G_in.copy()
            if G_in.vertices() != range(n):
                to = G_in.relabel(return_map=True)
                frm = {}
                for v in to.iterkeys():
                    frm[to[v]] = v
                if first:
                    partition = [[to[v] for v in cell] for cell in partn]
            else:
                to = range(n)
                frm = to
            if sparse:
                G = SparseGraph(n)
            else:
                G = DenseGraph(n)
            if G_in.is_directed():
                for i,j in G_in.edge_iterator(labels=False):
                    G.add_arc(i,j)
            else:
                for i,j in G_in.edge_iterator(labels=False):
                    G.add_arc(i,j)
                    G.add_arc(j,i)
        elif isinstance(G_in, CGraph):
            G = <CGraph> G_in
            if n == -1:
                n = G.num_verts
            elif n != G.num_verts:
                return False
            to = {}
            for a in G.verts(): to[a]=a
            frm = to
        else:
            raise TypeError("G must be a Sage graph.")
        if first: frm1=frm;to1=to
        else: frm2=frm;to2=to
        GS.G = G
        GS.directed = 1 if dig else 0
        GS.use_indicator = 1 if use_indicator_function else 0
        which_G += 1

    if n == 0:
        return {}

    part = <int **> sage_malloc((len(partn)+1) * sizeof(int *))
    ordering = <int *> sage_malloc(n * sizeof(int))
    if part is NULL or ordering is NULL:
        if part is not NULL: sage_free(part)
        if ordering is not NULL: sage_free(ordering)
        raise MemoryError
    for i from 0 <= i < len(partn):
        part[i] = <int *> sage_malloc((len(partn[i])+1) * sizeof(int))
        if part[i] is NULL:
            for j from 0 <= j < i:
                sage_free(part[j])
            sage_free(part)
            raise MemoryError
        for j from 0 <= j < len(partn[i]):
            part[i][j] = to1[partn[i][j]]
        part[i][len(partn[i])] = -1
    part[len(partn)] = NULL
    for i from 0 <= i < n:
        ordering[i] = to2[ordering2[i]]

    GS1.scratch = <int *> sage_malloc((3*n+1) * sizeof(int))
    GS2.scratch = <int *> sage_malloc((3*n+1) * sizeof(int))
    if GS1.scratch is NULL or GS2.scratch is NULL:
        if GS1.scratch is not NULL: sage_free(GS1.scratch)
        if GS2.scratch is not NULL: sage_free(GS2.scratch)
        for j from 0 <= j < len(partition):
            sage_free(part[j])
        sage_free(part)
        raise MemoryError

    output = double_coset(GS1, GS2, part, ordering, n, &all_children_are_equivalent, &refine_by_degree, &compare_graphs)

    for i from 0 <= i < len(partn):
        sage_free(part[i])
    sage_free(part)
    sage_free(ordering)
    sage_free(GS1.scratch)
    sage_free(GS2.scratch)

    if output is NULL:
        return False
    else:
        output_py = dict([[frm1[i], frm2[output[i]]] for i from 0 <= i < n])
        sage_free(output)
        return output_py

def search_tree(G_in, partition, lab=True, dig=False, dict_rep=False, certify=False,
                    verbosity=0, use_indicator_function=True, sparse=True,
                    base=False, order=False):
    """
    Compute canonical labels and automorphism groups of graphs.

    INPUT:
    G_in -- a Sage graph
    partition -- a list of lists representing a partition of the vertices
    lab -- if True, compute and return the canonical label in addition to the
        automorphism group.
    dig -- set to True for digraphs and graphs with loops.  If True, does not
        use optimizations based on Lemma 2.25 in [1] that are valid only for
        simple graphs.
    dict_rep -- if True, return a dictionary with keys the vertices of the
        input graph G_in and values elements of the set the permutation group
        acts on.  (The point is that graphs are arbitrarily labelled, often
        0..n-1, and permutation groups always act on 1..n.  This dictionary
        maps vertex labels (such as 0..n-1) to the domain of the permutations.)
    certify -- if True, return the permutation from G to its canonical label.
    verbosity -- currently ignored
    use_indicator_function -- option to turn off indicator function
        (True is generally faster)
    sparse -- whether to use sparse or dense representation of the graph
        (ignored if G is already a CGraph - see sage.graphs.base)
    base -- whether to return the first sequence of split vertices (used in
        computing the order of the group)
    order -- whether to return the order of the automorphism group

    OUTPUT:
    Depends on the options. If more than one thing is returned, they are in a
    tuple in the following order:
        list of generators in list-permutation format -- always
        dict -- if dict_rep
        graph -- if lab
        dict -- if certify
        list -- if base
        integer -- if order

    DOCTEST:
        sage: st = sage.groups.perm_gps.partn_ref.refinement_graphs.search_tree
        sage: from sage.graphs.base.dense_graph import DenseGraph
        sage: from sage.graphs.base.sparse_graph import SparseGraph

    Graphs on zero vertices:
        sage: G = Graph()
        sage: st(G, [[]], order=True)
        ([], Graph on 0 vertices, 1)

    Graphs on one vertex:
        sage: G = Graph(1)
        sage: st(G, [[0]], order=True)
        ([], Graph on 1 vertex, 1)

    Graphs on two vertices:
        sage: G = Graph(2)
        sage: st(G, [[0,1]], order=True)
        ([[1, 0]], Graph on 2 vertices, 2)
        sage: st(G, [[0],[1]], order=True)
        ([], Graph on 2 vertices, 1)
        sage: G.add_edge(0,1)
        sage: st(G, [[0,1]], order=True)
        ([[1, 0]], Graph on 2 vertices, 2)
        sage: st(G, [[0],[1]], order=True)
        ([], Graph on 2 vertices, 1)

    Graphs on three vertices:
        sage: G = Graph(3)
        sage: st(G, [[0,1,2]], order=True)
        ([[0, 2, 1], [1, 0, 2]], Graph on 3 vertices, 6)
        sage: st(G, [[0],[1,2]], order=True)
        ([[0, 2, 1]], Graph on 3 vertices, 2)
        sage: st(G, [[0],[1],[2]], order=True)
        ([], Graph on 3 vertices, 1)
        sage: G.add_edge(0,1)
        sage: st(G, [range(3)], order=True)
        ([[1, 0, 2]], Graph on 3 vertices, 2)
        sage: st(G, [[0],[1,2]], order=True)
        ([], Graph on 3 vertices, 1)
        sage: st(G, [[0,1],[2]], order=True)
        ([[1, 0, 2]], Graph on 3 vertices, 2)

    The Dodecahedron has automorphism group of size 120:
        sage: G = graphs.DodecahedralGraph()
        sage: Pi = [range(20)]
        sage: st(G, Pi, order=True)[2]
        120

    The three-cube has automorphism group of size 48:
        sage: G = graphs.CubeGraph(3)
        sage: G.relabel()
        sage: Pi = [G.vertices()]
        sage: st(G, Pi, order=True)[2]
        48

    We obtain the same output using different types of Sage graphs:
        sage: G = graphs.DodecahedralGraph()
        sage: GD = DenseGraph(20)
        sage: GS = SparseGraph(20)
        sage: for i,j,_ in G.edge_iterator():
        ...    GD.add_arc(i,j); GD.add_arc(j,i)
        ...    GS.add_arc(i,j); GS.add_arc(j,i)
        sage: Pi=[range(20)]
        sage: a,b = st(G, Pi)
        sage: asp,bsp = st(GS, Pi)
        sage: ade,bde = st(GD, Pi)
        sage: bsg = Graph()
        sage: bdg = Graph()
        sage: for i in range(20):
        ...    for j in range(20):
        ...        if bsp.has_arc(i,j):
        ...            bsg.add_edge(i,j)
        ...        if bde.has_arc(i,j):
        ...            bdg.add_edge(i,j)
        sage: print a, b.graph6_string()
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]] S?[PG__OQ@?_?_?P?CO?_?AE?EC?Ac?@O
        sage: a == asp
        True
        sage: a == ade
        True
        sage: b == bsg
        True
        sage: b == bdg
        True

    Cubes!
        sage: C = graphs.CubeGraph(1)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        2
        sage: C = graphs.CubeGraph(2)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        8
        sage: C = graphs.CubeGraph(3)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        48
        sage: C = graphs.CubeGraph(4)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        384
        sage: C = graphs.CubeGraph(5)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        3840
        sage: C = graphs.CubeGraph(6)
        sage: gens, order = st(C, [C.vertices()], lab=False, order=True); order
        46080

    One can also turn off the indicator function (note- this will take longer)
        sage: D1 = DiGraph({0:[2],2:[0],1:[1]}, loops=True)
        sage: D2 = DiGraph({1:[2],2:[1],0:[0]}, loops=True)
        sage: a,b = st(D1, [D1.vertices()], dig=True, use_indicator_function=False)
        sage: c,d = st(D2, [D2.vertices()], dig=True, use_indicator_function=False)
        sage: b==d
        True

    This example is due to Chris Godsil:
        sage: HS = graphs.HoffmanSingletonGraph()
        sage: clqs = (HS.complement()).cliques()
        sage: alqs = [Set(c) for c in clqs if len(c) == 15]
        sage: Y = Graph([alqs, lambda s,t: len(s.intersection(t))==0])
        sage: Y0,Y1 = Y.connected_components_subgraphs()
        sage: st(Y0, [Y0.vertices()])[1] == st(Y1, [Y1.vertices()])[1]
        True
        sage: st(Y0, [Y0.vertices()])[1] == st(HS, [HS.vertices()])[1]
        True
        sage: st(HS, [HS.vertices()])[1] == st(Y1, [Y1.vertices()])[1]
        True

    Certain border cases need to be tested as well:
        sage: G = Graph('Fll^G')
        sage: a,b,c = st(G, [range(G.num_verts())], order=True); b
        Graph on 7 vertices
        sage: c
        48
        sage: G = Graph(21)
        sage: st(G, [range(G.num_verts())], order=True)[2] == factorial(21)
        True

        sage: G = Graph('^????????????????????{??N??@w??FaGa?PCO@CP?AGa?_QO?Q@G?CcA??cc????Bo????{????F_')
        sage: perm = {3:15, 15:3}
        sage: H = G.relabel(perm, inplace=False)
        sage: st(G, [range(G.num_verts())])[1] == st(H, [range(H.num_verts())])[1]
        True

        sage: from sage.graphs.graph import graph_isom_equivalent_non_multi_graph
        sage: G = Graph(multiedges=True, sparse=True)
        sage: G.add_edge(('a', 'b'))
        sage: G.add_edge(('a', 'b'))
        sage: G.add_edge(('a', 'b'))
        sage: G, Pi = graph_isom_equivalent_non_multi_graph(G, [['a','b']])
        sage: s,b = st(G, Pi, lab=False, dict_rep=True)
        sage: sorted(b.items())
        [(('o', 'a'), 5), (('o', 'b'), 1), (('x', 0), 2), (('x', 1), 3), (('x', 2), 4)]

        sage: st(Graph(':Dkw'), [range(5)], lab=False, dig=True)
        [[4, 1, 2, 3, 0], [0, 2, 1, 3, 4]]

    """
    cdef CGraph G
    cdef int i, j, n
    cdef Integer I
    cdef aut_gp_and_can_lab_return *output
    cdef int **part
    from sage.graphs.graph import GenericGraph, Graph, DiGraph
    if isinstance(G_in, GenericGraph):
        n = G_in.num_verts()
        G_in = G_in.copy()
        if G_in.vertices() != range(n):
            to = G_in.relabel(return_map=True)
            frm = {}
            for v in to.iterkeys():
                frm[to[v]] = v
            partition = [[to[v] for v in cell] for cell in partition]
        else:
            to = dict(enumerate(range(n)))
            frm = to
        if sparse:
            G = SparseGraph(n)
        else:
            G = DenseGraph(n)
        if G_in.is_directed():
            for i,j in G_in.edge_iterator(labels=False):
                G.add_arc(i,j)
        else:
            for i,j in G_in.edge_iterator(labels=False):
                G.add_arc(i,j)
                G.add_arc(j,i)
    elif isinstance(G_in, CGraph):
        G = <CGraph> G_in
        n = G.num_verts
        to = {}
        for a in G.verts(): to[a]=a
        frm = to
    else:
        raise TypeError("G must be a Sage graph.")

    part = <int **> sage_malloc((len(partition)+1) * sizeof(int *))
    if part is NULL:
        raise MemoryError
    for i from 0 <= i < len(partition):
        part[i] = <int *> sage_malloc((len(partition[i])+1) * sizeof(int))
        if part[i] is NULL:
            for j from 0 <= j < i:
                sage_free(part[j])
            sage_free(part)
            raise MemoryError
        for j from 0 <= j < len(partition[i]):
            part[i][j] = partition[i][j]
        part[i][len(partition[i])] = -1
    part[len(partition)] = NULL

    cdef GraphStruct GS = GraphStruct()
    GS.G = G
    GS.directed = 1 if dig else 0
    GS.use_indicator = 1 if use_indicator_function else 0
    GS.scratch = <int *> sage_malloc( (3*G.num_verts + 1) * sizeof(int) )
    if GS.scratch is NULL:
        for j from 0 <= j < len(partition):
            sage_free(part[j])
        sage_free(part)
        raise MemoryError

    output = get_aut_gp_and_can_lab(GS, part, G.num_verts, &all_children_are_equivalent, &refine_by_degree, &compare_graphs, lab, base, order)
    sage_free( GS.scratch )
    # prepare output
    list_of_gens = []
    for i from 0 <= i < output.num_gens:
        list_of_gens.append([output.generators[j+i*G.num_verts] for j from 0 <= j < G.num_verts])
    return_tuple = [list_of_gens]
    if dict_rep:
        ddd = {}
        for v in frm.iterkeys():
            ddd[frm[v]] = v if v != 0 else n
        return_tuple.append(ddd)
    if lab:
        if isinstance(G_in, GenericGraph):
            G_C = G_in.copy()
            G_C.relabel([output.relabeling[i] for i from 0 <= i < n])
        else:
            if isinstance(G, SparseGraph):
                G_C = SparseGraph(n)
            else:
                G_C = DenseGraph(n)
            for i from 0 <= i < n:
                for j in G.out_neighbors(i):
                    G_C.add_arc(output.relabeling[i],output.relabeling[j])
        return_tuple.append(G_C)
    if certify:
        dd = {}
        for i from 0 <= i < G.num_verts:
            dd[frm[i]] = output.relabeling[i]
        return_tuple.append(dd)
    if base:
        return_tuple.append([output.base[i] for i from 0 <= i < output.base_size])
    if order:
        I = Integer()
        mpz_set(I.value, output.order)
        return_tuple.append(I)
    for i from 0 <= i < len(partition):
        sage_free(part[i])
    sage_free(part)
    mpz_clear(output.order)
    sage_free(output.generators)
    if base:
        sage_free(output.base)
    if lab:
        sage_free(output.relabeling)
    sage_free(output)
    if len(return_tuple) == 1:
        return return_tuple[0]
    else:
        return tuple(return_tuple)

cdef int refine_by_degree(PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len):
    r"""
    Refines the input partition by checking degrees of vertices to the given
    cells.

    INPUT:
    PS -- a partition stack, whose finest partition is the partition to be
        refined.
    S -- a graph struct object, which contains scratch space, the graph in
        question, and some flags.
    cells_to_refine_by -- a list of pointers to cells to check degrees against
        in refining the other cells (updated in place). Must be allocated to
        length at least the degree of PS, since the array may grow
    ctrb_len -- how many cells in cells_to_refine_by

    OUTPUT:

    An integer invariant under the orbits of $S_n$.  That is, if $\gamma$ is a
    permutation of the vertices, then
    $$ I(G, PS, cells_to_refine_by) = I( \gamma(G), \gamma(PS), \gamma(cells_to_refine_by) ) .$$

    """
    cdef GraphStruct GS = <GraphStruct> S
    cdef CGraph G = GS.G
    cdef int current_cell_against = 0
    cdef int current_cell, i, r
    cdef int first_largest_subcell
    cdef int invariant = 1
    cdef int max_degree
    cdef int *degrees = GS.scratch # length 3n+1
    cdef bint necessary_to_split_cell
    cdef int against_index
    while not PS_is_discrete(PS) and current_cell_against < ctrb_len:
        invariant += 1
        current_cell = 0
        while current_cell < PS.degree:
            invariant += 50
            i = current_cell
            necessary_to_split_cell = 0
            max_degree = 0
            while 1:
                degrees[i-current_cell] = degree(PS, G, i, cells_to_refine_by[current_cell_against], 0)
                if degrees[i-current_cell] != degrees[0]:
                    necessary_to_split_cell = 1
                if degrees[i-current_cell] > max_degree:
                    max_degree = degrees[i-current_cell]
                i += 1
                if PS.levels[i-1] <= PS.depth:
                    break
            # now, i points to the next cell (before refinement)
            if necessary_to_split_cell:
                invariant += 10
                first_largest_subcell = sort_by_function(PS, current_cell, degrees)
                invariant += first_largest_subcell + max_degree
                against_index = current_cell_against
                while against_index < ctrb_len:
                    if cells_to_refine_by[against_index] == current_cell:
                        cells_to_refine_by[against_index] = first_largest_subcell
                        break
                    against_index += 1
                r = current_cell
                while 1:
                    if r == current_cell or PS.levels[r-1] == PS.depth:
                        if r != first_largest_subcell:
                            cells_to_refine_by[ctrb_len] = r
                            ctrb_len += 1
                    r += 1
                    if r >= i:
                        break
                invariant += (i - current_cell)
            current_cell = i
        if GS.directed:
            # if we are looking at a digraph, also compute
            # the reverse degrees and sort by them
            current_cell = 0
            while current_cell < PS.degree: # current_cell is still a valid cell
                invariant += 20
                i = current_cell
                necessary_to_split_cell = 0
                max_degree = 0
                while 1:
                    degrees[i-current_cell] = degree(PS, G, i, cells_to_refine_by[current_cell_against], 1)
                    if degrees[i-current_cell] != degrees[0]:
                        necessary_to_split_cell = 1
                    if degrees[i-current_cell] > max_degree:
                        max_degree = degrees[i-current_cell]
                    i += 1
                    if PS.levels[i-1] <= PS.depth:
                        break
                # now, i points to the next cell (before refinement)
                if necessary_to_split_cell:
                    invariant += 7
                    first_largest_subcell = sort_by_function(PS, current_cell, degrees)
                    invariant += first_largest_subcell + max_degree
                    against_index = current_cell_against
                    while against_index < ctrb_len:
                        if cells_to_refine_by[against_index] == current_cell:
                            cells_to_refine_by[against_index] = first_largest_subcell
                            break
                        against_index += 1
                    against_index = ctrb_len
                    r = current_cell
                    while 1:
                        if r == current_cell or PS.levels[r-1] == PS.depth:
                            if r != first_largest_subcell:
                                cells_to_refine_by[against_index] = r
                                against_index += 1
                                ctrb_len += 1
                        r += 1
                        if r >= i:
                            break
                    invariant += (i - current_cell)
                current_cell = i
        current_cell_against += 1
    if GS.use_indicator:
        return invariant
    else:
        return 0

cdef int compare_graphs(int *gamma_1, int *gamma_2, object S1, object S2):
    r"""
    Compare gamma_1(S1) and gamma_2(S2).

    Return return -1 if gamma_1(S1) < gamma_2(S2), 0 if gamma_1(S1) ==
    gamma_2(S2), 1 if gamma_1(S1) > gamma_2(S2).  (Just like the python
    \code{cmp}) function.

    INPUT:
    gamma_1, gamma_2 -- list permutations (inverse)
    S1, S2 -- graph struct objects

    """
    cdef int i, j
    cdef GraphStruct GS1 = <GraphStruct> S1
    cdef GraphStruct GS2 = <GraphStruct> S2
    cdef CGraph G1 = GS1.G
    cdef CGraph G2 = GS2.G
    for i from 0 <= i < G1.num_verts:
        for j from 0 <= j < G1.num_verts:
            if G1.has_arc_unsafe(gamma_1[i], gamma_1[j]):
                if not G2.has_arc_unsafe(gamma_2[i], gamma_2[j]):
                    return 1
            elif G2.has_arc_unsafe(gamma_2[i], gamma_2[j]):
                return -1
    return 0

cdef bint all_children_are_equivalent(PartitionStack *PS, object S):
    """
    Return True if every refinement of the current partition results in the
    same structure.

    WARNING:
    Converse does not hold in general!  See Lemma 2.25 of [1] for details.

    INPUT:
    PS -- the partition stack to be checked
    S -- a graph struct object
    """
    cdef GraphStruct GS = <GraphStruct> S
    if GS.directed:
        return 0
    cdef CGraph G = GS.G
    cdef int i, n = PS.degree
    cdef bint in_cell = 0
    cdef int nontrivial_cells = 0
    cdef int total_cells = PS_num_cells(PS)
    if n <= total_cells + 4:
        return 1
    for i from 0 <= i < n-1:
        if PS.levels[i] <= PS.depth:
            if in_cell:
                nontrivial_cells += 1
            in_cell = 0
        else:
            in_cell = 1
    if in_cell:
        nontrivial_cells += 1
    if n == total_cells + nontrivial_cells:
        return 1
    if n == total_cells + nontrivial_cells + 1:
        return 1
    return 0

cdef inline int degree(PartitionStack *PS, CGraph G, int entry, int cell_index, bint reverse):
    """
    Returns the number of edges from the vertex corresponding to entry to
    vertices in the cell corresponding to cell_index.

    INPUT:
    PS -- the partition stack to be checked
    S -- a graph struct object
    entry -- the position of the vertex in question in the entries of PS
    cell_index -- the starting position of the cell in question in the entries
        of PS
    reverse -- whether to check for arcs in the other direction
    """
    cdef int num_arcs = 0
    entry = PS.entries[entry]
    if not reverse:
        while 1:
            if G.has_arc_unsafe(PS.entries[cell_index], entry):
                num_arcs += 1
            if PS.levels[cell_index] > PS.depth:
                cell_index += 1
            else:
                break
    else:
        while 1:
            if G.has_arc_unsafe(entry, PS.entries[cell_index]):
                num_arcs += 1
            if PS.levels[cell_index] > PS.depth:
                cell_index += 1
            else:
                break
    return num_arcs

cdef inline int sort_by_function(PartitionStack *PS, int start, int *degrees):
    """
    A simple counting sort, given the degrees of vertices to a certain cell.

    INPUT:
    PS -- the partition stack to be checked
    start -- beginning index of the cell to be sorted
    degrees -- the values to be sorted by

    """
    cdef int n = PS.degree
    cdef int i, j, max, max_location
    cdef int *counts = degrees + n, *output = degrees + 2*n + 1
    for i from 0 <= i <= n:
        counts[i] = 0
    i = 0
    while PS.levels[i+start] > PS.depth:
        counts[degrees[i]] += 1
        i += 1
    counts[degrees[i]] += 1
    # i+start is the right endpoint of the cell now
    max = counts[0]
    max_location = 0
    for j from 0 < j <= n:
        if counts[j] > max:
            max = counts[j]
            max_location = j
        counts[j] += counts[j - 1]
    for j from i >= j >= 0:
        counts[degrees[j]] -= 1
        output[counts[degrees[j]]] = PS.entries[start+j]
    max_location = counts[max_location]+start
    for j from 0 <= j <= i:
        PS.entries[start+j] = output[j]
    j = 1
    while j <= n and counts[j] <= i:
        if counts[j] > 0:
            PS.levels[start + counts[j] - 1] = PS.depth
        PS_move_min_to_front(PS, start + counts[j-1], start + counts[j] - 1)
        j += 1
    return max_location


def all_labeled_graphs(n):
    """
    Returns all labeled graphs on n vertices {0,1,...,n-1}. Used in
    classifying isomorphism types (naive approach), and more importantly
    in benchmarking the search algorithm.

    EXAMPLE::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import all_labeled_graphs
        sage: st = sage.groups.perm_gps.partn_ref.refinement_graphs.search_tree
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in [1..5]:
        ...    Glist[n] = all_labeled_graphs(n)
        ...    Giso[n] = []
        ...    for g in Glist[n]:
        ...        a, b = st(g, [range(n)])
        ...        inn = False
        ...        for gi in Giso[n]:
        ...            if b == gi:
        ...                inn = True
        ...        if not inn:
        ...            Giso[n].append(b)
        sage: for n in Giso:
        ...    print n, len(Giso[n])
        1 1
        2 2
        3 4
        4 11
        5 34

    """
    from sage.graphs.graph import Graph
    TE = []
    for i in range(n):
        for j in range(i):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = Graph(n)
        b = Integer(i).binary()
        b = '0'*(m-len(b)) + b
        for i in range(m):
            if int(b[i]):
                G.add_edge(TE[i])
        Glist.append(G)
    return Glist


def random_tests(t=10.0, n_max=60, perms_per_graph=10):
    """
    Tests to make sure that C(gamma(G)) == C(G) for random permutations gamma
    and random graphs G, and that isomorphic returns an isomorphism.

    INPUT:
    t -- run tests for approximately this many seconds
    n_max -- test graphs with at most this many vertices
    perms_per_graph -- test each graph with this many random permutations

    DISCUSSION:

    Until t seconds have elapsed, this code generates a random graph G on at
    most n_max vertices.  The density of edges is chosen randomly between 0
    and 1.

    For each graph G generated, we uniformly generate perms_per_graph random
    permutations and verify that the canonical labels of G and the image of G
    under the generated permutation are equal, and that the isomorphic function
    returns an isomorphism.

    TESTS:
        sage: import sage.groups.perm_gps.partn_ref.refinement_graphs
        sage: sage.groups.perm_gps.partn_ref.refinement_graphs.random_tests()
        All passed: ... random tests on ... graphs.
        sage: sage.groups.perm_gps.partn_ref.refinement_graphs.random_tests(180.0, 200, 30) # long time
        All passed: ... random tests on ... graphs.

    """
    from sage.misc.misc import walltime
    from sage.misc.prandom import random, randint
    from sage.graphs.graph_generators import GraphGenerators, DiGraphGenerators
    from sage.combinat.permutation import Permutations
    cdef int i, j, num_tests = 0, num_graphs = 0
    GG = GraphGenerators()
    DGG = DiGraphGenerators()
    t_0 = walltime()
    while walltime(t_0) < t:
        p = random()
        n = randint(1, n_max)
        S = Permutations(n)

        G = GG.RandomGNP(n, p)
        H = G.copy()
        for i from 0 <= i < perms_per_graph:
            G = H.copy()
            G1 = search_tree(G, [G.vertices()])[1]
            perm = list(S.random_element())
            perm = [perm[j]-1 for j from 0 <= j < n]
            G.relabel(perm)
            G2 = search_tree(G, [G.vertices()])[1]
            if G1 != G2:
                print "search_tree FAILURE: graph6-"
                print H.graph6_string()
                print perm
                return
            isom = isomorphic(G, H, [range(n)], range(n), 0, 1)
            if not isom or G.relabel(isom, inplace=False) != H:
                print "isom FAILURE: graph6-"
                print H.graph6_string()
                print perm
                return

        D = DGG.RandomDirectedGNP(n, p)
        D.allow_loops(True)
        for i from 0 <= i < n:
            if random() < p:
                D.add_edge(i,i)
        E = D.copy()
        for i from 0 <= i < perms_per_graph:
            D = E.copy()
            D1 = search_tree(D, [D.vertices()], dig=True)[1]
            perm = list(S.random_element())
            perm = [perm[j]-1 for j from 0 <= j < n]
            D.relabel(perm)
            D2 = search_tree(D, [D.vertices()], dig=True)[1]
            if D1 != D2:
                print "search_tree FAILURE: dig6-"
                print E.dig6_string()
                print perm
                return
            isom = isomorphic(D, E, [range(n)], range(n), 1, 1)
            if not isom or D.relabel(isom, inplace=False) != E:
                print "isom FAILURE: dig6-"
                print E.dig6_string()
                print perm
                print isom
                return
        num_tests += 4*perms_per_graph
        num_graphs += 2
    print "All passed: %d random tests on %d graphs."%(num_tests, num_graphs)

def orbit_partition(gamma, list_perm=False):
    r"""
    Assuming that G is a graph on vertices 0,1,...,n-1, and gamma is an
    element of SymmetricGroup(n), returns the partition of the vertex
    set determined by the orbits of gamma, considered as action on the
    set 1,2,...,n where we take 0 = n. In other words, returns the
    partition determined by a cyclic representation of gamma.

    INPUT:


    -  ``list_perm`` - if True, assumes
       ``gamma`` is a list representing the map
       `i \mapsto ``gamma``[i]`.


    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import orbit_partition
        sage: G = graphs.PetersenGraph()
        sage: S = SymmetricGroup(10)
        sage: gamma = S('(10,1,2,3,4)(5,6,7)(8,9)')
        sage: orbit_partition(gamma)
        [[1, 2, 3, 4, 0], [5, 6, 7], [8, 9]]
        sage: gamma = S('(10,5)(1,6)(2,7)(3,8)(4,9)')
        sage: orbit_partition(gamma)
        [[1, 6], [2, 7], [3, 8], [4, 9], [5, 0]]
    """
    if list_perm:
        n = len(gamma)
        seen = [1] + [0]*(n-1)
        i = 0
        p = 0
        partition = [[0]]
        while sum(seen) < n:
            if gamma[i] != partition[p][0]:
                partition[p].append(gamma[i])
                i = gamma[i]
                seen[i] = 1
            else:
                for j in range(n):
                    if seen[j]==0:
                        i = j
                        break
                partition.append([i])
                p += 1
                seen[i] = 1
        return partition
    else:
        n = len(gamma.list())
        l = []
        for i in range(1,n+1):
            orb = gamma.orbit(i)
            if orb not in l: l.append(orb)
        for i in l:
            for j in range(len(i)):
                if i[j] == n:
                    i[j] = 0
        return l

def perm_group_elt(lperm):
    """
    Given a list permutation of the set 0, 1, ..., n-1, returns the
    corresponding PermutationGroupElement where we take 0 = n.

    EXAMPLE::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import perm_group_elt
        sage: perm_group_elt([0,2,1])
        (1,2)
        sage: perm_group_elt([1,2,0])
        (1,2,3)
    """
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    n = len(lperm)
    S = SymmetricGroup(n)
    Part = orbit_partition(lperm, list_perm=True)
    gens = []
    for z in Part:
        if len(z) > 1:
            if 0 in z:
                zed = z.index(0)
                generator = z[:zed] + [n] + z[zed+1:]
                gens.append(tuple(generator))
            else:
                gens.append(tuple(z))
    E = S(gens)
    return E

def coarsest_equitable_refinement(CGraph G, list partition, bint directed):
    """
    Returns the coarsest equitable refinement of ``partition`` for ``G``.

    This is a helper function for the graph function of the same name.

    DOCTEST (More thorough testing in ``sage/graphs/graph.py``)::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import coarsest_equitable_refinement
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: coarsest_equitable_refinement(SparseGraph(7), [[0], [1,2,3,4], [5,6]], 0)
        [[0], [1, 2, 3, 4], [5, 6]]

    """
    cdef int i, j = 0, k = 0, n = G.num_verts

    # set up partition stack and graph struct
    cdef PartitionStack *nu = PS_new(n, 0)
    for cell in partition:
        for i in cell:
            nu.entries[j] = i
            nu.levels[j] = n
            j += 1
        nu.levels[j-1] = 0
        PS_move_min_to_front(nu, k, j-1)
        k = j

    cdef GraphStruct GS = GraphStruct()
    GS.G = G
    GS.scratch = <int *> sage_malloc((3*n+1) * sizeof(int))
    if GS.scratch is NULL:
        PS_clear(nu)
        raise MemoryError
    GS.directed = directed
    GS.use_indicator = 0

    # set up cells to refine by
    cdef int num_cells = len(partition)
    cdef int *alpha = <int *>sage_malloc(n * sizeof(int))
    if alpha is NULL:
        PS_clear(nu)
        sage_free(GS.scratch)
        raise MemoryError
    j = 0
    for i from 0 <= i < num_cells:
        alpha[i] = j
        j += len(partition[i])

    # refine, and get the result
    refine_by_degree(nu, GS, alpha, num_cells)

    eq_part = []
    cell = []
    for i from 0 <= i < n:
        cell.append(nu.entries[i])
        if nu.levels[i] <= 0:
            eq_part.append(cell)
            cell = []

    PS_clear(nu)
    sage_free(GS.scratch)
    sage_free(alpha)

    return eq_part


def get_orbits(list gens, int n):
    """
    Compute orbits given a list of generators of a permutation group, in list
    format.

    This is a helper function for automorphism groups of graphs.

    DOCTEST (More thorough testing in ``sage/graphs/graph.py``)::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import get_orbits
        sage: get_orbits([[1,2,3,0,4,5], [0,1,2,3,5,4]], 6)
        [[0, 1, 2, 3], [4, 5]]

    """
    cdef int i, j
    if len(gens) == 0:
        return [[i] for i from 0 <= i < n]

    cdef OrbitPartition *OP = OP_new(n)
    cdef int *perm_ints = <int *> sage_malloc(n * sizeof(int))
    if perm_ints is NULL:
        OP_dealloc(OP)
        raise MemoryError

    for gen in gens:
        for i from 0 <= i < n:
            perm_ints[i] = gen[i]
        OP_merge_list_perm(OP, perm_ints)

    orbit_dict = {}
    for i from 0 <= i < n:
        j = OP_find(OP, i)
        if orbit_dict.has_key(j):
            orbit_dict[j].append(i)
        else:
            orbit_dict[j] = [i]

    OP_dealloc(OP)
    sage_free(perm_ints)

    return orbit_dict.values()







