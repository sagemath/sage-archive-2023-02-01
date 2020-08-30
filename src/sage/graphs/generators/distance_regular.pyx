r"""
Database of distance regular graphs

In this module we construct several distance regular graphs
and group them in a function that maps intersection arrays
to graphs.

For a survey on distance-regular graph see [BCN1989]_ or [VDKT2016]_.

EXAMPLES::

   sage: G = graphs.cocliques_HoffmannSingleton()
   sage: G.is_distance_regular()
   True

AUTHORS:

- Ivo Maffei (2020-07-28): initial version

"""

# ****************************************************************************
#       Copyright (C) 2020 Ivo Maffei <ivomaffei@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.graphs.graph import Graph
from sage.libs.gap.libgap import libgap
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import Matrix
import itertools
from cysignals.signals cimport sig_check

def cocliques_HoffmannSingleton():
    r"""
    Return the graph obtained from the cocliques of the Hoffmann-Singleton graph.

    This is a distance-regular graph with intersection array
    `[15, 14, 10, 3; 1, 5, 12, 15]`.

    EXAMPLES::

        sage: G = graphs.cocliques_HoffmannSingleton()
        sage: G.is_distance_regular(True)
        ([15, 14, 10, 3, None], [None, 1, 5, 12, 15])

    REFERENCES:

    The construction of this graph can be found in [BCN1989]_ p. 392.
    """
    from sage.graphs.graph_generators import GraphGenerators

    D = GraphGenerators.HoffmanSingletonGraph()
    DC = D.complement()

    cocliques = [frozenset(c) for c in DC.cliques_maximum()]  # 100 of this

    edges = []
    for c1, c2 in itertools.combinations(cocliques, 2):
        if len(c1.intersection(c2)) == 8:
            edges.append((c1, c2))

    G = Graph(edges, format="list_of_edges")
    return G

def locally_GQ42_distance_transitive_graph():
    r"""
    Return the unique amply regular graph with `\mu = 6` which is locally
    a generalised quadrangle.

    This graph is distance-regular with intersection array
    `[45, 32, 12, 1; 1, 6, 32, 45]`.

    This graph is also distance-transitive.

    EXAMPLES::

        sage: G = graphs.locally_GQ42_distance_transitive_graph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([45, 32, 12, 1, None], [None, 1, 6, 32, 45])

    REFERENCES:

    A description of this graph can be found in [BCN1989]_ p.399.
    This construction is due to Dima Pasechnik.
    """
    H = libgap.AtlasGroup("3^2.U4(3).D8", libgap.NrMovedPoints, 756)
    Ns = H.NormalSubgroups()
    for N in Ns:
        if len(N.GeneratorsSmallest()) == 7:  # there is only one
            break

    G = Graph(libgap.Orbit(N, [1, 9], libgap.OnSets), format='list_of_edges')
    G.name("locally GQ(4,2) distance transitive graph")
    return G


def ConwaySmith_for_3S7():
    r"""
    Return the Conway-Smith graph related to `3 Sym(7)`.

    This is a distance-regular graph with intersection array
    `[10, 6, 4, 1; 1, 2, 6, 10]`.

    EXAMPLES::

        sage: G = graphs.ConwaySmith_for_3S7()
        sage: G.is_distance_regular(True)
        ([10, 6, 4, 1, None], [None, 1, 2, 6, 10])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 399.
    """
    from sage.rings.number_field.number_field import CyclotomicField

    F = CyclotomicField(3)
    w = F.gen()

    V= VectorSpace(GF(4), 6)
    z2 = GF(4)('z2') # GF(4) = {0, 1, z2, z2+1}

    W = V.span([(0,0,1,1,1,1), (0,1,0,1,z2,z2+1), (1,0,0,1,z2+1,z2)])
    # we only need the 45 vectors with 2 zero entries
    # we also embed everything into CC

    K = []
    for v in W:
        # check zero entries
        zeros = 0
        for x in v:
            if x.is_zero():
                zeros += 1

        if zeros == 2:
            # send to F and in K
            # z2 -> w
            # z2+1 -> w^2
            vv = []  # new vector
            for x in v:
                if x == z2:
                    vv.append(w)
                elif x == z2 + 1:
                    vv.append(w**2)
                else:
                    vv.append(int(x))

            # now vv is the new vector in F
            vv = vector(F, vv)
            K.append(vv)

    # we need to add other vectors
    for i in range(6):
        # create e_i
        ei = [0, 0, 0, 0, 0, 0]
        ei[i] = 1
        ei = vector(F, ei)

        K.append(2 * ei)
        K.append(2 * w * ei)
        K.append(2 * w**2 * ei)
    # now K is all the 63 vertices

    for v in K:
        v.set_immutable()

    def has_edge(u, v):
        return sum(u[i].conjugate() * v[i] for i in range(6)) == 2

    G = Graph()
    for Ki, Kj in itertools.combinations(K, 2):
        if has_edge(Ki, Kj):
            G.add_edge((Ki, Kj))

    G.name("Conway-Smith graph for 3S7")
    return G

def graph_3O73():
    r"""
    Return the graph related to the group `3 O(7,3)`.

    This graph is distance-regular with intersection array
    `[117, 80, 24, 1; 1, 12, 80, 117]`.

    The graph is also distance transitive with "3.O(7,3)" as automorphism
    group

    EXAMPLES::

        sage: G = graphs.graph_3O73()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([117, 80, 24, 1, None], [None, 1, 12, 80, 117])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 400.
    """
    group = libgap.AtlasGroup("3.O7(3)", libgap.NrMovedPoints, 1134)
    G = Graph(libgap.Orbit(group, [1, 3], libgap.OnSets), format='list_of_edges')
    G.name("Distance transitive graph with automorphism group 3.O_7(3)")
    return G

def FosterGraph3S6():
    r"""
    Return the Foster graph for `3.Sym(6)`.

    This graph is distance-regular with intersection array
    `[6, 4, 2, 1; 1, 1, 4, 6]`.

    The graph is also distance transitive.

    EXAMPLES::

        sage: G = graphs.FosterGraph3S6()
        sage: G.is_distance_regular(True)
        ([6, 4, 2, 1, None], [None, 1, 1, 4, 6])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 397.
    """

    a = libgap.eval(("(2,6)(3,5)(4,11)(7,17)(8,16)(9,14)(13,22)(15,25)"
                    "(18,29)(19,28)(20,21)(24,30)(26,35)(27,33)(31,39)"
                     "(34,38)(36,43)(37,40)(42,44)"))
    b = libgap.eval(("(1,2,7,12,4)(3,8,18,20,10)(5,9,19,21,11)(6,13,17,26,15)"
                     "(14,23,28,31,24)(16,22,29,36,27)(25,32,35,42,34)"
                     "(30,37,39,44,38)(33,40,43,45,41)"))

    group = libgap.Group(a,b)

    G = Graph(group.Orbit([1, 7], libgap.OnSets), format='list_of_edges')
    G.name("Foster graph for 3.Sym(6) graph")
    return G

def J2Graph():
    r"""
    Return the distance-transitive graph with automorphism group `J_2`.

    EXAMPLES::

        sage: G = graphs.J2Graph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True) # optional - internet gap_packages
        ([10, 8, 8, 2, None], [None, 1, 1, 4, 5])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 408.
    """
    group = libgap.AtlasGroup("J2", libgap.NrMovedPoints, 315)
    G = Graph(group.Orbit([1, 9], libgap.OnSets), format='list_of_edges')
    G.name("J_2 graph")
    return G

def IvanovIvanovFaradjevGraph():
    r"""
    Return the IvanovIvanovFaradjev graph.

    The graph is distance-transitive with automorphism group `3.M_{22}`.

    EXAMPLES::

        sage: G = graphs.IvanovIvanovFaradjevGraph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([7, 6, 4, 4, 4, 1, 1, 1, None], [None, 1, 1, 1, 2, 4, 4, 6, 7])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 369.
    """

    group = libgap.AtlasGroup("3.M22", libgap.NrMovedPoints, 990)
    graph = Graph(group.Orbit([1, 22], libgap.OnSets), format='list_of_edges')

    graph.name("Ivanov-Ivanov-Faradjev Graph")
    return graph

def LargeWittGraph():
    r"""
    Return the large Witt graph.

    This is a distance-regular graph with intersection array
    `[30,28,24;1,3,15]`.

    EXAMPLES::

        sage: g = graphs.LargeWittGraph()
        sage: g.is_distance_regular(True)
        ([30, 28, 24, None], [None, 1, 3, 15])

    REFERENCES:

    A description of this graph can be found in
    [BCN1989]_ p. 366.
    This construction is taken from
    http://mathworld.wolfram.com/LargeWittGraph.html
    """
    from sage.coding import codes_catalog as codes
    import itertools

    C = codes.GolayCode(GF(2), extended=True)
    vertices = [c for c in C if c.hamming_weight() == 8]

    edges = []
    for v, w in itertools.combinations(vertices, 2):
        if not set(v.support()).intersection(w.support()):
            edges.append((v, w))

    W = Graph(edges, format='list_of_edges')
    W.name("Large Witt graph")
    return W

def TruncatedWittGraph():
    r"""
    Return the truncated Witt graph.

    This builds the large Witt graph, then removes
    all vertices whose codeword start with a 1.

    The graph is distance-regular with intersection array
    `[15,14,12;1,1,9]`.

    EXAMPLES::

         sage: G = graphs.TruncatedWittGraph()
         sage: G.is_distance_regular(True)
         ([15, 14, 12, None], [None, 1, 1, 9])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 367.
    """
    # get large witt graph and remove all vertices which start with a 1
    G = LargeWittGraph()
    G.delete_vertices(filter(lambda x : x[0] == 1, G.vertices()))

    G.name("Truncated Witt graph")
    return G

def DoublyTruncatedWittGraph():
    r"""
    Return the doubly truncated Witt graph.

    This builds the truncated Witt graph, then removes
    all vertices whose codeword start with a 1.

    The graph is distance-regular with intersection array
    `[7,6,4,4;1,1,1,6]`.

    EXAMPLES::

         sage: G = graphs.DoublyTruncatedWittGraph()
         sage: G.is_distance_regular(True)
         ([7, 6, 4, 4, None], [None, 1, 1, 1, 6])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 368.
    """

    G = TruncatedWittGraph()
    G.delete_vertices(filter(lambda x : x[1] == 1, G.vertices()))

    G.name("Doubly Truncated Witt graph")
    return G

def distance_3_doubly_truncated_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[9, 8, 6, 3; 1, 1, 3, 8]`.

    EXAMPLES::

        sage: G = graphs.distance_3_doubly_truncated_Golay_code_graph()
        sage: G.is_distance_regular(True)
        ([9, 8, 6, 3, None], [None, 1, 1, 3, 8])

    ALGORITHM:

    Compute the binary Golay code and truncate it twice. Compute its coset graph.
    Take a vertex and compute the set of vertices at distance 3
    from the vertex choosen. This set constitutes the set of vertices of our
    distance-regular graph. Moreover we have an edge `(u,v)` if the coset graph
    contains such edge.

    REFERENCES:

    Description and construction of this graph are taken from [BCN1989]_ p. 364.
    """
    from sage.coding import codes_catalog as codes

    G = codes.GolayCode(GF(2),extended=False).punctured([0,1]).cosetGraph()
    v = G.vertices(sort=False)[0]
    it = G.breadth_first_search(v, distance=3, report_distance=True)
    vertices = [w for (w,d) in it if d == 3]

    edges =[(a ,b) for a, b in itertools.combinations(vertices, 2)
            if G.has_edge((a, b))]

    H = Graph(edges, format='list_of_edges')
    return H

def shortened_00_11_binary_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[21, 20, 16, 6, 2, 1; 1, 2, 6, 16, 20, 21]`.

    EXAMPLES::

        sage: G = graphs.shortened_00_11_binary_Golay_code_graph() # long time (25 s)
        sage: G.is_distance_regular(True) # long time
        ([21, 20, 16, 6, 2, 1, None], [None, 1, 2, 6, 16, 20, 21])

    ALGORITHM:

    Compute the binary Golay code. Compute the subcode whose codewords start
    with 00 or 11. Remove the first two entries from all codewords of the newly
    found linear code and compute its coset graph.

    REFERENCES:

    Description and construction of this graph can be found in [BCN1989]_ p. 365.
    """
    from sage.coding import codes_catalog as codes
    from sage.coding.linear_code import LinearCode

    code = codes.GolayCode(GF(2), False)
    C_basis = code.basis()

    # Now special shortening
    v = C_basis[0] + C_basis[1] # v has 11 at the start
    C_basis = C_basis[2:]
    C_basis.append(v)
    C_basis = list(map(lambda x: x[2:], C_basis))

    code = LinearCode(Matrix(GF(2), C_basis))

    G = code.cosetGraph()
    G.name("Shortened 00 11 binary Golay code")
    return G

def shortened_000_111_extended_binary_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[21, 20, 16, 9, 2, 1; 1, 2, 3, 16, 20, 21]`.

    EXAMPLES::

        sage: G = graphs.shortened_000_111_extended_binary_Golay_code_graph() # long time (2 min)
        sage: G.is_distance_regular(True)  # long time
        ([21, 20, 16, 9, 2, 1, None], [None, 1, 2, 3, 16, 20, 21])

    ALGORITHM:

    Compute the extended binary Golay code. Compute its subcode whose codewords
    start with 000 or 111. Remove the first 3 entries from all the codewords
    from the new linear code and compute its coset graph.

    REFERENCES:

    Description and construction of this graph can be found in [BCN1989]_ p. 365.
    """
    from sage.coding import codes_catalog as codes
    from sage.coding.linear_code import LinearCode

    code = codes.GolayCode(GF(2))
    C_basis = code.basis()

    # now special shortening
    v = C_basis[0] + C_basis[1] + C_basis[2] # v has 111 at the start
    C_basis = C_basis[3:]
    C_basis.append(v)
    C_basis = list(map(lambda x: x[3:], C_basis))

    code = LinearCode(Matrix(GF(2), C_basis))

    G = code.cosetGraph()
    G.name("Shortened 000 111 extended binary Golay code")
    return G

def LintSchrijverGraph():
    r"""
    Return the van Lint-Schrijver graph.

    The graph is distance-regular with intersection array
    `[6, 5, 5, 4; 1, 1, 2, 6]`.

    EXAMPLES::

         sage: G = graphs.LintSchrijverGraph()
         sage: G.is_distance_regular(True)
         ([6, 5, 5, 4, None], [None, 1, 1, 2, 6])

    REFERENCES:

    For a description of this graph see [BCN1989]_ p. 373.
    """
    from sage.coding.linear_code import LinearCode

    one = vector(GF(3), [1, 1, 1, 1, 1, 1])
    G = LinearCode(Matrix(GF(3), one)).cosetGraph()

    vertices = [v for v in G.vertices() if v.dot_product(one) in {1, 2}]
    edges = [(v, w) for v, w in itertools.combinations(vertices, 2)
             if G.has_edge((v, w))]

    H = Graph(edges, format='list_of_edges')
    H.name("Linst-Schrijver graph")
    return H

def LeonardGraph():
    r"""
    Return the Leonard graph.

    The graph is distance-regular with intersection array
    `[12, 11, 10, 7; 1, 2, 5, 12]`.

    EXAMPLES::

         sage: G = graphs.LeonardGraph()
         sage: G.is_distance_regular(True)
         ([12, 11, 10, 7, None], [None, 1, 2, 5, 12])

    REFERENCES:

    For a description of this graph see [BCN1989]_ p. 371.
    """
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix

    M = hadamard_matrix(12)
    edges = []
    for i, j, k, l in itertools.product(range(12), repeat=4):
        if i == k or j == l: continue
        if M[i, j] * M[i, l] * M[k, j] * M[k, l] == -1:
            edges.append(((i, j), (k, l)))

    D = Graph(edges, format="list_of_edges")
    blocks = [frozenset(cl) for cl in D.cliques_maximum()]

    edges = [(p, b) for b in blocks for p in b]
    G = Graph(edges, format="list_of_edges")
    return G

def UstimenkoGraph(const int m, const int q):
    r"""
    Return the Ustimenko graph with parameters `(m, q)`.

    This is the distance 1 or 2 graph of the dual polar graph `C_{m-1}(q)`.
    The graph is distance-regular with classical with parameters
    `(d,q^2, qbinom(3,1,q) -1, qbinom(m+1,1,q) -1)`

    INPUT:

    - ``m, q`` -- integers; ``q`` must be a prime power and ``m > 1``.

    EXAMPLES::

        sage: G = graphs.UstimenkoGraph(4, 2)
        sage: G.is_distance_regular(True)
        ([70, 32, None], [None, 1, 35])

    REFERENCES:

    See [BCN1989]_ p. 279 or [VDKT2016]_ p. 22.

    TESTS::

        sage: G = graphs.UstimenkoGraph(5, 2)
        sage: G.order()
        2295
        sage: G.is_distance_regular(True)
        ([310, 224, None], [None, 1, 35])
        sage: G = graphs.UstimenkoGraph(4,3)
        sage: G.is_distance_regular(True)
        ([390, 243, None], [None, 1, 130])
    """
    from sage.graphs.graph_generators import graphs

    G = graphs.SymplecticDualPolarGraph(2*m - 2, q)

    edgesToAdd = []
    for v in G:
        for w in G.neighbor_iterator(v):
            for u in G.neighbor_iterator(w):
                sig_check()
                if u != v and not G.has_edge(u, v):
                    # then u,v are at distance 2
                    edgesToAdd.append((u, v))

    G.add_edges(edgesToAdd)
    G.name(f"Ustimenko graph ({m}, {q})")
    return G

def BilinearFormsGraph(const int d, const int e, const int q):
    r"""
    Return a bilinear forms graph with the given parameters.

    This build a graph whose vertices are all ``d``x``e`` matrices over
    ``GF(q)``. Two vertices are adjacent if the difference of the two
    matrices has rank 1.

    The graph is distance-regular with classical parameters
    `(\min(d, e), q, q-1 , q^{\max(d, e)}-1)`.

    INPUT:

    - ``d, e`` -- integers; dimension of the matrices
    - ``q`` -- integer; a prime power

    EXAMPLES::

        sage: G = graphs.BilinearFormsGraph(3, 3, 2)
        sage: G.is_distance_regular(True)
        ([49, 36, 16, None], [None, 1, 6, 28])
        sage: G = graphs.BilinearFormsGraph(3,3,3)  # not tested (20 s)
        sage: G.order()  # not tested (due to above)
        19683
        sage: G = graphs.BilinearFormsGraph(3, 4, 2)
        sage: G.is_distance_regular(True)
        ([105, 84, 48, None], [None, 1, 6, 28])

    REFERENCES:

    See [BCN1989]_ pp. 280-282 for a rather detailed discussion, otherwise
    see [VDKT2016]_ p. 21.

    TESTS::

        sage: G = graphs.BilinearFormsGraph(2,3,2)
        sage: G.is_distance_regular(True)
        ([21, 12, None], [None, 1, 6])
        sage: H = graphs.BilinearFormsGraph(3,2,2)
        sage: H.is_isomorphic(G)
        True
        sage: G = graphs.BilinearFormsGraph(5, 1, 3)
        sage: K = graphs.CompleteGraph(G.order())
        sage: K.is_isomorphic(G)
        True
    """
    from sage.combinat.integer_vector import IntegerVectors

    Fq = GF(q)
    Fqelems = list(Fq)
    matricesOverq = IntegerVectors(k=d*e, max_part=q-1)

    rank1Matrices = []
    for u in VectorSpace(Fq, d):
        if u.is_zero() or not u[u.support()[0]].is_one():
            continue

        for v in VectorSpace(Fq, e):
            if v.is_zero():
                continue

            sig_check()
            M = [0] * (d*e)
            for row in range(d):
                for col in range(e):
                    M[e*row + col] = u[row] * v[col]

            rank1Matrices.append(M)

    edges = []
    for m1 in matricesOverq:
        m1 = tuple(map(lambda x: Fqelems[x], m1))
        for m2 in rank1Matrices:
            sig_check()
            m3 = tuple([m1[i] + m2[i] for i in range(d*e)])
            edges.append((m1, m3))

    G = Graph(edges, format='list_of_edges')
    G.name("Bilinear forms graphs over F_%d with parameters (%d, %d)"%(q, d, e))
    return G

def AlternatingFormsGraph(const int n, const int q):
    r"""
    Return the alternating forms graph with the given parameters.

    This construct a graph whose vertices are all ``n``x``n`` skew-symmetric
    matrices over ``GF(q)`` with zero diagonal. Two vertices are adjacent
    if and only if the difference of the two matrices has rank 2.

    This grap is distance-regular with classical parameters
    `(\lfloor \frac n 2 \rfloor,  q^2, q^2 - 1, q^{2 \lceil \frac n 2 \rceil -1})`.

    INPUT:

    - ``n`` -- integer
    - ``q`` -- a prime power

    EXAMPLES::

        sage: G = graphs.AlternatingFormsGraph(5,2)
        sage: G.is_distance_regular(True)
        ([155, 112, None], [None, 1, 20])

    REFERENCES:

    See [BCN1989]_ pp. 282-284 for a rather detailed discussion, otherwise
    see [VDKT2016]_ p. 22.

    TESTS::

         sage: G = graphs.AlternatingFormsGraph(6,2)  # not tested (2 min)
         sage: G.order()  # not tested (because of above)
         32768
         sage: G.is_distance_regular(True)  # not tested (33 min)
         ([651, 560, 256, None], [None, 1, 20, 336])
         sage: G = graphs.AlternatingFormsGraph(4, 3)
         sage: G.is_distance_regular(True)
         ([260, 162, None], [None, 1, 90])
    """
    # n x n zero-diagonal skew-symmetric matrix
    # can be represented by the upper triangular entries
    # there are n*(n-1) // 2 of them
    size = (n * (n-1)) // 2
    V = VectorSpace(GF(q), size)

    # construct all rank 2 matrices
    rank2Matrices = set()
    Vn = VectorSpace(GF(q), n)
    basis = set(Vn.basis())
    e = [Vn([0]*i + [1] + [0]*(n - i - 1)) for i in range(n)]
    for v in e:
        v.set_immutable()

    scalars = [x for x in GF(q) if not x.is_zero()]
    Vseen = set()
    for v in Vn:
        if v.is_zero() or not v[v.support()[0]].is_one():
            continue
        v.set_immutable()
        # remove from basis e_i s.t. (v[i-1] =) v_i != 0
        i = v.support()[0]
        Ubasis = basis.difference([e[i]])

        for u in Vn.span_of_basis(Ubasis):
            sig_check()
            if u.is_zero() or not u[u.support()[0]].is_one():
                continue
            u.set_immutable()
            if u in Vseen:
                continue

            M = []
            for row in range(n - 1):
                upperRow = [0] * (n - 1 - row)
                for col in range(row + 1, n):
                    upperRow[col - row - 1] = v[row]*u[col] - u[row]*v[col]
                M += upperRow

            for scalar in scalars:
                N = tuple(map(lambda x: scalar * x, M))
                rank2Matrices.add(N)

        Vseen.add(v)

    # now we have all matrices of rank 2
    edges = []
    for m1 in V:
        t1 = tuple(m1)
        for m2 in rank2Matrices:
            sig_check()
            t3 = tuple([t1[i] + m2[i] for i in range(size)])
            edges.append((t1, t3))

    G = Graph(edges, format='list_of_edges')
    G.name("Alternating forms graph on (F_%d)^%d"%(q, n))
    return G

def HermitianFormsGraph(const int n, const int r):
    r"""
    Return the Hermitian froms graph with the given parameters.

    We build a graph whose vertices are all ``n``x``n`` Hermitian matrices
    over ``GF(r^2)``. Two  vertices are adjacent if the difference of the two
    vertices has rank 1.

    This graph is distance-regular with classical parameters
    `(n, - r, - r - 1, - (- r)^d - 1)`.

    INPUT:

    - ``n`` -- integer
    - ``r`` -- a prime power

    EXAMPLES::

        sage: G = graphs.HermitianFormsGraph(2, 2)
        sage: G.is_distance_regular(True)
        ([5, 4, None], [None, 1, 2])
        sage: G = graphs.HermitianFormsGraph(3, 3)  # not tested (2 min)
        sage: G.order()  # not tested (bacuase of the above)
        19683

    REFERENCES:

    See [BCN1989]_ p. 285 or [VDKT2016]_ p. 22.

    TESTS::

         sage: G = graphs.HermitianFormsGraph(3, 2)
         sage: G.is_distance_regular(True)
         ([21, 20, 16, None], [None, 1, 2, 12])
         sage: G = graphs.HermitianFormsGraph(2, 3)
         sage: G.is_distance_regular(True)
         ([20, 18, None], [None, 1, 6])
    """
    q = r * r
    Fr = GF(r)
    Fq = GF(q)
    i = Fq.gen()
    ir = i**r

    toR = {(a + i*b): (a + ir*b) for a, b in itertools.product(Fr, repeat=2)}

    def build_mat(v, w):
        # get upper diagonal entries
        res = []
        used_v = 0
        used_w = 0
        for row in range(n):
            res += [v[used_v]] + [v[used_v + 1 + j] + i * w[used_w + j]
                                  for j in range(n - 1 - row)]
            used_v += n - row
            used_w += n - 1 - row

        return tuple(res)

    # produce all rank1 matrices
    rank1Matrices = []
    for w1 in VectorSpace(Fr, n):
        if not w1.is_zero():
            # build matrix
            nonZero = 0
            while nonZero < n and w1[nonZero] == 0:
                nonZero += 1

            for w2 in VectorSpace(Fr, n - nonZero - 1):
                # get upper triangular entries
                sig_check()

                v = [w1[nonZero]] + \
                    [w1[nonZero + 1 + j] + i * w2[j]
                     for j in range(n - nonZero - 1)]

                res = []
                for row in range(nonZero):
                    res += [0] * (n - row)

                res += v

                for row in range(1, n - nonZero):
                    factor = toR[v[row]] / v[0]
                    res += list(map(lambda x: factor * x, v[row:]))

                rank1Matrices.append(res)

    Vs = VectorSpace(Fr, (n * (n+1)) // 2)
    Va = VectorSpace(Fr, (n * (n-1)) // 2)

    edges = []
    for a, b in itertools.product(Vs, Va):
        M = build_mat(a, b)
        for R in rank1Matrices:
            N = tuple([M[i] + R[i] for i in range((n * (n+1)) // 2)])
            edges.append((M, N))

    G = Graph(edges, format='list_of_edges')
    G.name(f"Hermitian forms graph on (F_{q})^{n}")
    return G
