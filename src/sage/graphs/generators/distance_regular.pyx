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
    import itertools

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

        sage: G = graphs.locally_GQ42_distance_transitive_graph()  # optional - internet
        sage: G.is_distance_regular(True)  # optional - internet
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
    import itertools

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

        sage: G = graphs.graph_3O73()  # optional - internet
        sage: G.is_distance_regular(True)  # optional - internet
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

        sage: G = graphs.J2Graph()  # optional - internet
        sage: G.is_distance_regular(True) # optional - internet
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

        sage: G = graphs.IvanovIvanovFaradjevGraph()  # optional - internet
        sage: G.is_distance_regular(True)  # optional - internet
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
