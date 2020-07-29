r"""
<<<<<<< HEAD
Database of distance regular graphs
=======
Dabase of distance regular graphs
>>>>>>> moved file to generators; fixed tests

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

<<<<<<< HEAD
    This is a distance-regular graph with intersection array
=======
    This is a distance-regular graph with intersecion array
>>>>>>> moved file to generators; fixed tests
    `[15, 14, 10, 3; 1, 5, 12, 15]`.

    EXAMPLES::

        sage: G = graphs.cocliques_HoffmannSingleton()
        sage: G.is_distance_regular(True)
        ([15, 14, 10, 3, None], [None, 1, 5, 12, 15])

    REFERENCES:

    The construction of this graph can be found in [BCN1989]_ p. 392.
    """
    from sage.graphs.graph_generators import GraphGenerators
<<<<<<< HEAD
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
=======
    D = GraphGenerators.HoffmanSingletonGraph()
    DC = D.complement()

    cocliques = DC.cliques_maximum()  # 100 of this

    edges = []
    for i in range(100):
        sC = frozenset(cocliques[i])
        for j in range(i+1,100):
            if len(sC.intersection(cocliques[j])) == 8:
                sC2 = frozenset(cocliques[j])
                edges.append( (sC,sC2) )

    G = Graph(edges,format="list_of_edges")
    return G


def locally_GQ42_graph():
    r"""
    Return the unique amply regular graph which is locally a generalised
    quadrangle.
>>>>>>> moved file to generators; fixed tests

    This graph is distance-regular with intersection array
    `[45, 32, 12, 1; 1, 6, 32, 45]`.

    This graph is also distance-transitive.

    EXAMPLES::

<<<<<<< HEAD
        sage: G = graphs.locally_GQ42_distance_transitive_graph()  # optional - internet
        sage: G.is_distance_regular(True)  # optional - internet
        ([45, 32, 12, 1, None], [None, 1, 6, 32, 45])

    REFERENCES:

    A description of this graph can be found in [BCN1989]_ p.399.
    This construction is due to Dima Pasechnik.
    """
    H = libgap.AtlasGroup("3^2.U4(3).D8", libgap.NrMovedPoints, 756)
=======
        sage: G = graphs.locally_GQ42_graph()  # optional - gap_packages
        sage: G.is_distance_regular(True)  # optional - gap_packages
        ([45, 32, 12, 1, None], [None, 1, 6, 32, 45])

    .. NOTE::

        This function needs the GAP's package AtlasRep [WPNBBAtl]_.
        Install it via ``sage -i gap_packages``.
    """
    H = libgap.AtlasGroup("3^2.U4(3).D8",libgap.NrMovedPoints,756)
>>>>>>> moved file to generators; fixed tests
    Ns = H.NormalSubgroups()
    for N in Ns:
        if len(N.GeneratorsSmallest()) == 7:  # there is only one
            break

<<<<<<< HEAD
    G = Graph(libgap.Orbit(N, [1, 9], libgap.OnSets), format='list_of_edges')
    G.name("locally GQ(4,2) distance transitive graph")
=======
    G = Graph(libgap.Orbit(N,[1,9],libgap.OnSets), format='list_of_edges')
    G.name("locally GQ(4,2) graph")
>>>>>>> moved file to generators; fixed tests
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
<<<<<<< HEAD

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 399.
    """
    from sage.rings.number_field.number_field import CyclotomicField
    import itertools
=======
    """
    from sage.rings.number_field.number_field import CyclotomicField
>>>>>>> moved file to generators; fixed tests

    F = CyclotomicField(3)
    w = F.gen()

    V= VectorSpace(GF(4), 6)
<<<<<<< HEAD
    z2 = GF(4)('z2') # GF(4) = {0, 1, z2, z2+1}
=======
    z2 = GF(4)('z2') # GF(4) = {0,1,z2, z2+1}
>>>>>>> moved file to generators; fixed tests

    W = V.span([(0,0,1,1,1,1), (0,1,0,1,z2,z2+1), (1,0,0,1,z2+1,z2)])
    # we only need the 45 vectors with 2 zero entries
    # we also embed everything into CC

    K = []
    for v in W:
<<<<<<< HEAD
        # check zero entries
        zeros = 0
        for x in v:
            if x.is_zero():
=======
        #check zero entries
        zeros = 0
        for x in v:
            if x == 0:
>>>>>>> moved file to generators; fixed tests
                zeros += 1

        if zeros == 2:
            # send to F and in K
            # z2 -> w
            # z2+1 -> w^2
            vv = []  # new vector
            for x in v:
                if x == z2:
                    vv.append(w)
<<<<<<< HEAD
                elif x == z2 + 1:
                    vv.append(w**2)
                else:
                    vv.append(int(x))
=======
                elif x == z2+1:
                    vv.append(w**2)
                else:
                    vv.append(int(x))  # this is weirdly needed for some reason
>>>>>>> moved file to generators; fixed tests

            # now vv is the new vector in F
            vv = vector(F, vv)
            K.append(vv)

    # we need to add other vectors
    for i in range(6):
<<<<<<< HEAD
        # create e_i
        ei = [0, 0, 0, 0, 0, 0]
=======
        #create e_i
        ei = [0]*6
>>>>>>> moved file to generators; fixed tests
        ei[i] = 1
        ei = vector(F, ei)

        K.append(2 * ei)
        K.append(2 * w * ei)
        K.append(2 * w**2 * ei)
    # now K is all the 63 vertices

<<<<<<< HEAD
    for v in K:
        v.set_immutable()

    def has_edge(u, v):
        return sum(u[i].conjugate() * v[i] for i in range(6)) == 2

    G = Graph()
    for Ki, Kj in itertools.combinations(K, 2):
        if has_edge(Ki, Kj):
            G.add_edge((Ki, Kj))
=======
    def has_edge(u,v):
        com = 0
        for i in range(6):
            com += u[i].conjugate() * v[i]

        if com == 2:
            return True
        return False

    G = Graph()
    length = len(K)
    for i in range(length):
        K[i].set_immutable()
        for j in range(i+1, length):
            if has_edge(K[i], K[j]):
                K[j].set_immutable()
                G.add_edge((K[i], K[j]))
>>>>>>> moved file to generators; fixed tests

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

<<<<<<< HEAD
        sage: G = graphs.graph_3O73()  # optional - internet
        sage: G.is_distance_regular(True)  # optional - internet
        ([117, 80, 24, 1, None], [None, 1, 12, 80, 117])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 400.
    """
    group = libgap.AtlasGroup("3.O7(3)", libgap.NrMovedPoints, 1134)
=======
        sage: G = graphs.graph_3O73()  # optional - gap_packages
        sage: G.is_distance_regular(True)  # optional - gap_packages
        ([117, 80, 24, 1, None], [None, 1, 12, 80, 117])

    .. NOTE::

        This function needs the GAP's package AtlasRep [WPNBBAtl]_.
        Install it via ``sage -i gap_packages``.
    """
    group = libgap.AtlasGroup("3.O7(3)",libgap.NrMovedPoints,1134)
>>>>>>> moved file to generators; fixed tests
    G = Graph(libgap.Orbit(group, [1, 3], libgap.OnSets), format='list_of_edges')
    G.name("Distance transitive graph with automorphism group 3.O_7(3)")
    return G
