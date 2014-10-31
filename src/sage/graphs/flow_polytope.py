r"""
Procedure to compute the flow polytope of a directed graph.

Flow polytopes of a directed graph is a polytope formed by
assigning a nonnegative flow to each of edges of the graph such that
the flow is conserved on internal vertices, and there is a unit of
flow entering the sources and leaving the sinks.

The faces and volume of these polytopes are of interest. Examples of
these polytopes are the Chan-Robbins-Yuen polytope and the
Pitman-Stanley polytope.
"""
from sage.rings.integer_ring import ZZ


def flow_polytope(G):
    """
    Return the flow polytope of `G`

    EXAMPLES:

    A commutative square::

        sage: from sage.graphs.flow_polytope import flow_polytope
        sage: G = DiGraph({1:[2,3],2:[4],3:[4]})
        sage: fl = flow_polytope(G); fl
        A 1-dimensional polyhedron in QQ^4 defined as the convex hull
        of 2 vertices
        sage: fl.vertices()
        (A vertex at (0, 1, 0, 1), A vertex at (1, 0, 1, 0))

    A tournament on 4 vertices::

        sage: H = digraphs.TransitiveTournament(4)
        sage: fl = flow_polytope(H); fl
        A 3-dimensional polyhedron in QQ^6 defined as the convex hull
        of 4 vertices
        sage: fl.vertices()
        (A vertex at (0, 0, 1, 0, 0, 0),
         A vertex at (0, 1, 0, 0, 0, 1),
         A vertex at (1, 0, 0, 0, 1, 0),
         A vertex at (1, 0, 0, 1, 0, 1))
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    ineqs = [[0] + [ZZ(j==u) for j in G.edges()]
             for u in G.edges()]

    eqs = []
    for u in G:
        ins = G.incoming_edges(u)
        outs = G.outgoing_edges(u)
        eq = [ZZ(j in ins) - ZZ(j in outs) for j in G.edges()]

        if len(ins) == 0:  # sources
            eq = [1] + eq
        elif len(outs) == 0:  # sinks
            eq = [-1] + eq
        else:
            eq = [0] + eq
        eqs.append(eq)

    return Polyhedron(ieqs=ineqs, eqns=eqs)
