r"""
Tutte polynomial

This module implements a deletion-contraction algorithm for computing as described in the
paper 

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`tutte_polynomial` | Computes the Tutte polynomial of the input graph

Author:

- mhansen (06-2013), Implemented the algorithm.
- Jernej Azarija (06-2013), Tweaked the code, merged it into Sage and documented it.

Definition
-----------

Given a graph `G`, with `n` vertices and `m` edges and `k(G)` connected components we define 
the Tutte polynomial of `G` as

.. MATH::

    \sum (x-1) ^(k(H) - c) (y-1)^{k(H) - |E(H)|-n}
    
where the sum ranges over all induced subgraphs `H` of `G`,


..[Gordon10]  Computing Tutte Polynomials. Gary Haggard, David J. Pearce and
    Gordon Royle. In ACM Transactions on Mathematical Software, Volume
    37(3), article 24, 2010. Preprint:
    http://homepages.ecs.vuw.ac.nz/~djp/files/TOMS10.pdf

Functions
---------
"""


from contextlib import contextmanager
from sage.misc.lazy_attribute import lazy_attribute
from sage.graphs.graph import Graph
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
R = ZZ['x,y']
x,y = R.gen(0),R.gen(1)

######################
# Graph Modification #
######################

@contextmanager
def removed_multiedge(G, edge, multiplicity):
    """
    A context manager which removes an edge with multiplicity from the
    graph and restores it upon exiting.
    """
    for i in range(multiplicity):
        G.delete_edge(edge)
    try:
        yield
    finally:
        for i in range(multiplicity):
            G.add_edge(edge)

@contextmanager
def removed_edge(G, edge):
    """
    A context manager which removes an edge from the graph and
    restores it upon exiting.
    """
    G.delete_edge(edge)
    try:
        yield
    finally:
        G.add_edge(edge)

@contextmanager
def contracted_edge(G, edge):
    """
    Delete the first vertex in the edge, and make all the edges that
    went from it go to the second vertex.
    """
    v1, v2 = edge[:2]

    v1_edges = G.edges_incident(v1)
    G.delete_vertex(v1)
    added_edges = []

    for start, end, label in v1_edges:
        other_vertex = start if start != v1 else end
        edge = (other_vertex, v2, label)
        G.add_edge(edge)
        added_edges.append(edge)

    try:
        yield
    finally:
        for edge in added_edges:
            G.delete_edge(edge)
        for edge in v1_edges:
            G.add_edge(edge)

@contextmanager
def removed_loops(G):
    """
    A context manager which removes all the loops in the graph `G`.
    It yields a list of the the loops, and restores the loops upon
    exiting.
    """
    loops = G.loops()
    for edge in loops:
        G.delete_edge(edge)
    try:
        yield loops
    finally:
        for edge in loops:
            G.add_edge(edge)


def underlying_graph(G):
    """
    Given a graph `G` with multi-edges, returns a graph where all the
    multi-edges are replaced with a single edge.
    """
    g = Graph()
    g.allow_loops(True)
    for edge in set(G.edges(labels=False)):
        g.add_edge(edge)
    return g

def edge_multiplicites(G):
    """
    Returns the a dictionary of multiplicites of the edges in the
    graph.
    """
    d = {}
    for edge in G.edges(labels=False):
        d[edge] = d.setdefault(edge, 0) + 1
    return d

########
# Ears #
########

class Ear(object):
    """
    An ear is a sequence of vertices
    """
    def __init__(self, graph, end_points, interior, is_cycle):
        self.end_points = end_points
        self.interior = interior
        self.is_cycle = is_cycle
        self.graph = graph

    @property
    def s(self):
        """
        Returns the number of distinct edges in this ear
        """
        return len(self.interior) + 1

    @property
    def vertices(self):
        """
        Returns the vertices of this ear.
        """
        return sorted(self.end_points + self.interior)

    @lazy_attribute
    def edges(self):
        """
        Returns the edges in this ear.
        """
        return self.graph.edges_incident(vertices=self.interior, labels=False)

    @staticmethod
    def find_ear(g):
        """
        Finds the first ear in a graph.
        """
        degree_two_vertices = [v for v, degree in g.degree_iterator(labels=True) if degree == 2]
        subgraph = g.subgraph(degree_two_vertices)
        for component in subgraph.connected_components():
            edges = g.edges_incident(vertices=component, labels=False)
            all_vertices = list(sorted(set(sum(edges, tuple()))))
            if len(all_vertices) < 3:
                continue
            end_points = [v for v in all_vertices if v not in component]
            if not end_points:
                end_points = [component[0]]

            ear_is_cycle = end_points[0] == end_points[-1]

            if ear_is_cycle:
                for e in end_points:
                    if e in component:
                        component.remove(e)

            return Ear(g, end_points, component, ear_is_cycle)

    @contextmanager
    def removed_from(self, G):
        """
        A context manager which removes the ear from the graph `G`.
        """
        deleted_edges = []
        for edge in G.edges_incident(vertices=self.interior, labels=False):
            G.delete_edge(edge)
            deleted_edges.append(edge)
        for v in self.interior:
            G.delete_vertex(v)

        try:
            yield
        finally:
            for edge in deleted_edges:
                G.add_edge(edge)

##################
# Edge Selection #
##################
class EdgeSelection(object):
    pass

class VertexOrder(EdgeSelection):
    def __init__(self, order):
        self.order = list(order)
        self.inverse_order = dict(map(reversed, enumerate(order)))

    def __call__(self, graph):
        for v in self.order:
            edges = graph.edges_incident([v], labels=False)
            if edges:
                edges.sort(key=lambda x: self.inverse_order[x[0] if x[0] != v else x[1]])
                return edges[0]
        raise RuntimeError, "no edges left to select"

class MinimizeSingleDegree(EdgeSelection):
    def __call__(self, graph):
        degrees = list(graph.degree_iterator(labels=True))
        degrees.sort(key=lambda x: x[1]) #Sort by degree
        for v, degree in degrees:
            for e in graph.edges_incident([v], labels=False):
                return e
        raise RuntimeError, "no edges left to select"

class MinimizeDegree(EdgeSelection):
    def __call__(self, graph):
        degrees = dict(graph.degree_iterator(labels=True))
        edges = graph.edges(labels=False)
        edges.sort(key=lambda x: degrees[x[0]]+degrees[x[1]]) #Sort by degree
        for e in edges:
            return e
        raise RuntimeError, "no edges left to select"

class MaximizeDegree(EdgeSelection):
    def __call__(self, graph):
        degrees = dict(graph.degree_iterator(labels=True))
        edges = graph.edges(labels=False)
        edges.sort(key=lambda x: degrees[x[0]]+degrees[x[1]]) #Sort by degree
        for e in reversed(edges):
            return e
        raise RuntimeError, "no edges left to select"

###########
# Caching #
###########

def _cache_key(G):
    return tuple(sorted(G.canonical_label().edges(labels=False)))

def _cached(func):
    _cache = {}
    def wrapper(G, **kwds):
        key = _cache_key(G)
        if key in _cache:
            return _cache[key]
        result = func(G, **kwds)
        _cache[key] = result
        return result
    wrapper.cache = _cache
    wrapper.original_func = func
    return wrapper

####################
# Tutte Polynomial #
####################

@_cached
def tutte_polynomial(G, initial_call=True, edge_selector=None):
    """
    Returns the Tutte polynomial of the graph `G`. 

    INPUT:
    
    - ``edge_selector`` (method) This (optional) argument allows the user
    to specify his own heuristic for selecting edges used in the deletion 
    contraction recurrence.

    EXAMPLES:

    The Tutte polynomial of any tree of order `n` is `x^{n-1}`::

        sage: all(T.tutte_polynomial() == x**9 for T in graphs.trees(10))
        True

    The Tutte polynomial of the Petersen graph is::

        sage: P = graphs.PetersenGraph()
        sage: P.tutte_polynomial()
        x^9 + 6*x^8 + 21*x^7 + 56*x^6 + 12*x^5*y + y^6 + 114*x^5 + 70*x^4*y + 30*x^3*y^2 + 15*x^2*y^3 + 10*x*y^4 + 9*y^5 + 170*x^4 + 170*x^3*y + 105*x^2*y^2 + 65*x*y^3 + 35*y^4 + 180*x^3 + 240*x^2*y + 171*x*y^2 + 75*y^3 + 120*x^2 + 168*x*y + 84*y^2 + 36*x + 36*y


    The Tutte polynomial of `G``evaluated at (1,1) is the number of spanning trees of `G`::

        sage: G = graphs.RandomGNP(10,0.6)
        sage: G.tutte_polynomial()(1,1) == G.spanning_trees_count()
        True

    Given that `T(x,y)` is the Tutte polynomial of a graph `G` with `n` vertices and `c` connected components, then 
    (-1)^{n-c} x^k T(1-x,0) is the chromatic polynomial of `G`.

        sage: G = graphs.OctahedralGraph() 
        sage: T = G.tutte_polynomial()
        sage: (-1)^5*x*T(1-x,0).factor() # comparison (==) of this with G.chrompoly() wont work???!
        (x - 2)*(x - 1)*(x^3 - 9*x^2 + 29*x - 32)*x
        sage: G.chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^3 - 9*x^2 + 29*x - 32)


    """
    if G.num_edges() == 0:
        return x**0 * y**0

    if initial_call is True:
        G = G.copy()
        G.allow_loops(True)
        G.allow_multiple_edges(True)

    if edge_selector is None:
        edge_selector = MinimizeSingleDegree()

    def recursive_tp(graph=None):
        """
        The recursive call -- used so that we don't have to specify
        the same arguments everywhere.
        """
        if graph is None:
            graph = G
        return tutte_polynomial(graph, initial_call=False, edge_selector=edge_selector)

    #Remove loops
    with removed_loops(G) as loops:
        if loops:
            return y**len(loops) * recursive_tp()

    uG = underlying_graph(G)
    em = edge_multiplicites(G)
    d = em.values()

    def yy(start, end):
        return sum(y**i for i in range(start, end+1))

    #Lemma 1
    if G.is_forest():
        return prod(x + yy(1, d_i-1) for d_i in d)

    #Handle disconnected components
    if not G.is_connected():
        return prod([recursive_tp(G.subgraph(block)) for block in G.connected_components()])

    #Theorem 1: from Haggard, Pearce, Royle 2008
    blocks, cut_vertices = G.blocks_and_cut_vertices()
    if len(blocks) > 1:
        return prod([recursive_tp(G.subgraph(block)) for block in blocks])

    components = G.connected_components_number()
    edge = edge_selector(G)

    with removed_edge(G, edge):
        if G.connected_components_number() > components:
            with contracted_edge(G, edge):
                return x*recursive_tp()

    ##################################
    # We are in the biconnected case #
    ##################################

    # Theorem 4: from Haggard, Pearce, and Royle 
    # Note that the formula at
    # http://homepages.ecs.vuw.ac.nz/~djp/files/TOMS10.pdf is slightly
    # incorrect.  The initial sum should only go to n-2 instead of n
    # (allowing for the last part of the recursion).  Additionally, the first  
    # operand of the final product should be (x+y^{1...(d_n+d_{n-1}-1)}) instead of just
    # (x+y^(d_n+d_{n-1}-1)
    if uG.num_verts() == uG.num_edges():  #G is a multi-cycle
        n = len(d)
        result = 0
        for i in range(n-2):
            term =(prod((x+yy(1, d_j-1)) for d_j in d[i+1:])*
                   prod((yy(0, d_k-1)) for d_k in d[:i]))
            result += term
        #The last part of the recursion
        result += (x+yy(1, d[-1] + d[-2] - 1))*prod(yy(0, d_i-1) for d_i in d[:-2])
        return result

    # Theorem 3 from Haggard, Pearce, and Royle, adapted to multi-eaars
    ear = Ear.find_ear(uG)
    if ear is not None:
        if (ear.is_cycle and ear.vertices == G.vertices()): #The graph is an ear (cycle)
            #We should never be in this case since we check for multi-cycles above
            return  y + sum(x**i for i in range(1, ear.s))
        else:
            with ear.removed_from(G):
                #result = sum(x^i for i in range(ear.s)) #single ear case
                result = sum((prod(x+yy(1, em[e]-1) for e in ear.edges[i+1:])*
                              prod(yy(0,em[e]-1) for e in ear.edges[:i])) 
                             for i in range(len(ear.edges)))
                result *= recursive_tp()

                with contracted_edge(G, [ear.end_points[0], ear.end_points[-1]]):
                    result += prod(yy(0, em[e]-1) for e in ear.edges)*recursive_tp()
            
            return result

    #Theorem 2
    if len(em) == 1: #the graph is just a multiedge
        return x + sum(y**i for i in range(1, em[edge]))
    else:
        with removed_multiedge(G, edge, em[edge]):
            result = recursive_tp()
            with contracted_edge(G, edge):
                result += sum(y**i for i in range(em[edge]))*recursive_tp()
        return result

