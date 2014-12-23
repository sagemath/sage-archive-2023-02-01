r"""

Interface for graph automorphisms with bliss.

Implemented functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Returns the automorphism group of the given (di)graph
    :meth:`canonical_form` | Computes a canonical certificate for the given (di) graph. 
    :meth:`is_isomorpic` | Tests whether the passed (di) graphs are isomorphic. 


AUTHORS:

    - Jernej Azarija 

"""
include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

cdef extern from "graph.hh" namespace "bliss":

    cdef cppclass Stats:
        
        pass
    
    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int, 
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int, 
                    const unsigned int*), void*)

    cdef cppclass Digraph(AbstractGraph):

        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int, 
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int, 
                    const unsigned int*), void*)
        unsigned int get_hash()


cdef void add_gen(void *user_param, unsigned int n, const unsigned int *aut):

    covered = cur = 0
    perm = []       
    done = [False]*n

    gens, d2 = <object> <PyObject *> user_param
    
    while covered < n:
        while cur < n and done[cur]:
            cur+=1

        cycle = [d2[aut[cur]]]
        sec = aut[aut[cur]]
        done[aut[cur]] = done[cur] = True

        covered+=1

        while d2[sec] != cycle[0]:
            cycle.append(d2[sec])
            done[sec] = True
            sec = aut[sec]
            covered+=1
        perm+=[tuple(cycle)]
    gens += [perm]

# The function accepts a graph, a coloring  of its vertices (possibly None)
# and two empty dictionaries. The entries of the dicitionary are later set
# to record the labeling of our graph. They are taken as arguments to avoid
# technicalities of returning Python objects in Cython functions.
cdef Graph *bliss_graph(G, partition, vert2int, int2vert):    

    cdef Graph *g = new Graph(G.order())

    if g == NULL:
        raise MemoryError("Allocation Failed")

    for i,v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v 

    for x,y in G.edges(labels=False):
       g.add_edge(vert2int[x],vert2int[y])     

    if partition:
        for i in xrange(1,len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

# This is the same function as bliss_graph with the only exception
# being that it returns a digraph. This code duplication is of course
# ugly and I am open for suggestions (that do not waste time)            
cdef Digraph *bliss_digraph(G, partition, vert2int, int2vert):    

    cdef Digraph *g = new Digraph(G.order()) 

    for i,v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v 
    
    if g == NULL:
        raise MemoryError("Allocation Failed")

    for x,y in G.edges(labels=False):
        g.add_edge(vert2int[x],vert2int[y])     

    if partition:
        for i in xrange(1,len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

def automorphism_group(G, partition=None):
    """
    Computes the automorphism group of ``G`` subject to the coloring ``partition.`` 

    INPUT:

    - ``G`` -- A graph
    - ``partition`` -- A partition of the vertices of ``G`` into color classes. 
        Defaults to ``none``.

    TESTS::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                        # optional - bliss
        sage: automorphism_group(G).is_isomorphic(G.automorphism_group())       # optional - bliss
        True                                                                    # optional - bliss

        sage: G = graphs.HeawoodGraph()                                         # optional - bliss
        sage: p = G.bipartite_sets()                                            # optional - bliss
        sage: A = G.automorphism_group(partition=[list(p[0]), list(p[1])])      # optional - bliss
        sage: automorphism_group(G, partition=p).is_isomorphic(A)               # optional - bliss
        True                                                                    # optional - bliss
    """

    cv = 0 
    n = G.order()
    vert2int = {}
    int2vert = {}

    cdef Graph *g = NULL
    cdef Digraph *d = NULL
    cdef Stats s

    gens = [] 
    data = (gens, int2vert)        

    if G.is_directed(): 
        d = bliss_digraph(G, partition, vert2int, int2vert)
        d.find_automorphisms(s, add_gen, <PyObject *> data)
    else:
        g = bliss_graph(G, partition, vert2int, int2vert) 
        g.find_automorphisms(s, add_gen, <PyObject *> data)
    del g 
    del d

    from sage.groups.perm_gps.permgroup import PermutationGroup

    return PermutationGroup(gens,domain=G)

cdef void empty_hook(void *user_param , unsigned int n, const unsigned int *aut):
    return

def canonical_form(G, partition=None, return_graph=False, certify=False):
    """
    Returns a Python object such that two graphs ``G`` and ``H`` are 
    isomorphic if and only if canonical_form(G) == canonical_form(H).

    INPUT:

    - ``G`` -- A graph or digraph.
    - ``partition`` -- A partition of the vertices of ``G`` into color classes. 
        Defaults to ``none``.
    - ``return_graph`` -- If set to True, canonical_form returns the canonical graph 
        of G. 
    - ``certify`` -- If set to True returns the labeling of G into a canonical graph.

    TESTS::
        sage: from sage.graphs.bliss import canonical_form                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(G)                                             # optional - bliss
        [(2, 0), (2, 1), (3, 0), (4, 1), (4, 6), (5, 3), (5, 4), (6, 0), 
        (7, 1), (7, 3), (8, 2), (8, 5), (9, 6), (9, 7), (9, 8)]             # optional - bliss
        
        sage: P = graphs.GeneralizedPetersenGraph(5,2)                      # optional - bliss
        sage: Q = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(P) == canonical_form(Q)                        # optional - bliss
        True
    """
    cdef Graph *g = NULL
    cdef Digraph *d = NULL
    cdef const unsigned int *aut
    cdef Stats s

    vert2int = {}
    
    if G.is_directed():
        d = bliss_digraph(G, partition, vert2int, {}) 
        aut = d.canonical_form(s, empty_hook, NULL)
    else:
        g = bliss_graph(G, partition, vert2int, {})
        aut = g.canonical_form(s, empty_hook, NULL)     
     
    edges = []
    for x,y in G.edges(labels=False):
        e,f = aut[ vert2int[x] ], aut[ vert2int[y] ]
        edges.append( (e,f) if e > f else (f,e))

    del g
    del d 

    if return_graph:
        if G.is_directed():
            from sage.graphs.graph import DiGraph
            G = DiGraph(edges)
        else:
            from sage.graphs.graph import Graph
            G = Graph(edges)

        return G, vert2int if certify else G

    if certify:
        return sorted(edges),vert2int

    return sorted(edges)

# FIXME if cert=False then the best way to do this would actually
# be to compute the two canonical forms of G and H and use bliss::Graph::cmp
# which I somehow cannot make it work? Can someone look into that?
def is_isomorphic(G,H, cert=False):
    """
    Tests ``G`` and ``H`` for (di) graph isomorphism.

    INPUT:

    - ``G`` -- A graph
    - ``cert`` -- If set to True returns the actual isomorphism between 
        the vertices of ``G`` and ``H`` 

    TESTS::

        sage: from sage.graphs.bliss import is_isomorphic           # optional - bliss
        sage: G1 = graphs.RandomGNP(10,0.5)                         # optional - bliss
        sage: G2 = graphs.RandomGNP(10,0.8)                         # optional - bliss
        sage: is_isomorphic(G1,G2)
        False

        sage: G = graphs.CubeGraph(3)                               # optional - bliss
        sage: H = G.copy()                                          # optional - bliss
        sage: H.relabel()                                           # optional - bliss
        sage: is_isomorphic(G,H,cert=True)                          # optional - bliss
        {'000': 3, '001': 2, '010': 0, '011': 1, '100': 6, '101': 7, '110': 5, '111': 4}

    """
   

    c1,lab1 = canonical_form(G, certify=True)
    c2,lab2 = canonical_form(H, certify=True)
   
    if c1 == c2:
        if cert:
            lab2inv = { lab2[key]:key for key in lab2}
            return { v: lab2inv[lab1[v]] for v in G}
        else:
            return True
    return False  
