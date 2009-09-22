r"""
Cliquer : routines for finding cliques in graphs

This module defines functions based on Cliquer, an exact branch-and-bound algorithm developed by Patric R. J. Ostergard and written by Sampo Niskanen.

AUTHORS:

- Nathann Cohen (2009-08-14): Initial version

Functions
------------------------------

"""
def max_clique(graph):
    """
    Returns the vertex set of a maximum complete subgraph.

    Currently only implemented for undirected graphs. Use
    to_undirected to convert a digraph to an undirected graph.

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: max_clique(C)
          [7, 9]
    """
    graph,d = graph.relabel(inplace=False, return_map=True)
    d_inv = {}
    for v in d:
        d_inv[d[v]] = v

    cdef graph_t *g
    g=graph_new(graph.order())
    buf="p edges %d %d" %(graph.order(),graph.size())
    parse_input(buf,g)

    for e in graph.edge_iterator():
          (u,v,w)=e
          buf=' e %d %d' %(u+1,v+1)
          parse_input(buf,g)

    cdef int* list
    cdef int size
    size=sage_clique_max(g, &list)
    b = []
    cdef int i
    for i in range(size):
        b.append(list[i])
    return list_composition(b,d_inv)


# computes all the maximum clique of a graph and return its list

def all_max_clique(graph):
    """
    Returns the vertex sets of *ALL* the maximum complete subgraphs.

    Currently only implemented for undirected graphs. Use
    to_undirected to convert a digraph to an undirected graph.

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: all_max_clique(C)
          [[2, 7], [7, 9], [6, 8], [6, 9], [0, 4], [4, 9], [5, 7], [0, 5], [5, 8], [3, 4], [2, 3], [3, 8], [1, 6], [0, 1], [1, 2]]
    """

    graph,d = graph.relabel(inplace=False, return_map=True)
    d_inv = {}
    for v in d:
        d_inv[d[v]] = v

    cdef graph_t *g
    g=graph_new(graph.order())
    buf="p edges %d %d" %(graph.order(),graph.size())
    parse_input(buf,g)

    for e in graph.edge_iterator():
          (u,v,w)=e
          buf=' e %d %d' %(u+1,v+1)
          parse_input(buf,g)

    cdef int* list
    cdef int size
    size=sage_all_clique_max(g, &list)
    b = []
    c=[]
    cdef int i
    for i in range(size):
        if(list[i]!=-1):
                c.append(list[i])
        else:
                b.append(list_composition(c,d_inv))
                c=[]
    return b


#computes the clique number of a graph

def clique_number(graph):
    """
    Returns the size of the largest clique of the graph (clique
    number).

    Currently only implemented for undirected graphs. Use
    to_undirected to convert a digraph to an undirected graph.

    EXAMPLES::

        sage: C = Graph('DJ{')
        sage: clique_number(C)
        4
        sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
        sage: G.show(figsize=[2,2])
        sage: clique_number(G)
        3
    """
    graph=graph.relabel(inplace=False)
    cdef graph_t *g
    g=graph_new(graph.order())
    buf="p edges %d %d" %(graph.order(),graph.size())
    parse_input(buf,g)

    for e in graph.edge_iterator():
          (u,v,w)=e
          buf=' e %d %d' %(u+1,v+1)
          parse_input(buf,g)
    return sage_clique_number(g)


def list_composition(a,b):
    """
    Composes a list ``a`` with a map ``b``.

    EXAMPLE::

        sage: from sage.graphs.cliquer import list_composition
        sage: list_composition([1,3,'a'], {'a':'b', 1:2, 2:3, 3:4})
        [2, 4, 'b']

    """
    value=[]
    for i in a:
        value.append(b[i])
    return value


