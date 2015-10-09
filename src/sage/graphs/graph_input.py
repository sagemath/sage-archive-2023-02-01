r"""
Functions for reading/building graphs/digraphs.

This module gathers functions needed to build a graph from any other data.

Note that because they are called by the constructors of :class:`Graph` and
:class:`DiGraph`, most of these functions modify a graph inplace.

{INDEX_OF_FUNCTIONS}

Functions
---------

"""

def from_graph6(G, g6_string):
    r"""
    Fill ``G`` with the data of a graph6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a graph6 string

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_graph6
        sage: g = Graph()
        sage: from_graph6(g, 'IheA@GUAo')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True

    """
    from generic_graph_pyx import length_and_string_from_graph6, binary_string_from_graph6

    if not isinstance(g6_string, str):
        raise ValueError('If input format is graph6, then g6_string must be a string.')
    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    ss = g6_string[:n]
    n, s = length_and_string_from_graph6(ss)
    m = binary_string_from_graph6(s, n)
    expected = n*(n-1)/2 + (6 - n*(n-1)/2)%6
    if len(m) > expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
    elif len(m) < expected:
        raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
    G.add_vertices(range(n))
    k = 0
    for i in xrange(n):
        for j in xrange(i):
            if m[k] == '1':
                G._backend.add_edge(i, j, None, False)
            k += 1

def from_sparse6(G, g6_string):
    r"""
    Fill ``G`` with the data of a sparse6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a sparse6 string

    EXAMPLE::

        sage: from sage.graphs.graph_input import from_sparse6
        sage: g = Graph()
        sage: from_sparse6(g, ':I`ES@obGkqegW~')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from generic_graph_pyx import length_and_string_from_graph6, int_to_binary_string
    from math import ceil, floor
    from sage.misc.functional import log
    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    s = g6_string[:n]
    n, s = length_and_string_from_graph6(s[1:])
    if n == 0:
        edges = []
    else:
        k = int(ceil(log(n,2)))
        ords = [ord(i) for i in s]
        if any(o > 126 or o < 63 for o in ords):
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        bits = ''.join([int_to_binary_string(o-63).zfill(6) for o in ords])
        b = []
        x = []
        for i in xrange(int(floor(len(bits)/(k+1)))):
            b.append(int(bits[(k+1)*i:(k+1)*i+1],2))
            x.append(int(bits[(k+1)*i+1:(k+1)*i+k+1],2))
        v = 0
        edges = []
        for i in xrange(len(b)):
            if b[i] == 1:
                v += 1
            if x[i] > v:
                v = x[i]
            else:
                if v < n:
                    edges.append((x[i],v))
    G.add_vertices(range(n))
    G.add_edges(edges)

from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__]))
