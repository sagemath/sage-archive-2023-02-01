from sage.graphs.graph import Graph
from sage.graphs import graph
from sage.graphs.base.dense_graph import DenseGraph
def EgawaGraph(int p, int s):
	r"""
    Returns the Egawa graph with parameters `p, s`.

    Egawa graphs are a peculiar family of graphs devised by Yoshimi Egawa in 
    https://doi.org/10.1016/0097-3165(81)90007-8
    The Shrikhande graph is a special case of this family of graphs, with parameters
    (1,0). All the graphs in this family are not recognizable by 1-WL 
    (Weisfeiler Lehamn algorithm of the first order) and 2-WL,
    that is their orbits are not correctly returned by k-WL for k lower than 3.
    
    Furthermore, all the graphs in this family are distance-regular, but they are
    not distance-transitive if `p` != 0.
    
    Any Egawa graph with parameters (0, s) is isomorphic to the Hamming graph
    with parameters (s, 4)
    
    EXAMPLES:

    Every Egawa graph is distance regular.  ::

        sage: g = graphs.EgawaGraph(2, 3)
        sage: g.is_distance_regular()
        True

    An Egawa graph with parameters (0,s) is isomorphic to the Hamming graph
    with parameters (s, 4).  ::

        sage: g = graphs.EgawaGraph(0, 4)
        sage: g.is_isomorphic(graphs.HammingGraph(4,4))
        True
    """
	from sage.graphs.generators.basic import CompleteGraph
	from itertools import product, chain, repeat
	g = Graph(name="I Graph with parameters " + str(p) + "," + str(s))
	X = CompleteGraph(4)
	Y = Graph('O?Wse@UgqqT_LUebWkbT_')
	g.add_vertices(product(*chain(repeat(Y, p), repeat(X,s))))
	cdef int i,j
	cdef set used = set()
	for v in g:
		for i in range(p):
			for el in Y:
				if el == v[i]: continue
				if Y.has_edge(v[i], el):
					u = v[:i] + (el,) + v[i+1:]
					if (u,v) not in used:
						g.add_edge(v,u)
						used.add((v,u))
		for j in range(s):
			i = p+j
			for el in X:
				if el == v[i]: continue
				u = v[:i] + (el,) + v[i+1:]
				if (u,v) not in used:
					g.add_edge(v,u)
					used.add((v,u))
	return g

def HammingGraph(int n, int q, list X=[]):
	r"""
    Returns the Hamming graph with parameters `n, q` over set `X`.

    Hamming graphs are graphs over the cartesian product of n copies of X,
    where q = |X|, where the vertices, labelled with the corresponding tuple in X^n,
    are connected if the Hamming distance between their labels is 1.
    All Hamming graphs are regular, vertex-transitive and distance-regular.
    
    Hamming graphs with parameters (1,q) represent the complete graph with
    q vertices over the set X.
    
    EXAMPLES:

    Every Hamming graph is distance-regular, regular and vertex-transitive.  ::

        sage: g = graphs.HammingGraph(3, 7)
        sage: g.is_distance_regular()
        True
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True

    An Egawa graph with parameters (1,q) is isomorphic to the Complete graph
    with parameter q.  ::

        sage: g = graphs.HammingGraph(1, 23)
        sage: g.is_isomorphic(graphs.CompleteGraph(23))
        True
    """
	from itertools import product, repeat
	if not X:
		X = list(range(q))
	if q != len(X):
		raise ValueError("q must be the cardinality of X")
	g = Graph(name="Hamming Graph with parameters " + str(n) + "," + str(q))
	g.add_vertices(product(*repeat(X, n)))
	cdef set used = set()
	cdef int i
	for v in g:
		for i in range(n):
			for el in X:
				if el == v[i]: continue
				u = v[:i] + (el,) + v[i+1:]
				if (u,v) not in used:
					g.add_edge(v,u)
					used.add((v,u))
	return g
		
