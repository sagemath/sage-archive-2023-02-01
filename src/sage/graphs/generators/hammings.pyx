from sage.graphs.graph import Graph
from sage.graphs import graph
from sage.graphs.base.dense_graph import DenseGraph
def IGraph(int p, int s):
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
	from itertools import product, repeat
	if not X:
		X = list(range(q))
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
		
