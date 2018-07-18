from sage.graphs.graph import DiGraph
from sage.graphs import graph
def CaiFurerImmermanGraph(k):
	from itertools import repeat as rep, chain, combinations
	G = DiGraph()
	G.add_vertices(enumerate(rep('a', k)))
	G.add_vertices(enumerate(rep('b', k)))
	powerset = list(chain.from_iterable(combinations(range(k), r) for r in range(k+1) if r % 2 == 0))
	G.add_vertices(powerset)
	G.add_edges(chain.from_iterable([(s,(i,'a')) for i in s] for s in powerset))
	G.add_edges(chain.from_iterable([(s,(i,'b')) for i in range(k) if i not in s] for s in powerset))
	return G
	
