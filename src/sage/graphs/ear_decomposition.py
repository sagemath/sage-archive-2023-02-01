# -*- coding: utf-8 -*-
"""
This module implements the function for computing the Ear Decomposition
of undirected connected graphs.

Input:     

- ``graph`` -- The undirected connected graph for which ear decomposition needs to be computed

Output:

- A nested list representing the cycle and chains of input graph G.

Example:
	
	sage: from sage.graphs.ear_decomposition import ear_decomposition
	sage: g = Graph([[0, 1],[0, 7],[0, 3],[0, 11],[2, 3],[1, 2],[1, 12],[2, 12],[3, 4],[3, 6],[4, 5],[4, 7],[4, 6],[5, 6],[7, 10],[7, 8],[7, 11],[8, 9],[8, 10],[8, 11],[9, 11]])
	sage: ear_decomposition(g)
	[[0, 11, 7, 4, 3, 2, 1, 0],
	 [0, 3],
	 [0, 7],
	 [1, 12, 2],
	 [3, 6, 4],
	 [4, 5, 6],
	 [7, 10, 8, 11],
	 [7, 8],
	 [11, 9, 8]]

	
	If input is one edge
	sage: g = Graph([[1,2]])
	sage: ear_decomposition(g)
	ValueError: Please input a undirected connected graph with number of vertices > 2

	

	If input graph is not connected
	sage: g = Graph([[1,2], [3,4]])
	sage: ear_decomposition(g)
	ValueError: Graph must be undirected connected

	

	sage: G.allow_loops(True)
	sage: G.allow_multiple_edges(True)
	sage: for i in range(20):
	....:	u = randint(1, 7)
	....:	v = randint(1, 7)
	....:	G.add_edge(u, v)
	sage: G
	Looped multi-graph on 7 vertices
	sage: H = copy(G)
	sage: H
	Looped multi-graph on 7 vertices
	sage: H = copy(G)
	sage: H.allow_loops(False)
	sage: H.allow_multiple_edges(False)
	sage: H
	Graph on 7 vertices
	sage: ear_decomposition(H) == ear_decomposition(G)
	[[1, 3, 2, 1], [1, 6, 5, 4, 3], [1, 7, 6], [2, 5], [2, 7], [3, 7], [4, 7], [5, 7]]
	[[1, 3, 2, 1], [1, 6, 5, 4, 3], [1, 7, 6], [2, 5], [2, 7], [3, 7], [4, 7], [5, 7]]
	True
	
	

	sage: g = Graph([['0', '1'],['0', '7'],['0', '3'],['0', '11'],['2', '3'],['1', '2'],['1', '12'],['2', '12'],['3', '4'],['3', '6'],['4', '5'],['
	....: 4', '7'],['4', '6'],['5', '6'],['7', '10'],['7', '8'],['7', '11'],['8', '9'],['8', '10'],['8', '11'],['9', '11']])
	sage: ear_decomposition(g)
	[['11', '7', '4', '3', '2', '1', '0', '11'], ['11', '8', '10', '7'], ['11', '9', '8'], ['0', '7'], ['0', '3'], ['1', '12', '2'], ['3', '6', '4'], ['4', '5', '6'], ['7', '8']]

"""
from sage.graphs.graph import Graph
class Search:
	"""
	graph                   : Input graph
	time_of_visit           : List to store the order in which dfs visits vertices.
	visited_vertice         : Boolean dict to mark vertices as visited or unvisited during Dfs traversal in graph.
	traverse_visited_vertice: Boolean dict to mark vertices as visited or unvisited in Dfs tree traversal.
	parent                  : Dict to store parent vertex of all the visited vertices.
	value                   : List to store visit_time of vertices in Dfs traversal.
	chains                  : List to store all the chains and cycles of the input graph G.
	DFS()					: Function that performs depth first search on input graph G.
	"""
	
	def __init__(self, graph, vertices):
		self.graph = graph
		self.time_of_visit = []
		self.visited_vertice = {}
		self.traverse_visited_vertice = {}
		self.parent = {}
		self.value = {}
		self.count = 0
		self.chains = []
		
		#iniializing
		for i in vertices:
			self.visited_vertice[i] = False
			self.traverse_visited_vertice[i] = False

	def DFS(self, v):
		"""
			v = the current vertex need to expand.
		"""
		#make v are visited, update it's time of visited and value
		self.visited_vertice[v] = True
		self.time_of_visit.append(v)
		self.value[v] = self.count
		self.count += 1
		
		#Traverse though all the neighbor vertices of v
		for neighbor in self.graph.neighbors(v):
			#if any neighbor is not visited, enter
			if(not self.visited_vertice[neighbor]):
				#Upated neighbor's parent as v and expand neighbor
				self.parent[neighbor] = v
				self.DFS(neighbor)

	def traverse(self,start, end):
		#Make the firt end of non-tree edge visited
		self.traverse_visited_vertice[start] = True
		pointer = end
		chain = []
		chain.append(start)

		#Traverse DFS Tree of G and print all the not visited vertices.
		#Appending all the vertices in chain
		while True:
			if(self.traverse_visited_vertice[pointer]):
				chain.append(pointer)
				break
			if(pointer==start):
				chain.append(pointer)
				break
			chain.append(pointer)
			self.traverse_visited_vertice[pointer] = True
			pointer = self.parent[pointer]
		self.chains.append(chain)


def ear_decomposition(graph):
	n = graph.num_verts()
	m = graph.num_edges()
    
	if graph.is_directed():
		raise ValueError("Graph must be undirected")

	if not graph.is_connected():
		raise ValueError("Graph must be undirected connected")
	
	if(n<3):
		raise ValueError("Please input a undirected connected graph with number of vertices > 2")

	vertices = graph.get_vertices().keys()
	search = Search(graph, vertices)
	search.parent[vertices[0]] = -1
	
	#start the depth first search from first vertex
	search.DFS(vertices[0])

	#Traverse all the non Tree edges, according to depth first traversal
	for i in range(n):
		for neighbor in graph.neighbors(search.time_of_visit[i]):
			if(search.value[search.time_of_visit[i]] < search.value[neighbor] and search.time_of_visit[i] != search.parent[neighbor]):
				search.traverse(search.time_of_visit[i],neighbor)

	return search.chains
