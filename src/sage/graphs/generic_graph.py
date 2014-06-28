r"""
Generic graphs

This module implements the base class for graphs and digraphs, and methods that
can be applied on both. Here is what it can do:

**Basic Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.networkx_graph` | Creates a new NetworkX graph from the Sage graph
    :meth:`~GenericGraph.to_dictionary` | Creates a dictionary encoding the graph.
    :meth:`~GenericGraph.copy` | Return a copy of the graph.
    :meth:`~GenericGraph.adjacency_matrix` | Returns the adjacency matrix of the (di)graph.
    :meth:`~GenericGraph.incidence_matrix` | Returns an incidence matrix of the (di)graph
    :meth:`~GenericGraph.distance_matrix` | Returns the distance matrix of the (strongly) connected (di)graph
    :meth:`~GenericGraph.weighted_adjacency_matrix` | Returns the weighted adjacency matrix of the graph
    :meth:`~GenericGraph.kirchhoff_matrix` | Returns the Kirchhoff matrix (a.k.a. the Laplacian) of the graph.
    :meth:`~GenericGraph.get_boundary` | Returns the boundary of the (di)graph.
    :meth:`~GenericGraph.set_boundary` | Sets the boundary of the (di)graph.
    :meth:`~GenericGraph.has_loops` | Returns whether there are loops in the (di)graph.
    :meth:`~GenericGraph.allows_loops` | Returns whether loops are permitted in the (di)graph.
    :meth:`~GenericGraph.allow_loops` | Changes whether loops are permitted in the (di)graph.
    :meth:`~GenericGraph.loops` | Returns any loops in the (di)graph.
    :meth:`~GenericGraph.has_multiple_edges` | Returns whether there are multiple edges in the (di)graph.
    :meth:`~GenericGraph.allows_multiple_edges` | Returns whether multiple edges are permitted in the (di)graph.
    :meth:`~GenericGraph.allow_multiple_edges` | Changes whether multiple edges are permitted in the (di)graph.
    :meth:`~GenericGraph.multiple_edges` | Returns any multiple edges in the (di)graph.
    :meth:`~GenericGraph.name` | Returns or sets the graph's name.
    :meth:`~GenericGraph.weighted` | Whether the (di)graph is to be considered as a weighted (di)graph.
    :meth:`~GenericGraph.antisymmetric` | Tests whether the graph is antisymmetric
    :meth:`~GenericGraph.density` | Returns the density
    :meth:`~GenericGraph.order` | Returns the number of vertices.
    :meth:`~GenericGraph.size` | Returns the number of edges.
    :meth:`~GenericGraph.add_vertex` | Creates an isolated vertex.
    :meth:`~GenericGraph.add_vertices` | Add vertices to the (di)graph from an iterable container
    :meth:`~GenericGraph.delete_vertex` | Deletes a vertex, removing all incident edges.
    :meth:`~GenericGraph.delete_vertices` | Remove vertices from the (di)graph taken from an iterable container of vertices.
    :meth:`~GenericGraph.has_vertex` | Return True if vertex is one of the vertices of this graph.
    :meth:`~GenericGraph.random_vertex` | Returns a random vertex of self.
    :meth:`~GenericGraph.random_edge` | Returns a random edge of self.
    :meth:`~GenericGraph.vertex_boundary` | Returns a list of all vertices in the external boundary of vertices1, intersected with vertices2.
    :meth:`~GenericGraph.set_vertices` | Associate arbitrary objects with each vertex
    :meth:`~GenericGraph.set_vertex` | Associate an arbitrary object with a vertex.
    :meth:`~GenericGraph.get_vertex` | Retrieve the object associated with a given vertex.
    :meth:`~GenericGraph.get_vertices` | Return a dictionary of the objects associated to each vertex.
    :meth:`~GenericGraph.loop_vertices` | Returns a list of vertices with loops.
    :meth:`~GenericGraph.vertex_iterator` | Returns an iterator over the vertices.
    :meth:`~GenericGraph.neighbor_iterator` | Return an iterator over neighbors of vertex.
    :meth:`~GenericGraph.vertices` | Return a list of the vertices.
    :meth:`~GenericGraph.neighbors` | Return a list of neighbors (in and out if directed) of vertex.
    :meth:`~GenericGraph.merge_vertices` | Merge vertices.
    :meth:`~GenericGraph.add_edge` | Adds an edge from u and v.
    :meth:`~GenericGraph.add_edges` | Add edges from an iterable container.
    :meth:`~GenericGraph.subdivide_edge` | Subdivides an edge `k` times.
    :meth:`~GenericGraph.subdivide_edges` | Subdivides k times edges from an iterable container.
    :meth:`~GenericGraph.delete_edge` | Delete the edge from u to v
    :meth:`~GenericGraph.delete_edges` | Delete edges from an iterable container.
    :meth:`~GenericGraph.delete_multiedge` | Deletes all edges from u and v.
    :meth:`~GenericGraph.set_edge_label` | Set the edge label of a given edge.
    :meth:`~GenericGraph.has_edge` | Returns True if (u, v) is an edge, False otherwise.
    :meth:`~GenericGraph.edges` | Return a list of edges.
    :meth:`~GenericGraph.edge_boundary` | Returns a list of edges `(u,v,l)` with `u` in ``vertices1``
    :meth:`~GenericGraph.edge_iterator` | Returns an iterator over edges.
    :meth:`~GenericGraph.edges_incident` | Returns incident edges to some vertices.
    :meth:`~GenericGraph.edge_label` | Returns the label of an edge.
    :meth:`~GenericGraph.edge_labels` | Returns a list of edge labels.
    :meth:`~GenericGraph.remove_multiple_edges` | Removes all multiple edges, retaining one edge for each.
    :meth:`~GenericGraph.remove_loops` | Removes loops on vertices in vertices. If vertices is None, removes all loops.
    :meth:`~GenericGraph.loop_edges` | Returns a list of all loops in the graph.
    :meth:`~GenericGraph.number_of_loops` | Returns the number of edges that are loops.
    :meth:`~GenericGraph.clear` | Empties the graph of vertices and edges and removes name, boundary, associated objects, and position information.
    :meth:`~GenericGraph.degree` | Gives the degree (in + out for digraphs) of a vertex or of vertices.
    :meth:`~GenericGraph.average_degree` | Returns the average degree of the graph.
    :meth:`~GenericGraph.degree_histogram` | Returns a list, whose ith entry is the frequency of degree i.
    :meth:`~GenericGraph.degree_iterator` | Returns an iterator over the degrees of the (di)graph.
    :meth:`~GenericGraph.degree_sequence` | Return the degree sequence of this (di)graph.
    :meth:`~GenericGraph.random_subgraph` | Return a random subgraph that contains each vertex with prob. p.
    :meth:`~GenericGraph.add_cycle` | Adds a cycle to the graph with the given vertices.
    :meth:`~GenericGraph.add_path` | Adds a cycle to the graph with the given vertices.
    :meth:`~GenericGraph.complement` | Returns the complement of the (di)graph.
    :meth:`~GenericGraph.line_graph` | Returns the line graph of the (di)graph.
    :meth:`~GenericGraph.to_simple` | Returns a simple version of itself (i.e., undirected and loops and multiple edges are removed).
    :meth:`~GenericGraph.disjoint_union` | Returns the disjoint union of self and other.
    :meth:`~GenericGraph.union` | Returns the union of self and other.
    :meth:`~GenericGraph.relabel` | Relabels the vertices of ``self``
    :meth:`~GenericGraph.degree_to_cell` | Returns the number of edges from vertex to an edge in cell.
    :meth:`~GenericGraph.subgraph` | Returns the subgraph containing the given vertices and edges.
    :meth:`~GenericGraph.is_subgraph` | Tests whether self is a subgraph of other.

**Graph products:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.cartesian_product` | Returns the Cartesian product of self and other.
    :meth:`~GenericGraph.tensor_product` | Returns the tensor product, also called the categorical product, of self and other.
    :meth:`~GenericGraph.lexicographic_product` | Returns the lexicographic product of self and other.
    :meth:`~GenericGraph.strong_product` | Returns the strong product of self and other.
    :meth:`~GenericGraph.disjunctive_product` | Returns the disjunctive product of self and other.

**Paths and cycles:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.eulerian_orientation` | Returns a DiGraph which is an Eulerian orientation of the current graph.
    :meth:`~GenericGraph.eulerian_circuit` | Return a list of edges forming an eulerian circuit if one exists.
    :meth:`~GenericGraph.cycle_basis` | Returns a list of cycles which form a basis of the cycle space of ``self``.
    :meth:`~GenericGraph.interior_paths` | Returns an exhaustive list of paths (also lists) through only interior vertices from vertex start to vertex end in the (di)graph.
    :meth:`~GenericGraph.all_paths` | Returns a list of all paths (also lists) between a pair of vertices in the (di)graph.
    :meth:`~GenericGraph.triangles_count` | Returns the number of triangles in the (di)graph.

**Linear algebra:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.spectrum` | Returns a list of the eigenvalues of the adjacency matrix.
    :meth:`~GenericGraph.eigenvectors` | Returns the *right* eigenvectors of the adjacency matrix of the graph.
    :meth:`~GenericGraph.eigenspaces` | Returns the *right* eigenspaces of the adjacency matrix of the graph.

**Some metrics:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.cluster_triangles` | Returns the number of triangles for nbunch of vertices as a dictionary keyed by vertex.
    :meth:`~GenericGraph.clustering_average` | Returns the average clustering coefficient.
    :meth:`~GenericGraph.clustering_coeff` | Returns the clustering coefficient for each vertex in nbunch
    :meth:`~GenericGraph.cluster_transitivity` | Returns the transitivity (fraction of transitive triangles) of the graph.
    :meth:`~GenericGraph.szeged_index` | Returns the Szeged index of the graph.


**Automorphism group:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.coarsest_equitable_refinement` | Returns the coarsest partition which is finer than the input partition, and equitable with respect to self.
    :meth:`~GenericGraph.automorphism_group` | Returns the largest subgroup of the automorphism group of the (di)graph whose orbit partition is finer than the partition given.
    :meth:`~GenericGraph.is_vertex_transitive` | Returns whether the automorphism group of self is transitive within the partition provided
    :meth:`~GenericGraph.is_isomorphic` | Tests for isomorphism between self and other.
    :meth:`~GenericGraph.canonical_label` | Returns the unique graph on `\{0,1,...,n-1\}` ( ``n = self.order()`` ) which 1) is isomorphic to self 2) is invariant in the isomorphism class.

**Graph properties:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.is_eulerian` | Return true if the graph has a (closed) tour that visits each edge exactly once.
    :meth:`~GenericGraph.is_planar` | Tests whether the graph is planar.
    :meth:`~GenericGraph.is_circular_planar` | Tests whether the graph is circular planar (outerplanar)
    :meth:`~GenericGraph.is_regular` | Return ``True`` if this graph is (`k`-)regular.
    :meth:`~GenericGraph.is_chordal` | Tests whether the given graph is chordal.
    :meth:`~GenericGraph.is_circulant` | Tests whether the graph is a circulant graph.
    :meth:`~GenericGraph.is_interval` | Check whether self is an interval graph
    :meth:`~GenericGraph.is_gallai_tree` | Returns whether the current graph is a Gallai tree.
    :meth:`~GenericGraph.is_clique` | Tests whether a set of vertices is a clique
    :meth:`~GenericGraph.is_independent_set` | Tests whether a set of vertices is an independent set
    :meth:`~GenericGraph.is_transitively_reduced` | Tests whether the digraph is transitively reduced.
    :meth:`~GenericGraph.is_equitable` | Checks whether the given partition is equitable with respect to self.

**Traversals:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.breadth_first_search` | Returns an iterator over the vertices in a breadth-first ordering.
    :meth:`~GenericGraph.depth_first_search` | Returns an iterator over the vertices in a depth-first ordering.
    :meth:`~GenericGraph.lex_BFS` | Performs a Lex BFS on the graph.

**Distances:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.distance` | Returns the (directed) distance from u to v in the (di)graph
    :meth:`~GenericGraph.distance_all_pairs` | Returns the distances between all pairs of vertices.
    :meth:`~GenericGraph.distances_distribution` | Returns the distances distribution of the (di)graph in a dictionary.
    :meth:`~GenericGraph.eccentricity` | Return the eccentricity of vertex (or vertices) v.
    :meth:`~GenericGraph.radius` | Returns the radius of the (di)graph.
    :meth:`~GenericGraph.center` | Returns the set of vertices in the center of the graph
    :meth:`~GenericGraph.diameter` | Returns the largest distance between any two vertices.
    :meth:`~GenericGraph.distance_graph` | Returns the graph on the same vertex set as the original graph but vertices are adjacent in the returned graph if and only if they are at specified distances in the original graph.
    :meth:`~GenericGraph.girth` | Computes the girth of the graph.
    :meth:`~GenericGraph.periphery` | Returns the set of vertices in the periphery
    :meth:`~GenericGraph.shortest_path` | Returns a list of vertices representing some shortest path from `u` to `v`
    :meth:`~GenericGraph.shortest_path_length` | Returns the minimal length of paths from u to v
    :meth:`~GenericGraph.shortest_paths` | Returns a dictionary associating to each vertex v a shortest path from u to v, if it exists.
    :meth:`~GenericGraph.shortest_path_lengths` | Returns a dictionary of shortest path lengths keyed by targets that are connected by a path from u.
    :meth:`~GenericGraph.shortest_path_all_pairs` | Computes a shortest path between each pair of vertices.
    :meth:`~GenericGraph.wiener_index` | Returns the Wiener index of the graph.
    :meth:`~GenericGraph.average_distance` | Returns the average distance between vertices of the graph.


**Flows, connectivity, trees:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.is_connected` | Tests whether the (di)graph is connected.
    :meth:`~GenericGraph.connected_components` | Returns the list of connected components
    :meth:`~GenericGraph.connected_components_number` | Returns the number of connected components.
    :meth:`~GenericGraph.connected_components_subgraphs` | Returns a list of connected components as graph objects.
    :meth:`~GenericGraph.connected_component_containing_vertex` | Returns a list of the vertices connected to vertex.
    :meth:`~GenericGraph.blocks_and_cut_vertices` | Computes the blocks and cut vertices of the graph.
    :meth:`~GenericGraph.blocks_and_cuts_tree` | Computes the blocks-and-cuts tree of the graph.
    :meth:`~GenericGraph.is_cut_edge` | Returns True if the input edge is a cut-edge or a bridge.
    :meth:`~GenericGraph.is_cut_vertex` | Returns True if the input vertex is a cut-vertex.
    :meth:`~GenericGraph.edge_cut` | Returns a minimum edge cut between vertices `s` and `t`
    :meth:`~GenericGraph.vertex_cut` | Returns a minimum vertex cut between non-adjacent vertices `s` and `t`
    :meth:`~GenericGraph.flow` | Returns a maximum flow in the graph from ``x`` to ``y``
    :meth:`~GenericGraph.edge_disjoint_paths` | Returns a list of edge-disjoint paths between two vertices
    :meth:`~GenericGraph.vertex_disjoint_paths` | Returns a list of vertex-disjoint paths between two vertices as given by Menger's theorem.
    :meth:`~GenericGraph.edge_connectivity` | Returns the edge connectivity of the graph.
    :meth:`~GenericGraph.vertex_connectivity` | Returns the vertex connectivity of the graph.
    :meth:`~GenericGraph.transitive_closure` | Computes the transitive closure of a graph and returns it.
    :meth:`~GenericGraph.transitive_reduction` | Returns a transitive reduction of a graph.
    :meth:`~GenericGraph.min_spanning_tree` | Returns the edges of a minimum spanning tree.
    :meth:`~GenericGraph.spanning_trees_count` | Returns the number of spanning trees in a graph.

**Plot/embedding-related methods:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.set_embedding` | Sets a combinatorial embedding dictionary to ``_embedding`` attribute.
    :meth:`~GenericGraph.get_embedding` | Returns the attribute _embedding if it exists.
    :meth:`~GenericGraph.faces` | Returns the faces of an embedded graph.
    :meth:`~GenericGraph.get_pos` | Returns the position dictionary
    :meth:`~GenericGraph.set_pos` | Sets the position dictionary.
    :meth:`~GenericGraph.set_planar_positions` | Compute a planar layout for self using Schnyder's algorithm
    :meth:`~GenericGraph.layout_planar` | Uses Schnyder's algorithm to compute a planar layout for self.
    :meth:`~GenericGraph.is_drawn_free_of_edge_crossings` | Tests whether the position dictionary gives a planar embedding.
    :meth:`~GenericGraph.latex_options` | Returns an instance of :class:`~sage.graphs.graph_latex.GraphLatex` for the graph.
    :meth:`~GenericGraph.set_latex_options` | Sets multiple options for rendering a graph with LaTeX.
    :meth:`~GenericGraph.layout` | Returns a layout for the vertices of this graph.
    :meth:`~GenericGraph.layout_spring` | Computes a spring layout for this graph
    :meth:`~GenericGraph.layout_ranked` | Computes a ranked layout for this graph
    :meth:`~GenericGraph.layout_extend_randomly` | Extends randomly a partial layout
    :meth:`~GenericGraph.layout_circular` | Computes a circular layout for this graph
    :meth:`~GenericGraph.layout_tree` | Computes an ordered tree layout for this graph, which should be a tree (no non-oriented cycles).
    :meth:`~GenericGraph.layout_graphviz` | Calls ``graphviz`` to compute a layout of the vertices of this graph.
    :meth:`~GenericGraph.graphplot` | Returns a GraphPlot object.
    :meth:`~GenericGraph.plot` | Returns a graphics object representing the (di)graph.
    :meth:`~GenericGraph.show` | Shows the (di)graph.
    :meth:`~GenericGraph.plot3d` | Plot a graph in three dimensions.
    :meth:`~GenericGraph.show3d` | Plots the graph using Tachyon, and shows the resulting plot.
    :meth:`~GenericGraph.graphviz_string` | Returns a representation in the dot language.
    :meth:`~GenericGraph.graphviz_to_file_named` | Write a representation in the dot in a file.

**Algorithmically hard stuff:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraph.steiner_tree` | Returns a tree of minimum weight connecting the given set of vertices.
    :meth:`~GenericGraph.edge_disjoint_spanning_trees` | Returns the desired number of edge-disjoint spanning trees/arborescences.
    :meth:`~GenericGraph.feedback_vertex_set` | Computes the minimum feedback vertex set of a (di)graph.
    :meth:`~GenericGraph.multiway_cut` | Returns a minimum edge multiway cut
    :meth:`~GenericGraph.max_cut` | Returns a maximum edge cut of the graph.
    :meth:`~GenericGraph.longest_path` | Returns a longest path of ``self``.
    :meth:`~GenericGraph.traveling_salesman_problem` | Solves the traveling salesman problem (TSP)
    :meth:`~GenericGraph.is_hamiltonian` | Tests whether the current graph is Hamiltonian.
    :meth:`~GenericGraph.hamiltonian_cycle` | Returns a Hamiltonian cycle/circuit of the current graph/digraph
    :meth:`~GenericGraph.multicommodity_flow` | Solves a multicommodity flow problem.
    :meth:`~GenericGraph.disjoint_routed_paths` | Returns a set of disjoint routed paths.
    :meth:`~GenericGraph.dominating_set` | Returns a minimum dominating set of the graph
    :meth:`~GenericGraph.subgraph_search` | Returns a copy of ``G`` in ``self``.
    :meth:`~GenericGraph.subgraph_search_count` | Returns the number of labelled occurences of ``G`` in ``self``.
    :meth:`~GenericGraph.subgraph_search_iterator` | Returns an iterator over the labelled copies of ``G`` in ``self``.
    :meth:`~GenericGraph.characteristic_polynomial` | Returns the characteristic polynomial of the adjacency matrix of the (di)graph.
    :meth:`~GenericGraph.genus` | Returns the minimal genus of the graph.

Methods
-------
"""

from sage.misc.decorators import options
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import random
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from generic_graph_pyx import GenericGraph_pyx, spring_layout_fast
from sage.graphs.dot2tex_utils import assert_have_dot2tex
from sage.misc.superseded import deprecation, deprecated_function_alias

class GenericGraph(GenericGraph_pyx):
    """
    Base class for graphs and digraphs.
    """

    # Nice defaults for plotting arrays of graphs (see sage.misc.functional.show)
    graphics_array_defaults =  {'layout': 'circular', 'vertex_size':50, 'vertex_labels':False, 'graph_border':True}

    def __init__(self):
        r"""
        Every graph carries a dictionary of options, which is set
        here to ``None``.  Some options are added to the global
        :data:`sage.misc.latex.latex` instance which will insure
        that if LaTeX is used to render the graph,
        then the right packages are loaded and MathJax reacts
        properly.

        Most other initialization is done in the directed
        and undirected subclasses.

        TESTS::

            sage: g = Graph()
            sage: g
            Graph on 0 vertices
        """
        self._latex_opts = None

    def __add__(self, other_graph):
        """
        Returns the disjoint union of self and other.

        If there are common vertices to both, they will be renamed.

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G + H; J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4, 5, 6]
        """
        if isinstance(other_graph, GenericGraph):
            return self.disjoint_union(other_graph, verbose_relabel=False)

    def __eq__(self, other):
        """
        Comparison of self and other. For equality, must be in the same
        class, have the same settings for loops and multiedges, output the
        same vertex list (in order) and the same adjacency matrix.

        Note that this is _not_ an isomorphism test.

        EXAMPLES::

            sage: G = graphs.EmptyGraph()
            sage: H = Graph()
            sage: G == H
            True
            sage: G.to_directed() == H.to_directed()
            True
            sage: G = graphs.RandomGNP(8,.9999)
            sage: H = graphs.CompleteGraph(8)
            sage: G == H # most often true
            True
            sage: G = Graph( {0:[1,2,3,4,5,6,7]} )
            sage: H = Graph( {1:[0], 2:[0], 3:[0], 4:[0], 5:[0], 6:[0], 7:[0]} )
            sage: G == H
            True
            sage: G.allow_loops(True)
            sage: G == H
            False
            sage: G = graphs.RandomGNP(9,.3).to_directed()
            sage: H = graphs.RandomGNP(9,.3).to_directed()
            sage: G == H # most often false
            False
            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edge(0,1)
            sage: H = copy(G)
            sage: H.add_edge(0,1)
            sage: G == H
            False

        Note that graphs must be considered weighted, or Sage will not pay
        attention to edge label data in equality testing::

            sage: foo = Graph(sparse=True)
            sage: foo.add_edges([(0, 1, 1), (0, 2, 2)])
            sage: bar = Graph(sparse=True)
            sage: bar.add_edges([(0, 1, 2), (0, 2, 1)])
            sage: foo == bar
            True
            sage: foo.weighted(True)
            sage: foo == bar
            False
            sage: bar.weighted(True)
            sage: foo == bar
            False

        """
        # inputs must be (di)graphs:
        if not isinstance(other, GenericGraph):
            raise TypeError("cannot compare graph to non-graph (%s)"%str(other))
        from sage.graphs.all import Graph
        g1_is_graph = isinstance(self, Graph) # otherwise, DiGraph
        g2_is_graph = isinstance(other, Graph) # otherwise, DiGraph

        if (g1_is_graph != g2_is_graph):
            return False
        if self.allows_multiple_edges() != other.allows_multiple_edges():
            return False
        if self.allows_loops() != other.allows_loops():
            return False
        if self.order() != other.order():
            return False
        if self.size() != other.size():
            return False
        if any(x not in other for x in self):
            return False
        if self._weighted != other._weighted:
            return False

        verts = self.vertices()
        # Finally, we are prepared to check edges:
        if not self.allows_multiple_edges():
            for i in verts:
                for j in verts:
                    if self.has_edge(i,j) != other.has_edge(i,j):
                        return False
                    if self.has_edge(i,j) and self._weighted and other._weighted:
                        if self.edge_label(i,j) != other.edge_label(i,j):
                            return False
            return True
        else:
            for i in verts:
                for j in verts:
                    if self.has_edge(i, j):
                        edges1 = self.edge_label(i, j)
                    else:
                        edges1 = []
                    if other.has_edge(i, j):
                        edges2 = other.edge_label(i, j)
                    else:
                        edges2 = []
                    if len(edges1) != len(edges2):
                        return False
                    if sorted(edges1) != sorted(edges2) and self._weighted and other._weighted:
                        return False
            return True

    @cached_method
    def __hash__(self):
        """
        Only immutable graphs are hashable.

        The hash value of an immutable graph relies on the tuple
        of vertices and the tuple of edges. The resulting value
        is cached.

        EXAMPLE::

            sage: G = graphs.PetersenGraph()
            sage: {G:1}[G]
            Traceback (most recent call last):
            ...
            TypeError: This graph is mutable, and thus not hashable. Create
            an immutable copy by `g.copy(immutable=True)`
            sage: G_imm = Graph(G, data_structure="static_sparse")
            sage: G_imm == G
            True
            sage: {G_imm:1}[G_imm]  # indirect doctest
            1
            sage: G_imm.__hash__() is G_imm.__hash__()
            True

        """
        if getattr(self, "_immutable", False):
            return hash((tuple(self.vertices()), tuple(self.edges())))
        raise TypeError("This graph is mutable, and thus not hashable. "
                        "Create an immutable copy by `g.copy(immutable=True)`")

    def __mul__(self, n):
        """
        Returns the sum of a graph with itself n times.

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = G*3; H
            Cycle graph disjoint_union Cycle graph disjoint_union Cycle graph: Graph on 9 vertices
            sage: H.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8]
            sage: H = G*1; H
            Cycle graph: Graph on 3 vertices
        """
        if isinstance(n, (int, long, Integer)):
            if n < 1:
                raise TypeError('multiplication of a graph and a nonpositive integer is not defined')
            if n == 1:
                from copy import copy
                return copy(self)
            return sum([self]*(n-1), self)
        else:
            raise TypeError('multiplication of a graph and something other than an integer is not defined')

    def __ne__(self, other):
        """
        Tests for inequality, complement of __eq__.

        EXAMPLES::

            sage: g = Graph()
            sage: g2 = copy(g)
            sage: g == g
            True
            sage: g != g
            False
            sage: g2 == g
            True
            sage: g2 != g
            False
            sage: g is g
            True
            sage: g2 is g
            False
        """
        return (not (self == other))

    def __rmul__(self, n):
        """
        Returns the sum of a graph with itself n times.

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = int(3)*G; H
            Cycle graph disjoint_union Cycle graph disjoint_union Cycle graph: Graph on 9 vertices
            sage: H.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8]
        """
        return self*n

    def __str__(self):
        """
        str(G) returns the name of the graph, unless the name is the empty
        string, in which case it returns the default representation.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: str(G)
            'Petersen graph'
        """
        if self.name():
            return self.name()
        else:
            return repr(self)

    def _bit_vector(self):
        """
        Returns a string representing the edges of the (simple) graph for
        graph6 and dig6 strings.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G._bit_vector()
            '101001100110000010000001001000010110000010110'
            sage: len([a for a in G._bit_vector() if a == '1'])
            15
            sage: G.num_edges()
            15
        """
        self._scream_if_not_simple()
        vertices = self.vertices()
        n = len(vertices)
        if self._directed:
            total_length = n*n
            bit = lambda x,y : x*n + y
        else:
            total_length = (n*(n - 1))//2
            n_ch_2 = lambda b : int(b*(b-1))//2
            bit = lambda x,y : n_ch_2(max([x,y])) + min([x,y])
        bit_vector = set()
        for u,v,_ in self.edge_iterator():
            bit_vector.add(bit(vertices.index(u), vertices.index(v)))
        bit_vector = sorted(bit_vector)
        s = []
        j = 0
        for i in bit_vector:
            s.append( '0'*(i - j) + '1' )
            j = i + 1
        s = "".join(s)
        s += '0'*(total_length-len(s))
        return s

    def _latex_(self):
        r"""

        Returns a string to render the graph using LaTeX.

        To adjust the string, use the
        :meth:`set_latex_options` method to set options,
        or call the :meth:`latex_options` method to
        get a :class:`~sage.graphs.graph_latex.GraphLatex`
        object that may be used to also customize the
        output produced here.  Possible options are documented at
        :meth:`sage.graphs.graph_latex.GraphLatex.set_option`.

        EXAMPLES::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on TeX installation
            sage: g = graphs.CompleteGraph(2)
            sage: print g._latex_()
            \begin{tikzpicture}
            %
            \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
            %
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
            \definecolor{clv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v1}{rgb}{0.0,0.0,0.0}
            %
            \Vertex[style={minimum size=1.0cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$0$},x=5.0cm,y=5.0cm]{v0}
            \Vertex[style={minimum size=1.0cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$1$},x=0.0cm,y=0.0cm]{v1}
            %
            \Edge[lw=0.1cm,style={color=cv0v1,},](v0)(v1)
            %
            \end{tikzpicture}
        """
        from sage.graphs.graph_latex import setup_latex_preamble
        setup_latex_preamble()

        return self.latex_options().latex()

    def _matrix_(self, R=None):
        """
        Returns the adjacency matrix of the graph over the specified ring.

        EXAMPLES::

            sage: G = graphs.CompleteBipartiteGraph(2,3)
            sage: m = matrix(G); m.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Integer Ring
            sage: m
            [0 0 1 1 1]
            [0 0 1 1 1]
            [1 1 0 0 0]
            [1 1 0 0 0]
            [1 1 0 0 0]
            sage: G._matrix_()
            [0 0 1 1 1]
            [0 0 1 1 1]
            [1 1 0 0 0]
            [1 1 0 0 0]
            [1 1 0 0 0]
            sage: factor(m.charpoly())
            x^3 * (x^2 - 6)
        """
        if R is None:
            return self.am()
        else:
            return self.am().change_ring(R)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G._repr_()
            'Petersen graph: Graph on 10 vertices'
        """
        name = ""
        if self.allows_loops():
            name += "looped "
        if self.allows_multiple_edges():
            name += "multi-"
        if self._directed:
            name += "di"
        name += "graph on %d vert"%self.order()
        if self.order() == 1:
            name += "ex"
        else:
            name += "ices"
        name = name.capitalize()
        if self.name() != '':
            name = self.name() + ": " + name
        return name

    ### Formats

    def __copy__(self, implementation='c_graph', data_structure=None,
                 sparse=None, immutable=None):
        """
        Return a copy of the graph.

        INPUT:

         - ``implementation`` - string (default: 'c_graph') the implementation
           goes here.  Current options are only 'networkx' or 'c_graph'.

         - ``sparse`` (boolean) -- ``sparse=True`` is an alias for
           ``data_structure="sparse"``, and ``sparse=False`` is an alias for
           ``data_structure="dense"``. Only used when
           ``implementation='c_graph'`` and ``data_structure=None``.

         - ``data_structure`` -- one of ``"sparse"``, ``"static_sparse"``, or
           ``"dense"``. See the documentation of :class:`Graph` or
           :class:`DiGraph`. Only used when ``implementation='c_graph'``.

         - ``immutable`` (boolean) -- whether to create a mutable/immutable
           copy. Only used when ``implementation='c_graph'`` and
           ``data_structure=None``.

           * ``immutable=None`` (default) means that the graph and its copy will
             behave the same way.

           * ``immutable=True`` is a shortcut for
             ``data_structure='static_sparse'`` and ``implementation='c_graph'``

           * ``immutable=False`` sets ``implementation`` to ``'c_graph'``. When
             ``immutable=False`` is used to copy an immutable graph, the data
             structure used is ``"sparse"`` unless anything else is specified.

        .. NOTE::

            If the graph uses
            :class:`~sage.graphs.base.static_sparse_backend.StaticSparseBackend`
            and the ``_immutable`` flag, then ``self`` is returned rather than a
            copy (unless one of the optional arguments is used).

        OUTPUT:

        A Graph object.

        .. warning::

           Please use this method only if you need to copy but change the
           underlying implementation.  Otherwise simply do ``copy(g)``
           instead of ``g.copy()``.

        EXAMPLES::

            sage: g=Graph({0:[0,1,1,2]},loops=True,multiedges=True,sparse=True)
            sage: g==copy(g)
            True
            sage: g=DiGraph({0:[0,1,1,2],1:[0,1]},loops=True,multiedges=True,sparse=True)
            sage: g==copy(g)
            True

        Note that vertex associations are also kept::

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: T = graphs.TetrahedralGraph()
            sage: T.set_vertices(d)
            sage: T2 = copy(T)
            sage: T2.get_vertex(0)
            Dodecahedron: Graph on 20 vertices

        Notice that the copy is at least as deep as the objects::

            sage: T2.get_vertex(0) is T.get_vertex(0)
            False

        Examples of the keywords in use::

            sage: G = graphs.CompleteGraph(19)
            sage: H = G.copy(implementation='c_graph')
            sage: H == G; H is G
            True
            False
            sage: G1 = G.copy(sparse=True)
            sage: G1==G
            True
            sage: G1 is G
            False
            sage: G2 = copy(G)
            sage: G2 is G
            False

        TESTS:

        We make copies of the ``_pos`` and ``_boundary`` attributes::

            sage: g = graphs.PathGraph(3)
            sage: h = copy(g)
            sage: h._pos is g._pos
            False
            sage: h._boundary is g._boundary
            False

        We make sure that one can make immutable copies by providing the
        ``data_structure`` optional argument, and that copying an immutable graph
        returns the graph::

            sage: G = graphs.PetersenGraph()
            sage: hash(G)
            Traceback (most recent call last):
            ...
            TypeError: This graph is mutable, and thus not hashable. Create an
            immutable copy by `g.copy(immutable=True)`
            sage: g = G.copy(immutable=True)
            sage: hash(g)    # random
            1833517720
            sage: g==G
            True
            sage: g is copy(g) is g.copy()
            True

        ``immutable=True`` is a short-cut for ``data_structure='static_sparse'``::

            sage: g is g.copy(data_structure='static_sparse') is g.copy(immutable=True)
            True

        If a graph pretends to be immutable, but does not use the static sparse
        backend, then the copy is not identic with the graph, even though it is
        considered to be hashable::

            sage: P = Poset(([1,2,3,4], [[1,3],[1,4],[2,3]]), linear_extension=True, facade = False)
            sage: H = P.hasse_diagram()
            sage: H._immutable = True
            sage: hash(H)   # random
            -1843552882
            sage: copy(H) is H
            False

        TESTS:

        Bad input::

            sage: G.copy(data_structure="sparse", sparse=False)
            Traceback (most recent call last):
            ...
            ValueError: You cannot define 'immutable' or 'sparse' when 'data_structure' has a value.
            sage: G.copy(data_structure="sparse", immutable=True)
            Traceback (most recent call last):
            ...
            ValueError: You cannot define 'immutable' or 'sparse' when 'data_structure' has a value.
            sage: G.copy(immutable=True, sparse=False)
            Traceback (most recent call last):
            ...
            ValueError: There is no dense immutable backend at the moment.

        Which backend ?::

            sage: G.copy(data_structure="sparse")._backend
            <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
            sage: G.copy(data_structure="dense")._backend
            <class 'sage.graphs.base.dense_graph.DenseGraphBackend'>
            sage: G.copy(data_structure="static_sparse")._backend
            <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
            sage: G.copy(immutable=True)._backend
            <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
            sage: G.copy(immutable=True, sparse=True)._backend
            <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
            sage: G.copy(immutable=False, sparse=True)._backend
            <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
            sage: G.copy(immutable=False, sparse=False)._backend
            <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
            sage: Graph(implementation="networkx").copy(implementation='c_graph')._backend
            <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>

        Fake immutable graphs::

            sage: G._immutable = True
            sage: G.copy()._backend
            <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
        """
        # Which data structure should be used ?
        if implementation != 'c_graph':
            # We do not care about the value of data_structure. But let's check
            # the user did not define too much.
            if data_structure is not None or immutable is not None or sparse is not None:
                raise ValueError("'data_structure' 'immutable' and 'sparse' can"
                                 " only be defined when 'implementation'='c_graph'")
        elif data_structure is not None:
            # data_structure is already defined so there is nothing left to do
            # here ! Did the user try to define too much ?
            if immutable is not None or sparse is not None:
                raise ValueError("You cannot define 'immutable' or 'sparse' "
                                 "when 'data_structure' has a value.")
        # At this point :
        # - implementation is 'c_graph'
        # - data_structure is None.
        elif immutable is True:
            data_structure = 'static_sparse'
            if sparse is False:
                raise ValueError("There is no dense immutable backend at the moment.")
        elif immutable is False:
            # If the users requests a mutable graph and input is immutable, we
            # chose the 'sparse' cgraph backend. Unless the user explicitly
            # asked for something different.
            if getattr(self, '_immutable', False):
                data_structure = 'dense' if sparse is False else 'sparse'
        elif sparse is True:
            data_structure = "sparse"
        elif sparse is False:
            data_structure = "dense"

        # Immutable copy of an immutable graph ? return self !
        if getattr(self, '_immutable', False):
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            if (isinstance(self._backend, StaticSparseBackend) and
                implementation=='c_graph' and
                (data_structure=='static_sparse' or data_structure is None)):
                return self

        if data_structure is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            if isinstance(self._backend, DenseGraphBackend):
                data_structure = "dense"
            else:
                data_structure = "sparse"

        from copy import copy
        G = self.__class__(self, name=self.name(), pos=copy(self._pos), boundary=copy(self._boundary), implementation=implementation, data_structure=data_structure)

        attributes_to_copy = ('_assoc', '_embedding')
        for attr in attributes_to_copy:
            if hasattr(self, attr):
                copy_attr = {}
                old_attr = getattr(self, attr)
                if isinstance(old_attr, dict):
                    for v,value in old_attr.iteritems():
                        try:
                            copy_attr[v] = value.copy()
                        except AttributeError:
                            from copy import copy
                            copy_attr[v] = copy(value)
                    setattr(G, attr, copy_attr)
                else:
                    setattr(G, attr, copy(old_attr))

        G._weighted = self._weighted
        return G

    copy = __copy__

    def _scream_if_not_simple(self, allow_loops=False, allow_multiple_edges=False):
        r"""
        Raises an exception if the graph is not simple.

        This function is called by some functions of the Graph library when they
        have been written for simple graphs only (i.e. neither loops nor multiple
        edges). It raises an exception inviting the user to convert the graph to
        a simple graph first, before calling the function again.

        Note that this function does not check the existence of loops or
        multiple edges, which would take linear time: it merely checks that the
        graph *does not allow* multiple edges nor loops, which takes a constant
        time.

        INPUT:

        - ``allow_loops`` (boolean) -- whether to tolerate loops. Set to
          ``False`` by default.

        - ``allow_multiple_edges`` (boolean) -- whether to tolerate multiple
          edges. Set to ``False`` by default.

        .. SEEALSO::

            * :meth:`allow_loops`
            * :meth:`allow_multiple_edges`

        EXAMPLES::

            sage: g = graphs.PetersenGraph()

        No scream::

            sage: from itertools import product
            sage: for p,q in product((True,False),repeat=2):
            ....:     g.allow_loops(p)
            ....:     g.allow_multiple_edges(q)
            ....:     g._scream_if_not_simple(p,q)

        A lot of them::

            sage: g.allow_loops(True); g.allow_multiple_edges(True)
            sage: g._scream_if_not_simple()
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with multiedges/loops. Perhaps this method can be updated to handle them, but in the meantime if you want to use it please disallow multiedges/loops using allow_multiple_edges() and allow_loops().
            sage: g._scream_if_not_simple(allow_loops=True)
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with multiedges. Perhaps this method can be updated to handle them, but in the meantime if you want to use it please disallow multiedges using allow_multiple_edges().
            sage: g._scream_if_not_simple(allow_multiple_edges=True)
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with loops. Perhaps this method can be updated to handle them, but in the meantime if you want to use it please disallow loops using allow_loops().

            sage: g.allow_loops(True); g.allow_multiple_edges(False)
            sage: g._scream_if_not_simple()
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with loops. Perhaps this method can be updated to handle them, but in the meantime if you want to use it please disallow loops using allow_loops().
            sage: g._scream_if_not_simple(allow_multiple_edges=True)
            Traceback (most recent call last):
            ...
            ValueError: This method is not known to work on graphs with loops. Perhaps this method can be updated to handle them, but in the meantime if you want to use it please disallow loops using allow_loops().
        """
        pb_with_loops = not allow_loops and self.allows_loops()
        pb_with_multiple_edges = not allow_multiple_edges and self.allows_multiple_edges()
        if pb_with_loops or pb_with_multiple_edges:
            if pb_with_loops and pb_with_multiple_edges:
                name = "multiedges/loops"
                functions = "allow_multiple_edges() and allow_loops()"
            elif pb_with_loops:
                name = "loops"
                functions = "allow_loops()"
            elif pb_with_multiple_edges:
                name = "multiedges"
                functions = "allow_multiple_edges()"
            msg = ("This method is not known to work on graphs with "+name+". "+
                   "Perhaps this method can be updated to handle them, but in the "+
                   "meantime if you want to use it please disallow "+name+" using "+
                   functions+".")
            raise ValueError(msg)

    def networkx_graph(self, copy=True):
        """
        Creates a new NetworkX graph from the Sage graph.

        INPUT:


        -  ``copy`` - if False, and the underlying
           implementation is a NetworkX graph, then the actual object itself
           is returned.


        EXAMPLES::

            sage: G = graphs.TetrahedralGraph()
            sage: N = G.networkx_graph()
            sage: type(N)
            <class 'networkx.classes.graph.Graph'>

        ::

            sage: G = graphs.TetrahedralGraph()
            sage: G = Graph(G, implementation='networkx')
            sage: N = G.networkx_graph()
            sage: G._backend._nxg is N
            False

        ::

            sage: G = Graph(graphs.TetrahedralGraph(), implementation='networkx')
            sage: N = G.networkx_graph(copy=False)
            sage: G._backend._nxg is N
            True
        """
        try:
            if copy:
                return self._backend._nxg.copy()
            else:
                return self._backend._nxg
        except Exception:
            import networkx
            if self._directed and self.allows_multiple_edges():
                class_type = networkx.MultiDiGraph
            elif self._directed:
                class_type = networkx.DiGraph
            elif self.allows_multiple_edges():
                class_type = networkx.MultiGraph
            else:
                class_type = networkx.Graph
            N = class_type(selfloops=self.allows_loops(), multiedges=self.allows_multiple_edges(),
                           name=self.name())
            N.add_nodes_from(self.vertices())
            for u,v,l in self.edges():
                if l is None:
                    N.add_edge(u,v)
                else:
                    from networkx import NetworkXError
                    try:
                        N.add_edge(u,v,l)
                    except (TypeError, ValueError, NetworkXError):
                        N.add_edge(u,v,weight=l)
            return N

    def to_dictionary(self, edge_labels=False, multiple_edges=False):
        r"""
        Returns the graph as a dictionary.

        INPUT:

        - ``edge_labels`` (boolean) -- whether to include edge labels in the
          output.

        - ``multiple_edges`` (boolean) -- whether to include multiple edges in
          the output.

        OUTPUT:

        The output depends on the input:

        * If ``edge_labels == False`` and ``multiple_edges == False``, the
          output is a dictionary associating to each vertex the list of its
          neighbors.

        * If ``edge_labels == False`` and ``multiple_edges == True``, the output
          is a dictionary the same as previously with one difference: the
          neighbors are listed with multiplicity.

        * If ``edge_labels == True`` and ``multiple_edges == False``, the output
          is a dictionary associating to each vertex `u` [a dictionary
          associating to each vertex `v` incident to `u` the label of edge
          `(u,v)`].

        * If ``edge_labels == True`` and ``multiple_edges == True``, the output
          is a dictionary associating to each vertex `u` [a dictionary
          associating to each vertex `v` incident to `u` [the list of labels of
          all edges between `u` and `v`]].

        .. NOTE::

          When used on directed graphs, the explanations above can be understood
          by replacing the word "neigbours" by "out-neighbors"

        EXAMPLES::

            sage: g = graphs.PetersenGraph().to_dictionary()
            sage: [(key, sorted(g[key])) for key in g]
            [(0, [1, 4, 5]),
             (1, [0, 2, 6]),
             (2, [1, 3, 7]),
             (3, [2, 4, 8]),
             (4, [0, 3, 9]),
             (5, [0, 7, 8]),
             (6, [1, 8, 9]),
             (7, [2, 5, 9]),
             (8, [3, 5, 6]),
             (9, [4, 6, 7])]
            sage: graphs.PetersenGraph().to_dictionary(multiple_edges=True)
            {0: [1, 4, 5], 1: [0, 2, 6],
             2: [1, 3, 7], 3: [2, 4, 8],
             4: [0, 3, 9], 5: [0, 7, 8],
             6: [1, 8, 9], 7: [2, 5, 9],
             8: [3, 5, 6], 9: [4, 6, 7]}
            sage: graphs.PetersenGraph().to_dictionary(edge_labels=True)
            {0: {1: None, 4: None, 5: None},
             1: {0: None, 2: None, 6: None},
             2: {1: None, 3: None, 7: None},
             3: {8: None, 2: None, 4: None},
             4: {0: None, 9: None, 3: None},
             5: {0: None, 8: None, 7: None},
             6: {8: None, 1: None, 9: None},
             7: {9: None, 2: None, 5: None},
             8: {3: None, 5: None, 6: None},
             9: {4: None, 6: None, 7: None}}
            sage: graphs.PetersenGraph().to_dictionary(edge_labels=True,multiple_edges=True)
            {0: {1: [None], 4: [None], 5: [None]},
             1: {0: [None], 2: [None], 6: [None]},
             2: {1: [None], 3: [None], 7: [None]},
             3: {8: [None], 2: [None], 4: [None]},
             4: {0: [None], 9: [None], 3: [None]},
             5: {0: [None], 8: [None], 7: [None]},
             6: {8: [None], 1: [None], 9: [None]},
             7: {9: [None], 2: [None], 5: [None]},
             8: {3: [None], 5: [None], 6: [None]},
             9: {4: [None], 6: [None], 7: [None]}}
        """

        # Returning the resuls as a dictionary of lists
        #
        # dictionary :
        # {vertex : [list of (out-)neighbors]}

        if not edge_labels and not multiple_edges:
            d = {}

            if self.is_directed():
                for u in self:
                    d[u]=self.neighbors_out(u)
            else:
                for u in self:
                    d[u]=self.neighbors(u)


        # Returning the result as a dictionary of lists
        #
        # dictionary :
        # {vertex : [list of (out-)neighbors, with multiplicity]}
        elif not edge_labels and multiple_edges:
            d={v:[] for v in self}

            if self.is_directed():
                for u,v in self.edge_iterator(labels = False):
                    d[u].append(v)

            else:
                for u,v in self.edge_iterator(labels = False):
                    d[u].append(v)
                    d[v].append(u)

        # Returning the result as a dictionary of dictionaries
        #
        # Each vertex is associated with the dictionary associating to each of
        # its neighbors the corresponding edge label.
        #
        # dictionary :
        # {v : dictionary                          }
        #      {neighbor u of v : label of edge u,v}

        elif edge_labels and not multiple_edges:
            d={v:{} for v in self}

            if self.is_directed():
                for u,v,l in self.edge_iterator():
                    d[u][v] = l

            else:
                for u,v,l in self.edge_iterator():
                    d[u][v] = l
                    d[v][u] = l

        # Returning the result as a dictionary of dictionaries
        #
        # Each vertex is associated with the dictionary associating to each of
        # its neighbors the list of edge labels between the two vertices
        #
        # dictionary :
        # {v : dictionary                                          }
        #      {neighbor u of v : [labels of edges between u and v]}

        elif edge_labels and multiple_edges:
            d={v:{} for v in self}

            if self.is_directed():
                for u,v,l in self.edge_iterator():
                    if not v in d[u]:
                        d[u][v] = []
                    d[u][v].append(l)

            else:
                for u,v,l in self.edge_iterator():
                    if not v in d[u]:
                        d[u][v] = []
                        d[v][u] = []

                    d[u][v].append(l)
                    d[v][u].append(l)

        return d

    def adjacency_matrix(self, sparse=None, boundary_first=False):
        """
        Returns the adjacency matrix of the (di)graph.

        Each vertex is represented by its position in the list returned by the
        vertices() function.

        The matrix returned is over the integers. If a different ring is
        desired, use either the change_ring function or the matrix
        function.

        INPUT:

        -  ``sparse`` - whether to represent with a sparse
           matrix

        -  ``boundary_first`` - whether to represent the
           boundary vertices in the upper left block

        EXAMPLES::

            sage: G = graphs.CubeGraph(4)
            sage: G.adjacency_matrix()
            [0 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0]
            [1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0]
            [0 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0]
            [1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0]
            [0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0]
            [0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]
            [0 1 0 0 0 0 0 0 1 0 0 1 0 1 0 0]
            [0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1]
            [0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 1]
            [0 0 0 0 0 0 0 1 0 0 0 1 0 1 1 0]

        ::

            sage: matrix(GF(2),G) # matrix over GF(2)
            [0 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0]
            [1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0]
            [0 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0]
            [1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0]
            [0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0]
            [0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]
            [0 1 0 0 0 0 0 0 1 0 0 1 0 1 0 0]
            [0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1]
            [0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 1]
            [0 0 0 0 0 0 0 1 0 0 0 1 0 1 1 0]

        ::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.adjacency_matrix()
            [0 1 1 1 0 0]
            [1 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [1 0 0 0 0 1]
            [0 1 0 0 0 0]

        TESTS::

            sage: graphs.CubeGraph(8).adjacency_matrix().parent()
            Full MatrixSpace of 256 by 256 dense matrices over Integer Ring
            sage: graphs.CubeGraph(9).adjacency_matrix().parent()
            Full MatrixSpace of 512 by 512 sparse matrices over Integer Ring
            sage: Graph([(i,i+1) for i in range(500)]+[(0,1),], multiedges=True).adjacency_matrix().parent()
            Full MatrixSpace of 501 by 501 dense matrices over Integer Ring
        """
        n = self.order()
        if sparse is None:
            sparse=True
            if self.has_multiple_edges() or n <= 256 or self.density() > 0.05:
                sparse=False
        verts = self.vertices(boundary_first=boundary_first)
        new_indices = dict((v,i) for i,v in enumerate(verts))
        D = {}
        directed = self._directed
        multiple_edges = self.allows_multiple_edges()
        for i,j,l in self.edge_iterator():
            i = new_indices[i]
            j = new_indices[j]
            if multiple_edges and (i,j) in D:
                D[(i,j)] += 1
                if not directed and i != j:
                    D[(j,i)] += 1
            else:
                D[(i,j)] = 1
                if not directed and i != j:
                    D[(j,i)] = 1
        from sage.matrix.constructor import matrix
        M = matrix(ZZ, n, n, D, sparse=sparse)
        return M

    am = adjacency_matrix # shorter call makes life easier

    def incidence_matrix(self, sparse=True):
        """
        Returns the incidence matrix of the (di)graph.

        Each row is a vertex, and each column is an edge. Note that in the case
        of graphs, there is a choice of orientation for each edge.

        EXAMPLES::

            sage: G = graphs.CubeGraph(3)
            sage: G.incidence_matrix()
            [-1 -1 -1  0  0  0  0  0  0  0  0  0]
            [ 0  0  1 -1 -1  0  0  0  0  0  0  0]
            [ 0  1  0  0  0 -1 -1  0  0  0  0  0]
            [ 0  0  0  0  1  0  1 -1  0  0  0  0]
            [ 1  0  0  0  0  0  0  0 -1 -1  0  0]
            [ 0  0  0  1  0  0  0  0  0  1 -1  0]
            [ 0  0  0  0  0  1  0  0  1  0  0 -1]
            [ 0  0  0  0  0  0  0  1  0  0  1  1]


        A well known result states that the product of the incidence matrix
        with its transpose is in fact the Kirchhoff matrix::

            sage: G = graphs.PetersenGraph()
            sage: G.incidence_matrix()*G.incidence_matrix().transpose() == G.kirchhoff_matrix()
            True

        ::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incidence_matrix()
            [-1 -1 -1  0  0  0  0  0  1  1]
            [ 0  0  1 -1  0  0  0  1 -1  0]
            [ 0  1  0  1 -1  0  0  0  0  0]
            [ 1  0  0  0  1 -1  0  0  0  0]
            [ 0  0  0  0  0  1 -1  0  0 -1]
            [ 0  0  0  0  0  0  1 -1  0  0]
        """
        from sage.matrix.constructor import matrix
        from copy import copy
        n = self.order()
        verts = self.vertices()
        d = [0]*n
        cols = []
        if self._directed:
            for i, j, l in self.edge_iterator():
                col = copy(d)
                i = verts.index(i)
                j = verts.index(j)
                col[i] = -1
                col[j] = 1
                cols.append(col)
        else:
            for i, j, l in self.edge_iterator():
                col = copy(d)
                i,j = (i,j) if i <= j else (j,i)
                i = verts.index(i)
                j = verts.index(j)
                col[i] = -1
                col[j] = 1
                cols.append(col)
        cols.sort()
        return matrix(cols, sparse=sparse).transpose()

    def distance_matrix(self):
        """
        Returns the distance matrix of the (strongly) connected (di)graph.

        The distance matrix of a (strongly) connected (di)graph is a matrix whose
        rows and columns are indexed with the vertices of the (di) graph. The
        intersection of a row and column contains the respective distance between
        the vertices indexed at these position.

        .. WARNING::

            The ordering of vertices in the matrix has no reason to correspond
            to the order of vertices in
            :meth:`~sage.graphs.generic_graph.GenericGraph.vertices`. In
            particular, if two integers `i,j` are vertices of a graph `G` with
            distance matrix ``M``, then ``M[i][i]`` is not necessarily the
            distance between vertices `i` and `j`.

        EXAMPLES::

            sage: G = graphs.CubeGraph(3)
            sage: G.distance_matrix()
            [0 1 1 2 1 2 2 3]
            [1 0 2 1 2 1 3 2]
            [1 2 0 1 2 3 1 2]
            [2 1 1 0 3 2 2 1]
            [1 2 2 3 0 1 1 2]
            [2 1 3 2 1 0 2 1]
            [2 3 1 2 1 2 0 1]
            [3 2 2 1 2 1 1 0]

        The well known result of Graham and Pollak states that the determinant of
        the distance matrix of any tree of order n is (-1)^{n-1}(n-1)2^{n-2} ::

            sage: all(T.distance_matrix().det() == (-1)^9*(9)*2^8 for T in graphs.trees(10))
            True

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.distance_all_pairs`
              -- computes the distance between any two vertices.
        """

        from sage.matrix.constructor import matrix

        n = self.order()
        ret = matrix(n,n)
        V = self.vertices()

        dist = self.distance_all_pairs()

        for i in xrange(n):
            for j in xrange(i+1,n):
                d = (dist[V[i]])[V[j]]
                if d > n :
                    raise ValueError("Input (di)graph must be (strongly) connected.")
                ret[i,j] = ret[j,i] = d

        return ret

    def weighted_adjacency_matrix(self, sparse=True, boundary_first=False):
        """
        Returns the weighted adjacency matrix of the graph.

        Each vertex is represented by its position in the list returned by the
        vertices() function.

        EXAMPLES::

            sage: G = Graph(sparse=True, weighted=True)
            sage: G.add_edges([(0,1,1),(1,2,2),(0,2,3),(0,3,4)])
            sage: M = G.weighted_adjacency_matrix(); M
            [0 1 3 4]
            [1 0 2 0]
            [3 2 0 0]
            [4 0 0 0]
            sage: H = Graph(data=M, format='weighted_adjacency_matrix', sparse=True)
            sage: H == G
            True

        The following doctest verifies that \#4888 is fixed::

            sage: G = DiGraph({0:{}, 1:{0:1}, 2:{0:1}}, weighted = True,sparse=True)
            sage: G.weighted_adjacency_matrix()
            [0 0 0]
            [1 0 0]
            [1 0 0]

        """
        if self.has_multiple_edges():
            raise NotImplementedError("don't know how to represent weights for a multigraph")

        verts = self.vertices(boundary_first=boundary_first)
        new_indices = dict((v,i) for i,v in enumerate(verts))

        D = {}
        if self._directed:
            for i,j,l in self.edge_iterator():
                i = new_indices[i]
                j = new_indices[j]
                D[(i,j)] = l
        else:
            for i,j,l in self.edge_iterator():
                i = new_indices[i]
                j = new_indices[j]
                D[(i,j)] = l
                D[(j,i)] = l
        from sage.matrix.constructor import matrix
        M = matrix(self.num_verts(), D, sparse=sparse)
        return M

    def kirchhoff_matrix(self, weighted=None, indegree=True, normalized=False, **kwds):
        """
        Returns the Kirchhoff matrix (a.k.a. the Laplacian) of the graph.

        The Kirchhoff matrix is defined to be `D - M`, where `D` is
        the diagonal degree matrix (each diagonal entry is the degree
        of the corresponding vertex), and `M` is the adjacency matrix.
        If ``normalized`` is ``True``, then the returned matrix is
        `D^{-1/2}(D-M)D^{-1/2}`.

        ( In the special case of DiGraphs, `D` is defined as the diagonal
        in-degree matrix or diagonal out-degree matrix according to the
        value of ``indegree``)

        INPUT:

        - ``weighted`` -- Binary variable :
            - If ``True``, the weighted adjacency matrix is used for `M`,
              and the diagonal matrix `D` takes into account the weight of edges
              (replace in the definition "degree" by "sum of the incident edges" ).
            - Else, each edge is assumed to have weight 1.

            Default is to take weights into consideration if and only if the graph is
            weighted.

        - ``indegree`` -- Binary variable  :
            - If ``True``, each diagonal entry of `D` is equal to the
              in-degree of the corresponding vertex.
            - Else, each diagonal entry of `D` is equal to the
              out-degree of the corresponding vertex.

              By default, ``indegree`` is set to ``True``

            ( This variable only matters when the graph is a digraph )

        - ``normalized`` -- Binary variable :

            - If ``True``, the returned matrix is
              `D^{-1/2}(D-M)D^{-1/2}`, a normalized version of the
              Laplacian matrix.
              (More accurately, the normalizing matrix used is equal to `D^{-1/2}`
              only for non-isolated vertices.  If vertex `i` is isolated, then
              diagonal entry `i` in the matrix is 1, rather than a division by
              zero.)
            - Else, the matrix `D-M` is returned

        Note that any additional keywords will be passed on to either
        the ``adjacency_matrix`` or ``weighted_adjacency_matrix`` method.

        AUTHORS:

        - Tom Boothby
        - Jason Grout

        EXAMPLES::

            sage: G = Graph(sparse=True)
            sage: G.add_edges([(0,1,1),(1,2,2),(0,2,3),(0,3,4)])
            sage: M = G.kirchhoff_matrix(weighted=True); M
            [ 8 -1 -3 -4]
            [-1  3 -2  0]
            [-3 -2  5  0]
            [-4  0  0  4]
            sage: M = G.kirchhoff_matrix(); M
            [ 3 -1 -1 -1]
            [-1  2 -1  0]
            [-1 -1  2  0]
            [-1  0  0  1]
            sage: G.set_boundary([2,3])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: M = G.kirchhoff_matrix(weighted=True, boundary_first=True); M
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            [ 5  0 -3 -2]
            [ 0  4 -4  0]
            [-3 -4  8 -1]
            [-2  0 -1  3]
            sage: M = G.kirchhoff_matrix(boundary_first=True); M
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            [ 2  0 -1 -1]
            [ 0  1 -1  0]
            [-1 -1  3 -1]
            [-1  0 -1  2]
            sage: M = G.laplacian_matrix(boundary_first=True); M
            [ 2  0 -1 -1]
            [ 0  1 -1  0]
            [-1 -1  3 -1]
            [-1  0 -1  2]
            sage: M = G.laplacian_matrix(boundary_first=True, sparse=False); M
            [ 2  0 -1 -1]
            [ 0  1 -1  0]
            [-1 -1  3 -1]
            [-1  0 -1  2]
            sage: M = G.laplacian_matrix(normalized=True); M
            [                   1 -1/6*sqrt(3)*sqrt(2) -1/6*sqrt(3)*sqrt(2)         -1/3*sqrt(3)]
            [-1/6*sqrt(3)*sqrt(2)                    1                 -1/2                    0]
            [-1/6*sqrt(3)*sqrt(2)                 -1/2                    1                    0]
            [        -1/3*sqrt(3)                    0                    0                    1]

            sage: Graph({0:[],1:[2]}).laplacian_matrix(normalized=True)
            [ 0  0  0]
            [ 0  1 -1]
            [ 0 -1  1]

        A weighted directed graph with loops, changing the variable ``indegree`` ::

            sage: G = DiGraph({1:{1:2,2:3}, 2:{1:4}}, weighted=True,sparse=True)
            sage: G.laplacian_matrix()
            [ 4 -3]
            [-4  3]

        ::

            sage: G = DiGraph({1:{1:2,2:3}, 2:{1:4}}, weighted=True,sparse=True)
            sage: G.laplacian_matrix(indegree=False)
            [ 3 -3]
            [-4  4]
        """
        from sage.matrix.constructor import diagonal_matrix
        from sage.functions.all import sqrt

        if weighted is None:
            weighted = self._weighted

        if weighted:
            M = self.weighted_adjacency_matrix(**kwds)
        else:
            M = self.adjacency_matrix(**kwds)

        D = M.parent(0)

        if M.is_sparse():
            row_sums = {}
            if indegree:
                for (i,j), entry in M.dict().iteritems():
                    row_sums[j] = row_sums.get(j, 0) + entry
            else:
                for (i,j), entry in M.dict().iteritems():
                    row_sums[i] = row_sums.get(i, 0) + entry


            for i in range(M.nrows()):
                D[i,i] += row_sums.get(i, 0)

        else:
            if indegree:
                col_sums=[sum(v) for v in M.columns()]
                for i in range(M.nrows()):
                    D[i,i] += col_sums[i]
            else:
                row_sums=[sum(v) for v in M.rows()]
                for i in range(M.nrows()):
                    D[i,i] += row_sums[i]

        if normalized:
            Dsqrt = diagonal_matrix([1/sqrt(D[i,i]) if D[i,i]>0 else 1 \
                                     for i in range(D.nrows())])
            return Dsqrt*(D-M)*Dsqrt
        else:
            return D-M

    laplacian_matrix = kirchhoff_matrix

    ### Attributes

    def get_boundary(self):
        """
        Returns the boundary of the (di)graph.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.set_boundary([0,1,2,3,4])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: G.get_boundary()
            [0, 1, 2, 3, 4]
        """
        from sage.misc.superseded import deprecation
        deprecation(15494, "The boundary parameter is deprecated and will soon disappear.")

        return self._boundary

    def set_boundary(self, boundary):
        """
        Sets the boundary of the (di)graph.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.set_boundary([0,1,2,3,4])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: G.get_boundary()
            [0, 1, 2, 3, 4]
            sage: G.set_boundary((1..4))
            sage: G.get_boundary()
            [1, 2, 3, 4]
        """
        from sage.misc.superseded import deprecation
        deprecation(15494, "The boundary parameter is deprecated and will soon disappear.")

        if isinstance(boundary,list):
            self._boundary = boundary
        else:
            self._boundary = list(boundary)

    def set_embedding(self, embedding):
        """
        Sets a combinatorial embedding dictionary to ``_embedding`` attribute.

        Dictionary is organized with vertex labels as keys and a list of
        each vertex's neighbors in clockwise order.

        Dictionary is error-checked for validity.

        INPUT:

        -  ``embedding`` - a dictionary

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.set_embedding({0: [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]})
            sage: G.set_embedding({'s': [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]})
            Traceback (most recent call last):
            ...
            ValueError: embedding is not valid for Petersen graph
        """
        if self._check_embedding_validity(embedding):
            self._embedding = embedding
        else:
            raise ValueError('embedding is not valid for %s'%self)

    def get_embedding(self):
        """
        Returns the attribute _embedding if it exists.

        ``_embedding`` is a dictionary organized with vertex labels as keys and a
        list of each vertex's neighbors in clockwise order.

        Error-checked to insure valid embedding is returned.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.genus()
            1
            sage: G.get_embedding()
            {0: [1, 4, 5], 1: [0, 2, 6], 2: [1, 3, 7], 3: [2, 4, 8], 4: [0, 3, 9], 5: [0, 7, 8], 6: [1, 9, 8], 7: [2, 5, 9], 8: [3, 6, 5], 9: [4, 6, 7]}
        """
        if self._check_embedding_validity():
            return self._embedding
        else:
            raise ValueError('%s has been modified and the embedding is no longer valid'%self)

    def _check_embedding_validity(self, embedding=None):
        """
        Checks whether an _embedding attribute is well defined.

        If the ``_embedding`` attribute exists, it is checked for
        accuracy. Returns True if everything is okay, False otherwise.

        If embedding=None will test the attribute _embedding.

        EXAMPLES::

            sage: d = {0: [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]}
            sage: G = graphs.PetersenGraph()
            sage: G._check_embedding_validity(d)
            True

        TESTS::

            sage: G.check_embedding_validity(d)
            doctest:...: DeprecationWarning: check_embedding_validity is deprecated. Please use _check_embedding_validity instead.
            See http://trac.sagemath.org/15551 for details.
            True

        """
        if embedding is None:
            embedding = getattr(self, '_embedding', None)
        if embedding is None:
            return False
        if len(embedding) != self.order():
            return False
        if self._directed:
            connected = lambda u,v : self.has_edge(u,v) or self.has_edge(v,u)
        else:
            connected = lambda u,v : self.has_edge(u,v)
        for v in embedding:
            if not self.has_vertex(v):
                return False
            if len(embedding[v]) != len(self.neighbors(v)):
                return False
            for u in embedding[v]:
                if not connected(v,u):
                    return False
        return True

    check_embedding_validity = deprecated_function_alias(15551, _check_embedding_validity)

    def has_loops(self):
        """
        Returns whether there are loops in the (di)graph.

        EXAMPLES::

            sage: G = Graph(loops=True); G
            Looped graph on 0 vertices
            sage: G.has_loops()
            False
            sage: G.allows_loops()
            True
            sage: G.add_edge((0,0))
            sage: G.has_loops()
            True
            sage: G.loops()
            [(0, 0, None)]
            sage: G.allow_loops(False); G
            Graph on 1 vertex
            sage: G.has_loops()
            False
            sage: G.edges()
            []

            sage: D = DiGraph(loops=True); D
            Looped digraph on 0 vertices
            sage: D.has_loops()
            False
            sage: D.allows_loops()
            True
            sage: D.add_edge((0,0))
            sage: D.has_loops()
            True
            sage: D.loops()
            [(0, 0, None)]
            sage: D.allow_loops(False); D
            Digraph on 1 vertex
            sage: D.has_loops()
            False
            sage: D.edges()
            []
        """
        if self.allows_loops():
            for v in self:
                if self.has_edge(v,v):
                    return True
        return False

    def allows_loops(self):
        """
        Returns whether loops are permitted in the (di)graph.

        EXAMPLES::

            sage: G = Graph(loops=True); G
            Looped graph on 0 vertices
            sage: G.has_loops()
            False
            sage: G.allows_loops()
            True
            sage: G.add_edge((0,0))
            sage: G.has_loops()
            True
            sage: G.loops()
            [(0, 0, None)]
            sage: G.allow_loops(False); G
            Graph on 1 vertex
            sage: G.has_loops()
            False
            sage: G.edges()
            []

            sage: D = DiGraph(loops=True); D
            Looped digraph on 0 vertices
            sage: D.has_loops()
            False
            sage: D.allows_loops()
            True
            sage: D.add_edge((0,0))
            sage: D.has_loops()
            True
            sage: D.loops()
            [(0, 0, None)]
            sage: D.allow_loops(False); D
            Digraph on 1 vertex
            sage: D.has_loops()
            False
            sage: D.edges()
            []
        """
        return self._backend.loops(None)

    def allow_loops(self, new, check=True):
        """
        Changes whether loops are permitted in the (di)graph.

        INPUT:

        - ``new`` - boolean.

        EXAMPLES::

            sage: G = Graph(loops=True); G
            Looped graph on 0 vertices
            sage: G.has_loops()
            False
            sage: G.allows_loops()
            True
            sage: G.add_edge((0,0))
            sage: G.has_loops()
            True
            sage: G.loops()
            [(0, 0, None)]
            sage: G.allow_loops(False); G
            Graph on 1 vertex
            sage: G.has_loops()
            False
            sage: G.edges()
            []

            sage: D = DiGraph(loops=True); D
            Looped digraph on 0 vertices
            sage: D.has_loops()
            False
            sage: D.allows_loops()
            True
            sage: D.add_edge((0,0))
            sage: D.has_loops()
            True
            sage: D.loops()
            [(0, 0, None)]
            sage: D.allow_loops(False); D
            Digraph on 1 vertex
            sage: D.has_loops()
            False
            sage: D.edges()
            []
        """
        if new is False and check:
            self.remove_loops()
        self._backend.loops(new)

    def loops(self, labels=True):
        """
        Returns any loops in the (di)graph.

        INPUT:

        - ``new`` -- deprecated

        - ``labels`` -- whether returned edges have labels ((u,v,l)) or not ((u,v)).

        EXAMPLES::

            sage: G = Graph(loops=True); G
            Looped graph on 0 vertices
            sage: G.has_loops()
            False
            sage: G.allows_loops()
            True
            sage: G.add_edge((0,0))
            sage: G.has_loops()
            True
            sage: G.loops()
            [(0, 0, None)]
            sage: G.allow_loops(False); G
            Graph on 1 vertex
            sage: G.has_loops()
            False
            sage: G.edges()
            []

            sage: D = DiGraph(loops=True); D
            Looped digraph on 0 vertices
            sage: D.has_loops()
            False
            sage: D.allows_loops()
            True
            sage: D.add_edge((0,0))
            sage: D.has_loops()
            True
            sage: D.loops()
            [(0, 0, None)]
            sage: D.allow_loops(False); D
            Digraph on 1 vertex
            sage: D.has_loops()
            False
            sage: D.edges()
            []

            sage: G = graphs.PetersenGraph()
            sage: G.loops()
            []

        """
        loops = []
        for v in self:
            loops += self.edge_boundary([v], [v], labels)
        return loops

    def has_multiple_edges(self, to_undirected=False):
        """
        Returns whether there are multiple edges in the (di)graph.

        INPUT:

        - ``to_undirected`` -- (default: False) If True, runs the test on the undirected version of a DiGraph.
          Otherwise, treats DiGraph edges (u,v) and (v,u) as unique individual edges.

        EXAMPLES::

            sage: G = Graph(multiedges=True,sparse=True); G
            Multi-graph on 0 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.allows_multiple_edges()
            True
            sage: G.add_edges([(0,1)]*3)
            sage: G.has_multiple_edges()
            True
            sage: G.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: G.allow_multiple_edges(False); G
            Graph on 2 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.edges()
            [(0, 1, None)]

            sage: D = DiGraph(multiedges=True,sparse=True); D
            Multi-digraph on 0 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.allows_multiple_edges()
            True
            sage: D.add_edges([(0,1)]*3)
            sage: D.has_multiple_edges()
            True
            sage: D.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: D.allow_multiple_edges(False); D
            Digraph on 2 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.edges()
            [(0, 1, None)]

            sage: G = DiGraph({1:{2: 'h'}, 2:{1:'g'}},sparse=True)
            sage: G.has_multiple_edges()
            False
            sage: G.has_multiple_edges(to_undirected=True)
            True
            sage: G.multiple_edges()
            []
            sage: G.multiple_edges(to_undirected=True)
            [(1, 2, 'h'), (2, 1, 'g')]
        """
        if self.allows_multiple_edges() or (self._directed and to_undirected):
            if self._directed:
                for u in self:
                    s = set()
                    for a,b,c in self.outgoing_edge_iterator(u):
                        if b in s:
                            return True
                        s.add(b)
                    if to_undirected:
                        for a,b,c in self.incoming_edge_iterator(u):
                            if a in s:
                                return True
                            s.add(a)
            else:
                for u in self:
                    s = set()
                    for a,b,c in self.edge_iterator(u):
                        if a is u:
                            if b in s:
                                return True
                            s.add(b)
                        if b is u:
                            if a in s:
                                return True
                            s.add(a)
        return False

    def allows_multiple_edges(self):
        """
        Returns whether multiple edges are permitted in the (di)graph.

        EXAMPLES::

            sage: G = Graph(multiedges=True,sparse=True); G
            Multi-graph on 0 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.allows_multiple_edges()
            True
            sage: G.add_edges([(0,1)]*3)
            sage: G.has_multiple_edges()
            True
            sage: G.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: G.allow_multiple_edges(False); G
            Graph on 2 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.edges()
            [(0, 1, None)]

            sage: D = DiGraph(multiedges=True,sparse=True); D
            Multi-digraph on 0 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.allows_multiple_edges()
            True
            sage: D.add_edges([(0,1)]*3)
            sage: D.has_multiple_edges()
            True
            sage: D.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: D.allow_multiple_edges(False); D
            Digraph on 2 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.edges()
            [(0, 1, None)]
        """
        return self._backend.multiple_edges(None)

    def allow_multiple_edges(self, new, check=True):
        """
        Changes whether multiple edges are permitted in the (di)graph.

        INPUT:

        - ``new`` - boolean.

        EXAMPLES::

            sage: G = Graph(multiedges=True,sparse=True); G
            Multi-graph on 0 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.allows_multiple_edges()
            True
            sage: G.add_edges([(0,1)]*3)
            sage: G.has_multiple_edges()
            True
            sage: G.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: G.allow_multiple_edges(False); G
            Graph on 2 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.edges()
            [(0, 1, None)]

            sage: D = DiGraph(multiedges=True,sparse=True); D
            Multi-digraph on 0 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.allows_multiple_edges()
            True
            sage: D.add_edges([(0,1)]*3)
            sage: D.has_multiple_edges()
            True
            sage: D.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: D.allow_multiple_edges(False); D
            Digraph on 2 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.edges()
            [(0, 1, None)]
        """
        seen = set()

        # TODO: this should be much faster for c_graphs, but for now we just do this
        if self.allows_multiple_edges() and new is False and check:
            for u,v,l in self.multiple_edges():
                if (u,v) in seen:
                    self.delete_edge(u,v,l)
                else:
                    seen.add((u,v))

        self._backend.multiple_edges(new)

    def multiple_edges(self, to_undirected=False, labels=True):
        """
        Returns any multiple edges in the (di)graph.

        EXAMPLES::

            sage: G = Graph(multiedges=True,sparse=True); G
            Multi-graph on 0 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.allows_multiple_edges()
            True
            sage: G.add_edges([(0,1)]*3)
            sage: G.has_multiple_edges()
            True
            sage: G.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: G.allow_multiple_edges(False); G
            Graph on 2 vertices
            sage: G.has_multiple_edges()
            False
            sage: G.edges()
            [(0, 1, None)]

            sage: D = DiGraph(multiedges=True,sparse=True); D
            Multi-digraph on 0 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.allows_multiple_edges()
            True
            sage: D.add_edges([(0,1)]*3)
            sage: D.has_multiple_edges()
            True
            sage: D.multiple_edges()
            [(0, 1, None), (0, 1, None), (0, 1, None)]
            sage: D.allow_multiple_edges(False); D
            Digraph on 2 vertices
            sage: D.has_multiple_edges()
            False
            sage: D.edges()
            [(0, 1, None)]

            sage: G = DiGraph({1:{2: 'h'}, 2:{1:'g'}},sparse=True)
            sage: G.has_multiple_edges()
            False
            sage: G.has_multiple_edges(to_undirected=True)
            True
            sage: G.multiple_edges()
            []
            sage: G.multiple_edges(to_undirected=True)
            [(1, 2, 'h'), (2, 1, 'g')]
        """
        from sage.misc.superseded import deprecation
        multi_edges = []
        if self._directed and not to_undirected:
            for v in self:
                for u in self.neighbor_in_iterator(v):
                    edges = self.edge_boundary([u], [v], labels)
                    if len(edges) > 1:
                        multi_edges += edges
        else:
            to_undirected *= self._directed
            for v in self:
                for u in self.neighbor_iterator(v):
                    if hash(u) >= hash(v):
                        edges = self.edge_boundary([v], [u], labels)
                        if to_undirected:
                            edges += self.edge_boundary([u],[v], labels)
                        if len(edges) > 1:
                            multi_edges += edges
        return multi_edges

    def name(self, new=None):
        """
        Returns or sets the graph's name.

        INPUT:

        - ``new`` - if not None, then this becomes the new name of the (di)graph.
          (if new == '', removes any name)

        EXAMPLES::

            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G
            Graph on 10 vertices
            sage: G.name("Petersen Graph"); G
            Petersen Graph: Graph on 10 vertices
            sage: G.name(new=""); G
            Graph on 10 vertices
            sage: G.name()
            ''

        Name of an immutable graph :trac:`15681` ::

            sage: g = graphs.PetersenGraph()
            sage: gi = g.copy(immutable=True)
            sage: gi.name()
            'Petersen graph'
            sage: gi.name("Hey")
            Traceback (most recent call last):
            ...
            NotImplementedError: An immutable graph does not change name
        """
        if new is None:
            return getattr(self, '_name', "")

        if getattr(self, '_immutable', False):
            raise NotImplementedError("An immutable graph does not change name")

        self._name = str(new)

    def get_pos(self, dim = 2):
        """
        Returns the position dictionary, a dictionary specifying the
        coordinates of each vertex.

        EXAMPLES: By default, the position of a graph is None::

            sage: G = Graph()
            sage: G.get_pos()
            sage: G.get_pos() is None
            True
            sage: P = G.plot(save_pos=True)
            sage: G.get_pos()
            {}

        Some of the named graphs come with a pre-specified positioning::

            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: (...e-17, 1.0),
             ...
             9: (0.475..., 0.154...)}
        """
        if dim == 2:
            return self._pos
        elif dim == 3:
            return getattr(self, "_pos3d", None)
        else:
            raise ValueError("dim must be 2 or 3")

    def _check_pos_validity(self, pos=None, dim = 2):
        r"""
        Checks whether pos specifies two (resp. 3) coordinates for every vertex (and no more vertices).

        INPUT:

            - ``pos`` - a position dictionary for a set of vertices
            - ``dim`` - 2 or 3 (default: 3

        OUTPUT:

        If ``pos`` is ``None`` then the position dictionary of ``self`` is
        investigated, otherwise the position dictionary provided in  ``pos`` is
        investigated.  The function returns ``True`` if the dictionary is of the
        correct form for ``self``.

        EXAMPLES::

            sage: p = {0: [1, 5], 1: [0, 2], 2: [1, 3], 3: [8, 2], 4: [0, 9], 5: [0, 8], 6: [8, 1], 7: [9, 5], 8: [3, 5], 9: [6, 7]}
            sage: G = graphs.PetersenGraph()
            sage: G._check_pos_validity(p)
            True

        TESTS::

            sage: G.check_pos_validity(p)
            doctest:...: DeprecationWarning: check_pos_validity is deprecated. Please use _check_pos_validity instead.
            See http://trac.sagemath.org/15551 for details.
            True
        """
        if pos is None:
            pos = self.get_pos(dim = dim)
        if pos is None:
            return False
        if len(pos) != self.order():
            return False
        for v in pos:
            if not self.has_vertex(v):
                return False
            if len(pos[v]) != dim:
                return False
        return True

    check_pos_validity = deprecated_function_alias(15551,_check_pos_validity)

    def set_pos(self, pos, dim = 2):
        """
        Sets the position dictionary, a dictionary specifying the
        coordinates of each vertex.

        EXAMPLES: Note that set_pos will allow you to do ridiculous things,
        which will not blow up until plotting::

            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: (..., ...),
             ...
             9: (..., ...)}

        ::

            sage: G.set_pos('spam')
            sage: P = G.plot()
            Traceback (most recent call last):
            ...
            TypeError: string indices must be integers, not str
        """
        if dim == 2:
            self._pos = pos
        elif dim == 3:
            self._pos3d = pos
        else:
            raise ValueError("dim must be 2 or 3")

    def weighted(self, new=None):
        """
        Whether the (di)graph is to be considered as a weighted (di)graph.

        INPUT:

        - ``new`` (optional bool): If it is provided, then the weightedness
          flag is set accordingly. This is not allowed for immutable graphs.

        .. NOTE::

            Changing the weightedness flag changes the ``==``-class of
            a graph and is thus not allowed for immutable graphs.

            Edge weightings can still exist for (di)graphs ``G`` where
            ``G.weighted()`` is ``False``.

        EXAMPLES:

        Here we have two graphs with different labels, but ``weighted()`` is
        ``False`` for both, so we just check for the presence of edges::

            sage: G = Graph({0:{1:'a'}}, sparse=True)
            sage: H = Graph({0:{1:'b'}}, sparse=True)
            sage: G == H
            True

        Now one is weighted and the other is not, and thus the graphs are
        not equal::

            sage: G.weighted(True)
            sage: H.weighted()
            False
            sage: G == H
            False

        However, if both are weighted, then we finally compare 'a' to 'b'::

            sage: H.weighted(True)
            sage: G == H
            False

        TESTS:

        Ensure that :trac:`10490` is fixed: allows a weighted graph to be
        set as unweighted. ::

            sage: G = Graph({1:{2:3}})
            sage: G.weighted()
            False
            sage: G.weighted('a')
            sage: G.weighted(True)
            sage: G.weighted()
            True
            sage: G.weighted('a')
            sage: G.weighted()
            True
            sage: G.weighted(False)
            sage: G.weighted()
            False
            sage: G.weighted('a')
            sage: G.weighted()
            False
            sage: G.weighted(True)
            sage: G.weighted()
            True

        Ensure that graphs using the static sparse backend can not be mutated
        using this method, as fixed in :trac:`15278`::

            sage: G = graphs.PetersenGraph()
            sage: G.weighted()
            False
            sage: H = copy(G)
            sage: H == G
            True
            sage: H.weighted(True)
            sage: H == G
            False
            sage: G_imm = Graph(G, data_structure="static_sparse")
            sage: G_imm == G
            True
            sage: G_imm.weighted()
            False
            sage: G_imm.weighted(True)
            Traceback (most recent call last):
            ...
            TypeError: This graph is immutable and can thus not be changed.
            Create a mutable copy, e.g., by `g.copy(immutable=False)`
            sage: G_mut = G_imm.copy(immutable=False)
            sage: G_mut == G_imm
            True
            sage: G_mut.weighted(True)
            sage: G_mut == G_imm
            False
            sage: G_mut == H
            True

        """
        if new is not None:
            if getattr(self, '_immutable', False):
                raise TypeError("This graph is immutable and can thus not be changed. "
                                "Create a mutable copy, e.g., by `g.copy(immutable=False)`")
            if new in [True, False]:
                self._weighted = new
        else:
            return self._weighted

    ### Properties

    def antisymmetric(self):
        r"""
        Tests whether the graph is antisymmetric.

        A graph represents an antisymmetric relation if there being a path
        from a vertex x to a vertex y implies that there is not a path from
        y to x unless x=y.

        A directed acyclic graph is antisymmetric. An undirected graph is
        never antisymmetric unless it is just a union of isolated
        vertices.

        ::

            sage: graphs.RandomGNP(20,0.5).antisymmetric()
            False
            sage: digraphs.RandomDirectedGNR(20,0.5).antisymmetric()
            True
        """
        if not self._directed:
            if self.size()-len(self.loop_edges())>0:
                return False
            else:
                return True
        from copy import copy
        g = copy(self)
        g.allow_multiple_edges(False)
        g.allow_loops(False)
        g = g.transitive_closure()
        gpaths = g.edges(labels=False)
        for e in gpaths:
            if (e[1],e[0]) in gpaths:
                return False
        return True

    def density(self):
        """
        Returns the density (number of edges divided by number of possible
        edges).

        In the case of a multigraph, raises an error, since there is an
        infinite number of possible edges.

        EXAMPLES::

            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G.density()
            1/3
            sage: G = Graph({0:[1,2], 1:[0] }); G.density()
            2/3
            sage: G = DiGraph({0:[1,2], 1:[0] }); G.density()
            1/2

        Note that there are more possible edges on a looped graph::

            sage: G.allow_loops(True)
            sage: G.density()
            1/3
        """
        if self.has_multiple_edges():
            raise TypeError("Density is not well-defined for multigraphs.")
        n = self.order()
        if self.allows_loops():
            if n == 0:
                return Rational(0)
            if self._directed:
                return Rational(self.size())/Rational(n**2)
            else:
                return Rational(self.size())/Rational((n**2 + n)/2)
        else:
            if n < 2:
                return Rational(0)
            if self._directed:
                return Rational(self.size())/Rational((n**2 - n))
            else:
                return Rational(self.size())/Rational((n**2 - n)/2)

    def is_eulerian(self, path=False):
        r"""
        Return true if the graph has a (closed) tour that visits each edge exactly
        once.

        INPUT:

        - ``path`` -- by default this function finds if the graph contains a closed
          tour visiting each edge once, i.e. an eulerian cycle. If you want to test
          the existence of an eulerian path, set this argument to ``True``. Graphs
          with this property are sometimes called semi-eulerian.

        OUTPUT:

        ``True`` or ``False`` for the closed tour case. For an open tour search
        (``path``=``True``) the function returns ``False`` if the graph is not
        semi-eulerian, or a tuple (u, v) in the other case. This tuple defines the
        edge that would make the graph eulerian, i.e. close an existing open tour.
        This edge may or may not be already present in the graph.

        EXAMPLES::

            sage: graphs.CompleteGraph(4).is_eulerian()
            False
            sage: graphs.CycleGraph(4).is_eulerian()
            True
            sage: g = DiGraph({0:[1,2], 1:[2]}); g.is_eulerian()
            False
            sage: g = DiGraph({0:[2], 1:[3], 2:[0,1], 3:[2]}); g.is_eulerian()
            True
            sage: g = DiGraph({0:[1], 1:[2], 2:[0], 3:[]}); g.is_eulerian()
            True
            sage: g = Graph([(1,2), (2,3), (3,1), (4,5), (5,6), (6,4)]); g.is_eulerian()
            False

        ::

            sage: g = DiGraph({0: [1]}); g.is_eulerian(path=True)
            (1, 0)
            sage: graphs.CycleGraph(4).is_eulerian(path=True)
            False
            sage: g = DiGraph({0: [1], 1: [2,3], 2: [4]}); g.is_eulerian(path=True)
            False

        ::

            sage: g = Graph({0:[1,2,3], 1:[2,3], 2:[3,4], 3:[4]}, multiedges=True)
            sage: g.is_eulerian()
            False
            sage: e = g.is_eulerian(path=True); e
            (0, 1)
            sage: g.add_edge(e)
            sage: g.is_eulerian(path=False)
            True
            sage: g.is_eulerian(path=True)
            False

        TESTS::

            sage: g = Graph({0:[], 1:[], 2:[], 3:[]}); g.is_eulerian()
            True
        """

        # unconnected graph can still be eulerian if all components
        # up to one doesn't contain any edge
        nontrivial_components = 0
        for cc in self.connected_components():
            if len(cc) > 1:
                nontrivial_components += 1
            if nontrivial_components > 1:
                return False

        uv = [None, None]
        if self._directed:
            for v in self.vertex_iterator():
                # loops don't matter since they count in both the in and out degree.
                if self.in_degree(v) != self.out_degree(v):
                    if path:
                        diff = self.out_degree(v) - self.in_degree(v)
                        if abs(diff) > 1:
                            return False
                        else:
                            # if there was another vertex with the same sign of difference...
                            if uv[(diff+1)//2] is not None:
                                return False # ... the graph is not semi-eulerian
                            else:
                                uv[(diff+1)//2] = v
                    else:
                        return False
        else:
            for v in self.vertex_iterator():
                # loops don't matter since they add an even number to the degree
                if self.degree(v) % 2 != 0:
                    if not path:
                        return False
                    else:
                        if uv[0] is None or uv[1] is None:
                            uv[0 if uv[0] is None else 1] = v
                        else:
                            return False

        if path and (uv[0] is None or uv[1] is None):
            return False

        return True if not path else tuple(uv)

    def order(self):
        """
        Returns the number of vertices. Note that len(G) returns the number
        of vertices in G also.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.order()
            10

        ::

            sage: G = graphs.TetrahedralGraph()
            sage: len(G)
            4
        """
        return self._backend.num_verts()

    __len__ = order

    num_verts = order

    def size(self):
        """
        Returns the number of edges.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.size()
            15
        """
        return self._backend.num_edges(self._directed)

    num_edges = size

    ### Orientations

    def eulerian_orientation(self):
        r"""
        Returns a DiGraph which is an Eulerian orientation of the current graph.

        An Eulerian graph being a graph such that any vertex has an even degree,
        an Eulerian orientation of a graph is an orientation of its edges such
        that each vertex `v` verifies `d^+(v)=d^-(v)=d(v)/2`, where `d^+` and
        `d^-` respectively represent the out-degree and the in-degree of a vertex.

        If the graph is not Eulerian, the orientation verifies for any vertex `v`
        that `| d^+(v)-d^-(v) | \leq 1`.

        ALGORITHM:

        This algorithm is a random walk through the edges of the graph, which
        orients the edges according to the walk. When a vertex is reached which
        has no non-oriented edge ( this vertex must have odd degree ), the
        walk resumes at another vertex of odd degree, if any.

        This algorithm has complexity `O(m)`, where `m` is the number of edges
        in the graph.

        EXAMPLES:

        The CubeGraph with parameter 4, which is regular of even degree, has an
        Eulerian orientation such that `d^+=d^-`::

            sage: g=graphs.CubeGraph(4)
            sage: g.degree()
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
            sage: o=g.eulerian_orientation()
            sage: o.in_degree()
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            sage: o.out_degree()
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

        Secondly, the Petersen Graph, which is 3 regular has an orientation
        such that the difference between `d^+` and `d^-` is at most 1::

            sage: g=graphs.PetersenGraph()
            sage: o=g.eulerian_orientation()
            sage: o.in_degree()
            [2, 2, 2, 2, 2, 1, 1, 1, 1, 1]
            sage: o.out_degree()
            [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
        """
        from copy import copy
        g=copy(self)
        from sage.graphs.digraph import DiGraph
        d=DiGraph()
        d.add_vertices(g.vertex_iterator())


        # list of vertices of odd degree
        from itertools import izip
        odd=[x for (x,deg) in izip(g.vertex_iterator(),g.degree_iterator()) if deg%2==1]

        # Picks the first vertex, which is preferably an odd one
        if len(odd)>0:
            v=odd.pop()
        else:
            v=g.edge_iterator(labels=None).next()[0]
            odd.append(v)
        # Stops when there is no edge left
        while True:

            # If there is an edge adjacent to the current one
            if g.degree(v)>0:
                e = g.edge_iterator(v).next()
                g.delete_edge(e)
                if e[0]!=v:
                    e=(e[1],e[0],e[2])
                d.add_edge(e)
                v=e[1]

            # The current vertex is isolated
            else:
                odd.remove(v)

                # jumps to another odd vertex if possible
                if len(odd)>0:
                    v=odd.pop()
                # Else jumps to an ever vertex which is not isolated
                elif g.size()>0:
                    v=g.edge_iterator().next()[0]
                    odd.append(v)
                # If there is none, we are done !
                else:
                    return d

    def eulerian_circuit(self, return_vertices=False, labels=True, path=False):
        r"""
        Return a list of edges forming an eulerian circuit if one exists.
        Otherwise return False.

        This is implemented using Hierholzer's algorithm.

        INPUT:

        -  ``return_vertices`` -- (default: ``False``) optionally provide a list of
           vertices for the path

        -  ``labels`` -- (default: ``True``) whether to return edges with labels
           (3-tuples)

        -  ``path`` -- (default: ``False``) find an eulerian path instead

        OUTPUT:

        either ([edges], [vertices]) or [edges] of an Eulerian circuit (or path)

        EXAMPLES::

            sage: g=graphs.CycleGraph(5);
            sage: g.eulerian_circuit()
            [(0, 4, None), (4, 3, None), (3, 2, None), (2, 1, None), (1, 0, None)]
            sage: g.eulerian_circuit(labels=False)
            [(0, 4), (4, 3), (3, 2), (2, 1), (1, 0)]

        ::

            sage: g = graphs.CompleteGraph(7)
            sage: edges, vertices = g.eulerian_circuit(return_vertices=True)
            sage: vertices
            [0, 6, 5, 4, 6, 3, 5, 2, 4, 3, 2, 6, 1, 5, 0, 4, 1, 3, 0, 2, 1, 0]

        ::

            sage: graphs.CompleteGraph(4).eulerian_circuit()
            False

        A disconnected graph can be eulerian::

            sage: g = Graph({0: [], 1: [2], 2: [3], 3: [1], 4: []})
            sage: g.eulerian_circuit(labels=False)
            [(1, 3), (3, 2), (2, 1)]

        ::

            sage: g = DiGraph({0: [1], 1: [2, 4], 2:[3], 3:[1]})
            sage: g.eulerian_circuit(labels=False, path=True)
            [(0, 1), (1, 2), (2, 3), (3, 1), (1, 4)]

        ::

            sage: g = Graph({0:[1,2,3], 1:[2,3], 2:[3,4], 3:[4]})
            sage: g.is_eulerian(path=True)
            (0, 1)
            sage: g.eulerian_circuit(labels=False, path=True)
            [(1, 3), (3, 4), (4, 2), (2, 3), (3, 0), (0, 2), (2, 1), (1, 0)]

        TESTS::

            sage: Graph({'H': ['G','L','L','D'], 'L': ['G','D']}).eulerian_circuit(labels=False)
            [('H', 'D'), ('D', 'L'), ('L', 'G'), ('G', 'H'), ('H', 'L'), ('L', 'H')]
            sage: Graph({0: [0, 1, 1, 1, 1]}).eulerian_circuit(labels=False)
            [(0, 1), (1, 0), (0, 1), (1, 0), (0, 0)]
        """
        # trivial case
        if self.order() == 0:
            return ([], []) if return_vertices else []

        # check if the graph has proper properties to be eulerian
        edge = self.is_eulerian(path=path)
        if not edge:
            return False
        if path:
            start_vertex = edge[0]

        edges = []
        vertices = []

        # we'll remove edges as we go, so let's preserve the graph structure
        if self.is_directed():
            g = self.reverse()  # so the output will be in the proper order
        else:
            from copy import copy
            g = copy(self)

        if not path:
            # get the first vertex with degree>0
            start_vertex = None
            for v in g.vertex_iterator():
                if g.degree(v) != 0:
                    start_vertex = v
                    break

        # (where to return?, what was the way?)
        stack = [ (start_vertex, None) ]

        while len(stack) != 0:
            v, e = stack.pop()

            degr = g.out_degree(v) if self.is_directed() else g.degree(v)
            if degr == 0:
                vertices.append(v)
                if e is not None:
                    edges.append(e if labels else (e[0], e[1]))
            else:
                if self.is_directed():
                    next_edge = g.outgoing_edge_iterator(v).next()
                else:
                    next_edge = g.edge_iterator(v).next()

                if next_edge[0] == v:  # in the undirected case we want to
                                       # save the direction of traversal
                    next_edge_new = (next_edge[1], next_edge[0], next_edge[2])
                else:
                    next_edge_new = next_edge
                next_vertex = next_edge_new[0]

                stack.append((v, e))
                stack.append((next_vertex, next_edge_new))

                g.delete_edge(next_edge)

        if return_vertices:
            return edges, vertices
        else:
            return edges

    def min_spanning_tree(self,
                          weight_function=lambda e: 1,
                          algorithm="Kruskal",
                          starting_vertex=None,
                          check=False):
        r"""
        Returns the edges of a minimum spanning tree.

        INPUT:

        - ``weight_function`` -- A function that takes an edge and returns a
          numeric weight. Defaults to assigning each edge a weight of 1.

        - ``algorithm`` -- The algorithm to use in computing a minimum spanning
          tree of ``G``. The default is to use Kruskal's algorithm. The
          following algorithms are supported:

          - ``"Kruskal"`` -- Kruskal's algorithm.

          - ``"Prim_fringe"`` -- a variant of Prim's algorithm.
            ``"Prim_fringe"`` ignores the labels on the edges.

          - ``"Prim_edge"`` -- a variant of Prim's algorithm.

          - ``NetworkX`` -- Uses NetworkX's minimum spanning tree
            implementation.

        - ``starting_vertex`` -- The vertex from which to begin the search
          for a minimum spanning tree.

        - ``check`` -- Boolean; default: ``False``. Whether to first perform
          sanity checks on the input graph ``G``. If appropriate, ``check``
          is passed on to any minimum spanning tree functions that are
          invoked from the current method. See the documentation of the
          corresponding functions for details on what sort of sanity checks
          will be performed.

        OUTPUT:

        The edges of a minimum spanning tree of ``G``, if one exists, otherwise
        returns the empty list.

        .. seealso::

            - :func:`sage.graphs.spanning_tree.kruskal`

        EXAMPLES:

        Kruskal's algorithm::

            sage: g = graphs.CompleteGraph(5)
            sage: len(g.min_spanning_tree())
            4
            sage: weight = lambda e: 1 / ((e[0] + 1) * (e[1] + 1))
            sage: g.min_spanning_tree(weight_function=weight)
            [(3, 4, None), (2, 4, None), (1, 4, None), (0, 4, None)]
            sage: g = graphs.PetersenGraph()
            sage: g.allow_multiple_edges(True)
            sage: g.weighted(True)
            sage: g.add_edges(g.edges())
            sage: g.min_spanning_tree()
            [(0, 1, None), (0, 4, None), (0, 5, None), (1, 2, None), (1, 6, None), (2, 3, None), (2, 7, None), (3, 8, None), (4, 9, None)]

        Prim's algorithm::

            sage: g = graphs.CompleteGraph(5)
            sage: g.min_spanning_tree(algorithm='Prim_edge', starting_vertex=2, weight_function=weight)
            [(2, 4, None), (3, 4, None), (1, 4, None), (0, 4, None)]
            sage: g.min_spanning_tree(algorithm='Prim_fringe', starting_vertex=2, weight_function=weight)
            [(2, 4), (4, 3), (4, 1), (4, 0)]
        """
        if algorithm == "Kruskal":
            from spanning_tree import kruskal
            return kruskal(self, wfunction=weight_function, check=check)
        elif algorithm == "Prim_fringe":
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            tree = set([v])
            edges = []
            # Initialize fringe_list with v's neighbors. Fringe_list
            # contains fringe_vertex: (vertex_in_tree, weight) for each
            # fringe vertex.
            fringe_list = dict([u, (weight_function((v, u)), v)] for u in self[v])
            cmp_fun = lambda x: fringe_list[x]
            for i in range(self.order() - 1):
                # find the smallest-weight fringe vertex
                u = min(fringe_list, key=cmp_fun)
                edges.append((fringe_list[u][1], u))
                tree.add(u)
                fringe_list.pop(u)
                # update fringe list
                for neighbor in [v for v in self[u] if v not in tree]:
                    w = weight_function((u, neighbor))
                    if neighbor not in fringe_list or fringe_list[neighbor][0] > w:
                        fringe_list[neighbor] = (w, u)
            return edges
        elif algorithm == "Prim_edge":
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            sorted_edges = sorted(self.edges(), key=weight_function)
            tree = set([v])
            edges = []
            for _ in range(self.order() - 1):
                # Find a minimum-weight edge connecting a vertex in the tree
                # to something outside the tree. Remove the edges between
                # tree vertices for efficiency.
                i = 0
                while True:
                    e = sorted_edges[i]
                    v0, v1 = e[0], e[1]
                    if v0 in tree:
                        del sorted_edges[i]
                        if v1 in tree:
                            continue
                        edges.append(e)
                        tree.add(v1)
                        break
                    elif v1 in tree:
                        del sorted_edges[i]
                        edges.append(e)
                        tree.add(v0)
                        break
                    else:
                        i += 1
            return edges
        elif algorithm == "NetworkX":
            import networkx
            G = networkx.Graph([(u, v, dict(weight=weight_function((u, v)))) for u, v, l in self.edge_iterator()])
            return list(networkx.mst(G))
        else:
            raise NotImplementedError("Minimum Spanning Tree algorithm '%s' is not implemented." % algorithm)

    def spanning_trees_count(self, root_vertex=None):
        """
        Returns the number of spanning trees in a graph.

        In the case of a digraph, counts the number of spanning out-trees rooted
        in ``root_vertex``.  Default is to set first vertex as root.

        This computation uses Kirchhoff's Matrix Tree Theorem [1] to calculate
        the number of spanning trees. For complete graphs on `n` vertices the
        result can also be reached using Cayley's formula: the number of
        spanning trees are `n^(n-2)`.

        For digraphs, the augmented Kirchhoff Matrix as defined in [2] is
        used for calculations. Here the result is the number of out-trees
        rooted at a specific vertex.

        INPUT:

        - ``root_vertex`` -- integer (default: the first vertex) This is the vertex
          that will be used as root for all spanning out-trees if the graph
          is a directed graph.  This argument is ignored if the graph is not a digraph.

        .. SEEALSO::

            :meth:`~sage.graphs.graph.Graph.spanning_trees` -- enumerates all
            spanning trees of a graph.


        REFERENCES:

        - [1] http://mathworld.wolfram.com/MatrixTreeTheorem.html

        - [2] Lih-Hsing Hsu, Cheng-Kuan Lin, "Graph Theory and
          Interconnection Networks"

        AUTHORS:

        - Anders Jonsson (2009-10-10)

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.spanning_trees_count()
            2000

        ::

            sage: n = 11
            sage: G = graphs.CompleteGraph(n)
            sage: ST = G.spanning_trees_count()
            sage: ST == n^(n-2)
            True

        ::

            sage: M=matrix(3,3,[0,1,0,0,0,1,1,1,0])
            sage: D=DiGraph(M)
            sage: D.spanning_trees_count()
            1
            sage: D.spanning_trees_count(0)
            1
            sage: D.spanning_trees_count(2)
            2

        """

        if self.order() == 0:
            return 0

        if not self.is_directed():
            M = self.kirchhoff_matrix()
            M.subdivide(1,1)
            M2 = M.subdivision(1,1)
            return M2.determinant()
        else:
            if root_vertex is None:
                root_vertex=self.vertex_iterator().next()
            if root_vertex not in self.vertices():
                raise ValueError("Vertex (%s) not in the graph."%root_vertex)

            M=self.kirchhoff_matrix()

            index=self.vertices().index(root_vertex)
            M[index,index]+=1
            return abs(M.determinant())

    def cycle_basis(self, output='vertex'):
        r"""
        Returns a list of cycles which form a basis of the cycle space
        of ``self``.

        A basis of cycles of a graph is a minimal collection of cycles
        (considered as sets of edges) such that the edge set of any
        cycle in the graph can be written as a `Z/2Z` sum of the
        cycles in the basis.

        INPUT:

        - ``output`` (``'vertex'`` (default) or ``'edge'``) -- whether
          every cycle is given as a list of vertices or a list of
          edges.

        OUTPUT:

        A list of lists, each of them representing the vertices (or
        the edges) of a cycle in a basis.

        ALGORITHM:

        Uses the NetworkX library for graphs without multiple edges.

        Otherwise, by the standard algorithm using a spanning tree.

        EXAMPLE:

        A cycle basis in Petersen's Graph ::

            sage: g = graphs.PetersenGraph()
            sage: g.cycle_basis()
            [[1, 2, 7, 5, 0], [8, 3, 2, 7, 5], [4, 3, 2, 7, 5, 0], [4, 9, 7, 5, 0], [8, 6, 9, 7, 5], [1, 6, 9, 7, 5, 0]]

        One can also get the result as a list of lists of edges::

            sage: g.cycle_basis(output='edge')
            [[(1, 2, None), (2, 7, None), (7, 5, None), (5, 0, None),
            (0, 1, None)], [(8, 3, None), (3, 2, None), (2, 7, None),
            (7, 5, None), (5, 8, None)], [(4, 3, None), (3, 2, None),
            (2, 7, None), (7, 5, None), (5, 0, None), (0, 4, None)],
            [(4, 9, None), (9, 7, None), (7, 5, None), (5, 0, None),
            (0, 4, None)], [(8, 6, None), (6, 9, None), (9, 7, None),
            (7, 5, None), (5, 8, None)], [(1, 6, None), (6, 9, None),
            (9, 7, None), (7, 5, None), (5, 0, None), (0, 1, None)]]

        Checking the given cycles are algebraically free::

            sage: g = graphs.RandomGNP(30,.4)
            sage: basis = g.cycle_basis()

        Building the space of (directed) edges over `Z/2Z`. On the way,
        building a dictionary associating an unique vector to each
        undirected edge::

            sage: m = g.size()
            sage: edge_space = VectorSpace(FiniteField(2),m)
            sage: edge_vector = dict( zip( g.edges(labels = False), edge_space.basis() ) )
            sage: for (u,v),vec in edge_vector.items():
            ...      edge_vector[(v,u)] = vec

        Defining a lambda function associating a vector to the
        vertices of a cycle::

            sage: vertices_to_edges = lambda x : zip( x, x[1:] + [x[0]] )
            sage: cycle_to_vector = lambda x : sum( edge_vector[e] for e in vertices_to_edges(x) )

        Finally checking the cycles are a free set::

            sage: basis_as_vectors = map( cycle_to_vector, basis )
            sage: edge_space.span(basis_as_vectors).rank() == len(basis)
            True

        For undirected graphs with multiple edges::

            sage: G = Graph([(0,2,'a'),(0,2,'b'),(0,1,'c'),(1,2,'d')])
            sage: G.cycle_basis()
            [[0, 2], [2, 1, 0]]
            sage: G.cycle_basis(output='edge')
            [[(0, 2, 'a'), (2, 0, 'b')], [(2, 1, 'd'), (1, 0, 'c'),
            (0, 2, 'a')]]

        Disconnected graph::

            sage: G.add_cycle(["Hey", "Wuuhuu", "Really ?"])
            sage: G.cycle_basis()
            [[0, 2], [2, 1, 0], ['Really ?', 'Hey', 'Wuuhuu']]
            sage: G.cycle_basis(output='edge')
            [[(0, 2, 'a'), (2, 0, 'b')], [(2, 1, 'd'), (1, 0, 'c'), (0, 2, 'a')],
            [('Really ?', 'Hey', None), ('Hey', 'Wuuhuu', None), ('Wuuhuu', 'Really ?', None)]]

        Graph that allows multiple edges but does not contain any::

            sage: G = graphs.CycleGraph(3)
            sage: G.allow_multiple_edges(True)
            sage: G.cycle_basis()
            [[2, 1, 0]]

        Not yet implemented for directed graphs with multiple edges::

            sage: G = DiGraph([(0,2,'a'),(0,2,'b'),(0,1,'c'),(1,2,'d')])
            sage: G.cycle_basis()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented for directed graphs
            with multiple edges
        """
        if not output in ['vertex', 'edge']:
            raise ValueError('output must be either vertex or edge')

        if self.allows_multiple_edges():
            if self.is_directed():
                raise NotImplementedError('not implemented for directed '
                                          'graphs with multiple edges')

            if not self.is_connected():
                return sum([g.cycle_basis(output=output)
                            for g in self.connected_components_subgraphs()],
                           [])

            T = self.min_spanning_tree()
            return [self.subgraph(edges=T + [e]).is_forest(certificate=True,
                                                           output=output)[1]
                    for e in self.edges() if not e in T]

        # second case: there are no multiple edges
        import networkx
        cycle_basis_v = networkx.cycle_basis(self.networkx_graph(copy=False))
        if output == 'vertex':
            return cycle_basis_v

        def vertices_to_edges(x):
            return [(u[0], u[1], self.edge_label(u[0], u[1]))
                    for u in zip(x, x[1:] + [x[0]])]
        return map(vertices_to_edges, cycle_basis_v)


    ### Planarity

    def is_planar(self, on_embedding=None, kuratowski=False, set_embedding=False, set_pos=False):
        """
        Returns True if the graph is planar, and False otherwise. This
        wraps the reference implementation provided by John Boyer of the
        linear time planarity algorithm by edge addition due to Boyer
        Myrvold. (See reference code in graphs.planarity).

        Note - the argument on_embedding takes precedence over
        set_embedding. This means that only the on_embedding
        combinatorial embedding will be tested for planarity and no
        _embedding attribute will be set as a result of this function
        call, unless on_embedding is None.

        REFERENCE:

        - [1] John M. Boyer and Wendy J. Myrvold, On the Cutting Edge:
          Simplified O(n) Planarity by Edge Addition. Journal of Graph
          Algorithms and Applications, Vol. 8, No. 3, pp. 241-273,
          2004.

        INPUT:


        -  ``kuratowski`` - returns a tuple with boolean as
           first entry. If the graph is nonplanar, will return the Kuratowski
           subgraph or minor as the second tuple entry. If the graph is
           planar, returns None as the second entry.

        -  ``on_embedding`` - the embedding dictionary to test
           planarity on. (i.e.: will return True or False only for the given
           embedding.)

        -  ``set_embedding`` - whether or not to set the
           instance field variable that contains a combinatorial embedding
           (clockwise ordering of neighbors at each vertex). This value will
           only be set if a planar embedding is found. It is stored as a
           Python dict: v1: [n1,n2,n3] where v1 is a vertex and n1,n2,n3 are
           its neighbors.

        -  ``set_pos`` - whether or not to set the position
           dictionary (for plotting) to reflect the combinatorial embedding.
           Note that this value will default to False if set_emb is set to
           False. Also, the position dictionary will only be updated if a
           planar embedding is found.


        EXAMPLES::

            sage: g = graphs.CubeGraph(4)
            sage: g.is_planar()
            False

        ::

            sage: g = graphs.CircularLadderGraph(4)
            sage: g.is_planar(set_embedding=True)
            True
            sage: g.get_embedding()
            {0: [1, 4, 3],
             1: [2, 5, 0],
             2: [3, 6, 1],
             3: [0, 7, 2],
             4: [0, 5, 7],
             5: [1, 6, 4],
             6: [2, 7, 5],
             7: [4, 6, 3]}

        ::

            sage: g = graphs.PetersenGraph()
            sage: (g.is_planar(kuratowski=True))[1].adjacency_matrix()
            [0 1 0 0 0 1 0 0 0]
            [1 0 1 0 0 0 1 0 0]
            [0 1 0 1 0 0 0 1 0]
            [0 0 1 0 0 0 0 0 1]
            [0 0 0 0 0 0 1 1 0]
            [1 0 0 0 0 0 0 1 1]
            [0 1 0 0 1 0 0 0 1]
            [0 0 1 0 1 1 0 0 0]
            [0 0 0 1 0 1 1 0 0]

        ::

            sage: k43 = graphs.CompleteBipartiteGraph(4,3)
            sage: result = k43.is_planar(kuratowski=True); result
            (False, Graph on 6 vertices)
            sage: result[1].is_isomorphic(graphs.CompleteBipartiteGraph(3,3))
            True

        Multi-edged and looped graphs are partially supported::

            sage: G = Graph({0:[1,1]}, multiedges=True)
            sage: G.is_planar()
            True
            sage: G.is_planar(on_embedding={})
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute with embeddings of multiple-edged or looped graphs.
            sage: G.is_planar(set_pos=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute with embeddings of multiple-edged or looped graphs.
            sage: G.is_planar(set_embedding=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute with embeddings of multiple-edged or looped graphs.
            sage: G.is_planar(kuratowski=True)
            (True, None)

        ::

            sage: G = graphs.CompleteGraph(5)
            sage: G = Graph(G, multiedges=True)
            sage: G.add_edge(0,1)
            sage: G.is_planar()
            False
            sage: b,k = G.is_planar(kuratowski=True)
            sage: b
            False
            sage: k.vertices()
            [0, 1, 2, 3, 4]

        """
        if self.has_multiple_edges() or self.has_loops():
            if set_embedding or (on_embedding is not None) or set_pos:
                raise NotImplementedError("Cannot compute with embeddings of multiple-edged or looped graphs.")
            else:
                return self.to_simple().is_planar(kuratowski=kuratowski)
        if on_embedding:
            if self._check_embedding_validity(on_embedding):
                return (0 == self.genus(minimal=False,set_embedding=False,on_embedding=on_embedding))
            else:
                raise ValueError('on_embedding is not a valid embedding for %s'%self)
        else:
            from sage.graphs.planarity import is_planar
            G = self.to_undirected()
            planar = is_planar(G,kuratowski=kuratowski,set_pos=set_pos,set_embedding=set_embedding)
            if kuratowski:
                bool_result = planar[0]
            else:
                bool_result = planar
            if bool_result:
                if set_pos:
                    self._pos = G._pos
                if set_embedding:
                    self._embedding = G._embedding
            return planar

    def is_circular_planar(self, on_embedding=None, kuratowski=False,
                           set_embedding=True, boundary = None,
                           ordered=False, set_pos=False):
        """
        Tests whether the graph is circular planar (outerplanar)

        A graph is circular planar if it has a planar embedding in which all
        vertices can be drawn in order on a circle. This method can also be used
        to check the existence of a planar embedding in which the vertices of a
        specific set (the *boundary*) can be drawn on a circle, all other
        vertices being drawn inside of the circle. An order can be defined on
        the vertices of the boundary in order to define how they are to appear
        on the circle.

        INPUT:

        - ``kuratowski`` (boolean) - if set to True, returns a tuple with
           boolean first entry and the Kuratowski subgraph or minor as the
           second entry (see OUTPUT below). It is set to ``False`` by default.

        - ``on_embedding`` (boolean) - the embedding dictionary to test
           planarity on. (i.e.: will return ``True`` or ``False`` only for the
           given embedding). It is set to ``False`` by default.

        - ``set_embedding`` (boolean) - whether or not to set the instance field
           variable that contains a combinatorial embedding (clockwise ordering
           of neighbors at each vertex). This value will only be set if a
           circular planar embedding is found. It is stored as a Python dict:
           ``v1: [n1,n2,n3]`` where ``v1`` is a vertex and ``n1,n2,n3`` are its
           neighbors. It is set to ``True`` by default.

        - ``boundary`` - a set of vertices that are required to be drawn on the
          circle, all others being drawn inside of it. It is set to ``None`` by
          default, meaning that *all* vertices should be drawn on the boundary.

        - ``ordered`` (boolean) - whether or not to consider the order of the
           boundary. It is set to ``False`` by default, and required
           ``boundary`` to be defined.

        - ``set_pos`` - whether or not to set the position dictionary (for
           plotting) to reflect the combinatorial embedding.  Note that this
           value will default to False if set_emb is set to False. Also, the
           position dictionary will only be updated if a circular planar
           embedding is found.

        OUTPUT:

        The method returns ``True`` if the graph is circular planar, and
        ``False`` if it is not.

        If ``kuratowski`` is set to ``True``, then this function will return a
        tuple, whose first entry is a boolean and whose second entry is the
        Kuratowski subgraph or minor isolated by the Boyer-Myrvold
        algorithm. Note that this graph might contain a vertex or edges that
        were not in the initial graph. These would be elements referred to below
        as parts of the wheel and the star, which were added to the graph to
        require that the boundary can be drawn on the boundary of a disc, with
        all other vertices drawn inside (and no edge crossings). For more
        information, see [Kirkman]_.

        ALGORITHM:

        This is a linear time algorithm to test for circular planarity. It
        relies on the edge-addition planarity algorithm due to
        Boyer-Myrvold. We accomplish linear time for circular planarity by
        modifying the graph before running the general planarity
        algorithm.

        REFERENCE:

        .. [BM04] John M. Boyer and Wendy J. Myrvold, On the Cutting Edge:
          Simplified O(n) Planarity by Edge Addition. Journal of Graph
          Algorithms and Applications, Vol. 8, No. 3, pp. 241-273,
          2004.

        .. [Kirkman] Kirkman, Emily A. O(n) Circular Planarity Testing.
          [Online] Available: soon!

        EXAMPLES::

            sage: g439 = Graph({1:[5,7], 2:[5,6], 3:[6,7], 4:[5,6,7]})
            sage: g439.set_boundary([1,2,3,4])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: g439.show()
            sage: g439.is_circular_planar(boundary = [1,2,3,4])
            False
            sage: g439.is_circular_planar(kuratowski=True, boundary = [1,2,3,4])
            (False, Graph on 8 vertices)
            sage: g439.is_circular_planar(kuratowski=True, boundary = [1,2,3])
            (True, None)
            sage: g439.get_embedding()
            {1: [7, 5],
            2: [5, 6],
            3: [6, 7],
            4: [7, 6, 5],
            5: [1, 4, 2],
            6: [2, 4, 3],
            7: [3, 4, 1]}

        Order matters::

            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.is_circular_planar(boundary = [0,1,2,3])
            True
            sage: K23.is_circular_planar(ordered=True, boundary = [0,1,2,3])
            False

        With a different order::

            sage: K23.is_circular_planar(set_embedding=True, boundary = [0,2,1,3])
            True
        """
        if ordered and boundary is None:
            raise ValueError("boundary must be set when ordered is True.")

        if boundary is None:
            boundary = self

        # A local copy of self
        from sage.graphs.planarity import is_planar
        graph = self.to_undirected()
        if hasattr(graph, '_embedding'):
            del(graph._embedding)

        # Adds a new vertex to the graph
        extra = 0
        while graph.has_vertex(extra):
            extra=extra+1
        graph.add_vertex(extra)

        # Adds edges from the new vertex to all vertices of the boundary
        for vertex in boundary:
            graph.add_edge(vertex,extra)

        extra_edges = []

        # When ordered is True, we need a way to make sure that the ordering is
        # respected.
        if ordered:

            # We add edges between consecutive vertices of the boundary (only
            # len(boundary)-1 are actually sufficient)
            for u,v in zip(boundary[:-1],boundary[1:]):
                if not graph.has_edge(u,v):
                    extra_edges.append((u,v))

            graph.add_edges(extra_edges)

        result = is_planar(graph,kuratowski=kuratowski,set_embedding=set_embedding,circular=True)

        if kuratowski:
            bool_result = result[0]
        else:
            bool_result = result

        if bool_result:
            graph.delete_vertex(extra)
            graph.delete_edges(extra_edges)

            if hasattr(graph,'_embedding'):
                # strip the embedding to fit original graph
                for u,v in extra_edges:
                    graph._embedding[u].pop(graph._embedding[u].index(v))
                    graph._embedding[v].pop(graph._embedding[v].index(u))
                for w in boundary:
                    graph._embedding[w].pop(graph._embedding[w].index(extra))

                if set_embedding:
                    self._embedding = graph._embedding.copy()

            if (set_pos and set_embedding):
                self.set_planar_positions()

        del graph
        return result

    # TODO: rename into _layout_planar
    def set_planar_positions(self, test = False, **layout_options):
        """
        Compute a planar layout for self using Schnyder's algorithm,
        and save it as default layout.

        EXAMPLES::

            sage: g = graphs.CycleGraph(7)
            sage: g.set_planar_positions(test=True)
            True

        This method is deprecated since Sage-4.4.1.alpha2. Please use instead:

            sage: g.layout(layout = "planar", save_pos = True)
            {0: [1, 1], 1: [2, 2], 2: [3, 2], 3: [1, 4], 4: [5, 1], 5: [0, 5], 6: [1, 0]}
        """
        self.layout(layout = "planar", save_pos = True, test = test, **layout_options)
        if test:    # Optional error-checking, ( looking for edge-crossings O(n^2) ).
            return self.is_drawn_free_of_edge_crossings() # returns true if tests pass
        else:
            return

    def layout_planar(self, set_embedding=False, on_embedding=None, external_face=None, test=False, circular=False, **options):
        """
        Uses Schnyder's algorithm to compute a planar layout for self,
        raising an error if self is not planar.

        INPUT:

        -  ``set_embedding`` - if True, sets the combinatorial
           embedding used (see self.get_embedding())

        -  ``on_embedding`` - dict: provide a combinatorial
           embedding

        -  ``external_face`` - ignored

        -  ``test`` - if True, perform sanity tests along the way

        -  ``circular`` - ignored


        EXAMPLES::

            sage: g = graphs.PathGraph(10)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.BalancedTree(3,4)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.CycleGraph(7)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.CompleteGraph(5)
            sage: g.set_planar_positions(test=True,set_embedding=True)
            Traceback (most recent call last):
            ...
            ValueError: Complete graph is not a planar graph
        """
        from sage.graphs.schnyder import _triangulate, _normal_label, _realizer, _compute_coordinates

        G = self.to_undirected()
        try:
            G._embedding = self._embedding
        except AttributeError:
            pass
        embedding_copy = None
        if set_embedding:
            if not (G.is_planar(set_embedding=True)):
                raise ValueError('%s is not a planar graph'%self)
            embedding_copy = G._embedding
        else:
            if on_embedding is not None:
                if G._check_embedding_validity(on_embedding):
                    if not (G.is_planar(on_embedding=on_embedding)):
                        raise ValueError('provided embedding is not a planar embedding for %s'%self )
                else:
                    raise ValueError('provided embedding is not a valid embedding for %s. Try putting set_embedding=True'%self)
            else:
                if hasattr(G,'_embedding'):
                    if G._check_embedding_validity():
                        if not (G.is_planar(on_embedding=G._embedding)):
                            raise ValueError('%s has nonplanar _embedding attribute.  Try putting set_embedding=True'%self)
                        embedding_copy = G._embedding
                    else:
                        raise ValueError('provided embedding is not a valid embedding for %s. Try putting set_embedding=True'%self)
                else:
                    G.is_planar(set_embedding=True)

        # The following is what was breaking the code.  It is where we were specifying the external
        #       face ahead of time.  This is definitely a TODO:
        #
        # Running is_planar(set_embedding=True) has set attribute self._embedding
        #if external_face is None:
        #    faces = faces( self, self._embedding )
        #    faces.sort(key=len)
        #    external_face = faces[-1]

        #n = len(external_face)
        #other_added_edges = []
        #if n > 3:
        #    v1, v2, v3 = external_face[0][0], external_face[int(n/3)][0], external_face[int(2*n/3)][0]
        #    if not self.has_edge( (v1,v2) ):
        #        self.add_edge( (v1, v2) )
        #        other_added_edges.append( (v1, v2) )
        #    if not self.has_edge( (v2,v3) ):
        #        self.add_edge( (v2, v3) )
        #        other_added_edges.append( (v2, v3) )
        #    if not self.has_edge( (v3,v1) ):
        #        self.add_edge( (v3, v1) )
        #        other_added_edges.append( (v3, v1) )
        #    if not self.is_planar(set_embedding=True): # get new combinatorial embedding (with added edges)
        #        raise ValueError('modified graph %s is not planar.  Try specifying an external face'%self)

        # Triangulate the graph
        extra_edges = _triangulate( G, G._embedding)

        # Optional error-checking
        if test:
            G.is_planar(set_embedding=True) # to get new embedding
            test_faces = G.faces(G._embedding)
            for face in test_faces:
                if len(face) != 3:
                    raise RuntimeError('BUG: Triangulation returned face: %s'%face)

        G.is_planar(set_embedding=True)
        faces = G.faces(G._embedding)
        # Assign a normal label to the graph
        label = _normal_label( G, G._embedding, faces[0] )

        # Get dictionary of tree nodes from the realizer
        tree_nodes = _realizer( G, label)

        # Compute the coordinates and store in position dictionary (attr self._pos)
        _compute_coordinates( G, tree_nodes )

        # Delete all the edges added to the graph
        #G.delete_edges( extra_edges )
        #self.delete_edges( other_added_edges )

        if embedding_copy is not None:
            self._embedding = embedding_copy

        return G._pos

    def is_drawn_free_of_edge_crossings(self):
        """
        Returns True is the position dictionary for this graph is set and
        that position dictionary gives a planar embedding.

        This simply checks all pairs of edges that don't share a vertex to
        make sure that they don't intersect.

        .. note::

           This function require that _pos attribute is set. (Returns
           False otherwise.)

        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: D.set_planar_positions()
            sage: D.is_drawn_free_of_edge_crossings()
            True
        """
        if self._pos is None:
            return False

        G = self.to_undirected()
        for edge1 in G.edges(labels = False):
            for edge2 in G.edges(labels = False):
                if edge1[0] == edge2[0] or edge1[0] == edge2[1] or edge1[1] == edge2[0] or edge1[1] == edge2[1]:
                    continue
                p1, p2 = self._pos[edge1[0]], self._pos[edge1[1]]
                dy = Rational(p2[1] - p1[1])
                dx = Rational(p2[0] - p1[0])
                q1, q2 = self._pos[edge2[0]], self._pos[edge2[1]]
                db = Rational(q2[1] - q1[1])
                da = Rational(q2[0] - q1[0])
                if(da * dy == db * dx):
                    if dx != 0:
                        t1 = Rational(q1[0] - p1[0])/dx
                        t2 = Rational(q2[0] - p1[0])/dx
                        if (0 <= t1 and t1 <= 1) or (0 <= t2 and t2 <= 1):
                            if p1[1] + t1 * dy == q1[1] or p1[1] + t2 * dy == q2[1]:
                                return False
                    else:
                        t1 = Rational(q1[1] - p1[1])/dy
                        t2 = Rational(q2[1] - p1[1])/dy
                        if (0 <= t1 and t1 <= 1) or (0 <= t2 and t2 <= 1):
                            if p1[0] + t1 * dx == q1[0] or p1[0] + t2 * dx == q2[0]:
                                return False
                else:
                    s = (dx * Rational(q1[1] - p1[1]) + dy * Rational(p1[0] - q1[0])) / (da * dy - db * dx)
                    t = (da * Rational(p1[1] - q1[1]) + db * Rational(q1[0] - p1[0])) / (db * dx - da * dy)

                    if s >= 0 and s <= 1 and t >= 0 and t <= 1:
                        print 'fail on', p1, p2, ' : ',q1, q2
                        print edge1, edge2
                        return False
        return True

    def genus(self, set_embedding=True, on_embedding=None, minimal=True, maximal=False, circular=False, ordered=True):
        """
        Returns the minimal genus of the graph.

        The genus of a compact surface is the number of handles it
        has. The genus of a graph is the minimal genus of the surface
        it can be embedded into.

        Note - This function uses Euler's formula and thus it is necessary
        to consider only connected graphs.

        INPUT:

        -  ``set_embedding (boolean)`` - whether or not to
           store an embedding attribute of the computed (minimal) genus of the
           graph. (Default is True).

        -  ``on_embedding (dict)`` - a combinatorial embedding
           to compute the genus of the graph on. Note that this must be a
           valid embedding for the graph. The dictionary structure is given
           by: vertex1: [neighbor1, neighbor2, neighbor3], vertex2: [neighbor]
           where there is a key for each vertex in the graph and a (clockwise)
           ordered list of each vertex's neighbors as values. on_embedding
           takes precedence over a stored _embedding attribute if minimal is
           set to False. Note that as a shortcut, the user can enter
           on_embedding=True to compute the genus on the current _embedding
           attribute. (see eg's.)

        -  ``minimal (boolean)`` - whether or not to compute
           the minimal genus of the graph (i.e., testing all embeddings). If
           minimal is False, then either maximal must be True or on_embedding
           must not be None. If on_embedding is not None, it will take
           priority over minimal. Similarly, if maximal is True, it will take
           priority over minimal.

        -  ``maximal (boolean)`` - whether or not to compute
           the maximal genus of the graph (i.e., testing all embeddings). If
           maximal is False, then either minimal must be True or on_embedding
           must not be None. If on_embedding is not None, it will take
           priority over maximal. However, maximal takes priority over the
           default minimal.

        -  ``circular (boolean)`` - whether or not to compute
           the genus preserving a planar embedding of the boundary. (Default
           is False). If circular is True, on_embedding is not a valid
           option.

        -  ``ordered (boolean)`` - if circular is True, then
           whether or not the boundary order may be permuted. (Default is
           True, which means the boundary order is preserved.)

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.genus() # tests for minimal genus by default
            1
            sage: g.genus(on_embedding=True, maximal=True) # on_embedding overrides minimal and maximal arguments
            1
            sage: g.genus(maximal=True) # setting maximal to True overrides default minimal=True
            3
            sage: g.genus(on_embedding=g.get_embedding()) # can also send a valid combinatorial embedding dict
            3
            sage: (graphs.CubeGraph(3)).genus()
            0
            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.genus()
            0
            sage: K33 = graphs.CompleteBipartiteGraph(3,3)
            sage: K33.genus()
            1

        Using the circular argument, we can compute the minimal genus
        preserving a planar, ordered boundary::

            sage: cube = graphs.CubeGraph(2)
            sage: cube.set_boundary(['01','10'])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: cube.genus()
            0
            sage: cube.is_circular_planar()
            True
            sage: cube.genus(circular=True)
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            0
            sage: cube.genus(circular=True, on_embedding=True)
            0
            sage: cube.genus(circular=True, maximal=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute the maximal genus of a genus respecting a boundary.

        Note: not everything works for multigraphs, looped graphs or digraphs.  But the
        minimal genus is ultimately computable for every connected graph -- but the
        embedding we obtain for the simple graph can't be easily converted to an
        embedding of a non-simple graph.  Also, the maximal genus of a multigraph does
        not trivially correspond to that of its simple graph. ::

            sage: G = DiGraph({ 0 : [0,1,1,1], 1 : [2,2,3,3], 2 : [1,3,3], 3:[0,3]})
            sage: G.genus()
            Traceback (most recent call last):
            ...
            NotImplementedError: Can't work with embeddings of non-simple graphs
            sage: G.to_simple().genus()
            0
            sage: G.genus(set_embedding=False)
            0
            sage: G.genus(maximal=True, set_embedding=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: Can't compute the maximal genus of a graph with loops or multiple edges


        We break graphs with cut vertices into their blocks, which greatly speeds up
        computation of minimal genus.  This is not implemented for maximal genus. ::

            sage: K5 = graphs.CompleteGraph(5)
            sage: G = K5.copy()
            sage: s = 4
            sage: for i in range(1,100):
            ...       k = K5.relabel(range(s,s+5),inplace=False)
            ...       G.add_edges(k.edges())
            ...       s += 4
            ...
            sage: G.genus()
            100

        """
        if not self.is_connected():
            raise TypeError("Graph must be connected to use Euler's Formula to compute minimal genus.")

        G = self.to_simple()
        verts = G.order()
        edges = G.size()

        if maximal:
            minimal = False

        if circular:
            if maximal:
                raise NotImplementedError("Cannot compute the maximal genus of a genus respecting a boundary.")
            boundary = G.get_boundary()
            if hasattr(G, '_embedding'):
                del(G._embedding)

            extra = 0
            while G.has_vertex(extra):
                extra=extra+1
            G.add_vertex(extra)
            verts += 1

            for vertex in boundary:
                G.add_edge(vertex,extra)

            extra_edges = []
            if ordered: # WHEEL
                for i in range(len(boundary)-1):
                    if not G.has_edge(boundary[i],boundary[i+1]):
                        G.add_edge(boundary[i],boundary[i+1])
                        extra_edges.append((boundary[i],boundary[i+1]))
                if not G.has_edge(boundary[-1],boundary[0]):
                    G.add_edge(boundary[-1],boundary[0])
                    extra_edges.append((boundary[-1],boundary[0]))
                # else STAR (empty list of extra edges)

            edges = G.size()

        if on_embedding is not None:
            if self.has_loops() or self.is_directed() or self.has_multiple_edges():
                raise NotImplementedError("Can't work with embeddings of non-simple graphs")
            if on_embedding: #i.e., if on_embedding True (returns False if on_embedding is of type dict)
                try:
                    faces = len(self.faces(self._embedding))
                except AttributeError:
                    raise AttributeError('graph must have attribute _embedding set to compute current (embedded) genus')
                return (2-verts+edges-faces) // 2
            else: # compute genus on the provided dict
                faces = len(self.faces(on_embedding))
                return (2-verts+edges-faces) // 2
        else: # then compute either maximal or minimal genus of all embeddings
            import genus

            if set_embedding:
                if self.has_loops() or self.is_directed() or self.has_multiple_edges():
                    raise NotImplementedError("Can't work with embeddings of non-simple graphs")
                if minimal:
                    B,C = G.blocks_and_cut_vertices()
                    embedding = {}
                    g = 0
                    for block in B:
                        H = G.subgraph(block)
                        g += genus.simple_connected_graph_genus(H, set_embedding = True, check = False, minimal = True)
                        emb = H.get_embedding()
                        for v in emb:
                            if v in embedding:
                                embedding[v] += emb[v]
                            else:
                                embedding[v] = emb[v]
                    self._embedding = embedding
                else:
                    g = genus.simple_connected_graph_genus(G, set_embedding = True, check = False, minimal = minimal)
                    self._embedding = G._embedding
                return g
            else:
                if maximal and (self.has_multiple_edges() or self.has_loops()):
                    raise NotImplementedError("Can't compute the maximal genus of a graph with loops or multiple edges")
                if minimal:
                    B,C = G.blocks_and_cut_vertices()
                    g = 0
                    for block in B:
                        H = G.subgraph(block)
                        g += genus.simple_connected_graph_genus(H, set_embedding = False, check = False, minimal = True)
                    return g
                else:
                    return genus.simple_connected_graph_genus(G, set_embedding = False, check=False, minimal=minimal)

    def faces(self, embedding = None):
        """
        Return the faces of an embedded graph.

        A combinatorial embedding of a graph is a clockwise ordering of the
        neighbors of each vertex. From this information one can define the faces
        of the embedding, which is what this method returns.

        INPUT:

        - ``embedding`` - a combinatorial embedding dictionary. Format:
          ``{v1:[v2,v3], v2:[v1], v3:[v1]}`` (clockwise ordering of neighbors at
          each vertex). If set to ``None`` (default) the method will use the
          embedding stored as ``self._embedding``. If none is stored, the method
          will compute the set of faces from the embedding returned by
          :meth:`is_planar` (if the graph is, of course, planar).

        .. NOTE::

            ``embedding`` is an ordered list based on the hash order of the
            vertices of graph. To avoid confusion, it might be best to set the
            rot_sys based on a 'nice_copy' of the graph.

        .. SEEALSO::

            * :meth:`set_embedding`
            * :meth:`get_embedding`
            * :meth:`is_planar`

        EXAMPLES::

            sage: T = graphs.TetrahedralGraph()
            sage: T.faces({0: [1, 3, 2], 1: [0, 2, 3], 2: [0, 3, 1], 3: [0, 1, 2]})
            [[(0, 1), (1, 2), (2, 0)],
             [(3, 2), (2, 1), (1, 3)],
             [(2, 3), (3, 0), (0, 2)],
             [(0, 3), (3, 1), (1, 0)]]

        With no embedding provided::

            sage: graphs.TetrahedralGraph().faces()
            [[(0, 1), (1, 2), (2, 0)],
             [(3, 2), (2, 1), (1, 3)],
             [(2, 3), (3, 0), (0, 2)],
             [(0, 3), (3, 1), (1, 0)]]

        With no embedding provided (non-planar graph)::

            sage: graphs.PetersenGraph().faces()
            Traceback (most recent call last):
            ...
            ValueError: No embedding is provided and the graph is not planar.

        TESTS:

        :trac:`15551` deprecated the ``trace_faces`` name::

            sage: T.trace_faces({0: [1, 3, 2], 1: [0, 2, 3], 2: [0, 3, 1], 3: [0, 1, 2]})
            doctest:...: DeprecationWarning: trace_faces is deprecated. Please use faces instead.
            See http://trac.sagemath.org/15551 for details.
            [[(0, 1), (1, 2), (2, 0)], [(3, 2), (2, 1), (1, 3)], [(2, 3), (3, 0), (0, 2)], [(0, 3), (3, 1), (1, 0)]]

        """
        # Which embedding should we use ?
        if embedding is None:
            # Is self._embedding available ?
            if self._check_embedding_validity():
                embedding = self._embedding
            else:
                if self.is_planar(set_embedding=True):
                    embedding = self._embedding
                    self._embedding = None
                else:
                    raise ValueError("No embedding is provided and the graph is not planar.")

        from sage.sets.set import Set

        # Establish set of possible edges
        edgeset = Set([])
        for edge in self.to_undirected().edges():
            edgeset = edgeset.union(Set([(edge[0],edge[1]),(edge[1],edge[0])]))

        # Storage for face paths
        faces = []
        path = []
        for edge in edgeset:
            path.append(edge)
            edgeset -= Set([edge])
            break  # (Only one iteration)

        # Trace faces
        while (len(edgeset) > 0):
            neighbors = embedding[path[-1][-1]]
            next_node = neighbors[(neighbors.index(path[-1][-2])+1)%(len(neighbors))]
            tup = (path[-1][-1],next_node)
            if tup == path[0]:
                faces.append(path)
                path = []
                for edge in edgeset:
                    path.append(edge)
                    edgeset -= Set([edge])
                    break  # (Only one iteration)
            else:
                path.append(tup)
                edgeset -= Set([tup])
        if (len(path) != 0): faces.append(path)
        return faces

    trace_faces = deprecated_function_alias(15551, faces)

    ### Connectivity

    def is_connected(self):
        """
        Indicates whether the (di)graph is connected. Note that in a graph,
        path connected is equivalent to connected.

        EXAMPLES::

            sage: G = Graph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: G.is_connected()
            False
            sage: G.add_edge(0,3)
            sage: G.is_connected()
            True
            sage: D = DiGraph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: D.is_connected()
            False
            sage: D.add_edge(0,3)
            sage: D.is_connected()
            True
            sage: D = DiGraph({1:[0], 2:[0]})
            sage: D.is_connected()
            True
        """
        if self.order() == 0:
            return True

        try:
            return self._backend.is_connected()
        except AttributeError:
            v = self.vertex_iterator().next()
            conn_verts = list(self.depth_first_search(v, ignore_direction=True))
            return len(conn_verts) == self.num_verts()

    def connected_components(self):
        """
        Returns the list of connected components.

        Returns a list of lists of vertices, each list representing a
        connected component. The list is ordered from largest to smallest
        component.

        EXAMPLES::

            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]
        """
        seen = set()
        components = []
        for v in self:
            if v not in seen:
                c = self.connected_component_containing_vertex(v)
                seen.update(c)
                components.append(c)
        components.sort(lambda comp1, comp2: cmp(len(comp2), len(comp1)))
        return components

    def connected_components_number(self):
        """
        Returns the number of connected components.

        EXAMPLES::

            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components_number()
            2
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components_number()
            2
        """
        return len(self.connected_components())

    def connected_components_subgraphs(self):
        """
        Returns a list of connected components as graph objects.

        EXAMPLES::

            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = G.connected_components_subgraphs()
            sage: graphs_list.show_graphs(L)
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = D.connected_components_subgraphs()
            sage: graphs_list.show_graphs(L)
        """
        cc = self.connected_components()
        list = []
        for c in cc:
            list.append(self.subgraph(c, inplace=False))
        return list

    def connected_component_containing_vertex(self, vertex):
        """
        Returns a list of the vertices connected to vertex.

        EXAMPLES::

            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_component_containing_vertex(0)
            [0, 1, 2, 3]
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_component_containing_vertex(0)
            [0, 1, 2, 3]
        """
        try:
            c = list(self._backend.depth_first_search(vertex, ignore_direction=True))
        except AttributeError:
            c = list(self.depth_first_search(vertex, ignore_direction=True))

        c.sort()
        return c

    def blocks_and_cut_vertices(self):
        """
        Computes the blocks and cut vertices of the graph.

        In the case of a digraph, this computation is done on the underlying
        graph.

        A cut vertex is one whose deletion increases the number of
        connected components. A block is a maximal induced subgraph which
        itself has no cut vertices. Two distinct blocks cannot overlap in
        more than a single cut vertex.

        OUTPUT: ``( B, C )``, where ``B`` is a list of blocks- each is
        a list of vertices and the blocks are the corresponding induced
        subgraphs-and ``C`` is a list of cut vertices.

        ALGORITHM:

          We implement the algorithm proposed by Tarjan in [Tarjan72]_. The
          original version is recursive. We emulate the recursion using a stack.

        .. SEEALSO:: :meth:`blocks_and_cuts_tree`

        EXAMPLES::

            sage: graphs.PetersenGraph().blocks_and_cut_vertices()
            ([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]], [])
            sage: graphs.PathGraph(6).blocks_and_cut_vertices()
            ([[4, 5], [3, 4], [2, 3], [1, 2], [0, 1]], [1, 2, 3, 4])
            sage: graphs.CycleGraph(7).blocks_and_cut_vertices()
            ([[0, 1, 2, 3, 4, 5, 6]], [])
            sage: graphs.KrackhardtKiteGraph().blocks_and_cut_vertices()
            ([[8, 9], [7, 8], [0, 1, 2, 3, 4, 5, 6, 7]], [7, 8])
            sage: G=Graph()  # make a bowtie graph where 0 is a cut vertex
            sage: G.add_vertices(range(5))
            sage: G.add_edges([(0,1),(0,2),(0,3),(0,4),(1,2),(3,4)])
            sage: G.blocks_and_cut_vertices()
            ([[0, 1, 2], [0, 3, 4]], [0])
            sage: graphs.StarGraph(3).blocks_and_cut_vertices()
            ([[0, 1], [0, 2], [0, 3]], [0])

        TESTS::

            sage: Graph(0).blocks_and_cut_vertices()
            ([], [])
            sage: Graph(1).blocks_and_cut_vertices()
            ([[0]], [])
            sage: Graph(2).blocks_and_cut_vertices()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...

        REFERENCE:

        .. [Tarjan72] R.E. Tarjan. Depth-First Search and Linear Graph
          Algorithms. SIAM J. Comput. 1(2): 146-160 (1972).
        """
        if not self: # empty graph
            return [],[]

        start = self.vertex_iterator().next() # source

        if len(self) == 1: # only one vertex
            return [[start]],[]

        if not self.is_connected():
            raise NotImplementedError("Blocks and cut vertices is currently only implemented for connected graphs.")

        # Each vertex is number with an integer from 1...|V(G)|, corresponding
        # to the order in which it is discovered during the DFS.
        number = {}
        num = 1

        # Associates to each vertex v the smallest number of a vertex that can
        # be reached from v in the orientation of the graph that the algorithm
        # creates.
        low_point = {}

        # Associates to each vertex an iterator over its neighbors
        neighbors = {}

        blocks = []
        cut_vertices = set()

        stack = [start]
        edge_stack = []
        start_already_seen = False

        while stack:
            v = stack[-1]

            # The first time we meet v
            if not v in number:
                # We number the vertices in the order they are reached during
                # DFS
                number[v] = num
                neighbors[v] = self.neighbor_iterator(v)
                low_point[v] = num
                num += 1

            try:
                # We consider the next of its neighbors
                w = neighbors[v].next()

                # If we never met w before, we remember the direction of edge
                # vw, and add w to the stack.
                if not w in number:
                    edge_stack.append( (v,w) )
                    stack.append(w)

                # If w is an ancestor of v in the DFS tree, we remember the
                # direction of edge vw
                elif number[w]<number[v]:
                    edge_stack.append( (v,w) )
                    low_point[v] = min(low_point[v], number[w])

            # We went through all of v's neighbors
            except StopIteration:
                # We trackback, so w takes the value of v and we pop the stack
                w = stack.pop()

                # Test termination of the algorithm
                if not stack:
                    break

                v = stack[-1]

                # Propagating the information : low_point[v] indicates the
                # smallest vertex (the vertex x with smallest number[x]) that
                # can be reached from v
                low_point[v] = min(low_point[v], low_point[w])

                # The situation in which there is no path from w to an ancestor
                # of v : we have identified a new biconnected component
                if low_point[w] >= number[v]:
                    new_block = set()
                    nw = number[w]
                    u1,u2 = edge_stack.pop()
                    while number[u1] >= nw:
                        new_block.add(u1)
                        u1,u2 = edge_stack.pop()
                    new_block.add(u1)
                    blocks.append(sorted(list(new_block)))

                    # We update the set of cut vertices.
                    #
                    # If v is start, then we add it only if it belongs to
                    # several blocks.
                    if (not v is start) or start_already_seen:
                        cut_vertices.add(v)
                    else:
                        start_already_seen = True

        return blocks,sorted(list(cut_vertices))

    def blocks_and_cuts_tree(self):
        """
        Returns the blocks-and-cuts tree of ``self``.

        This new graph has two different kinds of vertices, some representing
        the blocks (type B) and some other the cut vertices of the graph
        ``self`` (type C).

        There is an edge between a vertex `u` of type B and a vertex `v` of type
        C if the cut-vertex corresponding to `v` is in the block corresponding
        to `u`.

        The resulting graph is a tree, with the additional characteristic
        property that the distance between two leaves is even.

        When ``self`` is biconnected, the tree is reduced to a single node of
        type `B`.

        .. SEEALSO:: :meth:`blocks_and_cut_vertices`

        EXAMPLES:

            sage: T = graphs.KrackhardtKiteGraph().blocks_and_cuts_tree(); T
            Graph on 5 vertices
            sage: T.is_isomorphic(graphs.PathGraph(5))
            True

        The distance between two leaves is even::

            sage: T = graphs.RandomTree(40).blocks_and_cuts_tree()
            sage: T.is_tree()
            True
            sage: leaves = [v for v in T if T.degree(v) == 1]
            sage: all(T.distance(u,v) % 2 == 0 for u in leaves for v in leaves)
            True

        The tree of a biconnected graph has a single vertex, of type `B`::

            sage: T = graphs.PetersenGraph().blocks_and_cuts_tree()
            sage: T.vertices()
            [('B', (0, 1, 2, 3, 4, 5, 6, 7, 8, 9))]

        REFERENCES:

        .. [HarPri] F. Harary and G. Prins. The block-cutpoint-tree of
           a graph. Publ. Math. Debrecen 13 1966 103-107.
        .. [Gallai] T. Gallai, Elementare Relationen bezueglich der
           Glieder und trennenden Punkte von Graphen, Magyar
           Tud. Akad. Mat. Kutato Int. Kozl. 9 (1964) 235-236
        """
        from sage.graphs.graph import Graph
        B, C = self.blocks_and_cut_vertices()
        B = map(tuple, B)
        G = Graph()
        for bloc in B:
            G.add_vertex(('B', bloc))
            for c in bloc:
                if c in C:
                    G.add_edge(('B', bloc), ('C', c))
        return G

    def is_cut_edge(self, u, v=None, label=None):
        """
        Returns True if the input edge is a cut-edge or a bridge.

        A cut edge (or bridge) is an edge that when removed increases
        the number of connected components.  This function works with
        simple graphs as well as graphs with loops and multiedges.  In
        a digraph, a cut edge is an edge that when removed increases
        the number of (weakly) connected components.

        INPUT: The following forms are accepted

        - G.is_cut_edge( 1, 2 )

        - G.is_cut_edge( (1, 2) )

        - G.is_cut_edge( 1, 2, 'label' )

        - G.is_cut_edge( (1, 2, 'label') )

        OUTPUT:

        - Returns True if (u,v) is a cut edge, False otherwise

        EXAMPLES::

            sage: G = graphs.CompleteGraph(4)
            sage: G.is_cut_edge(0,2)
            False

            sage: G = graphs.CompleteGraph(4)
            sage: G.add_edge((0,5,'silly'))
            sage: G.is_cut_edge((0,5,'silly'))
            True

            sage: G = Graph([[0,1],[0,2],[3,4],[4,5],[3,5]])
            sage: G.is_cut_edge((0,1))
            True

            sage: G = Graph([[0,1],[0,2],[1,1]])
            sage: G.allow_loops(True)
            sage: G.is_cut_edge((1,1))
            False

            sage: G = digraphs.Circuit(5)
            sage: G.is_cut_edge((0,1))
            False

            sage: G = graphs.CompleteGraph(6)
            sage: G.is_cut_edge((0,7))
            Traceback (most recent call last):
            ...
            ValueError: edge not in graph
        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except ValueError:
                    u, v = u
                    label = None

        if not self.has_edge(u,v):
            raise ValueError('edge not in graph')

        # If edge (u,v) is a pending edge, it is also a cut-edge
        if self.degree(u) == 1 or self.degree(v) == 1:
            return True
        elif self.allows_multiple_edges():
            # If we have two or more edges between u and v, it is not a cut-edge
            if len([(uu,vv) for uu,vv,ll in self.edges_incident(u) if uu == v or vv == v]) > 1:
                return False

        self.delete_edge(u,v,label)
        if self.is_directed():
            # (u,v) is a cut-edge if u is not in the connected
            # component containing v of self-(u,v)
            sol = not u in self.connected_component_containing_vertex(v)
        else:
            # (u,v) is a cut-edge if there is no path from u to v in
            # self-(u,v)
            sol = not self.distance(u,v) < self.order()

        self.add_edge(u,v,label)
        return sol


    def is_cut_vertex(self, u, weak=False):
        r"""
        Returns True if the input vertex is a cut-vertex.

        A vertex is a cut-vertex if its removal from the (di)graph increases the
        number of (strongly) connected components. Isolated vertices or leafs
        are not cut-vertices. This function works with simple graphs as well as
        graphs with loops and multiple edges.

        INPUT:

        - ``u`` -- a vertex

        - ``weak`` -- (default: ``False``) boolean set to `True` if the
          connectivity of directed graphs is to be taken in the weak sense, that
          is ignoring edges orientations.

        OUTPUT:

        Returns True if ``u`` is a cut-vertex, and False otherwise.

        EXAMPLES:

        Giving a LollipopGraph(4,2), that is a complete graph with 4 vertices with a pending edge::

            sage: G = graphs.LollipopGraph(4,2)
            sage: G.is_cut_vertex(0)
            False
            sage: G.is_cut_vertex(3)
            True

        Comparing the weak and strong connectivity of a digraph::

            sage: D = digraphs.Circuit(6)
            sage: D.is_strongly_connected()
            True
            sage: D.is_cut_vertex(2)
            True
            sage: D.is_cut_vertex(2, weak=True)
            False

        Giving a vertex that is not in the graph::

            sage: G = graphs.CompleteGraph(6)
            sage: G.is_cut_vertex(7)
            Traceback (most recent call last):
            ...
            ValueError: The input vertex is not in the vertex set.
        """
        if not u in self:
            raise ValueError('The input vertex is not in the vertex set.')

        # Initialization
        if not self.is_directed() or weak:
            # Weak connectivity

            if self.degree(u) < 2:
                # An isolated or a leaf vertex is not a cut vertex
                return False

            neighbors_func = [self.neighbor_iterator]
            start = self.neighbor_iterator(u).next()
            CC = set(self.vertex_iterator())

        else:
            # Strong connectivity for digraphs

            if self.out_degree(u) == 0 or self.in_degree(u) == 0:
                # A vertex without in or out neighbors is not a cut vertex
                return False

            # We consider only the strongly connected component containing u
            CC = set(self.strongly_connected_component_containing_vertex(u))

            # We perform two DFS starting from an out neighbor of u and avoiding
            # u. The first DFS follows the edges directions, and the second is
            # in the reverse order. If both allow to reach all neighbors of u,
            # then u is not a cut vertex
            neighbors_func = [self.neighbor_out_iterator, self.neighbor_in_iterator]
            start = self.neighbor_out_iterator(u).next()

        CC.discard(u)
        CC.discard(start)
        for neighbors in neighbors_func:

            # We perform a DFS starting from a neighbor of u and avoiding u
            queue = [start]
            seen = set(queue)
            targets = set(self.neighbor_iterator(u))&CC
            targets.discard(start)
            while queue and targets:
                v = queue.pop()
                for w in neighbors(v):
                    if not w in seen and w in CC:
                        seen.add(w)
                        queue.append(w)
                        targets.discard(w)

            # If some neighbors cannot be reached, u is a cut vertex.
            if targets:
                return True

        return False


    def steiner_tree(self,vertices, weighted = False, solver = None, verbose = 0):
        r"""
        Returns a tree of minimum weight connecting the given
        set of vertices.

        Definition :

        Computing a minimum spanning tree in a graph can be done in `n
        \log(n)` time (and in linear time if all weights are equal) where
        `n = V + E`. On the other hand, if one is given a large (possibly
        weighted) graph and a subset of its vertices, it is NP-Hard to
        find a tree of minimum weight connecting the given set of
        vertices, which is then called a Steiner Tree.

        `Wikipedia article on Steiner Trees
        <http://en.wikipedia.org/wiki/Steiner_tree_problem>`_.

        INPUT:

        - ``vertices`` -- the vertices to be connected by the Steiner
          Tree.

        - ``weighted`` (boolean) -- Whether to consider the graph as
          weighted, and use each edge's label as a weight, considering
          ``None`` as a weight of `1`. If ``weighted=False`` (default)
          all edges are considered to have a weight of `1`.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.


        .. NOTE::

            * This problem being defined on undirected graphs, the
              orientation is not considered if the current graph is
              actually a digraph.

            * The graph is assumed not to have multiple edges.

        ALGORITHM:

        Solved through Linear Programming.

        COMPLEXITY:

        NP-Hard.

        Note that this algorithm first checks whether the given
        set of vertices induces a connected graph, returning one of its
        spanning trees if ``weighted`` is set to ``False``, and thus
        answering very quickly in some cases

        EXAMPLES:

        The Steiner Tree of the first 5 vertices in a random graph is,
        of course, always a tree ::

            sage: g = graphs.RandomGNP(30,.5)
            sage: st = g.steiner_tree(g.vertices()[:5])
            sage: st.is_tree()
            True

        And all the 5 vertices are contained in this tree ::

            sage: all([v in st for v in g.vertices()[:5] ])
            True

        An exception is raised when the problem is impossible, i.e.
        if the given vertices are not all included in the same
        connected component ::

            sage: g = 2 * graphs.PetersenGraph()
            sage: st = g.steiner_tree([5,15])
            Traceback (most recent call last):
            ...
            EmptySetError: The given vertices do not all belong to the same connected component. This problem has no solution !
        """
        self._scream_if_not_simple(allow_loops=True)

        if self.is_directed():
            from sage.graphs.all import Graph
            g = Graph(self)
        else:
            g = self

        # Can the problem be solved ? Are all the vertices in the same
        # connected component ?
        cc = g.connected_component_containing_vertex(vertices[0])
        if not all([v in cc for v in vertices]):
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("The given vertices do not all belong to the same connected component. This problem has no solution !")

        # Can it be solved using the min spanning tree algorithm ?
        if not weighted:
            gg = g.subgraph(vertices)
            if gg.is_connected():
                st = g.subgraph(edges = gg.min_spanning_tree())
                st.delete_vertices([v for v in g if st.degree(v) == 0])
                return st

        # Then, LP formulation
        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(maximization = False, solver = solver)

        # Reorder an edge
        R = lambda x_y: x_y if x_y[0] < x_y[1] else (x_y[1], x_y[0])

        # edges used in the Steiner Tree
        edges = p.new_variable(binary=True)

        # relaxed edges to test for acyclicity
        r_edges = p.new_variable(nonnegative=True)

        # Whether a vertex is in the Steiner Tree
        vertex = p.new_variable(binary = True)
        for v in g:
            for e in g.edges_incident(v, labels=False):
                p.add_constraint(vertex[v] - edges[R(e)], min = 0)

        # We must have the given vertices in our tree
        for v in vertices:
            p.add_constraint(p.sum([edges[R(e)] for e in g.edges_incident(v,labels=False)]), min=1)

        # The number of edges is equal to the number of vertices in our tree minus 1
        p.add_constraint(p.sum([vertex[v] for v in g]) - p.sum([edges[R(e)] for e in g.edges(labels=None)]), max = 1, min = 1)

        # There are no cycles in our graph

        for u,v in g.edges(labels = False):
            p.add_constraint( r_edges[(u,v)]+ r_edges[(v,u)] - edges[R((u,v))] , min = 0 )

        eps = 1/(5*Integer(g.order()))
        for v in g:
            p.add_constraint(p.sum([r_edges[(u,v)] for u in g.neighbors(v)]), max = 1-eps)


        # Objective
        if weighted:
            w = lambda x_y: g.edge_label(x_y[0], x_y[1]) if g.edge_label(x_y[0], x_y[1]) is not None else 1
        else:
            w = lambda x_y: 1

        p.set_objective(p.sum([w(e)*edges[R(e)] for e in g.edges(labels = False)]))

        p.solve(log = verbose)

        edges = p.get_values(edges)

        st =  g.subgraph(edges=[e for e in g.edges(labels = False) if edges[R(e)] == 1])
        st.delete_vertices([v for v in g if st.degree(v) == 0])
        return st

    def edge_disjoint_spanning_trees(self,k, root=None, solver = None, verbose = 0):
        r"""
        Returns the desired number of edge-disjoint spanning
        trees/arborescences.

        INPUT:

        - ``k`` (integer) -- the required number of edge-disjoint
          spanning trees/arborescences

        - ``root`` (vertex) -- root of the disjoint arborescences
          when the graph is directed.
          If set to ``None``, the first vertex in the graph is picked.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        ALGORITHM:

        Mixed Integer Linear Program. The formulation can be found
        in [LPForm]_.

        There are at least two possible rewritings of this method
        which do not use Linear Programming:

            * The algorithm presented in the paper entitled "A short
              proof of the tree-packing theorem", by Thomas Kaiser
              [KaisPacking]_.

            * The implementation of a Matroid class and of the Matroid
              Union Theorem (see section 42.3 of [SchrijverCombOpt]_),
              applied to the cycle Matroid (see chapter 51 of
              [SchrijverCombOpt]_).

        EXAMPLES:

        The Petersen Graph does have a spanning tree (it is connected)::

            sage: g = graphs.PetersenGraph()
            sage: [T] = g.edge_disjoint_spanning_trees(1)
            sage: T.is_tree()
            True

        Though, it does not have 2 edge-disjoint trees (as it has less
        than `2(|V|-1)` edges)::

            sage: g.edge_disjoint_spanning_trees(2)
            Traceback (most recent call last):
            ...
            EmptySetError: This graph does not contain the required number of trees/arborescences !

        By Edmond's theorem, a graph which is `k`-connected always has `k` edge-disjoint
        arborescences, regardless of the root we pick::

            sage: g = digraphs.RandomDirectedGNP(28,.3) # reduced from 30 to 28, cf. #9584
            sage: k = Integer(g.edge_connectivity())
            sage: arborescences = g.edge_disjoint_spanning_trees(k)  # long time (up to 15s on sage.math, 2011)
            sage: all([a.is_directed_acyclic() for a in arborescences])  # long time
            True
            sage: all([a.is_connected() for a in arborescences])  # long time
            True

        In the undirected case, we can only ensure half of it::

            sage: g = graphs.RandomGNP(30,.3)
            sage: k = floor(Integer(g.edge_connectivity())/2)
            sage: trees = g.edge_disjoint_spanning_trees(k)
            sage: all([t.is_tree() for t in trees])
            True

        REFERENCES:

        .. [LPForm] Nathann Cohen,
          Several Graph problems and their Linear Program formulations,
          http://hal.archives-ouvertes.fr/inria-00504914/en

        .. [KaisPacking] Thomas Kaiser
          A short proof of the tree-packing theorem
          :arxiv:`0911.2809`

        .. [SchrijverCombOpt] Alexander Schrijver
          Combinatorial optimization: polyhedra and efficiency
          2003
        """

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = MixedIntegerLinearProgram(solver = solver)
        p.set_objective(None)

        # The colors we can use
        colors = range(0,k)

        # edges[j,e] is equal to one if and only if edge e belongs to color j
        edges = p.new_variable(binary=True)

        if root is None:
            root = self.vertex_iterator().next()

        # r_edges is a relaxed variable grater than edges. It is used to
        # check the presence of cycles
        r_edges = p.new_variable(nonnegative=True)

        epsilon = 1/(3*(Integer(self.order())))

        if self.is_directed():
            # Does nothing to an edge.. Useful when out of "if self.directed"
            S = lambda x_y: x_y

            # An edge belongs to at most arborescence
            for e in self.edges(labels=False):
                p.add_constraint(p.sum([edges[j,e] for j in colors]), max=1)


            for j in colors:
                # each color class has self.order()-1 edges
                p.add_constraint(p.sum([edges[j,e] for e in self.edges(labels=None)]), min=self.order()-1)

                # Each vertex different from the root has indegree equals to one
                for v in self.vertices():
                    if v is not root:
                        p.add_constraint(p.sum([edges[j,e] for e in self.incoming_edges(v, labels=None)]), max=1, min=1)
                    else:
                        p.add_constraint(p.sum([edges[j,e] for e in self.incoming_edges(v, labels=None)]), max=0, min=0)

                # r_edges is larger than edges
                for u,v in self.edges(labels=None):
                    if self.has_edge(v,u):
                        if v<u:
                            p.add_constraint(r_edges[j,(u,v)] + r_edges[j,(v, u)] - edges[j,(u,v)] - edges[j,(v,u)], min=0)
                    else:
                        p.add_constraint(r_edges[j,(u,v)] + r_edges[j,(v, u)] - edges[j,(u,v)], min=0)

                from sage.graphs.digraph import DiGraph
                D = DiGraph()
                D.add_vertices(self.vertices())
                D.set_pos(self.get_pos())
                classes = [D.copy() for j in colors]

        else:

            # Sort an edge
            S = lambda x_y: x_y if x_y[0] < x_y[1] else (x_y[1], x_y[0])

            # An edge belongs to at most one arborescence
            for e in self.edges(labels=False):
                p.add_constraint(p.sum([edges[j,S(e)] for j in colors]), max=1)


            for j in colors:
                # each color class has self.order()-1 edges
                p.add_constraint(p.sum([edges[j,S(e)] for e in self.edges(labels=None)]), min=self.order()-1)

                # Each vertex is in the tree
                for v in self.vertices():
                    p.add_constraint(p.sum([edges[j,S(e)] for e in self.edges_incident(v, labels=None)]), min=1)

                # r_edges is larger than edges
                for u,v in self.edges(labels=None):
                    p.add_constraint(r_edges[j,(u,v)] + r_edges[j,(v, u)] - edges[j,S((u,v))], min=0)

                from sage.graphs.graph import Graph
                D = Graph()
                D.add_vertices(self.vertices())
                D.set_pos(self.get_pos())
                classes = [D.copy() for j in colors]

        # no cycles
        for j in colors:
            for v in self:
                p.add_constraint(p.sum(r_edges[j,(u,v)] for u in self.neighbors(v)), max=1-epsilon)
        try:
            p.solve(log = verbose)

        except MIPSolverException:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("This graph does not contain the required number of trees/arborescences !")

        edges = p.get_values(edges)

        for j,g in enumerate(classes):
            for e in self.edges(labels=False):
                if edges[j,S(e)] == 1:
                    g.add_edge(e)
            if len(list(g.breadth_first_search(root))) != self.order():
                raise RuntimeError("The computation seems to have gone wrong somewhere..."+
                                   "This is probably because of the value of epsilon, but"+
                                   " in any case please report this bug, with the graph "+
                                   "that produced it ! ;-)")

        return classes

    def edge_cut(self, s, t, value_only=True, use_edge_labels=False, vertices=False, method="FF", solver=None, verbose=0):
        r"""
        Returns a minimum edge cut between vertices `s` and `t`
        represented by a list of edges.

        A minimum edge cut between two vertices `s` and `t` of self
        is a set `A` of edges of minimum weight such that the graph
        obtained by removing `A` from self is disconnected. For more
        information, see the
        `Wikipedia article on cuts
        <http://en.wikipedia.org/wiki/Cut_%28graph_theory%29>`_.

        INPUT:

        - ``s`` -- source vertex

        - ``t`` -- sink vertex

        - ``value_only`` -- boolean (default: ``True``). When set to
          ``True``, only the weight of a minimum cut is returned.
          Otherwise, a list of edges of a minimum cut is also returned.

        - ``use_edge_labels`` -- boolean (default: ``False``). When set to
          ``True``, computes a weighted minimum cut where each edge has
          a weight defined by its label (if an edge has no label, `1`
          is assumed). Otherwise, each edge has weight `1`.

        - ``vertices`` -- boolean (default: ``False``). When set to
          ``True``, returns a list of edges in the edge cut and the
          two sets of vertices that are disconnected by the cut.

          Note: ``vertices=True`` implies ``value_only=False``.

        - ``method`` -- There are currently two different
          implementations of this method :

              * If ``method = "FF"`` (default), a Python
                implementation of the Ford-Fulkerson algorithm is
                used.

              * If ``method = "LP"``, the flow problem is solved using
                Linear Programming.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        .. NOTE::

           The use of Linear Programming for non-integer problems may
           possibly mean the presence of a (slight) numerical noise.

        OUTPUT:

        Real number or tuple, depending on the given arguments
        (examples are given below).

        EXAMPLES:

        A basic application in the Pappus graph::

           sage: g = graphs.PappusGraph()
           sage: g.edge_cut(1, 2, value_only=True)
           3

        Or on Petersen's graph, with the corresponding bipartition of
        the vertex set::

           sage: g = graphs.PetersenGraph()
           sage: g.edge_cut(0, 3, vertices=True)
           [3, [(0, 1, None), (0, 4, None), (0, 5, None)], [[0], [1, 2, 3, 4, 5, 6, 7, 8, 9]]]

        If the graph is a path with randomly weighted edges::

           sage: g = graphs.PathGraph(15)
           sage: for (u,v) in g.edge_iterator(labels=None):
           ...      g.set_edge_label(u,v,random())

        The edge cut between the two ends is the edge of minimum weight::

           sage: minimum = min([l for u,v,l in g.edge_iterator()])
           sage: minimum == g.edge_cut(0, 14, use_edge_labels=True)
           True
           sage: [value,[e]] = g.edge_cut(0, 14, use_edge_labels=True, value_only=False)
           sage: g.edge_label(e[0],e[1]) == minimum
           True

        The two sides of the edge cut are obviously shorter paths::

           sage: value,edges,[set1,set2] = g.edge_cut(0, 14, use_edge_labels=True, vertices=True)
           sage: g.subgraph(set1).is_isomorphic(graphs.PathGraph(len(set1)))
           True
           sage: g.subgraph(set2).is_isomorphic(graphs.PathGraph(len(set2)))
           True
           sage: len(set1) + len(set2) == g.order()
           True

        TESTS:

        If method is set to an exotic value::

           sage: g = graphs.PetersenGraph()
           sage: g.edge_cut(0,1, method="Divination")
           Traceback (most recent call last):
           ...
           ValueError: The method argument has to be equal to either "FF" or "LP"

        Same result for both methods::

           sage: g = graphs.RandomGNP(20,.3)
           sage: for u,v in g.edges(labels=False):
           ...      g.set_edge_label(u,v,round(random(),5))
           sage: g.edge_cut(0,1, method="FF") == g.edge_cut(0,1,method="LP")
           True

        Rounded return value when using the LP method::

           sage: g = graphs.PappusGraph()
           sage: g.edge_cut(1, 2, value_only=True, method = "LP")
           3
        """
        self._scream_if_not_simple(allow_loops=True)
        if vertices:
            value_only = False

        if use_edge_labels:
            weight = lambda x: x if (x!={} and x is not None) else 1
        else:
            weight = lambda x: 1

        if method == "FF":
            if value_only:
                return self.flow(s,t,value_only=value_only,use_edge_labels=use_edge_labels, method=method)

            flow_value, flow_graph = self.flow(s,t,value_only=value_only,use_edge_labels=use_edge_labels, method=method)
            g = self.copy()
            for u,v,l in flow_graph.edge_iterator():
                if (not use_edge_labels or
                    (weight(g.edge_label(u,v)) == weight(l))):
                    g.delete_edge(u,v)

            return_value = [flow_value]

            reachable_from_s = list(g.breadth_first_search(s))

            return_value.append(self.edge_boundary(reachable_from_s))

            if vertices:
                return_value.append([reachable_from_s,list(set(self.vertices())-set(reachable_from_s))])

            return return_value

        if method != "LP":
            raise ValueError("The method argument has to be equal to either \"FF\" or \"LP\"")

        from sage.numerical.mip import MixedIntegerLinearProgram
        g = self
        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        b = p.new_variable(binary = True)
        v = p.new_variable(binary = True)

        # Some vertices belong to part 1, others to part 0
        p.add_constraint(v[s], min=0, max=0)
        p.add_constraint(v[t], min=1, max=1)

        if g.is_directed():

            # we minimize the number of edges
            p.set_objective(p.sum([weight(w) * b[x,y] for (x,y,w) in g.edges()]))

            # Adjacent vertices can belong to different parts only if the
            # edge that connects them is part of the cut
            for (x,y) in g.edges(labels=None):
                p.add_constraint(v[x] + b[x,y] - v[y], min=0)

        else:
            # we minimize the number of edges
            p.set_objective(p.sum([weight(w) * b[min(x,y),max(x,y)] for (x,y,w) in g.edges()]))
            # Adjacent vertices can belong to different parts only if the
            # edge that connects them is part of the cut
            for (x,y) in g.edges(labels=None):
                p.add_constraint(v[x] + b[min(x,y),max(x,y)] - v[y], min=0)
                p.add_constraint(v[y] + b[min(x,y),max(x,y)] - v[x], min=0)

        if value_only:
            if use_edge_labels:
                return p.solve(objective_only=True, log=verbose)
            else:
                return Integer(round(p.solve(objective_only=True, log=verbose)))
        else:
            obj = p.solve(log=verbose)

            if use_edge_labels is False:
                obj = Integer(round(obj))

            b = p.get_values(b)
            answer = [obj]
            if g.is_directed():
                answer.append([(x,y) for (x,y) in g.edges(labels=None) if b[x,y] == 1])
            else:
                answer.append([(x,y) for (x,y) in g.edges(labels=None) if b[min(x,y),max(x,y)] == 1])

            if vertices:
                v = p.get_values(v)
                l0 = []
                l1 = []
                for x in g.vertex_iterator():
                    if x in v and v[x] == 1:
                        l1.append(x)
                    else:
                        l0.append(x)
                answer.append([l0, l1])
            return tuple(answer)

    def vertex_cut(self, s, t, value_only=True, vertices=False, solver=None, verbose=0):
        r"""
        Returns a minimum vertex cut between non-adjacent vertices `s` and `t`
        represented by a list of vertices.

        A vertex cut between two non-adjacent vertices is a set `U`
        of vertices of self such that the graph obtained by removing
        `U` from self is disconnected. For more information, see the
        `Wikipedia article on cuts
        <http://en.wikipedia.org/wiki/Cut_%28graph_theory%29>`_.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``). When set to
          ``True``, only the size of the minimum cut is returned.

        - ``vertices`` -- boolean (default: ``False``). When set to
          ``True``, also returns the two sets of vertices that
          are disconnected by the cut. Implies ``value_only``
          set to False.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        Real number or tuple, depending on the given arguments
        (examples are given below).

        EXAMPLE:

        A basic application in the Pappus graph::

           sage: g = graphs.PappusGraph()
           sage: g.vertex_cut(1, 16, value_only=True)
           3

        In the bipartite complete graph `K_{2,8}`, a cut between the two
        vertices in the size `2` part consists of the other `8` vertices::

           sage: g = graphs.CompleteBipartiteGraph(2, 8)
           sage: [value, vertices] = g.vertex_cut(0, 1, value_only=False)
           sage: print value
           8
           sage: vertices == range(2,10)
           True

        Clearly, in this case the two sides of the cut are singletons ::

           sage: [value, vertices, [set1, set2]] = g.vertex_cut(0,1, vertices=True)
           sage: len(set1) == 1
           True
           sage: len(set2) == 1
           True
        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        g = self
        if g.has_edge(s,t):
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("There can be no vertex cut between adjacent vertices !")
        if vertices:
            value_only = False

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        b = p.new_variable(binary=True)
        v = p.new_variable(binary=True)

        # Some vertices belong to part 1, some others to part 0
        p.add_constraint(v[s] == 0)
        p.add_constraint(v[t] == 1)

        # b indicates whether the vertices belong to the cut
        p.add_constraint(b[s] ==0)
        p.add_constraint(b[t] ==0)

        if g.is_directed():

            p.set_objective(p.sum([b[x] for x in g.vertices()]))

            # adjacent vertices belong to the same part except if one of them
            # belongs to the cut
            for (x,y) in g.edges(labels=None):
                p.add_constraint(v[x] + b[y] - v[y], min=0)

        else:
            p.set_objective(p.sum([b[x] for x in g.vertices()]))
            # adjacent vertices belong to the same part except if one of them
            # belongs to the cut
            for (x,y) in g.edges(labels=None):
                p.add_constraint(v[x] + b[y] >= v[y])
                p.add_constraint(v[y] + b[x] >= v[x])

        if value_only:
            return Integer(round(p.solve(objective_only=True, log=verbose)))
        else:
            obj = Integer(round(p.solve(log=verbose)))
            b = p.get_values(b)
            answer = [obj,[x for x in g if b[x] == 1]]
            if vertices:
                v = p.get_values(v)
                l0 = []
                l1 = []
                for x in g.vertex_iterator():
                    # if the vertex is not in the cut
                    if not (x in b and b[x] == 1):
                        if (x in v and v[x] == 1):
                            l1.append(x)
                        else:
                            l0.append(x)
                answer.append([l0, l1])
            return tuple(answer)


    def multiway_cut(self, vertices, value_only = False, use_edge_labels = False, solver = None, verbose = 0):
        r"""
        Returns a minimum edge multiway cut corresponding to the
        given set of vertices
        ( cf. http://www.d.kth.se/~viggo/wwwcompendium/node92.html )
        represented by a list of edges.

        A multiway cut for a vertex set `S` in a graph or a digraph
        `G` is a set `C` of edges such that any two vertices `u,v`
        in `S` are disconnected when removing the edges from `C` from `G`.

        Such a cut is said to be minimum when its cardinality
        (or weight) is minimum.

        INPUT:

        - ``vertices`` (iterable)-- the set of vertices

        - ``value_only`` (boolean)

            - When set to ``True``, only the value of a minimum
              multiway cut is returned.

            - When set to ``False`` (default), the list of edges
              is returned

        - ``use_edge_labels`` (boolean)
            - When set to ``True``, computes a weighted minimum cut
              where each edge has a weight defined by its label. ( if
              an edge has no label, `1` is assumed )

            - when set to ``False`` (default), each edge has weight `1`.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        Of course, a multiway cut between two vertices correspond
        to a minimum edge cut ::

            sage: g = graphs.PetersenGraph()
            sage: g.edge_cut(0,3) == g.multiway_cut([0,3], value_only = True)
            True

        As Petersen's graph is `3`-regular, a minimum multiway cut
        between three vertices contains at most `2\times 3` edges
        (which could correspond to the neighborhood of 2
        vertices)::

            sage: g.multiway_cut([0,3,9], value_only = True) == 2*3
            True

        In this case, though, the vertices are an independent set.
        If we pick instead vertices `0,9,` and `7`, we can save `4`
        edges in the multiway cut ::

            sage: g.multiway_cut([0,7,9], value_only = True) == 2*3 - 1
            True

        This example, though, does not work in the directed case anymore,
        as it is not possible in Petersen's graph to mutualise edges ::

            sage: g = DiGraph(g)
            sage: g.multiway_cut([0,7,9], value_only = True) == 3*3
            True

        Of course, a multiway cut between the whole vertex set
        contains all the edges of the graph::

            sage: C = g.multiway_cut(g.vertices())
            sage: set(C) == set(g.edges())
            True
        """
        self._scream_if_not_simple(allow_loops=True)
        from sage.numerical.mip import MixedIntegerLinearProgram
        from itertools import combinations, chain

        p = MixedIntegerLinearProgram(maximization = False, solver= solver)

        # height[c,v] represents the height of vertex v for commodity c
        height = p.new_variable(nonnegative=True)

        # cut[e] represents whether e is in the cut
        cut = p.new_variable(binary = True)

        # Reorder
        R = lambda x,y : (x,y) if x<y else (y,x)

        # Weight function
        if use_edge_labels:
            w = lambda l : l if l is not None else 1
        else:
            w = lambda l : 1

        if self.is_directed():

            p.set_objective( p.sum([ w(l) * cut[u,v] for u,v,l in self.edge_iterator() ]) )

            for s,t in chain( combinations(vertices,2), map(lambda x_y: (x_y[1],x_y[0]), combinations(vertices,2))) :
                # For each commodity, the source is at height 0
                # and the destination is at height 1
                p.add_constraint( height[(s,t),s], min = 0, max = 0)
                p.add_constraint( height[(s,t),t], min = 1, max = 1)

                # given a commodity (s,t), the height of two adjacent vertices u,v
                # can differ of at most the value of the edge between them
                for u,v in self.edges(labels = False):
                    p.add_constraint( height[(s,t),u] - height[(s,t),v] - cut[u,v], max = 0)

        else:

            p.set_objective( p.sum([ w(l) * cut[R(u,v)] for u,v,l in self.edge_iterator() ]) )

            for s,t in combinations(vertices,2):
                # For each commodity, the source is at height 0
                # and the destination is at height 1
                p.add_constraint( height[(s,t),s], min = 0, max = 0)
                p.add_constraint( height[(s,t),t], min = 1, max = 1)

                # given a commodity (s,t), the height of two adjacent vertices u,v
                # can differ of at most the value of the edge between them
                for u,v in self.edges(labels = False):
                    p.add_constraint( height[(s,t),u] - height[(s,t),v] - cut[R(u,v)], max = 0)
                    p.add_constraint( height[(s,t),v] - height[(s,t),u] - cut[R(u,v)], max = 0)

        if value_only:
            if use_edge_labels:
                return p.solve(objective_only = True, log = verbose)
            else:
                return Integer(round(p.solve(objective_only = True, log = verbose)))

        p.solve(log = verbose)

        cut = p.get_values(cut)

        if self.is_directed():
            return [x for x in self.edge_iterator() if cut[x[0], x[1]] == 1]

        return [x for x in self.edge_iterator() if cut[R(x[0], x[1])] == 1]


    def max_cut(self, value_only=True, use_edge_labels=False, vertices=False, solver=None, verbose=0):
        r"""
        Returns a maximum edge cut of the graph. For more information, see the
        `Wikipedia article on cuts
        <http://en.wikipedia.org/wiki/Cut_%28graph_theory%29>`_.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``)

          - When set to ``True`` (default), only the value is returned.

          - When set to ``False``, both the value and a maximum edge cut
            are returned.

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - When set to ``True``, computes a maximum weighted cut
            where each edge has a weight defined by its label. (If
            an edge has no label, `1` is assumed.)

          - When set to ``False``, each edge has weight `1`.

        - ``vertices`` -- boolean (default: ``False``)

          - When set to ``True``, also returns the two sets of
            vertices that are disconnected by the cut. This implies
            ``value_only=False``.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLE:

        Quite obviously, the max cut of a bipartite graph
        is the number of edges, and the two sets of vertices
        are the the two sides ::

            sage: g = graphs.CompleteBipartiteGraph(5,6)
            sage: [ value, edges, [ setA, setB ]] = g.max_cut(vertices=True)
            sage: value == 5*6
            True
            sage: bsetA, bsetB  = map(list,g.bipartite_sets())
            sage: (bsetA == setA and bsetB == setB ) or ((bsetA == setB and bsetB == setA ))
            True

        The max cut of a Petersen graph::

           sage: g=graphs.PetersenGraph()
           sage: g.max_cut()
           12

        """
        self._scream_if_not_simple(allow_loops=True)
        g=self

        if vertices:
            value_only=False

        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            weight = lambda x: x if x in RR else 1
        else:
            weight = lambda x: 1

        if g.is_directed():
            reorder_edge = lambda x,y : (x,y)
        else:
            reorder_edge = lambda x,y : (x,y) if x<= y else (y,x)

        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=True, solver=solver)

        in_set = p.new_variable(binary = True)
        in_cut = p.new_variable(binary = True)

        # A vertex has to be in some set
        for v in g:
            p.add_constraint(in_set[0,v]+in_set[1,v],max=1,min=1)

        # There is no empty set
        p.add_constraint(p.sum([in_set[1,v] for v in g]),min=1)
        p.add_constraint(p.sum([in_set[0,v] for v in g]),min=1)

        if g.is_directed():
            # There is no edge from set 0 to set 1 which
            # is not in the cut
            # Besides, an edge can only be in the cut if its vertices
            # belong to different sets
            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u] + in_set[1,v] - in_cut[(u,v)], max = 1)
                p.add_constraint(in_set[0,u] + in_set[0,v] + in_cut[(u,v)], max = 2)
                p.add_constraint(in_set[1,u] + in_set[1,v] + in_cut[(u,v)], max = 2)
        else:

            # Two adjacent vertices are in different sets if and only if
            # the edge between them is in the cut

            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u] + in_set[1,v] - in_cut[reorder_edge(u,v)], max=1)
                p.add_constraint(in_set[1,u] + in_set[0,v] - in_cut[reorder_edge(u,v)], max=1)
                p.add_constraint(in_set[0,u] + in_set[0,v] + in_cut[reorder_edge(u,v)], max = 2)
                p.add_constraint(in_set[1,u] + in_set[1,v] + in_cut[reorder_edge(u,v)], max = 2)

        p.set_objective(p.sum([weight(l ) * in_cut[reorder_edge(u,v)] for (u,v,l ) in g.edge_iterator()]))

        if value_only:
            obj = p.solve(objective_only=True, log=verbose)
            return obj if use_edge_labels else Integer(round(obj))
        else:
            obj = p.solve(log=verbose)

            if use_edge_labels:
                obj = Integer(round(obj))

            val = [obj]

            in_cut = p.get_values(in_cut)
            in_set = p.get_values(in_set)

            edges = []
            for (u,v,l) in g.edge_iterator():
                if in_cut[reorder_edge(u,v)] == 1:
                    edges.append((u,v,l))

            val.append(edges)

            if vertices:
                a = []
                b = []
                for v in g:
                    if in_set[0,v] == 1:
                        a.append(v)
                    else:
                        b.append(v)
                val.append([a,b])

            return val

    def longest_path(self, s=None, t=None, use_edge_labels=False, algorithm="MILP", solver=None, verbose=0):
        r"""
        Returns a longest path of ``self``.

        INPUT:

        - ``s`` (vertex) -- forces the source of the path (the method then
          returns the longest path starting at ``s``). The argument is set to
          ``None`` by default, which means that no constraint is set upon the
          first vertex in the path.

        - ``t`` (vertex) -- forces the destination of the path (the method then
          returns the longest path ending at ``t``). The argument is set to
          ``None`` by default, which means that no constraint is set upon the
          last vertex in the path.

        - ``use_edge_labels`` (boolean) -- whether the labels on the edges are
          to be considered as weights (a label set to ``None`` or ``{}`` being
          considered as a weight of `1`). Set to ``False`` by default.

        - ``algorithm`` -- one of ``"MILP"`` (default) or ``"backtrack"``. Two
          remarks on this respect:

              * While the MILP formulation returns an exact answer, the
                backtrack algorithm is a randomized heuristic.

              * As the backtrack algorithm does not support edge weighting,
                setting ``use_edge_labels=True`` will force the use of the MILP
                algorithm.

        - ``solver`` -- (default: ``None``) Specify the Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        .. NOTE::

            The length of a path is assumed to be the number of its edges, or
            the sum of their labels.

        OUTPUT:

        A subgraph of ``self`` corresponding to a (directed if ``self`` is
        directed) longest path. If ``use_edge_labels == True``, a pair ``weight,
        path`` is returned.

        ALGORITHM:

        Mixed Integer Linear Programming. (This problem is known to be NP-Hard).

        EXAMPLES:

        Petersen's graph being hypohamiltonian, it has a longest path
        of length `n-2`::

            sage: g = graphs.PetersenGraph()
            sage: lp = g.longest_path()
            sage: lp.order() >= g.order() - 2
            True

        The heuristic totally agrees::

            sage: g = graphs.PetersenGraph()
            sage: g.longest_path(algorithm="backtrack").edges()
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 4, None), (4, 9, None), (5, 7, None), (5, 8, None), (6, 8, None), (6, 9, None)]

        Let us compute longest paths on random graphs with random weights. Each
        time, we ensure the resulting graph is indeed a path::

            sage: for i in range(20):
            ...       g = graphs.RandomGNP(15, 0.3)
            ...       for u, v in g.edges(labels=False):
            ...           g.set_edge_label(u, v, random())
            ...       lp = g.longest_path()
            ...       if (not lp.is_forest() or
            ...           not max(lp.degree()) <= 2 or
            ...           not lp.is_connected()):
            ...           print("Error!")
            ...           break

        TESTS:

        The argument ``algorithm`` must be either ``'backtrack'`` or
        ``'MILP'``::

            sage: graphs.PetersenGraph().longest_path(algorithm="abc")
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be either 'backtrack' or 'MILP'

        Disconnected graphs not weighted::

            sage: g1 = graphs.PetersenGraph()
            sage: g2 = 2 * g1
            sage: lp1 = g1.longest_path()
            sage: lp2 = g2.longest_path()
            sage: len(lp1) == len(lp2)
            True

        Disconnected graphs weighted::

            sage: g1 = graphs.PetersenGraph()
            sage: for u,v in g.edges(labels=False):
            ...       g.set_edge_label(u, v, random())
            sage: g2 = 2 * g1
            sage: lp1 = g1.longest_path(use_edge_labels=True)
            sage: lp2 = g2.longest_path(use_edge_labels=True)
            sage: lp1[0] == lp2[0]
            True

        Empty graphs::

            sage: Graph().longest_path()
            Graph on 0 vertices
            sage: Graph().longest_path(use_edge_labels=True)
            [0, Graph on 0 vertices]
            sage: graphs.EmptyGraph().longest_path()
            Graph on 0 vertices
            sage: graphs.EmptyGraph().longest_path(use_edge_labels=True)
            [0, Graph on 0 vertices]

        Trivial graphs::

            sage: G = Graph()
            sage: G.add_vertex(0)
            sage: G.longest_path()
            Graph on 0 vertices
            sage: G.longest_path(use_edge_labels=True)
            [0, Graph on 0 vertices]
            sage: graphs.CompleteGraph(1).longest_path()
            Graph on 0 vertices
            sage: graphs.CompleteGraph(1).longest_path(use_edge_labels=True)
            [0, Graph on 0 vertices]

        Random test for digraphs::

            sage: for i in range(20):
            ...       g = digraphs.RandomDirectedGNP(15, 0.3)
            ...       for u, v in g.edges(labels=False):
            ...           g.set_edge_label(u, v, random())
            ...       lp = g.longest_path()
            ...       if (not lp.is_directed_acyclic() or
            ...           not max(lp.out_degree()) <= 1 or
            ...           not max(lp.in_degree()) <= 1 or
            ...           not lp.is_connected()):
            ...           print("Error!")
            ...           print g.edges()
            ...           break

        :trac:`13019`::

            sage: g = graphs.CompleteGraph(5).to_directed()
            sage: g.longest_path(s=1,t=2)
            Subgraph of (Complete graph): Digraph on 5 vertices

        :trac:`14412`::

            sage: l = [(0, 1), (0, 3), (2, 0)]
            sage: G = DiGraph(l)
            sage: G.longest_path().edges()
            [(0, 1, None), (2, 0, None)]
        """
        self._scream_if_not_simple()

        if use_edge_labels:
            algorithm = "MILP"
        if algorithm not in ("backtrack", "MILP"):
            raise ValueError("algorithm must be either 'backtrack' or 'MILP'")

        # Quick improvement
        if not self.is_connected():
            if use_edge_labels:
                return max(g.longest_path(s=s, t=t,
                                          use_edge_labels=use_edge_labels,
                                          algorithm=algorithm)
                           for g in self.connected_components_subgraphs())
            else:
                return max((g.longest_path(s=s, t=t,
                                           use_edge_labels=use_edge_labels,
                                           algorithm=algorithm)
                            for g in self.connected_components_subgraphs()),
                           key=lambda x: x.order())

        # Stupid cases
        # - Graph having <= 1 vertex.
        #
        # - The source has outdegree 0 in a directed graph, or
        #   degree 0, or is not a vertex of the graph.
        #
        # - The destination has indegree 0 in a directed graph, or
        #   degree 0, or is not a vertex of the graph.
        #
        # - Both s and t are specified, but there is no path between
        #   the two in a directed graph (the graph is connected).
        if (self.order() <= 1 or
            (s is not None and (
                    (s not in self) or
                    (self._directed and self.out_degree(s) == 0) or
                    (not self._directed and self.degree(s) == 0))) or
            (t is not None and (
                    (t not in self) or
                    (self._directed and self.in_degree(t) == 0) or
                    (not self._directed and self.degree(t) == 0))) or
            (self._directed and (s is not None) and (t is not None) and
             len(self.shortest_path(s, t)) == 0)):
            if self._directed:
                from sage.graphs.all import DiGraph
                return [0, DiGraph()] if use_edge_labels else DiGraph()
            from sage.graphs.all import Graph
            return [0, Graph()] if use_edge_labels else Graph()

        # Calling the backtrack heuristic if asked
        if algorithm == "backtrack":
            from sage.graphs.generic_graph_pyx import find_hamiltonian as fh
            x = fh(self, find_path=True)[1]
            return self.subgraph(vertices=x, edges=zip(x[:-1], x[1:]))

        ##################
        # LP Formulation #
        ##################

        # Epsilon... Must be less than 1/(n+1), but we want to avoid
        # numerical problems...
        epsilon = 1/(2*float(self.order()))

        # Associating a weight to a label
        if use_edge_labels:
            weight = lambda x: x if (x is not None and x != {}) else 1
        else:
            weight = lambda x: 1

        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(solver=solver)

        # edge_used[(u,v)] == 1 if (u,v) is used
        edge_used = p.new_variable(binary=True)

        # relaxed version of the previous variable, to prevent the
        # creation of cycles
        r_edge_used = p.new_variable(nonnegative=True)

        # vertex_used[v] == 1 if vertex v is used
        vertex_used = p.new_variable(binary=True)

        if self._directed:

            # if edge uv is used, vu can not be
            for u, v in self.edges(labels=False):
                if self.has_edge(v, u):
                    p.add_constraint(edge_used[(u,v)] + edge_used[(v,u)] <= 1)

            # A vertex is used if one of its incident edges is
            for u,v in self.edges(labels = False):
                p.add_constraint(vertex_used[v] >= edge_used[(u,v)])
                p.add_constraint(vertex_used[u] >= edge_used[(u,v)])

            # A path is a tree. If n vertices are used, at most n-1 edges are
            p.add_constraint(
                  p.sum(vertex_used[v] for v in self)
                - p.sum(edge_used[e] for e in self.edges(labels=False))
                  == 1)

            # A vertex has at most one incoming used edge and at most
            # one outgoing used edge
            for v in self:
                p.add_constraint(
                    p.sum(edge_used[(u,v)] for u in self.neighbors_in(v)) <= 1)
                p.add_constraint(
                    p.sum(edge_used[(v,u)] for u in self.neighbors_out(v)) <= 1)

            # r_edge_used is "more" than edge_used, though it ignores
            # the direction
            for u, v in self.edges(labels=False):
                p.add_constraint(r_edge_used[(u,v)] + r_edge_used[(v,u)]
                                 >= edge_used[(u,v)])

            # No cycles
            for v in self:
                p.add_constraint(
                    p.sum(r_edge_used[(u,v)] for u in self.neighbors(v))
                    <= 1-epsilon)

            # Enforcing the source if asked.. If s is set, it has no
            # incoming edge and exactly one son
            if s is not None:
                p.add_constraint(
                    p.sum(edge_used[(u,s)] for u in self.neighbors_in(s)),
                    max=0, min=0)
                p.add_constraint(
                    p.sum(edge_used[(s,u)] for u in self.neighbors_out(s)),
                    min=1, max=1)

            # Enforcing the destination if asked.. If t is set, it has
            # no outgoing edge and exactly one parent
            if t is not None:
                p.add_constraint(
                    p.sum(edge_used[(u,t)] for u in self.neighbors_in(t)),
                    min=1, max=1)
                p.add_constraint(
                    p.sum(edge_used[(t,u)] for u in self.neighbors_out(t)),
                    max=0, min=0)

            # Defining the objective
            p.set_objective(
                p.sum(weight(l) * edge_used[(u,v)] for u, v, l in self.edges()))
        else:
            # f_edge_used calls edge_used through reordering u and v
            # to avoid having two different variables
            f_edge_used = lambda u, v: edge_used[
                (u,v) if hash(u) < hash(v) else (v,u)]
            # A vertex is used if one of its incident edges is
            for v in self:
                for u in self.neighbors(v):
                    p.add_constraint(vertex_used[v] - f_edge_used(u,v), min=0)
            # A path is a tree. If n vertices are used, at most n-1 edges are
            p.add_constraint(
                p.sum(vertex_used[v] for v in self)
                - p.sum(f_edge_used(u,v) for u, v in self.edges(labels=False)),
                min=1, max=1)
            # A vertex has at most two incident edges used
            for v in self:
                p.add_constraint(
                    p.sum(f_edge_used(u,v) for u in self.neighbors(v)), max=2)
            # r_edge_used is "more" than edge_used
            for u, v in self.edges(labels=False):
                p.add_constraint(r_edge_used[(u,v)]
                                 + r_edge_used[(v,u)]
                                 - f_edge_used(u,v),
                                 min=0)
            # No cycles
            for v in self:
                p.add_constraint(
                    p.sum(r_edge_used[(u,v)] for u in self.neighbors(v)),
                    max=1-epsilon)
            # Enforcing the destination if asked.. If s or t are set,
            # they have exactly one incident edge
            if s is not None:
                p.add_constraint(
                    p.sum(f_edge_used(s,u) for u in self.neighbors(s)),
                    max=1, min=1)
            if t is not None:
                p.add_constraint(
                    p.sum(f_edge_used(t,u) for u in self.neighbors(t)),
                    max=1, min=1)
            # Defining the objective
            p.set_objective(p.sum(weight(l) * f_edge_used(u,v)
                                for u, v, l in self.edges()))
        # Computing the result. No exception has to be raised, as this
        # problem always has a solution (there is at least one edge,
        # and a path from s to t if they are specified).
        p.solve(log=verbose)
        edge_used = p.get_values(edge_used)
        vertex_used = p.get_values(vertex_used)
        if self._directed:
            g = self.subgraph(
                vertices=(v for v in self if vertex_used[v] == 1),
                edges=((u,v,l) for u, v, l in self.edges()
                       if edge_used[(u,v)] == 1))
        else:
            g = self.subgraph(
                vertices=(v for v in self if vertex_used[v] == 1),
                edges=((u,v,l) for u, v, l in self.edges()
                       if f_edge_used(u,v) == 1))
        if use_edge_labels:
            return sum(map(weight, g.edge_labels())), g
        else:
            return g

    def traveling_salesman_problem(self, use_edge_labels = False, solver = None, constraint_generation = None, verbose = 0, verbose_constraints = False):
        r"""
        Solves the traveling salesman problem (TSP)

        Given a graph (resp. a digraph) `G` with weighted edges, the traveling
        salesman problem consists in finding a Hamiltonian cycle (resp. circuit)
        of the graph of minimum cost.

        This TSP is one of the most famous NP-Complete problems, this function
        can thus be expected to take some time before returning its result.

        INPUT:

        - ``use_edge_labels`` (boolean) -- whether to consider the weights of
          the edges.

              - If set to ``False`` (default), all edges are assumed to weight
                `1`

              - If set to ``True``, the weights are taken into account, and the
                circuit returned is the one minimizing the sum of the weights
                (an edge with no label is assumed to have weight `1`).

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``constraint_generation`` (boolean) -- whether to use constraint
          generation when solving the Mixed Integer Linear Program.

          When ``constraint_generation = None``, constraint generation is used
          whenever the graph has a density larger than 70%.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``verbose_constraints`` -- whether to display which constraints are
          being generated.

        OUTPUT:

        A solution to the TSP, as a ``Graph`` object whose vertex set is `V(G)`,
        and whose edges are only those of the solution.

        ALGORITHM:

        This optimization problem is solved through the use of Linear
        Programming.

        .. NOTE::

            This function is correctly defined for both graph and digraphs. In
            the second case, the returned cycle is a circuit of optimal cost.

        EXAMPLES:

        The Heawood graph is known to be Hamiltonian::

            sage: g = graphs.HeawoodGraph()
            sage: tsp = g.traveling_salesman_problem()
            sage: tsp
            TSP from Heawood graph: Graph on 14 vertices

        The solution to the TSP has to be connected ::

            sage: tsp.is_connected()
            True

        It must also be a `2`-regular graph::

            sage: tsp.is_regular(k=2)
            True

        And obviously it is a subgraph of the Heawood graph::

            sage: tsp.is_subgraph(g, induced=False)
            True

        On the other hand, the Petersen Graph is known not to
        be Hamiltonian::

            sage: g = graphs.PetersenGraph()
            sage: tsp = g.traveling_salesman_problem()
            Traceback (most recent call last):
            ...
            EmptySetError: The given graph is not hamiltonian

        One easy way to change it is obviously to add to this graph the edges
        corresponding to a Hamiltonian cycle. If we do this by setting the cost
        of these new edges to `2`, while the others are set to `1`, we notice
        that not all the edges we added are used in the optimal solution ::

            sage: for u, v in g.edges(labels = None):
            ...      g.set_edge_label(u,v,1)

            sage: cycle = graphs.CycleGraph(10)
            sage: for u,v in cycle.edges(labels = None):
            ...      if not g.has_edge(u,v):
            ...          g.add_edge(u,v)
            ...      g.set_edge_label(u,v,2)

            sage: tsp = g.traveling_salesman_problem(use_edge_labels = True)
            sage: sum( tsp.edge_labels() ) < 2*10
            True

        If we pick `1/2` instead of `2` as a cost for these new edges, they
        clearly become the optimal solution::

            sage: for u,v in cycle.edges(labels = None):
            ...      g.set_edge_label(u,v,1/2)

            sage: tsp = g.traveling_salesman_problem(use_edge_labels = True)
            sage: sum( tsp.edge_labels() ) == (1/2)*10
            True

        TESTS:

        Comparing the results returned according to the value of
        ``constraint_generation``. First, for graphs::

            sage: from operator import itemgetter
            sage: n = 20
            sage: for i in range(20):
            ...       g = Graph()
            ...       g.allow_multiple_edges(False)
            ...       for u,v in graphs.RandomGNP(n,.2).edges(labels = False):
            ...            g.add_edge(u,v,round(random(),5))
            ...       for u,v in graphs.CycleGraph(n).edges(labels = False):
            ...            if not g.has_edge(u,v):
            ...                g.add_edge(u,v,round(random(),5))
            ...       v1 = g.traveling_salesman_problem(constraint_generation = False, use_edge_labels = True)
            ...       v2 = g.traveling_salesman_problem(use_edge_labels = True)
            ...       c1 = sum(map(itemgetter(2), v1.edges()))
            ...       c2 = sum(map(itemgetter(2), v2.edges()))
            ...       if c1 != c2:
            ...           print "Error !",c1,c2
            ...           break

        Then for digraphs::

            sage: from operator import itemgetter
            sage: set_random_seed(0)
            sage: n = 20
            sage: for i in range(20):
            ...       g = DiGraph()
            ...       g.allow_multiple_edges(False)
            ...       for u,v in digraphs.RandomDirectedGNP(n,.2).edges(labels = False):
            ...            g.add_edge(u,v,round(random(),5))
            ...       for u,v in digraphs.Circuit(n).edges(labels = False):
            ...            if not g.has_edge(u,v):
            ...                g.add_edge(u,v,round(random(),5))
            ...       v2 = g.traveling_salesman_problem(use_edge_labels = True)
            ...       v1 = g.traveling_salesman_problem(constraint_generation = False, use_edge_labels = True)
            ...       c1 = sum(map(itemgetter(2), v1.edges()))
            ...       c2 = sum(map(itemgetter(2), v2.edges()))
            ...       if c1 != c2:
            ...           print "Error !",c1,c2
            ...           print "With constraint generation :",c2
            ...           print "Without constraint generation :",c1
            ...           break

        Simple tests for multiple edges and loops::

            sage: G = DiGraph(multiedges=True, loops=True)
            sage: G.is_hamiltonian()
            Traceback (most recent call last):
            ...
            ValueError: The traveling salesman problem (or finding Hamiltonian cycle) is not well defined for graphs with less than 2 vertices
            sage: G.add_vertex(0)
            sage: G.is_hamiltonian()
            Traceback (most recent call last):
            ...
            ValueError: The traveling salesman problem (or finding Hamiltonian cycle) is not well defined for graphs with less than 2 vertices
            sage: G.add_edge(0,0,1)
            sage: G.add_edge(0,0,2)
            sage: tsp = G.traveling_salesman_problem(use_edge_labels=True)
            Traceback (most recent call last):
            ...
            ValueError: The traveling salesman problem (or finding Hamiltonian cycle) is not well defined for graphs with less than 2 vertices
            sage: G.add_vertex(1)
            sage: G.is_hamiltonian()
            False
            sage: G.add_edge(0,1,2)
            sage: G.add_edge(0,1,3)
            sage: G.add_edge(1,1,1)
            sage: G.add_edge(1,0,2)
            sage: G.is_hamiltonian()
            True
            sage: tsp = G.traveling_salesman_problem(use_edge_labels=True)
            sage: sum(tsp.edge_labels())
            4

        Graphs on 2 vertices::

            sage: Graph([(0,1),(0,1)],multiedges=True).is_hamiltonian()
            True
            sage: DiGraph([(0,1),(0,1)],multiedges=True).is_hamiltonian()
            False
            sage: DiGraph([(0,1),(1,0)],multiedges=True).is_hamiltonian()
            True
            sage: G = DiGraph(loops=True)
            sage: G.add_edges([(0,0),(0,1),(1,1),(1,0)])  # i.e. complete digraph with loops
            sage: G.is_hamiltonian()
            True
            sage: G.remove_loops()
            sage: G.is_hamiltonian()
            True
            sage: G.allow_loops(False)
            sage: G.is_hamiltonian()
            True

        Check that weight 0 edges are handled correctly (see :trac:`16214`)::

            sage: G = Graph([(0,1,1),(0,2,0),(0,3,1),(1,2,1),(1,3,0),(2,3,1)])
            sage: tsp = G.traveling_salesman_problem(use_edge_labels=True)
            sage: sum(tsp.edge_labels())
            2
        """
        from sage.categories.sets_cat import EmptySetError

        weight = lambda l : 1 if l is None else l

        ########################
        # 0 or 1 vertex graphs #
        ########################

        if self.order() < 2:
            raise ValueError("The traveling salesman problem (or finding "+
                 "Hamiltonian cycle) is not well defined for graphs with "+
                 "less than 2 vertices")

        #####################
        # 2-vertices graphs #
        #####################

        if self.order() == 2:
            u,v = self.vertices()
            if self.is_directed():
                if self.has_edge(u,v) and self.has_edge(v,u):
                    if self.has_multiple_edges():
                        edges = [(u,v,min(self.edge_label(u,v), key=weight)),
                                 (v,u,min(self.edge_label(v,u), key=weight))]
                    else:
                        edges = [(u,v,self.edge_label(u,v)),
                                 (v,u,self.edge_label(v,u))]
                    answer = self.subgraph(edges = edges)
                    answer.set_pos(self.get_pos())
                    answer.name("TSP from "+self.name())
                    return answer
            else:
                if self.has_multiple_edges() and len(self.edge_label(u,v)) > 1:
                    edges = self.edges()
                    edges.sort(key=weight)
                    answer = self.subgraph(edges = edges[:2])
                    answer.set_pos(self.get_pos())
                    answer.name("TSP from "+self.name())
                    return answer

            raise EmptySetError("The given graph is not hamiltonian")

        ################################
        # Quick checks of connectivity #
        ################################

        # TODO : Improve it by checking vertex-connectivity instead of
        # edge-connectivity.... But calling the vertex_connectivity (which
        # builds a LP) is way too slow. These tests only run traversals written
        # in Cython --> hence FAST

        if self.is_directed():
            if not self.is_strongly_connected():
                raise EmptySetError("The given graph is not hamiltonian")

        else:
            # Checks whether the graph is 2-connected
            if not self.strong_orientation().is_strongly_connected():
                raise EmptySetError("The given graph is not hamiltonian")

        ############################
        # Deal with multiple edges #
        ############################

        if self.has_loops() or self.has_multiple_edges():
            g = self.copy(immutable=False)
        else:
            g = self

        if g.has_multiple_edges():
            multi = g.multiple_edges()
            g.delete_edges(multi)
            g.allow_multiple_edges(False)
            if use_edge_labels:
                e = {}

                # The weight of an edge is the minimum over the weights of the parallel edges
                for u,v,l in multi:
                    if (u,v) in e:
                        e[u,v] = weight(l) if weight(l) < e[u,v] else e[u,v]
                    else:
                        e[u,v] = l

                g.add_edges([(u,v) for (u,v),l in e.iteritems()])

            else:
                g.add_edges(multi)

        if g.has_loops():
            g.remove_loops()
            g.allow_loops(False)

        if constraint_generation is None:
            if g.density() > .7:
                constraint_generation = False
            else:
                constraint_generation = True

        from sage.numerical.mip import MixedIntegerLinearProgram
        from sage.numerical.mip import MIPSolverException

        ####################################################
        # Constraint-generation formulation of the problem #
        ####################################################

        if constraint_generation:

            p = MixedIntegerLinearProgram(maximization = False,
                                          solver = solver,
                                          constraint_generation = True)


            # Directed Case #
            #################
            if g.is_directed():

                from sage.graphs.all import DiGraph
                b = p.new_variable(binary = True)

                # Objective function
                if use_edge_labels:
                    p.set_objective(p.sum([ weight(l)*b[u,v]
                                          for u,v,l in g.edges()]))

                # All the vertices have in-degree 1 and out-degree 1
                for v in g:
                    p.add_constraint(p.sum([b[u,v] for u in g.neighbors_in(v)]),
                                     min = 1,
                                     max = 1)
                    p.add_constraint(p.sum([b[v,u] for u in g.neighbors_out(v)]),
                                     min = 1,
                                     max = 1)

                # Initial Solve
                try:
                    p.solve(log = verbose)
                except MIPSolverException:
                    raise EmptySetError("The given graph is not hamiltonian")

                while True:
                    # We build the DiGraph representing the current solution
                    h = DiGraph()
                    for u,v,l in g.edges():
                        if p.get_values(b[u,v]) == 1:
                            h.add_edge(u,v,l)

                    # If there is only one circuit, we are done !
                    cc = h.connected_components()
                    if len(cc) == 1:
                        break

                    # Adding the corresponding constraint
                    if verbose_constraints:
                        print "Adding a constraint on set",cc[0]


                    p.add_constraint(p.sum(b[u,v] for u,v in
                                         g.edge_boundary(cc[0], labels = False)),
                                     min = 1)

                    try:
                        p.solve(log = verbose)
                    except MIPSolverException:
                        raise EmptySetError("The given graph is not hamiltonian")

            # Undirected Case #
            ###################
            else:

                from sage.graphs.all import Graph
                b = p.new_variable(binary = True)
                B = lambda u,v : b[frozenset((u,v))]

                # Objective function
                if use_edge_labels:
                    p.set_objective(p.sum([ weight(l)*B(u,v)
                                          for u,v,l in g.edges()]) )

                # All the vertices have degree 2
                for v in g:
                    p.add_constraint(p.sum([ B(u,v) for u in g.neighbors(v)]),
                                     min = 2,
                                     max = 2)

                # Initial Solve
                try:
                    p.solve(log = verbose)
                except MIPSolverException:
                    raise EmptySetError("The given graph is not hamiltonian")

                while True:
                    # We build the DiGraph representing the current solution
                    h = Graph()
                    for u,v,l in g.edges():
                        if p.get_values(B(u,v)) == 1:
                            h.add_edge(u,v,l)

                    # If there is only one circuit, we are done !
                    cc = h.connected_components()
                    if len(cc) == 1:
                        break

                    # Adding the corresponding constraint
                    if verbose_constraints:
                        print "Adding a constraint on set",cc[0]

                    p.add_constraint(p.sum(B(u,v) for u,v in
                                         g.edge_boundary(cc[0], labels = False)),
                                     min = 2)

                    try:
                        p.solve(log = verbose)
                    except MIPSolverException:
                        raise EmptySetError("The given graph is not hamiltonian")

            # We can now return the TSP !
            answer = self.subgraph(edges = h.edges())
            answer.set_pos(self.get_pos())
            answer.name("TSP from "+g.name())
            return answer

        #################################################
        # ILP formulation without constraint generation #
        #################################################

        p = MixedIntegerLinearProgram(maximization = False, solver = solver)

        f = p.new_variable(binary=True)
        r = p.new_variable(nonnegative=True)

        eps = 1/(2*Integer(g.order()))
        x = g.vertex_iterator().next()


        if g.is_directed():

            # returns the variable corresponding to arc u,v
            E = lambda u,v : f[(u,v)]

            # All the vertices have in-degree 1 and out-degree 1
            for v in g:
                p.add_constraint(p.sum([ f[(u,v)] for u in g.neighbors_in(v)]),
                                 min = 1,
                                 max = 1)

                p.add_constraint(p.sum([ f[(v,u)] for u in g.neighbors_out(v)]),
                                 min = 1,
                                 max = 1)

            # r is greater than f
            for u,v in g.edges(labels = None):
                if g.has_edge(v,u):
                    if u < v:
                        p.add_constraint( r[(u,v)] + r[(v,u)]- f[(u,v)] - f[(v,u)], min = 0)

                        # no 2-cycles
                        p.add_constraint( f[(u,v)] + f[(v,u)], max = 1)

                else:
                    p.add_constraint( r[(u,v)] + r[(v,u)] - f[(u,v)], min = 0)

            # defining the answer when g is directed
            from sage.graphs.all import DiGraph
            tsp = DiGraph()
        else:

            # reorders the edge as they can appear in the two different ways
            R = lambda x,y : frozenset((x,y))

            # returns the variable corresponding to arc u,v
            E = lambda u,v : f[R(u,v)]

            # All the vertices have degree 2
            for v in g:
                p.add_constraint(p.sum([ E(u,v) for u in g.neighbors(v)]),
                                 min = 2,
                                 max = 2)

            # r is greater than f
            for u,v in g.edges(labels = None):
                p.add_constraint( r[(u,v)] + r[(v,u)] - E(u,v), min = 0)

            from sage.graphs.all import Graph

            # defining the answer when g is not directed
            tsp = Graph()

        # no cycle which does not contain x
        for v in g:
            if v != x:
                p.add_constraint(p.sum([ r[(u,v)] for u in g.neighbors(v)]),max = 1-eps)

        if use_edge_labels:
            p.set_objective(p.sum([ weight(l)*E(u,v) for u,v,l in g.edges()]) )
        else:
            p.set_objective(None)

        try:
            obj = p.solve(log = verbose)
            f = p.get_values(f)
            tsp.add_vertices(g.vertices())
            tsp.set_pos(g.get_pos())
            tsp.name("TSP from "+g.name())
            tsp.add_edges([(u,v,l) for u,v,l in g.edges() if E(u,v) == 1])

            return tsp

        except MIPSolverException:
            raise EmptySetError("The given graph is not Hamiltonian")


    def hamiltonian_cycle(self, algorithm='tsp' ):
        r"""
        Returns a Hamiltonian cycle/circuit of the current graph/digraph

        A graph (resp. digraph) is said to be Hamiltonian
        if it contains as a subgraph a cycle (resp. a circuit)
        going through all the vertices.

        Computing a Hamiltonian cycle/circuit being NP-Complete,
        this algorithm could run for some time depending on
        the instance.

        ALGORITHM:

        See ``Graph.traveling_salesman_problem`` for 'tsp' algorithm and
        ``find_hamiltonian`` from ``sage.graphs.generic_graph_pyx``
        for 'backtrack' algorithm.

        INPUT:

            - ``algorithm`` - one of 'tsp' or 'backtrack'.

        OUTPUT:

        If using the 'tsp' algorithm, returns a Hamiltonian cycle/circuit if it
        exists; otherwise, raises a ``EmptySetError`` exception. If using the
        'backtrack' algorithm, returns a pair (B,P). If B is True then P is a
        Hamiltonian cycle and if B is False, P is a longest path found by the
        algorithm. Observe that if B is False, the graph may still be Hamiltonian.
        The 'backtrack' algorithm is only implemented for undirected
        graphs.

        .. WARNING::

            The 'backtrack' algorithm may loop endlessly on graphs
            with vertices of degree 1.

        NOTE:

        This function, as ``is_hamiltonian``, computes a Hamiltonian
        cycle if it exists: the user should *NOT* test for
        Hamiltonicity using ``is_hamiltonian`` before calling this
        function, as it would result in computing it twice.

        The backtrack algorithm is only implemented for undirected graphs.

        EXAMPLES:

        The Heawood Graph is known to be Hamiltonian ::

            sage: g = graphs.HeawoodGraph()
            sage: g.hamiltonian_cycle()
            TSP from Heawood graph: Graph on 14 vertices

        The Petersen Graph, though, is not ::

            sage: g = graphs.PetersenGraph()
            sage: g.hamiltonian_cycle()
            Traceback (most recent call last):
            ...
            EmptySetError: The given graph is not hamiltonian

        Now, using the backtrack algorithm in the Heawood graph ::

            sage: G=graphs.HeawoodGraph()
            sage: G.hamiltonian_cycle(algorithm='backtrack')
            (True, [11, 10, 1, 2, 3, 4, 9, 8, 7, 6, 5, 0, 13, 12])

        And now in the Petersen graph ::

            sage: G=graphs.PetersenGraph()
            sage: G.hamiltonian_cycle(algorithm='backtrack')
            (False, [6, 8, 5, 0, 1, 2, 7, 9, 4, 3])

        Finally, we test the algorithm in a cube graph, which is Hamiltonian ::

            sage: G=graphs.CubeGraph(3)
            sage: G.hamiltonian_cycle(algorithm='backtrack')
            (True, ['010', '110', '100', '000', '001', '101', '111', '011'])

        """
        if algorithm=='tsp':
            from sage.numerical.mip import MIPSolverException

            try:
                return self.traveling_salesman_problem(use_edge_labels = False)
            except MIPSolverException:
                from sage.categories.sets_cat import EmptySetError
                raise EmptySetError("The given graph is not Hamiltonian")
        elif algorithm=='backtrack':
            from sage.graphs.generic_graph_pyx import find_hamiltonian as fh
            return fh( self )

        else:
            raise ValueError("``algorithm`` (%s) should be 'tsp' or 'backtrack'."%(algorithm))

    def feedback_vertex_set(self, value_only=False, solver=None, verbose=0, constraint_generation = True):
        r"""
        Computes the minimum feedback vertex set of a (di)graph.

        The minimum feedback vertex set of a (di)graph is a set of vertices that
        intersect all of its cycles.  Equivalently, a minimum feedback vertex
        set of a (di)graph is a set `S` of vertices such that the digraph `G-S`
        is acyclic. For more information, see :wikipedia:`Feedback_vertex_set`.

        INPUT:

        - ``value_only`` -- boolean (default: ``False``)

          - When set to ``True``, only the minimum cardinal of a minimum vertex
            set is returned.

          - When set to ``False``, the ``Set`` of vertices of a minimal feedback
            vertex set is returned.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used,
          see the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``constraint_generation`` (boolean) -- whether to use constraint
          generation when solving the Mixed Integer Linear Program (default:
          ``True``).

        ALGORITHMS:

        (Constraints generation)

        When the parameter ``constraint_generation`` is enabled (default) the
        following MILP formulation is used to solve the problem:

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_{v}\\
            \mbox{Such that : }&\\
            &\forall C\text{ circuits }\subseteq G, \sum_{v\in C}b_{v}\geq 1\\

        As the number of circuits contained in a graph is exponential, this LP
        is solved through constraint generation. This means that the solver is
        sequentially asked to solve the problem, knowing only a portion of the
        circuits contained in `G`, each time adding to the list of its
        constraints the circuit which its last answer had left intact.

        (Another formulation based on an ordering of the vertices)

        When the graph is directed, a second (and very slow) formulation is
        available, which should only be used to check the result of the first
        implementation in case of doubt.

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\\
            &\forall (u,v)\in G, d_u-d_v+nb_u+nb_v\geq 0\\
            &\forall u\in G, 0\leq d_u\leq |G|\\

        A brief explanation:

        An acyclic digraph can be seen as a poset, and every poset has a linear
        extension. This means that in any acyclic digraph the vertices can be
        ordered with a total order `<` in such a way that if `(u,v)\in G`, then
        `u<v`.  Thus, this linear program is built in order to assign to each
        vertex `v` a number `d_v\in [0,\dots,n-1]` such that if there exists an
        edge `(u,v)\in G` then either `d_v<d_u` or one of `u` or `v` is removed.
        The number of vertices removed is then minimized, which is the
        objective.

        EXAMPLES:

        The necessary example::

            sage: g = graphs.PetersenGraph()
            sage: fvs = g.feedback_vertex_set()
            sage: len(fvs)
            3
            sage: g.delete_vertices(fvs)
            sage: g.is_forest()
            True

        In a digraph built from a graph, any edge is replaced by arcs going in
        the two opposite directions, thus creating a cycle of length two.
        Hence, to remove all the cycles from the graph, each edge must see one
        of its neighbors removed: a feedback vertex set is in this situation a
        vertex cover::

            sage: cycle = graphs.CycleGraph(5)
            sage: dcycle = DiGraph(cycle)
            sage: cycle.vertex_cover(value_only=True)
            3
            sage: feedback = dcycle.feedback_vertex_set()
            sage: len(feedback)
            3
            sage: (u,v,l) = cycle.edge_iterator().next()
            sage: u in feedback or v in feedback
            True

        For a circuit, the minimum feedback arc set is clearly `1`::

            sage: circuit = digraphs.Circuit(5)
            sage: circuit.feedback_vertex_set(value_only=True) == 1
            True

        TESTS:

        Comparing with/without constraint generation::

            sage: g = digraphs.RandomDirectedGNP(10,.3)
            sage: x = g.feedback_vertex_set(value_only = True)
            sage: y = g.feedback_vertex_set(value_only = True,
            ....:          constraint_generation = False)
            sage: x == y
            True

        Bad algorithm::

            sage: g = graphs.PetersenGraph()
            sage: g.feedback_vertex_set(constraint_generation = False)
            Traceback (most recent call last):
            ...
            ValueError: The only implementation available for undirected graphs is with constraint_generation set to True.
        """
        if not constraint_generation and not self.is_directed():
            raise ValueError("The only implementation available for "
                             "undirected graphs is with constraint_generation "
                             "set to True.")

        # It would be a pity to start a LP if the graph is already acyclic
        if ((not self.is_directed() and self.is_forest()) or
            (    self.is_directed() and self.is_directed_acyclic())):
            if value_only:
                return 0
            return []

        from sage.numerical.mip import MixedIntegerLinearProgram

        ########################################
        # Constraint Generation Implementation #
        ########################################
        if constraint_generation:

            p = MixedIntegerLinearProgram(constraint_generation = True,
                                          maximization = False)

            # A variable for each vertex
            b = p.new_variable(binary = True)

            # Variables are binary, and their coefficient in the objective is 1

            p.set_objective(p.sum( b[v] for v in self))

            p.solve(log = verbose)

            # For as long as we do not break because the digraph is
            # acyclic....
            while True:

                # Building the graph without the edges removed by the LP
                h = self.subgraph(vertices =
                                  [v for v in self if p.get_values(b[v]) == 0])

                # Is the graph acyclic ?
                if self.is_directed():
                    isok, certificate = h.is_directed_acyclic(certificate = True)
                else:
                    isok, certificate = h.is_forest(certificate = True)

                # If so, we are done !
                if isok:
                    break

                if verbose:
                    print "Adding a constraint on circuit: ",certificate

                # There is a circuit left. Let's add the corresponding
                # constraint !

                p.add_constraint(p.sum(b[v] for v in certificate), min = 1)

                obj = p.solve(log = verbose)

            if value_only:
                return obj

            else:

                # listing the edges contained in the MFVS
                return [v for v in self if p.get_values(b[v]) == 1]

        else:

        ######################################
        # Ordering-based MILP Implementation #
        ######################################

            p = MixedIntegerLinearProgram(maximization = False, solver = solver)

            b = p.new_variable(binary = True)
            d = p.new_variable(integer = True)
            n = self.order()

            # The removed vertices cover all the back arcs ( third condition )
            for (u,v) in self.edges(labels = None):
                p.add_constraint(d[u]-d[v]+n*(b[u]+b[v]), min = 1)

            for u in self:
                p.add_constraint(d[u], max = n)

            p.set_objective(p.sum([b[v] for v in self]))

            if value_only:
                return Integer(round(p.solve(objective_only = True, log = verbose)))
            else:
                p.solve(log=verbose)
                b_sol = p.get_values(b)

                return [v for v in self if b_sol[v] == 1]

    def flow(self, x, y, value_only=True, integer=False, use_edge_labels=True, vertex_bound=False, method = None, solver=None, verbose=0):
        r"""
        Returns a maximum flow in the graph from ``x`` to ``y``
        represented by an optimal valuation of the edges. For more
        information, see the
        `Wikipedia article on maximum flow
        <http://en.wikipedia.org/wiki/Max_flow>`_.

        As an optimization problem, is can be expressed this way :

        .. MATH::

            \mbox{Maximize : }&\sum_{e\in G.edges()} w_e b_e\\
            \mbox{Such that : }&\forall v \in G, \sum_{(u,v)\in G.edges()} b_{(u,v)}\leq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``x`` -- Source vertex

        - ``y`` -- Sink vertex

        - ``value_only`` -- boolean (default: ``True``)

          - When set to ``True``, only the value of a maximal
            flow is returned.

          - When set to ``False``, is returned a pair whose first element
            is the value of the maximum flow, and whose second value is
            a flow graph (a copy of the current graph, such that each edge
            has the flow using it as a label, the edges without flow being
            omitted).

        - ``integer`` -- boolean (default: ``False``)

          - When set to ``True``, computes an optimal solution under the
            constraint that the flow going through an edge has to be an
            integer.

        - ``use_edge_labels`` -- boolean (default: ``True``)

          - When set to ``True``, computes a maximum flow
            where each edge has a capacity defined by its label. (If
            an edge has no label, `1` is assumed.)

          - When set to ``False``, each edge has capacity `1`.

        - ``vertex_bound`` -- boolean (default: ``False``)

          - When set to ``True``, sets the maximum flow leaving
            a vertex different from `x` to `1` (useful for vertex
            connectivity parameters).

        - ``method`` -- There are currently two different
          implementations of this method :

              * If ``method = "FF"``, a Python implementation of the
                Ford-Fulkerson algorithm is used (only available when
                ``vertex_bound = False``)

              * If ``method = "LP"``, the flow problem is solved using
                Linear Programming.

              * If ``method = None`` (default), the Ford-Fulkerson
                implementation is used iif ``vertex_bound = False``.

        - ``solver`` -- Specify a Linear Program solver to be used.
          If set to ``None``, the default one is used.  function of
          ``MixedIntegerLinearProgram``. See the documentation of
          ``MixedIntegerLinearProgram.solve`` for more information.

          Only useful when LP is used to solve the flow problem.

        - ``verbose`` (integer) -- sets the level of verbosity. Set to 0
          by default (quiet).

          Only useful when LP is used to solve the flow problem.

        .. NOTE::

           Even though the two different implementations are meant to
           return the same Flow values, they can not be expected to
           return the same Flow graphs.

           Besides, the use of Linear Programming may possibly mean a
           (slight) numerical noise.

        EXAMPLES:

        Two basic applications of the flow method for the ``PappusGraph`` and the
        ``ButterflyGraph`` with parameter `2` ::

           sage: g=graphs.PappusGraph()
           sage: g.flow(1,2)
           3

        ::

           sage: b=digraphs.ButterflyGraph(2)
           sage: b.flow(('00',1),('00',2))
           1

        The flow method can be used to compute a matching in a bipartite graph
        by linking a source `s` to all the vertices of the first set and linking
        a sink `t` to all the vertices of the second set, then computing
        a maximum `s-t` flow ::

            sage: g = DiGraph()
            sage: g.add_edges([('s',i) for i in range(4)])
            sage: g.add_edges([(i,4+j) for i in range(4) for j in range(4)])
            sage: g.add_edges([(4+i,'t') for i in range(4)])
            sage: [cardinal, flow_graph] = g.flow('s','t',integer=True,value_only=False)
            sage: flow_graph.delete_vertices(['s','t'])
            sage: len(flow_graph.edges())
            4

        TESTS:

        An exception if raised when forcing "FF" with ``vertex_bound = True``::

            sage: g = graphs.PetersenGraph()
            sage: g.flow(0,1,vertex_bound = True, method = "FF")
            Traceback (most recent call last):
            ...
            ValueError: This method does not support both vertex_bound=True and method="FF".

        Or if the method is different from the expected values::

            sage: g.flow(0,1, method="Divination")
            Traceback (most recent call last):
            ...
            ValueError: The method argument has to be equal to either "FF", "LP" or None

        The two methods are indeed returning the same results (possibly with
        some numerical noise, cf. :trac:`12362`)::

           sage: g = graphs.RandomGNP(20,.3)
           sage: for u,v in g.edges(labels=False):
           ...      g.set_edge_label(u,v,round(random(),5))
           sage: flow_ff = g.flow(0,1, method="FF")
           sage: flow_lp = g.flow(0,1,method="LP")
           sage: abs(flow_ff-flow_lp) < 0.01
           True
        """
        self._scream_if_not_simple(allow_loops=True)
        if vertex_bound and method == "FF":
            raise ValueError("This method does not support both vertex_bound=True and method=\"FF\".")

        if (method == "FF" or
            (method is None and not vertex_bound)):
            return self._ford_fulkerson(x,y, value_only=value_only, integer=integer, use_edge_labels=use_edge_labels)

        if method != "LP" and not method is None:
            raise ValueError("The method argument has to be equal to either \"FF\", \"LP\" or None")


        from sage.numerical.mip import MixedIntegerLinearProgram
        g=self
        p=MixedIntegerLinearProgram(maximization=True, solver = solver)
        flow=p.new_variable(nonnegative=True)

        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            capacity=lambda x: x if x in RR else 1
        else:
            capacity=lambda x: 1

        if g.is_directed():
            # This function return the balance of flow at X
            flow_sum=lambda X: p.sum([flow[(X,v)] for (u,v) in g.outgoing_edges([X],labels=None)])-p.sum([flow[(u,X)] for (u,v) in g.incoming_edges([X],labels=None)])

            # The flow leaving x
            flow_leaving = lambda X : p.sum([flow[(uu,vv)] for (uu,vv) in g.outgoing_edges([X],labels=None)])

            # The flow to be considered when defining the capacity contraints
            capacity_sum = lambda u,v : flow[(u,v)]

        else:
            # This function return the balance of flow at X
            flow_sum=lambda X:p.sum([flow[(X,v)]-flow[(v,X)] for v in g[X]])

            # The flow leaving x
            flow_leaving = lambda X : p.sum([flow[(X,vv)] for vv in g[X]])

            # The flow to be considered when defining the capacity contraints
            capacity_sum = lambda u,v : flow[(u,v)] + flow[(v,u)]

        # Maximizes the flow leaving x
        p.set_objective(flow_sum(x))

        # Elsewhere, the flow is equal to 0
        for v in g:
            if v!=x and v!=y:
                p.add_constraint(flow_sum(v),min=0,max=0)

        # Capacity constraints
        for (u,v,w) in g.edges():
            p.add_constraint(capacity_sum(u,v),max=capacity(w))

        # No vertex except the sources can send more than 1
        if vertex_bound:
            for v in g:
                if v!=x and v!=y:
                    p.add_constraint(flow_leaving(v),max=1)

        if integer:
            p.set_integer(flow)


        if value_only:
            return p.solve(objective_only=True, log = verbose)

        obj=p.solve(log = verbose)

        if integer or use_edge_labels is False:
            obj = Integer(round(obj))

        flow=p.get_values(flow)
        # Builds a clean flow Draph
        flow_graph = g._build_flow_graph(flow, integer=integer)

        # Which could be a Graph
        if not self.is_directed():
            from sage.graphs.graph import Graph
            flow_graph = Graph(flow_graph)

        return [obj,flow_graph]

    def _ford_fulkerson(self, s, t, use_edge_labels = False, integer = False, value_only = True):
        r"""
        Python implementation of the Ford-Fulkerson algorithm.

        This method is a Python implementation of the Ford-Fulkerson
        max-flow algorithm, which is (slightly) faster than the LP
        implementation.

        INPUT:

        - ``s`` -- Source vertex

        - ``t`` -- Sink vertex

        - ``value_only`` -- boolean (default: ``True``)

          - When set to ``True``, only the value of a maximal
            flow is returned.

          - When set to ``False``, is returned a pair whose first element
            is the value of the maximum flow, and whose second value is
            a flow graph (a copy of the current graph, such that each edge
            has the flow using it as a label, the edges without flow being
            omitted).

        - ``integer`` -- boolean (default: ``False``)

          - When set to ``True``, computes an optimal solution under the
            constraint that the flow going through an edge has to be an
            integer.

        - ``use_edge_labels`` -- boolean (default: ``True``)

          - When set to ``True``, computes a maximum flow
            where each edge has a capacity defined by its label. (If
            an edge has no label, `1` is assumed.)

          - When set to ``False``, each edge has capacity `1`.

        EXAMPLES:

        Two basic applications of the flow method for the ``PappusGraph`` and the
        ``ButterflyGraph`` with parameter `2` ::

           sage: g=graphs.PappusGraph()
           sage: g._ford_fulkerson(1,2)
           3

        ::

           sage: b=digraphs.ButterflyGraph(2)
           sage: b._ford_fulkerson(('00',1),('00',2))
           1

        The flow method can be used to compute a matching in a bipartite graph
        by linking a source `s` to all the vertices of the first set and linking
        a sink `t` to all the vertices of the second set, then computing
        a maximum `s-t` flow ::

            sage: g = DiGraph()
            sage: g.add_edges([('s',i) for i in range(4)])
            sage: g.add_edges([(i,4+j) for i in range(4) for j in range(4)])
            sage: g.add_edges([(4+i,'t') for i in range(4)])
            sage: [cardinal, flow_graph] = g._ford_fulkerson('s','t',integer=True,value_only=False)
            sage: flow_graph.delete_vertices(['s','t'])
            sage: len(flow_graph.edges(labels=None))
            4
        """
        from sage.graphs.digraph import DiGraph
        from sage.functions.other import floor

        # Whether we should consider the edges labeled
        if use_edge_labels:
            l_capacity=lambda x: 1 if (x is None or x == {}) else (floor(x) if integer else x)
        else:
            l_capacity=lambda x: 1

        directed = self.is_directed()

        # Associated to each edge (u,v) of the flow graph its capacity
        capacity = {}
        # Associates to each edge (u,v) of the graph the (directed)
        # flow going through it
        flow = {}

        # Residual graph. Only contains edge on which some flow can be
        # sent. This can happen both when the flow going through the
        # current edge is strictly less than its capacity, or when
        # there exists a back arc with non-null flow
        residual = DiGraph()

        # Initializing the variables
        if directed:
            for u,v,l in self.edge_iterator():
                if l_capacity(l) > 0:
                    capacity[(u,v)] = l_capacity(l) + capacity.get((u,v),0)
                    capacity[(v,u)] = capacity.get((v,u),0)
                    residual.add_edge(u,v)
                    flow[(u,v)] = 0
                    flow[(v,u)] = 0
        else:
            for u,v,l in self.edge_iterator():
                if l_capacity(l) > 0:
                    capacity[(u,v)] = l_capacity(l) + capacity.get((u,v),0)
                    capacity[(v,u)] = l_capacity(l) + capacity.get((v,u),0)
                    residual.add_edge(u,v)
                    residual.add_edge(v,u)
                    flow[(u,v)] = 0
                    flow[(v,u)] = 0

        # Reqrites a path as a list of edges :
        # ex : [0,1,2,3,4,5] becomes [(0,1), (1,2), (2,3), (3,4), (4,5)]
        path_to_edges = lambda P : zip(P[:-1],P[1:])

        # Rewrites a path as a list of edges labeled with their
        # available capacity
        path_to_labelled_edges = lambda P : map(lambda x_y: (x_y[0], x_y[1], capacity[(x_y[0], x_y[1])] - flow[(x_y[0], x_y[1])] + flow[(x_y[1], x_y[0])]),
                                                path_to_edges(P))

        # Total flow going from s to t
        flow_intensity = 0

        while True:

            # If there is a shortest path from s to t
            path = residual.shortest_path(s,t)
            if not path:
                break

            # We are rewriting the shortest path as a sequence of
            # edges whose labels are their available capacities
            edges = path_to_labelled_edges(path)

            # minimum capacity available on the whole path
            epsilon = min(map( lambda x : x[2], edges))

            flow_intensity = flow_intensity + epsilon

            # Updating variables
            for uu,vv,ll in edges:

                # The flow on the back arc
                other = flow[(vv,uu)]
                flow[(uu,vv)] = flow[(uu,vv)] + max(0,epsilon-other)
                flow[(vv,uu)] = other - min(other, epsilon)

                # If the current edge is fully used, we do not need it
                # anymore in the residual graph
                if capacity[(uu,vv)] - flow[(uu,vv)] + flow[(vv,uu)] == 0:
                    residual.delete_edge(uu,vv)

                # If the back arc does not exist, it now does as the
                # edge (uu,vv) has a flow of at least epsilon>0
                if not residual.has_edge(vv,uu):
                    residual.add_edge(vv,uu,epsilon)

        if value_only:
            return flow_intensity

        # Building and returning the flow graph
        g = DiGraph()
        g.add_edges([(x,y,l) for ((x,y),l) in flow.iteritems() if l > 0])
        g.set_pos(self.get_pos())

        return flow_intensity, g

    def multicommodity_flow(self, terminals, integer=True, use_edge_labels=False,vertex_bound=False, solver=None, verbose=0):
        r"""
        Solves a multicommodity flow problem.

        In the multicommodity flow problem, we are given a set of pairs
        `(s_i, t_i)`, called terminals meaning that `s_i` is willing
        some flow to `t_i`.

        Even though it is a natural generalisation of the flow problem
        this version of it is NP-Complete to solve when the flows
        are required to be integer.

        For more information, see the
        :wikipedia:`Wikipedia page on multicommodity flows
        <Multi-commodity_flow_problem>`.

        INPUT:

        - ``terminals`` -- a list of pairs `(s_i, t_i)` or triples
          `(s_i, t_i, w_i)` representing a flow from `s_i` to `t_i`
          of intensity `w_i`. When the pairs are of size `2`, a intensity
          of `1` is assumed.

        - ``integer`` (boolean) -- whether to require an integer multicommodity
          flow

        - ``use_edge_labels`` (boolean) -- whether to consider the label of edges
          as numerical values representing a capacity. If set to ``False``, a capacity
          of `1` is assumed

        - ``vertex_bound`` (boolean) -- whether to require that a vertex can stand at most
          `1` commodity of flow through it of intensity `1`. Terminals can obviously
          still send or receive several units of flow even though vertex_bound is set
          to ``True``, as this parameter is meant to represent topological properties.

        - ``solver`` -- Specify a Linear Program solver to be used.
          If set to ``None``, the default one is used.
          function of ``MixedIntegerLinearProgram``. See the documentation  of ``MixedIntegerLinearProgram.solve``
          for more informations.

        - ``verbose`` (integer) -- sets the level of verbosity. Set to 0
          by default (quiet).

        ALGORITHM:

        (Mixed Integer) Linear Program, depending on the value of ``integer``.

        EXAMPLE:

        An easy way to obtain a satisfiable multiflow is to compute
        a matching in a graph, and to consider the paired vertices
        as terminals ::

            sage: g = graphs.PetersenGraph()
            sage: matching = [(u,v) for u,v,_ in g.matching()]
            sage: h = g.multicommodity_flow(matching)
            sage: len(h)
            5

        We could also have considered ``g`` as symmetric and computed
        the multiflow in this version instead. In this case, however
        edges can be used in both directions at the same time::

            sage: h = DiGraph(g).multicommodity_flow(matching)
            sage: len(h)
            5

        An exception is raised when the problem has no solution ::

            sage: h = g.multicommodity_flow([(u,v,3) for u,v in matching])
            Traceback (most recent call last):
            ...
            EmptySetError: The multiflow problem has no solution
        """
        self._scream_if_not_simple(allow_loops=True)
        from sage.numerical.mip import MixedIntegerLinearProgram
        g=self
        p=MixedIntegerLinearProgram(maximization=True, solver = solver)

        # Adding the intensity if not present
        terminals = [(x if len(x) == 3 else (x[0],x[1],1)) for x in terminals]

        # defining the set of terminals
        set_terminals = set([])
        for s,t,_ in terminals:
            set_terminals.add(s)
            set_terminals.add(t)

        # flow[i,(u,v)] is the flow of commodity i going from u to v
        flow=p.new_variable(nonnegative=True)

        # Whether to use edge labels
        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            capacity=lambda x: x if x in RR else 1
        else:
            capacity=lambda x: 1

        if g.is_directed():
            # This function return the balance of flow at X
            flow_sum=lambda i,X: p.sum([flow[i,(X,v)] for (u,v) in g.outgoing_edges([X],labels=None)])-p.sum([flow[i,(u,X)] for (u,v) in g.incoming_edges([X],labels=None)])

            # The flow leaving x
            flow_leaving = lambda i,X : p.sum([flow[i,(uu,vv)] for (uu,vv) in g.outgoing_edges([X],labels=None)])

            # the flow to consider when defining the capacity contraints
            capacity_sum = lambda i,u,v : flow[i,(u,v)]

        else:
            # This function return the balance of flow at X
            flow_sum=lambda i,X:p.sum([flow[i,(X,v)]-flow[i,(v,X)] for v in g[X]])

            # The flow leaving x
            flow_leaving = lambda i, X : p.sum([flow[i,(X,vv)] for vv in g[X]])

            # the flow to consider when defining the capacity contraints
            capacity_sum = lambda i,u,v : flow[i,(u,v)] + flow[i,(v,u)]


        # Flow constraints
        for i,(s,t,l) in enumerate(terminals):
            for v in g:
                if v == s:
                    p.add_constraint(flow_sum(i,v),min=l,max=l)
                elif v == t:
                    p.add_constraint(flow_sum(i,v),min=-l,max=-l)
                else:
                    p.add_constraint(flow_sum(i,v),min=0,max=0)

        # Capacity constraints
        for (u,v,w) in g.edges():
            p.add_constraint(p.sum([capacity_sum(i,u,v) for i in range(len(terminals))]),max=capacity(w))

        if vertex_bound:

            # Any vertex
            for v in g.vertices():

                # which is an endpoint
                if v in set_terminals:
                    for i,(s,t,_) in enumerate(terminals):

                        # only tolerates the commodities of which it is an endpoint
                        if not (v==s or v==t):
                            p.add_constraint(flow_leaving(i,v), max = 0)

                # which is not an endpoint
                else:
                    # can stand at most 1 unit of flow through itself
                    p.add_constraint(p.sum([flow_leaving(i,v) for i in range(len(terminals))]), max = 1)

        p.set_objective(None)

        if integer:
            p.set_integer(flow)

        from sage.numerical.mip import MIPSolverException

        try:
            obj=p.solve(log = verbose)
        except MIPSolverException:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("The multiflow problem has no solution")

        flow=p.get_values(flow)

        # building clean flow digraphs
        flow_graphs = [g._build_flow_graph({e:f for (ii,e),f in flow.iteritems() if ii == i}, integer=integer)
                       for i in range(len(terminals))]

        # which could be .. graphs !
        if not self.is_directed():
            from sage.graphs.graph import Graph
            flow_graphs = map(Graph, flow_graphs)

        return flow_graphs

    def _build_flow_graph(self, flow, integer):
        r"""
        Builds a "clean" flow graph

        It build it, then looks for circuits and removes them

        INPUT:

        - ``flow`` -- a dictionary associating positive numerical values
          to edges

        - ``integer`` (boolean) -- whether the values from ``flow`` are the solution
          of an integer flow. In this case, a value of less than .5 is assumed to be 0


        EXAMPLE:

        This method is tested in ``flow`` and ``multicommodity_flow``::

            sage: g = Graph()
            sage: g.add_edge(0,1)
            sage: f = g._build_flow_graph({(0,1):1}, True)

        The method removes zero-cost flow cycles and updates the values accordingly::

            sage: g = digraphs.DeBruijn(2,3)
            sage: flow = {('001','010'):1,('010','100'):1,('010','101'):1,('101','010'):1}
            sage: flow_graph = g._build_flow_graph(flow,True)
            sage: flow_graph.edges()
            [('001', '010', 1), ('010', '100', 1)]
            sage: flow = {('001','010'):2,('010','101'):3,('101','011'):2,('101','010'):1}
            sage: flow_graph = g._build_flow_graph(flow,True)
            sage: flow_graph.edges()
            [('001', '010', 2), ('010', '101', 2), ('101', '011', 2)]

        Isolated zero-cost flow cycles are also removed::

            sage: g = digraphs.DeBruijn(2,3)
            sage: flow = {('000','001'):1,('010','101'):1,('101','010'):1}
            sage: flow_graph = g._build_flow_graph(flow,True)
            sage: flow_graph.edges()
            [('000', '001', 1)]
        """

        from sage.graphs.digraph import DiGraph
        g = DiGraph()

        # add significant edges
        for (u,v),l in flow.iteritems():
            if l > 0 and not (integer and l<.5):
                g.add_edge(u,v,l)

        # stupid way to find Cycles. Will be fixed by #8932
        # for any vertex, for any of its in-neighbors, tried to find a cycle
        for v in g:
            for u in g.neighbor_in_iterator(v):

                # the edge from u to v could have been removed in a previous iteration
                if not g.has_edge(u,v):
                    break
                sp = g.shortest_path(v,u)
                if sp != []:

                    #find the minimm value along the cycle.
                    m = g.edge_label(u,v)
                    for i in range(len(sp)-1):
                        m = min(m,g.edge_label(sp[i],sp[i+1]))

                    # removes it from all the edges of the cycle
                    sp.append(v)
                    for i in range(len(sp)-1):
                        l = g.edge_label(sp[i],sp[i+1]) - m

                        # an edge with flow 0 is removed
                        if l == 0:
                            g.delete_edge(sp[i],sp[i+1])
                        else:
                            g.set_edge_label(sp[i],sp[i+1],l)

        # if integer is set, round values and deletes zeroes
        if integer:
            for (u,v,l) in g.edges():
                if l<.5:
                    g.delete_edge(u,v)
                else:
                    g.set_edge_label(u,v, int(round(l)))

        # returning a graph with the same embedding, the corresponding name, etc ...
        h = self.subgraph(edges=[])
        h.delete_vertices([v for v in self if (v not in g) or g.degree(v)==0])
        h.add_edges(g.edges())

        return h

    def disjoint_routed_paths(self,pairs, solver=None, verbose=0):
        r"""
        Returns a set of disjoint routed paths.

        Given a set of pairs `(s_i,t_i)`, a set
        of disjoint routed paths is a set of
        `s_i-t_i` paths which can interset at their endpoints
        and are vertex-disjoint otherwise.

        INPUT:

        - ``pairs`` -- list of pairs of vertices

        - ``solver`` -- Specify a Linear Program solver to be used.
          If set to ``None``, the default one is used.
          function of ``MixedIntegerLinearProgram``. See the documentation  of ``MixedIntegerLinearProgram.solve``
          for more informations.

        - ``verbose`` (integer) -- sets the level of verbosity. Set to `0`
          by default (quiet).

        EXAMPLE:

        Given a grid, finding two vertex-disjoint
        paths, the first one from the top-left corner
        to the bottom-left corner, and the second from
        the top-right corner to the bottom-right corner
        is easy ::

            sage: g = graphs.GridGraph([5,5])
            sage: p1,p2 = g.disjoint_routed_paths( [((0,0), (0,4)), ((4,4), (4,0))])

        Though there is obviously no solution to the problem
        in which each corner is sending information to the opposite
        one::

            sage: g = graphs.GridGraph([5,5])
            sage: p1,p2 = g.disjoint_routed_paths( [((0,0), (4,4)), ((0,4), (4,0))])
            Traceback (most recent call last):
            ...
            EmptySetError: The disjoint routed paths do not exist.
        """
        from sage.categories.sets_cat import EmptySetError
        try:
            return self.multicommodity_flow(pairs, vertex_bound = True, solver=solver, verbose=verbose)
        except EmptySetError:
            raise EmptySetError("The disjoint routed paths do not exist.")

    def edge_disjoint_paths(self, s, t, method = "FF"):
        r"""
        Returns a list of edge-disjoint paths between two
        vertices as given by Menger's theorem.

        The edge version of Menger's theorem asserts that the size
        of the minimum edge cut between two vertices `s` and`t`
        (the minimum number of edges whose removal disconnects `s`
        and `t`) is equal to the maximum number of pairwise
        edge-independent paths from `s` to `t`.

        This function returns a list of such paths.

        INPUT:

        - ``method`` -- There are currently two different
          implementations of this method :

              * If ``method = "FF"`` (default), a Python implementation of the
                Ford-Fulkerson algorithm is used.

              * If ``method = "LP"``, the flow problem is solved using
                Linear Programming.

        .. NOTE::

            This function is topological: it does not take the eventual
            weights of the edges into account.

        EXAMPLE:

        In a complete bipartite graph ::

            sage: g = graphs.CompleteBipartiteGraph(2,3)
            sage: g.edge_disjoint_paths(0,1)
            [[0, 2, 1], [0, 3, 1], [0, 4, 1]]
        """

        [obj, flow_graph] = self.flow(s,t,value_only=False, integer=True, use_edge_labels=False, method=method)

        paths = []

        while True:
            path = flow_graph.shortest_path(s,t)
            if not path:
                break
            v = s
            edges = []
            for w in path:
                edges.append((v,w))
                v=w
            flow_graph.delete_edges(edges)
            paths.append(path)

        return paths

    def vertex_disjoint_paths(self, s, t):
        r"""
        Returns a list of vertex-disjoint paths between two
        vertices as given by Menger's theorem.

        The vertex version of Menger's theorem asserts that the size
        of the minimum vertex cut between two vertices `s` and`t`
        (the minimum number of vertices whose removal disconnects `s`
        and `t`) is equal to the maximum number of pairwise
        vertex-independent paths from `s` to `t`.

        This function returns a list of such paths.

        EXAMPLE:

        In a complete bipartite graph ::

            sage: g = graphs.CompleteBipartiteGraph(2,3)
            sage: g.vertex_disjoint_paths(0,1)
            [[0, 2, 1], [0, 3, 1], [0, 4, 1]]
        """

        [obj, flow_graph] = self.flow(s,t,value_only=False, integer=True, use_edge_labels=False, vertex_bound=True)

        paths = []

        while True:
            path = flow_graph.shortest_path(s,t)
            if not path:
                break
            flow_graph.delete_vertices(path[1:-1])
            paths.append(path)

        return paths

    def dominating_set(self, independent=False, value_only=False, solver=None, verbose=0):
        r"""
        Returns a minimum dominating set of the graph
        represented by the list of its vertices. For more information, see the
        `Wikipedia article on dominating sets
        <http://en.wikipedia.org/wiki/Dominating_set>`_.

        A minimum dominating set `S` of a graph `G` is
        a set of its vertices of minimal cardinality such
        that any vertex of `G` is in `S` or has one of its neighbors
        in `S`.

        As an optimization problem, it can be expressed as:

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\forall v \in G, b_v+\sum_{(u,v)\in G.edges()} b_u\geq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``independent`` -- boolean (default: ``False``). If
          ``independent=True``, computes a minimum independent dominating set.

        - ``value_only`` -- boolean (default: ``False``)

          - If ``True``, only the cardinality of a minimum
            dominating set is returned.

          - If ``False`` (default), a minimum dominating set
            is returned as the list of its vertices.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        A basic illustration on a ``PappusGraph``::

           sage: g=graphs.PappusGraph()
           sage: g.dominating_set(value_only=True)
           5

        If we build a graph from two disjoint stars, then link their centers
        we will find a difference between the cardinality of an independent set
        and a stable independent set::

           sage: g = 2 * graphs.StarGraph(5)
           sage: g.add_edge(0,6)
           sage: len(g.dominating_set())
           2
           sage: len(g.dominating_set(independent=True))
           6

        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        g=self
        p=MixedIntegerLinearProgram(maximization=False, solver=solver)
        b=p.new_variable(binary=True)

        # For any vertex v, one of its neighbors or v itself is in
        # the minimum dominating set
        for v in g.vertices():
            p.add_constraint(b[v]+p.sum([b[u] for u in g.neighbors(v)]),min=1)


        if independent:
            # no two adjacent vertices are in the set
            for (u,v) in g.edges(labels=None):
                p.add_constraint(b[u]+b[v],max=1)

        # Minimizes the number of vertices used
        p.set_objective(p.sum([b[v] for v in g.vertices()]))

        if value_only:
            return Integer(round(p.solve(objective_only=True, log=verbose)))
        else:
            p.solve(log=verbose)
            b=p.get_values(b)
            return [v for v in g.vertices() if b[v]==1]

    def edge_connectivity(self, value_only=True, use_edge_labels=False, vertices=False, solver=None, verbose=0):
        r"""
        Returns the edge connectivity of the graph. For more information, see
        the
        `Wikipedia article on connectivity
        <http://en.wikipedia.org/wiki/Connectivity_(graph_theory)>`_.

        .. NOTE::

            When the graph is a directed graph, this method actually computes
            the *strong* connectivity, (i.e. a directed graph is strongly
            `k`-connected if there are `k` disjoint paths between any two
            vertices `u, v`). If you do not want to consider strong
            connectivity, the best is probably to convert your ``DiGraph``
            object to a ``Graph`` object, and compute the connectivity of this
            other graph.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``)

          - When set to ``True`` (default), only the value is returned.

          - When set to ``False``, both the value and a minimum edge cut
            are returned.

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - When set to ``True``, computes a weighted minimum cut
            where each edge has a weight defined by its label. (If
            an edge has no label, `1` is assumed.)

          - When set to ``False``, each edge has weight `1`.

        - ``vertices`` -- boolean (default: ``False``)

          - When set to ``True``, also returns the two sets of
            vertices that are disconnected by the cut. Implies
            ``value_only=False``.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        A basic application on the PappusGraph::

           sage: g = graphs.PappusGraph()
           sage: g.edge_connectivity()
           3

        The edge connectivity of a complete graph ( and of a random graph )
        is its minimum degree, and one of the two parts of the bipartition
        is reduced to only one vertex. The cutedges isomorphic to a
        Star graph::

           sage: g = graphs.CompleteGraph(5)
           sage: [ value, edges, [ setA, setB ]] = g.edge_connectivity(vertices=True)
           sage: print value
           4
           sage: len(setA) == 1 or len(setB) == 1
           True
           sage: cut = Graph()
           sage: cut.add_edges(edges)
           sage: cut.is_isomorphic(graphs.StarGraph(4))
           True

        Even if obviously in any graph we know that the edge connectivity
        is less than the minimum degree of the graph::

           sage: g = graphs.RandomGNP(10,.3)
           sage: min(g.degree()) >= g.edge_connectivity()
           True

        If we build a tree then assign to its edges a random value, the
        minimum cut will be the edge with minimum value::

           sage: g = graphs.RandomGNP(15,.5)
           sage: tree = Graph()
           sage: tree.add_edges(g.min_spanning_tree())
           sage: for u,v in tree.edge_iterator(labels=None):
           ...        tree.set_edge_label(u,v,random())
           sage: minimum = min([l for u,v,l in tree.edge_iterator()])
           sage: [value, [(u,v,l)]] = tree.edge_connectivity(value_only=False, use_edge_labels=True)
           sage: l == minimum
           True

        When ``value_only = True``, this function is optimized for small
        connectivity values and does not need to build a linear program.

        It is the case for connected graphs which are not
        connected ::

           sage: g = 2 * graphs.PetersenGraph()
           sage: g.edge_connectivity()
           0.0

        Or if they are just 1-connected ::

           sage: g = graphs.PathGraph(10)
           sage: g.edge_connectivity()
           1.0

        For directed graphs, the strong connectivity is tested
        through the dedicated function ::

           sage: g = digraphs.ButterflyGraph(3)
           sage: g.edge_connectivity()
           0.0
        """
        self._scream_if_not_simple(allow_loops=True)
        g=self

        if vertices:
            value_only=False

        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            weight=lambda x: x if x in RR else 1
        else:
            weight=lambda x: 1


        # Better methods for small connectivity tests,
        # when one is not interested in cuts...
        if value_only and not use_edge_labels:

            if self.is_directed():
                if not self.is_strongly_connected():
                    return 0.0

            else:
                if not self.is_connected():
                    return 0.0

                h = self.strong_orientation()
                if not h.is_strongly_connected():
                    return 1.0


        if g.is_directed():
            reorder_edge = lambda x,y : (x,y)
        else:
            reorder_edge = lambda x,y : (x,y) if x<= y else (y,x)

        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)

        in_set = p.new_variable(binary = True)
        in_cut = p.new_variable(binary = True)

        # A vertex has to be in some set
        for v in g:
            p.add_constraint(in_set[0,v]+in_set[1,v],max=1,min=1)

        # There is no empty set
        p.add_constraint(p.sum([in_set[1,v] for v in g]),min=1)
        p.add_constraint(p.sum([in_set[0,v] for v in g]),min=1)

        if g.is_directed():
            # There is no edge from set 0 to set 1 which
            # is not in the cut
            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u] + in_set[1,v] - in_cut[(u,v)], max = 1)
        else:

            # Two adjacent vertices are in different sets if and only if
            # the edge between them is in the cut
            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u]+in_set[1,v]-in_cut[reorder_edge(u,v)],max=1)
                p.add_constraint(in_set[1,u]+in_set[0,v]-in_cut[reorder_edge(u,v)],max=1)

        p.set_objective(p.sum([weight(l ) * in_cut[reorder_edge(u,v)] for (u,v,l) in g.edge_iterator()]))

        obj = p.solve(objective_only=value_only, log=verbose)

        if use_edge_labels is False:
            obj = Integer(round(obj))

        if value_only:
            return obj

        else:
            val = [obj]

            in_cut = p.get_values(in_cut)
            in_set = p.get_values(in_set)

            edges = []
            for (u,v,l) in g.edge_iterator():
                if in_cut[reorder_edge(u,v)] == 1:
                    edges.append((u,v,l))

            val.append(edges)

            if vertices:
                a = []
                b = []
                for v in g:
                    if in_set[0,v] == 1:
                        a.append(v)
                    else:
                        b.append(v)
                val.append([a,b])

            return val

    def vertex_connectivity(self, value_only=True, sets=False, solver=None, verbose=0):
        r"""
        Returns the vertex connectivity of the graph. For more information,
        see the
        `Wikipedia article on connectivity
        <http://en.wikipedia.org/wiki/Connectivity_(graph_theory)>`_.

        .. NOTE::

            * When the graph is a directed graph, this method actually computes
              the *strong* connectivity, (i.e. a directed graph is strongly
              `k`-connected if there are `k` disjoint paths between any two
              vertices `u, v`). If you do not want to consider strong
              connectivity, the best is probably to convert your ``DiGraph``
              object to a ``Graph`` object, and compute the connectivity of this
              other graph.

            * By convention, a complete graph on `n` vertices is `n-1`
              connected. In this case, no certificate can be given as there is
              no pair of vertices split by a cut of size `k-1`. For this reason,
              the certificates returned in this situation are empty.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``)

          - When set to ``True`` (default), only the value is returned.

          - When set to ``False`` , both the value and a minimum vertex cut are
            returned.

        - ``sets`` -- boolean (default: ``False``)

          - When set to ``True``, also returns the two sets of
            vertices that are disconnected by the cut.
            Implies ``value_only=False``

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        A basic application on a ``PappusGraph``::

           sage: g=graphs.PappusGraph()
           sage: g.vertex_connectivity()
           3

        In a grid, the vertex connectivity is equal to the
        minimum degree, in which case one of the two sets is
        of cardinality `1`::

           sage: g = graphs.GridGraph([ 3,3 ])
           sage: [value, cut, [ setA, setB ]] = g.vertex_connectivity(sets=True)
           sage: len(setA) == 1 or len(setB) == 1
           True

        A vertex cut in a tree is any internal vertex::

           sage: g = graphs.RandomGNP(15,.5)
           sage: tree = Graph()
           sage: tree.add_edges(g.min_spanning_tree())
           sage: [val, [cut_vertex]] = tree.vertex_connectivity(value_only=False)
           sage: tree.degree(cut_vertex) > 1
           True

        When ``value_only = True``, this function is optimized for small
        connectivity values and does not need to build a linear program.

        It is the case for connected graphs which are not connected::

           sage: g = 2 * graphs.PetersenGraph()
           sage: g.vertex_connectivity()
           0

        Or if they are just 1-connected::

           sage: g = graphs.PathGraph(10)
           sage: g.vertex_connectivity()
           1

        For directed graphs, the strong connectivity is tested
        through the dedicated function::

           sage: g = digraphs.ButterflyGraph(3)
           sage: g.vertex_connectivity()
           0

        A complete graph on `10` vertices is `9`-connected::

           sage: g = graphs.CompleteGraph(10)
           sage: g.vertex_connectivity()
           9

        A complete digraph on `10` vertices is `9`-connected::

           sage: g = DiGraph(graphs.CompleteGraph(10))
           sage: g.vertex_connectivity()
           9
        """
        g=self

        if sets:
            value_only=False

        # When the graph is complete, the MILP below is infeasible.
        if ((not g.is_directed() and g.is_clique())
            or
            (g.is_directed() and g.size() >= (g.order()-1)*g.order() and
             ((not g.allows_loops() and not g.allows_multiple_edges())
              or
              all(g.has_edge(u,v) for u in g for v in g if v != u)))):
            if value_only:
                return g.order()-1
            elif not sets:
                return g.order()-1, []
            else:
                return g.order()-1, [], [[],[]]

        if value_only:
            if self.is_directed():
                if not self.is_strongly_connected():
                    return 0

            else:
                if not self.is_connected():
                    return 0

                if len(self.blocks_and_cut_vertices()[0]) > 1:
                    return 1

        if g.is_directed():
            reorder_edge = lambda x,y : (x,y)
        else:
            reorder_edge = lambda x,y : (x,y) if x<= y else (y,x)

        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)

        # Sets 0 and 2 are "real" sets while set 1 represents the cut
        in_set = p.new_variable(binary=True)

        # A vertex has to be in some set
        for v in g:
            p.add_constraint(in_set[0,v]+in_set[1,v]+in_set[2,v],max=1,min=1)

        # There is no empty set
        p.add_constraint(p.sum([in_set[0,v] for v in g]),min=1)
        p.add_constraint(p.sum([in_set[2,v] for v in g]),min=1)

        if g.is_directed():
            # There is no edge from set 0 to set 1 which
            # is not in the cut
            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u] + in_set[2,v], max = 1)
        else:

            # Two adjacent vertices are in different sets if and only if
            # the edge between them is in the cut

            for (u,v) in g.edge_iterator(labels=None):
                p.add_constraint(in_set[0,u]+in_set[2,v],max=1)
                p.add_constraint(in_set[2,u]+in_set[0,v],max=1)

        p.set_binary(in_set)

        p.set_objective(p.sum([in_set[1,v] for v in g]))

        if value_only:
            return Integer(round(p.solve(objective_only=True, log=verbose)))
        else:
            val = [Integer(round(p.solve(log=verbose)))]

            in_set = p.get_values(in_set)

            cut = []
            a = []
            b = []

            for v in g:
                if in_set[0,v] == 1:
                    a.append(v)
                elif in_set[1,v]==1:
                    cut.append(v)
                else:
                    b.append(v)

            val.append(cut)

            if sets:
                val.append([a,b])

            return val

    ### Vertex handlers

    def add_vertex(self, name=None):
        """
        Creates an isolated vertex. If the vertex already exists, then
        nothing is done.

        INPUT:

        -  ``name`` - Name of the new vertex. If no name is
           specified, then the vertex will be represented by the least integer
           not already representing a vertex. Name must be an immutable
           object, and cannot be None.

        As it is implemented now, if a graph `G` has a large number
        of vertices with numeric labels, then G.add_vertex() could
        potentially be slow, if name is None.

        OUTPUT:

        If ``name``=``None``, the new vertex name is returned. ``None`` otherwise.

        EXAMPLES::

            sage: G = Graph(); G.add_vertex(); G
            0
            Graph on 1 vertex

        ::

            sage: D = DiGraph(); D.add_vertex(); D
            0
            Digraph on 1 vertex

        """
        return self._backend.add_vertex(name)

    def add_vertices(self, vertices):
        """
        Add vertices to the (di)graph from an iterable container of
        vertices. Vertices that already exist in the graph will not be
        added again.

        INPUT:

        - ``vertices``: iterator of vertex labels. A new label is created, used and returned in
          the output list for all ``None`` values in ``vertices``.

        OUTPUT:

        Generated names of new vertices if there is at least one ``None`` value
        present in ``vertices``. ``None`` otherwise.

        EXAMPLES::

            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7,8], 6: [8,9], 7: [9]}
            sage: G = Graph(d)
            sage: G.add_vertices([10,11,12])
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            sage: G.add_vertices(graphs.CycleGraph(25).vertices())
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

        ::

            sage: G = Graph()
            sage: G.add_vertices([1,2,3])
            sage: G.add_vertices([4,None,None,5])
            [0, 6]

        """
        return self._backend.add_vertices(vertices)

    def delete_vertex(self, vertex, in_order=False):
        """
        Deletes vertex, removing all incident edges. Deleting a
        non-existent vertex will raise an exception.

        INPUT:


        -  ``in_order`` - (default False) If True, this
           deletes the ith vertex in the sorted list of vertices, i.e.
           G.vertices()[i]


        EXAMPLES::

            sage: G = Graph(graphs.WheelGraph(9))
            sage: G.delete_vertex(0); G.show()

        ::

            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertex(0); D
            Digraph on 5 vertices
            sage: D.vertices()
            [1, 2, 3, 4, 5]
            sage: D.delete_vertex(0)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (0) not in the graph.

        ::

            sage: G = graphs.CompleteGraph(4).line_graph(labels=False)
            sage: G.vertices()
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: G.delete_vertex(0, in_order=True)
            sage: G.vertices()
            [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: G = graphs.PathGraph(5)
            sage: G.set_vertices({0: 'no delete', 1: 'delete'})
            sage: G.set_boundary([1,2])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: G.delete_vertex(1)
            sage: G.get_vertices()
            {0: 'no delete', 2: None, 3: None, 4: None}
            sage: G.get_boundary()
            [2]
            sage: G.get_pos()
            {0: (0, 0), 2: (2, 0), 3: (3, 0), 4: (4, 0)}
        """
        if in_order:
            vertex = self.vertices()[vertex]
        if vertex not in self:
            raise RuntimeError("Vertex (%s) not in the graph."%str(vertex))

        self._backend.del_vertex(vertex)
        attributes_to_update = ('_pos', '_assoc', '_embedding')
        for attr in attributes_to_update:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                getattr(self, attr).pop(vertex, None)
        self._boundary = [v for v in self._boundary if v != vertex]

    def delete_vertices(self, vertices):
        """
        Remove vertices from the (di)graph taken from an iterable container
        of vertices. Deleting a non-existent vertex will raise an
        exception.

        EXAMPLES::

            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertices([1,2,3,4,5]); D
            Digraph on 1 vertex
            sage: D.vertices()
            [0]
            sage: D.delete_vertices([1])
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (1) not in the graph.

        """
        for vertex in vertices:
            if vertex not in self:
                raise RuntimeError("Vertex (%s) not in the graph."%str(vertex))

        self._backend.del_vertices(vertices)
        attributes_to_update = ('_pos', '_assoc', '_embedding')
        for attr in attributes_to_update:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                attr_dict = getattr(self, attr)
                for vertex in vertices:
                    attr_dict.pop(vertex, None)

        self._boundary = [v for v in self._boundary if v not in vertices]

    def has_vertex(self, vertex):
        """
        Return True if vertex is one of the vertices of this graph.

        INPUT:


        -  ``vertex`` - an integer


        OUTPUT:


        -  ``bool`` - True or False


        EXAMPLES::

            sage: g = Graph({0:[1,2,3], 2:[4]}); g
            Graph on 5 vertices
            sage: 2 in g
            True
            sage: 10 in g
            False
            sage: graphs.PetersenGraph().has_vertex(99)
            False
        """
        try:
            hash(vertex)
        except Exception:
            return False
        return self._backend.has_vertex(vertex)

    __contains__ = has_vertex

    def random_vertex(self, **kwds):
        r"""
        Returns a random vertex of self.

        INPUT:

        - ``**kwds`` - arguments to be passed down to the
          ``vertex_iterator`` method.

        EXAMPLE:

        The returned value is a vertex of self::

            sage: g = graphs.PetersenGraph()
            sage: v = g.random_vertex()
            sage: v in g
            True
        """
        from sage.misc.prandom import randint
        it = self.vertex_iterator(**kwds)
        for i in xrange(0, randint(0, self.order() - 1)):
            it.next()
        return it.next()

    def random_edge(self,**kwds):
        r"""
        Returns a random edge of self.

        INPUT:

        - ``**kwds`` - arguments to be passed down to the
          ``edge_iterator`` method.

        EXAMPLE:

        The returned value is an edge of self::

            sage: g = graphs.PetersenGraph()
            sage: u,v = g.random_edge(labels=False)
            sage: g.has_edge(u,v)
            True

        As the ``edges()`` method would, this function returns
        by default a triple ``(u,v,l)`` of values, in which
        ``l`` is the label of edge `(u,v)`::

            sage: g.random_edge()
            (3, 4, None)
        """
        from sage.misc.prandom import randint
        it = self.edge_iterator(**kwds)
        for i in xrange(0, randint(0, self.size() - 1)):
            it.next()
        return it.next()

    def vertex_boundary(self, vertices1, vertices2=None):
        """
        Returns a list of all vertices in the external boundary of
        vertices1, intersected with vertices2. If vertices2 is None, then
        vertices2 is the complement of vertices1. This is much faster if
        vertices1 is smaller than vertices2.

        The external boundary of a set of vertices is the union of the
        neighborhoods of each vertex in the set. Note that in this
        implementation, since vertices2 defaults to the complement of
        vertices1, if a vertex `v` has a loop, then
        vertex_boundary(v) will not contain `v`.

        In a digraph, the external boundary of a vertex v are those
        vertices u with an arc (v, u).

        EXAMPLES::

            sage: G = graphs.CubeGraph(4)
            sage: l = ['0111', '0000', '0001', '0011', '0010', '0101', '0100', '1111', '1101', '1011', '1001']
            sage: G.vertex_boundary(['0000', '1111'], l)
            ['0111', '0001', '0010', '0100', '1101', '1011']

        ::

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: D.vertex_boundary([0])
            [1, 2]
        """
        vertices1 = [v for v in vertices1 if v in self]
        output = set()
        if self._directed:
            for v in vertices1:
                output.update(self.neighbor_out_iterator(v))
        else:
            for v in vertices1:
                output.update(self.neighbor_iterator(v))
        if vertices2 is not None:
            output.intersection_update(vertices2)
        return list(output)

    def set_vertices(self, vertex_dict):
        """
        Associate arbitrary objects with each vertex, via an association
        dictionary.

        INPUT:


        -  ``vertex_dict`` - the association dictionary


        EXAMPLES::

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.set_vertices(d)
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices
        """
        if hasattr(self, '_assoc') is False:
            self._assoc = {}

        self._assoc.update(vertex_dict)

    def set_vertex(self, vertex, object):
        """
        Associate an arbitrary object with a vertex.

        INPUT:


        -  ``vertex`` - which vertex

        -  ``object`` - object to associate to vertex


        EXAMPLES::

            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.set_vertex(1, graphs.FlowerSnark())
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices
        """
        if hasattr(self, '_assoc') is False:
            self._assoc = {}

        self._assoc[vertex] = object

    def get_vertex(self, vertex):
        """
        Retrieve the object associated with a given vertex.

        INPUT:


        -  ``vertex`` - the given vertex


        EXAMPLES::

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.set_vertices(d)
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices
        """
        if hasattr(self, '_assoc') is False:
            return None

        return self._assoc.get(vertex, None)

    def get_vertices(self, verts=None):
        """
        Return a dictionary of the objects associated to each vertex.

        INPUT:


        -  ``verts`` - iterable container of vertices


        EXAMPLES::

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: T = graphs.TetrahedralGraph()
            sage: T.set_vertices(d)
            sage: T.get_vertices([1,2])
            {1: Flower Snark: Graph on 20 vertices,
             2: Moebius-Kantor Graph: Graph on 16 vertices}
        """
        if verts is None:
            verts = self.vertices()

        if hasattr(self, '_assoc') is False:
            return dict.fromkeys(verts, None)

        output = {}

        for v in verts:
            output[v] = self._assoc.get(v, None)

        return output

    def loop_vertices(self):
        """
        Returns a list of vertices with loops.

        EXAMPLES::

            sage: G = Graph({0 : [0], 1: [1,2,3], 2: [3]}, loops=True)
            sage: G.loop_vertices()
            [0, 1]
        """
        if self.allows_loops():
            return [v for v in self if self.has_edge(v,v)]
        else:
            return []

    def vertex_iterator(self, vertices=None):
        """
        Returns an iterator over the given vertices.

        Returns False if not given a vertex, sequence, iterator or None. None is
        equivalent to a list of every vertex. Note that ``for v in G`` syntax is
        allowed.

        INPUT:

        -  ``vertices`` - iterated vertices are these
           intersected with the vertices of the (di)graph

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: for v in P.vertex_iterator():
            ...    print v
            ...
            0
            1
            2
            ...
            8
            9

        ::

            sage: G = graphs.TetrahedralGraph()
            sage: for i in G:
            ...    print i
            0
            1
            2
            3

        Note that since the intersection option is available, the
        vertex_iterator() function is sub-optimal, speed-wise, but note the
        following optimization::

            sage: timeit V = P.vertices()                   # not tested
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage: timeit V = list(P.vertex_iterator())      # not tested
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage: timeit V = list(P._nxg.adj.iterkeys())    # not tested
            100000 loops, best of 3: 3.45 [micro]s per loop

        In other words, if you want a fast vertex iterator, call the
        dictionary directly.
        """
        return self._backend.iterator_verts(vertices)

    __iter__ = vertex_iterator

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors of vertex.

        EXAMPLES::

            sage: G = graphs.CubeGraph(3)
            sage: for i in G.neighbor_iterator('010'):
            ...    print i
            011
            000
            110
            sage: D = G.to_directed()
            sage: for i in D.neighbor_iterator('010'):
            ...    print i
            011
            000
            110

        ::

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: list(D.neighbor_iterator(0))
            [1, 2, 3]
        """
        return self._backend.iterator_nbrs(vertex)

    def vertices(self, key=None, boundary_first=False):
        r"""
        Return a list of the vertices.

        INPUT:

        - ``key`` - default: ``None`` - a function that takes
          a vertex as its one argument and returns a value that
          can be used for comparisons in the sorting algorithm.

        - ``boundary_first`` - default:  ``False`` - if ``True``,
          return the boundary vertices first.

        OUTPUT:

        The vertices of the list.

        .. warning::

            There is always an attempt to sort the list before
            returning the result. However, since any object may
            be a vertex, there is no guarantee that any two
            vertices will be comparable.  With default objects
            for vertices (all integers), or when all the vertices
            are of the same simple type, then there should not be
            a problem with how the vertices will be sorted. However,
            if you need to guarantee a total order for the sort,
            use the ``key`` argument, as illustrated in the
            examples below.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        If you do not care about sorted output and you are
        concerned about the time taken to sort, consider the
        following alternatives.  The moral is: if you want a
        fast vertex iterator, call the dictionary directly. ::

            sage: timeit V = P.vertices()                     # not tested
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage: timeit V = list(P.vertex_iterator())        # not tested
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage: timeit V = list(P._nxg.adj.iterkeys())      # not tested
            100000 loops, best of 3: 3.45 [micro]s per loop

        We illustrate various ways to use a ``key`` to sort the list::

            sage: H=graphs.HanoiTowerGraph(3,3,labels=False)
            sage: H.vertices()
            [0, 1, 2, 3, 4, ... 22, 23, 24, 25, 26]
            sage: H.vertices(key=lambda x: -x)
            [26, 25, 24, 23, 22, ... 4, 3, 2, 1, 0]

        ::

            sage: G=graphs.HanoiTowerGraph(3,3)
            sage: G.vertices()
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), ... (2, 2, 1), (2, 2, 2)]
            sage: G.vertices(key = lambda x: (x[1], x[2], x[0]))
            [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 0, 1), ... (1, 2, 2), (2, 2, 2)]

        The discriminant of a polynomial is a function that returns an integer.
        We build a graph whose vertices are polynomials, and use the discriminant
        function to provide an ordering.  Note that since functions are first-class
        objects in Python, we can specify precisely the function from the Sage library
        that we wish to use as the key. ::

            sage: t = polygen(QQ, 't')
            sage: K = Graph({5*t:[t^2], t^2:[t^2+2], t^2+2:[4*t^2-6], 4*t^2-6:[5*t]})
            sage: dsc = sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint.discriminant
            sage: verts = K.vertices(key=dsc)
            sage: verts
            [t^2 + 2, t^2, 5*t, 4*t^2 - 6]
            sage: [x.discriminant() for x in verts]
            [-8, 0, 1, 96]

        If boundary vertices are requested first, then they are sorted
        separately from the remainder (which are also sorted). ::

            sage: P = graphs.PetersenGraph()
            sage: P.set_boundary((5..9))
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: P.vertices(boundary_first=True)
            [5, 6, 7, 8, 9, 0, 1, 2, 3, 4]
            sage: P.vertices(boundary_first=True, key=lambda x: -x)
            [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        """
        if not boundary_first:
            return sorted(list(self.vertex_iterator()), key=key)

        from sage.misc.superseded import deprecation
        deprecation(15494, "The boundary parameter is deprecated and will soon disappear.")
        bdy_verts = []
        int_verts = []
        for v in self.vertex_iterator():
            if v in self._boundary:
                bdy_verts.append(v)
            else:
                int_verts.append(v)
        return sorted(bdy_verts, key=key) + sorted(int_verts, key=key)

    def neighbors(self, vertex):
        """
        Return a list of neighbors (in and out if directed) of vertex.

        G[vertex] also works.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: sorted(P.neighbors(3))
            [2, 4, 8]
            sage: sorted(P[4])
            [0, 3, 9]
        """
        return list(self.neighbor_iterator(vertex))

    __getitem__ = neighbors


    def merge_vertices(self,vertices):
        r"""
        Merge vertices.

        This function replaces a set `S` of vertices by a single vertex
        `v_{new}`, such that the edge `uv_{new}` exists if and only if
        `\exists v'\in S: (u,v')\in G`.

        The new vertex is named after the first vertex in the list
        given in argument. If this first name is None, a new vertex
        is created.

        In the case of multigraphs, the multiplicity is preserved.

        INPUT:

        - ``vertices`` -- the set of vertices to be merged

        EXAMPLE::

            sage: g=graphs.CycleGraph(3)
            sage: g.merge_vertices([0,1])
            sage: g.edges()
            [(0, 2, None)]

            sage: # With a Multigraph :
            sage: g=graphs.CycleGraph(3)
            sage: g.allow_multiple_edges(True)
            sage: g.merge_vertices([0,1])
            sage: g.edges(labels=False)
            [(0, 2), (0, 2)]

            sage: P=graphs.PetersenGraph()
            sage: P.merge_vertices([5,7])
            sage: P.vertices()
            [0, 1, 2, 3, 4, 5, 6, 8, 9]

            sage: g=graphs.CycleGraph(5)
            sage: g.vertices()
            [0, 1, 2, 3, 4]
            sage: g.merge_vertices([None, 1, 3])
            sage: g.edges(labels=False)
            [(0, 4), (0, 5), (2, 5), (4, 5)]

        """

        if len(vertices) > 0 and vertices[0] is None:
            vertices[0] = self.add_vertex()

        if self.is_directed():
            out_edges=self.edge_boundary(vertices)
            in_edges=self.edge_boundary([v for v in self if not v in vertices])
            self.delete_vertices(vertices[1:])
            self.add_edges([(vertices[0],v,l) for (u,v,l) in out_edges if u!=vertices[0]])
            self.add_edges([(v,vertices[0],l) for (v,u,l) in in_edges if u!=vertices[0]])
        else:
            edges=self.edge_boundary(vertices)
            self.delete_vertices(vertices[1:])
            add_edges=[]
            for (u,v,l) in edges:
                if v in vertices and v != vertices[0]:
                    add_edges.append((vertices[0],u,l))
                if u in vertices and u!=vertices[0]:
                    add_edges.append((vertices[0],v,l))
            self.add_edges(add_edges)

    ### Edge handlers

    def add_edge(self, u, v=None, label=None):
        """
        Adds an edge from u and v.

        INPUT: The following forms are all accepted:

        - G.add_edge( 1, 2 )
        - G.add_edge( (1, 2) )
        - G.add_edges( [ (1, 2) ])
        - G.add_edge( 1, 2, 'label' )
        - G.add_edge( (1, 2, 'label') )
        - G.add_edges( [ (1, 2, 'label') ] )

        WARNING: The following intuitive input results in nonintuitive
        output::

            sage: G = Graph()
            sage: G.add_edge((1,2), 'label')
            sage: G.networkx_graph().adj           # random output order
            {'label': {(1, 2): None}, (1, 2): {'label': None}}

        Use one of these instead::

            sage: G = Graph()
            sage: G.add_edge((1,2), label="label")
            sage: G.networkx_graph().adj           # random output order
            {1: {2: 'label'}, 2: {1: 'label'}}

        ::

            sage: G = Graph()
            sage: G.add_edge(1,2,'label')
            sage: G.networkx_graph().adj           # random output order
            {1: {2: 'label'}, 2: {1: 'label'}}

        The following syntax is supported, but note that you must use
        the ``label`` keyword::

            sage: G = Graph()
            sage: G.add_edge((1,2), label='label')
            sage: G.edges()
            [(1, 2, 'label')]
            sage: G = Graph()
            sage: G.add_edge((1,2), 'label')
            sage: G.edges()
            [('label', (1, 2), None)]

        Vertex name cannot be None, so::

            sage: G = Graph()
            sage: G.add_edge(None, 4)
            sage: G.vertices()
            [0, 4]
        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    try:
                        u, v = u
                    except Exception:
                        pass
        else:
            if v is None:
                try:
                    u, v = u
                except Exception:
                    pass

        self._backend.add_edge(u, v, label, self._directed)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.

        All elements of ``edges`` must follow the same format, i.e. have the
        same length.

        EXAMPLES::

            sage: G = graphs.DodecahedralGraph()
            sage: H = Graph()
            sage: H.add_edges( G.edge_iterator() ); H
            Graph on 20 vertices
            sage: G = graphs.DodecahedralGraph().to_directed()
            sage: H = DiGraph()
            sage: H.add_edges( G.edge_iterator() ); H
            Digraph on 20 vertices
            sage: H.add_edges(iter([]))

            sage: H = Graph()
            sage: H.add_edges([(0,1),(0,2)])
            sage: H.edges()
            [(0, 1, None), (0, 2, None)]
        """
        it = iter(edges)

        try:
            e0 = it.next()
        except StopIteration:
            return

        if len(e0) == 3:
            self._backend.add_edge(e0[0], e0[1], e0[2], self._directed)
            for u,v,label in it:
                self._backend.add_edge(u, v, label, self._directed)
        else:
            self._backend.add_edge(e0[0], e0[1], None, self._directed)
            for u,v in it:
                self._backend.add_edge(u, v, None, self._directed)

    def subdivide_edge(self, *args):
        """
        Subdivides an edge `k` times.

        INPUT:

        The following forms are all accepted to subdivide `8` times
        the edge between vertices `1` and `2` labeled with
        ``"my_label"``.

        - ``G.subdivide_edge( 1, 2, 8 )``
        - ``G.subdivide_edge( (1, 2), 8 )``
        - ``G.subdivide_edge( (1, 2, "my_label"), 8 )``

        .. NOTE::

            * If the given edge is labelled with `l`, all the edges
              created by the subdivision will have the same label.

            * If no label is given, the label used will be the one
              returned by the method :meth:`edge_label` on the pair
              ``u,v``

        EXAMPLE:

        Subdividing `5` times an edge in a path of length
        `3` makes it a path of length `8`::

            sage: g = graphs.PathGraph(3)
            sage: edge = g.edges()[0]
            sage: g.subdivide_edge(edge, 5)
            sage: g.is_isomorphic(graphs.PathGraph(8))
            True

        Subdividing a labelled edge in two ways ::

            sage: g = Graph()
            sage: g.add_edge(0,1,"label1")
            sage: g.add_edge(1,2,"label2")
            sage: print sorted(g.edges())
            [(0, 1, 'label1'), (1, 2, 'label2')]

        Specifying the label::

            sage: g.subdivide_edge(0,1,"label1", 3)
            sage: print sorted(g.edges())
            [(0, 3, 'label1'), (1, 2, 'label2'), (1, 5, 'label1'), (3, 4, 'label1'), (4, 5, 'label1')]

        The lazy way::

            sage: g.subdivide_edge(1,2,"label2", 5)
            sage: print sorted(g.edges())
            [(0, 3, 'label1'), (1, 5, 'label1'), (1, 6, 'label2'), (2, 10, 'label2'), (3, 4, 'label1'), (4, 5, 'label1'), (6, 7, 'label2'), (7, 8, 'label2'), (8, 9, 'label2'), (9, 10, 'label2')]

        If too many arguments are given, an exception is raised ::

            sage: g.subdivide_edge(0,1,1,1,1,1,1,1,1,1,1)
            Traceback (most recent call last):
            ...
            ValueError: This method takes at most 4 arguments !

        The same goes when the given edge does not exist::

            sage: g.subdivide_edge(0,1,"fake_label",5)
            Traceback (most recent call last):
            ...
            ValueError: The given edge does not exist.

        .. SEEALSO::

            - :meth:`subdivide_edges` -- subdivides multiples edges at a time

        TESTS:

        :trac:`15895` is fixed::

            sage: F = graphs.PathGraph(3)
            sage: S = 'S'; F.add_vertex(S)
            sage: F.add_edge(S,0)
            sage: F2 = Graph(F)
            sage: F2.subdivide_edges(list(F2.edges(labels=False)),2)
            sage: 0 in F2.degree()
            False
        """
        if len(args) == 2:
            edge, k = args

            if len(edge) == 2:
                u, v = edge
                l = self.edge_label(u, v)
            elif len(edge) == 3:
                u, v, l = edge

        elif len(args) == 3:
            u, v, k = args
            l = self.edge_label(u, v)

        elif len(args) == 4:
            u, v, l, k = args

        else:
            raise ValueError("This method takes at most 4 arguments !")

        if not self.has_edge(u, v, l):
            raise ValueError("The given edge does not exist.")

        new_verts = [self.add_vertex() for i in xrange(k)]

        self.delete_edge(u, v, l)

        edges = []
        for uu in new_verts + [v]:
            edges.append((u, uu, l))
            u = uu

        self.add_edges(edges)

    def subdivide_edges(self, edges, k):
        """
        Subdivides k times edges from an iterable container.

        For more information on the behaviour of this method, please
        refer to the documentation of :meth:`subdivide_edge`.

        INPUT:

        - ``edges`` -- a list of edges

        - ``k`` (integer) -- common length of the subdivisions


        .. NOTE::

            If a given edge is labelled with `l`, all the edges
            created by its subdivision will have the same label.

        EXAMPLE:

        If we are given the disjoint union of several paths::

            sage: paths = [2,5,9]
            sage: paths = map(graphs.PathGraph, paths)
            sage: g = Graph()
            sage: for P in paths:
            ...     g = g + P

        ... subdividing edges in each of them will only change their
        lengths::

            sage: edges = [P.edges()[0] for P in g.connected_components_subgraphs()]
            sage: k = 6
            sage: g.subdivide_edges(edges, k)

        Let us check this by creating the graph we expect to have built
        through subdivision::

            sage: paths2 = [2+k, 5+k, 9+k]
            sage: paths2 = map(graphs.PathGraph, paths2)
            sage: g2 = Graph()
            sage: for P in paths2:
            ...     g2 = g2 + P
            sage: g.is_isomorphic(g2)
            True

        .. SEEALSO::

        - :meth:`subdivide_edge` -- subdivides one edge
        """
        for e in edges:
            self.subdivide_edge(e, k)

    def delete_edge(self, u, v=None, label=None):
        r"""
        Delete the edge from u to v, returning silently if vertices or edge
        does not exist.

        INPUT: The following forms are all accepted:

        - G.delete_edge( 1, 2 )
        - G.delete_edge( (1, 2) )
        - G.delete_edges( [ (1, 2) ] )
        - G.delete_edge( 1, 2, 'label' )
        - G.delete_edge( (1, 2, 'label') )
        - G.delete_edges( [ (1, 2, 'label') ] )

        EXAMPLES::

            sage: G = graphs.CompleteGraph(19).copy(implementation='c_graph')
            sage: G.size()
            171
            sage: G.delete_edge( 1, 2 )
            sage: G.delete_edge( (3, 4) )
            sage: G.delete_edges( [ (5, 6), (7, 8) ] )
            sage: G.size()
            167

        Note that NetworkX accidentally deletes these edges, even though the
        labels do not match up::

            sage: N = graphs.CompleteGraph(19).copy(implementation='networkx')
            sage: N.size()
            171
            sage: N.delete_edge( 1, 2 )
            sage: N.delete_edge( (3, 4) )
            sage: N.delete_edges( [ (5, 6), (7, 8) ] )
            sage: N.size()
            167
            sage: N.delete_edge( 9, 10, 'label' )
            sage: N.delete_edge( (11, 12, 'label') )
            sage: N.delete_edges( [ (13, 14, 'label') ] )
            sage: N.size()
            167
            sage: N.has_edge( (11, 12) )
            True

        However, CGraph backends handle things properly::

            sage: G.delete_edge( 9, 10, 'label' )
            sage: G.delete_edge( (11, 12, 'label') )
            sage: G.delete_edges( [ (13, 14, 'label') ] )
            sage: G.size()
            167

        ::

            sage: C = graphs.CompleteGraph(19).to_directed(sparse=True)
            sage: C.size()
            342
            sage: C.delete_edge( 1, 2 )
            sage: C.delete_edge( (3, 4) )
            sage: C.delete_edges( [ (5, 6), (7, 8) ] )

            sage: D = graphs.CompleteGraph(19).to_directed(sparse=True, implementation='networkx')
            sage: D.size()
            342
            sage: D.delete_edge( 1, 2 )
            sage: D.delete_edge( (3, 4) )
            sage: D.delete_edges( [ (5, 6), (7, 8) ] )
            sage: D.delete_edge( 9, 10, 'label' )
            sage: D.delete_edge( (11, 12, 'label') )
            sage: D.delete_edges( [ (13, 14, 'label') ] )
            sage: D.size()
            338
            sage: D.has_edge( (11, 12) )
            True

        ::

            sage: C.delete_edge( 9, 10, 'label' )
            sage: C.delete_edge( (11, 12, 'label') )
            sage: C.delete_edges( [ (13, 14, 'label') ] )
            sage: C.size() # correct!
            338
            sage: C.has_edge( (11, 12) ) # correct!
            True

        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    u, v = u
                    label = None
        self._backend.del_edge(u, v, label, self._directed)

    def delete_edges(self, edges):
        """
        Delete edges from an iterable container.

        EXAMPLES::

            sage: K12 = graphs.CompleteGraph(12)
            sage: K4 = graphs.CompleteGraph(4)
            sage: K12.size()
            66
            sage: K12.delete_edges(K4.edge_iterator())
            sage: K12.size()
            60

        ::

            sage: K12 = graphs.CompleteGraph(12).to_directed()
            sage: K4 = graphs.CompleteGraph(4).to_directed()
            sage: K12.size()
            132
            sage: K12.delete_edges(K4.edge_iterator())
            sage: K12.size()
            120
        """
        for e in edges:
            self.delete_edge(e)

    def delete_multiedge(self, u, v):
        """
        Deletes all edges from u and v.

        EXAMPLES::

            sage: G = Graph(multiedges=True,sparse=True)
            sage: G.add_edges([(0,1), (0,1), (0,1), (1,2), (2,3)])
            sage: G.edges()
            [(0, 1, None), (0, 1, None), (0, 1, None), (1, 2, None), (2, 3, None)]
            sage: G.delete_multiedge( 0, 1 )
            sage: G.edges()
            [(1, 2, None), (2, 3, None)]

        ::

            sage: D = DiGraph(multiedges=True,sparse=True)
            sage: D.add_edges([(0,1,1), (0,1,2), (0,1,3), (1,0,None), (1,2,None), (2,3,None)])
            sage: D.edges()
            [(0, 1, 1), (0, 1, 2), (0, 1, 3), (1, 0, None), (1, 2, None), (2, 3, None)]
            sage: D.delete_multiedge( 0, 1 )
            sage: D.edges()
            [(1, 0, None), (1, 2, None), (2, 3, None)]
        """
        if self.allows_multiple_edges():
            for l in self.edge_label(u, v):
                self.delete_edge(u, v, l)
        else:
            self.delete_edge(u, v)

    def set_edge_label(self, u, v, l):
        """
        Set the edge label of a given edge.

        .. note::

           There can be only one edge from u to v for this to make
           sense. Otherwise, an error is raised.

        INPUT:


        -  ``u, v`` - the vertices (and direction if digraph)
           of the edge

        -  ``l`` - the new label


        EXAMPLES::

            sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] }, sparse=True)
            sage: SD.set_edge_label(1, 18, 'discrete')
            sage: SD.set_edge_label(4, 7, 'discrete')
            sage: SD.set_edge_label(2, 5, 'h = 0')
            sage: SD.set_edge_label(7, 18, 'h = 0')
            sage: SD.set_edge_label(7, 10, 'aut')
            sage: SD.set_edge_label(8, 10, 'aut')
            sage: SD.set_edge_label(8, 9, 'label')
            sage: SD.set_edge_label(8, 6, 'no label')
            sage: SD.set_edge_label(13, 17, 'k > h')
            sage: SD.set_edge_label(13, 14, 'k = h')
            sage: SD.set_edge_label(17, 15, 'v_k finite')
            sage: SD.set_edge_label(14, 15, 'v_k m.c.r.')
            sage: posn = {1:[ 3,-3],  2:[0,2],  3:[0, 13],  4:[3,9],  5:[3,3],  6:[16, 13], 7:[6,1],  8:[6,6],  9:[6,11], 10:[9,1], 11:[10,6], 12:[13,6], 13:[16,2], 14:[10,-6], 15:[0,-10], 16:[14,-6], 17:[16,-10], 18:[6,-4]}
            sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True).show() # long time

        ::

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.edges()
                [(0, 1, '(0,1)'),
                 (0, 5, '(0,5)'),
                 (0, 13, '(0,13)'),
                 ...
                 (11, 12, '(11,12)'),
                 (12, 13, '(12,13)')]

        ::

            sage: g = Graph({0: [0,1,1,2]}, loops=True, multiedges=True, sparse=True)
            sage: g.set_edge_label(0,0,'test')
            sage: g.edges()
            [(0, 0, 'test'), (0, 1, None), (0, 1, None), (0, 2, None)]
            sage: g.add_edge(0,0,'test2')
            sage: g.set_edge_label(0,0,'test3')
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot set edge label, since there are multiple edges from 0 to 0.

        ::

            sage: dg = DiGraph({0 : [1], 1 : [0]}, sparse=True)
            sage: dg.set_edge_label(0,1,5)
            sage: dg.set_edge_label(1,0,9)
            sage: dg.outgoing_edges(1)
            [(1, 0, 9)]
            sage: dg.incoming_edges(1)
            [(0, 1, 5)]
            sage: dg.outgoing_edges(0)
            [(0, 1, 5)]
            sage: dg.incoming_edges(0)
            [(1, 0, 9)]

        ::

            sage: G = Graph({0:{1:1}}, sparse=True)
            sage: G.num_edges()
            1
            sage: G.set_edge_label(0,1,1)
            sage: G.num_edges()
            1
        """
        if self.allows_multiple_edges():
            if len(self.edge_label(u, v)) > 1:
                raise RuntimeError("Cannot set edge label, since there are multiple edges from %s to %s."%(u,v))
        self._backend.set_edge_label(u, v, l, self._directed)

    def has_edge(self, u, v=None, label=None):
        r"""
        Returns True if (u, v) is an edge, False otherwise.

        INPUT: The following forms are accepted by NetworkX:

        - G.has_edge( 1, 2 )
        - G.has_edge( (1, 2) )
        - G.has_edge( 1, 2, 'label' )
        - G.has_edge( (1, 2, 'label') )

        EXAMPLES::

            sage: graphs.EmptyGraph().has_edge(9,2)
            False
            sage: DiGraph().has_edge(9,2)
            False
            sage: G = Graph(sparse=True)
            sage: G.add_edge(0,1,"label")
            sage: G.has_edge(0,1,"different label")
            False
            sage: G.has_edge(0,1,"label")
            True
        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    u, v = u
                    label = None
        return self._backend.has_edge(u, v, label)

    def edges(self, labels=True, sort=True, key=None):
        r"""
        Return a list of edges.

        Each edge is a triple (u,v,l) where u and v are vertices and l is a
        label. If the parameter ``labels`` is False then a list of couple (u,v)
        is returned where u and v are vertices.

        INPUT:

        - ``labels`` - default: ``True`` - if ``False``, each
          edge is simply a pair (u,v) of vertices.

        - ``sort`` - default: ``True`` - if ``True``, edges are
          sorted according to the default ordering.

        - ``key`` - default: ``None`` - a function  takes an edge
          (a pair or a triple, according to the ``labels`` keyword)
          as its one argument and returns a value that can be used
          for comparisons in the sorting algorithm.

        OUTPUT: A list of tuples. It is safe to change the returned list.

        .. warning::

            Since any object may be a vertex, there is no guarantee
            that any two vertices will be comparable, and thus no
            guarantee how two edges may compare.  With default
            objects for vertices (all integers), or when all the
            vertices are of the same simple type, then there should
            not be a problem with how the vertices will be sorted.
            However, if you need to guarantee a total order for
            the sorting of the edges, use the ``key`` argument,
            as illustrated in the examples below.

        EXAMPLES::

            sage: graphs.DodecahedralGraph().edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 2, None), (1, 8, None), (2, 3, None), (2, 6, None), (3, 4, None), (3, 19, None), (4, 5, None), (4, 17, None), (5, 6, None), (5, 15, None), (6, 7, None), (7, 8, None), (7, 14, None), (8, 9, None), (9, 10, None), (9, 13, None), (10, 11, None), (11, 12, None), (11, 18, None), (12, 13, None), (12, 16, None), (13, 14, None), (14, 15, None), (15, 16, None), (16, 17, None), (17, 18, None), (18, 19, None)]

        ::

            sage: graphs.DodecahedralGraph().edges(labels=False)
            [(0, 1), (0, 10), (0, 19), (1, 2), (1, 8), (2, 3), (2, 6), (3, 4), (3, 19), (4, 5), (4, 17), (5, 6), (5, 15), (6, 7), (7, 8), (7, 14), (8, 9), (9, 10), (9, 13), (10, 11), (11, 12), (11, 18), (12, 13), (12, 16), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19)]

        ::

            sage: D = graphs.DodecahedralGraph().to_directed()
            sage: D.edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 0, None), (1, 2, None), (1, 8, None), (2, 1, None), (2, 3, None), (2, 6, None), (3, 2, None), (3, 4, None), (3, 19, None), (4, 3, None), (4, 5, None), (4, 17, None), (5, 4, None), (5, 6, None), (5, 15, None), (6, 2, None), (6, 5, None), (6, 7, None), (7, 6, None), (7, 8, None), (7, 14, None), (8, 1, None), (8, 7, None), (8, 9, None), (9, 8, None), (9, 10, None), (9, 13, None), (10, 0, None), (10, 9, None), (10, 11, None), (11, 10, None), (11, 12, None), (11, 18, None), (12, 11, None), (12, 13, None), (12, 16, None), (13, 9, None), (13, 12, None), (13, 14, None), (14, 7, None), (14, 13, None), (14, 15, None), (15, 5, None), (15, 14, None), (15, 16, None), (16, 12, None), (16, 15, None), (16, 17, None), (17, 4, None), (17, 16, None), (17, 18, None), (18, 11, None), (18, 17, None), (18, 19, None), (19, 0, None), (19, 3, None), (19, 18, None)]
            sage: D.edges(labels = False)
            [(0, 1), (0, 10), (0, 19), (1, 0), (1, 2), (1, 8), (2, 1), (2, 3), (2, 6), (3, 2), (3, 4), (3, 19), (4, 3), (4, 5), (4, 17), (5, 4), (5, 6), (5, 15), (6, 2), (6, 5), (6, 7), (7, 6), (7, 8), (7, 14), (8, 1), (8, 7), (8, 9), (9, 8), (9, 10), (9, 13), (10, 0), (10, 9), (10, 11), (11, 10), (11, 12), (11, 18), (12, 11), (12, 13), (12, 16), (13, 9), (13, 12), (13, 14), (14, 7), (14, 13), (14, 15), (15, 5), (15, 14), (15, 16), (16, 12), (16, 15), (16, 17), (17, 4), (17, 16), (17, 18), (18, 11), (18, 17), (18, 19), (19, 0), (19, 3), (19, 18)]

        The default is to sort the returned list in the default fashion, as in the above examples.
        this can be overridden by specifying a key function. This first example just ignores
        the labels in the third component of the triple.  ::

            sage: G=graphs.CycleGraph(5)
            sage: G.edges(key = lambda x: (x[1],-x[0]))
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 4, None), (0, 4, None)]

        We set the labels to characters and then perform a default sort
        followed by a sort according to the labels. ::

            sage: G=graphs.CycleGraph(5)
            sage: for e in G.edges():
            ...     G.set_edge_label(e[0], e[1], chr(ord('A')+e[0]+5*e[1]))
            sage: G.edges(sort=True)
            [(0, 1, 'F'), (0, 4, 'U'), (1, 2, 'L'), (2, 3, 'R'), (3, 4, 'X')]
            sage: G.edges(key=lambda x: x[2])
            [(0, 1, 'F'), (1, 2, 'L'), (2, 3, 'R'), (0, 4, 'U'), (3, 4, 'X')]

        TESTS:

        It is an error to turn off sorting while providing a key function for sorting. ::

            sage: P=graphs.PetersenGraph()
            sage: P.edges(sort=False, key=lambda x: x)
            Traceback (most recent call last):
            ...
            ValueError: sort keyword is False, yet a key function is given
        """
        if not(sort) and key:
            raise ValueError('sort keyword is False, yet a key function is given')
        L = list(self.edge_iterator(labels=labels))
        if sort:
            L.sort(key=key)
        return L

    def edge_boundary(self, vertices1, vertices2=None, labels=True, sort=True):
        """
        Returns a list of edges `(u,v,l)` with `u` in ``vertices1``
        and `v` in ``vertices2``. If ``vertices2`` is ``None``, then
        it is set to the complement of ``vertices1``.

        In a digraph, the external boundary of a vertex `v` are those
        vertices `u` with an arc `(v, u)`.

        INPUT:


        -  ``labels`` - if ``False``, each edge is a tuple `(u,v)` of
           vertices.


        EXAMPLES::

            sage: K = graphs.CompleteBipartiteGraph(9,3)
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            27

        Note that the edge boundary preserves direction::

            sage: K = graphs.CompleteBipartiteGraph(9,3).to_directed()
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            54

        ::

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: D.edge_boundary([0])
            [(0, 1, None), (0, 2, None)]
            sage: D.edge_boundary([0], labels=False)
            [(0, 1), (0, 2)]

        TESTS::

            sage: G = graphs.DiamondGraph()
            sage: G.edge_boundary([0,1])
            [(0, 2, None), (1, 2, None), (1, 3, None)]
            sage: G.edge_boundary([0], [0])
            []
            sage: G.edge_boundary([2], [0])
            [(0, 2, None)]
        """
        vertices1 = set([v for v in vertices1 if v in self])
        if self._directed:
            if vertices2 is not None:
                vertices2 = set([v for v in vertices2 if v in self])
                output = [e for e in self.outgoing_edge_iterator(vertices1,labels=labels)
                            if e[1] in vertices2]
            else:
                output = [e for e in self.outgoing_edge_iterator(vertices1,labels=labels)
                            if e[1] not in vertices1]
        else:
            if vertices2 is not None:
                vertices2 = set([v for v in vertices2 if v in self])
                output = [e for e in self.edge_iterator(vertices1,labels=labels)
                            if (e[0] in vertices1 and e[1] in vertices2) or
                            (e[1] in vertices1 and e[0] in vertices2)]
            else:
                output = [e for e in self.edge_iterator(vertices1,labels=labels)
                            if e[1] not in vertices1 or e[0] not in vertices1]
        if sort:
            output.sort()
        return output

    def edge_iterator(self, vertices=None, labels=True, ignore_direction=False):
        """
        Returns an iterator over edges.

        The iterator returned is over the edges incident with any vertex given
        in the parameter ``vertices``. If the graph is directed, iterates over
        edges going out only. If vertices is None, then returns an iterator over
        all edges. If self is directed, returns outgoing edges only.

        INPUT:

        - ``vertices`` - (default: None) a vertex, a list of vertices or None

        -  ``labels`` - if False, each edge is a tuple (u,v) of
           vertices.

        -  ``ignore_direction`` - bool (default: False) - only applies
           to directed graphs. If True, searches across edges in either
           direction.


        EXAMPLES::

            sage: for i in graphs.PetersenGraph().edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 4, None)
            (0, 5, None)
            sage: D = DiGraph( { 0 : [1,2], 1: [0] } )
            sage: for i in D.edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 2, None)

        ::

            sage: G = graphs.TetrahedralGraph()
            sage: list(G.edge_iterator(labels=False))
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

        ::

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.edge_iterator(0))
            []
            sage: list(D.edge_iterator(0, ignore_direction=True))
            [(1, 0, None), (2, 0, None)]
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]

        if ignore_direction and self._directed:
            from itertools import chain
            return chain(self._backend.iterator_out_edges(vertices, labels),
                         self._backend.iterator_in_edges(vertices, labels))
        elif self._directed:
            return self._backend.iterator_out_edges(vertices, labels)
        else:
            return self._backend.iterator_edges(vertices, labels)

    def edges_incident(self, vertices=None, labels=True, sort=True):
        """
        Returns incident edges to some vertices.

        If ``vertices` is a vertex, then it returns the list of edges incident to
        that vertex. If ``vertices`` is a list of vertices then it returns the
        list of all edges adjacent to those vertices. If ``vertices``
        is None, returns a list of all edges in graph. For digraphs, only
        lists outward edges.

        INPUT:

        - ``vertices`` - object (default: None) - a vertex, a list of vertices
          or None.

        -  ``labels`` - bool (default: True) - if False, each edge is a tuple
           (u,v) of vertices.

        - ``sort`` - bool (default: True) - if True the returned list is sorted.


        EXAMPLES::

            sage: graphs.PetersenGraph().edges_incident([0,9], labels=False)
            [(0, 1), (0, 4), (0, 5), (4, 9), (6, 9), (7, 9)]
            sage: D = DiGraph({0:[1]})
            sage: D.edges_incident([0])
            [(0, 1, None)]
            sage: D.edges_incident([1])
            []

        TESTS::

            sage: G = Graph({0:[0]}, loops=True)  # ticket 9581
            sage: G.edges_incident(0)
            [(0, 0, None)]
        """
        if vertices in self:
            vertices = [vertices]

        if sort:
            return sorted(self.edge_iterator(vertices=vertices,labels=labels))
        return list(self.edge_iterator(vertices=vertices,labels=labels))

    def edge_label(self, u, v=None):
        """
        Returns the label of an edge. Note that if the graph allows
        multiple edges, then a list of labels on the edge is returned.

        EXAMPLES::

            sage: G = Graph({0 : {1 : 'edgelabel'}}, sparse=True)
            sage: G.edges(labels=False)
            [(0, 1)]
            sage: G.edge_label( 0, 1 )
            'edgelabel'
            sage: D = DiGraph({0 : {1 : 'edgelabel'}}, sparse=True)
            sage: D.edges(labels=False)
            [(0, 1)]
            sage: D.edge_label( 0, 1 )
            'edgelabel'

        ::

            sage: G = Graph(multiedges=True, sparse=True)
            sage: [G.add_edge(0,1,i) for i in range(1,6)]
            [None, None, None, None, None]
            sage: sorted(G.edge_label(0,1))
            [1, 2, 3, 4, 5]

        TESTS::

            sage: G = Graph()
            sage: G.add_edge(0,1,[7])
            sage: G.add_edge(0,2,[7])
            sage: G.edge_label(0,1)[0] += 1
            sage: G.edges()
            [(0, 1, [8]), (0, 2, [7])]

        """
        return self._backend.get_edge_label(u,v)

    def edge_labels(self):
        """
        Returns a list of edge labels.

        EXAMPLES::

            sage: G = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, sparse=True)
            sage: G.edge_labels()
            ['x', 'z', 'a', 'out']
            sage: G = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, sparse=True)
            sage: G.edge_labels()
            ['x', 'z', 'a', 'out']
        """
        labels = []
        for u,v,l in self.edges():
            labels.append(l)
        return labels

    def remove_multiple_edges(self):
        """
        Removes all multiple edges, retaining one edge for each.

        EXAMPLES::

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges( [ (0,1), (0,1), (0,1), (0,1), (1,2) ] )
            sage: G.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]

        ::

            sage: G.remove_multiple_edges()
            sage: G.edges(labels=False)
            [(0, 1), (1, 2)]

        ::

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: D.add_edges( [ (0,1,1), (0,1,2), (0,1,3), (0,1,4), (1,2,None) ] )
            sage: D.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            sage: D.remove_multiple_edges()
            sage: D.edges(labels=False)
            [(0, 1), (1, 2)]
        """
        if self.allows_multiple_edges():
            if self._directed:
                for v in self:
                    for u in self.neighbor_in_iterator(v):
                        edges = self.edge_boundary([u], [v])
                        if len(edges) > 1:
                            self.delete_edges(edges[1:])
            else:
                for v in self:
                    for u in self.neighbor_iterator(v):
                        edges = self.edge_boundary([v], [u])
                        if len(edges) > 1:
                            self.delete_edges(edges[1:])

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes
        all loops.

        EXAMPLE

        ::

            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.remove_loops()
            sage: G.edges(labels=False)
            [(2, 3)]
            sage: G.allows_loops()
            True
            sage: G.has_loops()
            False

            sage: D = DiGraph(4, loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.remove_loops()
            sage: D.edges(labels=False)
            [(2, 3)]
            sage: D.allows_loops()
            True
            sage: D.has_loops()
            False
        """
        if vertices is None:
            vertices = self
        for v in vertices:
            if self.has_edge(v,v):
                self.delete_multiedge(v,v)

    def loop_edges(self):
        """
        Returns a list of all loops in the graph.

        EXAMPLES::

            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

        ::

            sage: D = DiGraph(4, loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

        ::

            sage: G = Graph(4, loops=True, multiedges=True, sparse=True)
            sage: G.add_edges([(i,i) for i in range(4)])
            sage: G.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]
        """
        if self.allows_multiple_edges():
            return [(v,v,l) for v in self.loop_vertices() for l in self.edge_label(v,v)]
        else:
            return [(v,v,self.edge_label(v,v)) for v in self.loop_vertices()]

    def number_of_loops(self):
        """
        Returns the number of edges that are loops.

        EXAMPLES::

            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.number_of_loops()
            4

        ::

            sage: D = DiGraph(4, loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.number_of_loops()
            4
        """
        return len(self.loop_edges())

    ### Modifications

    def clear(self):
        """
        Empties the graph of vertices and edges and removes name, boundary,
        associated objects, and position information.

        EXAMPLES::

            sage: G=graphs.CycleGraph(4); G.set_vertices({0:'vertex0'})
            sage: G.order(); G.size()
            4
            4
            sage: len(G._pos)
            4
            sage: G.name()
            'Cycle graph'
            sage: G.get_vertex(0)
            'vertex0'
            sage: H = G.copy(implementation='c_graph', sparse=True)
            sage: H.clear()
            sage: H.order(); H.size()
            0
            0
            sage: len(H._pos)
            0
            sage: H.name()
            ''
            sage: H.get_vertex(0)
            sage: H = G.copy(implementation='c_graph', sparse=False)
            sage: H.clear()
            sage: H.order(); H.size()
            0
            0
            sage: len(H._pos)
            0
            sage: H.name()
            ''
            sage: H.get_vertex(0)
            sage: H = G.copy(implementation='networkx')
            sage: H.clear()
            sage: H.order(); H.size()
            0
            0
            sage: len(H._pos)
            0
            sage: H.name()
            ''
            sage: H.get_vertex(0)
        """
        self.name('')
        self.delete_vertices(self.vertices())

    ### Degree functions

    def degree(self, vertices=None, labels=False):
        """
        Gives the degree (in + out for digraphs) of a vertex or of
        vertices.

        INPUT:


        -  ``vertices`` - If vertices is a single vertex,
           returns the number of neighbors of vertex. If vertices is an
           iterable container of vertices, returns a list of degrees. If
           vertices is None, same as listing all vertices.

        -  ``labels`` - see OUTPUT


        OUTPUT: Single vertex- an integer. Multiple vertices- a list of
        integers. If labels is True, then returns a dictionary mapping each
        vertex to its degree.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.degree(5)
            3

        ::

            sage: K = graphs.CompleteGraph(9)
            sage: K.degree()
            [8, 8, 8, 8, 8, 8, 8, 8, 8]

        ::

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.degree(vertices = [0,1,2], labels=True)
            {0: 5, 1: 4, 2: 3}
            sage: D.degree()
            [5, 4, 3, 3, 3, 2]
        """
        if labels:
            return dict(self.degree_iterator(vertices,labels))
        elif vertices in self and not labels:
            return self.degree_iterator(vertices,labels).next()
        else:
            return list(self.degree_iterator(vertices,labels))

    def average_degree(self):
        r"""
        Returns the average degree of the graph.

        The average degree of a graph `G=(V,E)` is equal to
        ``\frac {2|E|}{|V|}``.

        EXAMPLES:

        The average degree of a regular graph is equal to the
        degree of any vertex::

            sage: g = graphs.CompleteGraph(5)
            sage: g.average_degree() == 4
            True

        The average degree of a tree is always strictly less than
        `2`::

           sage: g = graphs.RandomGNP(20,.5)
           sage: tree = Graph()
           sage: tree.add_edges(g.min_spanning_tree())
           sage: tree.average_degree() < 2
           True

        For any graph, it is equal to ``\frac {2|E|}{|V|}``::

            sage: g = graphs.RandomGNP(50,.8)
            sage: g.average_degree() == 2*g.size()/g.order()
            True
        """

        return 2*Integer(self.size())/Integer(self.order())

    def degree_histogram(self):
        """
        Returns a list, whose ith entry is the frequency of degree i.

        EXAMPLES::

            sage: G = graphs.Grid2dGraph(9,12)
            sage: G.degree_histogram()
            [0, 0, 4, 34, 70]

        ::

            sage: G = graphs.Grid2dGraph(9,12).to_directed()
            sage: G.degree_histogram()
            [0, 0, 0, 0, 4, 0, 34, 0, 70]
        """
        degree_sequence = self.degree()
        dmax = max(degree_sequence) + 1
        frequency = [0]*dmax
        for d in degree_sequence:
            frequency[d] += 1
        return frequency

    def degree_iterator(self, vertices=None, labels=False):
        """
        Returns an iterator over the degrees of the (di)graph.

        In the case of a digraph, the degree is defined as the sum of the
        in-degree and the out-degree, i.e. the total number of edges incident to
        a given vertex.

        INPUT:

        - ``labels`` (boolean) -- if set to ``False`` (default) the method
          returns an iterator over degrees. Otherwise it returns an iterator
          over tuples (vertex, degree).

        -  ``vertices`` - if specified, restrict to this
           subset.


        EXAMPLES::

            sage: G = graphs.Grid2dGraph(3,4)
            sage: for i in G.degree_iterator():
            ...    print i
            3
            4
            2
            ...
            2
            4
            sage: for i in G.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 4)
            ((0, 0), 2)
            ...
            ((0, 3), 2)
            ((1, 1), 4)

        ::

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.degree_iterator():
            ...    print i
            6
            6
            ...
            4
            6
            sage: for i in D.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 6)
            ((1, 2), 6)
            ...
            ((0, 3), 4)
            ((1, 1), 6)
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        if labels:
            filter_ = lambda v, self: (v, self._backend.degree(v, self._directed))
        else:
            filter_ = lambda v, self: self._backend.degree(v, self._directed)
        for v in vertices:
            yield filter_(v, self)

    def degree_sequence(self):
        r"""
        Return the degree sequence of this (di)graph.

        EXAMPLES:

        The degree sequence of an undirected graph::

            sage: g = Graph({1: [2, 5], 2: [1, 5, 3, 4], 3: [2, 5], 4: [3], 5: [2, 3]})
            sage: g.degree_sequence()
            [4, 3, 3, 2, 2]

        The degree sequence of a digraph::

            sage: g = DiGraph({1: [2, 5, 6], 2: [3, 6], 3: [4, 6], 4: [6], 5: [4, 6]})
            sage: g.degree_sequence()
            [5, 3, 3, 3, 3, 3]

        Degree sequences of some common graphs::

            sage: graphs.PetersenGraph().degree_sequence()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: graphs.HouseGraph().degree_sequence()
            [3, 3, 2, 2, 2]
            sage: graphs.FlowerSnark().degree_sequence()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        """
        return sorted(self.degree_iterator(), reverse=True)

    def is_regular(self, k = None):
        """
        Return ``True`` if this graph is (`k`-)regular.

        INPUT:

        - ``k`` (default: ``None``) - the degree of regularity to
          check for

        EXAMPLES::

            sage: G = graphs.HoffmanSingletonGraph()
            sage: G.is_regular()
            True
            sage: G.is_regular(9)
            False

        So the Hoffman-Singleton graph is regular, but not
        9-regular.  In fact, we can now find the degree easily as
        follows::

            sage: G.degree_iterator().next()
            7

        The house graph is not regular::

            sage: graphs.HouseGraph().is_regular()
            False

        A graph without vertices is `k`-regular for every `k`::

            sage: Graph().is_regular()
            True
        """
        if self.order() == 0:
            return True

        deg_it = self.degree_iterator()
        if k is None:
            k = deg_it.next()

        for d in deg_it:
            if d != k:
                return False

        return True

    ### Substructures

    def subgraph(self, vertices=None, edges=None, inplace=False,
                       vertex_property=None, edge_property=None, algorithm=None):
        """
        Returns the subgraph containing the given vertices and edges.

        If either vertices or edges are not specified, they are assumed to be
        all vertices or edges. If edges are not specified, returns the subgraph
        induced by the vertices.

        INPUT:


        -  ``inplace`` - Using inplace is True will simply
           delete the extra vertices and edges from the current graph. This
           will modify the graph.

        -  ``vertices`` - Vertices can be a single vertex or an
           iterable container of vertices, e.g. a list, set, graph, file or
           numeric array. If not passed, defaults to the entire graph.

        -  ``edges`` - As with vertices, edges can be a single
           edge or an iterable container of edges (e.g., a list, set, file,
           numeric array, etc.). If not edges are not specified, then all
           edges are assumed and the returned graph is an induced subgraph. In
           the case of multiple edges, specifying an edge as (u,v) means to
           keep all edges (u,v), regardless of the label.

        -  ``vertex_property`` - If specified, this is
           expected to be a function on vertices, which is intersected with
           the vertices specified, if any are.

        -  ``edge_property`` - If specified, this is expected
           to be a function on edges, which is intersected with the edges
           specified, if any are.

        - ``algorithm`` - If ``algorithm=delete`` or ``inplace=True``,
          then the graph is constructed by deleting edges and
          vertices.  If ``add``, then the graph is constructed by
          building a new graph from the appropriate vertices and
          edges.  If not specified, then the algorithm is chosen based
          on the number of vertices in the subgraph.


        EXAMPLES::

            sage: G = graphs.CompleteGraph(9)
            sage: H = G.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: J = G.subgraph(edges=[(0,1)])
            sage: J.edges(labels=False)
            [(0, 1)]
            sage: J.vertices()==G.vertices()
            True
            sage: G.subgraph([0,1,2], inplace=True); G
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G.subgraph()==G
            True

        ::

            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: H = D.subgraph(edges=[(0,1), (0,2)])
            sage: H.edges(labels=False)
            [(0, 1), (0, 2)]
            sage: H.vertices()==D.vertices()
            True
            sage: D
            Complete graph: Digraph on 9 vertices
            sage: D.subgraph([0,1,2], inplace=True); D
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: D.subgraph()==D
            True

        A more complicated example involving multiple edges and labels.

        ::

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: G.subgraph(edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 1, 'c'), (0, 2, 'd')]
            sage: J = G.subgraph(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: J.edges()
            [(0, 1, 'a')]
            sage: J.vertices()
            [0, 1]
            sage: G.subgraph(vertices=G.vertices())==G
            True

        ::

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: D.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: D.subgraph(edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 2, 'd')]
            sage: H = D.subgraph(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: H.edges()
            [(0, 1, 'a')]
            sage: H.vertices()
            [0, 1]

        Using the property arguments::

            sage: P = graphs.PetersenGraph()
            sage: S = P.subgraph(vertex_property = lambda v : v%2 == 0)
            sage: S.vertices()
            [0, 2, 4, 6, 8]

        ::

            sage: C = graphs.CubeGraph(2)
            sage: S = C.subgraph(edge_property=(lambda e: e[0][0] == e[1][0]))
            sage: C.edges()
            [('00', '01', None), ('00', '10', None), ('01', '11', None), ('10', '11', None)]
            sage: S.edges()
            [('00', '01', None), ('10', '11', None)]


        The algorithm is not specified, then a reasonable choice is made for speed.

        ::

            sage: g=graphs.PathGraph(1000)
            sage: g.subgraph(range(10)) # uses the 'add' algorithm
            Subgraph of (Path Graph): Graph on 10 vertices



        TESTS: The appropriate properties are preserved.

        ::

            sage: g = graphs.PathGraph(10)
            sage: g.is_planar(set_embedding=True)
            True
            sage: g.set_vertices(dict((v, 'v%d'%v) for v in g.vertices()))
            sage: h = g.subgraph([3..5])
            sage: h.get_pos().keys()
            [3, 4, 5]
            sage: h.get_vertices()
            {3: 'v3', 4: 'v4', 5: 'v5'}
        """
        if vertices is None:
            vertices=self.vertices()
        elif vertices in self:
            vertices=[vertices]
        else:
            vertices=list(vertices)

        if vertex_property is not None:
            vertices = [v for v in vertices if vertex_property(v)]

        if algorithm is not None and algorithm not in ("delete", "add"):
            raise ValueError('algorithm should be None, "delete", or "add"')

        if inplace or len(vertices)>0.05*self.order() or algorithm=="delete":
            return self._subgraph_by_deleting(vertices=vertices, edges=edges,
                                              inplace=inplace,
                                              edge_property=edge_property)
        else:
            return self._subgraph_by_adding(vertices=vertices, edges=edges,
                                            edge_property=edge_property)

    def _subgraph_by_adding(self, vertices=None, edges=None, edge_property=None):
        """
        Returns the subgraph containing the given vertices and edges.
        The edges also satisfy the edge_property, if it is not None.
        The subgraph is created by creating a new empty graph and
        adding the necessary vertices, edges, and other properties.

        INPUT:

        -  ``vertices`` - Vertices is a list of vertices

        - ``edges`` - Edges can be a single edge or an iterable
           container of edges (e.g., a list, set, file, numeric array,
           etc.). If not edges are not specified, then all edges are
           assumed and the returned graph is an induced subgraph. In
           the case of multiple edges, specifying an edge as (u,v)
           means to keep all edges (u,v), regardless of the label.

        -  ``edge_property`` - If specified, this is expected
           to be a function on edges, which is intersected with the edges
           specified, if any are.


        EXAMPLES::

            sage: G = graphs.CompleteGraph(9)
            sage: H = G._subgraph_by_adding([0,1,2]); H
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: J = G._subgraph_by_adding(vertices=G.vertices(), edges=[(0,1)])
            sage: J.edges(labels=False)
            [(0, 1)]
            sage: J.vertices()==G.vertices()
            True
            sage: G._subgraph_by_adding(vertices=G.vertices())==G
            True

        ::

            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D._subgraph_by_adding([0,1,2]); H
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: H = D._subgraph_by_adding(vertices=D.vertices(), edges=[(0,1), (0,2)])
            sage: H.edges(labels=False)
            [(0, 1), (0, 2)]
            sage: H.vertices()==D.vertices()
            True
            sage: D
            Complete graph: Digraph on 9 vertices
            sage: D._subgraph_by_adding(D.vertices())==D
            True

        A more complicated example involving multiple edges and labels.

        ::

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: G._subgraph_by_adding(G.vertices(), edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 1, 'c'), (0, 2, 'd')]
            sage: J = G._subgraph_by_adding(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: J.edges()
            [(0, 1, 'a')]
            sage: J.vertices()
            [0, 1]
            sage: G._subgraph_by_adding(vertices=G.vertices())==G
            True

        ::

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: D.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: D._subgraph_by_adding(vertices=D.vertices(), edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 2, 'd')]
            sage: H = D._subgraph_by_adding(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: H.edges()
            [(0, 1, 'a')]
            sage: H.vertices()
            [0, 1]

        Using the property arguments::

            sage: C = graphs.CubeGraph(2)
            sage: S = C._subgraph_by_adding(vertices=C.vertices(), edge_property=(lambda e: e[0][0] == e[1][0]))
            sage: C.edges()
            [('00', '01', None), ('00', '10', None), ('01', '11', None), ('10', '11', None)]
            sage: S.edges()
            [('00', '01', None), ('10', '11', None)]

        TESTS: Properties of the graph are preserved.

        ::

            sage: g = graphs.PathGraph(10)
            sage: g.is_planar(set_embedding=True)
            True
            sage: g.set_vertices(dict((v, 'v%d'%v) for v in g.vertices()))
            sage: h = g._subgraph_by_adding([3..5])
            sage: h.get_pos().keys()
            [3, 4, 5]
            sage: h.get_vertices()
            {3: 'v3', 4: 'v4', 5: 'v5'}
        """
        G = self.__class__(weighted=self._weighted, loops=self.allows_loops(),
                           multiedges= self.allows_multiple_edges())
        G.name("Subgraph of (%s)"%self.name())
        G.add_vertices(vertices)
        if edges is not None:
            if G._directed:
                edges_graph = (e for e in self.edge_iterator(vertices) if e[1] in vertices)
                edges_to_keep_labeled = [e for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [e for e in edges if len(e)==2]
            else:
                edges_graph = (sorted(e[0:2])+[e[2]] for e in self.edge_iterator(vertices) if e[0] in vertices and e[1] in vertices)
                edges_to_keep_labeled = [sorted(e[0:2])+[e[2]] for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [sorted(e) for e in edges if len(e)==2]
            edges_to_keep = [tuple(e) for e in edges_graph if e in edges_to_keep_labeled
                             or e[0:2] in edges_to_keep_unlabeled]
        else:
            edges_to_keep=[e for e in self.edge_iterator(vertices) if e[0] in vertices and e[1] in vertices]

        if edge_property is not None:
            edges_to_keep = [e for e in edges_to_keep if edge_property(e)]
        G.add_edges(edges_to_keep)

        attributes_to_update = ('_pos', '_assoc')
        for attr in attributes_to_update:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                value = dict([(v, getattr(self, attr).get(v, None)) for v in G])
                setattr(G, attr,value)

        G._boundary = [v for v in self._boundary if v in G]

        return G

    def _subgraph_by_deleting(self, vertices=None, edges=None, inplace=False,
                              edge_property=None):
        """
        Returns the subgraph containing the given vertices and edges.
        The edges also satisfy the edge_property, if it is not None.
        The subgraph is created by creating deleting things that are
        not needed.

        INPUT:

        -  ``vertices`` - Vertices is a list of vertices

        - ``edges`` - Edges can be a single edge or an iterable
           container of edges (e.g., a list, set, file, numeric array,
           etc.). If not edges are not specified, then all edges are
           assumed and the returned graph is an induced subgraph. In
           the case of multiple edges, specifying an edge as (u,v)
           means to keep all edges (u,v), regardless of the label.

        -  ``edge_property`` - If specified, this is expected
           to be a function on edges, which is intersected with the edges
           specified, if any are.

        -  ``inplace`` - Using inplace is True will simply
           delete the extra vertices and edges from the current graph. This
           will modify the graph.


        EXAMPLES::

            sage: G = graphs.CompleteGraph(9)
            sage: H = G._subgraph_by_deleting([0,1,2]); H
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: J = G._subgraph_by_deleting(vertices=G.vertices(), edges=[(0,1)])
            sage: J.edges(labels=False)
            [(0, 1)]
            sage: J.vertices()==G.vertices()
            True
            sage: G._subgraph_by_deleting([0,1,2], inplace=True); G
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G._subgraph_by_deleting(vertices=G.vertices())==G
            True

        ::

            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D._subgraph_by_deleting([0,1,2]); H
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: H = D._subgraph_by_deleting(vertices=D.vertices(), edges=[(0,1), (0,2)])
            sage: H.edges(labels=False)
            [(0, 1), (0, 2)]
            sage: H.vertices()==D.vertices()
            True
            sage: D
            Complete graph: Digraph on 9 vertices
            sage: D._subgraph_by_deleting([0,1,2], inplace=True); D
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: D._subgraph_by_deleting(D.vertices())==D
            True

        A more complicated example involving multiple edges and labels.

        ::

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: G._subgraph_by_deleting(G.vertices(), edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 1, 'c'), (0, 2, 'd')]
            sage: J = G._subgraph_by_deleting(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: J.edges()
            [(0, 1, 'a')]
            sage: J.vertices()
            [0, 1]
            sage: G._subgraph_by_deleting(vertices=G.vertices())==G
            True

        ::

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: D.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: D._subgraph_by_deleting(vertices=D.vertices(), edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 2, 'd')]
            sage: H = D._subgraph_by_deleting(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: H.edges()
            [(0, 1, 'a')]
            sage: H.vertices()
            [0, 1]

        Using the property arguments::

            sage: C = graphs.CubeGraph(2)
            sage: S = C._subgraph_by_deleting(vertices=C.vertices(), edge_property=(lambda e: e[0][0] == e[1][0]))
            sage: C.edges()
            [('00', '01', None), ('00', '10', None), ('01', '11', None), ('10', '11', None)]
            sage: S.edges()
            [('00', '01', None), ('10', '11', None)]

        TESTS: Properties of the graph are preserved.

        ::

            sage: g = graphs.PathGraph(10)
            sage: g.is_planar(set_embedding=True)
            True
            sage: g.set_vertices(dict((v, 'v%d'%v) for v in g.vertices()))
            sage: h = g._subgraph_by_deleting([3..5])
            sage: h.get_pos().keys()
            [3, 4, 5]
            sage: h.get_vertices()
            {3: 'v3', 4: 'v4', 5: 'v5'}
        """
        if inplace:
            G = self
        else:
            G = self.copy()
        G.name("Subgraph of (%s)"%self.name())

        G.delete_vertices([v for v in G if v not in vertices])

        edges_to_delete=[]
        if edges is not None:
            if G._directed:
                edges_graph = G.edge_iterator()
                edges_to_keep_labeled = [e for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [e for e in edges if len(e)==2]
            else:
                edges_graph = [sorted(e[0:2])+[e[2]] for e in G.edge_iterator()]
                edges_to_keep_labeled = [sorted(e[0:2])+[e[2]] for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [sorted(e) for e in edges if len(e)==2]
            for e in edges_graph:
                if e not in edges_to_keep_labeled and e[0:2] not in edges_to_keep_unlabeled:
                    edges_to_delete.append(tuple(e))
        if edge_property is not None:
            for e in G.edge_iterator():
                if not edge_property(e):
                    # We might get duplicate edges, but this does
                    # handle the case of multiple edges.
                    edges_to_delete.append(e)

        G.delete_edges(edges_to_delete)
        if not inplace:
            return G

    def subgraph_search(self, G, induced=False):
        r"""
        Returns a copy of ``G`` in ``self``.

        INPUT:

        - ``G`` -- the graph whose copy we are looking for in ``self``.

        - ``induced`` -- boolean (default: ``False``). Whether or not to
          search for an induced copy of ``G`` in ``self``.

        OUTPUT:

        - If ``induced=False``, return a copy of ``G`` in this graph.
          Otherwise, return an induced copy of ``G`` in ``self``. If ``G``
          is the empty graph, return the empty graph since it is a subgraph
          of every graph. Now suppose ``G`` is not the empty graph. If there
          is no copy (induced or otherwise) of ``G`` in ``self``, we return
          ``None``.

        .. NOTE::

            This method also works on digraphs.

        .. SEEALSO::

            - :meth:`~GenericGraph.subgraph_search_count` -- Counts the number
              of copies of a graph `H` inside of a graph `G`

            - :meth:`~GenericGraph.subgraph_search_iterator` -- Iterate on the
              copies of a graph `H` inside of a graph `G`

        ALGORITHM:

        Brute-force search.

        EXAMPLES:

        The Petersen graph contains the path graph `P_5`::

             sage: g = graphs.PetersenGraph()
             sage: h1 = g.subgraph_search(graphs.PathGraph(5)); h1
             Subgraph of (Petersen graph): Graph on 5 vertices
             sage: h1.vertices(); h1.edges(labels=False)
             [0, 1, 2, 3, 4]
             [(0, 1), (1, 2), (2, 3), (3, 4)]
             sage: I1 = g.subgraph_search(graphs.PathGraph(5), induced=True); I1
             Subgraph of (Petersen graph): Graph on 5 vertices
             sage: I1.vertices(); I1.edges(labels=False)
             [0, 1, 2, 3, 8]
             [(0, 1), (1, 2), (2, 3), (3, 8)]

        It also contains the claw `K_{1,3}`::

             sage: h2 = g.subgraph_search(graphs.ClawGraph()); h2
             Subgraph of (Petersen graph): Graph on 4 vertices
             sage: h2.vertices(); h2.edges(labels=False)
             [0, 1, 4, 5]
             [(0, 1), (0, 4), (0, 5)]
             sage: I2 = g.subgraph_search(graphs.ClawGraph(), induced=True); I2
             Subgraph of (Petersen graph): Graph on 4 vertices
             sage: I2.vertices(); I2.edges(labels=False)
             [0, 1, 4, 5]
             [(0, 1), (0, 4), (0, 5)]

        Of course the induced copies are isomorphic to the graphs we were
        looking for::

             sage: I1.is_isomorphic(graphs.PathGraph(5))
             True
             sage: I2.is_isomorphic(graphs.ClawGraph())
             True

        However, the Petersen graph does not contain a subgraph isomorphic to
        `K_3`::

             sage: g.subgraph_search(graphs.CompleteGraph(3)) is None
             True

        Nor does it contain a nonempty induced subgraph isomorphic to `P_6`::

             sage: g.subgraph_search(graphs.PathGraph(6), induced=True) is None
             True

        The empty graph is a subgraph of every graph::

             sage: g.subgraph_search(graphs.EmptyGraph())
             Graph on 0 vertices
             sage: g.subgraph_search(graphs.EmptyGraph(), induced=True)
             Graph on 0 vertices

        The subgraph may just have edges missing::

            sage: k3=graphs.CompleteGraph(3); p3=graphs.PathGraph(3)
            sage: k3.relabel(list('abc'))
            sage: s=k3.subgraph_search(p3)
            sage: s.edges(labels=False)
            [('a', 'b'), ('b', 'c')]

        Of course, `P_3` is not an induced subgraph of `K_3`, though::

            sage: k3=graphs.CompleteGraph(3); p3=graphs.PathGraph(3)
            sage: k3.relabel(list('abc'))
            sage: k3.subgraph_search(p3, induced=True) is None
            True

        TESTS:

        Inside of a small graph (:trac:`13906`)::

            sage: Graph(5).subgraph_search(Graph(1))
            Graph on 1 vertex
        """
        from sage.graphs.generic_graph_pyx import SubgraphSearch
        from sage.graphs.graph_generators import GraphGenerators

        if G.order() == 0:
            return GraphGenerators().EmptyGraph()

        # SubgraphSearch assumes the graph we are searching for has order at least 2.
        if G.order() == 1:
            if self.order() >= 1:
                import graph
                return graph.Graph({ self.vertices()[0]:[]})
            else:
                return None

        S = SubgraphSearch(self, G, induced = induced)

        for g in S:
            if induced:
                return self.subgraph(g)
            else:
                Gcopy=G.copy()
                Gcopy.relabel(g)
                return self.subgraph(vertices=Gcopy.vertices(), edges=Gcopy.edges())

        return None

    def subgraph_search_count(self, G, induced=False):
        r"""
        Returns the number of labelled occurences of ``G`` in ``self``.

        INPUT:

        - ``G`` -- the graph whose copies we are looking for in
          ``self``.

        - ``induced`` -- boolean (default: ``False``). Whether or not
          to count induced copies of ``G`` in ``self``.

        ALGORITHM:

        Brute-force search.

        .. NOTE::

            This method also works on digraphs.

        .. SEEALSO::

            - :meth:`~GenericGraph.subgraph_search` -- finds an subgraph
              isomorphic to `H` inside of a graph `G`

            - :meth:`~GenericGraph.subgraph_search_iterator` -- Iterate on the
              copies of a graph `H` inside of a graph `G`

        EXAMPLES:

        Counting the number of paths `P_5` in a PetersenGraph::

            sage: g = graphs.PetersenGraph()
            sage: g.subgraph_search_count(graphs.PathGraph(5))
            240

        Requiring these subgraphs be induced::

            sage: g.subgraph_search_count(graphs.PathGraph(5), induced = True)
            120

        If we define the graph `T_k` (the transitive tournament on `k`
        vertices) as the graph on `\{0, ..., k-1\}` such that `ij \in
        T_k` iif `i<j`, how many directed triangles can be found in
        `T_5` ? The answer is of course `0` ::

             sage: T5 = DiGraph()
             sage: T5.add_edges([(i,j) for i in xrange(5) for j in xrange(i+1, 5)])
             sage: T5.subgraph_search_count(digraphs.Circuit(3))
             0

        If we count instead the number of `T_3` in `T_5`, we expect
        the answer to be `{5 \choose 3}`::

             sage: T3 = T5.subgraph([0,1,2])
             sage: T5.subgraph_search_count(T3)
             10
             sage: binomial(5,3)
             10

        The empty graph is a subgraph of every graph::

            sage: g.subgraph_search_count(graphs.EmptyGraph())
            1

        TESTS:

        Inside of a small graph (:trac:`13906`)::

            sage: Graph(5).subgraph_search_count(Graph(1))
            5
        """
        from sage.graphs.generic_graph_pyx import SubgraphSearch

        if G.order() == 0:
            return 1

        if self.order() == 0:
            return 0

        if G.order() == 1:
            return self.order()

        S = SubgraphSearch(self, G, induced = induced)

        return S.cardinality()

    def subgraph_search_iterator(self, G, induced=False):
        r"""
        Returns an iterator over the labelled copies of ``G`` in ``self``.

        INPUT:

        - ``G`` -- the graph whose copies we are looking for in
          ``self``.

        - ``induced`` -- boolean (default: ``False``). Whether or not
          to iterate over the induced copies of ``G`` in ``self``.

        ALGORITHM:

        Brute-force search.

        OUTPUT:

            Iterator over the labelled copies of ``G`` in ``self``, as
            *lists*. For each value `(v_1, v_2, ..., v_k)` returned,
            the first vertex of `G` is associated with `v_1`, the
            second with `v_2`, etc ...

        .. NOTE::

            This method also works on digraphs.

        .. SEEALSO::

            - :meth:`~GenericGraph.subgraph_search` -- finds an subgraph
              isomorphic to `H` inside of a graph `G`

            - :meth:`~GenericGraph.subgraph_search_count` -- Counts the number
              of copies of a graph `H` inside of a graph `G`

        EXAMPLE:

        Iterating through all the labelled `P_3` of `P_5`::

            sage: g = graphs.PathGraph(5)
            sage: for p in g.subgraph_search_iterator(graphs.PathGraph(3)):
            ...      print p
            [0, 1, 2]
            [1, 2, 3]
            [2, 1, 0]
            [2, 3, 4]
            [3, 2, 1]
            [4, 3, 2]

        TESTS:

        Inside of a small graph (:trac:`13906`)::

            sage: list(Graph(5).subgraph_search_iterator(Graph(1)))
            [Graph on 1 vertex, Graph on 1 vertex, Graph on 1 vertex, Graph on 1 vertex, Graph on 1 vertex]
        """

        if G.order() == 0:
            from sage.graphs.graph_generators import GraphGenerators
            return [GraphGenerators().EmptyGraph()]

        elif self.order() == 0:
            return []

        elif G.order() == 1:
            import graph
            return iter([graph.Graph({v:[]}) for v in self.vertices()])
        else:
            from sage.graphs.generic_graph_pyx import SubgraphSearch
            return SubgraphSearch(self, G, induced = induced)

    def random_subgraph(self, p, inplace=False):
        """
        Return a random subgraph that contains each vertex with prob. p.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.random_subgraph(.25)
            Subgraph of (Petersen graph): Graph on 4 vertices
        """
        vertices = []
        p = float(p)
        for v in self:
            if random() < p:
                vertices.append(v)
        return self.subgraph(vertices=vertices, inplace=inplace)

    def is_chordal(self, certificate = False, algorithm = "B"):
        r"""
        Tests whether the given graph is chordal.

        A Graph `G` is said to be chordal if it contains no induced hole (a
        cycle of length at least 4).

        Alternatively, chordality can be defined using a Perfect Elimination
        Order :

        A Perfect Elimination Order of a graph `G` is an ordering `v_1,...,v_n`
        of its vertex set such that for all `i`, the neighbors of `v_i` whose
        index is greater that `i` induce a complete subgraph in `G`. Hence, the
        graph `G` can be totally erased by successively removing vertices whose
        neighborhood is a clique (also called *simplicial* vertices)
        [Fulkerson65]_.

        (It can be seen that if `G` contains an induced hole, then it can not
        have a perfect elimination order. Indeed, if we write `h_1,...,h_k` the
        `k` vertices of such a hole, then the first of those vertices to be
        removed would have two non-adjacent neighbors in the graph.)

        A Graph is then chordal if and only if it has a Perfect Elimination
        Order.

        INPUT:

        - ``certificate`` (boolean) -- Whether to return a certificate.

            * If ``certificate = False`` (default), returns ``True`` or
              ``False`` accordingly.

            * If ``certificate = True``, returns :

                * ``(True, peo)`` when the graph is chordal, where ``peo`` is a
                  perfect elimination order of its vertices.

                * ``(False, Hole)`` when the graph is not chordal, where
                  ``Hole`` (a ``Graph`` object) is an induced subgraph of
                  ``self`` isomorphic to a hole.

        - ``algorithm`` -- Two algorithms are available for this method (see
          next section), which can be selected by setting ``algorithm = "A"`` or
          ``algorithm = "B"`` (default). While they will agree on whether the
          given graph is chordal, they can not be expected to return the same
          certificates.

        ALGORITHM:

        This algorithm works through computing a Lex BFS on the graph, then
        checking whether the order is a Perfect Elimination Order by computing
        for each vertex `v` the subgraph induces by its non-deleted neighbors,
        then testing whether this graph is complete.

        This problem can be solved in `O(m)` [Rose75]_ ( where `m` is the number
        of edges in the graph ) but this implementation is not linear because of
        the complexity of Lex BFS.

        .. NOTE::

            Because of a past bug (#11735, #11961), the first implementation
            (algorithm A) of this method sometimes returned as certificates
            subgraphs which were **not** holes. Since then, this bug has been
            fixed and the values are now double-checked before being returned,
            so that the algorithm only returns correct values or raises an
            exception. In the case where an exception is raised, the user is
            advised to switch to the other algorithm. And to **please** report
            the bug :-)

        EXAMPLES:

        The lexicographic product of a Path and a Complete Graph
        is chordal ::

            sage: g = graphs.PathGraph(5).lexicographic_product(graphs.CompleteGraph(3))
            sage: g.is_chordal()
            True

        The same goes with the product of a random lobster
        ( which is a tree ) and a Complete Graph ::

            sage: g = graphs.RandomLobster(10,.5,.5).lexicographic_product(graphs.CompleteGraph(3))
            sage: g.is_chordal()
            True

        The disjoint union of chordal graphs is still chordal::

            sage: (2*g).is_chordal()
            True

        Let us check the certificate given by Sage is indeed a perfect elimintion order::

            sage: (_, peo) = g.is_chordal(certificate = True)
            sage: for v in peo:
            ...       if not g.subgraph(g.neighbors(v)).is_clique():
            ...            print "This should never happen !"
            ...       g.delete_vertex(v)
            sage: print "Everything is fine !"
            Everything is fine !

        Of course, the Petersen Graph is not chordal as it has girth 5 ::

            sage: g = graphs.PetersenGraph()
            sage: g.girth()
            5
            sage: g.is_chordal()
            False

        We can even obtain such a cycle as a certificate ::

            sage: (_, hole) = g.is_chordal(certificate = True)
            sage: hole
            Subgraph of (Petersen graph): Graph on 5 vertices
            sage: hole.is_isomorphic(graphs.CycleGraph(5))
            True

        TESTS:

        This shouldn't raise exceptions (:trac:`10899`)::

            sage: Graph(1).is_chordal()
            True
            sage: for g in graphs(5):
            ....:     _ = g.is_chordal()

        REFERENCES:

        .. [Rose75] Rose, D.J. and Tarjan, R.E.,
          Algorithmic aspects of vertex elimination,
          Proceedings of seventh annual ACM symposium on Theory of computing
          Page 254, ACM 1975

        .. [Fulkerson65] Fulkerson, D.R. and Gross, OA
          Incidence matrices and interval graphs
          Pacific J. Math 1965
          Vol. 15, number 3, pages 835--855

        TESTS:

        Trac Ticket #11735::

           sage: g = Graph({3:[2,1,4],2:[1],4:[1],5:[2,1,4]})
           sage: _, g1 = g.is_chordal(certificate=True); g1.is_chordal()
           False
           sage: g1.is_isomorphic(graphs.CycleGraph(g1.order()))
           True
        """
        self._scream_if_not_simple()

        # If the graph is not connected, we are computing the result on each component
        if not self.is_connected():

            # If the user wants a certificate, we had no choice but to
            # collect the perfect elimination orders... But we return
            # a hole immediately if we find any !
            if certificate:
                peo = []
                for gg in self.connected_components_subgraphs():

                    b, certif = gg.is_chordal(certificate = True)
                    if not b:
                        return False, certif
                    else:
                        peo.extend(certif)

                return True, peo

            # One line if no certificate is requested
            else:
                return all( gg.is_chordal() for gg in self.connected_components_subgraphs() )

        hole = None
        g = self.copy()

        # Algorithms
        #
        # They find the perfect elimination ordering or produce a hole

        if algorithm == "A":

            peo,t_peo = self.lex_BFS(tree=True)
            peo.reverse()

            # Iteratively removing vertices and checking everything is fine.
            for v in peo:

                if t_peo.out_degree(v) == 0:
                    g.delete_vertex(v)
                    continue

                x = t_peo.neighbor_out_iterator(v).next()
                S = self.neighbors(x)+[x]

                if not frozenset(g.neighbors(v)).issubset(S):

                    # Do we need to return a hole ?
                    if certificate:

                        # In this case, let us take two nonadjacent neighbors of v
                        # In this case, let us take two nonadjacent neighbors of
                        # v. In order to do so, we pick a vertex y which is a
                        # neighbor of v but is not adjacent to x, which we know
                        # exists by the test written two lines above.

                        for y in g.neighbors(v):
                            if y not in S:
                                break

                        g.delete_vertices([vv for vv in g.neighbors(v) if vv != y and vv != x])
                        g.delete_vertex(v)

                        # Our hole is v + (a shortest path between x and y not
                        # containing v or any of its neighbors).

                        hole = self.subgraph([v] + g.shortest_path(x,y))

                        # End of the algorithm
                        break
                    else:
                        return False

                g.delete_vertex(v)

        elif algorithm == "B":

            peo,t_peo = self.lex_BFS(reverse=True, tree=True)

            # Remembering the (closed) neighborhoods of each vertex
            neighbors_subsets = dict([(v,self.neighbors(v)+[v]) for v in g])
            pos_in_peo = dict(zip(peo, range(self.order())))

            # Iteratively removing vertices and checking everything is fine.
            for v in reversed(peo):

                if (t_peo.out_degree(v)>0 and
                    not frozenset([v1 for v1 in g.neighbors(v) if pos_in_peo[v1] > pos_in_peo[v]]).issubset(
                        neighbors_subsets[t_peo.neighbor_out_iterator(v).next()])):

                    # Do we need to return a hole ?
                    if certificate:

                        # In this case, let us take two nonadjacent neighbors of
                        # v. In order to do so, we pick a vertex y which is a
                        # neighbor of v but is not adjacent to x, which we know
                        # exists by the test written two lines above.
                        max_tup = (-1, 0)
                        nb1 = [u for u in g.neighbors(v) if pos_in_peo[u] > pos_in_peo[v]]
                        for xi in nb1:
                            for yi in nb1:
                                if not yi in neighbors_subsets[xi]:
                                    new_tup = (pos_in_peo[xi], pos_in_peo[yi])
                                    if max_tup < new_tup:
                                        max_tup = new_tup
                                        x, y = xi, yi

                        # Our hole is v + (a shortest path between x and y not
                        # containing v or any of its neighbors).

                        #g.delete_vertices([vv for vv in g.vertices() if pos_in_peo[vv] < pos_in_peo[v]])

                        g.delete_vertices([vv for vv in g.neighbors(v) if vv != y and vv != x])
                        g.delete_vertex(v)

                        hole = self.subgraph([v] + g.shortest_path(x,y))

                        # End of the algorithm
                        break
                    else:
                        return False


        # Returning values
        # ----------------

        # 1- The graph is not chordal

        if not hole is None:
            # There was a bug there once, so it's better to check the
            # answer is valid, especally when it is so cheap ;-)

            if hole.order() <= 3 or not hole.is_regular(k=2):
                raise RuntimeError("the graph is not chordal, and something went wrong in the computation of the certificate. Please report this bug, providing the graph if possible!")

            return (False, hole)


        # 2- The graph is chordal
        if certificate:
            return True, peo

        else:
            return True

    def is_circulant(self, certificate = False):
        r"""
        Tests whether the graph is circulant.

        For more information on circulant graphs, see the
        :wikipedia:`Wikipedia page on circulant graphs
        <Circulant_graph>`.

        INPUT:

        - ``certificate`` (boolean) -- whether to return a certificate for
          yes-answers. See OUTPUT section. Set to ``False`` by default.

        OUTPUT:

        When ``certificate`` is set to ``False`` (default) this method only
        returns ``True`` or ``False`` answers. When ``certificate`` is set to
        ``True``, the method either returns ``(False, None)`` or ``(True,
        lists_of_parameters)`` each element of ``lists_of_parameters`` can be
        used to define the graph as a circulant graph.

        See the documentation of
        :func:`~sage.graphs.graph_generators.GraphGenerators.CirculantGraph` and
        :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.Circulant` for
        more information, and the examples below.

        .. SEEALSO::

            :meth:`~sage.graphs.graph_generators.GraphGenerators.CirculantGraph`
            -- a constructor for circulant graphs.

        EXAMPLES:

        The Petersen graph is not a circulant graph::

            sage: g = graphs.PetersenGraph()
            sage: g.is_circulant()
            False

        A cycle is obviously a circulant graph, but several sets of parameters
        can be used to define it::

            sage: g = graphs.CycleGraph(5)
            sage: g.is_circulant(certificate = True)
            (True, [(5, [1, 4]), (5, [2, 3])])

        The same goes for directed graphs::

            sage: g = digraphs.Circuit(5)
            sage: g.is_circulant(certificate = True)
            (True, [(5, [1]), (5, [3]), (5, [2]), (5, [4])])

        With this information, it is very easy to create (and plot) all possible
        drawings of a circulant graph::

            sage: g = graphs.CirculantGraph(13, [2, 3, 10, 11])
            sage: for param in g.is_circulant(certificate = True)[1]:
            ...      graphs.CirculantGraph(*param)
            Circulant graph ([2, 3, 10, 11]): Graph on 13 vertices
            Circulant graph ([1, 5, 8, 12]): Graph on 13 vertices
            Circulant graph ([4, 6, 7, 9]): Graph on 13 vertices

        TESTS::

            sage: digraphs.DeBruijn(3,1).is_circulant(certificate = True)
            (True, [(3, [0, 1, 2])])
            sage: Graph(1).is_circulant(certificate = True)
            (True, (1, []))
            sage: Graph(0).is_circulant(certificate = True)
            (True, (0, []))
            sage: Graph([(0,0)]).is_circulant(certificate = True)
            (True, (1, [0]))
        """
        self._scream_if_not_simple(allow_loops=True)
        # Stupid cases
        if self.order() <= 1:
            if certificate:
                return (True,(self.order(),[0] if self.size() else []))
            else:
                return True

        certif_list = []

        # The automorphism group, the translation between the vertices of self
        # and 1..n, and the orbits.
        ag, orbits = self.automorphism_group([self.vertices()],
                          order=False,
                          return_group=True,
                          orbits=True)

        # Not transitive ? Not a circulant graph !
        if len(orbits) != 1:
            return (False, None) if certificate else False

        # We go through all conjugacy classes of the automorphism
        # group, and only keep the cycles of length n
        for e in ag.conjugacy_classes_representatives():
            cycles = e.cycle_tuples()

            # If the automorphism is not the identity and has exactly one
            # cycle that contains all vertices.
            if ((not cycles) or
                len(cycles[0]) != self.order()):
                continue

            # From now on the graph is a circulant graph !

            if not certificate:
                return True

            # We build the list of integers defining the circulant graph, and
            # add it to the list.
            parameters = []
            cycle = cycles[0]
            u = cycle[0]
            integers = [i for i,v in enumerate(cycle) if self.has_edge(u,v)]
            certif_list.append((self.order(),integers))

        if not certificate:
            return False
        else:
            if certif_list:
                return (True, certif_list)
            else:
                return (False, None)

    def is_interval(self, certificate = False):
        r"""
        Check whether self is an interval graph

        INPUT:

        - ``certificate`` (boolean) -- The function returns ``True``
          or ``False`` according to the graph, when ``certificate =
          False`` (default). When ``certificate = True`` and the graph
          is an interval graph, a dictionary whose keys are the
          vertices and values are pairs of integers are returned
          instead of ``True``. They correspond to an embedding of the
          interval graph, each vertex being represented by an interval
          going from the first of the two values to the second.

        ALGORITHM:

        Through the use of PQ-Trees

        AUTHOR:

        Nathann Cohen (implementation)

        EXAMPLES:

        A Petersen Graph is not chordal, nor car it be an interval
        graph ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_interval()
            False

        Though we can build intervals from the corresponding random
        generator::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_interval()
            True

        This method can also return, given an interval graph, a
        possible embedding (we can actually compute all of them
        through the PQ-Tree structures)::

            sage: g = Graph(':S__@_@A_@AB_@AC_@ACD_@ACDE_ACDEF_ACDEFG_ACDEGH_ACDEGHI_ACDEGHIJ_ACDEGIJK_ACDEGIJKL_ACDEGIJKLMaCEGIJKNaCEGIJKNaCGIJKNPaCIP', loops=False, multiedges=False)
            sage: d = g.is_interval(certificate = True)
            sage: print d                                    # not tested
            {0: (0, 20), 1: (1, 9), 2: (2, 36), 3: (3, 5), 4: (4, 38), 5: (6, 21), 6: (7, 27), 7: (8, 12), 8: (10, 29), 9: (11, 16), 10: (13, 39), 11: (14, 31), 12: (15, 32), 13: (17, 23), 14: (18, 22), 15: (19, 33), 16: (24, 25), 17: (26, 35), 18: (28, 30), 19: (34, 37)}

        From this embedding, we can clearly build an interval graph
        isomorphic to the previous one::

            sage: g2 = graphs.IntervalGraph(d.values())
            sage: g2.is_isomorphic(g)
            True

        .. SEEALSO::

        - :mod:`Interval Graph Recognition <sage.graphs.pq_trees>`.

        - :meth:`PQ <sage.graphs.pq_trees.PQ>`
          -- Implementation of PQ-Trees.

        """
        self._scream_if_not_simple()

        # An interval graph first is a chordal graph. Without this,
        # there is no telling how we should find its maximal cliques,
        # by the way :-)

        if not self.is_chordal():
            return False

        # First, we need to gather the list of maximal cliques, which
        # is easy as the graph is chordal

        cliques = []

        # As we will be deleting vertices ...
        g = self.copy()

        for cc in self.connected_components_subgraphs():

            # We pick a perfect elimination order for every connected
            # component. We will then iteratively take the last vertex
            # in the order (a simplicial vertex) and consider the
            # clique it forms with its neighbors. If we do not have an
            # inclusion-wise larger clique in our list, we add it !

            peo = cc.lex_BFS()



            while peo:
                v = peo.pop()
                clique = frozenset( [v] + cc.neighbors(v))
                cc.delete_vertex(v)

                if not any([clique.issubset(c) for c in cliques]):
                    cliques.append(clique)

        from sage.graphs.pq_trees import reorder_sets

        try:
            ordered_sets = reorder_sets(cliques)
            if not certificate:
                return True

        except ValueError:
            return False

        # We are now listing the maximal cliques in the given order,
        # and keeping track of the vertices appearing/disappearing

        current = set([])
        beg = {}
        end = {}

        i = 0

        ordered_sets.append([])
        for S in map(set,ordered_sets):
            for v in current-S:
                end[v] = i
                i = i + 1

            for v in S-current:
                beg[v] = i
                i = i + 1

            current = S


        return dict([(v, (beg[v], end[v])) for v in self])


    def is_gallai_tree(self):
        r"""
        Returns whether the current graph is a Gallai tree.

        A graph is a Gallai tree if and only if it is
        connected and its `2`-connected components are all
        isomorphic to complete graphs or odd cycles.

        A connected graph is not degree-choosable if and
        only if it is a Gallai tree [erdos1978choos]_.

        REFERENCES:

        .. [erdos1978choos] Erdos, P. and Rubin, A.L. and Taylor, H.
          Proc. West Coast Conf. on Combinatorics
          Graph Theory and Computing, Congressus Numerantium
          vol 26, pages 125--157, 1979

        EXAMPLES:

        A complete graph is, or course, a Gallai Tree::

            sage: g = graphs.CompleteGraph(15)
            sage: g.is_gallai_tree()
            True

        The Petersen Graph is not::

            sage: g = graphs.PetersenGraph()
            sage: g.is_gallai_tree()
            False

        A Graph built from vertex-disjoint complete graphs
        linked by one edge to a special vertex `-1` is a
        ''star-shaped'' Gallai tree ::

            sage: g = 8 * graphs.CompleteGraph(6)
            sage: g.add_edges([(-1,c[0]) for c in g.connected_components()])
            sage: g.is_gallai_tree()
            True
        """
        self._scream_if_not_simple()
        if not self.is_connected():
            return False

        for c in self.blocks_and_cut_vertices()[0]:
            gg = self.subgraph(c)
            #                    is it an odd cycle ?              a complete graph ?
            if not ( (len(c)%2 == 1 and gg.size() == len(c)+1) or gg.is_clique() ):
                return False

        return True

    def is_clique(self, vertices=None, directed_clique=False):
        """
        Tests whether a set of vertices is a clique

        A clique is a set of vertices such that there is an edge between any two
        vertices.

        INPUT:

        -  ``vertices`` - Vertices can be a single vertex or an
           iterable container of vertices, e.g. a list, set, graph, file or
           numeric array. If not passed, defaults to the entire graph.

        -  ``directed_clique`` - (default False) If set to
           False, only consider the underlying undirected graph. If set to
           True and the graph is directed, only return True if all possible
           edges in _both_ directions exist.

        EXAMPLES::

            sage: g = graphs.CompleteGraph(4)
            sage: g.is_clique([1,2,3])
            True
            sage: g.is_clique()
            True
            sage: h = graphs.CycleGraph(4)
            sage: h.is_clique([1,2])
            True
            sage: h.is_clique([1,2,3])
            False
            sage: h.is_clique()
            False
            sage: i = graphs.CompleteGraph(4).to_directed()
            sage: i.delete_edge([0,1])
            sage: i.is_clique()
            True
            sage: i.is_clique(directed_clique=True)
            False
        """
        if directed_clique and self._directed:
            subgraph=self.subgraph(vertices)
            subgraph.allow_loops(False)
            subgraph.allow_multiple_edges(False)
            n=subgraph.order()
            return subgraph.size()==n*(n-1)
        else:
            if vertices is None:
                subgraph = self
            else:
                subgraph=self.subgraph(vertices)

            if self._directed:
                subgraph = subgraph.to_simple()

            n=subgraph.order()
            return subgraph.size()==n*(n-1)/2

    def is_independent_set(self, vertices=None):
        """
        Returns True if the set ``vertices`` is an independent
        set, False if not. An independent set is a set of vertices such
        that there is no edge between any two vertices.

        INPUT:

        -  ``vertices`` - Vertices can be a single vertex or an
           iterable container of vertices, e.g. a list, set, graph, file or
           numeric array. If not passed, defaults to the entire graph.

        EXAMPLES::

            sage: graphs.CycleGraph(4).is_independent_set([1,3])
            True
            sage: graphs.CycleGraph(4).is_independent_set([1,2,3])
            False
        """
        return self.subgraph(vertices).size()==0

    def is_subgraph(self, other, induced=True):
        """
        Tests whether ``self`` is a subgraph of ``other``.

        .. WARNING::

            Please note that this method does not check whether ``self``
            contains a subgraph *isomorphic* to ``other``, but only if it
            directly contains it as a subgraph !

            By default ``induced`` is ``True`` for backwards compatibility.

        INPUT:

        - ``induced`` - boolean (default: ``True``) If set to ``True`` tests
          whether ``self`` is an *induced* subgraph of ``other`` that is if
          the vertices of ``self`` are also vertices of ``other``, and the
          edges of  ``self`` are equal to the edges of ``other`` between the
          vertices contained in ``self`.
          If set to ``False`` tests whether ``self`` is a subgraph of ``other``
          that is if all vertices of ``self`` are also in ``other`` and all
          edges of ``self`` are also in ``other``.

        OUTPUT:

        boolean -- ``True`` iff ``self`` is a (possibly induced) subgraph of ``other``.

        .. SEEALSO::

            If you are interested in the (possibly induced) subgraphs isomorphic
            to ``self`` in ``other``, you are looking for the following methods:

            - :meth:`~GenericGraph.subgraph_search` -- finds a subgraph
              isomorphic to `G` inside of a `self`

            - :meth:`~GenericGraph.subgraph_search_count` -- Counts the number
              of such copies.

            - :meth:`~GenericGraph.subgraph_search_iterator` -- Iterate over all
              the copies of `G` contained in `self`.

        EXAMPLES:

            sage: P = graphs.PetersenGraph()
            sage: G = P.subgraph(range(6))
            sage: G.is_subgraph(P)
            True

            sage: H=graphs.CycleGraph(5)
            sage: G=graphs.PathGraph(5)
            sage: G.is_subgraph(H)
            False
            sage: G.is_subgraph(H, induced=False)
            True
            sage: H.is_subgraph(G, induced=False)
            False

        TESTS:

        Raise an error when self and other are of different types::

            sage: Graph([(0,1)]).is_subgraph( DiGraph([(0,1)]) )
            Traceback (most recent call last):
            ...
            ValueError: The input parameter must be a Graph.
            sage: DiGraph([(0,1)]).is_subgraph( Graph([(0,1)]) )
            Traceback (most recent call last):
            ...
            ValueError: The input parameter must be a DiGraph.
        """
        from sage.graphs.graph import Graph
        from sage.graphs.digraph import DiGraph
        if isinstance(self,Graph) and not isinstance(other,Graph):
            raise ValueError('The input parameter must be a Graph.')

        if isinstance(self,DiGraph) and not isinstance(other,DiGraph):
            raise ValueError('The input parameter must be a DiGraph.')

        if self.num_verts() > other.num_verts():
            return False

        if any(v not in other for v in self.vertex_iterator()):
            return False

        if induced:
            return other.subgraph(self.vertices()) == self
        else:
            self._scream_if_not_simple(allow_loops=True)
            return all(other.has_edge(e) for e in self.edge_iterator())

    ### Cluster

    def cluster_triangles(self, nbunch=None, with_labels=False):
        r"""
        Returns the number of triangles for the set `nbunch` of vertices
        as a dictionary keyed by vertex.

        See also section "Clustering" in chapter "Algorithms" of [HSSNX]_.

        INPUT:

        -  ``nbunch`` - The vertices to inspect. If ``nbunch=None``,
           returns data for all vertices in the graph.

        REFERENCE:

        .. [HSSNX] Aric Hagberg, Dan Schult and Pieter Swart. NetworkX
          documentation. [Online] Available:
          http://networkx.github.io/documentation/latest/reference/index.html

        EXAMPLES::

            sage: (graphs.FruchtGraph()).cluster_triangles().values()
            [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0]
            sage: (graphs.FruchtGraph()).cluster_triangles()
            {0: 1, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 0, 9: 1, 10: 1, 11: 0}
            sage: (graphs.FruchtGraph()).cluster_triangles(nbunch=[0,1,2])
            {0: 1, 1: 1, 2: 0}
        """
        import networkx
        return networkx.triangles(self.networkx_graph(copy=False), nbunch)

    def clustering_average(self):
        r"""
        Returns the average clustering coefficient.

        The clustering coefficient of a node `i` is the fraction
        of existing triangles containing node `i` and all
        possible triangles containing `i`: `c_i = T(i) / \binom {k_i} 2`
        where `T(i)` is the number of existing triangles through `i`, and
        `k_i` is the degree of vertex `i`.

        A coefficient for the whole graph is the average of the `c_i`.

        See also section "Clustering" in chapter "Algorithms" of [HSSNX]_.

        EXAMPLES::

            sage: (graphs.FruchtGraph()).clustering_average()
            0.25
        """
        import networkx
        return networkx.average_clustering(self.networkx_graph(copy=False))

    def clustering_coeff(self, nodes=None, weight=False, return_vertex_weights=True):
        r"""
        Returns the clustering coefficient for each vertex in ``nodes`` as
        a dictionary keyed by vertex.

        For an unweighted graph, the clustering coefficient of a node `i`
        is the fraction of existing triangles containing node `i` and all
        possible triangles containing `i`: `c_i = T(i) / \binom {k_i} 2`
        where `T(i)` is the number of existing triangles through `i`, and
        `k_i` is the degree of vertex `i`.

        For weighted graphs the clustering is defined as the geometric
        average of the subgraph edge weights, normalized by the
        maximum weight in the network.

        The value of `c_i` is assigned `0` if `k_i < 2`.

        See also section "Clustering" in chapter "Algorithms" of [HSSNX]_.

        INPUT:

        - ``nodes`` - the vertices to inspect (default ``None``, returns data
          on all vertices in graph)

        - ``weight`` - string or boolean (default is ``False``). If it is
          a string it used the indicated edge property as weight.
          ``weight = True`` is equivalent to ``weight = 'weight'``

        - ``return_vertex_weights`` is a boolean ensuring backwards
          compatibility with deprecated features of NetworkX 1.2. It
          should be set to ``False`` for all production code.

        EXAMPLES::

            sage: (graphs.FruchtGraph()).clustering_coeff().values()
            [0.3333333333333333, 0.3333333333333333, 0.0, 0.3333333333333333,
             0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
             0.3333333333333333, 0.0, 0.3333333333333333, 0.3333333333333333,
             0.0]
            sage: (graphs.FruchtGraph()).clustering_coeff()
            {0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0,
             3: 0.3333333333333333, 4: 0.3333333333333333,
             5: 0.3333333333333333, 6: 0.3333333333333333,
             7: 0.3333333333333333, 8: 0.0, 9: 0.3333333333333333,
             10: 0.3333333333333333, 11: 0.0}

            sage: (graphs.FruchtGraph()).clustering_coeff(weight=True,
            ...     return_vertex_weights=False)
            {0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0,
            3: 0.3333333333333333, 4: 0.3333333333333333,
            5: 0.3333333333333333, 6: 0.3333333333333333,
            7: 0.3333333333333333, 8: 0.0, 9: 0.3333333333333333,
            10: 0.3333333333333333, 11: 0.0}
            sage: (graphs.FruchtGraph()).clustering_coeff(nodes=[0,1,2])
            {0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0}

            sage: (graphs.FruchtGraph()).clustering_coeff(nodes=[0,1,2],
            ...     weight=True, return_vertex_weights=False)
            {0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0}

        TESTS:

        Doctests that demonstrate the deprecation of the two-dictionary
        return value due to the NetworkX API change after 1.2. The
        return_vertex_weights keyword is provided with a default value
        of True for backwards compatibility with older versions of Sage.
        When the deprecation period has expired and the keyword is
        removed, these doctests should be removed as well. ::

            sage: (graphs.FruchtGraph()).clustering_coeff(weight=True,
            ...     return_vertex_weights=True)
            doctest:...: DeprecationWarning: The option 'return_vertex_weights'
            has been deprecated. Only offered for backwards compatibility with
            NetworkX 1.2.
            See http://trac.sagemath.org/12806 for details.
            ({0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0,
              3: 0.3333333333333333, 4: 0.3333333333333333,
              5: 0.3333333333333333, 6: 0.3333333333333333,
              7: 0.3333333333333333, 8: 0.0, 9: 0.3333333333333333,
              10: 0.3333333333333333, 11: 0.0}, {0: 0.08333333333333333,
              1: 0.08333333333333333, 2: 0.08333333333333333,
              3: 0.08333333333333333, 4: 0.08333333333333333,
              5: 0.08333333333333333, 6: 0.08333333333333333,
              7: 0.08333333333333333, 8: 0.08333333333333333,
              9: 0.08333333333333333, 10: 0.08333333333333333,
              11: 0.08333333333333333})

            sage: (graphs.FruchtGraph()).clustering_coeff(nodes=[0, 1, 2],
            ...     weight=True, return_vertex_weights=True)
            ({0: 0.3333333333333333, 1: 0.3333333333333333, 2: 0.0},
             {0: 0.3333333333333333, 1: 0.3333333333333333,
             2: 0.3333333333333333})
        """
        import networkx
        if weight and return_vertex_weights:
            # Running in compatibility mode with deprecated NetworkX 1.2 features
            # All this code should be removed when the deprecation warning expires
            from sage.misc.superseded import deprecation
            deprecation(12806, "The option 'return_vertex_weights' has been " +\
                        "deprecated. Only offered for backwards" +\
                        " compatibility with NetworkX 1.2.")
            G = self.networkx_graph(copy=False)
            if G.is_directed():
                raise NetworkXError("Clustering algorithms are not defined for directed graphs.")
            clusterc={}
            weights={}
            for v,d,t in networkx.cluster._triangles_and_degree_iter(G,nodes):
                weights[v]=float(d*(d-1))
                if t==0:
                    clusterc[v]=0.0
                else:
                    clusterc[v]=t/float(d*(d-1))
            scale=1./sum(weights.itervalues())
            for v,w in weights.iteritems():
                weights[v]=w*scale
            return clusterc,weights

        else:
            return networkx.clustering(self.networkx_graph(copy=False), nodes, weight=weight)

    def cluster_transitivity(self):
        r"""
        Returns the transitivity (fraction of transitive triangles) of the
        graph.

        Transitivity is the fraction of all existing triangles and all
        connected triples (triads), `T = 3\times\text{triangles}
        / \text{triads}`.

        See also section "Clustering" in chapter "Algorithms" of [HSSNX]_.

        EXAMPLES::

            sage: (graphs.FruchtGraph()).cluster_transitivity()
            0.25
        """
        import networkx
        return networkx.transitivity(self.networkx_graph(copy=False))

    ### Distance

    def distance(self, u, v, by_weight=False):
        """
        Returns the (directed) distance from u to v in the (di)graph, i.e.
        the length of the shortest path from u to v.

        INPUT:

        -  ``by_weight`` - if ``False``, uses a breadth first
           search. If True, takes edge weightings into account, using
           Dijkstra's algorithm.

        EXAMPLES::

            sage: G = graphs.CycleGraph(9)
            sage: G.distance(0,1)
            1
            sage: G.distance(0,4)
            4
            sage: G.distance(0,5)
            4
            sage: G = Graph( {0:[], 1:[]} )
            sage: G.distance(0,1)
            +Infinity
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse = True)
            sage: G.plot(edge_labels=True).show() # long time
            sage: G.distance(0, 3)
            2
            sage: G.distance(0, 3, by_weight=True)
            3

        """
        return self.shortest_path_length(u, v, by_weight = by_weight)

    def distance_all_pairs(self, algorithm = "auto"):
        r"""
        Returns the distances between all pairs of vertices.

        INPUT:

            - ``"algorithm"`` (string) -- two algorithms are available

                * ``algorithm = "BFS"`` in which case the distances are computed
                  through `n` different breadth-first-search.

                * ``algorithm = "Floyd-Warshall"``, in which case the
                  Floyd-Warshall algorithm is used.

                * ``algorithm = "auto"``, in which case the Floyd-Warshall
                  algorithm is used for graphs on less than 20 vertices, and BFS
                  otherwise.

                The default is ``algorithm = "BFS"``.

        OUTPUT:

        A doubly indexed dictionary

        .. NOTE::

           There is a Cython version of this method that is usually
           much faster for large graphs, as most of the time is
           actually spent building the final double
           dictionary. Everything on the subject is to be found in the
           :mod:`~sage.graphs.distances_all_pairs` module.

        EXAMPLE:

        The Petersen Graph::

            sage: g = graphs.PetersenGraph()
            sage: print g.distance_all_pairs()
            {0: {0: 0, 1: 1, 2: 2, 3: 2, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2, 9: 2}, 1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 1, 7: 2, 8: 2, 9: 2}, 2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 2, 5: 2, 6: 2, 7: 1, 8: 2, 9: 2}, 3: {0: 2, 1: 2, 2: 1, 3: 0, 4: 1, 5: 2, 6: 2, 7: 2, 8: 1, 9: 2}, 4: {0: 1, 1: 2, 2: 2, 3: 1, 4: 0, 5: 2, 6: 2, 7: 2, 8: 2, 9: 1}, 5: {0: 1, 1: 2, 2: 2, 3: 2, 4: 2, 5: 0, 6: 2, 7: 1, 8: 1, 9: 2}, 6: {0: 2, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 0, 7: 2, 8: 1, 9: 1}, 7: {0: 2, 1: 2, 2: 1, 3: 2, 4: 2, 5: 1, 6: 2, 7: 0, 8: 2, 9: 1}, 8: {0: 2, 1: 2, 2: 2, 3: 1, 4: 2, 5: 1, 6: 1, 7: 2, 8: 0, 9: 2}, 9: {0: 2, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 2, 9: 0}}

        Testing on Random Graphs::

            sage: g = graphs.RandomGNP(20,.3)
            sage: distances = g.distance_all_pairs()
            sage: all([g.distance(0,v) == distances[0][v] for v in g])
            True

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.distance_matrix`
        """
        if algorithm == "auto":
            if self.order() <= 20:
                algorithm = "Floyd-Warshall"
            else:
                algorithm = "BFS"

        if algorithm == "BFS":
            from sage.graphs.distances_all_pairs import distances_all_pairs
            return distances_all_pairs(self)

        elif algorithm == "Floyd-Warshall":
            from sage.graphs.distances_all_pairs import floyd_warshall
            return floyd_warshall(self,paths = False, distances = True)

        else:
            raise ValueError("The algorithm keyword can be equal to either \"BFS\" or \"Floyd-Warshall\" or \"auto\"")


    def eccentricity(self, v=None, dist_dict=None, with_labels=False):
        """
        Return the eccentricity of vertex (or vertices) v.

        The eccentricity of a vertex is the maximum distance to any other
        vertex.

        INPUT:


        -  ``v`` - either a single vertex or a list of
           vertices. If it is not specified, then it is taken to be all
           vertices.

        -  ``dist_dict`` - optional, a dict of dicts of
           distance.

        -  ``with_labels`` - Whether to return a list or a
           dict.


        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.eccentricity()
            [4, 4, 4, 4, 4, 3, 3, 2, 3, 4]
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: G.eccentricity(7)
            2
            sage: G.eccentricity([7,8,9])
            [3, 4, 2]
            sage: G.eccentricity([7,8,9], with_labels=True) == {8: 3, 9: 4, 7: 2}
            True
            sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
            sage: G.eccentricity()
            [+Infinity, +Infinity, +Infinity]
            sage: G = Graph({0:[]})
            sage: G.eccentricity(with_labels=True)
            {0: 0}
            sage: G = Graph({0:[], 1:[]})
            sage: G.eccentricity(with_labels=True)
            {0: +Infinity, 1: +Infinity}
        """
        if v is None:
            if dist_dict is None:
                from sage.graphs.distances_all_pairs import eccentricity

                if with_labels:
                    return dict(zip(self.vertices(), eccentricity(self)))
                else:
                    return eccentricity(self)

            v = self.vertices()
        elif not isinstance(v, list):
            v = [v]
        e = {}
        infinite = False
        for u in v:
            if dist_dict is None:
                length = self.shortest_path_lengths(u)
            else:
                length = dist_dict[u]
            if len(length) != self.num_verts():
                infinite = True
                break
            e[u] = max(length.values())
        if infinite:
            from sage.rings.infinity import Infinity
            for u in v:
                e[u] = Infinity
        if with_labels:
            return e
        else:
            if len(e)==1: return e.values()[0] # return single value
            return e.values()

    def radius(self):
        """
        Returns the radius of the (di)graph.

        The radius is defined to be the minimum eccentricity of any vertex,
        where the eccentricity is the maximum distance to any other
        vertex.

        EXAMPLES: The more symmetric a graph is, the smaller (diameter -
        radius) is.

        ::

            sage: G = graphs.BarbellGraph(9, 3)
            sage: G.radius()
            3
            sage: G.diameter()
            6

        ::

            sage: G = graphs.OctahedralGraph()
            sage: G.radius()
            2
            sage: G.diameter()
            2

        TEST::

            sage: g = Graph()
            sage: g.radius()
            Traceback (most recent call last):
            ...
            ValueError: This method has no meaning on empty graphs.
        """
        if self.order() == 0:
            raise ValueError("This method has no meaning on empty graphs.")

        return min(self.eccentricity())

    def center(self):
        """
        Returns the set of vertices in the center, i.e. whose eccentricity
        is equal to the radius of the (di)graph.

        In other words, the center is the set of vertices achieving the
        minimum eccentricity.

        EXAMPLES::

            sage: G = graphs.DiamondGraph()
            sage: G.center()
            [1, 2]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.center()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.center()
            [0]
            sage: G = Graph()
            sage: G.center()
            []
            sage: G.add_vertex()
            0
            sage: G.center()
            [0]
        """
        e = self.eccentricity(with_labels=True)
        try:
            r = min(e.values())
        except Exception:
            return []
        return [v for v in e if e[v]==r]

    def diameter(self):
        """
        Returns the largest distance between any two vertices. Returns
        Infinity if the (di)graph is not connected.

        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.diameter()
            2
            sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
            sage: G.diameter()
            +Infinity

        Although max( ) is usually defined as -Infinity, since the diameter
        will never be negative, we define it to be zero::

            sage: G = graphs.EmptyGraph()
            sage: G.diameter()
            0

        """

        if self.order() > 0:
            return max(self.eccentricity())
        else:
            return 0

    def distance_graph(self, dist):
        r"""
        Returns the graph on the same vertex set as
        the original graph but vertices are adjacent
        in the returned graph if and only if they are
        at specified distances in the original graph.

        INPUT:

        - ``dist`` is a nonnegative integer or
          a list of nonnegative integers.
          ``Infinity`` may be used here to describe
          vertex pairs in separate components.

        OUTPUT:

        The returned value is an undirected graph.  The
        vertex set is identical to the calling graph, but edges
        of the returned graph join vertices whose distance in
        the calling graph are present in the input ``dist``.
        Loops will only be present if distance 0 is included.  If
        the original graph has a position dictionary specifying
        locations of vertices for plotting, then this information
        is copied over to the distance graph.  In some instances
        this layout may not be the best, and might even be confusing
        when edges run on top of each other due to symmetries
        chosen for the layout.

        EXAMPLES::

            sage: G = graphs.CompleteGraph(3)
            sage: H = G.cartesian_product(graphs.CompleteGraph(2))
            sage: K = H.distance_graph(2)
            sage: K.am()
            [0 0 0 1 0 1]
            [0 0 1 0 1 0]
            [0 1 0 0 0 1]
            [1 0 0 0 1 0]
            [0 1 0 1 0 0]
            [1 0 1 0 0 0]

        To obtain the graph where vertices are adjacent if their
        distance apart is ``d`` or less use a ``range()`` command
        to create the input, using ``d+1`` as the input to ``range``.
        Notice that this will include distance 0 and hence place a loop
        at each vertex.  To avoid this, use ``range(1,d+1)``. ::

            sage: G = graphs.OddGraph(4)
            sage: d = G.diameter()
            sage: n = G.num_verts()
            sage: H = G.distance_graph(range(d+1))
            sage: H.is_isomorphic(graphs.CompleteGraph(n))
            False
            sage: H = G.distance_graph(range(1,d+1))
            sage: H.is_isomorphic(graphs.CompleteGraph(n))
            True

        A complete collection of distance graphs will have
        adjacency matrices that sum to the matrix of all ones. ::

            sage: P = graphs.PathGraph(20)
            sage: all_ones = sum([P.distance_graph(i).am() for i in range(20)])
            sage: all_ones == matrix(ZZ, 20, 20, [1]*400)
            True

        Four-bit strings differing in one bit is the same as
        four-bit strings differing in three bits.  ::

            sage: G = graphs.CubeGraph(4)
            sage: H = G.distance_graph(3)
            sage: G.is_isomorphic(H)
            True

        The graph of eight-bit strings, adjacent if different
        in an odd number of bits.  ::

            sage: G = graphs.CubeGraph(8) # long time
            sage: H = G.distance_graph([1,3,5,7]) # long time
            sage: degrees = [0]*sum([binomial(8,j) for j in [1,3,5,7]]) # long time
            sage: degrees.append(2^8) # long time
            sage: degrees == H.degree_histogram() # long time
            True

        An example of using ``Infinity`` as the distance in
        a graph that is not connected. ::

            sage: G = graphs.CompleteGraph(3)
            sage: H = G.disjoint_union(graphs.CompleteGraph(2))
            sage: L = H.distance_graph(Infinity)
            sage: L.am()
            [0 0 0 1 1]
            [0 0 0 1 1]
            [0 0 0 1 1]
            [1 1 1 0 0]
            [1 1 1 0 0]

        TESTS:

        Empty input, or unachievable distances silently yield empty graphs. ::

            sage: G = graphs.CompleteGraph(5)
            sage: G.distance_graph([]).num_edges()
            0
            sage: G = graphs.CompleteGraph(5)
            sage: G.distance_graph(23).num_edges()
            0

        It is an error to provide a distance that is not an integer type. ::

            sage: G = graphs.CompleteGraph(5)
            sage: G.distance_graph('junk')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=junk) to an integer

        It is an error to provide a negative distance. ::

            sage: G = graphs.CompleteGraph(5)
            sage: G.distance_graph(-3)
            Traceback (most recent call last):
            ...
            ValueError: Distance graph for a negative distance (d=-3) is not defined

        AUTHOR:

        Rob Beezer, 2009-11-25
        """
        from sage.rings.infinity import Infinity
        from copy import copy
        # If input is not a list, make a list with this single object
        if not isinstance(dist, list):
            dist = [dist]
        # Create a list of positive integer (or infinite) distances
        distances = []
        for d in dist:
            if d == Infinity:
                distances.append(d)
            else:
                dint = ZZ(d)
                if dint < 0:
                    raise ValueError('Distance graph for a negative distance (d=%d) is not defined' % dint)
                distances.append(dint)
        # Build a graph on the same vertex set, with loops for distance 0
        vertices = {}
        for v in self.vertex_iterator():
            vertices[v] = {}
        positions = copy(self.get_pos())
        if ZZ(0) in distances:
            looped = True
        else:
            looped = False
        from sage.graphs.all import Graph
        D = Graph(vertices, pos=positions, multiedges=False, loops=looped)
        if len(distances) == 1:
            dstring = "distance " + str(distances[0])
        else:
            dstring = "distances " + str(sorted(distances))
        D.name("Distance graph for %s in " % dstring + self.name())

        # Create the appropriate edges
        d = self.distance_all_pairs()
        for u in self.vertex_iterator():
            for v in self.vertex_iterator():
                if d[u][v] in distances:
                    D.add_edge(u,v)
        return D

    def girth(self):
        """
        Computes the girth of the graph. For directed graphs, computes the
        girth of the undirected graph.

        The girth is the length of the shortest cycle in the graph. Graphs
        without cycles have infinite girth.

        EXAMPLES::

            sage: graphs.TetrahedralGraph().girth()
            3
            sage: graphs.CubeGraph(3).girth()
            4
            sage: graphs.PetersenGraph().girth()
            5
            sage: graphs.HeawoodGraph().girth()
            6
            sage: graphs.trees(9).next().girth()
            +Infinity


        .. SEEALSO::

            * :meth:`~sage.graphs.graph.Graph.odd_girth` -- computes
              the odd girth of a graph.

        TESTS:

        Prior to :trac:`12243`, the girth computation assumed
        vertices were integers (and failed).  The example below
        tests the computation for graphs with vertices that are
        not integers.  In this example the vertices are sets.  ::

            sage: G = graphs.OddGraph(3)
            sage: type(G.vertices()[0])
            <class 'sage.sets.set.Set_object_enumerated_with_category'>
            sage: G.girth()
            5

        Ticket :trac:`12355`::

            sage: H=Graph([(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (1, 6), (2, 5), (3, 4), (5, 6)])
            sage: H.girth()
            3

        Girth < 3 (see :trac:`12355`)::

           sage: g = graphs.PetersenGraph()
           sage: g.allow_multiple_edges(True)
           sage: g.allow_loops(True)
           sage: g.girth()
           5
           sage: g.add_edge(0,0)
           sage: g.girth()
           1
           sage: g.delete_edge(0,0)
           sage: g.add_edge(0,1)
           sage: g.girth()
           2
           sage: g.delete_edge(0,1)
           sage: g.girth()
           5
           sage: g = DiGraph(g)
           sage: g.girth()
           2
        """

        # Cases where girth <= 2
        if self.has_loops():
            return 1
        if self.is_directed():
            if any(self.has_edge(v,u) for u,v in self.edges(labels = False)):
                return 2
        else:
            if self.has_multiple_edges():
                return 2

        n = self.num_verts()
        best = n+1
        seen = {}
        for w in self.vertex_iterator():
            seen[w] = None
            span = set([w])
            depth = 1
            thisList = set([w])
            while 2*depth <= best and 3 < best:
                nextList = set()
                for v in thisList:
                    for u in self.neighbors(v):
                        if u in seen: continue
                        if not u in span:
                            span.add(u)
                            nextList.add(u)
                        else:
                            if u in thisList:
                                best = depth*2-1
                                break
                            if u in nextList:
                                best = depth*2
                    if best == 2*depth-1:
                        break
                thisList = nextList
                depth += 1
        if best == n+1:
            from sage.rings.infinity import Infinity
            return Infinity
        return best



    def periphery(self):
        """
        Returns the set of vertices in the periphery, i.e. whose
        eccentricity is equal to the diameter of the (di)graph.

        In other words, the periphery is the set of vertices achieving the
        maximum eccentricity.

        EXAMPLES::

            sage: G = graphs.DiamondGraph()
            sage: G.periphery()
            [0, 3]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.periphery()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.periphery()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: G = Graph()
            sage: G.periphery()
            []
            sage: G.add_vertex()
            0
            sage: G.periphery()
            [0]
        """
        e = self.eccentricity(with_labels=True)
        try:
            r = max(e.values())
        except Exception:
            return []
        return [v for v in e if e[v]==r]

    ### Paths

    def interior_paths(self, start, end):
        """
        Returns an exhaustive list of paths (also lists) through only
        interior vertices from vertex start to vertex end in the
        (di)graph.

        Note - start and end do not necessarily have to be boundary
        vertices.

        INPUT:


        -  ``start`` - the vertex of the graph to search for
           paths from

        -  ``end`` - the vertex of the graph to search for
           paths to


        EXAMPLES::

            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: sorted(eg1.all_paths(0,6))
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = copy(eg1)
            sage: eg2.set_boundary([0,1,3])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: sorted(eg2.interior_paths(0,6))
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            [[0, 2, 4, 5, 6]]
            sage: sorted(eg2.all_paths(0,6))
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg3 = graphs.PetersenGraph()
            sage: eg3.set_boundary([0,1,2,3,4])
            sage: sorted(eg3.all_paths(1,4))
            [[1, 0, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
            sage: sorted(eg3.interior_paths(1,4))
            [[1, 6, 8, 5, 7, 9, 4], [1, 6, 9, 4]]
            sage: dg = DiGraph({0:[1,3,4], 1:[3], 2:[0,3,4],4:[3]}, boundary=[4])
            sage: sorted(dg.all_paths(0,3))
            [[0, 1, 3], [0, 3], [0, 4, 3]]
            sage: sorted(dg.interior_paths(0,3))
            [[0, 1, 3], [0, 3]]
            sage: ug = dg.to_undirected()
            sage: sorted(ug.all_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 2, 4, 3], [0, 3], [0, 4, 2, 3], [0, 4, 3]]
            sage: sorted(ug.interior_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 3]]
        """
        from sage.misc.superseded import deprecation
        deprecation(15494, "The boundary parameter is deprecated and will soon disappear.")

        from copy import copy
        H = copy(self)
        for vertex in self.get_boundary():
            if (vertex != start and vertex != end):
                H.delete_edges(H.edges_incident(vertex))
        return H.all_paths(start, end)

    def all_paths(self, start, end):
        """
        Returns a list of all paths (also lists) between a pair of
        vertices (start, end) in the (di)graph. If ``start`` is the same
        vertex as ``end``, then ``[[start]]`` is returned -- a list
        containing the 1-vertex, 0-edge path "``start``".

        EXAMPLES::

            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: eg1.all_paths(0,6)
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = graphs.PetersenGraph()
            sage: sorted(eg2.all_paths(1,4))
            [[1, 0, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
            sage: dg = DiGraph({0:[1,3], 1:[3], 2:[0,3]})
            sage: sorted(dg.all_paths(0,3))
            [[0, 1, 3], [0, 3]]
            sage: ug = dg.to_undirected()
            sage: sorted(ug.all_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 3]]

        Starting and ending at the same vertex (see :trac:`13006`)::

            sage: graphs.CompleteGraph(4).all_paths(2,2)
            [[2]]
        """
        if self.is_directed():
            iterator=self.neighbor_out_iterator
        else:
            iterator=self.neighbor_iterator
        if start == end:
            return [[start]]
        all_paths = []      # list of
        act_path = []       # the current path
        act_path_iter = []  # the neighbor/successor-iterators of the current path
        done = False
        s=start
        while not done:
            if s==end:      # if path completes, add to list
                all_paths.append(act_path+[s])
            else:
                if s not in act_path:   # we want vertices just once in a path
                    act_path.append(s)  # extend current path
                    act_path_iter.append(iterator(s))  # save the state of the neighbor/successor-iterator of the current vertex
            s=None
            while (s is None) and not done:
                try:
                    s=act_path_iter[-1].next()  # try to get the next neighbor/successor, ...
                except (StopIteration):         # ... if there is none ...
                    act_path.pop()              # ... go one step back
                    act_path_iter.pop()
                if len(act_path)==0:            # there is no other vertex ...
                    done = True                 # ... so we are done
        return all_paths


    def triangles_count(self, algorithm='iter'):
        """
        Returns the number of triangles in the (di)graph.

        For digraphs, we count the number of directed circuit of length 3.

        INPUT:

        - ``algorithm`` -- (default: ``'matrix'``) specifies the algorithm to
          use among:

            - ``'matrix'`` uses the trace of the cube of the adjacency matrix.

            - ``'iter'`` iterates over the pairs of neighbors of each
              vertex. This is faster for sparse graphs.

        EXAMPLES:

        The Petersen graph is triangle free and thus::

            sage: G = graphs.PetersenGraph()
            sage: G.triangles_count()
            0

        Any triple of vertices in the complete graph induces a triangle so we have::

            sage: G = graphs.CompleteGraph(150)
            sage: G.triangles_count() == binomial(150,3)
            True

        The 2-dimensional DeBruijn graph of 2 symbols has 2 directed C3::

            sage: G = digraphs.DeBruijn(2,2)
            sage: G.triangles_count()
            2

        The directed n-cycle is trivially triangle free for n > 3::

            sage: G = digraphs.Circuit(10)
            sage: G.triangles_count()
            0

        TESTS:

        Comparison on algorithms::

            sage: for i in xrange(10): # long test
            ...       G = graphs.RandomBarabasiAlbert(50,2)
            ...       tm = G.triangles_count(algorithm='matrix')
            ...       te = G.triangles_count(algorithm='iter')
            ...       if tm!=te:
            ...          print "That's not good!"

        Asking for an unknown algorithm::

            sage: G = Graph()
            sage: G.triangles_count(algorithm='tip top')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'tip top' not yet implemented. Please contribute.

        """
        if self.is_directed():
            from sage.graphs.digraph_generators import digraphs
            return self.subgraph_search_count(digraphs.Circuit(3)) // 3

        else:
            if algorithm=='iter':
                from sage.combinat.combination import Combinations
                tr = 0
                ggnx = self.networkx_graph()
                for u in ggnx.nodes_iter():
                    tr += sum(ggnx.has_edge(v,w) for v,w in Combinations(ggnx.neighbors(u),2))
                return tr//3

            elif algorithm=='matrix':
                return (self.adjacency_matrix()**3).trace() // 6

            else:
                raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))

    def shortest_path(self, u, v, by_weight=False, bidirectional=True):
        """
        Returns a list of vertices representing some shortest path from u
        to v: if there is no path from u to v, the list is empty.

        INPUT:


        -  ``by_weight`` - if False, uses a breadth first
           search. If True, takes edge weightings into account, using
           Dijkstra's algorithm.

        -  ``bidirectional`` - if True, the algorithm will
           expand vertices from u and v at the same time, making two spheres
           of half the usual radius. This generally doubles the speed
           (consider the total volume in each case).


        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path(4, 9)
            [4, 17, 16, 12, 13, 9]
            sage: D.shortest_path(5, 5)
            [5]
            sage: D.delete_edges(D.edges_incident(13))
            sage: D.shortest_path(13, 4)
            []
            sage: G = Graph( { 0: [1], 1: [2], 2: [3], 3: [4], 4: [0] })
            sage: G.plot(edge_labels=True).show() # long time
            sage: G.shortest_path(0, 3)
            [0, 4, 3]
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse = True)
            sage: G.shortest_path(0, 3, by_weight=True)
            [0, 1, 2, 3]
        """ #         TODO- multiple edges??
        if u == v: # to avoid a NetworkX bug
            return [u]
        import networkx
        if by_weight:
            if bidirectional:
                try:
                    L = self._backend.bidirectional_dijkstra(u,v)
                except AttributeError:
                    try:
                        L = networkx.bidirectional_dijkstra(self.networkx_graph(copy=False), u, v)[1]
                    except Exception:
                        L = False
            else:
                L = networkx.dijkstra_path(self.networkx_graph(copy=False), u, v)
        else:
            if bidirectional:
                # If the graph is a C_graph, use shortest_path from its backend !
                try:
                    L = self._backend.shortest_path(u,v)
                except AttributeError:
                    L = networkx.shortest_path(self.networkx_graph(copy=False), u, v)
            else:
                try:
                    L = networkx.single_source_shortest_path(self.networkx_graph(copy=False), u)[v]
                except Exception:
                    L = False
        if L:
            return L
        else:
            return []

    def shortest_path_length(self, u, v, by_weight=False,
                                         bidirectional=True,
                                         weight_sum=None):
        """
        Returns the minimal length of paths from u to v.

        If there is no path from u to v, returns Infinity.

        INPUT:

        -  ``by_weight`` - if False, uses a breadth first
           search. If True, takes edge weightings into account, using
           Dijkstra's algorithm.

        -  ``bidirectional`` - if True, the algorithm will
           expand vertices from u and v at the same time, making two spheres
           of half the usual radius. This generally doubles the speed
           (consider the total volume in each case).

        -  ``weight_sum`` - if False, returns the number of
           edges in the path. If True, returns the sum of the weights of these
           edges. Default behavior is to have the same value as by_weight.

        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path_length(4, 9)
            5
            sage: D.shortest_path_length(5, 5)
            0
            sage: D.delete_edges(D.edges_incident(13))
            sage: D.shortest_path_length(13, 4)
            +Infinity
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse = True)
            sage: G.plot(edge_labels=True).show() # long time
            sage: G.shortest_path_length(0, 3)
            2
            sage: G.shortest_path_length(0, 3, by_weight=True)
            3
        """
        if weight_sum is None:
            weight_sum = by_weight
        path = self.shortest_path(u, v, by_weight, bidirectional)
        length = len(path) - 1
        if length == -1:
            from sage.rings.infinity import Infinity
            return Infinity
        if weight_sum:
            wt = 0
            for j in range(length):
                wt += self.edge_label(path[j], path[j+1])
            return wt
        else:
            return length

    def shortest_paths(self, u, by_weight=False, cutoff=None):
        """
        Returns a dictionary associating to each vertex v a shortest path from u
        to v, if it exists.

        INPUT:

        -  ``by_weight`` - if False, uses a breadth first
           search. If True, uses Dijkstra's algorithm to find the shortest
           paths by weight.

        -  ``cutoff`` - integer depth to stop search.

           (ignored if ``by_weight == True``)

        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_paths(0)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 19, 3], 4: [0, 19, 3, 4], 5: [0, 1, 2, 6, 5], 6: [0, 1, 2, 6], 7: [0, 1, 8, 7], 8: [0, 1, 8], 9: [0, 10, 9], 10: [0, 10], 11: [0, 10, 11], 12: [0, 10, 11, 12], 13: [0, 10, 9, 13], 14: [0, 1, 8, 7, 14], 15: [0, 19, 18, 17, 16, 15], 16: [0, 19, 18, 17, 16], 17: [0, 19, 18, 17], 18: [0, 19, 18], 19: [0, 19]}

        All these paths are obviously induced graphs::

            sage: all([D.subgraph(p).is_isomorphic(graphs.PathGraph(len(p)) )for p in D.shortest_paths(0).values()])
            True

        ::

            sage: D.shortest_paths(0, cutoff=2)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 19, 3], 8: [0, 1, 8], 9: [0, 10, 9], 10: [0, 10], 11: [0, 10, 11], 18: [0, 19, 18], 19: [0, 19]}
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True)
            sage: G.plot(edge_labels=True).show() # long time
            sage: G.shortest_paths(0, by_weight=True)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 1, 2, 3], 4: [0, 4]}
        """

        if by_weight:
            import networkx
            return networkx.single_source_dijkstra_path(self.networkx_graph(copy=False), u)
        else:
            try:
                return self._backend.shortest_path_all_vertices(u, cutoff)
            except AttributeError:
                return networkx.single_source_shortest_path(self.networkx_graph(copy=False), u, cutoff)

    def shortest_path_lengths(self, u, by_weight=False, weight_sums=None):
        """
        Returns a dictionary of shortest path lengths keyed by targets that
        are connected by a path from u.

        INPUT:


        -  ``by_weight`` - if False, uses a breadth first
           search. If True, takes edge weightings into account, using
           Dijkstra's algorithm.

        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path_lengths(0)
            {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4, 6: 3, 7: 3, 8: 2, 9: 2, 10: 1, 11: 2, 12: 3, 13: 3, 14: 4, 15: 5, 16: 4, 17: 3, 18: 2, 19: 1}
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show() # long time
            sage: G.shortest_path_lengths(0, by_weight=True)
            {0: 0, 1: 1, 2: 2, 3: 3, 4: 2}
        """
        paths = self.shortest_paths(u, by_weight)
        if by_weight:
            weights = {}
            for v in paths:
                wt = 0
                path = paths[v]
                for j in range(len(path) - 1):
                    wt += self.edge_label(path[j], path[j+1])
                weights[v] = wt
            return weights
        else:
            lengths = {}
            for v in paths:
                lengths[v] = len(paths[v]) - 1
            return lengths

    def shortest_path_all_pairs(self, by_weight=False, default_weight=1, algorithm = "auto"):
        """
        Computes a shortest path between each pair of vertices.

        INPUT:


        - ``by_weight`` - Whether to use the labels defined over the edges as
           weights. If ``False`` (default), the distance between `u` and `v` is
           the minimum number of edges of a path from `u` to `v`.

        - ``default_weight`` - (defaults to 1) The default weight to assign
           edges that don't have a weight (i.e., a label).

           Implies ``by_weight == True``.

        - ``algorithm`` -- four options :

           * ``"BFS"`` -- the computation is done through a BFS
             centered on each vertex successively. Only implemented
             when ``default_weight = 1`` and ``by_weight = False``.

           * ``"Floyd-Warshall-Cython"`` -- through the Cython implementation of
             the Floyd-Warshall algorithm.

           * ``"Floyd-Warshall-Python"`` -- through the Python implementation of
             the Floyd-Warshall algorithm.

           * ``"auto"`` -- use the fastest algorithm depending on the input
             (``"BFS"`` if possible, and ``"Floyd-Warshall-Python"`` otherwise)

             This is the default value.

        OUTPUT:

            A tuple ``(dist, pred)``. They are both dicts of dicts. The first
            indicates the length ``dist[u][v]`` of the shortest weighted path
            from `u` to `v`. The second is a compact representation of all the
            paths- it indicates the predecessor ``pred[u][v]`` of `v` in the
            shortest path from `u` to `v`.

        .. NOTE::

            Three different implementations are actually available through this method :

                * BFS (Cython)
                * Floyd-Warshall (Cython)
                * Floyd-Warshall (Python)

            The BFS algorithm is the fastest of the three, then comes the Cython
            implementation of Floyd-Warshall, and last the Python
            implementation. The first two implementations, however, only compute
            distances based on the topological distance (each edge is of weight
            1, or equivalently the length of a path is its number of
            edges). Besides, they do not deal with graphs larger than 65536
            vertices (which already represents 16GB of ram).

        .. NOTE::

           There is a Cython version of this method that is usually
           much faster for large graphs, as most of the time is
           actually spent building the final double
           dictionary. Everything on the subject is to be found in the
           :mod:`~sage.graphs.distances_all_pairs` module.

        EXAMPLES::

            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show() # long time
            sage: dist, pred = G.shortest_path_all_pairs(by_weight = True)
            sage: dist
            {0: {0: 0, 1: 1, 2: 2, 3: 3, 4: 2}, 1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 3}, 2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 3}, 3: {0: 3, 1: 2, 2: 1, 3: 0, 4: 2}, 4: {0: 2, 1: 3, 2: 3, 3: 2, 4: 0}}
            sage: pred
            {0: {0: None, 1: 0, 2: 1, 3: 2, 4: 0}, 1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0}, 2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3}, 3: {0: 1, 1: 2, 2: 3, 3: None, 4: 3}, 4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}}
            sage: pred[0]
            {0: None, 1: 0, 2: 1, 3: 2, 4: 0}

        So for example the shortest weighted path from `0` to `3` is obtained as
        follows. The predecessor of `3` is ``pred[0][3] == 2``, the predecessor
        of `2` is ``pred[0][2] == 1``, and the predecessor of `1` is
        ``pred[0][1] == 0``.

        ::

            sage: G = Graph( { 0: {1:None}, 1: {2:None}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.shortest_path_all_pairs()
            ({0: {0: 0, 1: 1, 2: 2, 3: 2, 4: 1},
            1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 2},
            2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 2},
            3: {0: 2, 1: 2, 2: 1, 3: 0, 4: 1},
            4: {0: 1, 1: 2, 2: 2, 3: 1, 4: 0}},
            {0: {0: None, 1: 0, 2: 1, 3: 4, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})
            sage: G.shortest_path_all_pairs(by_weight = True)
            ({0: {0: 0, 1: 1, 2: 2, 3: 3, 4: 2},
            1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 3},
            2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 3},
            3: {0: 3, 1: 2, 2: 1, 3: 0, 4: 2},
            4: {0: 2, 1: 3, 2: 3, 3: 2, 4: 0}},
            {0: {0: None, 1: 0, 2: 1, 3: 2, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 1, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})
            sage: G.shortest_path_all_pairs(default_weight=200)
            ({0: {0: 0, 1: 200, 2: 5, 3: 4, 4: 2},
            1: {0: 200, 1: 0, 2: 200, 3: 201, 4: 202},
            2: {0: 5, 1: 200, 2: 0, 3: 1, 4: 3},
            3: {0: 4, 1: 201, 2: 1, 3: 0, 4: 2},
            4: {0: 2, 1: 202, 2: 3, 3: 2, 4: 0}},
            {0: {0: None, 1: 0, 2: 3, 3: 4, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 4, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})

        Checking the distances are equal regardless of the algorithm used::

            sage: g = graphs.Grid2dGraph(5,5)
            sage: d1, _ = g.shortest_path_all_pairs(algorithm="BFS")
            sage: d2, _ = g.shortest_path_all_pairs(algorithm="Floyd-Warshall-Cython")
            sage: d3, _ = g.shortest_path_all_pairs(algorithm="Floyd-Warshall-Python")
            sage: d1 == d2 == d3
            True

        Checking a random path is valid ::

            sage: dist, path = g.shortest_path_all_pairs(algorithm="BFS")
            sage: u,v = g.random_vertex(), g.random_vertex()
            sage: p = [v]
            sage: while p[0] is not None:
            ...     p.insert(0,path[u][p[0]])
            sage: len(p) == dist[u][v] + 2
            True

        TESTS:

        Wrong name for ``algorithm``::

            sage: g.shortest_path_all_pairs(algorithm="Bob")
            Traceback (most recent call last):
            ...
            ValueError: The algorithm keyword can only be set to "auto", "BFS", "Floyd-Warshall-Python" or "Floyd-Warshall-Cython"
        """
        if default_weight != 1:
            by_weight = True

        if algorithm == "auto":
            if by_weight is False:
                algorithm = "BFS"
            else:
                algorithm = "Floyd-Warshall-Python"

        if algorithm == "BFS":
            from sage.graphs.distances_all_pairs import distances_and_predecessors_all_pairs
            return distances_and_predecessors_all_pairs(self)

        elif algorithm == "Floyd-Warshall-Cython":
            from sage.graphs.distances_all_pairs import floyd_warshall
            return floyd_warshall(self, distances = True)

        elif algorithm != "Floyd-Warshall-Python":
            raise ValueError("The algorithm keyword can only be set to "+
                             "\"auto\","+
                             " \"BFS\", "+
                             "\"Floyd-Warshall-Python\" or "+
                             "\"Floyd-Warshall-Cython\"")

        from sage.rings.infinity import Infinity
        dist = {}
        pred = {}
        verts = self.vertices()
        for u in verts:
            du = {}
            pu = {}
            for v in verts:
                if self.has_edge(u, v):
                    if by_weight is False:
                        du[v] = 1
                    else:
                        edge_label = self.edge_label(u, v)
                        if edge_label is None or edge_label == {}:
                            du[v] = default_weight
                        else:
                            du[v] = edge_label
                    pu[v] = u
                else:
                    du[v] = Infinity
                    pu[v] = None
            du[u] = 0
            dist[u] = du
            pred[u] = pu

        for w in verts:
            dw = dist[w]
            for u in verts:
                du = dist[u]
                for v in verts:
                    if du[v] > du[w] + dw[v]:
                        du[v] = du[w] + dw[v]
                        pred[u][v] = pred[w][v]

        return dist, pred

    def average_distance(self):
        r"""
        Returns the average distance between vertices of the graph.

        Formally, for a graph `G` this value is equal to
        `\frac 1 {n(n-1)} \sum_{u,v\in G} d(u,v)` where `d(u,v)`
        denotes the distance between vertices `u` and `v` and `n`
        is the number of vertices in `G`.

        EXAMPLE:

        From [GYLL93]_::

            sage: g=graphs.PathGraph(10)
            sage: w=lambda x: (x*(x*x -1)/6)/(x*(x-1)/2)
            sage: g.average_distance()==w(10)
            True

        REFERENCE:

        .. [GYLL93] I. Gutman, Y.-N. Yeh, S.-L. Lee, and Y.-L. Luo. Some recent
          results in the theory of the Wiener number. *Indian Journal of
          Chemistry*, 32A:651--661, 1993.

        TEST::

            sage: g = Graph()
            sage: g.average_distance()
            Traceback (most recent call last):
            ...
            ValueError: The graph must have at least two vertices for this value to be defined
        """
        if self.order() < 2:
            raise ValueError("The graph must have at least two vertices for this value to be defined")

        return Integer(self.wiener_index())/Integer((self.order()*(self.order()-1))/2)

    def szeged_index(self):
        r"""
        Returns the Szeged index of the graph.

        For any `uv\in E(G)`, let
        `N_u(uv) = \{w\in G:d(u,w)<d(v,w)\}, n_u(uv)=|N_u(uv)|`

        The Szeged index of a graph is then defined as [1]:
        `\sum_{uv \in E(G)}n_u(uv)\times n_v(uv)`

        EXAMPLE:

        True for any connected graph [1]::

            sage: g=graphs.PetersenGraph()
            sage: g.wiener_index()<= g.szeged_index()
            True

        True for all trees [1]::

            sage: g=Graph()
            sage: g.add_edges(graphs.CubeGraph(5).min_spanning_tree())
            sage: g.wiener_index() == g.szeged_index()
            True


        REFERENCE:

        [1] Klavzar S., Rajapakse A., Gutman I. (1996). The Szeged and the
        Wiener index of graphs. Applied Mathematics Letters, 9 (5), pp. 45-49.
        """
        distances=self.distance_all_pairs()
        s=0
        for (u,v) in self.edges(labels=None):
            du=distances[u]
            dv=distances[v]
            n1=n2=0
            for w in self:
                if du[w] < dv[w]:
                    n1+=1
                elif dv[w] < du[w]:
                    n2+=1
            s+=(n1*n2)
        return s

    ### Searches

    def breadth_first_search(self, start, ignore_direction=False,
                             distance=None, neighbors=None):
        """
        Returns an iterator over the vertices in a breadth-first ordering.

        INPUT:


        - ``start`` - vertex or list of vertices from which to start
          the traversal

        - ``ignore_direction`` - (default False) only applies to
          directed graphs. If True, searches across edges in either
          direction.

        - ``distance`` - the maximum distance from the ``start`` nodes
          to traverse.  The ``start`` nodes are distance zero from
          themselves.

        - ``neighbors`` - a function giving the neighbors of a vertex.
          The function should take a vertex and return a list of
          vertices.  For a graph, ``neighbors`` is by default the
          :meth:`.neighbors` function of the graph.  For a digraph,
          the ``neighbors`` function defaults to the
          :meth:`.successors` function of the graph.

        .. SEEALSO::

        - :meth:`breadth_first_search <sage.graphs.base.c_graph.CGraphBackend.breadth_first_search>`
          -- breadth-first search for fast compiled graphs.

        - :meth:`depth_first_search <sage.graphs.base.c_graph.CGraphBackend.depth_first_search>`
          -- depth-first search for fast compiled graphs.

        - :meth:`depth_first_search` -- depth-first search for generic graphs.

        EXAMPLES::

            sage: G = Graph( { 0: [1], 1: [2], 2: [3], 3: [4], 4: [0]} )
            sage: list(G.breadth_first_search(0))
            [0, 1, 4, 2, 3]

        By default, the edge direction of a digraph is respected, but this
        can be overridden by the ``ignore_direction`` parameter::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.breadth_first_search(0))
            [0, 1, 2, 3, 4, 5, 6, 7]
            sage: list(D.breadth_first_search(0, ignore_direction=True))
            [0, 1, 2, 3, 7, 4, 5, 6]

        You can specify a maximum distance in which to search.  A
        distance of zero returns the ``start`` vertices::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.breadth_first_search(0,distance=0))
            [0]
            sage: list(D.breadth_first_search(0,distance=1))
            [0, 1, 2, 3]

        Multiple starting vertices can be specified in a list::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.breadth_first_search([0]))
            [0, 1, 2, 3, 4, 5, 6, 7]
            sage: list(D.breadth_first_search([0,6]))
            [0, 6, 1, 2, 3, 7, 4, 5]
            sage: list(D.breadth_first_search([0,6],distance=0))
            [0, 6]
            sage: list(D.breadth_first_search([0,6],distance=1))
            [0, 6, 1, 2, 3, 7]
            sage: list(D.breadth_first_search(6,ignore_direction=True,distance=2))
            [6, 3, 7, 0, 5]

        More generally, you can specify a ``neighbors`` function.  For
        example, you can traverse the graph backwards by setting
        ``neighbors`` to be the :meth:`.neighbors_in` function of the graph::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.breadth_first_search(5,neighbors=D.neighbors_in, distance=2))
            [5, 1, 2, 0]
            sage: list(D.breadth_first_search(5,neighbors=D.neighbors_out, distance=2))
            [5, 7, 0]
            sage: list(D.breadth_first_search(5,neighbors=D.neighbors, distance=2))
            [5, 1, 2, 7, 0, 4, 6]


        TESTS::

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.breadth_first_search(0))
            [0]
            sage: list(D.breadth_first_search(0, ignore_direction=True))
            [0, 1, 2]

        """
        # Preferably use the Cython implementation
        if neighbors is None and not isinstance(start,list) and distance is None and hasattr(self._backend,"breadth_first_search"):
            for v in self._backend.breadth_first_search(start, ignore_direction = ignore_direction):
                yield v
        else:
            if neighbors is None:
                if not self._directed or ignore_direction:
                    neighbors=self.neighbor_iterator
                else:
                    neighbors=self.neighbor_out_iterator
            seen=set([])
            if isinstance(start, list):
                queue=[(v,0) for v in start]
            else:
                queue=[(start,0)]

            for v,d in queue:
                yield v
                seen.add(v)

            while len(queue)>0:
                v,d = queue.pop(0)
                if distance is None or d<distance:
                    for w in neighbors(v):
                        if w not in seen:
                            seen.add(w)
                            queue.append((w, d+1))
                            yield w

    def depth_first_search(self, start, ignore_direction=False,
                           distance=None, neighbors=None):
        """
        Returns an iterator over the vertices in a depth-first ordering.

        INPUT:


        - ``start`` - vertex or list of vertices from which to start
          the traversal

        - ``ignore_direction`` - (default False) only applies to
          directed graphs. If True, searches across edges in either
          direction.

        - ``distance`` - the maximum distance from the ``start`` nodes
          to traverse.  The ``start`` nodes are distance zero from
          themselves.

        - ``neighbors`` - a function giving the neighbors of a vertex.
          The function should take a vertex and return a list of
          vertices.  For a graph, ``neighbors`` is by default the
          :meth:`.neighbors` function of the graph.  For a digraph,
          the ``neighbors`` function defaults to the
          :meth:`.successors` function of the graph.

        .. SEEALSO::

        - :meth:`breadth_first_search`

        - :meth:`breadth_first_search <sage.graphs.base.c_graph.CGraphBackend.breadth_first_search>`
          -- breadth-first search for fast compiled graphs.

        - :meth:`depth_first_search <sage.graphs.base.c_graph.CGraphBackend.depth_first_search>`
          -- depth-first search for fast compiled graphs.

        EXAMPLES::

            sage: G = Graph( { 0: [1], 1: [2], 2: [3], 3: [4], 4: [0]} )
            sage: list(G.depth_first_search(0))
            [0, 4, 3, 2, 1]

        By default, the edge direction of a digraph is respected, but this
        can be overridden by the ``ignore_direction`` parameter::


            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.depth_first_search(0))
            [0, 3, 6, 7, 2, 5, 1, 4]
            sage: list(D.depth_first_search(0, ignore_direction=True))
            [0, 7, 6, 3, 5, 2, 1, 4]

        You can specify a maximum distance in which to search.  A
        distance of zero returns the ``start`` vertices::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.depth_first_search(0,distance=0))
            [0]
            sage: list(D.depth_first_search(0,distance=1))
            [0, 3, 2, 1]

        Multiple starting vertices can be specified in a list::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.depth_first_search([0]))
            [0, 3, 6, 7, 2, 5, 1, 4]
            sage: list(D.depth_first_search([0,6]))
            [0, 3, 6, 7, 2, 5, 1, 4]
            sage: list(D.depth_first_search([0,6],distance=0))
            [0, 6]
            sage: list(D.depth_first_search([0,6],distance=1))
            [0, 3, 2, 1, 6, 7]
            sage: list(D.depth_first_search(6,ignore_direction=True,distance=2))
            [6, 7, 5, 0, 3]

        More generally, you can specify a ``neighbors`` function.  For
        example, you can traverse the graph backwards by setting
        ``neighbors`` to be the :meth:`.neighbors_in` function of the graph::

            sage: D = DiGraph( { 0: [1,2,3], 1: [4,5], 2: [5], 3: [6], 5: [7], 6: [7], 7: [0]})
            sage: list(D.depth_first_search(5,neighbors=D.neighbors_in, distance=2))
            [5, 2, 0, 1]
            sage: list(D.depth_first_search(5,neighbors=D.neighbors_out, distance=2))
            [5, 7, 0]
            sage: list(D.depth_first_search(5,neighbors=D.neighbors, distance=2))
            [5, 7, 6, 0, 2, 1, 4]

        TESTS::

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.depth_first_search(0))
            [0]
            sage: list(D.depth_first_search(0, ignore_direction=True))
            [0, 2, 1]

        """
        # Preferably use the Cython implementation
        if neighbors is None and not isinstance(start,list) and  distance is None and hasattr(self._backend,"depth_first_search"):
            for v in self._backend.depth_first_search(start, ignore_direction = ignore_direction):
                yield v
        else:
            if neighbors is None:
                if not self._directed or ignore_direction:
                    neighbors=self.neighbor_iterator
                else:
                    neighbors=self.neighbor_out_iterator
            seen=set([])
            if isinstance(start, list):
                # Reverse the list so that the initial vertices come out in the same order
                queue=[(v,0) for v in reversed(start)]
            else:
                queue=[(start,0)]

            while len(queue)>0:
                v,d = queue.pop()
                if v not in seen:
                    yield v
                    seen.add(v)
                    if distance is None or d<distance:
                        for w in neighbors(v):
                            if w not in seen:
                                queue.append((w, d+1))

    def lex_BFS(self,reverse=False,tree=False, initial_vertex = None):
        r"""
        Performs a Lex BFS on the graph.

        A Lex BFS ( or Lexicographic Breadth-First Search ) is a Breadth
        First Search used for the recognition of Chordal Graphs. For more
        information, see the
        `Wikipedia article on Lex-BFS
        <http://en.wikipedia.org/wiki/Lexicographic_breadth-first_search>`_.

        INPUT:

        - ``reverse`` (boolean) -- whether to return the vertices
          in discovery order, or the reverse.

          ``False`` by default.

        - ``tree`` (boolean) -- whether to return the discovery
          directed tree (each vertex being linked to the one that
          saw it for the first time)

          ``False`` by default.

        - ``initial_vertex`` -- the first vertex to consider.

          ``None`` by default.

        ALGORITHM:

        This algorithm maintains for each vertex left in the graph
        a code corresponding to the vertices already removed. The
        vertex of maximal code ( according to the lexicographic
        order ) is then removed, and the codes are updated.

        This algorithm runs in time `O(n^2)` ( where `n` is the
        number of vertices in the graph ), which is not optimal.
        An optimal algorithm would run in time `O(m)` ( where `m`
        is the number of edges in the graph ), and require the use
        of a doubly-linked list which are not available in python
        and can not really be written efficiently. This could be
        done in Cython, though.

        EXAMPLE:

        A Lex BFS is obviously an ordering of the vertices::

            sage: g = graphs.PetersenGraph()
            sage: len(g.lex_BFS()) == g.order()
            True

        For a Chordal Graph, a reversed Lex BFS is a Perfect
        Elimination Order ::

            sage: g = graphs.PathGraph(3).lexicographic_product(graphs.CompleteGraph(2))
            sage: g.lex_BFS(reverse=True)
            [(2, 1), (2, 0), (1, 1), (1, 0), (0, 1), (0, 0)]


        And the vertices at the end of the tree of discovery are, for
        chordal graphs, simplicial vertices (their neighborhood is
        a complete graph)::

            sage: g = graphs.ClawGraph().lexicographic_product(graphs.CompleteGraph(2))
            sage: v = g.lex_BFS()[-1]
            sage: peo, tree = g.lex_BFS(initial_vertex = v,  tree=True)
            sage: leaves = [v for v in tree if tree.in_degree(v) ==0]
            sage: all([g.subgraph(g.neighbors(v)).is_clique() for v in leaves])
            True

        TESTS:

        There were some problems with the following call in the past (trac 10899) -- now
        it should be fine::

            sage: Graph(1).lex_BFS(tree=True)
            ([0], Digraph on 1 vertex)

        """
        id_inv = dict([(i,v) for (v,i) in zip(self.vertices(),range(self.order()))])
        code = [[] for i in range(self.order())]
        m = self.am()

        l = lambda x : code[x]
        vertices = set(range(self.order()))

        value = []
        pred = [-1]*self.order()

        add_element = (lambda y:value.append(id_inv[y])) if not reverse else (lambda y: value.insert(0,id_inv[y]))

        # Should we take care of the first vertex we pick ?
        first = True if initial_vertex is not None else False


        while vertices:

            if not first:
                v = max(vertices,key=l)
            else:
                v = self.vertices().index(initial_vertex)
                first = False

            vertices.remove(v)
            vector = m.column(v)
            for i in vertices:
                code[i].append(vector[i])
                if vector[i]:
                    pred[i] = v
            add_element(v)

        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            g.add_vertices(self.vertices())
            edges = [(id_inv[i], id_inv[pred[i]]) for i in range(self.order()) if pred[i]!=-1]
            g.add_edges(edges)
            return value, g

        else:
            return value

    ### Constructors

    def add_cycle(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices
        are already present, only the edges are added.

        For digraphs, adds the directed cycle, whose orientation is
        determined by the list. Adds edges (vertices[u], vertices[u+1]) and
        (vertices[-1], vertices[0]).

        INPUT:

        -  ``vertices`` -- a list of indices for the vertices of
           the cycle to be added.


        EXAMPLES::

            sage: G = Graph()
            sage: G.add_vertices(range(10)); G
            Graph on 10 vertices
            sage: show(G)
            sage: G.add_cycle(range(20)[10:20])
            sage: show(G)
            sage: G.add_cycle(range(10))
            sage: show(G)

        ::

            sage: D = DiGraph()
            sage: D.add_cycle(range(4))
            sage: D.edges()
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 0, None)]
        """
        self.add_path(vertices)
        self.add_edge(vertices[-1], vertices[0])

    def add_path(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices
        are already present, only the edges are added.

        For digraphs, adds the directed path vertices[0], ...,
        vertices[-1].

        INPUT:


        -  ``vertices`` - a list of indices for the vertices of
           the cycle to be added.


        EXAMPLES::

            sage: G = Graph()
            sage: G.add_vertices(range(10)); G
            Graph on 10 vertices
            sage: show(G)
            sage: G.add_path(range(20)[10:20])
            sage: show(G)
            sage: G.add_path(range(10))
            sage: show(G)

        ::

            sage: D = DiGraph()
            sage: D.add_path(range(4))
            sage: D.edges()
            [(0, 1, None), (1, 2, None), (2, 3, None)]
        """
        vert1 = vertices[0]
        for v in vertices[1:]:
            self.add_edge(vert1, v)
            vert1 = v

    def complement(self):
        """
        Returns the complement of the (di)graph.

        The complement of a graph has the same vertices, but exactly those
        edges that are not in the original graph. This is not well defined
        for graphs with multiple edges.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.plot() # long time
            sage: PC = P.complement()
            sage: PC.plot() # long time

        ::

            sage: graphs.TetrahedralGraph().complement().size()
            0
            sage: graphs.CycleGraph(4).complement().edges()
            [(0, 2, None), (1, 3, None)]
            sage: graphs.CycleGraph(4).complement()
            complement(Cycle graph): Graph on 4 vertices
            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edges([(0,1)]*3)
            sage: G.complement()
            Traceback (most recent call last):
            ...
            TypeError: complement not well defined for (di)graphs with multiple edges

        TESTS:

        We check that :trac:`15699` is fixed::

            sage: G = graphs.PathGraph(5).copy(immutable=True)
            sage: G.complement()
            complement(Path Graph): Graph on 5 vertices
        """
        if self.has_multiple_edges():
            raise TypeError('complement not well defined for (di)graphs with multiple edges')
        self._scream_if_not_simple()
        G = self.copy(immutable=False) # Make sure it's a mutable copy
        G.delete_edges(G.edges())
        G.name('complement(%s)'%self.name())
        for u in self:
            for v in self:
                if not self.has_edge(u,v):
                    G.add_edge(u,v)
        if getattr(self, '_immutable', False):
            return G.copy(immutable=True)
        return G

    def to_simple(self):
        """
        Returns a simple version of itself (i.e., undirected and loops and
        multiple edges are removed).

        EXAMPLES::

            sage: G = DiGraph(loops=True,multiedges=True,sparse=True)
            sage: G.add_edges( [ (0,0,None), (1,1,None), (2,2,None), (2,3,1), (2,3,2), (3,2,None) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (2, 3), (3, 2)]
            sage: H=G.to_simple()
            sage: H.edges(labels=False)
            [(2, 3)]
            sage: H.is_directed()
            False
            sage: H.allows_loops()
            False
            sage: H.allows_multiple_edges()
            False
        """
        g=self.to_undirected()
        g.allow_loops(False)
        g.allow_multiple_edges(False)
        return g

    def disjoint_union(self, other, verbose_relabel=True):
        """
        Returns the disjoint union of self and other.

        INPUT:

        - ``verbose_relabel`` - (defaults to True) If True, each
          vertex v in the first graph will be named '0,v' and each
          vertex u in the second graph will be named '1,u' in the
          final graph. If False, the vertices of the first graph and
          the second graph will be relabeled with consecutive
          integers.

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.union`

            * :meth:`~Graph.join`

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G.disjoint_union(H); J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3)]
            sage: J = G.disjoint_union(H, verbose_relabel=False); J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4, 5, 6]

        ::

            sage: G=Graph({'a': ['b']})
            sage: G.name("Custom path")
            sage: G.name()
            'Custom path'
            sage: H=graphs.CycleGraph(3)
            sage: J=G.disjoint_union(H); J
            Custom path disjoint_union Cycle graph: Graph on 5 vertices
            sage: J.vertices()
            [(0, 'a'), (0, 'b'), (1, 0), (1, 1), (1, 2)]

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('both arguments must be of the same class')

        if not verbose_relabel:
            r_self = {}; r_other = {}; i = 0
            for v in self:
                r_self[v] = i; i += 1
            for v in other:
                r_other[v] = i; i += 1
            G = self.relabel(r_self, inplace=False).union(other.relabel(r_other, inplace=False))
        else:
            r_self = dict([[v,(0,v)] for v in self])
            r_other = dict([[v,(1,v)] for v in other])
            G = self.relabel(r_self, inplace=False).union(other.relabel(r_other, inplace=False))

        G._name = '%s disjoint_union %s'%(self.name(), other.name())
        return G

    def union(self, other):
        """
        Returns the union of self and other.

        If the graphs have common vertices, the common vertices will be
        identified.

        If one of the two graphs allows loops (or multiple edges), the resulting
        graph will allow loops (or multiple edges).

        If both graphs are immutable, the resulting graph is immutable.

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.disjoint_union`

            * :meth:`~Graph.join`

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G.union(H); J
            Graph on 4 vertices
            sage: J.vertices()
            [0, 1, 2, 3]
            sage: J.edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]

        TESTS:

        Multiple edges and loops (:trac:`15627`)::

            sage: g = Graph(multiedges=True, loops=True)
            sage: g.add_edges(graphs.PetersenGraph().edges())
            sage: g.add_edges(graphs.PetersenGraph().edges())
            sage: g.add_edge(0,0)
            sage: g.add_edge(0,0,"Hey")
            sage: g.add_edge(0,9)
            sage: g.add_edge(0,9)
            sage: g.add_edge(0,9)
            sage: (2*g.size()) == (2*g).size()
            True

        Immutable input ? Immutable output (:trac:`15627`)::

            sage: g = g.copy(immutable=True)
            sage: (2*g)._backend
            <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('both arguments must be of the same class')

        multiedges = self.allows_multiple_edges() or other.allows_multiple_edges()
        loops      = self.allows_loops()          or other.allows_loops()

        if self._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph(multiedges=multiedges, loops=loops)
        else:
            from sage.graphs.all import Graph
            G = Graph(multiedges=multiedges, loops=loops)
        G.add_vertices(self.vertices())
        G.add_vertices(other.vertices())
        G.add_edges(self.edges())
        G.add_edges(other.edges())

        if getattr(self, "_immutable", False) and getattr(other, "_immutable", False):
            G = G.copy(immutable=True)

        return G

    def cartesian_product(self, other):
        r"""
        Returns the Cartesian product of self and other.

        The Cartesian product of `G` and `H` is the graph `L` with vertex set
        `V(L)` equal to the Cartesian product of the vertices `V(G)` and `V(H)`,
        and `((u,v), (w,x))` is an edge iff either - `(u, w)` is an edge of self
        and `v = x`, or - `(v, x)` is an edge of other and `u = w`.

        .. SEEALSO::

            - :meth:`~sage.graphs.graph_decompositions.graph_products.is_cartesian_product`
              -- factorization of graphs according to the cartesian product

            - :mod:`~sage.graphs.graph_decompositions.graph_products`
              -- a module on graph products.

        TESTS:

        Cartesian product of graphs::

            sage: G = Graph([(0,1),(1,2)])
            sage: H = Graph([('a','b')])
            sage: C1 = G.cartesian_product(H)
            sage: C1.edges(labels=None)
            [((0, 'a'), (0, 'b')), ((0, 'a'), (1, 'a')), ((0, 'b'), (1, 'b')), ((1, 'a'), (1, 'b')), ((1, 'a'), (2, 'a')), ((1, 'b'), (2, 'b')), ((2, 'a'), (2, 'b'))]
            sage: C2 = H.cartesian_product(G)
            sage: C1.is_isomorphic(C2)
            True

        Construction of a Toroidal grid::

            sage: A = graphs.CycleGraph(3)
            sage: B = graphs.CycleGraph(4)
            sage: T = A.cartesian_product(B)
            sage: T.is_isomorphic( graphs.ToroidalGrid2dGraph(3,4) )
            True

        Cartesian product of digraphs::

            sage: P = DiGraph([(0,1)])
            sage: B = digraphs.DeBruijn( ['a','b'], 2 )
            sage: Q = P.cartesian_product(B)
            sage: Q.edges(labels=None)
            [((0, 'aa'), (0, 'aa')), ((0, 'aa'), (0, 'ab')), ((0, 'aa'), (1, 'aa')), ((0, 'ab'), (0, 'ba')), ((0, 'ab'), (0, 'bb')), ((0, 'ab'), (1, 'ab')), ((0, 'ba'), (0, 'aa')), ((0, 'ba'), (0, 'ab')), ((0, 'ba'), (1, 'ba')), ((0, 'bb'), (0, 'ba')), ((0, 'bb'), (0, 'bb')), ((0, 'bb'), (1, 'bb')), ((1, 'aa'), (1, 'aa')), ((1, 'aa'), (1, 'ab')), ((1, 'ab'), (1, 'ba')), ((1, 'ab'), (1, 'bb')), ((1, 'ba'), (1, 'aa')), ((1, 'ba'), (1, 'ab')), ((1, 'bb'), (1, 'ba')), ((1, 'bb'), (1, 'bb'))]
            sage: Q.strongly_connected_components_digraph().num_verts()
            2
            sage: V = Q.strongly_connected_component_containing_vertex( (0, 'aa') )
            sage: B.is_isomorphic( Q.subgraph(V) )
            True
        """
        self._scream_if_not_simple(allow_loops=True)
        if self._directed and other._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph( loops = (self.has_loops() or other.has_loops()) )
        elif (not self._directed) and (not other._directed):
            from sage.graphs.all import Graph
            G = Graph()
        else:
            raise TypeError('the graphs should be both directed or both undirected')

        G.add_vertices( [(u,v) for u in self for v in other] )
        for u,w in self.edge_iterator(labels=None):
            for v in other:
                G.add_edge((u,v), (w,v))
        for v,x in other.edge_iterator(labels=None):
            for u in self:
                G.add_edge((u,v), (u,x))
        return G

    def tensor_product(self, other):
        r"""
        Returns the tensor product of self and other.

        The tensor product of `G` and `H` is the graph `L` with vertex set
        `V(L)` equal to the Cartesian product of the vertices `V(G)` and `V(H)`,
        and `((u,v), (w,x))` is an edge iff - `(u, w)` is an edge of self, and -
        `(v, x)` is an edge of other.

        The tensor product is also known as the categorical product and the
        kronecker product (refering to the kronecker matrix product). See
        :wikipedia:`Wikipedia article on the Kronecker product <Kronecker_product>`.

        EXAMPLES::

            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: T = C.tensor_product(Z); T
            Graph on 10 vertices
            sage: T.size()
            10
            sage: T.plot() # long time

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: T = D.tensor_product(P); T
            Graph on 200 vertices
            sage: T.size()
            900
            sage: T.plot() # long time

        TESTS:

        Tensor product of graphs::

            sage: G = Graph([(0,1), (1,2)])
            sage: H = Graph([('a','b')])
            sage: T = G.tensor_product(H)
            sage: T.edges(labels=None)
            [((0, 'a'), (1, 'b')), ((0, 'b'), (1, 'a')), ((1, 'a'), (2, 'b')), ((1, 'b'), (2, 'a'))]
            sage: T.is_isomorphic( H.tensor_product(G) )
            True

        Tensor product of digraphs::

            sage: I = DiGraph([(0,1), (1,2)])
            sage: J = DiGraph([('a','b')])
            sage: T = I.tensor_product(J)
            sage: T.edges(labels=None)
            [((0, 'a'), (1, 'b')), ((1, 'a'), (2, 'b'))]
            sage: T.is_isomorphic( J.tensor_product(I) )
            True

        The tensor product of two DeBruijn digraphs of same diameter is a DeBruijn digraph::

            sage: B1 = digraphs.DeBruijn(2, 3)
            sage: B2 = digraphs.DeBruijn(3, 3)
            sage: T = B1.tensor_product( B2 )
            sage: T.is_isomorphic( digraphs.DeBruijn( 2*3, 3) )
            True
        """
        self._scream_if_not_simple(allow_loops=True)
        if self._directed and other._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph( loops = (self.has_loops() or other.has_loops()) )
        elif (not self._directed) and (not other._directed):
            from sage.graphs.all import Graph
            G = Graph()
        else:
            raise TypeError('the graphs should be both directed or both undirected')
        G.add_vertices( [(u, v) for u in self for v in other] )
        for u, w in self.edge_iterator(labels=None):
            for v, x in other.edge_iterator(labels=None):
                G.add_edge((u, v), (w, x))
                if not G._directed:
                    G.add_edge((u, x), (w, v))
        return G

    categorical_product = tensor_product
    kronecker_product = tensor_product

    def lexicographic_product(self, other):
        r"""
        Returns the lexicographic product of self and other.

        The lexicographic product of `G` and `H` is the graph `L` with vertex
        set `V(L)=V(G)\times V(H)`, and `((u,v), (w,x))` is an edge iff :

        * `(u, w)` is an edge of `G`, or
        * `u = w` and `(v, x)` is an edge of `H`.

        EXAMPLES::

            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: L = C.lexicographic_product(Z); L
            Graph on 10 vertices
            sage: L.plot() # long time

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: L = D.lexicographic_product(P); L
            Graph on 200 vertices
            sage: L.plot() # long time

        TESTS:

        Lexicographic product of graphs::

            sage: G = Graph([(0,1), (1,2)])
            sage: H = Graph([('a','b')])
            sage: T = G.lexicographic_product(H)
            sage: T.edges(labels=None)
            [((0, 'a'), (0, 'b')), ((0, 'a'), (1, 'a')), ((0, 'a'), (1, 'b')), ((0, 'b'), (1, 'a')), ((0, 'b'), (1, 'b')), ((1, 'a'), (1, 'b')), ((1, 'a'), (2, 'a')), ((1, 'a'), (2, 'b')), ((1, 'b'), (2, 'a')), ((1, 'b'), (2, 'b')), ((2, 'a'), (2, 'b'))]
            sage: T.is_isomorphic( H.lexicographic_product(G) )
            False

        Lexicographic product of digraphs::

            sage: I = DiGraph([(0,1), (1,2)])
            sage: J = DiGraph([('a','b')])
            sage: T = I.lexicographic_product(J)
            sage: T.edges(labels=None)
            [((0, 'a'), (0, 'b')), ((0, 'a'), (1, 'a')), ((0, 'a'), (1, 'b')), ((0, 'b'), (1, 'a')), ((0, 'b'), (1, 'b')), ((1, 'a'), (1, 'b')), ((1, 'a'), (2, 'a')), ((1, 'a'), (2, 'b')), ((1, 'b'), (2, 'a')), ((1, 'b'), (2, 'b')), ((2, 'a'), (2, 'b'))]
            sage: T.is_isomorphic( J.lexicographic_product(I) )
            False
        """
        self._scream_if_not_simple(allow_loops=True)
        if self._directed and other._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph( loops = (self.has_loops() or other.has_loops()) )
        elif (not self._directed) and (not other._directed):
            from sage.graphs.all import Graph
            G = Graph()
        else:
            raise TypeError('the graphs should be both directed or both undirected')
        G.add_vertices( [(u,v) for u in self for v in other] )
        for u,w in self.edge_iterator(labels=None):
            for v in other:
                for x in other:
                    G.add_edge((u,v), (w,x))
        for u in self:
            for v,x in other.edge_iterator(labels=None):
                G.add_edge((u,v), (u,x))
        return G

    def strong_product(self, other):
        r"""
        Returns the strong product of self and other.

        The strong product of `G` and `H` is the graph `L` with vertex set
        `V(L)=V(G)\times V(H)`, and `((u,v), (w,x))` is an edge of `L` iff
        either :

        * `(u, w)` is an edge of `G` and `v = x`, or
        * `(v, x)` is an edge of `H` and `u = w`, or
        * `(u, w)` is an edge of `G` and `(v, x)` is an edge of `H`.

        In other words, the edges of the strong product is the union of the
        edges of the tensor and Cartesian products.

        EXAMPLES::

            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: S = C.strong_product(Z); S
            Graph on 10 vertices
            sage: S.plot() # long time

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: S = D.strong_product(P); S
            Graph on 200 vertices
            sage: S.plot() # long time

        TESTS:

        Strong product of graphs is commutative::

            sage: G = Graph([(0,1), (1,2)])
            sage: H = Graph([('a','b')])
            sage: T = G.strong_product(H)
            sage: T.is_isomorphic( H.strong_product(G) )
            True

        Strong product of digraphs is commutative::

            sage: I = DiGraph([(0,1), (1,2)])
            sage: J = DiGraph([('a','b')])
            sage: T = I.strong_product(J)
            sage: T.is_isomorphic( J.strong_product(I) )
            True

        Counting the edges (see :trac:`13699`)::

            sage: g = graphs.RandomGNP(5,.5)
            sage: gn,gm = g.order(), g.size()
            sage: h = graphs.RandomGNP(5,.5)
            sage: hn,hm = h.order(), h.size()
            sage: product_size = g.strong_product(h).size()
            sage: expected = gm*hn + hm*gn + 2*gm*hm
            sage: if product_size != expected:
            ...       print "Something is really wrong here...", product_size, "!=", expected
        """
        self._scream_if_not_simple(allow_loops=True)
        if self._directed and other._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph( loops = (self.has_loops() or other.has_loops()) )
        elif (not self._directed) and (not other._directed):
            from sage.graphs.all import Graph
            G = Graph()
        else:
            raise TypeError('the graphs should be both directed or both undirected')

        G.add_vertices( [(u,v) for u in self for v in other] )
        for u,w in self.edge_iterator(labels=None):
            for v in other:
                G.add_edge((u,v), (w,v))
            for v,x in other.edge_iterator(labels=None):
                G.add_edge((u,v), (w,x))
                if not self._directed:
                    G.add_edge((w,v), (u,x))
        for v,x in other.edge_iterator(labels=None):
            for u in self:
                G.add_edge((u,v), (u,x))
        return G

    def disjunctive_product(self, other):
        r"""
        Returns the disjunctive product of self and other.

        The disjunctive product of `G` and `H` is the graph `L` with vertex set
        `V(L)=V(G)\times V(H)`, and `((u,v), (w,x))` is an edge iff either :

        * `(u, w)` is an edge of `G`, or
        * `(v, x)` is an edge of `H`.

        EXAMPLES::

            sage: Z = graphs.CompleteGraph(2)
            sage: D = Z.disjunctive_product(Z); D
            Graph on 4 vertices
            sage: D.plot() # long time

        ::

            sage: C = graphs.CycleGraph(5)
            sage: D = C.disjunctive_product(Z); D
            Graph on 10 vertices
            sage: D.plot() # long time

        TESTS:

        Disjunctive product of graphs::

            sage: G = Graph([(0,1), (1,2)])
            sage: H = Graph([('a','b')])
            sage: T = G.disjunctive_product(H)
            sage: T.edges(labels=None)
            [((0, 'a'), (0, 'b')), ((0, 'a'), (1, 'a')), ((0, 'a'), (1, 'b')), ((0, 'a'), (2, 'b')), ((0, 'b'), (1, 'a')), ((0, 'b'), (1, 'b')), ((0, 'b'), (2, 'a')), ((1, 'a'), (1, 'b')), ((1, 'a'), (2, 'a')), ((1, 'a'), (2, 'b')), ((1, 'b'), (2, 'a')), ((1, 'b'), (2, 'b')), ((2, 'a'), (2, 'b'))]
            sage: T.is_isomorphic( H.disjunctive_product(G) )
            True

        Disjunctive product of digraphs::

            sage: I = DiGraph([(0,1), (1,2)])
            sage: J = DiGraph([('a','b')])
            sage: T = I.disjunctive_product(J)
            sage: T.edges(labels=None)
            [((0, 'a'), (0, 'b')), ((0, 'a'), (1, 'a')), ((0, 'a'), (1, 'b')), ((0, 'a'), (2, 'b')), ((0, 'b'), (1, 'a')), ((0, 'b'), (1, 'b')), ((1, 'a'), (0, 'b')), ((1, 'a'), (1, 'b')), ((1, 'a'), (2, 'a')), ((1, 'a'), (2, 'b')), ((1, 'b'), (2, 'a')), ((1, 'b'), (2, 'b')), ((2, 'a'), (0, 'b')), ((2, 'a'), (1, 'b')), ((2, 'a'), (2, 'b'))]
            sage: T.is_isomorphic( J.disjunctive_product(I) )
            True
        """
        self._scream_if_not_simple(allow_loops=True)
        if self._directed and other._directed:
            from sage.graphs.all import DiGraph
            G = DiGraph( loops = (self.has_loops() or other.has_loops()) )
        elif (not self._directed) and (not other._directed):
            from sage.graphs.all import Graph
            G = Graph()
        else:
            raise TypeError('the graphs should be both directed or both undirected')

        G.add_vertices( [(u,v) for u in self for v in other] )
        for u,w in self.edge_iterator(labels=None):
            for v in other:
                for x in other:
                    G.add_edge((u,v), (w,x))
        for v,x in other.edge_iterator(labels=None):
            for u in self:
                for w in self:
                    G.add_edge((u,v), (w,x))
        return G

    def transitive_closure(self):
        r"""
        Computes the transitive closure of a graph and returns it. The
        original graph is not modified.

        The transitive closure of a graph G has an edge (x,y) if and only
        if there is a path between x and y in G.

        The transitive closure of any strongly connected component of a
        graph is a complete graph. In particular, the transitive closure of
        a connected undirected graph is a complete graph. The transitive
        closure of a directed acyclic graph is a directed acyclic graph
        representing the full partial order.

        EXAMPLES::

            sage: g=graphs.PathGraph(4)
            sage: g.transitive_closure()
            Transitive closure of Path Graph: Graph on 4 vertices
            sage: g.transitive_closure()==graphs.CompleteGraph(4)
            True
            sage: g=DiGraph({0:[1,2], 1:[3], 2:[4,5]})
            sage: g.transitive_closure().edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 3), (2, 4), (2, 5)]

        """
        from copy import copy
        G = copy(self)
        G.name('Transitive closure of ' + self.name())
        for v in G:
            # todo optimization opportunity: we are adding edges that
            # are already in the graph and we are adding edges
            # one at a time.
            for e in G.breadth_first_search(v):
                G.add_edge((v,e))
        return G

    def transitive_reduction(self):
        r"""
        Returns a transitive reduction of a graph. The original graph is
        not modified.

        A transitive reduction H of G has a path from x to y if and only if
        there was a path from x to y in G. Deleting any edge of H destroys
        this property. A transitive reduction is not unique in general. A
        transitive reduction has the same transitive closure as the
        original graph.

        A transitive reduction of a complete graph is a tree. A transitive
        reduction of a tree is itself.

        EXAMPLES::

            sage: g=graphs.PathGraph(4)
            sage: g.transitive_reduction()==g
            True
            sage: g=graphs.CompleteGraph(5)
            sage: edges = g.transitive_reduction().edges(); len(edges)
            4
            sage: g=DiGraph({0:[1,2], 1:[2,3,4,5], 2:[4,5]})
            sage: g.transitive_reduction().size()
            5
        """
        from copy import copy
        from sage.rings.infinity import Infinity
        G = copy(self)
        for e in self.edge_iterator():
            # Try deleting the edge, see if we still have a path
            # between the vertices.
            G.delete_edge(e)
            if G.distance(e[0],e[1])==Infinity:
                # oops, we shouldn't have deleted it
                G.add_edge(e)
        return G

    def is_transitively_reduced(self):
        r"""
        Tests whether the digraph is transitively reduced.

        A digraph is transitively reduced if it is equal to its transitive
        reduction.

        EXAMPLES::

            sage: d = DiGraph({0:[1],1:[2],2:[3]})
            sage: d.is_transitively_reduced()
            True

            sage: d = DiGraph({0:[1,2],1:[2]})
            sage: d.is_transitively_reduced()
            False

            sage: d = DiGraph({0:[1,2],1:[2],2:[]})
            sage: d.is_transitively_reduced()
            False
        """
        from copy import copy
        from sage.rings.infinity import Infinity
        G = copy(self)
        for e in self.edge_iterator():
            G.delete_edge(e)
            if G.distance(e[0],e[1]) == Infinity:
                G.add_edge(e)
            else:
                return False
        return True


    ### Visualization

    def _color_by_label(self, format='hex', as_function=False, default_color = "black"):
        """
        Coloring of the edges according to their label for plotting

        INPUT:

         - ``format`` -- "rgbtuple", "hex", ``True`` (same as "hex"),
           or a function or dictionary assigning colors to labels
           (default: "hex")
         - ``default_color`` -- a color (as a string) or None (default: "black")
         - ``as_function`` -- boolean (default: ``False``)

        OUTPUT: A coloring of the edges of this graph.

        If ``as_function`` is ``True``, then the coloring is returned
        as a function assigning a color to each label. Otherwise (the
        default, for backward compatibility), the coloring is returned
        as a dictionary assigning to each color the list of the edges
        of the graph of that color.

        This is factored out from plot() for use in 3d plots, etc.

        If ``format`` is a function, then it is used directly as
        coloring. Otherwise, for each label a default color is chosen
        along a rainbow (see :func:`sage.plot.colors.rainbow`). If
        ``format`` is a dictionary, then the colors specified there
        override the default choices.

        EXAMPLES:

        We consider the Cayley graph of the symmetric group, whose
        edges are labelled by the numbers 1,2, and 3::

            sage: G = SymmetricGroup(4).cayley_graph()
            sage: set(G.edge_labels())
            set([1, 2, 3])

        We first request the coloring as a function::

            sage: f = G._color_by_label(as_function=True)
            sage: [f(1), f(2), f(3)]
            ['#00ff00', '#ff0000', '#0000ff']
            sage: f = G._color_by_label({1: "blue", 2: "red", 3: "green"}, as_function=True)
            sage: [f(1), f(2), f(3)]
            ['blue', 'red', 'green']
            sage: f = G._color_by_label({1: "red"}, as_function=True)
            sage: [f(1), f(2), f(3)]
            ['red', 'black', 'black']
            sage: f = G._color_by_label({1: "red"}, as_function=True, default_color = 'blue')
            sage: [f(1), f(2), f(3)]
            ['red', 'blue', 'blue']

        The default output is a dictionary assigning edges to colors::

            sage: G._color_by_label()
            {'#00ff00': [((1,4,3,2), (1,4,3), 1), ... ((1,2)(3,4), (3,4),     1)],
             '#ff0000': [((1,4,3,2), (1,4,2), 2), ... ((1,2)(3,4), (1,3,4,2), 2)],
             '#0000ff': [((1,4,3,2), (1,3,2), 3), ... ((1,2)(3,4), (1,2),     3)]}

            sage: G._color_by_label({1: "blue", 2: "red", 3: "green"})
            {'blue':    [((1,4,3,2), (1,4,3), 1), ... ((1,2)(3,4), (3,4),     1)],
             'green':   [((1,4,3,2), (1,3,2), 3), ... ((1,2)(3,4), (1,2),     3)],
             'red':     [((1,4,3,2), (1,4,2), 2), ... ((1,2)(3,4), (1,3,4,2), 2)]}

        TESTS:

        We check what happens when several labels have the same color::

            sage: result = G._color_by_label({1: "blue", 2: "blue", 3: "green"})
            sage: result.keys()
            ['blue', 'green']
            sage: len(result['blue'])
            48
            sage: len(result['green'])
            24
        """
        if format is True:
            format = "hex"
        if isinstance(format, str):
            # Find all labels; this is slower and huglier than:
            #    labels = set(edge[2] for edge in self.edge_iterator())
            # but works with non hashable edge labels, and keeps backward
            # compatibility for the label ordering.
            labels = []
            for edge in self.edge_iterator():
                label = edge[2]
                if label not in labels:
                    labels.append(label)

            from sage.plot.colors import rainbow
            colors = rainbow(len(labels), format=format)
            color_of_label = dict(zip(labels, colors))
            color_of_label = color_of_label.__getitem__
        elif isinstance(format, dict):
            color_of_label = lambda label: format.get(label, default_color)
        else:
            # This assumes that ``format`` is already a function
            color_of_label = format

        if as_function:
            return color_of_label

        edges_by_color = {}
        for edge in self.edge_iterator():
            color = color_of_label(edge[2])
            if color in edges_by_color:
                edges_by_color[color].append(edge)
            else:
                edges_by_color[color] = [edge]
        return edges_by_color

    def latex_options(self):
        r"""
        Returns an instance of
        :class:`~sage.graphs.graph_latex.GraphLatex` for the graph.

        Changes to this object will affect the LaTeX version of the graph.  For
        a full explanation of how to use LaTeX to render graphs, see the
        introduction to the :mod:`~sage.graphs.graph_latex` module.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts
            LaTeX options for Petersen graph: {}
            sage: opts.set_option('tkz_style', 'Classic')
            sage: opts
            LaTeX options for Petersen graph: {'tkz_style': 'Classic'}
        """
        if self._latex_opts is None:
            from sage.graphs.graph_latex import GraphLatex
            self._latex_opts = GraphLatex(self)
        return self._latex_opts

    def set_latex_options(self, **kwds):
        r"""
        Sets multiple options for rendering a graph with LaTeX.

        INPUTS:

        - ``kwds`` - any number of option/value pairs to set many graph
          latex options at once (a variable number, in any
          order). Existing values are overwritten, new values are
          added.  Existing values can be cleared by setting the value
          to ``None``.  Possible options are documented at
          :meth:`sage.graphs.graph_latex.GraphLatex.set_option`.

        This method is a convenience for setting the options of a graph
        directly on an instance of the graph.  For a full explanation of
        how to use LaTeX to render graphs, see the introduction to the
        :mod:`~sage.graphs.graph_latex` module.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.set_latex_options(tkz_style = 'Welsh')
            sage: opts = g.latex_options()
            sage: opts.get_option('tkz_style')
            'Welsh'
        """
        opts = self.latex_options()
        opts.set_options(**kwds)


    def layout(self, layout = None, pos = None, dim = 2, save_pos = False, **options):
        """
        Returns a layout for the vertices of this graph.

        INPUT:

         - layout -- one of "acyclic", "circular", "ranked", "graphviz", "planar", "spring", or "tree"

         - pos -- a dictionary of positions or None (the default)

         - save_pos -- a boolean

         - layout options -- (see below)

        If ``layout=algorithm`` is specified, this algorithm is used
        to compute the positions.

        Otherwise, if ``pos`` is specified, use the given positions.

        Otherwise, try to fetch previously computed and saved positions.

        Otherwise use the default layout (usually the spring layout)

        If ``save_pos = True``, the layout is saved for later use.

        EXAMPLES::

            sage: g = digraphs.ButterflyGraph(1)
            sage: g.layout()
            {('1', 1): [2.50..., -0.545...],
             ('0', 0): [2.22..., 0.832...],
             ('1', 0): [1.12..., -0.830...],
             ('0', 1): [0.833..., 0.543...]}

            sage: 1+1
            2
            sage: x = g.layout(layout = "acyclic_dummy", save_pos = True)
            sage: x =  {('1', 1): [41, 18], ('0', 0): [41, 90], ('1', 0): [140, 90], ('0', 1): [141, 18]}

            {('1', 1): [41, 18], ('0', 0): [41, 90], ('1', 0): [140, 90], ('0', 1): [141, 18]}


            sage: g.layout(dim = 3)
            {('1', 1): [1.07..., -0.260..., 0.927...],
             ('0', 0): [2.02..., 0.528..., 0.343...],
             ('1', 0): [0.674..., -0.528..., -0.343...],
             ('0', 1): [1.61..., 0.260..., -0.927...]}

        Here is the list of all the available layout options::

            sage: from sage.graphs.graph_plot import layout_options
            sage: for key, value in list(sorted(layout_options.iteritems())):
            ...      print "option", key, ":", value
            option by_component : Whether to do the spring layout by connected component -- a boolean.
            option dim : The dimension of the layout -- 2 or 3.
            option heights : A dictionary mapping heights to the list of vertices at this height.
            option iterations : The number of times to execute the spring layout algorithm.
            option layout : A layout algorithm -- one of : "acyclic", "circular" (plots the graph with vertices evenly distributed on a circle), "ranked", "graphviz", "planar", "spring" (traditional spring layout, using the graph's current positions as initial positions), or "tree" (the tree will be plotted in levels, depending on minimum distance for the root).
            option prog : Which graphviz layout program to use -- one of "circo", "dot", "fdp", "neato", or "twopi".
            option save_pos : Whether or not to save the computed position for the graph.
            option spring : Use spring layout to finalize the current layout.
            option tree_orientation : The direction of tree branches -- 'up', 'down', 'left' or 'right'.
            option tree_root : A vertex designation for drawing trees. A vertex of the tree to be used as the root for the ``layout='tree'`` option. If no root is specified, then one is chosen close to the center of the tree. Ignored unless ``layout='tree'``

        Some of them only apply to certain layout algorithms. For
        details, see :meth:`.layout_acyclic`, :meth:`.layout_planar`,
        :meth:`.layout_circular`, :meth:`.layout_spring`, ...

        ..warning: unknown optional arguments are silently ignored

        ..warning: ``graphviz`` and ``dot2tex`` are currently required
        to obtain a nice 'acyclic' layout. See
        :meth:`.layout_graphviz` for installation instructions.

        A subclass may implement another layout algorithm `blah`, by
        implementing a method ``.layout_blah``. It may override
        the default layout by overriding :meth:`.layout_default`, and
        similarly override the predefined layouts.

        TODO: use this feature for all the predefined graphs classes
        (like for the Petersen graph, ...), rather than systematically
        building the layout at construction time.
        """
        if layout is None:
            if pos is None:
                pos = self.get_pos(dim = dim)
            if pos is None:
                layout = 'default'

        if hasattr(self, "layout_%s"%layout):
            pos = getattr(self, "layout_%s"%layout)(dim = dim, **options)
        elif layout is not None:
            raise ValueError("unknown layout algorithm: %s"%layout)

        if len(pos) < self.order():
            pos = self.layout_extend_randomly(pos, dim = dim)

        if save_pos:
            self.set_pos(pos, dim = dim)
        return pos


    def layout_spring(self, by_component = True, **options):
        """
        Computes a spring layout for this graph

        INPUT:

         - ``iterations`` -- a positive integer
         - ``dim`` -- 2 or 3 (default: 2)

        OUTPUT: a dictionary mapping vertices to positions

        Returns a layout computed by randomly arranging the vertices
        along the given heights

        EXAMPLES::

            sage: g = graphs.LadderGraph(3) #TODO!!!!
            sage: g.layout_spring()
            {0: [1.28..., -0.943...],
             1: [1.57..., -0.101...],
             2: [1.83..., 0.747...],
             3: [0.531..., -0.757...],
             4: [0.795..., 0.108...],
             5: [1.08..., 0.946...]}
            sage: g = graphs.LadderGraph(7)
            sage: g.plot(layout = "spring")
        """
        return spring_layout_fast(self, by_component = by_component, **options)

    layout_default = layout_spring

#     if not isinstance(graph.get_pos(), dict):
#         if graph.is_planar():
#             graph.set_planar_positions()
#         else:
#             import sage.graphs.generic_graph_pyx as ggp
#             graph.set_pos(ggp.spring_layout_fast_split(graph, iterations=1000))

    def layout_ranked(self, heights = None, dim = 2, spring = False, **options):
        """
        Computes a ranked layout for this graph

        INPUT:

         - heights -- a dictionary mapping heights to the list of vertices at this height

        OUTPUT: a dictionary mapping vertices to positions

        Returns a layout computed by randomly arranging the vertices
        along the given heights

        EXAMPLES::

            sage: g = graphs.LadderGraph(3)
            sage: g.layout_ranked(heights = dict( (i,[i, i+3]) for i in range(3) ))
            {0: [0.668..., 0],
             1: [0.667..., 1],
             2: [0.677..., 2],
             3: [1.34..., 0],
             4: [1.33..., 1],
             5: [1.33..., 2]}
            sage: g = graphs.LadderGraph(7)
            sage: g.plot(layout = "ranked", heights = dict( (i,[i, i+7]) for i in range(7) ))
        """
        assert heights is not None

        from sage.misc.randstate import current_randstate
        random = current_randstate().python_random().random

        if self.order() == 0:
            return {}

        pos = {}
        mmax = max([len(ccc) for ccc in heights.values()])
        ymin = min(heights.keys())
        ymax = max(heights.keys())
        dist = (max(ymax-ymin, 1)) / (mmax+1.0)
        for height in heights:
            num_xs = len(heights[height])
            if num_xs == 0:
                continue
            j = (mmax - num_xs)/2.0
            for k in range(num_xs):
                pos[heights[height][k]] = [ dist*(j+k+1) + random()*(dist*0.03) for i in range(dim-1) ] + [height]
        if spring:
            # This does not work that well in 2d, since the vertices on
            # the same level are unlikely to cross. It is also hard to
            # set a good equilibrium distance (parameter k in
            # run_spring). If k<1, the layout gets squished
            # horizontally.  If k>1, then two adjacent vertices in
            # consecutive levels tend to be further away than desired.
            newpos = spring_layout_fast(self,
                                        vpos = pos,
                                        dim = dim,
                                        height = True,
                                        **options)
            # spring_layout_fast actually *does* touch the last coordinates
            # (conversion to floats + translation)
            # We restore back the original height.
            for x in self.vertices():
                newpos[x][dim-1] = pos[x][dim-1]
            pos = newpos
        return pos

    def layout_extend_randomly(self, pos, dim = 2):
        """
        Extends randomly a partial layout

        INPUT:

         - ``pos``: a dictionary mapping vertices to positions

        OUTPUT: a dictionary mapping vertices to positions

        The vertices not referenced in ``pos`` are assigned random
        positions within the box delimited by the other vertices.

        EXAMPLES::

            sage: H = digraphs.ButterflyGraph(1)
            sage: H.layout_extend_randomly({('0',0): (0,0), ('1',1): (1,1)})
            {('1', 1): (1, 1),
             ('0', 0): (0, 0),
             ('1', 0): [0.111..., 0.514...],
             ('0', 1): [0.0446..., 0.332...]}
        """
        assert dim == 2 # 3d not yet implemented
        from sage.misc.randstate import current_randstate
        random = current_randstate().python_random().random

        xmin, xmax,ymin, ymax = self._layout_bounding_box(pos)

        dx = xmax - xmin
        dy = ymax - ymin
        # Check each vertex position is in pos, add position
        # randomly within the plot range if none is defined
        for v in self:
            if not v in pos:
                pos[v] = [xmin + dx*random(), ymin + dy*random()]
        return pos


    def layout_circular(self, dim = 2, **options):
        """
        Computes a circular layout for this graph

        OUTPUT: a dictionary mapping vertices to positions

        EXAMPLES::

            sage: G = graphs.CirculantGraph(7,[1,3])
            sage: G.layout_circular()
            {0: [6.12...e-17, 1.0],
             1: [-0.78...,  0.62...],
             2: [-0.97..., -0.22...],
             3: [-0.43..., -0.90...],
             4: [0.43...,  -0.90...],
             5: [0.97...,  -0.22...],
             6: [0.78...,   0.62...]}
            sage: G.plot(layout = "circular")
        """
        assert dim == 2, "3D circular layout not implemented"
        from math import sin, cos, pi
        verts = self.vertices()
        n = len(verts)
        pos = {}
        for i in range(n):
            x = float(cos((pi/2) + ((2*pi)/n)*i))
            y = float(sin((pi/2) + ((2*pi)/n)*i))
            pos[verts[i]] = [x,y]
        return pos

    def layout_tree(self, tree_orientation = "down", tree_root = None, dim = 2, **options):
        r"""
        Computes an ordered tree layout for this graph, which should
        be a tree (no non-oriented cycles).

        INPUT:

        - ``tree_root`` -- the root vertex. By default ``None``. In
          this case, a vertex is chosen close to the center of the
          tree.

        - ``tree_orientation`` -- the direction in which the tree is
          growing, can be 'up', 'down', 'left' or 'right' (default is
          'down')

        OUTPUT: a dictionary mapping vertices to positions

        EXAMPLES::

            sage: T = graphs.RandomLobster(25, 0.3, 0.3)
            sage: T.show(layout='tree', tree_orientation='up')

            sage: G = graphs.HoffmanSingletonGraph()
            sage: T = Graph()
            sage: T.add_edges(G.min_spanning_tree(starting_vertex=0))
            sage: T.show(layout='tree', tree_root=0)

            sage: G = graphs.BalancedTree(2, 2)
            sage: G.layout_tree(tree_root = 0)
            {0: (1.5, 0),
             1: (2.5, -1),
             2: (0.5, -1),
             3: (3.0, -2),
             4: (2.0, -2),
             5: (1.0, -2),
             6: (0.0, -2)}

            sage: G = graphs.BalancedTree(2,4)
            sage: G.plot(layout="tree", tree_root = 0, tree_orientation = "up")

            sage: G = graphs.RandomTree(80)
            sage: G.plot(layout="tree", tree_orientation = "right")

        TESTS::

            sage: G = graphs.CycleGraph(3)
            sage: G.plot(layout='tree')
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot use tree layout on this graph: self.is_tree() returns False.
        """
        if not(dim == 2):
            raise ValueError('only implemented in 2D')

        from sage.graphs.all import Graph
        if not Graph(self).is_tree():
            raise RuntimeError("Cannot use tree layout on this graph: self.is_tree() returns False.")

        n = self.order()
        vertices = self.vertices()

        if tree_root is None:
            root = self.center()[0]
        else:
            root = tree_root

        pos = {}

        # The children and parent of each vertex
        children = {root:self.neighbors(root)}
        parent = {u:root for u in children[root]}

        # stack[i] is the list of children of stick[i] which have not been given
        # a position yet.
        stack = [list(children[root])]
        stick = [root]

        # obstruction[y] is the smallest value of x to which a vertex at height
        # y can be assigned. All vertices at height y which have already been
        # assigned are on the left of (x-1,y).
        obstruction = [0.0]*self.num_verts()

        if tree_orientation in ['down', 'left']:
            o = -1
        elif tree_orientation in ['up', 'right']:
            o = 1
        else:
            raise ValueError('orientation should be "up", "down", "left" or "right"')

        def slide(v, dx):
            """
            shift the vertex v and its descendants to the right by dx

            Precondition: v and its descendents have already had their
            positions computed.
            """
            level = [v]
            while level:
                nextlevel = []
                for u in level:
                    x, y = pos[u]
                    x += dx
                    obstruction[y] = max(x+1, obstruction[y])
                    pos[u] = x, y
                    nextlevel += children[u]

                level = nextlevel

        while stack:
            C = stack[-1]

            # If all the children of stick[-1] have been given a position
            if len(C) == 0:
                p = stick.pop()
                stack.pop()
                cp = children[p]
                y = o*len(stack)

                if len(cp) == 0:
                    # If p has no children, we draw it at the leftmost position
                    # which has not been forbidden
                    x = obstruction[y]
                    pos[p] = x, y
                else:
                    # If p has children, we put v on a vertical line going
                    # through the barycenter of its children
                    x = sum([pos[c][0] for c in cp])/len(cp)
                    pos[p] = x, y
                    ox = obstruction[y]
                    if x < ox:
                        slide(p, ox-x)
                        x = ox

                # If the vertex to the right of p has not children, we want it
                # at distance 1 from p
                obstruction[y] = x+1

            # Otherwise, we take one of the children and add it to the
            # stack. Note that this vertex is removed from the list C.
            else:
                t = C.pop()

                pt = parent[t]
                ct = [u for u in self.neighbors(t) if u != pt]
                children[t] = ct

                for c in ct:
                    parent[c] = t

                stack.append([c for c in ct])
                stick.append(t)

        if tree_orientation in ['right', 'left']:
            return {p:(py,px) for p,(px,py) in pos.iteritems()}

        return pos

    def layout_graphviz(self, dim = 2, prog = 'dot', **options):
        """
        Calls ``graphviz`` to compute a layout of the vertices of this graph.

        INPUT:

         - ``prog`` -- one of "dot", "neato", "twopi", "circo", or "fdp"

        EXAMPLES::

            sage: g = digraphs.ButterflyGraph(2)
            sage: g.layout_graphviz() # optional - dot2tex, graphviz
            {('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...],
             ('...', ...): [...,...]}
            sage: g.plot(layout = "graphviz") # optional - dot2tex, graphviz

        Note: the actual coordinates are not deterministic

        By default, an acyclic layout is computed using ``graphviz``'s
        ``dot`` layout program. One may specify an alternative layout
        program::

            sage: g.plot(layout = "graphviz", prog = "dot")   # optional - dot2tex, graphviz
            sage: g.plot(layout = "graphviz", prog = "neato") # optional - dot2tex, graphviz
            sage: g.plot(layout = "graphviz", prog = "twopi") # optional - dot2tex, graphviz
            sage: g.plot(layout = "graphviz", prog = "fdp")   # optional - dot2tex, graphviz
            sage: g = graphs.BalancedTree(5,2)
            sage: g.plot(layout = "graphviz", prog = "circo") # optional - dot2tex, graphviz

        TODO: put here some cool examples showcasing graphviz features.

        This requires ``graphviz`` and the ``dot2tex`` spkg. Here are
        some installation tips:

         - Install graphviz >= 2.14 so that the programs dot, neato, ...
           are in your path. The graphviz suite can be download from
           http://graphviz.org.

         - Download dot2tex-2.8.?.spkg from http://trac.sagemath.org/sage_trac/ticket/7004
           and install it with ``sage -i dot2tex-2.8.?.spkg``

        TODO: use the graphviz functionality of Networkx 1.0 once it
        will be merged into Sage.
        """
        assert_have_dot2tex()
        assert dim == 2, "3D graphviz layout not implemented"

        key = self._keys_for_vertices()
        key_to_vertex = dict( (key(v), v) for v in self )

        import dot2tex
        positions = dot2tex.dot2tex(self.graphviz_string(**options), format = "positions", prog = prog)

        return dict( (key_to_vertex[key], pos) for (key, pos) in positions.iteritems() )

    def _layout_bounding_box(self, pos):
        """
        INPUT:

         - pos -- a dictionary of positions

        Returns a bounding box around the specified positions

        EXAMPLES::

            sage: Graph()._layout_bounding_box( {} )
            [-1, 1, -1, 1]
            sage: Graph()._layout_bounding_box( {0: [3,5], 1: [2,7], 2: [-4,2] } )
            [-4, 3, 2, 7]
            sage: Graph()._layout_bounding_box( {0: [3,5], 1: [3.00000000001,4.999999999999999] } )
            [2, 4.00000000001000, 4.00000000000000, 6]
        """
        xs = [pos[v][0] for v in pos]
        ys = [pos[v][1] for v in pos]
        if len(xs) == 0:
            xmin = -1
            xmax =  1
            ymin = -1
            ymax =  1
        else:
            xmin = min(xs)
            xmax = max(xs)
            ymin = min(ys)
            ymax = max(ys)

        if xmax - xmin < 0.00000001:
            xmax += 1
            xmin -= 1

        if ymax - ymin < 0.00000001:
            ymax += 1
            ymin -= 1

        return [xmin, xmax, ymin, ymax]

    def graphplot(self, **options):
        """
        Returns a GraphPlot object.

        EXAMPLES:

        Creating a graphplot object uses the same options as graph.plot()::

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: g.set_boundary([0,1])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: GP = g.graphplot(edge_labels=True, color_by_label=True, edge_style='dashed')
            sage: GP.plot()

        We can modify the graphplot object.  Notice that the changes are cumulative::

            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            sage: GP.set_vertices(talk=True)
            sage: GP.plot()
        """
        from sage.graphs.graph_plot import GraphPlot
        return GraphPlot(graph=self, options=options)

    @options()
    def plot(self, **options):
        r"""
        Returns a graphics object representing the (di)graph.

        INPUT:

        - ``pos`` - an optional positioning dictionary

        - ``layout`` - what kind of layout to use, takes precedence
          over pos

           - 'circular' -- plots the graph with vertices evenly
             distributed on a circle

           - 'spring' - uses the traditional spring layout, using the
             graph's current positions as initial positions

           - 'tree' - the (di)graph must be a tree. One can specify
             the root of the tree using the keyword tree_root,
             otherwise a root will be selected at random. Then the
             tree will be plotted in levels, depending on minimum
             distance for the root.

        - ``vertex_labels`` - whether to print vertex labels

        - ``edge_labels`` - whether to print edge labels. By default,
          False, but if True, the result of str(l) is printed on the
          edge for each label l. Labels equal to None are not printed
          (to set edge labels, see set_edge_label).

        - ``vertex_size`` - size of vertices displayed

        - ``vertex_shape`` - the shape to draw the vertices (Not
          available for multiedge digraphs.)

        - ``graph_border`` - whether to include a box around the graph

        - ``vertex_colors`` - optional dictionary to specify vertex
          colors: each key is a color recognizable by matplotlib, and
          each corresponding entry is a list of vertices. If a vertex
          is not listed, it looks invisible on the resulting plot (it
          doesn't get drawn).

        - ``edge_colors`` - a dictionary specifying edge colors: each
          key is a color recognized by matplotlib, and each entry is a
          list of edges.

        - ``partition`` - a partition of the vertex set. if specified,
          plot will show each cell in a different color. vertex_colors
          takes precedence.

        - ``talk`` - if true, prints large vertices with white
          backgrounds so that labels are legible on slides

        - ``iterations`` - how many iterations of the spring layout
          algorithm to go through, if applicable

        - ``color_by_label`` - a boolean or dictionary or function (default: False)
           whether to color each edge with a different color according
           to its label; the colors are chosen along a rainbow, unless
           they are specified by a function or dictionary mapping
           labels to colors; this option is incompatible with
           ``edge_color`` and ``edge_colors``.

        - ``heights`` - if specified, this is a dictionary from a set
          of floating point heights to a set of vertices

        - ``edge_style`` - keyword arguments passed into the
          edge-drawing routine.  This currently only works for
          directed graphs, since we pass off the undirected graph to
          networkx

        - ``tree_root`` - a vertex of the tree to be used as the root
          for the layout="tree" option. If no root is specified, then one
          is chosen at random. Ignored unless layout='tree'.

        - ``tree_orientation`` - "up" or "down" (default is "down").
          If "up" (resp., "down"), then the root of the tree will
          appear on the bottom (resp., top) and the tree will grow
          upwards (resp. downwards). Ignored unless layout='tree'.

        - ``save_pos`` - save position computed during plotting

        .. NOTE::

            - See the documentation of the :mod:`sage.graphs.graph_plot` module
              for information and examples of how to define parameters that will
              be applied to **all** graph plots.

            - Default parameters for this method *and a specific graph* can also
              be set through the :class:`~sage.misc.decorators.options`
              mechanism. For more information on this different way to set
              default parameters, see the help of the :class:`options decorator
              <sage.misc.decorators.options>`.

            - See also the :mod:`sage.graphs.graph_latex` module for ways to use
              LaTeX to produce an image of a graph.

        EXAMPLES::

            sage: from sage.graphs.graph_plot import graphplot_options
            sage: list(sorted(graphplot_options.iteritems()))
            [...]

            sage: from math import sin, cos, pi
            sage: P = graphs.PetersenGraph()
            sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
            sage: pos_dict = {}
            sage: for i in range(5):
            ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: for i in range(10)[5:]:
            ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: pl = P.plot(pos=pos_dict, vertex_colors=d)
            sage: pl.show()

        ::

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()

        ::

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).show()

        ::

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} , sparse=True)
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.plot(edge_labels=True, layout='circular').show()

        ::

            sage: from sage.plot.colors import rainbow
            sage: C = graphs.CubeGraph(5)
            sage: R = rainbow(5)
            sage: edge_colors = {}
            sage: for i in range(5):
            ...    edge_colors[R[i]] = []
            sage: for u,v,l in C.edges():
            ...    for i in range(5):
            ...        if u[i] != v[i]:
            ...            edge_colors[R[i]].append((u,v,l))
            sage: C.plot(vertex_labels=False, vertex_size=0, edge_colors=edge_colors).show()

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: Pi = [[6,5,15,14,7],[16,13,8,2,4],[12,17,9,3,1],[0,19,18,10,11]]
            sage: D.show(partition=Pi)

        ::

            sage: G = graphs.PetersenGraph()
            sage: G.allow_loops(True)
            sage: G.add_edge(0,0)
            sage: G.show()

        ::

            sage: D = DiGraph({0:[0,1], 1:[2], 2:[3]}, loops=True)
            sage: D.show()
            sage: D.show(edge_colors={(0,1,0):[(0,1,None),(1,2,None)],(0,0,0):[(2,3,None)]})

        ::

            sage: pos = {0:[0.0, 1.5], 1:[-0.8, 0.3], 2:[-0.6, -0.8], 3:[0.6, -0.8], 4:[0.8, 0.3]}
            sage: g = Graph({0:[1], 1:[2], 2:[3], 3:[4], 4:[0]})
            sage: g.plot(pos=pos, layout='spring', iterations=0)

        ::

            sage: G = Graph()
            sage: P = G.plot()
            sage: P.axes()
            False
            sage: G = DiGraph()
            sage: P = G.plot()
            sage: P.axes()
            False

        ::

            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: (6.12..., 1.0...),
             1: (-0.95..., 0.30...),
             2: (-0.58..., -0.80...),
             3: (0.58..., -0.80...),
             4: (0.95..., 0.30...),
             5: (1.53..., 0.5...),
             6: (-0.47..., 0.15...),
             7: (-0.29..., -0.40...),
             8: (0.29..., -0.40...),
             9: (0.47..., 0.15...)}
            sage: P = G.plot(save_pos=True, layout='spring')

            The following illustrates the format of a position dictionary.

            sage: G.get_pos() # currently random across platforms, see #9593
            {0: [1.17..., -0.855...],
             1: [1.81..., -0.0990...],
             2: [1.35..., 0.184...],
             3: [1.51..., 0.644...],
             4: [2.00..., -0.507...],
             5: [0.597..., -0.236...],
             6: [2.04..., 0.687...],
             7: [1.46..., -0.473...],
             8: [0.902..., 0.773...],
             9: [2.48..., -0.119...]}

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})
            sage: t.set_edge_label(0,1,-7)
            sage: t.set_edge_label(0,5,3)
            sage: t.set_edge_label(0,5,99)
            sage: t.set_edge_label(1,2,1000)
            sage: t.set_edge_label(3,2,'spam')
            sage: t.set_edge_label(2,6,3/2)
            sage: t.set_edge_label(0,4,66)
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]}, edge_labels=True)

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(layout='tree')

        ::

            sage: t = DiGraph('JCC???@A??GO??CO??GO??')
            sage: t.plot(layout='tree', tree_root=0, tree_orientation="up")
            sage: D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sage: D.plot()

            sage: D = DiGraph(multiedges=True,sparse=True)
            sage: for i in range(5):
            ...     D.add_edge((i,i+1,'a'))
            ...     D.add_edge((i,i-1,'b'))
            sage: D.plot(edge_labels=True,edge_colors=D._color_by_label())
            sage: D.plot(edge_labels=True, color_by_label={'a':'blue', 'b':'red'}, edge_style='dashed')

            sage: g = Graph({}, loops=True, multiedges=True,sparse=True)
            sage: g.add_edges([(0,0,'a'),(0,0,'b'),(0,1,'c'),(0,1,'d'),
            ...     (0,1,'e'),(0,1,'f'),(0,1,'f'),(2,1,'g'),(2,2,'h')])
            sage: g.plot(edge_labels=True, color_by_label=True, edge_style='dashed')

        ::

            sage: S = SupersingularModule(389)
            sage: H = S.hecke_matrix(2)
            sage: D = DiGraph(H,sparse=True)
            sage: P = D.plot()

        ::

            sage: G=Graph({'a':['a','b','b','b','e'],'b':['c','d','e'],'c':['c','d','d','d'],'d':['e']},sparse=True)
            sage: G.show(pos={'a':[0,1],'b':[1,1],'c':[2,0],'d':[1,0],'e':[0,0]})

        TESTS::

            sage: G = DiGraph({0:{1:'a', 2:'a'}, 1:{0:'b'}, 2:{0:'c'}})
            sage: p = G.plot(edge_labels=True, color_by_label={'a':'yellow', 'b':'purple'}); p
            sage: sorted([x.options()['rgbcolor'] for x in p if isinstance(x, sage.plot.arrow.CurveArrow)])
            ['black', 'purple', 'yellow', 'yellow']
        """
        return self.graphplot(**options).plot()

    def show(self, method = "matplotlib", **kwds):
        """
        Shows the (di)graph.

        INPUT:

        - ``method`` --

            - If ``method="matplotlib"`` (default) then graph is drawn as a
              picture file, then displayed. In this situation, the method
              accepts any other option understood by
              :meth:`~sage.graphs.generic_graph.plot` (graph-specific) or by
              :meth:`sage.plot.graphics.Graphics.show`.

            - If ``method="js"`` the graph is displayed using the `d3.js
              <http://d3js.org/>`_ library in a browser. In this situation, the
              method accepts any other option understood by
              :meth:`sage.graphs.graph_plot_js.gen_html_code`. Depending on
              whether d3js optional package is installed or not, the javascript
              code will be used locally or fetched from d3js.org website by the
              browser. 

        This method accepts any other option understood by
        :meth:`~sage.graphs.generic_graph.GenericGraph.plot` (graph-specific) or
        by :meth:`sage.plot.graphics.Graphics.show`.

        .. NOTE::

            - See the documentation of the :mod:`sage.graphs.graph_plot` module
              for information on default arguments of this method.

            - For the javascript counterpart, refer to
              :mod:`sage.graphs.graph_plot_js`.

        EXAMPLES::

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()  # long time (3s on sage.math, 2011)
        """
        if method == "js":
            from sage.graphs.graph_plot_js import gen_html_code
            from sage.doctest import DOCTEST_MODE
            filename = gen_html_code(self, **kwds)

            if DOCTEST_MODE:
                return
            from sage.misc.viewer import browser
            import os
            os.system('%s %s 2>/dev/null 1>/dev/null &'% (browser(), filename))
            return

        from graph_plot import graphplot_options

        # This dictionary only contains the options that graphplot
        # understands. These options are removed from kwds at the same
        # time.
        plot_kwds = {k:kwds.pop(k) for k in graphplot_options if k in kwds}

        return self.graphplot(**plot_kwds).show(**kwds)

    def plot3d(self, bgcolor=(1,1,1),
                     vertex_colors=None, vertex_size=0.06, vertex_labels=False,
                     edge_colors=None, edge_size=0.02, edge_size2=0.0325,
                     pos3d=None, color_by_label=False,
                     engine='jmol', **kwds):
        r"""
        Plot a graph in three dimensions.

        See also the :mod:`sage.graphs.graph_latex` module for ways to use LaTeX
        to produce an image of a graph.

        INPUT:

        -  ``bgcolor`` - rgb tuple (default: (1,1,1))

        -  ``vertex_size`` - float (default: 0.06)

        -  ``vertex_labels`` -- a boolean (default: False)
           whether to display vertices using text labels instead of spheres

        -  ``vertex_colors`` - optional dictionary to specify
           vertex colors: each key is a color recognizable by tachyon (rgb
           tuple (default: (1,0,0))), and each corresponding entry is a list
           of vertices. If a vertex is not listed, it looks invisible on the
           resulting plot (it doesn't get drawn).

        -  ``edge_colors`` - a dictionary specifying edge
           colors: each key is a color recognized by tachyon ( default:
           (0,0,0) ), and each entry is a list of edges.

        - ``color_by_label`` - a boolean or dictionary or function (default: False)
           whether to color each edge with a different color according
           to its label; the colors are chosen along a rainbow, unless
           they are specified by a function or dictionary mapping
           labels to colors; this option is incompatible with
           ``edge_color`` and ``edge_colors``.

        -  ``edge_size`` - float (default: 0.02)

        -  ``edge_size2`` - float (default: 0.0325), used for
           Tachyon sleeves

        -  ``pos3d`` - a position dictionary for the vertices

        -  ``layout``, ``iterations``, ... - layout options; see :meth:`.layout`

        -  ``engine`` - which renderer to use. Options:

           -  ``'jmol'`` - default

           -  ``'tachyon'``

        -  ``xres`` - resolution

        -  ``yres`` - resolution

        -  ``**kwds`` - passed on to the rendering engine


        EXAMPLES::

            sage: G = graphs.CubeGraph(5)
            sage: G.plot3d(iterations=500, edge_size=None, vertex_size=0.04) # long time

        We plot a fairly complicated Cayley graph::

            sage: A5 = AlternatingGroup(5); A5
            Alternating group of order 5!/2 as a permutation group
            sage: G = A5.cayley_graph()
            sage: G.plot3d(vertex_size=0.03, edge_size=0.01, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, iterations=200) # long time

        Some Tachyon examples::

            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d(engine='tachyon')
            sage: P3D.show() # long time

        ::

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(engine='tachyon', vertex_colors={(0,0,1):G.vertices()}).show() # long time

        ::

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(engine='tachyon', edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)).show() # long time

        ::

            sage: K = graphs.CompleteGraph(3)
            sage: K.plot3d(engine='tachyon', edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}).show() # long time

        A directed version of the dodecahedron

        ::

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: D.plot3d().show() # long time

        ::

            sage: P = graphs.PetersenGraph().to_directed()
            sage: from sage.plot.colors import rainbow
            sage: edges = P.edges()
            sage: R = rainbow(len(edges), 'rgbtuple')
            sage: edge_colors = {}
            sage: for i in range(len(edges)):
            ...       edge_colors[R[i]] = [edges[i]]
            sage: P.plot3d(engine='tachyon', edge_colors=edge_colors).show() # long time


        ::

            sage: G=Graph({'a':['a','b','b','b','e'],'b':['c','d','e'],'c':['c','d','d','d'],'d':['e']},sparse=True)
            sage: G.show3d()
            Traceback (most recent call last):
            ...
            NotImplementedError: 3D plotting of multiple edges or loops not implemented.

        TESTS::

            sage: G = DiGraph({0:{1:'a', 2:'a'}, 1:{0:'b'}, 2:{0:'c'}})
            sage: p = G.plot3d(edge_labels=True, color_by_label={'a':'yellow', 'b':'cyan'})
            sage: s = p.x3d_str()

        This 3D plot contains four yellow objects (two cylinders and
        two cones), two black objects and 2 cyan objects::

            sage: s.count("Material diffuseColor='1.0 1.0 0.0'")
            4
            sage: s.count("Material diffuseColor='0.0 0.0 0.0'")
            2
            sage: s.count("Material diffuseColor='0.0 1.0 1.0'")
            2

        .. SEEALSO::

            - :meth:`plot`
            - :meth:`graphviz_string`
        """
        import graph_plot
        layout_options = dict( (key,kwds[key]) for key in kwds.keys() if key     in graph_plot.layout_options )
        kwds           = dict( (key,kwds[key]) for key in kwds.keys() if key not in graph_plot.layout_options )
        if pos3d is None:
            pos3d = self.layout(dim=3, **layout_options)

        if self.has_multiple_edges() or self.has_loops():
            raise NotImplementedError("3D plotting of multiple edges or loops not implemented.")
        if engine == 'jmol':
            from sage.plot.plot3d.all import sphere, line3d, arrow3d, text3d
            from sage.plot.plot3d.texture import Texture
            kwds.setdefault('aspect_ratio', [1,1,1])
            verts = self.vertices()

            if vertex_colors is None:
                vertex_colors = { (1,0,0) : verts }

            if color_by_label:
                if edge_colors is  None:
                        # do the coloring
                        edge_colors = self._color_by_label(format=color_by_label)
            elif edge_colors is None:
                edge_colors = { (0,0,0) : self.edges() }

            # by default turn off the frame
            if 'frame' not in kwds:
                kwds['frame'] = False
            # by default make the background given by bgcolor
            if 'background' not in kwds:
                kwds['background'] = bgcolor
            try:
                graphic = 0
                for color in vertex_colors:
                    texture = Texture(color=color, ambient=0.1, diffuse=0.9, specular=0.03)
                    for v in vertex_colors[color]:
                        if vertex_labels:
                            graphic += text3d(repr(v), pos3d[v])
                        else:
                            graphic += sphere(center=pos3d[v], size=vertex_size, texture=texture, **kwds)
                if self._directed:
                    for color in edge_colors:
                        for u, v, l in edge_colors[color]:
                            graphic += arrow3d(pos3d[u], pos3d[v], radius=edge_size, color=color, closed=False, **kwds)

                else:
                    for color in edge_colors:
                        texture = Texture(color=color, ambient=0.1, diffuse=0.9, specular=0.03)
                        for u, v, l in edge_colors[color]:
                            graphic += line3d([pos3d[u], pos3d[v]], radius=edge_size, texture=texture, closed=False, **kwds)

                return graphic

            except KeyError:
                raise KeyError("Oops! You haven't specified positions for all the vertices.")

        elif engine == 'tachyon':
            TT, pos3d = tachyon_vertex_plot(self, bgcolor=bgcolor, vertex_colors=vertex_colors,
                                            vertex_size=vertex_size, pos3d=pos3d, **kwds)
            edges = self.edges()

            if color_by_label:
                if edge_colors is  None:
                    # do the coloring
                    edge_colors = self._color_by_label(format=color_by_label)

            if edge_colors is None:
                edge_colors = { (0,0,0) : edges }

            i = 0

            for color in edge_colors:
                i += 1
                TT.texture('edge_color_%d'%i, ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=color)
                if self._directed:
                    for u,v,l in edge_colors[color]:
                        TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]),
                                      (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)
                        TT.fcylinder( (0.25*pos3d[u][0] + 0.75*pos3d[v][0],
                                       0.25*pos3d[u][1] + 0.75*pos3d[v][1],
                                       0.25*pos3d[u][2] + 0.75*pos3d[v][2],),
                                      (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size2,'edge_color_%d'%i)
                else:
                    for u, v, l in edge_colors[color]:
                        TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]), (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)

            return TT

        else:
            raise TypeError("Rendering engine (%s) not implemented."%engine)

    def show3d(self, bgcolor=(1,1,1), vertex_colors=None, vertex_size=0.06,
                     edge_colors=None, edge_size=0.02, edge_size2=0.0325,
                     pos3d=None, color_by_label=False,
                     engine='jmol', **kwds):
        """
        Plots the graph using Tachyon, and shows the resulting plot.

        INPUT:


        -  ``bgcolor`` - rgb tuple (default: (1,1,1))

        -  ``vertex_size`` - float (default: 0.06)

        -  ``vertex_colors`` - optional dictionary to specify
           vertex colors: each key is a color recognizable by tachyon (rgb
           tuple (default: (1,0,0))), and each corresponding entry is a list
           of vertices. If a vertex is not listed, it looks invisible on the
           resulting plot (it doesn't get drawn).

        -  ``edge_colors`` - a dictionary specifying edge
           colors: each key is a color recognized by tachyon ( default:
           (0,0,0) ), and each entry is a list of edges.

        -  ``edge_size`` - float (default: 0.02)

        -  ``edge_size2`` - float (default: 0.0325), used for
           Tachyon sleeves

        -  ``pos3d`` - a position dictionary for the vertices

        -  ``iterations`` - how many iterations of the spring
           layout algorithm to go through, if applicable

        -  ``engine`` - which renderer to use. Options:

        -  ``'jmol'`` - default 'tachyon'

        -  ``xres`` - resolution

        -  ``yres`` - resolution

        -  ``**kwds`` - passed on to the Tachyon command


        EXAMPLES::

            sage: G = graphs.CubeGraph(5)
            sage: G.show3d(iterations=500, edge_size=None, vertex_size=0.04) # long time

        We plot a fairly complicated Cayley graph::

            sage: A5 = AlternatingGroup(5); A5
            Alternating group of order 5!/2 as a permutation group
            sage: G = A5.cayley_graph()
            sage: G.show3d(vertex_size=0.03, edge_size=0.01, edge_size2=0.02, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, iterations=200) # long time

        Some Tachyon examples::

            sage: D = graphs.DodecahedralGraph()
            sage: D.show3d(engine='tachyon') # long time

        ::

            sage: G = graphs.PetersenGraph()
            sage: G.show3d(engine='tachyon', vertex_colors={(0,0,1):G.vertices()}) # long time

        ::

            sage: C = graphs.CubeGraph(4)
            sage: C.show3d(engine='tachyon', edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)) # long time

        ::

            sage: K = graphs.CompleteGraph(3)
            sage: K.show3d(engine='tachyon', edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}) # long time
        """
        self.plot3d(bgcolor=bgcolor, vertex_colors=vertex_colors,
                    edge_colors=edge_colors, vertex_size=vertex_size, engine=engine,
                    edge_size=edge_size, edge_size2=edge_size2, pos3d=pos3d,
                    color_by_label=color_by_label, **kwds).show()

    @cached_method
    def _keys_for_vertices(self):
        """
        Returns a function mapping each vertex to a unique identifier.

        The identifier is stable iff all vertex labels are unique. It
        is a string not starting with a number, as required by
        dot2tex.

        EXAMPLE::

            sage: g = graphs.Grid2dGraph(5,5)
            sage: g._keys_for_vertices()
            <function get_label at ...>
        """
        import hashlib
        label = dict()
        for i, v in enumerate(self.vertices()):
            label[v] = 'node_{0}'.format(i)
        def get_label(vertex):
            return label[vertex]
        return get_label

    ### String representation to be used by other programs
    @options(labels="string",
            vertex_labels=True,edge_labels=False,
            edge_color=None,edge_colors=None,
            edge_options = (),
            color_by_label=False)
    def graphviz_string(self, **options):
        r"""
        Returns a representation in the dot language.

        The dot language is a text based format for graphs. It is used
        by the software suite graphviz. The specifications of the
        language are available on the web (see the reference [dotspec]_).

        INPUT:

        - ``labels`` - "string" or "latex" (default: "string"). If labels is
          string latex command are not interpreted. This option stands for both
          vertex labels and edge labels.

        - ``vertex_labels`` - boolean (default: True) whether to add the labels
          on vertices.

        - ``edge_labels`` - boolean (default: False) whether to add
          the labels on edges.

        - ``edge_color`` - (default: None) specify a default color for the
          edges.

        - ``edge_colors`` - (default: None) a dictionary whose keys
          are colors and values are list of edges. The list of edges need not to
          be complete in which case the default color is used.

        - ``color_by_label`` - a boolean or dictionary or function (default: False)
           whether to color each edge with a different color according
           to its label; the colors are chosen along a rainbow, unless
           they are specified by a function or dictionary mapping
           labels to colors; this option is incompatible with
           ``edge_color`` and ``edge_colors``.

        - ``edge_options`` - a function (or tuple thereof) mapping
          edges to a dictionary of options for this edge.

        EXAMPLES::

            sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}},sparse=True)
            sage: print G.graphviz_string(edge_labels=True)
            graph {
              node_0  [label="0"];
              node_1  [label="1"];
              node_2  [label="2"];
              node_3  [label="3"];
            <BLANKLINE>
              node_0 -- node_1;
              node_0 -- node_2;
              node_1 -- node_2;
              node_2 -- node_3 [label="foo"];
            }

        A variant, with the labels in latex, for post-processing with ``dot2tex``::

            sage: print G.graphviz_string(edge_labels=True,labels = "latex")
            graph {
              node [shape="plaintext"];
              node_0  [label=" ", texlbl="$0$"];
              node_1  [label=" ", texlbl="$1$"];
              node_2  [label=" ", texlbl="$2$"];
              node_3  [label=" ", texlbl="$3$"];
            <BLANKLINE>
              node_0 -- node_1;
              node_0 -- node_2;
              node_1 -- node_2;
              node_2 -- node_3 [label=" ", texlbl="$\text{\texttt{foo}}$"];
            }

        Same, with a digraph and a color for edges::

            sage: G = DiGraph({0:{1:None,2:None}, 1:{2:None}, 2:{3:'foo'}, 3:{}} ,sparse=True)
            sage: print G.graphviz_string(edge_color="red")
            digraph {
              node_0  [label="0"];
              node_1  [label="1"];
              node_2  [label="2"];
              node_3  [label="3"];
            <BLANKLINE>
            edge [color="red"];
              node_0 -> node_1;
              node_0 -> node_2;
              node_1 -> node_2;
              node_2 -> node_3;
            }

        A digraph using latex labels for vertices and edges::

            sage: f(x) = -1/x
            sage: g(x) = 1/(x+1)
            sage: G = DiGraph()
            sage: G.add_edges([(i,f(i),f) for i in (1,2,1/2,1/4)])
            sage: G.add_edges([(i,g(i),g) for i in (1,2,1/2,1/4)])
            sage: print G.graphviz_string(labels="latex",edge_labels=True)
            digraph {
              node [shape="plaintext"];
              node_7  [label=" ", texlbl="$\frac{2}{3}$"];
              node_5  [label=" ", texlbl="$\frac{1}{3}$"];
              node_6  [label=" ", texlbl="$\frac{1}{2}$"];
              node_9  [label=" ", texlbl="$1$"];
              node_4  [label=" ", texlbl="$\frac{1}{4}$"];
              node_8  [label=" ", texlbl="$\frac{4}{5}$"];
              node_0  [label=" ", texlbl="$-4$"];
              node_10  [label=" ", texlbl="$2$"];
              node_1  [label=" ", texlbl="$-2$"];
              node_3  [label=" ", texlbl="$-\frac{1}{2}$"];
              node_2  [label=" ", texlbl="$-1$"];
            <BLANKLINE>
              node_6 -> node_1 [label=" ", texlbl="$x \ {\mapsto}\ -\frac{1}{x}$"];
              node_6 -> node_7 [label=" ", texlbl="$x \ {\mapsto}\ \frac{1}{x + 1}$"];
              node_9 -> node_2 [label=" ", texlbl="$x \ {\mapsto}\ -\frac{1}{x}$"];
              node_9 -> node_6 [label=" ", texlbl="$x \ {\mapsto}\ \frac{1}{x + 1}$"];
              node_4 -> node_0 [label=" ", texlbl="$x \ {\mapsto}\ -\frac{1}{x}$"];
              node_4 -> node_8 [label=" ", texlbl="$x \ {\mapsto}\ \frac{1}{x + 1}$"];
              node_10 -> node_3 [label=" ", texlbl="$x \ {\mapsto}\ -\frac{1}{x}$"];
              node_10 -> node_5 [label=" ", texlbl="$x \ {\mapsto}\ \frac{1}{x + 1}$"];
            }

            sage: print G.graphviz_string(labels="latex",color_by_label=True)
            digraph {
              node [shape="plaintext"];
              node_7  [label=" ", texlbl="$\frac{2}{3}$"];
              node_5  [label=" ", texlbl="$\frac{1}{3}$"];
              node_6  [label=" ", texlbl="$\frac{1}{2}$"];
              node_9  [label=" ", texlbl="$1$"];
              node_4  [label=" ", texlbl="$\frac{1}{4}$"];
              node_8  [label=" ", texlbl="$\frac{4}{5}$"];
              node_0  [label=" ", texlbl="$-4$"];
              node_10  [label=" ", texlbl="$2$"];
              node_1  [label=" ", texlbl="$-2$"];
              node_3  [label=" ", texlbl="$-\frac{1}{2}$"];
              node_2  [label=" ", texlbl="$-1$"];
            <BLANKLINE>
              node_6 -> node_1 [color = "#ff0000"];
              node_6 -> node_7 [color = "#00ffff"];
              node_9 -> node_2 [color = "#ff0000"];
              node_9 -> node_6 [color = "#00ffff"];
              node_4 -> node_0 [color = "#ff0000"];
              node_4 -> node_8 [color = "#00ffff"];
              node_10 -> node_3 [color = "#ff0000"];
              node_10 -> node_5 [color = "#00ffff"];
            }

            sage: print G.graphviz_string(labels="latex",color_by_label={ f: "red", g: "blue" })
            digraph {
              node [shape="plaintext"];
              node_7  [label=" ", texlbl="$\frac{2}{3}$"];
              node_5  [label=" ", texlbl="$\frac{1}{3}$"];
              node_6  [label=" ", texlbl="$\frac{1}{2}$"];
              node_9  [label=" ", texlbl="$1$"];
              node_4  [label=" ", texlbl="$\frac{1}{4}$"];
              node_8  [label=" ", texlbl="$\frac{4}{5}$"];
              node_0  [label=" ", texlbl="$-4$"];
              node_10  [label=" ", texlbl="$2$"];
              node_1  [label=" ", texlbl="$-2$"];
              node_3  [label=" ", texlbl="$-\frac{1}{2}$"];
              node_2  [label=" ", texlbl="$-1$"];
            <BLANKLINE>
              node_6 -> node_1 [color = "red"];
              node_6 -> node_7 [color = "blue"];
              node_9 -> node_2 [color = "red"];
              node_9 -> node_6 [color = "blue"];
              node_4 -> node_0 [color = "red"];
              node_4 -> node_8 [color = "blue"];
              node_10 -> node_3 [color = "red"];
              node_10 -> node_5 [color = "blue"];
            }

        Edge-specific options can also be specified by providing a
        function (or tuple thereof) which maps each edge to a
        dictionary of options. Valid options are "color", "backward"
        (a boolean), "dot" (a string containing a sequence of options
        in dot format), "label" (a string), "label_style" ("string" or
        "latex"), "edge_string" ("--" or "->"). Here we state that the
        graph should be laid out so that edges starting from ``1`` are
        going backward (e.g. going up instead of down)::

            sage: def edge_options((u,v,label)):
            ...       return { "backward": u == 1 }
            sage: print G.graphviz_string(edge_options = edge_options)
            digraph {
              node_7  [label="2/3"];
              node_5  [label="1/3"];
              node_6  [label="1/2"];
              node_9  [label="1"];
              node_4  [label="1/4"];
              node_8  [label="4/5"];
              node_0  [label="-4"];
              node_10  [label="2"];
              node_1  [label="-2"];
              node_3  [label="-1/2"];
              node_2  [label="-1"];
            <BLANKLINE>
              node_6 -> node_1;
              node_6 -> node_7;
              node_2 -> node_9 [dir=back];
              node_6 -> node_9 [dir=back];
              node_4 -> node_0;
              node_4 -> node_8;
              node_10 -> node_3;
              node_10 -> node_5;
            }

        We now test all options::

            sage: def edge_options((u,v,label)):
            ...       options = { "color": { f: "red", g: "blue" }[label] }
            ...       if (u,v) == (1/2, -2): options["label"]       = "coucou"; options["label_style"] = "string"
            ...       if (u,v) == (1/2,2/3): options["dot"]         = "x=1,y=2"
            ...       if (u,v) == (1,   -1): options["label_style"] = "latex"
            ...       if (u,v) == (1,  1/2): options["edge_string"] = "<-"
            ...       if (u,v) == (1/2,  1): options["backward"]    = True
            ...       return options
            sage: print G.graphviz_string(edge_options = edge_options)
            digraph {
              node_7  [label="2/3"];
              node_5  [label="1/3"];
              node_6  [label="1/2"];
              node_9  [label="1"];
              node_4  [label="1/4"];
              node_8  [label="4/5"];
              node_0  [label="-4"];
              node_10  [label="2"];
              node_1  [label="-2"];
              node_3  [label="-1/2"];
              node_2  [label="-1"];
            <BLANKLINE>
              node_6 -> node_1 [label="coucou", color = "red"];
              node_6 -> node_7 [x=1,y=2, color = "blue"];
              node_9 -> node_2 [label=" ", texlbl="$x \ {\mapsto}\ -\frac{1}{x}$", color = "red"];
              node_9 <- node_6 [color = "blue"];
              node_4 -> node_0 [color = "red"];
              node_4 -> node_8 [color = "blue"];
              node_10 -> node_3 [color = "red"];
              node_10 -> node_5 [color = "blue"];
            }

        TESTS:

        The following digraph has tuples as vertices::

            sage: print digraphs.ButterflyGraph(1).graphviz_string()
            digraph {
              node_3  [label="('1', 1)"];
              node_0  [label="('0', 0)"];
              node_2  [label="('1', 0)"];
              node_1  [label="('0', 1)"];
            <BLANKLINE>
              node_0 -> node_3;
              node_0 -> node_1;
              node_2 -> node_3;
              node_2 -> node_1;
            }

        The following digraph has vertices with newlines in their
        string representations::

            sage: m1 = matrix(3,3)
            sage: m2 = matrix(3,3, 1)
            sage: m1.set_immutable()
            sage: m2.set_immutable()
            sage: g = DiGraph({ m1: [m2] })
            sage: print g.graphviz_string()
            digraph {
              node_0  [label="[0 0 0]\n\
            [0 0 0]\n\
            [0 0 0]"];
              node_1  [label="[1 0 0]\n\
            [0 1 0]\n\
            [0 0 1]"];
            <BLANKLINE>
              node_0 -> node_1;
            }

        REFERENCES:

        .. [dotspec] http://www.graphviz.org/doc/info/lang.html

        """
        from sage.graphs.dot2tex_utils import quoted_latex, quoted_str
        
        if self.is_directed():
            graph_string = "digraph"
            default_edge_string = "->"
        else:
            graph_string = "graph"
            default_edge_string = "--"

        edge_option_functions = options['edge_options']
        if not isinstance(edge_option_functions, (tuple,list)):
            edge_option_functions = [edge_option_functions]
        else:
            edge_option_functions = list(edge_option_functions)

        if options['edge_color'] is not None:
            default_color = options['edge_color']
        else:
            default_color = None

        if options['color_by_label'] is not False:
            color_by_label = self._color_by_label(format = options['color_by_label'], as_function = True, default_color=default_color)
            edge_option_functions.append(lambda u_v_label: {"color": color_by_label(u_v_label[2])})
        elif options['edge_colors'] is not None:
            if not isinstance(options['edge_colors'],dict):
                raise ValueError("incorrect format for edge_colors")
            color_by_edge = {}
            for color in options['edge_colors'].keys():
                for edge in options['edge_colors'][color]:
                    assert isinstance(edge, tuple) and len(edge) >= 2 and len(edge) <= 3,\
                        "%s is not a valid format for edge"%(edge)
                    u = edge[0]
                    v = edge[1]
                    assert self.has_edge(*edge), "%s is not an edge"%(edge)
                    if len(edge) == 2:
                        if self.has_multiple_edges():
                            for label in self.edge_label(u,v):
                                color_by_edge[(u,v,label)] = color
                        else:
                            label = self.edge_label(u,v)
                            color_by_edge[(u,v,label)] = color
                    elif len(edge) == 3:
                        color_by_edge[edge] = color

            edge_option_functions.append(lambda edge: {"color": color_by_edge[edge]} if edge in color_by_edge else {})

        else:
            edges_by_color = []
            not_colored_edges = self.edge_iterator(labels=True)

        key = self._keys_for_vertices()

        s = '%s {\n' % graph_string
        if (options['vertex_labels'] and
            options['labels'] == "latex"): # not a perfect option name
            # TODO: why do we set this only for latex labels?
            s += '  node [shape="plaintext"];\n'
        for v in self.vertex_iterator():
            if not options['vertex_labels']:
                node_options = ""
            elif options['labels'] == "latex":
                node_options = " [label=\" \", texlbl=\"$%s$\"]"%quoted_latex(v)
            else:
                node_options = " [label=\"%s\"]" %quoted_str(v)

            s += '  %s %s;\n'%(key(v),node_options)

        s += "\n"
        if default_color is not None:
            s += 'edge [color="%s"];\n'%default_color

        for (u,v,label) in self.edge_iterator():
            edge_options = {
                'backward': False,
                'dot': None,
                'edge_string': default_edge_string,
                'color'   : default_color,
                'label'   : label,
                'label_style': options['labels'] if options['edge_labels'] else None
                }
            for f in edge_option_functions:
                edge_options.update(f((u,v,label)))

            dot_options = []

            if edge_options['dot'] is not None:
                assert isinstance(edge_options['dot'], str)
                dot_options.append(edge_options['dot'])

            label = edge_options['label']
            if label is not None and edge_options['label_style'] is not None:
                if edge_options['label_style'] == 'latex':
                    dot_options.append('label=" ", texlbl="$%s$"'%quoted_latex(label))
                else:
                    dot_options.append('label="%s"'% label)

            if edge_options['color'] != default_color:
                dot_options.append('color = "%s"'%edge_options['color'])

            if edge_options['backward']:
                v,u = u,v
                dot_options.append('dir=back')

            s+= '  %s %s %s' % (key(u), edge_options['edge_string'], key(v))
            if len(dot_options) > 0:
                s += " [" + ", ".join(dot_options)+"]"
            s+= ";\n"
        s += "}"

        return s

    def graphviz_to_file_named(self, filename, **options):
        r"""
        Write a representation in the dot in a file.

        The dot language is a plaintext format for graph structures. See the
        documentation of :meth:`.graphviz_string` for available options.

        INPUT:

        ``filename`` - the name of the file to write in

        ``options`` - options for the graphviz string

        EXAMPLES::

            sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}},sparse=True)
            sage: tempfile = os.path.join(SAGE_TMP, 'temp_graphviz')
            sage: G.graphviz_to_file_named(tempfile, edge_labels=True)
            sage: print open(tempfile).read()
            graph {
              node_0  [label="0"];
              node_1  [label="1"];
              node_2  [label="2"];
              node_3  [label="3"];
            <BLANKLINE>
              node_0 -- node_1;
              node_0 -- node_2;
              node_1 -- node_2;
              node_2 -- node_3 [label="foo"];
            }
        """
        return open(filename, 'wt').write(self.graphviz_string(**options))

    ### Spectrum

    def spectrum(self, laplacian=False):
        r"""
        Returns a list of the eigenvalues of the adjacency matrix.

        INPUT:

        -  ``laplacian`` - if ``True``, use the Laplacian matrix
           (see :meth:`kirchhoff_matrix`)

        OUTPUT:

        A list of the eigenvalues, including multiplicities, sorted
        with the largest eigenvalue first.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.spectrum()
            [3, 1, 1, 1, 1, 1, -2, -2, -2, -2]
            sage: P.spectrum(laplacian=True)
            [5, 5, 5, 5, 2, 2, 2, 2, 2, 0]
            sage: D = P.to_directed()
            sage: D.delete_edge(7,9)
            sage: D.spectrum()
            [2.9032119259..., 1, 1, 1, 1, 0.8060634335..., -1.7092753594..., -2, -2, -2]

        ::

            sage: C = graphs.CycleGraph(8)
            sage: C.spectrum()
            [2, 1.4142135623..., 1.4142135623..., 0, 0, -1.4142135623..., -1.4142135623..., -2]

        A digraph may have complex eigenvalues.  Previously, the complex parts
        of graph eigenvalues were being dropped.  For a 3-cycle, we have::

            sage: T = DiGraph({0:[1], 1:[2], 2:[0]})
            sage: T.spectrum()
            [1, -0.5000000000... + 0.8660254037...*I, -0.5000000000... - 0.8660254037...*I]

        TESTS:

        The Laplacian matrix of a graph is the negative of the adjacency matrix with the degree of each vertex on the diagonal.  So for a regular graph, if `\delta` is an eigenvalue of a regular graph of degree `r`, then `r-\delta` will be an eigenvalue of the Laplacian.  The Hoffman-Singleton graph is regular of degree 7, so the following will test both the Laplacian construction and the computation of eigenvalues. ::

            sage: H = graphs.HoffmanSingletonGraph()
            sage: evals = H.spectrum()
            sage: lap = map(lambda x : 7 - x, evals)
            sage: lap.sort(reverse=True)
            sage: lap == H.spectrum(laplacian=True)
            True
        """
        # Ideally the spectrum should return something like a Factorization object
        # containing each eigenvalue once, along with its multiplicity.
        # This function, returning a list. could then just be renamed "eigenvalues"
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.adjacency_matrix()
        evals = M.eigenvalues()
        evals.sort(reverse=True)
        return evals

    def characteristic_polynomial(self, var='x', laplacian=False):
        r"""
        Returns the characteristic polynomial of the adjacency matrix of
        the (di)graph.

        Let `G` be a (simple) graph with adjacency matrix `A`. Let `I` be the
        identity matrix of dimensions the same as `A`. The characteristic
        polynomial of `G` is defined as the determinant `\det(xI - A)`.

        .. note::

            ``characteristic_polynomial`` and ``charpoly`` are aliases and
            thus provide exactly the same method.

        INPUT:

        - ``x`` -- (default: ``'x'``) the variable of the characteristic
          polynomial.

        - ``laplacian`` -- (default: ``False``) if ``True``, use the
          Laplacian matrix.

        .. SEEALSO::

            - :meth:`kirchhoff_matrix`

            - :meth:`laplacian_matrix`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.characteristic_polynomial()
            x^10 - 15*x^8 + 75*x^6 - 24*x^5 - 165*x^4 + 120*x^3 + 120*x^2 - 160*x + 48
            sage: P.charpoly()
            x^10 - 15*x^8 + 75*x^6 - 24*x^5 - 165*x^4 + 120*x^3 + 120*x^2 - 160*x + 48
            sage: P.characteristic_polynomial(laplacian=True)
            x^10 - 30*x^9 + 390*x^8 - 2880*x^7 + 13305*x^6 -
            39882*x^5 + 77640*x^4 - 94800*x^3 + 66000*x^2 - 20000*x
        """
        if laplacian:
            return self.kirchhoff_matrix().charpoly(var=var)
        else:
            return self.adjacency_matrix().charpoly(var=var)

    # alias, consistent with linear algebra code
    charpoly = characteristic_polynomial

    def eigenvectors(self, laplacian=False):
        r"""
        Returns the *right* eigenvectors of the adjacency matrix of the graph.

        INPUT:

        -  ``laplacian`` - if True, use the Laplacian matrix
           (see :meth:`kirchhoff_matrix`)

        OUTPUT:

        A list of triples.  Each triple begins with an eigenvalue of
        the adjacency matrix of the graph.  This is followed by
        a list of eigenvectors for the eigenvalue, when the
        eigenvectors are placed on the right side of the matrix.
        Together, the eigenvectors form a basis for the eigenspace.
        The triple concludes with the algebraic multiplicity of
        the eigenvalue.

        For some graphs, the exact eigenspaces provided by
        :meth:`eigenspaces` provide additional insight into
        the structure of the eigenspaces.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.eigenvectors()
            [(3, [
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
            ], 1), (-2, [
            (1, 0, 0, 0, -1, -1, -1, 0, 1, 1),
            (0, 1, 0, 0, -1, 0, -2, -1, 1, 2),
            (0, 0, 1, 0, -1, 1, -1, -2, 0, 2),
            (0, 0, 0, 1, -1, 1, 0, -1, -1, 1)
            ], 4), (1, [
            (1, 0, 0, 0, 0, 1, -1, 0, 0, -1),
            (0, 1, 0, 0, 0, -1, 1, -1, 0, 0),
            (0, 0, 1, 0, 0, 0, -1, 1, -1, 0),
            (0, 0, 0, 1, 0, 0, 0, -1, 1, -1),
            (0, 0, 0, 0, 1, -1, 0, 0, -1, 1)
            ], 5)]

        Eigenspaces for the Laplacian should be identical since the
        Petersen graph is regular.  However, since the output also
        contains the eigenvalues, the two outputs are slightly
        different. ::

            sage: P.eigenvectors(laplacian=True)
            [(0, [
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
            ], 1), (5, [
            (1, 0, 0, 0, -1, -1, -1, 0, 1, 1),
            (0, 1, 0, 0, -1, 0, -2, -1, 1, 2),
            (0, 0, 1, 0, -1, 1, -1, -2, 0, 2),
            (0, 0, 0, 1, -1, 1, 0, -1, -1, 1)
            ], 4), (2, [
            (1, 0, 0, 0, 0, 1, -1, 0, 0, -1),
            (0, 1, 0, 0, 0, -1, 1, -1, 0, 0),
            (0, 0, 1, 0, 0, 0, -1, 1, -1, 0),
            (0, 0, 0, 1, 0, 0, 0, -1, 1, -1),
            (0, 0, 0, 0, 1, -1, 0, 0, -1, 1)
            ], 5)]

        ::

            sage: C = graphs.CycleGraph(8)
            sage: C.eigenvectors()
            [(2, [
            (1, 1, 1, 1, 1, 1, 1, 1)
            ], 1), (-2, [
            (1, -1, 1, -1, 1, -1, 1, -1)
            ], 1), (0, [
            (1, 0, -1, 0, 1, 0, -1, 0),
            (0, 1, 0, -1, 0, 1, 0, -1)
            ], 2), (-1.4142135623..., [(1, 0, -1, 1.4142135623..., -1, 0, 1, -1.4142135623...), (0, 1, -1.4142135623..., 1, 0, -1, 1.4142135623..., -1)], 2), (1.4142135623..., [(1, 0, -1, -1.4142135623..., -1, 0, 1, 1.4142135623...), (0, 1, 1.4142135623..., 1, 0, -1, -1.4142135623..., -1)], 2)]

        A digraph may have complex eigenvalues.  Previously, the complex parts
        of graph eigenvalues were being dropped.  For a 3-cycle, we have::

            sage: T = DiGraph({0:[1], 1:[2], 2:[0]})
            sage: T.eigenvectors()
            [(1, [
            (1, 1, 1)
            ], 1), (-0.5000000000... - 0.8660254037...*I, [(1, -0.5000000000... - 0.8660254037...*I, -0.5000000000... + 0.8660254037...*I)], 1), (-0.5000000000... + 0.8660254037...*I, [(1, -0.5000000000... + 0.8660254037...*I, -0.5000000000... - 0.8660254037...*I)], 1)]
        """
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.adjacency_matrix()
        return M.right_eigenvectors()

    def eigenspaces(self, laplacian=False):
        r"""
        Returns the *right* eigenspaces of the adjacency matrix of the graph.

        INPUT:

        -  ``laplacian`` - if True, use the Laplacian matrix
           (see :meth:`kirchhoff_matrix`)

        OUTPUT:

        A list of pairs.  Each pair is an eigenvalue of the
        adjacency matrix of the graph, followed by
        the vector space that is the eigenspace for that eigenvalue,
        when the eigenvectors are placed on the right of the matrix.

        For some graphs, some of the the eigenspaces are described
        exactly by vector spaces over a
        :func:`~sage.rings.number_field.number_field.NumberField`.
        For numerical eigenvectors use :meth:`eigenvectors`.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.eigenspaces()
            [
            (3, Vector space of degree 10 and dimension 1 over Rational Field
            User basis matrix:
            [1 1 1 1 1 1 1 1 1 1]),
            (-2, Vector space of degree 10 and dimension 4 over Rational Field
            User basis matrix:
            [ 1  0  0  0 -1 -1 -1  0  1  1]
            [ 0  1  0  0 -1  0 -2 -1  1  2]
            [ 0  0  1  0 -1  1 -1 -2  0  2]
            [ 0  0  0  1 -1  1  0 -1 -1  1]),
            (1, Vector space of degree 10 and dimension 5 over Rational Field
            User basis matrix:
            [ 1  0  0  0  0  1 -1  0  0 -1]
            [ 0  1  0  0  0 -1  1 -1  0  0]
            [ 0  0  1  0  0  0 -1  1 -1  0]
            [ 0  0  0  1  0  0  0 -1  1 -1]
            [ 0  0  0  0  1 -1  0  0 -1  1])
            ]

        Eigenspaces for the Laplacian should be identical since the
        Petersen graph is regular.  However, since the output also
        contains the eigenvalues, the two outputs are slightly
        different. ::

            sage: P.eigenspaces(laplacian=True)
            [
            (0, Vector space of degree 10 and dimension 1 over Rational Field
            User basis matrix:
            [1 1 1 1 1 1 1 1 1 1]),
            (5, Vector space of degree 10 and dimension 4 over Rational Field
            User basis matrix:
            [ 1  0  0  0 -1 -1 -1  0  1  1]
            [ 0  1  0  0 -1  0 -2 -1  1  2]
            [ 0  0  1  0 -1  1 -1 -2  0  2]
            [ 0  0  0  1 -1  1  0 -1 -1  1]),
            (2, Vector space of degree 10 and dimension 5 over Rational Field
            User basis matrix:
            [ 1  0  0  0  0  1 -1  0  0 -1]
            [ 0  1  0  0  0 -1  1 -1  0  0]
            [ 0  0  1  0  0  0 -1  1 -1  0]
            [ 0  0  0  1  0  0  0 -1  1 -1]
            [ 0  0  0  0  1 -1  0  0 -1  1])
            ]

        Notice how one eigenspace below is described with a square root of
        2.  For the two possible values (positive and negative) there is a
        corresponding eigenspace.  ::

            sage: C = graphs.CycleGraph(8)
            sage: C.eigenspaces()
            [
            (2, Vector space of degree 8 and dimension 1 over Rational Field
            User basis matrix:
            [1 1 1 1 1 1 1 1]),
            (-2, Vector space of degree 8 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -1  1 -1  1 -1  1 -1]),
            (0, Vector space of degree 8 and dimension 2 over Rational Field
            User basis matrix:
            [ 1  0 -1  0  1  0 -1  0]
            [ 0  1  0 -1  0  1  0 -1]),
            (a3, Vector space of degree 8 and dimension 2 over Number Field in a3 with defining polynomial x^2 - 2
            User basis matrix:
            [  1   0  -1 -a3  -1   0   1  a3]
            [  0   1  a3   1   0  -1 -a3  -1])
            ]

        A digraph may have complex eigenvalues and eigenvectors.
        For a 3-cycle, we have::

            sage: T = DiGraph({0:[1], 1:[2], 2:[0]})
            sage: T.eigenspaces()
            [
            (1, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [1 1 1]),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 + x + 1
            User basis matrix:
            [      1      a1 -a1 - 1])
            ]
        """
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.adjacency_matrix()
        # could pass format='all' to get QQbar eigenvalues and eigenspaces
        # which would be a change in default behavior
        return M.right_eigenspaces(format='galois', algebraic_multiplicity=False)

    ### Automorphism and isomorphism

    def relabel(self, perm=None, inplace=True, return_map=False, check_input = True, complete_partial_function = True):
        r"""
        Relabels the vertices of ``self``

        INPUT:

         - ``perm`` -- a function, dictionary, list, permutation, or
           ``None`` (default: ``None``)

         - ``inplace`` -- a boolean (default: ``True``)

         - ``return_map`` -- a boolean (default: ``False``)

         - ``check_input`` (boolean) -- whether to test input for
           correctness. *This can potentially be very time-consuming !*.

         - ``complete_partial_function`` (boolean) -- whether to automatically
           complete the permutation if some elements of the graph are not
           associated with any new name. In this case, those elements are not
           relabeled *This can potentially be very time-consuming !*.

        If ``perm`` is a function ``f``, then each vertex ``v`` is
        relabeled to ``f(v)``.

        If ``perm`` is a dictionary ``d``, then each vertex ``v`` (which should
        be a key of ``d``) is relabeled to ``d[v]``.  Similarly, if ``perm`` is
        a list or tuple ``l`` of length ``n``, then the first vertex returned by
        ``G.vertices()`` is relabeled to ``l[0]``, the second to ``l[1]``, ...

        If ``perm`` is a permutation, then each vertex ``v`` is
        relabeled to ``perm(v)``. Caveat: this assumes that the
        vertices are labelled `\{0,1,...,n-1\}`; since permutations
        act by default on the set `\{1,2,...,n\}`, this is achieved by
        identifying `n` and `0`.

        If ``perm`` is ``None``, the graph is relabeled to be on the
        vertices `\{0,1,...,n-1\}`.

        .. note:: at this point, only injective relabeling are supported.

        If ``inplace`` is ``True``, the graph is modified in place and
        ``None`` is returned. Otherwise a relabeled copy of the graph
        is returned.

        If ``return_map`` is ``True`` a dictionary representing the
        relabelling map is returned (incompatible with ``inplace==False``).

        EXAMPLES::

            sage: G = graphs.PathGraph(3)
            sage: G.am()
            [0 1 0]
            [1 0 1]
            [0 1 0]

        Relabeling using a dictionary. Note that the dictionary does not define
        the new label of vertex `0`::

            sage: G.relabel({1:2,2:1}, inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        This is because the method automatically "extends" the relabeling to the
        missing vertices (whose label will not change). Checking that all
        vertices have an image can require some time, and this feature can be
        disabled (at your own risk)::

            sage: G.relabel({1:2,2:1}, inplace=False, complete_partial_function = False).am()
            Traceback (most recent call last):
            ...
            KeyError: 0

        Relabeling using a list::

            sage: G.relabel([0,2,1], inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using a tuple::

            sage: G.relabel((0,2,1), inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using a Sage permutation::

            sage: G = graphs.PathGraph(3)
            sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            sage: S = SymmetricGroup(3)
            sage: gamma = S('(1,2)')
            sage: G.relabel(gamma, inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using an injective function::

            sage: G.edges()
            [(0, 1, None), (1, 2, None)]
            sage: H = G.relabel(lambda i: i+10, inplace=False)
            sage: H.vertices()
            [10, 11, 12]
            sage: H.edges()
            [(10, 11, None), (11, 12, None)]

        Relabeling using a non injective function has no meaning::

            sage: G.edges()
            [(0, 1, None), (1, 2, None)]
            sage: G.relabel(lambda i: 0, inplace=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: Non injective relabeling

        But this test can be disabled, which can lead to ... problems::

            sage: G.edges()
            [(0, 1, None), (1, 2, None)]
            sage: G.relabel(lambda i: 0, check_input = False)
            sage: G.edges()
            []

        Relabeling to simpler labels::

            sage: G = graphs.CubeGraph(3)
            sage: G.vertices()
            ['000', '001', '010', '011', '100', '101', '110', '111']
            sage: G.relabel()
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7]

        Recovering the relabeling with ``return_map``::

            sage: G = graphs.CubeGraph(3)
            sage: expecting = {'000': 0, '001': 1, '010': 2, '011': 3, '100': 4, '101': 5, '110': 6, '111': 7}
            sage: G.relabel(return_map=True) == expecting
            True

        ::

            sage: G = graphs.PathGraph(3)
            sage: G.relabel(lambda i: i+10, return_map=True)
            {0: 10, 1: 11, 2: 12}

        TESTS::

            sage: P = Graph(graphs.PetersenGraph())
            sage: P.delete_edge([0,1])
            sage: P.add_edge((4,5))
            sage: P.add_edge((2,6))
            sage: P.delete_vertices([0,1])
            sage: P.relabel()

        The attributes are properly updated too

        ::

            sage: G = graphs.PathGraph(5)
            sage: G.set_vertices({0: 'before', 1: 'delete', 2: 'after'})
            sage: G.set_boundary([1,2,3])
            doctest:...: DeprecationWarning: The boundary parameter is deprecated and will soon disappear.
            See http://trac.sagemath.org/15494 for details.
            sage: G.delete_vertex(1)
            sage: G.relabel()
            sage: G.get_vertices()
            {0: 'before', 1: 'after', 2: None, 3: None}
            sage: G.get_boundary()
            [1, 2]
            sage: G.get_pos()
            {0: (0, 0), 1: (2, 0), 2: (3, 0), 3: (4, 0)}

        Check that :trac:`12477` is fixed::

            sage: g = Graph({1:[2,3]})
            sage: rel = {1:'a', 2:'b'}
            sage: g.relabel(rel)
            sage: g.vertices()
            [3, 'a', 'b']
            sage: rel
            {1: 'a', 2: 'b'}

        Immutable graphs cannot be relabeled::

            sage: Graph(graphs.PetersenGraph(), immutable=True).relabel({})
            Traceback (most recent call last):
            ...
            ValueError: To relabel an immutable graph use inplace=False

        :trac:`16257`::

            sage: G = graphs.PetersenGraph()
            sage: G.relabel( [ i+1 for i in range(G.order()) ], inplace=True )
            sage: G.relabel( [ i+1 for i in range(G.order()) ], inplace=True )
        """
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement

        if not inplace:
            G = self.copy(immutable=False)
            perm2 = G.relabel(perm,
                              return_map= return_map,
                              check_input = check_input,
                              complete_partial_function = complete_partial_function)

            if getattr(self, "_immutable", False):
                G = self.__class__(G, immutable = True)

            if return_map:
                return G, perm2
            else:
                return G

        if getattr(self, "_immutable", False):
            raise ValueError("To relabel an immutable graph use inplace=False")

        # If perm is not a dictionary, we build one !

        if perm is None:

            # vertices() returns a sorted list:
            # this guarantees consistent relabeling
            perm = dict({v:i for i,v in enumerate(self.vertices())})
            complete_partial_function = False
            check_input = False

        elif isinstance(perm, dict):

            # If all vertices do not have a new label, the code will touch the
            # dictionary. Let us keep the one we received from the user clean !
            from copy import copy
            perm = copy(perm)

        elif isinstance(perm, (list, tuple)):
            perm = dict(zip(self.vertices(),perm))

        elif isinstance(perm, PermutationGroupElement):
            n = self.order()
            ddict = {}
            for i in xrange(1,n):
                ddict[i] = perm(i)%n
            if n > 0:
                ddict[0] = perm(n)%n
            perm = ddict

        elif callable(perm):
            perm = dict( [ i, perm(i) ] for i in self.vertices() )
            complete_partial_function = False

        else:
            raise TypeError("Type of perm is not supported for relabeling.")

        # Whether to complete the relabeling function if some vertices do not
        # appear in the permutation.
        if complete_partial_function:
            for v in self:
                if v not in perm:
                    perm[v] = v

        # Whether to check input
        if check_input:
            if len(set(perm.values())) < len(perm):
                raise NotImplementedError("Non injective relabeling")

            for v in perm.iterkeys():
                if v in self:
                    try:
                        hash(perm[v])
                    except TypeError:
                        raise ValueError("perm dictionary must be of the format {a:a1, b:b1, ...} where a,b,... are vertices and a1,b1,... are hashable")

        self._backend.relabel(perm, self._directed)

        attributes_to_update = ('_pos', '_assoc', '_embedding')
        for attr in attributes_to_update:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                new_attr = {}
                for v,value in getattr(self, attr).iteritems():
                    new_attr[perm[v]] = value

                setattr(self, attr, new_attr)

        self._boundary = [perm[v] for v in self._boundary]

        if return_map:
            return perm

    def degree_to_cell(self, vertex, cell):
        """
        Returns the number of edges from vertex to an edge in cell. In the
        case of a digraph, returns a tuple (in_degree, out_degree).

        EXAMPLES::

            sage: G = graphs.CubeGraph(3)
            sage: cell = G.vertices()[:3]
            sage: G.degree_to_cell('011', cell)
            2
            sage: G.degree_to_cell('111', cell)
            0

        ::

            sage: D = DiGraph({ 0:[1,2,3], 1:[3,4], 3:[4,5]})
            sage: cell = [0,1,2]
            sage: D.degree_to_cell(5, cell)
            (0, 0)
            sage: D.degree_to_cell(3, cell)
            (2, 0)
            sage: D.degree_to_cell(0, cell)
            (0, 2)
        """
        if self._directed:
            in_neighbors_in_cell = set([a for a,_,_ in self.incoming_edges(vertex)]) & set(cell)
            out_neighbors_in_cell = set([a for _,a,_ in self.outgoing_edges(vertex)]) & set(cell)
            return (len(in_neighbors_in_cell), len(out_neighbors_in_cell))
        else:
            neighbors_in_cell = set(self.neighbors(vertex)) & set(cell)
            return len(neighbors_in_cell)

    def is_equitable(self, partition, quotient_matrix=False):
        """
        Checks whether the given partition is equitable with respect to
        self.

        A partition is equitable with respect to a graph if for every pair
        of cells C1, C2 of the partition, the number of edges from a vertex
        of C1 to C2 is the same, over all vertices in C1.

        INPUT:


        -  ``partition`` - a list of lists

        -  ``quotient_matrix`` - (default False) if True, and
           the partition is equitable, returns a matrix over the integers
           whose rows and columns represent cells of the partition, and whose
           i,j entry is the number of vertices in cell j adjacent to each
           vertex in cell i (since the partition is equitable, this is well
           defined)


        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8],[7]])
            False
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8,7]])
            True
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8,7]], quotient_matrix=True)
            [1 2 0]
            [1 0 2]
            [0 2 1]

        ::

            sage: ss = (graphs.WheelGraph(6)).line_graph(labels=False)
            sage: prt = [[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]

        ::

            sage: ss.is_equitable(prt)
            Traceback (most recent call last):
            ...
            TypeError: Partition ([[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]) is not valid for this graph: vertices are incorrect.

        ::

            sage: ss = (graphs.WheelGraph(5)).line_graph(labels=False)
            sage: ss.is_equitable(prt)
            False
        """
        from sage.misc.flatten import flatten
        from sage.misc.misc import uniq
        if sorted(flatten(partition, max_level=1)) != self.vertices():
            raise TypeError("Partition (%s) is not valid for this graph: vertices are incorrect."%partition)
        if any(len(cell)==0 for cell in partition):
            raise TypeError("Partition (%s) is not valid for this graph: there is a cell of length 0."%partition)
        if quotient_matrix:
            from sage.matrix.constructor import Matrix
            from sage.rings.integer_ring import IntegerRing
            n = len(partition)
            M = Matrix(IntegerRing(), n)
            for i in xrange(n):
                for j in xrange(n):
                    cell_i = partition[i]
                    cell_j = partition[j]
                    degrees = [self.degree_to_cell(u, cell_j) for u in cell_i]
                    if len(uniq(degrees)) > 1:
                        return False
                    if self._directed:
                        M[i, j] = degrees[0][0]
                    else:
                        M[i, j] = degrees[0]
            return M
        else:
            for cell1 in partition:
                for cell2 in partition:
                    degrees = [self.degree_to_cell(u, cell2) for u in cell1]
                    if len(uniq(degrees)) > 1:
                        return False
            return True

    def coarsest_equitable_refinement(self, partition, sparse=True):
        """
        Returns the coarsest partition which is finer than the input
        partition, and equitable with respect to self.

        A partition is equitable with respect to a graph if for every pair
        of cells C1, C2 of the partition, the number of edges from a vertex
        of C1 to C2 is the same, over all vertices in C1.

        A partition P1 is finer than P2 (P2 is coarser than P1) if every
        cell of P1 is a subset of a cell of P2.

        INPUT:


        -  ``partition`` - a list of lists

        -  ``sparse`` - (default False) whether to use sparse
           or dense representation- for small graphs, use dense for speed


        EXAMPLES::

            sage: G = graphs.PetersenGraph()
            sage: G.coarsest_equitable_refinement([[0],range(1,10)])
            [[0], [2, 3, 6, 7, 8, 9], [1, 4, 5]]
            sage: G = graphs.CubeGraph(3)
            sage: verts = G.vertices()
            sage: Pi = [verts[:1], verts[1:]]
            sage: Pi
            [['000'], ['001', '010', '011', '100', '101', '110', '111']]
            sage: G.coarsest_equitable_refinement(Pi)
            [['000'], ['011', '101', '110'], ['111'], ['001', '010', '100']]

        Note that given an equitable partition, this function returns that
        partition::

            sage: P = graphs.PetersenGraph()
            sage: prt = [[0], [1, 4, 5], [2, 3, 6, 7, 8, 9]]
            sage: P.coarsest_equitable_refinement(prt)
            [[0], [1, 4, 5], [2, 3, 6, 7, 8, 9]]

        ::

            sage: ss = (graphs.WheelGraph(6)).line_graph(labels=False)
            sage: prt = [[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]
            sage: ss.coarsest_equitable_refinement(prt)
            Traceback (most recent call last):
            ...
            TypeError: Partition ([[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]) is not valid for this graph: vertices are incorrect.

        ::

            sage: ss = (graphs.WheelGraph(5)).line_graph(labels=False)
            sage: ss.coarsest_equitable_refinement(prt)
            [[(0, 1)], [(1, 2), (1, 4)], [(0, 3)], [(0, 2), (0, 4)], [(2, 3), (3, 4)]]

        ALGORITHM: Brendan D. McKay's Master's Thesis, University of
        Melbourne, 1976.
        """
        from sage.misc.flatten import flatten
        if sorted(flatten(partition, max_level=1)) != self.vertices():
            raise TypeError("Partition (%s) is not valid for this graph: vertices are incorrect."%partition)
        if any(len(cell)==0 for cell in partition):
            raise TypeError("Partition (%s) is not valid for this graph: there is a cell of length 0."%partition)
        if self.has_multiple_edges():
            raise TypeError("Refinement function does not support multiple edges.")
        from copy import copy
        G = copy(self)
        perm_to = G.relabel(return_map=True)
        partition = [[perm_to[b] for b in cell] for cell in partition]
        perm_from = {}
        for v in self:
            perm_from[perm_to[v]] = v
        n = G.num_verts()
        if sparse:
            from sage.graphs.base.sparse_graph import SparseGraph
            CG = SparseGraph(n)
        else:
            from sage.graphs.base.dense_graph import DenseGraph
            CG = DenseGraph(n)
        for i in range(n):
            for j in range(n):
                if G.has_edge(i,j):
                    CG.add_arc(i,j)

        from sage.groups.perm_gps.partn_ref.refinement_graphs import coarsest_equitable_refinement
        result = coarsest_equitable_refinement(CG, partition, G._directed)
        return [[perm_from[b] for b in cell] for cell in result]

    def automorphism_group(self, partition=None, verbosity=0,
                           edge_labels=False, order=False,
                           return_group=True, orbits=False):
        """
        Returns the largest subgroup of the automorphism group of the
        (di)graph whose orbit partition is finer than the partition given.
        If no partition is given, the unit partition is used and the entire
        automorphism group is given.

        INPUT:

        -  ``partition`` - default is the unit partition,
           otherwise computes the subgroup of the full automorphism group
           respecting the partition.

        -  ``edge_labels`` - default False, otherwise allows
           only permutations respecting edge labels.

        -  ``order`` - (default False) if True, compute the
           order of the automorphism group

        -  ``return_group`` - default True

        -  ``orbits`` - returns the orbits of the group acting
           on the vertices of the graph

        .. WARNING::

            Since :trac:`14319` the domain of the automorphism group is equal to
            the graph's vertex set, and the ``translation`` argument has become
            useless.

        OUTPUT: The order of the output is group, order, orbits. However, there
        are options to turn each of these on or off.

        EXAMPLES:

        Graphs::

            sage: graphs_query = GraphQuery(display_cols=['graph6'],num_vertices=4)
            sage: L = graphs_query.get_graphs_list()
            sage: graphs_list.show_graphs(L)
            sage: for g in L:
            ...    G = g.automorphism_group()
            ...    G.order(), G.gens()
            (24, [(2,3), (1,2), (0,1)])
            (4, [(2,3), (0,1)])
            (2, [(1,2)])
            (6, [(1,2), (0,1)])
            (6, [(2,3), (1,2)])
            (8, [(1,2), (0,1)(2,3)])
            (2, [(0,1)(2,3)])
            (2, [(1,2)])
            (8, [(2,3), (0,1), (0,2)(1,3)])
            (4, [(2,3), (0,1)])
            (24, [(2,3), (1,2), (0,1)])
            sage: C = graphs.CubeGraph(4)
            sage: G = C.automorphism_group()
            sage: M = G.character_table() # random order of rows, thus abs() below
            sage: QQ(M.determinant()).abs()
            712483534798848
            sage: G.order()
            384

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: G = D.automorphism_group()
            sage: A5 = AlternatingGroup(5)
            sage: Z2 = CyclicPermutationGroup(2)
            sage: H = A5.direct_product(Z2)[0] #see documentation for direct_product to explain the [0]
            sage: G.is_isomorphic(H)
            True

        Multigraphs::

            sage: G = Graph(multiedges=True,sparse=True)
            sage: G.add_edge(('a', 'b'))
            sage: G.add_edge(('a', 'b'))
            sage: G.add_edge(('a', 'b'))
            sage: G.automorphism_group()
            Permutation Group with generators [('a','b')]

        Digraphs::

            sage: D = DiGraph( { 0:[1], 1:[2], 2:[3], 3:[4], 4:[0] } )
            sage: D.automorphism_group()
            Permutation Group with generators [(0,1,2,3,4)]

        Edge labeled graphs::

            sage: G = Graph(sparse=True)
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: G.automorphism_group(edge_labels=True)
            Permutation Group with generators [(1,4)(2,3)]

        ::

            sage: G = Graph({0 : {1 : 7}})
            sage: G.automorphism_group(edge_labels=True)
            Permutation Group with generators [(0,1)]

            sage: foo = Graph(sparse=True)
            sage: bar = Graph(implementation='c_graph',sparse=True)
            sage: foo.add_edges([(0,1,1),(1,2,2), (2,3,3)])
            sage: bar.add_edges([(0,1,1),(1,2,2), (2,3,3)])
            sage: foo.automorphism_group(edge_labels=True)
            Permutation Group with generators [()]
            sage: foo.automorphism_group()
            Permutation Group with generators [(0,3)(1,2)]
            sage: bar.automorphism_group(edge_labels=True)
            Permutation Group with generators [()]

        You can also ask for just the order of the group::

            sage: G = graphs.PetersenGraph()
            sage: G.automorphism_group(return_group=False, order=True)
            120

        Or, just the orbits (note that each graph here is vertex transitive)

        ::

            sage: G = graphs.PetersenGraph()
            sage: G.automorphism_group(return_group=False, orbits=True)
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
            sage: G.automorphism_group(partition=[[0],range(1,10)], return_group=False, orbits=True)
            [[0], [2, 3, 6, 7, 8, 9], [1, 4, 5]]
            sage: C = graphs.CubeGraph(3)
            sage: C.automorphism_group(orbits=True, return_group=False)
            [['000', '001', '010', '011', '100', '101', '110', '111']]

        TESTS:

        We get a KeyError when given an invalid partition (trac #6087)::

            sage: g=graphs.CubeGraph(3)
            sage: g.relabel()
            sage: g.automorphism_group(partition=[[0,1,2],[3,4,5]])
            Traceback (most recent call last):
            ...
            KeyError: 6

        Labeled automorphism group::

            sage: digraphs.DeBruijn(3,2).automorphism_group()
            Permutation Group with generators [('01','02')('10','20')('11','22')('12','21'), ('00','11')('01','10')('02','12')('20','21')]
            sage: d = digraphs.DeBruijn(3,2)
            sage: d.allow_multiple_edges(True)
            sage: d.add_edge(d.edges()[0])
            sage: d.automorphism_group()
            Permutation Group with generators [('01','02')('10','20')('11','22')('12','21')]

        The labeling is correct::

            sage: g = graphs.PetersenGraph()
            sage: ag = g.automorphism_group()
            sage: for u,v in g.edges(labels = False):
            ...       if len(ag.orbit((u,v),action="OnPairs")) != 30:
            ...           print "ARggggggggggggg !!!"

        Empty group, correct domain::

            sage: Graph({'a':['a'], 'b':[]}).automorphism_group()
            Permutation Group with generators [()]
            sage: Graph({'a':['a'], 'b':[]}).automorphism_group().domain()
            {'a', 'b'}

        We can check that the subgroups are labelled correctly
        (:trac:`15656`)::

            sage: G1 = Graph(':H`ECw@HGXGAGUG`e')
            sage: G = G1.automorphism_group()
            sage: G.subgroups()
            [Subgroup of (Permutation Group with generators [(0,7)(1,4)(2,3)(6,8)]) generated by [()],
            Subgroup of (Permutation Group with generators [(0,7)(1,4)(2,3)(6,8)]) generated by [(0,7)(1,4)(2,3)(6,8)]]

        """
        from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        from sage.groups.perm_gps.permgroup import PermutationGroup
        dig = (self._directed or self.has_loops())
        if partition is None:
            partition = [self.vertices()]
        if edge_labels or self.has_multiple_edges():
            G, partition, relabeling = graph_isom_equivalent_non_edge_labeled_graph(self, partition, return_relabeling=True, ignore_edge_labels=(not edge_labels))
            G_vertices = sum(partition, [])
            G_to = {}
            for i in xrange(len(G_vertices)):
                G_to[G_vertices[i]] = i
            from sage.graphs.all import Graph, DiGraph
            DoDG = DiGraph if self._directed else Graph
            H = DoDG(len(G_vertices), implementation='c_graph', loops=G.allows_loops())
            HB = H._backend
            for u,v in G.edge_iterator(labels=False):
                u = G_to[u]; v = G_to[v]
                HB.add_edge(u,v,None,G._directed)
            GC = HB._cg
            partition = [[G_to[v] for v in cell] for cell in partition]
            A = search_tree(GC, partition, lab=False, dict_rep=True, dig=dig, verbosity=verbosity, order=order)
            if order:
                a,b,c = A
            else:
                a,b = A
            b_new = {}
            for v in G_to:
                b_new[v] = b[G_to[v]]
            b = b_new
            # b is a translation of the labellings
            acting_vertices = {}
            translation_d = {}
            m = G.order()
            for v in self:
                if b[relabeling[v]] == m:
                    translation_d[v] = self.order()
                    acting_vertices[v] = 0
                else:
                    translation_d[v] = b[relabeling[v]]
                    acting_vertices[v] = b[relabeling[v]]
            real_aut_gp = []
            n = self.order()
            for gen in a:
                gen_restr = [0]*n
                for v in self.vertex_iterator():
                    gen_restr[acting_vertices[v]] = gen[acting_vertices[v]]
                if gen_restr not in real_aut_gp:
                    real_aut_gp.append(gen_restr)
            id = range(n)
            if id in real_aut_gp:
                real_aut_gp.remove(id)
            a = real_aut_gp
            b = translation_d
        else:
            G_vertices = sum(partition, [])
            G_to = {}
            for i in xrange(len(G_vertices)):
                G_to[G_vertices[i]] = i
            from sage.graphs.all import Graph, DiGraph
            DoDG = DiGraph if self._directed else Graph
            H = DoDG(len(G_vertices), implementation='c_graph', loops=self.allows_loops())
            HB = H._backend
            for u,v in self.edge_iterator(labels=False):
                u = G_to[u]; v = G_to[v]
                HB.add_edge(u,v,None,self._directed)
            GC = HB._cg
            partition = [[G_to[v] for v in cell] for cell in partition]

            if return_group:
                A = search_tree(GC, partition, dict_rep=True, lab=False, dig=dig, verbosity=verbosity, order=order)
                if order:
                    a,b,c = A
                else:
                    a,b = A
                b_new = {}
                for v in G_to:
                    b_new[v] = b[G_to[v]]
                b = b_new
            else:
                a = search_tree(GC, partition, dict_rep=False, lab=False, dig=dig, verbosity=verbosity, order=order)
                if order:
                    a,c = a

        output = []
        if return_group:
            if len(a) != 0:
                # We translate the integer permutations into a collection of
                # cycles.
                from sage.combinat.permutation import Permutation
                gens = [Permutation([x+1 for x in aa]).to_cycles() for aa in a]

                # We relabel the cycles using the vertices' names instead of integers
                n = self.order()
                int_to_vertex = {((i+1) if i != n else 1):v for v,i in b.iteritems()}
                gens = [ [ tuple([int_to_vertex[i] for i in cycle]) for cycle in gen] for gen in gens]
                output.append(PermutationGroup(gens = gens, domain = int_to_vertex.values()))
            else:
                output.append(PermutationGroup([[]], domain = self.vertices()))
        if order:
            output.append(c)
        if orbits:
            G_from = {}
            for v in G_to:
                G_from[G_to[v]] = v
            from sage.groups.perm_gps.partn_ref.refinement_graphs import get_orbits
            output.append([[G_from[v] for v in W] for W in get_orbits(a, self.num_verts())])

        # A Python switch statement!
        return { 0: None,
                 1: output[0],
                 2: tuple(output),
                 3: tuple(output),
                 4: tuple(output)
               }[len(output)]

    def is_vertex_transitive(self, partition=None, verbosity=0,
                           edge_labels=False, order=False,
                           return_group=True, orbits=False):
        """
        Returns whether the automorphism group of self is transitive within
        the partition provided, by default the unit partition of the
        vertices of self (thus by default tests for vertex transitivity in
        the usual sense).

        EXAMPLES::

            sage: G = Graph({0:[1],1:[2]})
            sage: G.is_vertex_transitive()
            False
            sage: P = graphs.PetersenGraph()
            sage: P.is_vertex_transitive()
            True
            sage: D = graphs.DodecahedralGraph()
            sage: D.is_vertex_transitive()
            True
            sage: R = graphs.RandomGNP(2000, .01)
            sage: R.is_vertex_transitive()
            False
        """
        if partition is None:
            partition = [self.vertices()]

        for p in partition:
            if len(p) == 0:
                continue
            d = self.degree(p[0])
            if not all(self.degree(x) == d for x in p):
                return False

        new_partition = self.automorphism_group(partition,
                          verbosity=verbosity, edge_labels=edge_labels,
                          order=False, return_group=False, orbits=True)

        return (len(partition) == len(new_partition))

    def is_hamiltonian(self):
        r"""
        Tests whether the current graph is Hamiltonian.

        A graph (resp. digraph) is said to be Hamiltonian
        if it contains as a subgraph a cycle (resp. a circuit)
        going through all the vertices.

        Testing for Hamiltonicity being NP-Complete, this
        algorithm could run for some time depending on
        the instance.

        ALGORITHM:

        See ``Graph.traveling_salesman_problem``.

        OUTPUT:

        Returns ``True`` if a Hamiltonian cycle/circuit exists, and
        ``False`` otherwise.

        NOTE:

        This function, as ``hamiltonian_cycle`` and
        ``traveling_salesman_problem``, computes a Hamiltonian
        cycle if it exists: the user should *NOT* test for
        Hamiltonicity using ``is_hamiltonian`` before calling
        ``hamiltonian_cycle`` or ``traveling_salesman_problem``
        as it would result in computing it twice.

        EXAMPLES:

        The Heawood Graph is known to be Hamiltonian ::

            sage: g = graphs.HeawoodGraph()
            sage: g.is_hamiltonian()
            True

        The Petergraph, though, is not ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_hamiltonian()
            False

        TESTS:

        When no solver is installed, a
        ``OptionalPackageNotFoundError`` exception is raised::

            sage: from sage.misc.exceptions import OptionalPackageNotFoundError
            sage: try:
            ...       g = graphs.ChvatalGraph()
            ...       if not g.is_hamiltonian():
            ...          print "There is something wrong here !"
            ... except OptionalPackageNotFoundError:
            ...       pass

        :trac:`16210`::

            sage: g=graphs.CycleGraph(10)
            sage: g.allow_loops(True)
            sage: g.add_edge(0,0)
            sage: g.is_hamiltonian()
            True
        """
        from sage.categories.sets_cat import EmptySetError
        try:
            tsp = self.traveling_salesman_problem(use_edge_labels = False)
            return True

        except EmptySetError:
            return False

    def is_isomorphic(self, other, certify=False, verbosity=0, edge_labels=False):
        """
        Tests for isomorphism between self and other.

        INPUT:


        -  ``certify`` - if True, then output is (a,b), where a
           is a boolean and b is either a map or None.

        -  ``edge_labels`` - default False, otherwise allows
           only permutations respecting edge labels.


        EXAMPLES: Graphs::

            sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            sage: D = graphs.DodecahedralGraph()
            sage: E = copy(D)
            sage: gamma = SymmetricGroup(20).random_element()
            sage: E.relabel(gamma)
            sage: D.is_isomorphic(E)
            True

        ::

            sage: D = graphs.DodecahedralGraph()
            sage: S = SymmetricGroup(20)
            sage: gamma = S.random_element()
            sage: E = copy(D)
            sage: E.relabel(gamma)
            sage: a,b = D.is_isomorphic(E, certify=True); a
            True
            sage: from sage.plot.graphics import GraphicsArray
            sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
            sage: position_D = spring_layout_fast(D)
            sage: position_E = {}
            sage: for vert in position_D:
            ...    position_E[b[vert]] = position_D[vert]
            sage: GraphicsArray([D.plot(pos=position_D), E.plot(pos=position_E)]).show() # long time

        ::

            sage: g=graphs.HeawoodGraph()
            sage: g.is_isomorphic(g)
            True

        Multigraphs::

            sage: G = Graph(multiedges=True,sparse=True)
            sage: G.add_edge((0,1,1))
            sage: G.add_edge((0,1,2))
            sage: G.add_edge((0,1,3))
            sage: G.add_edge((0,1,4))
            sage: H = Graph(multiedges=True,sparse=True)
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: G.is_isomorphic(H)
            True

        Digraphs::

            sage: A = DiGraph( { 0 : [1,2] } )
            sage: B = DiGraph( { 1 : [0,2] } )
            sage: A.is_isomorphic(B, certify=True)
            (True, {0: 1, 1: 0, 2: 2})

        Edge labeled graphs::

            sage: G = Graph(sparse=True)
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: H = G.relabel([1,2,3,4,0], inplace=False)
            sage: G.is_isomorphic(H, edge_labels=True)
            True

        Edge labeled digraphs::

            sage: G = DiGraph()
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: H = G.relabel([1,2,3,4,0], inplace=False)
            sage: G.is_isomorphic(H, edge_labels=True)
            True
            sage: G.is_isomorphic(H, edge_labels=True, certify=True)
            (True, {0: 1, 1: 2, 2: 3, 3: 4, 4: 0})

        TESTS::

            sage: g1 = '~?A[~~{ACbCwV_~__OOcCW_fAA{CF{CCAAAC__bCCCwOOV___~____OOOOcCCCW___fAAAA'+\
            ...   '{CCCF{CCCCAAAAAC____bCCCCCwOOOOV_____~_O@ACG_@ACGOo@ACG?{?`A?GV_GO@AC}@?_OGC'+\
            ...   'C?_OI@?K?I@?_OM?_OGD?F_A@OGC@{A@?_OG?O@?gCA?@_GCA@O?B_@OGCA?BoA@?gC?@{A?GO`?'+\
            ...   '??_GO@AC??E?O`?CG??[?O`A?G??{?GO`A???|A?_GOC`AC@_OCGACEAGS?HA?_SA`aO@G?cOC_N'+\
            ...   'G_C@AOP?GnO@_GACOE?g?`OGACCOGaGOc?HA?`GORCG_AO@B?K@[`A?OCI@A@By?_K@?SCABA?H?'+\
            ...   'SA?a@GC`CH?Q?C_c?cGRC@G_AOCOa@Ax?QC?_GOo_CNg@A?oC@CaCGO@CGA_O`?GSGPAGOC_@OO_'+\
            ...   'aCHaG?cO@CB?_`Ax?GQC?_cAOCG^OGAC@_D?IGO`?D?O_I?HAOO`AGOHA?cC?oAO`AW_Q?HCACAC'+\
            ...   'GO`[_OCHA?_cCACG^O_@CAGO`A?GCOGc@?I?OQOC?IGC_o@CAGCCE?A@DBG_OA@C_CP?OG_VA_CO'+\
            ...   'G@D?_OA_DFgA@CO?aH?Ga@?a?_I?S@A@@Oa@?@P@GCO_AACO_a_?`K_GCQ@?cAOG_OGAwQ@?K?cC'+\
            ...   'GH?I?ABy@C?G_S@@GCA@C`?OI?_D?OP@G?IGGP@O_AGCP?aG?GCPAX?cA?OGSGCGCAGCJ`?oAGCC'+\
            ...   'HAA?A_CG^O@CAG_GCSCAGCCGOCG@OA_`?`?g_OACG_`CAGOAO_H?a_?`AXA?OGcAAOP?a@?CGVAC'+\
            ...   'OG@_AGG`OA_?O`|?Ga?COKAAGCA@O`A?a?S@?HCG`?_?gO`AGGaC?PCAOGI?A@GO`K_CQ@?GO_`O'+\
            ...   'GCAACGVAG@_COOCQ?g?I?O`ByC?G_P?O`A?H@G?_P?`OAGC?gD?_C@_GCAGDG_OA@CCPC?AOQ??g'+\
            ...   '_R@_AGCO____OCC_@OAbaOC?g@C_H?AOOC@?a`y?PC?G`@OOH??cOG_OOAG@_COAP?WA?_KAGC@C'+\
            ...   '_CQ@?HAACH??c@P?_AWGaC?P?gA_C_GAD?I?Awa?S@?K?`C_GAOGCS?@|?COGaA@CAAOQ?AGCAGO'+\
            ...   'ACOG@_G_aC@_G@CA@@AHA?OGc?WAAH@G?P?_?cH_`CAGOGACc@@GA?S?CGVCG@OA_CICAOOC?PO?'+\
            ...   'OG^OG_@CAC_cC?AOP?_OICG@?oAGCO_GO_GB@?_OG`AH?cA?OH?`P??cC_O?SCGR@O_AGCAI?Q?_'+\
            ...   'GGS?D?O`[OI?_D@@CCA?cCA_?_O`By?_PC?IGAGOQ?@A@?aO`A?Q@?K?__`_E?_GCA@CGO`C_GCQ'+\
            ...   '@A?gAOQ?@C?DCACGR@GCO_AGPA@@GAA?A_CO`Aw_I?S@?SCB@?OC_?_P@ACNgOC@A?aCGOCAGCA@'+\
            ...   'CA?H@GG_C@AOGa?OOG_O?g_OA?oDC_AO@GOCc?@P?_A@D??cC``O?cGAOGD?@OA_CAGCA?_cwKA?'+\
            ...   '`?OWGG?_PO?I?S?H@?^OGAC@_Aa@CAGC?a@?_Q?@H?_OCHA?OQA_P?_G_O?WA?_IG_Q?HC@A@ADC'+\
            ...   'A?AI?AC_?QAWOHA?cAGG_I?S?G_OG@GA?`[D?O_IA?`GGCS?OA_?c@?Q?^OAC@_G_Ca@CA@?OGCO'+\
            ...   'H@G@A@?GQC?_Q@GP?_OG?IGGB?OCGaG?cO@A__QGC?E?A@CH@G?GRAGOC_@GGOW@O?O_OGa?_c?G'+\
            ...   'V@CGA_OOaC?a_?a?A_CcC@?CNgA?oC@GGE@?_OH?a@?_?QA`A@?QC?_KGGO_OGCAa@?A?_KCGPC@'+\
            ...   'G_AOAGPGC?D@?a_A?@GGO`KH?Q?C_QGAA_?gOG_OA?_GG`AwH?SA?`?cAI?A@D?I?@?QA?`By?K@'+\
            ...   '?O`GGACA@CGCA@CC_?WO`?`A?OCH?`OCA@COG?I?oC@ACGPCG_AO@_aAA?Aa?g?GD@G?CO`AWOc?'+\
            ...   'HA?OcG_?g@OGCAAAOC@ACJ_`OGACAGCS?CAGI?A`@?OCACG^'
            sage: g2 = '~?A[??osR?WARSETCJ_QWASehOXQg`QwChK?qSeFQ_sTIaWIV?XIR?KAC?B?`?COCG?o?O_'+\
            ...   '@_?`??B?`?o@_O_WCOCHC@_?`W?E?AD_O?WCCeO?WCSEGAGAIaA@_?aw?OK?ER?`?@_HQXA?B@Q_'+\
            ...   'pA?a@Qg_`?o?h[?GOK@IR?@A?BEQcoCG?K\IB?GOCWiTC?GOKWIV??CGEKdH_H_?CB?`?DC??_WC'+\
            ...   'G?SO?AP?O_?g_?D_?`?C__?D_?`?CCo??@_O_XDC???WCGEGg_??a?`G_aa??E?AD_@cC??K?CJ?'+\
            ...   '@@K?O?WCCe?aa?G?KAIB?Gg_A?a?ag_@DC?OK?CV??EOO@?o?XK??GH`A?B?Qco?Gg`A?B@Q_o?C'+\
            ...   'SO`?P?hSO?@DCGOK?IV???K_`A@_HQWC??_cCG?KXIRG?@D?GO?WySEG?@D?GOCWiTCC??a_CGEK'+\
            ...   'DJ_@??K_@A@bHQWAW?@@K??_WCG?g_?CSO?A@_O_@P??Gg_?Ca?`?@P??Gg_?D_?`?C__?EOO?Ao'+\
            ...   '?O_AAW?@@K???WCGEPP??Gg_??B?`?pDC??aa??AGACaAIG?@DC??K?CJ?BGG?@cC??K?CJ?@@K?'+\
            ...   '?_e?G?KAAR?PP??Gg_A?B?a_oAIG?@DC?OCOCTC?Gg_?CSO@?o?P[??X@??K__A@_?qW??OR??GH'+\
            ...   '`A?B?Qco?Gg_?CSO`?@_hOW?AIG?@DCGOCOITC??PP??Gg`A@_@Qw??@cC??qACGE?dH_O?AAW?@'+\
            ...   '@GGO?WqSeO?AIG?@D?GO?WySEG?@DC??a_CGAKTIaA??PP??Gg@A@b@Qw?O?BGG?@c?GOKXIR?KA'+\
            ...   'C?H_?CCo?A@_O_?WCG@P??Gg_?CB?`?COCG@P??Gg_?Ca?`?E?AC?g_?CSO?Ao?O_@_?`@GG?@cC'+\
            ...   '??k?CG??WCGOR??GH_??B?`?o@_O`DC??aa???KACB?a?`AIG?@DC??COCHC@_?`AIG?@DC??K?C'+\
            ...   'J??o?O`cC??qA??E?AD_O?WC?OR??GH_A?B?_cq?B?_AIG?@DC?O?WCSEGAGA?Gg_?CSO@?P?PSO'+\
            ...   'OK?C?PP??Gg_A@_?aw?OK?C?X@??K__A@_?qWCG?K??GH_?CCo`?@_HQXA?B??AIG?@DCGO?WISE'+\
            ...   'GOCO??PP??Gg`A?a@Qg_`?o??@DC??aaCGE?DJ_@A@_??BGG?@cCGOK@IR?@A?BO?AAW?@@GGO?W'+\
            ...   'qSe?`?@g?@DC??a_CG?K\IB?GOCQ??PP??Gg@A?bDQg_@A@_O?AIG?@D?GOKWIV??CGE@??K__?E'+\
            ...   'O?`?pchK?_SA_OI@OGD?gCA_SA@OI?c@H?Q?c_H?QOC_HGAOCc?QOC_HGAOCc@GAQ?c@H?QD?gCA'+\
            ...   '_SA@OI@?gD?_SA_OKA_SA@OI@?gD?_SA_OI@OHI?c_H?QOC_HGAOCc@GAQ?eC_H?QOC_HGAOCc@G'+\
            ...   'AQ?c@XD?_SA_OI@OGD?gCA_SA@PKGO`A@ACGSGO`?`ACICGO_?ACGOcGO`?O`AC`ACHACGO???^?'+\
            ...   '????}Bw????Fo^???????Fo?}?????Bw?^?Bw?????GO`AO`AC`ACGACGOcGO`??aCGO_O`ADACG'+\
            ...   'OGO`A@ACGOA???@{?N_@{?????Fo?}????OFo????N_}????@{????Bw?OACGOgO`A@ACGSGO`?`'+\
            ...   'ACG?OaCGO_GO`AO`AC`ACGACGO_@G???Fo^?????}Bw????Fo??AC@{?????Fo?}?Fo?????^??A'+\
            ...   'OGO`AO`AC@ACGQCGO_GO`A?HAACGOgO`A@ACGOGO`A`ACG?GQ??^?Bw?????N_@{?????Fo?QC??'+\
            ...   'Fo^?????}????@{Fo???CHACGO_O`ACACGOgO`A@ACGO@AOcGO`?O`AC`ACGACGOcGO`?@GQFo??'+\
            ...   '??N_????^@{????Bw??`GRw?????N_@{?????Fo?}???HAO_OI@OGD?gCA_SA@OI@?gDK_??C@GA'+\
            ...   'Q?c@H?Q?c_H?QOC_HEW????????????????????????~~~~~'
            sage: G1 = Graph(g1)
            sage: G2 = Graph(g2)
            sage: G1.is_isomorphic(G2)
            True

        Ensure that isomorphic looped graphs with non-range vertex labels report
        correctly (:trac:`10814`, fixed by :trac:`8395`)::

            sage: G1 = Graph([(0,1), (1,1)])
            sage: G2 = Graph([(0,2), (2,2)])
            sage: G1.is_isomorphic(G2)
            True
            sage: G = Graph(multiedges = True, loops = True)
            sage: H = Graph(multiedges = True, loops = True)
            sage: G.add_edges([(0,1,0),(1,0,1),(1,1,2),(0,0,3)])
            sage: H.add_edges([(0,1,3),(1,0,2),(1,1,1),(0,0,0)])
            sage: G.is_isomorphic(H, certify=True)
            (True, {0: 0, 1: 1})
            sage: set_random_seed(0)
            sage: D = digraphs.RandomDirectedGNP(6, .2)
            sage: D.is_isomorphic(D, certify = True)
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5})
            sage: D.is_isomorphic(D,edge_labels=True, certify = True)
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5})

        Ensure that trac :trac:`11620` is fixed::

            sage: G1 = DiGraph([(0, 0, 'c'), (0, 4, 'b'), (0, 5, 'c'),
            ...   (0, 5, 't'), (1, 1, 'c'), (1, 3,'c'), (1, 3, 't'), (1, 5, 'b'),
            ...   (2, 2, 'c'), (2, 3, 'b'), (2, 4, 'c'),(2, 4, 't'), (3, 1, 't'),
            ...   (3, 2, 'b'), (3, 2, 'c'), (3, 4, 'c'), (4, 0,'b'), (4, 0, 'c'),
            ...   (4, 2, 't'), (4, 5, 'c'), (5, 0, 't'), (5, 1, 'b'), (5, 1, 'c'),
            ...   (5, 3, 'c')], loops=True, multiedges=True)
            sage: G2 = G1.relabel({0:4, 1:5, 2:3, 3:2, 4:1,5:0}, inplace=False)
            sage: G1.canonical_label(edge_labels=True) == G2.canonical_label(edge_labels=True)
            True
            sage: G1.is_isomorphic(G2,edge_labels=True)
            True

        Ensure that :trac:`13114` is fixed ::

            sage: g = Graph([(0, 0, 0), (0, 2, 0), (1, 1, 0), (1, 2, 0), (1, 2, 1), (2, 2, 0)])
            sage: gg = Graph([(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 2, 0), (2, 2, 0), (2, 2, 1)])
            sage: g.is_isomorphic(gg)
            False

        Ensure that trac:`14777` is fixed ::

            sage: g = Graph()
            sage: h = Graph()
            sage: g.is_isomorphic(h)
            True
        """
        from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic
        possible = True
        if self.num_verts() == 0 and other.num_verts() == 0:
            return True
        if self._directed != other._directed:
            possible = False
        if self.order() != other.order():
            possible = False
        if self.size() != other.size():
            possible = False
        if not possible and certify:
            return False, None
        elif not possible:
            return False
        self_vertices = self.vertices()
        other_vertices = other.vertices()
        if edge_labels or self.has_multiple_edges():
            if edge_labels and sorted(self.edge_labels()) != sorted(other.edge_labels()):
                return (False, None) if certify else False
            else:
                G, partition, relabeling, G_edge_labels = graph_isom_equivalent_non_edge_labeled_graph(self, return_relabeling=True, ignore_edge_labels=(not edge_labels), return_edge_labels=True)
                self_vertices = sum(partition,[])
                G2, partition2, relabeling2, G2_edge_labels = graph_isom_equivalent_non_edge_labeled_graph(other, return_relabeling=True, ignore_edge_labels=(not edge_labels), return_edge_labels=True)
                if map(len, partition) != map(len, partition2):
                    return (False, None) if certify else False
                multilabel = (lambda e:e) if edge_labels else (lambda e:map(lambda el: [None, el[1]], e))
                if map(multilabel, G_edge_labels) != map(multilabel, G2_edge_labels):
                    return (False, None) if certify else False
                partition2 = sum(partition2,[])
                other_vertices = partition2
        else:
            G = self; partition = [self_vertices]
            G2 = other; partition2 = other_vertices
        G_to = {}
        for i in xrange(len(self_vertices)):
            G_to[self_vertices[i]] = i
        from sage.graphs.all import Graph, DiGraph
        DoDG = DiGraph if self._directed else Graph
        H = DoDG(len(self_vertices), implementation='c_graph', loops=G.allows_loops())
        HB = H._backend
        for u,v in G.edge_iterator(labels=False):
            u = G_to[u]; v = G_to[v]
            HB.add_edge(u,v,None,G._directed)
        G = HB._cg
        partition = [[G_to[v] for v in cell] for cell in partition]
        GC = G
        G2_to = {}
        for i in xrange(len(other_vertices)):
            G2_to[other_vertices[i]] = i
        H2 = DoDG(len(other_vertices), implementation='c_graph', loops=G2.allows_loops())
        H2B = H2._backend
        for u,v in G2.edge_iterator(labels=False):
            u = G2_to[u]; v = G2_to[v]
            H2B.add_edge(u,v,None,G2._directed)
        G2 = H2B._cg
        partition2 = [G2_to[v] for v in partition2]
        GC2 = G2
        isom = isomorphic(GC, GC2, partition, partition2, (self._directed or self.has_loops()), 1)

        if not isom and certify:
            return False, None
        elif not isom:
            return False
        elif not certify:
            return True
        else:
            isom_trans = {}
            if edge_labels:
                relabeling2_inv = {}
                for x in relabeling2:
                    relabeling2_inv[relabeling2[x]] = x
                for v in self.vertices():
                    isom_trans[v] = relabeling2_inv[other_vertices[isom[G_to[relabeling[v]]]]]
            else:
                for v in self.vertices():
                    isom_trans[v] = other_vertices[isom[G_to[v]]]
            return True, isom_trans

    def canonical_label(self, partition=None, certify=False, verbosity=0, edge_labels=False):
        """
        Returns the unique graph on `\{0,1,...,n-1\}` ( ``n = self.order()`` ) which

        - is isomorphic to self,

        - is invariant in the isomorphism class.

        In other words, given two graphs ``G`` and ``H`` which are isomorphic,
        suppose ``G_c`` and ``H_c`` are the graphs returned by
        ``canonical_label``. Then the following hold:

        - ``G_c == H_c``

        - ``G_c.adjacency_matrix() == H_c.adjacency_matrix()``

        - ``G_c.graph6_string() == H_c.graph6_string()``

        INPUT:

        -  ``partition`` - if given, the canonical label with
           respect to this set partition will be computed. The default is the unit
           set partition.

        -  ``certify`` - if True, a dictionary mapping from the
           (di)graph to its canonical label will be given.

        -  ``verbosity`` - gets passed to nice: prints helpful
           output.

        -  ``edge_labels`` - default False, otherwise allows
           only permutations respecting edge labels.


        EXAMPLES::

            sage: D = graphs.DodecahedralGraph()
            sage: E = D.canonical_label(); E
            Dodecahedron: Graph on 20 vertices
            sage: D.canonical_label(certify=True)
            (Dodecahedron: Graph on 20 vertices, {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18})
            sage: D.is_isomorphic(E)
            True

        Multigraphs::

            sage: G = Graph(multiedges=True,sparse=True)
            sage: G.add_edge((0,1))
            sage: G.add_edge((0,1))
            sage: G.add_edge((0,1))
            sage: G.canonical_label()
            Multi-graph on 2 vertices
            sage: Graph('A?', implementation='c_graph').canonical_label()
            Graph on 2 vertices

        Digraphs::

            sage: P = graphs.PetersenGraph()
            sage: DP = P.to_directed()
            sage: DP.canonical_label().adjacency_matrix()
            [0 0 0 0 0 0 0 1 1 1]
            [0 0 0 0 1 0 1 0 0 1]
            [0 0 0 1 0 0 1 0 1 0]
            [0 0 1 0 0 1 0 0 0 1]
            [0 1 0 0 0 1 0 0 1 0]
            [0 0 0 1 1 0 0 1 0 0]
            [0 1 1 0 0 0 0 1 0 0]
            [1 0 0 0 0 1 1 0 0 0]
            [1 0 1 0 1 0 0 0 0 0]
            [1 1 0 1 0 0 0 0 0 0]

        Edge labeled graphs::

            sage: G = Graph(sparse=True)
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: G.canonical_label(edge_labels=True)
            Graph on 5 vertices
            sage: G.canonical_label(edge_labels=True,certify=True)
            (Graph on 5 vertices, {0: 4, 1: 3, 2: 0, 3: 1, 4: 2})
        """
        from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        from copy import copy

        dig = (self.has_loops() or self._directed)
        if partition is None:
            partition = [self.vertices()]
        if edge_labels or self.has_multiple_edges():
            G, partition, relabeling = graph_isom_equivalent_non_edge_labeled_graph(self, partition, return_relabeling=True)
            G_vertices = sum(partition, [])
            G_to = {}
            for i in xrange(len(G_vertices)):
                G_to[G_vertices[i]] = i
            from sage.graphs.all import Graph, DiGraph
            DoDG = DiGraph if self._directed else Graph
            H = DoDG(len(G_vertices), implementation='c_graph', loops=G.allows_loops())
            HB = H._backend
            for u,v in G.edge_iterator(labels=False):
                u = G_to[u]; v = G_to[v]
                HB.add_edge(u,v,None,G._directed)
            GC = HB._cg
            partition = [[G_to[v] for v in cell] for cell in partition]
            a,b,c = search_tree(GC, partition, certify=True, dig=dig, verbosity=verbosity)
            # c is a permutation to the canonical label of G, which depends only on isomorphism class of self.
            H = copy(self)
            c_new = {}
            for v in self.vertices():
                c_new[v] = c[G_to[relabeling[v]]]
            H.relabel(c_new)
            if certify:
                return H, c_new
            else:
                return H
        G_vertices = sum(partition, [])
        G_to = {}
        for i in xrange(len(G_vertices)):
            G_to[G_vertices[i]] = i
        from sage.graphs.all import Graph, DiGraph
        DoDG = DiGraph if self._directed else Graph
        H = DoDG(len(G_vertices), implementation='c_graph', loops=self.allows_loops())
        HB = H._backend
        for u,v in self.edge_iterator(labels=False):
            u = G_to[u]; v = G_to[v]
            HB.add_edge(u,v,None,self._directed)
        GC = HB._cg
        partition = [[G_to[v] for v in cell] for cell in partition]
        a,b,c = search_tree(GC, partition, certify=True, dig=dig, verbosity=verbosity)
        H = copy(self)
        c_new = {}
        for v in G_to:
            c_new[v] = c[G_to[v]]
        H.relabel(c_new)
        if certify:
            return H, c_new
        else:
            return H

import types

import sage.graphs.distances_all_pairs
GenericGraph.distances_distribution = types.MethodType(sage.graphs.distances_all_pairs.distances_distribution, None, GenericGraph)
GenericGraph.wiener_index = types.MethodType(sage.graphs.distances_all_pairs.wiener_index, None, GenericGraph)

# From Python modules
import sage.graphs.line_graph
GenericGraph.line_graph = sage.graphs.line_graph.line_graph

def tachyon_vertex_plot(g, bgcolor=(1,1,1),
                        vertex_colors=None,
                        vertex_size=0.06,
                        pos3d=None,
                        **kwds):
    """
    Helper function for plotting graphs in 3d with Tachyon. Returns a
    plot containing only the vertices, as well as the 3d position
    dictionary used for the plot.

    INPUT:
     - `pos3d` - a 3D layout of the vertices
     - various rendering options

    EXAMPLES::

        sage: G = graphs.TetrahedralGraph()
        sage: from sage.graphs.generic_graph import tachyon_vertex_plot
        sage: T,p = tachyon_vertex_plot(G, pos3d = G.layout(dim=3))
        sage: type(T)
        <class 'sage.plot.plot3d.tachyon.Tachyon'>
        sage: type(p)
        <type 'dict'>
    """
    assert pos3d is not None
    from math import sqrt
    from sage.plot.plot3d.tachyon import Tachyon

    c = [0,0,0]
    r = []
    verts = g.vertices()

    if vertex_colors is None:
        vertex_colors = { (1,0,0) : verts }
    try:
        for v in verts:
            c[0] += pos3d[v][0]
            c[1] += pos3d[v][1]
            c[2] += pos3d[v][2]
    except KeyError:
        raise KeyError("Oops! You haven't specified positions for all the vertices.")

    order = g.order()
    c[0] = c[0]/order
    c[1] = c[1]/order
    c[2] = c[2]/order
    for v in verts:
        pos3d[v][0] = pos3d[v][0] - c[0]
        pos3d[v][1] = pos3d[v][1] - c[1]
        pos3d[v][2] = pos3d[v][2] - c[2]
        r.append(abs(sqrt((pos3d[v][0])**2 + (pos3d[v][1])**2 + (pos3d[v][2])**2)))
    r = max(r)
    if r == 0:
        r = 1
    for v in verts:
        pos3d[v][0] = pos3d[v][0]/r
        pos3d[v][1] = pos3d[v][1]/r
        pos3d[v][2] = pos3d[v][2]/r
    TT = Tachyon(camera_center=(1.4,1.4,1.4), antialiasing=13, **kwds)
    TT.light((4,3,2), 0.02, (1,1,1))
    TT.texture('bg', ambient=1, diffuse=1, specular=0, opacity=1.0, color=bgcolor)
    TT.plane((-1.6,-1.6,-1.6), (1.6,1.6,1.6), 'bg')

    i = 0
    for color in vertex_colors:
        i += 1
        TT.texture('node_color_%d'%i, ambient=0.1, diffuse=0.9,
                   specular=0.03, opacity=1.0, color=color)
        for v in vertex_colors[color]:
            TT.sphere((pos3d[v][0],pos3d[v][1],pos3d[v][2]), vertex_size, 'node_color_%d'%i)

    return TT, pos3d

def graph_isom_equivalent_non_edge_labeled_graph(g, partition=None, standard_label=None, return_relabeling=False, return_edge_labels=False, inplace=False, ignore_edge_labels=False):
    """
    Helper function for canonical labeling of edge labeled (di)graphs.

    Translates to a bipartite incidence-structure type graph
    appropriate for computing canonical labels of edge labeled and/or multi-edge graphs.
    Note that this is actually computationally equivalent to
    implementing a change on an inner loop of the main algorithm-
    namely making the refinement procedure sort for each label.

    If the graph is a multigraph, it is translated to a non-multigraph,
    where each edge is labeled with a dictionary describing how many
    edges of each label were originally there. Then in either case we
    are working on a graph without multiple edges. At this point, we
    create another (bipartite) graph, whose left vertices are the
    original vertices of the graph, and whose right vertices represent
    the edges. We partition the left vertices as they were originally,
    and the right vertices by common labels: only automorphisms taking
    edges to like-labeled edges are allowed, and this additional
    partition information enforces this on the bipartite graph.

    INPUT:

    - ``g`` -- Graph or DiGraph
    - ``partition`` -- (default:None) if given, the partition of the vertices is as well relabeled
    - ``standard_label`` -- (default:None) the standard label is not considered to be changed
    - ``return_relabeling`` -- (defaut:False) if True, a dictionary containing the relabeling is returned
    - ``return_edge_labels`` -- (defaut:False) if True, the different edge_labels are returned (useful if inplace is True)
    - ``inplace`` -- (default:False) if True, g is modified, otherwise the result is returned. Note that attributes of g are *not* copied for speed issues, only edges and vertices.

    OUTPUT:

    - if not inplace: the unlabeled graph without multiple edges
    - the partition of the vertices
    - if return_relabeling: a dictionary containing the relabeling
    - if return_edge_labels: the list of (former) edge labels is returned

    EXAMPLES::

        sage: from sage.graphs.generic_graph import graph_isom_equivalent_non_edge_labeled_graph

        sage: G = Graph(multiedges=True,sparse=True)
        sage: G.add_edges( (0,1,i) for i in range(10) )
        sage: G.add_edge(1,2,'string')
        sage: G.add_edge(2,123)
        sage: g = graph_isom_equivalent_non_edge_labeled_graph(G, partition=[[0,123],[1,2]]); g
        [Graph on 6 vertices, [[0, 3], [1, 2], [4], [5]]]

        sage: g = graph_isom_equivalent_non_edge_labeled_graph(G); g
        [Graph on 6 vertices, [[0, 1, 2, 3], [4], [5]]]
        sage: g[0].edges()
        [(0, 4, None), (1, 4, None), (1, 5, None), (2, 3, None), (2, 5, None)]

        sage: g = graph_isom_equivalent_non_edge_labeled_graph(G,standard_label='string',return_edge_labels=True); g
        [Graph on 6 vertices, [[0, 1, 2, 3], [5], [4]], [[[None, 1]], [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1], [8, 1], [9, 1]], [['string', 1]]]]
        sage: g[0].edges()
        [(0, 4, None), (1, 2, None), (1, 4, None), (2, 5, None), (3, 5, None)]

        sage: graph_isom_equivalent_non_edge_labeled_graph(G,inplace=True)
        [[[0, 1, 2, 3], [4], [5]]]
        sage: G.edges()
        [(0, 4, None), (1, 4, None), (1, 5, None), (2, 3, None), (2, 5, None)]

    Ensure that #14108 is fixed::

        sage: G=DiGraph([[0,0],[0,0],[0,0],[1,1],[1,1],[1,1]])
        sage: H=DiGraph([[0,0],[0,0],[0,0],[0,0],[1,1],[1,1]])
        sage: G.is_isomorphic(H)
        False
        sage: H=DiGraph([[0,0],[0,0],[0,0],[0,0],[0,0],[1,1],[1,1]])
        sage: HH=DiGraph([[0,0],[0,0],[0,0],[0,0],[1,1],[1,1],[1,1]])
        sage: H.is_isomorphic(HH)
        False
        sage: H.is_isomorphic(HH, edge_labels=True)
        False

    """
    from copy import copy
    from sage.graphs.all import Graph, DiGraph

    g_has_multiple_edges = g.has_multiple_edges()

    if g_has_multiple_edges:
        if g._directed:
            G = DiGraph(loops=g.allows_loops(),sparse=True)
            edge_iter = g._backend.iterator_in_edges(g,True)
        else:
            G = Graph(loops=g.allows_loops(),sparse=True)
            edge_iter = g._backend.iterator_edges(g,True)
        for u,v,l in edge_iter:
            if ignore_edge_labels:
                l = None
            if not G.has_edge(u,v):
                G.add_edge(u,v,[[l,1]])
            else:
                label_list = copy( G.edge_label(u,v) )
                seen_label = False
                for i in xrange(len(label_list)):
                    if label_list[i][0] == l:
                        label_list[i][1] += 1
                        G.set_edge_label(u,v,label_list)
                        seen_label = True
                        break
                if not seen_label:
                    label_list.append([l,1])
                    label_list.sort()
                    G.set_edge_label(u,v,label_list)
        if G.order() < g.order():
            G.add_vertices(g)
        if inplace:
            g._backend = G._backend
    elif not inplace:
        G = copy( g )
    else:
        G = g

    G_order = G.order()
    V = range(G_order)
    if G.vertices() != V:
        relabel_dict = G.relabel(return_map=True)
    else:
        relabel_dict = dict( (i,i) for i in xrange(G_order) )
    if partition is None:
        partition = [V]
    else:
        partition = [ [ relabel_dict[i] for i in part ] for part in partition ]

    if G._directed:
        edge_iter = G._backend.iterator_in_edges(G,True)
    else:
        edge_iter = G._backend.iterator_edges(G,True)

    edges = [ edge for edge in edge_iter ]
    edge_labels = sorted([ label for v1,v2,label in edges if not label == standard_label])
    i = 1

    # edge_labels is sorted. We now remove values which are not unique
    while i < len(edge_labels):
        if edge_labels[i] == edge_labels[i-1]:
            edge_labels.pop(i)
        else:
            i += 1
    i = G_order
    edge_partition = [(el,[]) for el in edge_labels]

    if g_has_multiple_edges: standard_label = [[standard_label,1]]

    for u,v,l in edges:
        if not l == standard_label:
            for el, part in edge_partition:
                if el == l:
                    part.append(i)
                    break

            G._backend.add_edge(u,i,None,True)
            G._backend.add_edge(i,v,None,True)
            G.delete_edge(u,v)
            i += 1
        elif standard_label is not None:
            G._backend.set_edge_label(u,v,None,True)

    # Should we pay attention to edge labels ?
    if ignore_edge_labels:

        # If there are no multiple edges, we can just say that all edges are
        # equivalent to each other without any further consideration.
        if not g_has_multiple_edges:
            edge_partition = [el[1] for el in sorted(edge_partition)]
            edge_partition = [sum(edge_partition,[])]

        # An edge between u and v with label l and multiplicity k being encoded
        # as an uv edge with label [l,k], we must not assume that an edge with
        # multiplicity 2 is equivalent to a simple edge !
        # Hence, we still distinguish edges with different multiplicity
        if g_has_multiple_edges:

            # Compute the multiplicity the label
            multiplicity = lambda x : sum(map(lambda y:y[1],x))

            # Sort the edge according to their multiplicity
            edge_partition = sorted([[multiplicity(el),part] for el, part in sorted(edge_partition)])

            # Gather together the edges with same multiplicity
            i = 1
            while i < len(edge_partition):
                if edge_partition[i][0] == edge_partition[i-1][0]:
                    edge_partition[i-1][1].extend(edge_partition[i][1])
                    edge_partition.pop(i)
                else:
                    i += 1

            # now edge_partition has shape [[multiplicity, list_of_edges],
            # [multiplicity, liste of edges], ...], and we can flatted it to
            # [list of edges, list of edges, ...]
            edge_partition = [el[1] for el in sorted(edge_partition)]

            # Now the edges are partitionned according to the multiplicity they
            # represent, and edge labels are forgotten.

    else:
        edge_partition = [el[1] for el in sorted(edge_partition)]

    new_partition = [ part for part in partition + edge_partition if not part == [] ]

    return_data = []
    if not inplace:
        return_data.append( G )
    return_data.append( new_partition )
    if return_relabeling:
        return_data.append( relabel_dict )
    if return_edge_labels:
        return_data.append( edge_labels )
    return return_data
