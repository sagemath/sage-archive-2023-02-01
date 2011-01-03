r"""
Common Digraphs

Generators for common digraphs.

AUTHORS:

- Robert L. Miller (2006)
- Emily A. Kirkman (2006)
- Michael C. Yurko (2009)
"""

################################################################################
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################

import graph
from   math import sin, cos, pi
from sage.misc.randstate import current_randstate
from sage.graphs.digraph import DiGraph


class DiGraphGenerators():
    r"""
    A class consisting of constructors for several common digraphs,
    including orderly generation of isomorphism class representatives.

    A list of all graphs and graph structures in this database is
    available via tab completion. Type "digraphs." and then hit tab to
    see which graphs are available.

    The docstrings include educational information about each named
    digraph with the hopes that this class can be used as a reference.

    The constructors currently in this class include::

                Random Directed Graphs:
                    - RandomDirectedGN
                    - RandomDirectedGNC
                    - RandomDirectedGNR

                Families of Graphs:
                    - DeBruijn



    ORDERLY GENERATION: digraphs(vertices, property=lambda x: True,
    augment='edges', size=None)

    Accesses the generator of isomorphism class representatives.
    Iterates over distinct, exhaustive representatives.

    INPUT:


    -  ``vertices`` - natural number

    -  ``property`` - any property to be tested on digraphs
       before generation.

    -  ``augment`` - choices:

    -  ``'vertices'`` - augments by adding a vertex, and
       edges incident to that vertex. In this case, all digraphs on up to
       n=vertices are generated. If for any digraph G satisfying the
       property, every subgraph, obtained from G by deleting one vertex
       and only edges incident to that vertex, satisfies the property,
       then this will generate all digraphs with that property. If this
       does not hold, then all the digraphs generated will satisfy the
       property, but there will be some missing.

    -  ``'edges'`` - augments a fixed number of vertices by
       adding one edge In this case, all digraphs on exactly n=vertices
       are generated. If for any graph G satisfying the property, every
       subgraph, obtained from G by deleting one edge but not the vertices
       incident to that edge, satisfies the property, then this will
       generate all digraphs with that property. If this does not hold,
       then all the digraphs generated will satisfy the property, but
       there will be some missing.

    -  ``implementation`` - which underlying implementation to use (see DiGraph?)

    -  ``sparse`` - ignored if implementation is not ``c_graph``

    EXAMPLES: Print digraphs on 2 or less vertices.

    ::

        sage: for D in digraphs(2, augment='vertices'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Note that we can also get digraphs with underlying Cython implementation::

        sage: for D in digraphs(2, augment='vertices', implementation='c_graph'):
        ...    print D
        ...
        Digraph on 0 vertices
        Digraph on 1 vertex
        Digraph on 2 vertices
        Digraph on 2 vertices
        Digraph on 2 vertices

    Print digraphs on 3 vertices.

    ::

        sage: for D in digraphs(3):
        ...    print D
        Digraph on 3 vertices
        Digraph on 3 vertices
        ...
        Digraph on 3 vertices
        Digraph on 3 vertices

    Generate all digraphs with 4 vertices and 3 edges.

    ::

        sage: L = digraphs(4, size=3)
        sage: len(list(L))
        13

    Generate all digraphs with 4 vertices and up to 3 edges.

    ::

        sage: L = list(digraphs(4, lambda G: G.size() <= 3))
        sage: len(L)
        20
        sage: graphs_list.show_graphs(L)  # long time

    Generate all digraphs with degree at most 2, up to 5 vertices.

    ::

        sage: property = lambda G: ( max([G.degree(v) for v in G] + [0]) <= 2 )
        sage: L = list(digraphs(5, property, augment='vertices'))
        sage: len(L)
        75

    Generate digraphs on the fly: (see http://oeis.org/classic/A000273)

    ::

        sage: for i in range(0, 5):
        ...    print len(list(digraphs(i)))
        1
        1
        3
        16
        218

    REFERENCE:

    - Brendan D. McKay, Isomorph-Free Exhaustive generation.  Journal
      of Algorithms Volume 26, Issue 2, February 1998, pages 306-324.
    """

    def ButterflyGraph(self, n, vertices='strings'):
        """
        Returns a n-dimensional butterfly graph. The vertices consist of
        pairs (v,i), where v is an n-dimensional tuple (vector) with binary
        entries (or a string representation of such) and i is an integer in
        [0..n]. A directed edge goes from (v,i) to (w,i+1) if v and w are
        identical except for possibly v[i] != w[i].

        A butterfly graph has `(2^n)(n+1)` vertices and
        `n2^{n+1}` edges.

        INPUT:


        -  ``vertices`` - 'strings' (default) or 'vectors',
           specifying whether the vertices are zero-one strings or actually
           tuples over GF(2).


        EXAMPLES::

            sage: digraphs.ButterflyGraph(2).edges(labels=False)
            [(('00', 0), ('00', 1)),
            (('00', 0), ('10', 1)),
            (('00', 1), ('00', 2)),
            (('00', 1), ('01', 2)),
            (('01', 0), ('01', 1)),
            (('01', 0), ('11', 1)),
            (('01', 1), ('00', 2)),
            (('01', 1), ('01', 2)),
            (('10', 0), ('00', 1)),
            (('10', 0), ('10', 1)),
            (('10', 1), ('10', 2)),
            (('10', 1), ('11', 2)),
            (('11', 0), ('01', 1)),
            (('11', 0), ('11', 1)),
            (('11', 1), ('10', 2)),
            (('11', 1), ('11', 2))]
            sage: digraphs.ButterflyGraph(2,vertices='vectors').edges(labels=False)
            [(((0, 0), 0), ((0, 0), 1)),
            (((0, 0), 0), ((1, 0), 1)),
            (((0, 0), 1), ((0, 0), 2)),
            (((0, 0), 1), ((0, 1), 2)),
            (((0, 1), 0), ((0, 1), 1)),
            (((0, 1), 0), ((1, 1), 1)),
            (((0, 1), 1), ((0, 0), 2)),
            (((0, 1), 1), ((0, 1), 2)),
            (((1, 0), 0), ((0, 0), 1)),
            (((1, 0), 0), ((1, 0), 1)),
            (((1, 0), 1), ((1, 0), 2)),
            (((1, 0), 1), ((1, 1), 2)),
            (((1, 1), 0), ((0, 1), 1)),
            (((1, 1), 0), ((1, 1), 1)),
            (((1, 1), 1), ((1, 0), 2)),
            (((1, 1), 1), ((1, 1), 2))]
        """
        # We could switch to Sage integers to handle arbitrary n.
        if vertices=='strings':
            if n>=31:
                raise NotImplementedError, "vertices='strings' is only valid for n<=30."
            from sage.graphs.generic_graph_pyx import binary
            butterfly = {}
            for v in xrange(2**n):
                for i in range(n):
                    w = v
                    w ^= (1 << i)   # push 1 to the left by i and xor with w
                    bv = binary(v)
                    bw = binary(w)
                    # pad and reverse the strings
                    padded_bv = ('0'*(n-len(bv))+bv)[::-1]
                    padded_bw = ('0'*(n-len(bw))+bw)[::-1]
                    butterfly[(padded_bv,i)]=[(padded_bv,i+1), (padded_bw,i+1)]
        elif vertices=='vectors':
            from sage.modules.free_module import VectorSpace
            from sage.rings.finite_rings.constructor import FiniteField
            from copy import copy
            butterfly = {}
            for v in VectorSpace(FiniteField(2),n):
                for i in xrange(n):
                    w=copy(v)
                    w[i] += 1 # Flip the ith bit
                    # We must call tuple since vectors are mutable.  To obtain
                    # a vector from the tuple t, just call vector(t).
                    butterfly[(tuple(v),i)]=[(tuple(v),i+1), (tuple(w),i+1)]
        else:
            raise NotImplementedError, "vertices must be 'strings' or 'vectors'."
        return DiGraph(butterfly)

    def Circuit(self,n):
        r"""
        Returns the circuit on `n` vertices

        The circuit is an oriented ``CycleGraph``

        EXAMPLE:

        A circuit is the smallest strongly connected digraph::

            sage: circuit = digraphs.Circuit(15)
            sage: len(circuit.strongly_connected_components()) == 1
            True
        """
        if n<0:
            raise ValueError("The number of vertices must be a positive integer.")

        g = DiGraph()
        g.name("Circuit on "+str(n)+" vertices")

        if n==0:
            return g
        elif n == 1:
            g.allow_loops(True)
            g.add_edge(0,0)
            return g
        else:
            g.add_edges([(i,i+1) for i in xrange(n-1)])
            g.add_edge(n-1,0)
            return g

    def DeBruijn(self,k,n):
        r"""
        Returns the De Bruijn diraph with parameters `k,n`.

        The De Bruijn digraph with parameters `k,n` is built
        upon a set of vertices equal to the set of words of
        length `n` from a dictionary of `k` letters.

        In this digraph, there is an arc `w_1w_2` if `w_2`
        can be obtained from `w_1` by removing the leftmost
        letter and adding a new letter at its right end.
        For more information, see the
        `Wikipedia article on De Bruijn graph
        <http://en.wikipedia.org/wiki/De_Bruijn_graph>`_.

        INPUT:

        - ``k`` -- Two possibilities for this parameter :
              - an integer equal to the cardinality of the
                alphabet to use.
              - An iterable object to be used as the set
                of letters
        - ``n`` -- An integer equal to the length of words in
          the De Bruijn digraph.

        EXAMPLES::

            sage: db=digraphs.DeBruijn(2,2); db
            De Bruijn digraph (k=2, n=2): Looped digraph on 4 vertices
            sage: db.order()
            4
            sage: db.size()
            8

        TESTS::

            sage: digraphs.DeBruijn(5,0)
            De Bruijn digraph (k=5, n=0): Looped multi-digraph on 1 vertex
            sage: digraphs.DeBruijn(0,0)
            De Bruijn digraph (k=0, n=0): Looped multi-digraph on 0 vertices

        """
        from sage.combinat.words.words import Words
        from sage.rings.integer import Integer

        W = Words(range(k) if isinstance(k, Integer) else k, n)
        A = Words(range(k) if isinstance(k, Integer) else k, 1)
        g = DiGraph(loops=True)

        if n == 0 :
            g.allow_multiple_edges(True)
            v = W[0]
            for a in A:
                g.add_edge(v.string_rep(), v.string_rep(), a.string_rep())
        else:
            for w in W:
                ww = w[1:]
                for a in A:
                    g.add_edge(w.string_rep(), (ww*a).string_rep(), a.string_rep())

        g.name( "De Bruijn digraph (k=%s, n=%s)"%(k,n) )
        return g


    def RandomDirectedGN(self, n, kernel=lambda x:x, seed=None):
        """
        Returns a random GN (growing network) digraph with n vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The default attachment kernel is a linear function of
        degree. The digraph is always a tree, so in particular it is a
        directed acyclic graph.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``kernel`` - the attachment kernel

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGN(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (3, 1), (4, 0), (5, 0), (6, 1), (7, 0), (8, 3), (9, 0), (10, 8), (11, 3), (12, 9), (13, 8), (14, 0), (15, 11), (16, 11), (17, 5), (18, 11), (19, 6), (20, 5), (21, 14), (22, 5), (23, 18), (24, 11)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gn_graph(n, kernel, seed=seed))

    def RandomDirectedGNC(self, n, seed=None):
        """
        Returns a random GNC (growing network with copying) digraph with n
        vertices.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen with a
        preferential attachment model, i.e. probability is proportional to
        degree. The new vertex is also linked to all of the previously
        added vertex's successors.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``seed`` - for the random number generator


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNC(25)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Network Growth by
          Copying, Phys. Rev. E vol. 71 (2005), p. 036118.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gnc_graph(n, seed=seed))

    def RandomDirectedGNP(self, n, p):
        r"""
        Returns a random digraph on `n` nodes. Each edge is
        inserted independently with probability `p`.

        REFERENCES:

        - [1] P. Erdos and A. Renyi, On Random Graphs, Publ.  Math. 6,
          290 (1959).

        - [2] E. N. Gilbert, Random Graphs, Ann. Math.  Stat., 30,
          1141 (1959).

        PLOTTING: When plotting, this graph will use the default
        spring-layout algorithm, unless a position dictionary is
        specified.

        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNP(10, .2)
            sage: D.num_verts()
            10
            sage: D.edges(labels=False)
            [(0, 1), (0, 3), (0, 6), (0, 8), (1, 4), (3, 7), (4, 1), (4, 8), (5, 2), (5, 6), (5, 8), (6, 4), (7, 6), (8, 4), (8, 5), (8, 7), (8, 9), (9, 3), (9, 4), (9, 6)]
        """
        from sage.misc.prandom import random
        D = DiGraph(n)
        for i in xrange(n):
            for j in xrange(i):
                if random() < p:
                    D.add_edge(i,j)
            for j in xrange(i+1,n):
                if random() < p:
                    D.add_edge(i,j)
        return D

    def RandomDirectedGNR(self, n, p, seed=None):
        """
        Returns a random GNR (growing network with redirection) digraph
        with n vertices and redirection probability p.

        The digraph is constructed by adding vertices with a link to one
        previously added vertex. The vertex to link to is chosen uniformly.
        With probability p, the arc is instead redirected to the successor
        vertex. The digraph is always a tree.

        INPUT:


        -  ``n`` - number of vertices.

        -  ``p`` - redirection probability

        -  ``seed`` - for the random number generator.


        EXAMPLE::

            sage: D = digraphs.RandomDirectedGNR(25, .2)
            sage: D.edges(labels=False)
            [(1, 0), (2, 0), (2, 1), (3, 0), (4, 0), (4, 1), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (7, 0), (7, 1), (7, 4), (8, 0), (9, 0), (9, 8), (10, 0), (10, 1), (10, 2), (10, 5), (11, 0), (11, 8), (11, 9), (12, 0), (12, 8), (12, 9), (13, 0), (13, 1), (14, 0), (14, 8), (14, 9), (14, 12), (15, 0), (15, 8), (15, 9), (15, 12), (16, 0), (16, 1), (16, 4), (16, 7), (17, 0), (17, 8), (17, 9), (17, 12), (18, 0), (18, 8), (19, 0), (19, 1), (19, 4), (19, 7), (20, 0), (20, 1), (20, 4), (20, 7), (20, 16), (21, 0), (21, 8), (22, 0), (22, 1), (22, 4), (22, 7), (22, 19), (23, 0), (23, 8), (23, 9), (23, 12), (23, 14), (24, 0), (24, 8), (24, 9), (24, 12), (24, 15)]
            sage: D.show()  # long time

        REFERENCE:

        - [1] Krapivsky, P.L. and Redner, S. Organization of Growing
          Random Networks, Phys. Rev. E vol. 63 (2001), p. 066123.
        """
        if seed is None:
            seed = current_randstate().long_seed()
        import networkx
        return DiGraph(networkx.gnc_graph(n, seed=seed))

################################################################################
#   DiGraph Iterators
################################################################################

    def __call__(self, vertices, property=lambda x: True, augment='edges', size=None, implementation='c_graph', sparse=True):
        """
        Accesses the generator of isomorphism class representatives.
        Iterates over distinct, exhaustive representatives.

        INPUT:


        -  ``vertices`` - natural number

        -  ``property`` - any property to be tested on digraphs
           before generation.

        -  ``augment`` - choices:

        -  ``'vertices'`` - augments by adding a vertex, and
           edges incident to that vertex. In this case, all digraphs on up to
           n=vertices are generated. If for any digraph G satisfying the
           property, every subgraph, obtained from G by deleting one vertex
           and only edges incident to that vertex, satisfies the property,
           then this will generate all digraphs with that property. If this
           does not hold, then all the digraphs generated will satisfy the
           property, but there will be some missing.

        -  ``'edges'`` - augments a fixed number of vertices by
           adding one edge In this case, all digraphs on exactly n=vertices
           are generated. If for any graph G satisfying the property, every
           subgraph, obtained from G by deleting one edge but not the vertices
           incident to that edge, satisfies the property, then this will
           generate all digraphs with that property. If this does not hold,
           then all the digraphs generated will satisfy the property, but
           there will be some missing.

        -  ``implementation`` - which underlying implementation to use (see DiGraph?)

        -  ``sparse`` - ignored if implementation is not ``c_graph``

        EXAMPLES: Print digraphs on 2 or less vertices.

        ::

            sage: for D in digraphs(2, augment='vertices'):
            ...    print D
            ...
            Digraph on 0 vertices
            Digraph on 1 vertex
            Digraph on 2 vertices
            Digraph on 2 vertices
            Digraph on 2 vertices

        Print digraphs on 3 vertices.

        ::

            sage: for D in digraphs(3):
            ...    print D
            Digraph on 3 vertices
            Digraph on 3 vertices
            ...
            Digraph on 3 vertices
            Digraph on 3 vertices

        For more examples, see the class level documentation, or type

        ::

            sage: digraphs? # not tested

        REFERENCE:

        - Brendan D. McKay, Isomorph-Free Exhaustive generation.
          Journal of Algorithms Volume 26, Issue 2, February 1998,
          pages 306-324.
        """
        if size is not None:
            extra_property = lambda x: x.size() == size
        else:
            extra_property = lambda x: True
        if augment == 'vertices':
            from sage.graphs.graph_generators import canaug_traverse_vert
            g = DiGraph(implementation=implementation, sparse=sparse)
            for gg in canaug_traverse_vert(g, [], vertices, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        elif augment == 'edges':
            from sage.graphs.graph_generators import canaug_traverse_edge
            g = DiGraph(vertices, implementation=implementation, sparse=sparse)
            gens = []
            for i in range(vertices-1):
                gen = range(i)
                gen.append(i+1); gen.append(i)
                gen += range(i+2, vertices)
                gens.append(gen)
            for gg in canaug_traverse_edge(g, gens, property, dig=True, implementation=implementation, sparse=sparse):
                if extra_property(gg):
                    yield gg
        else:
            raise NotImplementedError()


# Easy access to the graph generators from the command line:
digraphs = DiGraphGenerators()




