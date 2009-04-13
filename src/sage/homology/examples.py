"""
Examples of simplicial complexes

AUTHORS:

- John H. Palmieri (2009-04)

This file constructs some examples of simplicial complexes.  There are
two main types: surfaces and examples related to graph theory.

For surfaces, this file introduces a class, SimplicialSurface, which
is based on SimplicialComplex but also includes a ``connected_sum``
method.  It also has routines for defining the 2-sphere, the n-sphere
for any n, the torus, the real projective plane, the Klein bottle, and
surfaces of arbitrary genus, all as simplicial complexes.

Aside from surfaces, this file also provides some functions for
constructing some other simplicial complexes: the simplicial complex
of not-i-connected graphs on n vertices, the matching complex on n
vertices, and the chessboard complex for an n by i chessboard.  These
provide examples of large simplicial complexes; for example,
not_i_connected_graphs(7,2) has over a million simplices.

All of these examples are accessible by typing
"simplicial_complexes.NAME", where "NAME" is the name of the example;
you can get a list by typing "simplicial_complexes." and hitting the
TAB key::

    Sphere
    Simplex
    Torus
    ProjectivePlane
    KleinBottle
    SurfaceOfGenus
    MooreSpace
    NotIConnectedGraphs
    MatchingComplex
    ChessboardComplex

See the documentation for ``simplicial_complexes`` and for each
particular type of example for full details.
"""

from sage.homology.simplicial_complex import SimplicialComplex, Simplex
from sage.sets.set import Set
from sage.misc.functional import is_even

def rename_vertex(n, keep, left = True):
    """
    Rename a vertex: the three vertices from the list 'keep' get
    relabeled 0, 1, 2, in order.  Any other vertex (e.g. 4) gets
    renamed to by prepending an 'L' or an 'R' (thus to either 'L4' or
    'R4'), depending on whether the argument left is True or False.

    INPUT:

    -  ``n`` - a 'vertex': either an integer or a string

    -  ``keep`` - a list of three vertices

    -  ``left`` - boolean (optional, default True)

    This is used by the ``connected_sum`` method for simplicial
    surfaces.

    EXAMPLES::

        sage: from sage.homology.examples import rename_vertex
        sage: rename_vertex(6, [5, 6, 7])
        1
        sage: rename_vertex(3, [5, 6, 7])
        'L3'
        sage: rename_vertex(3, [5, 6, 7], left=False)
        'R3'
    """
    lookup = dict(zip(keep, [0,1,2]))
    try:
        return lookup[n]
    except:
        if left:
            return "L" + str(n)
        else:
            return "R" + str(n)

class SimplicialSurface(SimplicialComplex):
    """
    A simplicial surface is a simplicial complex structure on a
    surface without boundary.  (The only difference between this and a
    SimplicialComplex is that there is a connected sum operation
    defined for these surfaces.)

    This can either take one argument, a simplicial complex, or two
    arguments, the usual arguments to define a simplicial complex.  It
    does almost no error checking to see whether the resulting complex
    actually represents a surface.

    INPUT:

    -  ``complex`` - either a SimplicialComplex or a vertex set
       suitable for defining a simplicial complex.

    - ``faces`` - (optional, default []).  If the argument complex is
       a simplicial complex, then the faces argument is ignored.  If
       complex is a vertex set, then faces is treated as the list of
       maximal faces to define the simplicial complex.

    OUTPUT: a simplicial complex representing a surface.

    EXAMPLES::

        sage: from sage.homology.examples import SimplicialSurface
        sage: sphere = SimplicialSurface(3, [[1,2,3], [0,2,3], [0,1,3], [0,1,2]])
        sage: sphere
        Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
        sage: sphere == loads(dumps(sphere))
        True
        """

    def __init__(self, complex, faces=[]):
        """
        See ``SimplicialSurface`` for full documentation.

        EXAMPLES::

            sage: from sage.homology.examples import SimplicialSurface
            sage: sphere = SimplicialSurface(3, [[1,2,3], [0,2,3], [0,1,3], [0,1,2]])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
        """
        if isinstance(complex, SimplicialComplex):
            S = complex
        else:
            S = SimplicialComplex(complex, faces)
        if (not S.is_pure()) or (S.dimension() != 2):
            raise ValueError, "This does not appear to be a surface."
        else:
            SimplicialComplex.__init__(self, S._vertex_set, S.facets)

    def connected_sum(self, other):
        """
        The connected sum of this simplicial surface with another one.

        The connected sum of two surfaces is defined as follows:
        remove a 2-simplex from each, producing a boundary for each
        surface consisting of three edges, and then glue the two
        surfaces together along this triangular boundary.

        INPUT:

        -  ``other`` - the other simplicial surface.

        OUTPUT: the connected sum (self # other) as a simplicial surface

        Algorithm: a facet is chosen from each surface, and removed.
        The vertices of these two facets are relabeled to (0,1,2).  Of
        the remaining vertices, the ones from the left-hand factor are
        renamed by prepending an "L", and similarly the remaining
        vertices in the right-hand factor are renamed by prepending an
        "R".

        EXAMPLES::

            sage: from sage.homology.examples import SimplicialSurface
            sage: S = SimplicialSurface(3, [[1,2,3], [0,2,3], [0,1,3], [0,1,2]])
            sage: S     # S is the sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: S.connected_sum(S)
            Simplicial complex with vertex set (0, 1, 2, 'L0', 'R0') and 6 facets
            sage: P = SimplicialSurface(5, [[0,1,2], [0,2,3], [0,1,5], [0,4,5], [0,3,4], [1,2,4], [1,3,4], [1,3,5], [2,3,5], [2,4,5]])
            sage: P     # P is the projective plane
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
            sage: P.connected_sum(P)    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets

        The notation '+' may be used for connected sum, also::

            sage: P + P    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets
            sage: (P + P).homology()
            {0: 0, 1: Z x C2, 2: 0}
        """
        # first find a 2-simplex to remove from each surface
        keep_left = self.facets[0]
        keep_right = other.facets[0]
        # construct the set of vertices:
        left = set(self.vertices()).difference(set(keep_left))
        right = set(other.vertices()).difference(set(keep_right))
        vertex_set = ([0,1,2] + ["L" + str(v) for v in left]
                      + ["R" + str(v) for v in right])
        # construct the set of facets:
        left = set(self.facets).difference(set([keep_left]))
        right = set(other.facets).difference(set([keep_right]))
        facet_set = ([[rename_vertex(v, keep=list(keep_left))
                       for v in face] for face in left]
                     + [[rename_vertex(v, keep=list(keep_right), left=False)
                         for v in face] for face in right])
        # return the new surface
        return SimplicialSurface(vertex_set, facet_set)

    __add__ = connected_sum

def matching(A, B):
    """
    List of maximal matchings between the sets A and B: a matching
    is a set of pairs (a,b) in A x B where each a, b appears in at
    most one pair.  A maximal matching is one which is maximal
    with respect to inclusion of subsets of A x B.

    INPUT:

    -  ``A``, ``B`` - list, tuple, or indeed anything which can be
       converted to a set.

    EXAMPLES::

        sage: from sage.homology.examples import matching
        sage: matching([1,2], [3,4])
        [set([(1, 3), (2, 4)]), set([(2, 3), (1, 4)])]
        sage: matching([0,2], [0])
        [set([(0, 0)]), set([(2, 0)])]
    """
    answer = []
    if len(A) == 0 or len(B) == 0:
        return [set([])]
    for v in A:
        for w in B:
            for M in matching(set(A).difference([v]), set(B).difference([w])):
                new = M.union([(v,w)])
                if new not in answer:
                    answer.append(new)
    return answer

class SimplicialComplexExamples():
    """
    Some examples of simplicial complexes.

    Here are the available examples; you can also type
    "simplicial_complexes."  and hit tab to get a list::

        Sphere
        Simplex
        Torus
        ProjectivePlane
        KleinBottle
        SurfaceOfGenus
        MooreSpace
        NotIConnectedGraphs
        MatchingComplex
        ChessboardComplex

    EXAMPLES::

        sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
        sage: S.homology()
        {0: 0, 1: 0, 2: Z}
        sage: simplicial_complexes.SurfaceOfGenus(3)
        Simplicial complex with 15 vertices and 38 facets
        sage: M4 = simplicial_complexes.MooreSpace(4)
        sage: M4.homology()
        {0: 0, 1: C4, 2: 0}
        sage: simplicial_complexes.MatchingComplex(6).homology()
        {0: 0, 1: Z^16, 2: 0}
    """

    def Sphere(self,n):
        """
        A minimal triangulation of the n-dimensional sphere.

        INPUT:

        -  ``n`` - positive integer

        EXAMPLES::

            sage: simplicial_complexes.Sphere(2)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: simplicial_complexes.Sphere(5).homology()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z}
            sage: [simplicial_complexes.Sphere(n).euler_characteristic() for n in range(6)]
            [2, 0, 2, 0, 2, 0]
            sage: [simplicial_complexes.Sphere(n).f_vector() for n in range(6)]
            [[1, 2],
             [1, 3, 3],
             [1, 4, 6, 4],
             [1, 5, 10, 10, 5],
             [1, 6, 15, 20, 15, 6],
             [1, 7, 21, 35, 35, 21, 7]]
        """
        S = Simplex(n+1)
        facets = S.faces()
        S = SimplicialComplex(n+1, facets)
        if n == 2:
            return(SimplicialSurface(S))
        else:
            return S


    def Simplex(self, n):
        """
        An `n`-dimensional simplex, as a simplicial complex.

        INPUT:

        -  ``n`` - a non-negative integer

        OUTPUT: the simplicial complex consisting of the `n`-simplex
        on vertices `(0, 1, ..., n)` and all of its faces.

        EXAMPLES::

            sage: simplicial_complexes.Simplex(3)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2, 3)}
            sage: simplicial_complexes.Simplex(5).euler_characteristic()
            1
        """
        S = Simplex(n)
        return SimplicialComplex(n, list(S))


    def Torus(self):
        """
        A minimal triangulation of the torus.

        EXAMPLES::

            sage: simplicial_complexes.Torus().homology(1)
            Z x Z
        """
        return SimplicialSurface(6, [[0,1,2], [1,2,4], [1,3,4], [1,3,6],
                                     [0,1,5], [1,5,6], [2,3,5], [2,4,5],
                                     [2,3,6], [0,2,6], [0,3,4], [0,3,5],
                                     [4,5,6], [0,4,6]])

    def ProjectivePlane(self):
        """
        A minimal triangulation of the real projective plane.

        EXAMPLES::

            sage: P = simplicial_complexes.ProjectivePlane()
            sage: P.cohomology()
            {0: 0, 1: 0, 2: C2}
            sage: P.cohomology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
            1: Vector space of dimension 1 over Finite Field of size 2,
            2: Vector space of dimension 1 over Finite Field of size 2}
        """
        return SimplicialSurface(5, [[0,1,2], [0,2,3], [0,1,5], [0,4,5],
                                     [0,3,4], [1,2,4], [1,3,4], [1,3,5],
                                     [2,3,5], [2,4,5]])

    def KleinBottle(self):
        """
        A triangulation of the Klein bottle, formed by taking the
        connected sum of a projective plane with itself.  (This is not
        a minimal triangulation.)

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle()
            Simplicial complex with 9 vertices and 18 facets
        """
        P = simplicial_complexes.ProjectivePlane()
        return P.connected_sum(P)

    def SurfaceOfGenus(self, g, orientable=True):
        """
        A surface of genus g.

        INPUT:

        -  ``g`` - a non-negative integer.  The desired genus

        -  ``orientable`` - boolean (optional, default True). If True,
           return an orientable surface, and if False, return a
           non-orientable surface.

        In the orientable case, return a sphere if `g` is zero, and
        otherwise return a `g`-fold connected sum of a torus with
        itself.

        In the non-orientable case, raise an error if `g` is zero.  If
        `g` is positive, return a `g`-fold connected sum of a
        projective plane with itself.

        EXAMPLES::

            sage: simplicial_complexes.SurfaceOfGenus(2)
            Simplicial complex with 11 vertices and 26 facets
            sage: simplicial_complexes.SurfaceOfGenus(1, orientable=False)
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
        """
        if g == 0:
            if not orientable:
                raise ValueError, "No non-orientable surface of genus zero."
            else:
                return simplicial_complexes.Sphere(2)
        if orientable:
            T = simplicial_complexes.Torus()
        else:
            T = simplicial_complexes.ProjectivePlane()
        S = T
        for i in range(g-1):
            S = S.connected_sum(T)
        return S

    def MooreSpace(self, q):
        """
        Triangulation of the mod q Moore space.

        INPUT:

        -  ``q`` - integer, at least 2

        This is a simplicial complex with simplices of dimension 0, 1,
        and 2, such that its reduced homology is isomorphic to
        `\\ZZ/q\\ZZ` in dimension 1, zero otherwise.

        If `q=2`, this is the projective plane.  If `q>2`, then
        construct it as follows: start with a triangle with vertices
        1, 2, 3.  We take a `3q`-gon forming a `q`-fold cover of the
        triangle, and we form the resulting complex as an
        identification space of the `3q`-gon.  To triangulate this
        identification space, put `q` vertices `A_0`, ..., `A_{q-1}`,
        in the interior, each of which is connected to 1, 2, 3 (two
        facets each: `[1, 2, A_i]`, `[2, 3, A_i]`).  Put `q` more
        vertices in the interior: `B_0`, ..., `B_{q-1}`, with facets
        `[3, 1, B_i]`, `[3, B_i, A_i]`, `[1, B_i, A_{i+1}]`, `[B_i,
        A_i, A_{i+1}]`.  Then triangulate the interior polygon with
        vertices `A_0`, `A_1`, ..., `A_{q-1}`.

        EXAMPLES::

            sage: simplicial_complexes.MooreSpace(3).homology()
            {0: 0, 1: C3, 2: 0}
            sage: simplicial_complexes.MooreSpace(4).suspension().homology()
            {0: 0, 1: 0, 2: C4, 3: 0}
            sage: simplicial_complexes.MooreSpace(8)
            Simplicial complex with 19 vertices and 54 facets
        """
        if q <= 1:
            raise ValueError, "The mod q Moore space is only defined if q is at least 2"
        if q == 2:
            return projective_plane()
        vertices = [1, 2, 3]
        facets = []
        for i in range(q):
            Ai = "A" + str(i)
            Aiplus = "A" + str((i+1)%q)
            Bi = "B" + str(i)
            vertices.append(Ai)
            vertices.append(Bi)
            facets.append([1, 2, Ai])
            facets.append([2, 3, Ai])
            facets.append([3, 1, Bi])
            facets.append([3, Bi, Ai])
            facets.append([1, Bi, Aiplus])
            facets.append([Bi, Ai, Aiplus])
        for i in range(1, q-1):
            Ai = "A" + str(i)
            Aiplus = "A" + str((i+1)%q)
            facets.append(["A0", Ai, Aiplus])
        return SimplicialComplex(vertices, facets)

    # examples from graph theory:

    def NotIConnectedGraphs(self, n, i):
        """
        The simplicial complex of all graphs on n vertices which are
        not i-connected.

        Fix an integer `n>0` and consider the set of graphs on `n`
        vertices.  View each graph as its set of edges, so it is a
        subset of a set of size `n` choose 2.  A graph is
        `i`-connected if, for any `j<i`, if any `j` vertices are
        removed along with the edges emanating from them, then the
        graph remains connected.  Now fix `i`: it is clear that if `G`
        is not `i`-connected, then the same is true for any graph
        obtained from `G` by deleting edges. Thus the set of all
        graphs which are not `i`-connected, viewed as a set of subsets
        of the `n` choose 2 possible edges, is closed under taking
        subsets, and thus forms a simplicial complex.  This function
        produces that simplicial complex.

        INPUT:

        -  ``n``, ``i`` - non-negative integers with `i` at most `n`

        See Dumas et al. for information on computing its homology by
        computer, and see Babson et al. for theory.  For example,
        Babson et al. show that when `i=2`, the reduced homology of
        this complex is nonzero only in dimension `2n-5`, where it is
        free abelian of rank `(n-2)!`.

        EXAMPLES::

            sage: simplicial_complexes.NotIConnectedGraphs(5,2).f_vector()
            [1, 10, 45, 120, 210, 240, 140, 20]
            sage: simplicial_complexes.NotIConnectedGraphs(5,2).homology(5).ngens()
            6

        REFERENCES:

        - Babson, Bjorner, Linusson, Shareshian, and Welker,
          "Complexes of not i-connected graphs," Topology 38 (1999),
          271-299

        - Dumas, Heckenbach, Saunders, Welker, "Computing simplicial
          homology based on efficient Smith normal form algorithms,"
          in "Algebra, geometry, and software systems" (2003),
          177-206.
        """
        G_list = range(1,n+1)
        G_vertices = Set(G_list)
        E_list = []
        for w in G_list:
            for v in range(1,w):
                E_list.append((v,w))
        E = Set(E_list)
        facets = []
        i_minus_one_sets = list(G_vertices.subsets(size=i-1))
        for A in i_minus_one_sets:
            G_minus_A = G_vertices.difference(A)
            for B in G_minus_A.subsets():
                if len(B) > 0 and len(B) < len(G_minus_A):
                    C = G_minus_A.difference(B)
                    facet = E
                    for v in B:
                        for w in C:
                            bad_edge = (min(v,w), max(v,w))
                            facet = facet.difference(Set([bad_edge]))
                    facets.append(facet)
        return SimplicialComplex(E_list, facets)

    def MatchingComplex(self, n):
        """
        The matching complex of graphs on n vertices.

        Fix an integer `n>0` and consider a set `V` of `n` vertices.
        A 'partial matching' on `V` is a graph formed by edges so that
        each vertex is in at most one edge.  If `G` is a partial
        matching, then so is any graph obtained by deleting edges from
        `G`.  Thus the set of all partial matchings on `n` vertices,
        viewed as a set of subsets of the `n` choose 2 possible edges,
        is closed under taking subsets, and thus forms a simplicial
        complex called the 'matching complex'.  This function produces
        that simplicial complex.

        INPUT:

        -  ``n`` - positive integer.

        See Dumas et al. for information on computing its homology by
        computer, and see Wachs for an expository article about the
        theory.  For example, the homology of these complexes seems to
        have only mod 3 torsion, and this has been proved for the
        bottom non-vanishing homology group for the matching complex
        `M_n`.

        EXAMPLES::

            sage: M = simplicial_complexes.MatchingComplex(7)
            sage: H = M.homology()
            sage: H
            {0: 0, 1: C3, 2: Z^20}
            sage: H[2].ngens()
            20
            sage: simplicial_complexes.MatchingComplex(8).homology(2)  # long time (a few seconds)
            Z^132

        REFERENCES:

        - Dumas, Heckenbach, Saunders, Welker, "Computing simplicial
          homology based on efficient Smith normal form algorithms,"
          in "Algebra, geometry, and software systems" (2003),
          177-206.

        - Wachs, "Topology of Matching, Chessboard and General Bounded
          Degree Graph Complexes" (Algebra Universalis Special Issue
          in Memory of Gian-Carlo Rota, Algebra Universalis, 49 (2003)
          345-385)
        """
        G_vertices = Set(range(1,n+1))
        E_list = []
        for w in G_vertices:
            for v in range(1,w):
                E_list.append((v,w))
        facets = []
        if is_even(n):
            half = int(n/2)
            half_n_sets = list(G_vertices.subsets(size=half))
        else:
            half = int((n-1)/2)
            half_n_sets = list(G_vertices.subsets(size=half))
        for X in half_n_sets:
            Xcomp = G_vertices.difference(X)
            if is_even(n):
                if 1 in X:
                    A = X
                    B = Xcomp
                else:
                    A = Xcomp
                    B = X
                for M in matching(A, B):
                    facet = []
                    for pair in M:
                        facet.append(tuple(sorted(pair)))
                        facets.append(facet)
            else:
                for w in Xcomp:
                    if 1 in X or (w == 1 and 2 in X):
                        A = X
                        B = Xcomp.difference([w])
                    else:
                        B = X
                        A = Xcomp.difference([w])
                    for M in matching(A, B):
                        facet = []
                        for pair in M:
                            facet.append(tuple(sorted(pair)))
                        facets.append(facet)
        return SimplicialComplex(E_list, facets)

    def ChessboardComplex(self, n, i):
        r"""
        The chessboard complex for an n by i chessboard.

        Fix integers `n, i > 0` and consider sets `V` of `n` vertices
        and `W` of `i` vertices.  A 'partial matching' between `V` and
        `W` is a graph formed by edges `(v,w)` with `v \in V` and `w
        \in W` so that each vertex is in at most one edge.  If `G` is
        a partial matching, then so is any graph obtained by deleting
        edges from `G`.  Thus the set of all partial matchings on `V`
        and `W`, viewed as a set of subsets of the `n+i` choose 2
        possible edges, is closed under taking subsets, and thus forms
        a simplicial complex called the 'chessboard complex'.  This
        function produces that simplicial complex.  (It is called the
        chessboard complex because such graphs also correspond to ways
        of placing rooks on an `n` by `i` chessboard so that none of
        them are attacking each other.)

        INPUT:

        -  ``n, i`` - positive integers.

        See Dumas et al. for information on computing its homology by
        computer, and see Wachs for an expository article about the
        theory.

        EXAMPLES::

            sage: C = simplicial_complexes.ChessboardComplex(5,5)
            sage: C.f_vector()
            [1, 25, 200, 600, 600, 120]
            sage: simplicial_complexes.ChessboardComplex(3,3).homology()
            {0: 0, 1: Z x Z x Z x Z, 2: 0}

        REFERENCES:

        - Dumas, Heckenbach, Saunders, Welker, "Computing simplicial
          homology based on efficient Smith normal form algorithms,"
          in "Algebra, geometry, and software systems" (2003),
          177-206.

        - Wachs, "Topology of Matching, Chessboard and General Bounded
          Degree Graph Complexes" (Algebra Universalis Special Issue
          in Memory of Gian-Carlo Rota, Algebra Universalis, 49 (2003)
          345-385)
        """
        A = range(n)
        B = range(i)
        E_dict = {}
        index = 0
        for v in A:
            for w in B:
                E_dict[(v,w)] = index
                index += 1
        E = range(n*i)
        facets = []
        for M in matching(A, B):
            facet = []
            for pair in M:
                facet.append(E_dict[pair])
            facets.append(facet)
        return SimplicialComplex(E, facets)

simplicial_complexes = SimplicialComplexExamples()
