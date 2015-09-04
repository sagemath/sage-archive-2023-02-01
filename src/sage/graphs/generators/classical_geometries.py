# -*- coding: utf-8 -*-
r"""
Families of graphs derived from classical geometries over finite fields

These include graphs of polar spaces, affine polar graphs, graphs
related to Hermitean unitals, graphs on nonistropic points, etc

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

###########################################################################
#
#           Copyright (C) 2015 Sagemath project
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################


from copy import copy
from math import sin, cos, pi
from sage.graphs.graph import Graph
from sage.graphs import graph

def SymplecticPolarGraph(d, q, algorithm=None):
    r"""
    Returns the Symplectic polar graph `Sp(d,q)`.

    The Symplectic Polar Graph `Sp(d,q)` is built from a projective space of dimension
    `d-1` over a field `F_q`, and a symplectic form `f`. Two vertices `u,v` are
    made adjacent if `f(u,v)=0`.

    See the page `on symplectic graphs on Andries Brouwer's website
    <http://www.win.tue.nl/~aeb/graphs/Sp.html>`_.

    INPUT:

    - ``d,q`` (integers) -- note that only even values of `d` are accepted by
      the function.

    - ``algorithm`` -- if set to 'gap' then the computation is carried via GAP
      library interface, computing totally singular subspaces, which is faster for `q>3`.
      Otherwise it is done directly.

    EXAMPLES:

    Computation of the spectrum of `Sp(6,2)`::

        sage: g = graphs.SymplecticPolarGraph(6,2)
        sage: g.is_strongly_regular(parameters=True)
        (63, 30, 13, 15)
        sage: set(g.spectrum()) == {-5, 3, 30}
        True

    The parameters of `Sp(4,q)` are the same as of `O(5,q)`, but they are
    not isomorphic if `q` is odd::

        sage: G = graphs.SymplecticiPolarGraph(4,3)
        sage: G.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)
        sage: O=graphs.OrthogonalPolarGraph(5,3)
        sage: O.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)
        sage: O.is_isomorphic(G)
        False
        sage: graphs.SymplecticPolarGraph(6,4,algorithm="gap").is_strongly_regular(parameters=True) # not tested (long time)
        (1365, 340, 83, 85)

    TESTS::

        sage: graphs.SymplecticPolarGraph(4,4,algorithm="gap").is_strongly_regular(parameters=True)
        (85, 20, 3, 5)
        sage: graphs.SymplecticPolarGraph(4,4).is_strongly_regular(parameters=True)
        (85, 20, 3, 5)
        sage: graphs.SymplecticPolarGraph(4,4,algorithm="blah")
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm!
    """
    if d < 1 or d%2 != 0:
        raise ValueError("d must be even and greater than 2")

    if algorithm == "gap":     # faster for larger (q>3)  fields
        from sage.libs.gap.libgap import libgap
        G = _polar_graph(d, q, libgap.SymplecticGroup(d, q))

    elif algorithm == None:    # faster for small (q<4) fields
        from sage.rings.finite_rings.constructor import FiniteField
        from sage.modules.free_module import VectorSpace
        from sage.schemes.projective.projective_space import ProjectiveSpace
        from sage.matrix.constructor import identity_matrix, block_matrix, zero_matrix

        F = FiniteField(q,"x")
        M = block_matrix(F, 2, 2,
                         [zero_matrix(F,d/2),
                          identity_matrix(F,d/2),
                          -identity_matrix(F,d/2),
                          zero_matrix(F,d/2)])

        V = VectorSpace(F,d)
        PV = list(ProjectiveSpace(d-1,F))
        G = Graph([[tuple(_) for _ in PV], lambda x,y:V(x)*(M*V(y)) == 0], loops = False)

    else:
        raise ValueError("unknown algorithm!")

    G.name("Symplectic Polar Graph Sp("+str(d)+","+str(q)+")")
    G.relabel()
    return G

from sage.misc.superseded import deprecated_function_alias
SymplecticGraph = deprecated_function_alias(19136, SymplecticPolarGraph)

def AffineOrthogonalPolarGraph(d,q,sign="+"):
    r"""
    Returns the affine polar graph `VO^+(d,q),VO^-(d,q)` or `VO(d,q)`.

    Affine Polar graphs are built from a `d`-dimensional vector space over
    `F_q`, and a quadratic form which is hyperbolic, elliptic or parabolic
    according to the value of ``sign``.

    Note that `VO^+(d,q),VO^-(d,q)` are strongly regular graphs, while `VO(d,q)`
    is not.

    For more information on Affine Polar graphs, see `Affine Polar
    Graphs page of Andries Brouwer's website
    <http://www.win.tue.nl/~aeb/graphs/VO.html>`_.

    INPUT:

    - ``d`` (integer) -- ``d`` must be even if ``sign != None``, and odd
      otherwise.

    - ``q`` (integer) -- a power of a prime number, as `F_q` must exist.

    - ``sign`` -- must be equal to ``"+"``, ``"-"``, or ``None`` to compute
      (respectively) `VO^+(d,q),VO^-(d,q)` or `VO(d,q)`. By default
      ``sign="+"``.

    .. NOTE::

        The graph `VO^\epsilon(d,q)` is the graph induced by the
        non-neighbors of a vertex in an :meth:`Orthogonal Polar Graph
        <OrthogonalPolarGraph>` `O^\epsilon(d+2,q)`.

    EXAMPLES:

    The :meth:`Brouwer-Haemers graph <BrouwerHaemersGraph>` is isomorphic to
    `VO^-(4,3)`::

        sage: g = graphs.AffineOrthogonalPolarGraph(4,3,"-")
        sage: g.is_isomorphic(graphs.BrouwerHaemersGraph())
        True

    Some examples from `Brouwer's table or strongly regular graphs
    <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`_::

        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"-"); g
        Affine Polar Graph VO^-(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 27, 10, 12)
        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"+"); g
        Affine Polar Graph VO^+(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)

    When ``sign is None``::

        sage: g = graphs.AffineOrthogonalPolarGraph(5,2,None); g
        Affine Polar Graph VO^-(5,2): Graph on 32 vertices
        sage: g.is_strongly_regular(parameters=True)
        False
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True
    """
    if sign in ["+","-"]:
        s = 1 if sign == "+" else -1
        if d%2 == 1:
            raise ValueError("d must be even when sign!=None")
    else:
        if d%2 == 0:
            raise ValueError("d must be odd when sign==None")
        s = 0

    from sage.interfaces.gap import gap
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module import VectorSpace
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap
    from itertools import combinations

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(s,d,q))['matrix'])
    F = libgap.GF(q).sage()
    V = list(VectorSpace(F,d))

    G = Graph()
    G.add_vertices([tuple(_) for _ in V])
    for x,y in combinations(V,2):
        if not (x-y)*M*(x-y):
            G.add_edge(tuple(x),tuple(y))

    G.name("Affine Polar Graph VO^"+str('+' if s == 1 else '-')+"("+str(d)+","+str(q)+")")
    G.relabel()
    return G

def OrthogonalPolarGraph(m, q, sign="+"):
    r"""
    Returns the Orthogonal Polar Graph `O^{\epsilon}(m,q)`.

    For more information on Orthogonal Polar graphs, see see the `page of
    Andries Brouwer's website <http://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m,q`` (integers) -- `q` must be a prime power.

    - ``sign`` -- ``"+"`` or ``"-"`` if `m` is even, ``"+"`` (default)
      otherwise.

    EXAMPLES::

        sage: G = graphs.OrthogonalPolarGraph(6,3,"+"); G
        Orthogonal Polar Graph O^+(6, 3): Graph on 130 vertices
        sage: G.is_strongly_regular(parameters=True)
        (130, 48, 20, 16)
        sage: G = graphs.OrthogonalPolarGraph(6,3,"-"); G
        Orthogonal Polar Graph O^-(6, 3): Graph on 112 vertices
        sage: G.is_strongly_regular(parameters=True)
        (112, 30, 2, 10)
        sage: G = graphs.OrthogonalPolarGraph(5,3); G
        Orthogonal Polar Graph O(5, 3): Graph on 40 vertices
        sage: G.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)

    TESTS::

        sage: G = graphs.OrthogonalPolarGraph(4,3,"")
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '-' or '+' when m is even
        sage: G = graphs.OrthogonalPolarGraph(5,3,"-")
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '' or '+' when m is odd
    """
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module_element import free_module_element as vector
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap
    from itertools import combinations

    if m % 2 == 0:
        if sign != "+" and sign != "-":
            raise ValueError("sign must be equal to either '-' or '+' when "
                             "m is even")
    else:
        if sign != "" and sign != "+":
            raise ValueError("sign must be equal to either '' or '+' when "
                             "m is odd")
        sign = ""

    e = {'+': 1,
         '-': -1,
         '' : 0}[sign]

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(e,m,q))['matrix'])
    Fq = libgap.GF(q).sage()
    PG = ProjectiveSpace(m - 1, Fq)
    m_over_two = m // 2

    def F(x):
        return x*M*x

    V = [x for x in PG if F(vector(x)) == 0]

    G = Graph([V,lambda x,y:F(vector(x)-vector(y))==0],loops=False)

    G.relabel()
    G.name("Orthogonal Polar Graph O" + ("^" + sign if sign else "") + str((m, q)))
    return G

def _polar_graph(m, q, g, intersection_size=None):
    r"""
    The helper function to build graphs `(D)U(m,q)` and `(D)Sp(m,q)`

    Building a graph on an orbit of a group `g` of `m\times m` matrices over `GF(q)` on
    the points (or subspaces of dimension ``m//2``) isotropic w.r.t. the form `F`
    left invariant by the group `g`.

    The only constraint is that the first ``m//2`` elements of the standard
    basis must generate a totally isotropic w.r.t. `F` subspace; this is the case with
    these groups coming from GAP; namely, `F` has the anti-diagonal all-1 matrix.

    INPUT:

    - ``m`` -- the dimension of the underlying vector space

    - ``q`` -- the size of the field

    - ``g`` -- the group acting

    - ``intersection_size`` -- if ``None``, build the graph on the isotropic points, with
      adjacency being orthogonality w.r.t. `F`. Otherwise, build the graph on the maximal
      totally isotropic subspaces, with adjacency specified by ``intersection_size`` being
      as given.

    TESTS::

        sage: from sage.graphs.generators.classical_geometries import _polar_graph
        sage: _polar_graph(4, 4, libgap.GeneralUnitaryGroup(4, 2))
        Graph on 45 vertices
        sage: _polar_graph(4, 4, libgap.GeneralUnitaryGroup(4, 2), intersection_size=1)
        Graph on 27 vertices
    """
    from sage.libs.gap.libgap import libgap
    from itertools import combinations
    W=libgap.FullRowSpace(libgap.GF(q), m)  # F_q^m
    B=libgap.Elements(libgap.Basis(W))      # the standard basis of W
    V = libgap.Orbit(g,B[0],libgap.OnLines) # orbit on isotropic points
    gp = libgap.Action(g,V,libgap.OnLines)  # make a permutation group
    s = libgap.Subspace(W,[B[i] for i in range(m//2)]) # a totally isotropic subspace
    # and the points there
    sp = [libgap.Elements(libgap.Basis(x))[0] for x in libgap.Elements(s.Subspaces(1))]
    h = libgap.Set(map(lambda x: libgap.Position(V, x), sp)) # indices of the points in s
    L = libgap.Orbit(gp, h, libgap.OnSets) # orbit on these subspaces
    if intersection_size == None:
        G = Graph()
        for x in L: # every pair of points in the subspace is adjacent to each other in G
            G.add_edges(combinations(x, 2))
        return G
    else:
        return Graph([L, lambda i,j: libgap.Size(libgap.Intersection(i,j))==intersection_size],
                        loops=False)

def UnitaryPolarGraph(m, q, algorithm="gap"):
    r"""
    Returns the Unitary Polar Graph `U(m,q)`.

    For more information on Unitary Polar graphs, see the `page of
    Andries Brouwer's website <http://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m,q`` (integers) -- `q` must be a prime power.

    - ``algorithm`` -- if set to 'gap' then the computation is carried via GAP
      library interface, computing totally singular subspaces, which is faster for
      large examples (especially with `q>2`). Otherwise it is done directly.

    EXAMPLES::

        sage: G = graphs.UnitaryPolarGraph(4,2); G
        Unitary Polar Graph U(4, 2); GQ(4, 2): Graph on 45 vertices
        sage: G.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)
        sage: graphs.UnitaryPolarGraph(5,2).is_strongly_regular(parameters=True)
        (165, 36, 3, 9)
        sage: graphs.UnitaryPolarGraph(6,2)    # not tested (long time)
        Unitary Polar Graph U(6, 2): Graph on 693 vertices

    TESTS::

        sage: graphs.UnitaryPolarGraph(4,3, algorithm="gap").is_strongly_regular(parameters=True)
        (280, 36, 8, 4)
        sage: graphs.UnitaryPolarGraph(4,3).is_strongly_regular(parameters=True) # not tested (long time)
        (280, 36, 8, 4)
        sage: graphs.UnitaryPolarGraph(4,3, algorithm="foo")
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm!
    """
    if algorithm == "gap":
        from sage.libs.gap.libgap import libgap
        G = _polar_graph(m, q**2, libgap.GeneralUnitaryGroup(m, q))

    elif algorithm == None: # slow on large examples
        from sage.schemes.projective.projective_space import ProjectiveSpace
        from sage.rings.finite_rings.constructor import FiniteField
        from sage.modules.free_module_element import free_module_element as vector
        from __builtin__ import sum as psum
        Fq = FiniteField(q**2, 'a')
        PG = ProjectiveSpace(m - 1, Fq)
        def P(xx,yy):
            x = vector(xx)
            y = vector(yy)
            return psum(map(lambda j: x[j]*y[m-1-j]**q, xrange(m)))==0  

        V = filter(lambda x: P(x,x), PG)
        G = Graph([V,lambda x,y:  # bottleneck is here, of course:
                     P(vector(x),vector(y))], loops=False)
    else:
        raise ValueError("unknown algorithm!")

    G.relabel()
    G.name("Unitary Polar Graph U" + str((m, q)))
    if m==4:
        G.name(G.name()+'; GQ'+str((q**2,q)))
    if m==5:
        G.name(G.name()+'; GQ'+str((q**2,q**3)))
    return G


def UnitaryDualPolarGraph(m, q):
    r"""
    Returns the Dual Unitary Polar Graph `U(m,q)`.

    For more information on Unitary Dual Polar graphs, see [BCN89]_ and
    Sect. 2.3.1 of [Co81]_.

    INPUT:

    - ``m,q`` (integers) -- `q` must be a prime power.

    EXAMPLES:

    The point graph of a generalized quadrangle of order (8,4)::

        sage: G = graphs.UnitaryDualPolarGraph(5,2); G   # long time
        Unitary Dual Polar Graph DU(5, 2); GQ(8, 4): Graph on 297 vertices
        sage: G.is_strongly_regular(parameters=True)     # long time
        (297, 40, 7, 5)

    Another way to get the  generalized quadrangle of order (2,4)::

        sage: G = graphs.UnitaryDualPolarGraph(4,2); G
        Unitary Dual Polar Graph DU(4, 2); GQ(2, 4): Graph on 27 vertices
        sage: G.is_isomorphic(graphs.OrthogonalPolarGraph(6,2,'-'))
        True

    A bigger graph::

        sage: G = graphs.UnitaryDualPolarGraph(6,2); G   # not tested (long time)
        Unitary Dual Polar Graph DU(6, 2): Graph on 891 vertices
        sage: G.is_distance_regular(parameters=True)     # not tested (long time)
        ([42, 40, 32, None], [None, 1, 5, 21])

    TESTS::

        sage: graphs.UnitaryDualPolarGraph(6,6)
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Error, <subfield> must be a prime or a finite field
    """
    from sage.libs.gap.libgap import libgap
    G = _polar_graph(m, q**2, libgap.GeneralUnitaryGroup(m, q),
            intersection_size=(q**(2*(m//2-1))-1)/(q**2-1))
    G.relabel()
    G.name("Unitary Dual Polar Graph DU" + str((m, q)))
    if m==4:
        G.name(G.name()+'; GQ'+str((q,q**2)))
    if m==5:
        G.name(G.name()+'; GQ'+str((q**3,q**2)))
    return G


def SymplecticDualPolarGraph(m, q):
    r"""
    Returns the Symplectic Dual Polar Graph `DSp(m,q)`.

    For more information on Symplectic Dual Polar graphs, see [BCN89]_ and
    Sect. 2.3.1 of [Co81]_.

    REFERENCE:

    .. [Co81] A. M. Cohen,
      `A synopsis of known distance-regular graphs with large diameters
      <http://persistent-identifier.org/?identifier=urn:nbn:nl:ui:18-6775>`_,
      Stichting Mathematisch Centrum, 1981.

    INPUT:

    - ``m,q`` (integers) -- `q` must be a prime power, and `m` must be even.

    EXAMPLES::

        sage: G = graphs.SymplecticDualPolarGraph(6,3); G       # not tested (long time)
        Symplectic Dual Polar Graph DSp(6, 3): Graph on 1120 vertices
        sage: G.is_distance_regular(parameters=True)            # not tested (long time)
        ([39, 36, 27, None], [None, 1, 4, 13])

    TESTS::

        sage: G = graphs.SymplecticDualPolarGraph(6,2); G
        Symplectic Dual Polar Graph DSp(6, 2): Graph on 135 vertices
        sage: G.is_distance_regular(parameters=True)
        ([14, 12, 8, None], [None, 1, 3, 7])
        sage: graphs.SymplecticDualPolarGraph(6,6)
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Error, <subfield> must be a prime or a finite field
    """
    from sage.libs.gap.libgap import libgap
    G = _polar_graph(m, q, libgap.SymplecticGroup(m, q),
             intersection_size=(q**(m/2-1)-1)/(q-1))

    G.relabel()
    G.name("Symplectic Dual Polar Graph DSp" + str((m, q)))
    if m==4:
        G.name(G.name()+'; GQ'+str((q,q)))
    return G

def TaylorTwographDescendantSRG(q, clique_partition=None):
    r"""
    constructing the descendant graph of the Taylor's two-graph for `U_3(q)`, `q` odd

    This is a strongly regular graph with parameters
    `(v,k,\lambda,\mu)=(q^3, (q^2+1)(q-1)/2, (q-1)^3/4-1, (q^2+1)(q-1)/4)`
    obtained as a two-graph descendant of the
    :func:`Taylor's two-graph <sage.combinat.designs.twographs.taylor_twograph>` `T`.
    This graph admits a partition into cliques of size `q`, which are useful in
    :func:`TaylorTwographSRG <sage.graphs.generators.classical_geometries.TaylorTwographSRG>`,
    a strongly regular graph on `q^3+1` vertices in the
    Seidel switching class of `T`, for which we need `(q^2+1)/2` cliques.
    The cliques are the `q^2` lines on `v_0` of the projective plane containing the unital
    for `U_3(q)`, and intersecting the unital (i.e. the vertices of the graph and the point
    we remove) in `q+1` points. This is all taken from §7E of [BvL84]_.

    INPUT:

    - ``q`` -- an odd prime power

    - ``clique_partition`` -- if ``True``, return `q^2-1` cliques of size `q`
      with empty pairwise intersection. (Removing all of them leaves a clique, too),
      and the point removed from the unital.

    EXAMPLES::

        sage: g=graphs.TaylorTwographDescendantSRG(3); g
        Taylor two-graph descendant SRG: Graph on 27 vertices
        sage: g.is_strongly_regular(parameters=True)
        (27, 10, 1, 5)
        sage: from sage.combinat.designs.twographs import taylor_twograph
        sage: T = taylor_twograph(3)                           # long time
        sage: g.is_isomorphic(T.descendant(T.ground_set()[1])) # long time
        True
        sage: g=graphs.TaylorTwographDescendantSRG(5)    # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)  # not tested (long time)
        (125, 52, 15, 26)

    TESTS::

        sage: g,l,_=graphs.TaylorTwographDescendantSRG(3,clique_partition=True)
        sage: all(map(lambda x: g.is_clique(x), l))
        True
        sage: graphs.TaylorTwographDescendantSRG(4)
        Traceback (most recent call last):
        ...
        ValueError: q must be a an odd prime power 
        sage: graphs.TaylorTwographDescendantSRG(6)
        Traceback (most recent call last):
        ...
        ValueError: q must be a an odd prime power 
    """
    from sage.rings.arith import is_prime_power
    p, k = is_prime_power(q,get_data=True)
    if k==0 or p==2:
       raise ValueError('q must be a an odd prime power')
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module_element import free_module_element as vector
    from sage.rings.finite_rings.integer_mod import mod
    from __builtin__ import sum
    Fq = FiniteField(q**2, 'a')
    PG = ProjectiveSpace(2, Fq)
    def S(xx,yy):
        x = vector(xx)
        y = vector(yy)
        return sum(map(lambda j: x[j]*y[2-j]**q, xrange(3)))

    V = filter(lambda x: S(x,x)==0, PG) # the points of the unital
    v0 = V[0]
    V.remove(v0)
    if mod(q,4)==1:
        G = Graph([V,lambda y,z: not (S(v0,y)*S(y,z)*S(z,v0)).is_square()], loops=False)
    else:
        G = Graph([V,lambda y,z:     (S(v0,y)*S(y,z)*S(z,v0)).is_square()], loops=False)
    G.name("Taylor two-graph descendant SRG")
    if clique_partition:
        lines = map(lambda x: filter(lambda t: vector((Fq.one(), x, Fq.zero()))*vector(t)==0, V),
                     filter(lambda z: z != 0, Fq))
        return (G, lines, v0)
    else:
        return G

def TaylorTwographSRG(q):
    r"""
    constructing a strongly regular graph from the Taylor's two-graph for `U_3(q)`, `q` odd

    This is a strongly regular graph with parameters
    `(v,k,\lambda,\mu)=(q^3+1, q(q^2+1)/2, (q^2+3)(q-1)/4, (q^2+1)(q+1)/4)`
    in the Seidel switching class of
    :func:`Taylor two-graph <sage.combinat.designs.twographs.taylor_twograph>`.
    Details are in §7E of [BvL84]_.

    .. SEEALSO::

        * :func:`TaylorTwographDescendantSRG <sage.graphs.generators.classical_geometries.TaylorTwographDescendantSRG>`

    EXAMPLES::

        sage: t=graphs.TaylorTwographSRG(3); t
        Taylor two-graph SRG: Graph on 28 vertices
        sage: t.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)

    """
    from sage.misc.flatten import flatten
    from sage.graphs.generators.classical_geometries import TaylorTwographDescendantSRG
    H, l, v0 = TaylorTwographDescendantSRG(q, clique_partition=True)
    G = H.union(Graph([[v0], lambda x, y: x != y])) # to make sure vertices are not relabeled
    G.seidel_switching(flatten(l[:(q**2+1)/2]))
    G.name("Taylor two-graph SRG")
    return G
