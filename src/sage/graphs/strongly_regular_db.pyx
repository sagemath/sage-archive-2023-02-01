# -*- coding: utf-8 -*-
r"""
Database of strongly regular graphs

This module manages a database associating to a set of four integers
`(v,k,\lambda,\mu)` a strongly regular graphs with these parameters, when one
exists.

Using Andries Brouwer's `database of strongly regular graphs
<http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__, it can also return
non-existence results. Note that some constructions are missing, and that some
strongly regular graphs that exist in the database cannot be automatically built
by Sage. Help us if you know any.

.. NOTE::

    Any missing/incorrect information in the database must be reported to
    `Andries E. Brouwer <http://www.win.tue.nl/~aeb/>`__ directly, in order to
    have a unique and updated source of information.

REFERENCES:

.. [BvL84] \A. Brouwer, J van Lint,
   Strongly regular graphs and partial geometries,
   Enumeration and design,
   (Waterloo, Ont., 1982) (1984): 85-122.
   http://oai.cwi.nl/oai/asset/1817/1817A.pdf

Functions
---------
"""
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown
from sage.arith.all import is_square, is_prime_power, divisors
from sage.misc.cachefunc import cached_function
from sage.combinat.designs.orthogonal_arrays import orthogonal_array
from sage.combinat.designs.bibd import balanced_incomplete_block_design
from sage.graphs.graph import Graph
from libc.math cimport sqrt, floor
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.coding.linear_code import LinearCode
from sage.rings.sum_of_squares cimport two_squares_c
from libc.stdint cimport uint_fast32_t

cdef dict _brouwer_database = None
_small_srg_database = None

@cached_function
def is_paley(int v,int k,int l,int mu):
    r"""
    Test whether some Paley graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_paley
        sage: t = is_paley(13,6,2,3); t
        (..., 13)
        sage: g = t[0](*t[1:]); g
        Paley graph with parameter 13: Graph on 13 vertices
        sage: g.is_strongly_regular(parameters=True)
        (13, 6, 2, 3)
        sage: t = is_paley(5,5,5,5); t
    """
    if (v%4 == 1 and is_prime_power(v) and
        k   == (v-1)/2 and
        l   == (v-5)/4 and
        mu  == (v-1)/4):
        from sage.graphs.generators.families import PaleyGraph
        return (PaleyGraph,v)

@cached_function
def is_mathon_PC_srg(int v,int k,int l,int mu):
    r"""
    Test whether some Mathon's Pseudocyclic s.r.g. is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    .. TODO::

        The current implementation only gives a subset of all possible graphs that can be
        obtained using this construction. A  full implementation should rely on a database
        of conference matrices (or, equivalently, on a database of s.r.g.'s with parameters
        `(4t+1,2t,t-1,t)`. Currently we make an extra assumtion that `4t+1` is a prime power.
        The first case where we miss a construction is `t=11`, where we could (recursively)
        use the graph for `t=1` to construct a graph on 83205 vertices.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_mathon_PC_srg
        sage: t = is_mathon_PC_srg(45,22,10,11); t
        (..., 1)
        sage: g = t[0](*t[1:]); g
        Mathon's PC SRG on 45 vertices: Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)
        (45, 22, 10, 11)

    TESTS::

        sage: t = is_mathon_PC_srg(5,5,5,5); t
        sage: mu = 1895  # t=5 case -- the construction cannot work
        sage: t = is_mathon_PC_srg(4*mu+1,2*mu,mu-1,mu); t
    """
    cdef int t
    if (v%4 == 1 and
        k   == (v-1)/2 and
        l   == (v-5)/4 and
        mu  == (v-1)/4):
        from sage.rings.integer_ring import ZZ
        K = ZZ['x']
        x = K.gen()
        rpoly = filter(lambda w: w[0]>0, (x*(4*x*(4*x-1)-1)-mu).roots())
        if rpoly != []:
            t = rpoly[0][0]
            if (is_prime_power(4*t-1) and
                is_prime_power(4*t+1)): # extra assumption in TODO!
                from sage.graphs.generators.families import \
                                    MathonPseudocyclicStronglyRegularGraph
                return (MathonPseudocyclicStronglyRegularGraph,t)

@cached_function
def is_orthogonal_array_block_graph(int v,int k,int l,int mu):
    r"""
    Test whether some (pseudo)Orthogonal Array graph is `(v,k,\lambda,\mu)`-strongly regular.

    We know how to construct graphs with parameters of an Orthogonal Array (`OA(m,n)`),
    also known as Latin squares graphs `L_m(n)`, in several cases where no orthogonal
    array is known, or even in some cases for which they are known not to exist.

    Such graphs are usually called pseudo-Latin squares graphs. Namely, Sage
    can construct a graph with parameters of an `OA(m,n)`-graph whenever there
    exists a skew-Hadamard matrix of order `n+1`, and `m=(n+1)/2` or `m=(n-1)/2`.
    The construction in the former case is due to Goethals-Seidel [BvL84]_, and in the
    latter case due to Pasechnik [Pa92]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_orthogonal_array_block_graph
        sage: t = is_orthogonal_array_block_graph(64, 35, 18, 20); t
        (..., 5, 8)
        sage: g = t[0](*t[1:]); g
        OA(5,8): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)
        sage: t=is_orthogonal_array_block_graph(225,98,43,42); t
        (..., 4)
        sage: g = t[0](*t[1:]); g
        Pasechnik Graph_4: Graph on 225 vertices
        sage: g.is_strongly_regular(parameters=True)
        (225, 98, 43, 42)
        sage: t=is_orthogonal_array_block_graph(225,112,55,56); t
        (..., 4)
        sage: g = t[0](*t[1:]); g
        skewhad^2_4: Graph on 225 vertices
        sage: g.is_strongly_regular(parameters=True)
        (225, 112, 55, 56)

        sage: t = is_orthogonal_array_block_graph(5,5,5,5); t

    REFERENCE:

    .. [Pa92] \D. V. Pasechnik,
      Skew-symmetric association schemes with two classes and strongly
      regular graphs of type `L_{2n-1}(4n- 1)`,
      Acta Applicandaie Math. 29(1992), 129-138
    """
    # notations from
    # http://www.win.tue.nl/~aeb/graphs/OA.html
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    try:
        m, n = latin_squares_graph_parameters(v,k,l,mu)
    except:
        return
    if orthogonal_array(m,n,existence=True):
        from sage.graphs.generators.intersection import OrthogonalArrayBlockGraph
        return (lambda m,n : OrthogonalArrayBlockGraph(m, n), m,n)

    elif n>2 and skew_hadamard_matrix(n+1, existence=True):
        if m==(n+1)/2:
            from sage.graphs.generators.families import SquaredSkewHadamardMatrixGraph as G
        elif m==(n-1)/2:
            from sage.graphs.generators.families import PasechnikGraph as G
        else:
            return
        return (G, (n+1)/4)

@cached_function
def is_johnson(int v,int k,int l,int mu):
    r"""
    Test whether some Johnson graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_johnson
        sage: t = is_johnson(10,6,3,4); t
        (..., 5)
        sage: g = t[0](*t[1:]); g
        Johnson graph with parameters 5,2: Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)
        (10, 6, 3, 4)

        sage: t = is_johnson(5,5,5,5); t
    """
    # Using notations of http://www.win.tue.nl/~aeb/graphs/Johnson.html
    #
    # J(n,m) has parameters v = m(m – 1)/2, k = 2(m – 2), λ = m – 2, μ = 4.
    m = l + 2
    if (mu == 4 and
        k  == 2*(m-2) and
        v  == m*(m-1)/2):
        from sage.graphs.generators.families import JohnsonGraph
        return (lambda m: JohnsonGraph(m,2), m)

@cached_function
def is_steiner(int v,int k,int l,int mu):
    r"""
    Test whether some Steiner graph is `(v,k,\lambda,\mu)`-strongly regular.

    A Steiner graph is the intersection graph of a Steiner set system. For more
    information, see http://www.win.tue.nl/~aeb/graphs/S.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_steiner
        sage: t = is_steiner(26,15,8,9); t
        (..., 13, 3)
        sage: g = t[0](*t[1:]); g
        Intersection Graph: Graph on 26 vertices
        sage: g.is_strongly_regular(parameters=True)
        (26, 15, 8, 9)

        sage: t = is_steiner(5,5,5,5); t
    """
    # Using notations from http://www.win.tue.nl/~aeb/graphs/S.html
    #
    # The block graph of a Steiner 2-design S(2,m,n) has parameters:
    # v = n(n-1)/m(m-1), k = m(n-m)/(m-1), λ = (m-1)^2 + (n-1)/(m–1)–2, μ = m^2.
    if mu <= 1 or not is_square(mu):
        return
    m = int(sqrt(mu))
    n = (k*(m-1))//m+m
    if (v == (n*(n-1))/(m*(m-1)) and
        k == m*(n-m)/(m-1) and
        l == (m-1)**2 + (n-1)/(m-1)-2 and
        balanced_incomplete_block_design(n,m,existence=True)):
        from sage.graphs.generators.intersection import IntersectionGraph
        return (lambda n,m: IntersectionGraph(map(frozenset,balanced_incomplete_block_design(n,m))),n,m)

@cached_function
def is_affine_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Affine Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/VO.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_affine_polar
        sage: t = is_affine_polar(81,32,13,12); t
        (..., 4, 3)
        sage: g = t[0](*t[1:]); g
        Affine Polar Graph VO^+(4,3): Graph on 81 vertices
        sage: g.is_strongly_regular(parameters=True)
        (81, 32, 13, 12)

        sage: t = is_affine_polar(5,5,5,5); t
    """
    # Using notations from http://www.win.tue.nl/~aeb/graphs/VO.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    if (not is_square(v) or
        not is_prime_power(v)):
        return
    prime,power = is_prime_power(v,get_data=True)
    if power%2:
        return
    for e in divisors(power/2):
        q = prime**(power//(2*e))
        assert v == q**(2*e)
        if (k == (q**(e-1) + 1)*(q**e-1) and
            l == q*(q**(e-2) + 1)*(q**(e-1)-1)+q-2 and
            mu== q**(e-1)*(q**(e-1) + 1)):
            from sage.graphs.generators.classical_geometries import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='+'),2*e,q)
        if (k == (q**(e-1) - 1)*(q**e+1) and
            l == q*(q**(e-2)- 1)*(q**(e-1)+1)+q-2 and
            mu== q**(e-1)*(q**(e-1) - 1)):
            from sage.graphs.generators.classical_geometries import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='-'),2*e,q)

@cached_function
def is_orthogonal_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Orthogonal Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_orthogonal_polar
        sage: t = is_orthogonal_polar(85, 20, 3, 5); t
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: g = t[0](*t[1:]); g
        Orthogonal Polar Graph O(5, 4): Graph on 85 vertices
        sage: g.is_strongly_regular(parameters=True)
        (85, 20, 3, 5)

        sage: t = is_orthogonal_polar(5,5,5,5); t

    TESTS:

    All of ``O(2m+1,q)``, ``O^+(2m,q)`` and ``O^-(2m,q)`` appear::

        sage: is_orthogonal_polar(85, 20, 3, 5)
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: is_orthogonal_polar(119,54,21,27)
        (<function OrthogonalPolarGraph at ...>, 8, 2, '-')
        sage: is_orthogonal_polar(130,48,20,16)
        (<function OrthogonalPolarGraph at ...>, 6, 3, '+')

    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q_pow_m_minus_one = -s-1 if abs(s) > r else r+1

    if is_prime_power(q_pow_m_minus_one):
        prime,power = is_prime_power(q_pow_m_minus_one,get_data=True)
        for d in divisors(power):
            q = prime**d
            m = (power//d)+1

            # O(2m+1,q)
            if (v == (q**(2*m)-1)/(q-1)              and
                k == q*(q**(2*m-2)-1)/(q-1)          and
                l == q**2*(q**(2*m-4)-1)/(q-1) + q-1 and
                mu== (q**(2*m-2)-1)/(q-1)):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m+1, q, "")

            # O^+(2m,q)
            if (v ==   (q**(2*m-1)-1)/(q-1) + q**(m-1)   and
                k == q*(q**(2*m-3)-1)/(q-1) + q**(m-1) and
                k == q**(2*m-3) + l + 1                  and
                mu== k/q):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "+")

            # O^+(2m+1,q)
            if (v ==   (q**(2*m-1)-1)/(q-1) - q**(m-1)   and
                k == q*(q**(2*m-3)-1)/(q-1) - q**(m-1) and
                k == q**(2*m-3) + l + 1                  and
                mu== k/q):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "-")

@cached_function
def is_goethals_seidel(int v,int k,int l,int mu):
    r"""
    Test whether some
    :func:`~sage.graphs.graph_generators.GraphGenerators.GoethalsSeidelGraph` graph is
    `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_goethals_seidel
        sage: t = is_goethals_seidel(28, 15, 6, 10); t
        [<function GoethalsSeidelGraph at ...>, 3, 3]
        sage: g = t[0](*t[1:]); g
        Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)

        sage: t = is_goethals_seidel(256, 135, 70, 72); t
        [<function GoethalsSeidelGraph at ...>, 2, 15]
        sage: g = t[0](*t[1:]); g
        Graph on 256 vertices
        sage: g.is_strongly_regular(parameters=True)
        (256, 135, 70, 72)

        sage: t = is_goethals_seidel(5,5,5,5); t

    TESTS::

        sage: for p in [(16, 9, 4, 6), (28, 15, 6, 10), (64, 35, 18, 20), (120, 63, 30, 36),
        ....:           (144, 77, 40, 42), (256, 135, 70, 72), (400, 209, 108, 110),
        ....:           (496, 255, 126, 136), (540, 275, 130, 150), (576, 299, 154, 156),
        ....:           (780, 399, 198, 210), (784, 405, 208, 210), (976, 495, 238, 264)]:
        ....:     print is_goethals_seidel(*p)
        [<function GoethalsSeidelGraph at ...>, 2, 3]
        [<function GoethalsSeidelGraph at ...>, 3, 3]
        [<function GoethalsSeidelGraph at ...>, 2, 7]
        [<function GoethalsSeidelGraph at ...>, 3, 7]
        [<function GoethalsSeidelGraph at ...>, 2, 11]
        [<function GoethalsSeidelGraph at ...>, 2, 15]
        [<function GoethalsSeidelGraph at ...>, 2, 19]
        [<function GoethalsSeidelGraph at ...>, 3, 15]
        [<function GoethalsSeidelGraph at ...>, 5, 11]
        [<function GoethalsSeidelGraph at ...>, 2, 23]
        [<function GoethalsSeidelGraph at ...>, 3, 19]
        [<function GoethalsSeidelGraph at ...>, 2, 27]
        [<function GoethalsSeidelGraph at ...>, 5, 15]
    """
    from sage.combinat.designs.bibd import balanced_incomplete_block_design
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix

    # here we guess the parameters v_bibd,k_bibd and r_bibd of the block design
    #
    # - the number of vertices v is equal to v_bibd*(r_bibd+1)
    # - the degree k of the graph is equal to k=(v+r_bibd-1)/2

    r_bibd = k - (v-1-k)
    v_bibd = v//(r_bibd+1)
    k_bibd = (v_bibd-1)/r_bibd + 1 if r_bibd>0 else -1

    if (v   == v_bibd*(r_bibd+1)                  and
        2*k == v+r_bibd-1                         and
        4*l == -2*v + 6*k -v_bibd -k_bibd         and
        hadamard_matrix(r_bibd+1, existence=True) and
        balanced_incomplete_block_design(v_bibd, k_bibd, existence = True)):
        from sage.graphs.generators.families import GoethalsSeidelGraph
        return [GoethalsSeidelGraph, k_bibd, r_bibd]

@cached_function
def is_NOodd(int v,int k,int l,int mu):
    r"""
    Test whether some NO^e(2n+1,q) graph is `(v,k,\lambda,\mu)`-strongly regular.

    Here `q>2`, for in the case `q=2` this graph is complete. For more
    information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`
    and Sect. 7.C of [BvL84]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NOodd
        sage: t = is_NOodd(120, 51, 18, 24); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '-')
        sage: g = t[0](*t[1:]); g
        NO^-(5, 4): Graph on 120 vertices
        sage: g.is_strongly_regular(parameters=True)
        (120, 51, 18, 24)

    TESTS:

    All of ``NO^+(2m+1,q)`` and ``NO^-(2m+1,q)`` appear::

        sage: t = is_NOodd(120, 51, 18, 24); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '-')
        sage: t = is_NOodd(136, 75, 42, 40); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '+')
        sage: t=is_NOodd(378, 260, 178, 180); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 7, 3, '+')
        sage: t=is_NOodd(45, 32, 22, 24); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 3, '+')
        sage: t=is_NOodd(351, 224, 142, 144); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 7, 3, '-')
        sage: t = is_NOodd(325, 144, 68, 60); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '+')
        sage: t = is_NOodd(300, 104, 28, 40); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '-')
        sage: t = is_NOodd(5,5,5,5); t
    """
    cdef int n, q
    r,s = eigenvalues(v,k,l,mu) # -eq^(n-1)-1 and eq^(n-1)(q-2)-1; q=3 is special case
    if r is None:
        return
    r += 1
    s += 1
    if abs(r)>abs(s):
        (r,s) = (s,r) # r=-eq^(n-1) s= eq^(n-1)(q-2)
    q = 2 - s/r
    p, t = is_prime_power(q, get_data=True)
    pp, kk = is_prime_power(abs(r), get_data=True)
    if p == pp and t != 0:
        n  = kk/t + 1
        e = 1 if v  == (q**n)*(q**n+1)/2 else -1
        if (v  == (q**n)*(q**n+e)/2                 and
            k  == (q**n-e)*(q**(n-1)+e)             and
            l  == 2*(q**(2*n-2)-1)+e*q**(n-1)*(q-1) and
            mu == 2*q**(n-1)*(q**(n-1)+e)):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n+1, q, '+' if e==1 else '-')

@cached_function
def is_NOperp_F5(int v,int k,int l,int mu):
    r"""
    Test whether some NO^e,perp(2n+1,5) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`
    and Sect. 7.D of [BvL84]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NOperp_F5
        sage: t = is_NOperp_F5(10, 3, 0, 1); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 3, 5, '-', 1)
        sage: g = t[0](*t[1:]); g
        NO^-,perp(3, 5): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)
        (10, 3, 0, 1)

    TESTS:

    All of ``NO^+,perp(2m+1,5)`` and ``NO^-,perp(2m+1,5)`` appear::

        sage: t = is_NOperp_F5(325, 60, 15, 10); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '+', 1)
        sage: t = is_NOperp_F5(300, 65, 10, 15); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '-', 1)
        sage: t = is_NOperp_F5(5,5,5,5); t
    """
    cdef int n
    r,s = eigenvalues(v,k,l,mu) # 2*e*5**(n-1), -e*5**(n-1); note exceptional case n=1
    if r is None:
        return
    if abs(r)<abs(s):
        (r,s) = (s,r)
    e = 1 if s<0 else -1
    p, n = is_prime_power(abs(s), get_data=True)
    if (5 == p and n != 0) or (abs(r)==2 and abs(s)==1):
        n += 1
        if (v  == (5**n)*(5**n+e)/2           and
            k  == (5**n-e)*5**(n-1)/2         and
            l  == 5**(n-1)*(5**(n-1)+e)/2     and
            mu == 5**(n-1)*(5**(n-1)-e)/2):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n+1, 5, '+' if e==1 else '-', 1)

@cached_function
def is_NO_F2(int v,int k,int l,int mu):
    r"""
    Test whether some NO^e,perp(2n,2) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NO_F2
        sage: t = is_NO_F2(10, 3, 0, 1); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 4, 2, '-')
        sage: g = t[0](*t[1:]); g
        NO^-(4, 2): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)
        (10, 3, 0, 1)

    TESTS:

    All of ``NO^+(2m,2)`` and ``NO^-(2m,2)`` appear::

        sage: t = is_NO_F2(36, 15, 6, 6); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 2, '-')
        sage: t = is_NO_F2(28, 15, 6, 10); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 2, '+')
        sage: t = is_NO_F2(5,5,5,5); t
    """
    cdef int n, e, p
    p, n = is_prime_power(k+1, get_data=True) # k+1==2**(2*n-2)
    if 2 == p and n != 0 and n % 2 == 0:
        n = (n+2)/2
        e = (2**(2*n-1)-v)/2**(n-1)
        if (abs(e) == 1                           and
            v  == 2**(2*n-1)-e*2**(n-1)           and
            k  == 2**(2*n-2)-1                    and
            l  == 2**(2*n-3)-2                    and
            mu == 2**(2*n-3)+e*2**(n-2)):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n, 2, '+' if e==1 else '-')

@cached_function
def is_NO_F3(int v,int k,int l,int mu):
    r"""
    Test whether some NO^e,perp(2n,3) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NO_F3
        sage: t = is_NO_F3(15, 6, 1, 3); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 4, 3, '-')
        sage: g = t[0](*t[1:]); g
        NO^-(4, 3): Graph on 15 vertices
        sage: g.is_strongly_regular(parameters=True)
        (15, 6, 1, 3)

    TESTS:

    All of ``NO^+(2m,3)`` and ``NO^-(2m,3)`` appear::

        sage: t = is_NO_F3(126, 45, 12, 18); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 3, '-')
        sage: t = is_NO_F3(117, 36, 15, 9); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 3, '+')
        sage: t = is_NO_F3(5,5,5,5); t
    """
    cdef int n, e, p
    r,s = eigenvalues(v,k,l,mu) # e*3**(n-1), -e*3**(n-2)
    if r is None:
        return
    if abs(r)<abs(s):
        (r,s) = (s,r)
    e = 1 if r>0 else -1
    p, n = is_prime_power(abs(r), get_data=True)
    if (3 == p and n != 0):
        n += 1
        if (v  == 3**(n-1)*(3**n-e)/2           and
            k  == 3**(n-1)*(3**(n-1)-e)/2           and
            l  == 3**(n-2)*(3**(n-1)+e)/2           and
            mu == 3**(n-1)*(3**(n-2)-e)/2):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n, 3, '+' if e==1 else '-')

@cached_function
def is_NU(int v,int k,int l,int mu):
    r"""
    Test whether some NU(n,q)-graph, is `(v,k,\lambda,\mu)`-strongly regular.

    Note that n>2; for n=2 there is no s.r.g. For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicUnitaryPolarGraph`
    and series C14 in [Hu75]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NU
        sage: t = is_NU(40, 27, 18, 18); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 4, 2)
        sage: g = t[0](*t[1:]); g
        NU(4, 2): Graph on 40 vertices
        sage: g.is_strongly_regular(parameters=True)
        (40, 27, 18, 18)

    TESTS::

        sage: t = is_NU(176, 135, 102, 108); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 5, 2)
        sage: t = is_NU(540, 224, 88, 96); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 4, 3)
        sage: t = is_NU(208, 75, 30, 25); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 3, 4)
        sage: t = is_NU(5,5,5,5); t
    """
    cdef int n, q, e               # special cases: n=3 or q=2
    r,s = eigenvalues(v,k,l,mu) #r,s = eq^{n-2} - 1, -e(q^2-q-1)q^{n-3} - 1, e=(-1)^n
    if r is None:
        return
    r += 1
    s += 1
    if abs(r)>abs(s):
        (r,s) = (s,r)
    p, t = is_prime_power(abs(r), get_data=True)
    if p==2: # it can be that q=2, then we'd have r>s now
        pp, kk = is_prime_power(abs(s), get_data=True)
        if pp==2 and kk>0:
            (r,s) = (s,r)
            p, t = is_prime_power(abs(r), get_data=True)
    if r==1:
        return
    kr = k/(r-1) # eq^{n-1}+1
    e = 1 if kr>0 else -1
    q = (kr-1)/r
    pp, kk = is_prime_power(q, get_data=True)
    if p == pp and kk != 0:
        n  = t/kk + 2
        if (v  == q**(n-1)*(q**n - e)/(q + 1)               and
            k  == (q**(n-1) + e)*(q**(n-2) - e)             and
            l  == q**(2*n-5)*(q+1) - e*q**(n-2)*(q-1) - 2  and
            mu == q**(n-3)*(q + 1)*(q**(n-2) - e)):
            from sage.graphs.generators.classical_geometries import NonisotropicUnitaryPolarGraph
            return (NonisotropicUnitaryPolarGraph, n, q)

@cached_function
def is_haemers(int v,int k,int l,int mu):
    r"""
    Test whether some HaemersGraph graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`~sage.graphs.graph_generators.GraphGenerators.HaemersGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_haemers
        sage: t = is_haemers(96, 19, 2, 4); t
        (<function HaemersGraph at ...>, 4)
        sage: g = t[0](*t[1:]); g
        Haemers(4): Graph on 96 vertices
        sage: g.is_strongly_regular(parameters=True)
        (96, 19, 2, 4)

    TESTS::

        sage: t = is_haemers(5,5,5,5); t
    """
    cdef int q, n, p
    p, n = is_prime_power(mu, get_data=True)
    q = mu
    if 2 == p and n != 0:
        if (v  == q**2*(q+2)     and
            k  == q*(q+1)-1      and
            l  == q-2):
            from sage.graphs.generators.classical_geometries import HaemersGraph
            return (HaemersGraph, q)

@cached_function
def is_cossidente_penttila(int v,int k,int l,int mu):
    r"""
    Test whether some CossidentePenttilaGraph graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`~sage.graphs.graph_generators.GraphGenerators.CossidentePenttilaGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_cossidente_penttila
        sage: t =  is_cossidente_penttila(378, 52, 1, 8); t
        (<function CossidentePenttilaGraph at ...>, 5)
        sage: g = t[0](*t[1:]); g                      # optional - gap_packages
        CossidentePenttila(5): Graph on 378 vertices
        sage: g.is_strongly_regular(parameters=True)   # optional - gap_packages
        (378, 52, 1, 8)

    TESTS::

        sage: t =  is_cossidente_penttila(56,10,0,2); t
        (<function CossidentePenttilaGraph at ...>, 3)
        sage: t =  is_cossidente_penttila(1376,150,2,18); t
        (<function CossidentePenttilaGraph at ...>, 7)
        sage: t = is_cossidente_penttila(5,5,5,5); t
    """
    cdef int q, n, p
    q = 2*l+3
    p, n = is_prime_power(q, get_data=True)
    if 2 < p and n != 0:
        if (v  == (q**3+1)*(q+1)/2     and
            k  == (q**2+1)*(q-1)/2      and
            mu  == (q-1)**2/2):
            from sage.graphs.generators.classical_geometries import CossidentePenttilaGraph
            return (CossidentePenttilaGraph, q)

@cached_function
def is_complete_multipartite(int v,int k,int l,int mu):
    r"""
    Test whether some complete multipartite graph is `(v,k,\lambda,\mu)`-strongly regular.

    Any complete multipartite graph with parts of the same size is strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_complete_multipartite
        sage: t = is_complete_multipartite(12,8,4,8); t
        (<cyfunction is_complete_multipartite.<locals>.CompleteMultipartiteSRG at ...>,
         3,
         4)
        sage: g = t[0](*t[1:]); g
        Multipartite Graph with set sizes [4, 4, 4]: Graph on 12 vertices
        sage: g.is_strongly_regular(parameters=True)
        (12, 8, 4, 8)

    TESTS::

        sage: t = is_complete_multipartite(5,5,5,5); t
        sage: t = is_complete_multipartite(11,8,4,8); t
        sage: t = is_complete_multipartite(20,16,12,16);
        sage: g = t[0](*t[1:]); g
        Multipartite Graph with set sizes [4, 4, 4, 4, 4]: Graph on 20 vertices
        sage: g.is_strongly_regular(parameters=True)
        (20, 16, 12, 16)
    """
    if v>k:
        r = v/(v-k) # number of parts (of size v-k each)
        if l==(v-k)*(r-2) and k==mu and v == r*(v-k):
            from sage.graphs.generators.basic import CompleteMultipartiteGraph
            def CompleteMultipartiteSRG(nparts, partsize):
                return CompleteMultipartiteGraph([partsize]*nparts)
            return (CompleteMultipartiteSRG, r, v-k)

@cached_function
def is_polhill(int v,int k,int l,int mu):
    r"""
    Test whether some graph from [Polhill09]_ is `(1024,k,\lambda,\mu)`-strongly regular.

    .. NOTE::

        This function does not actually explore *all* strongly regular graphs
        produced in [Polhill09]_, but only those on 1024 vertices.

        John Polhill offered his help if we attempt to write a code to guess,
        given `(v,k,\lambda,\mu)`, which of his construction must be applied to
        find the graph.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_polhill
        sage: t = is_polhill(1024, 231,  38,  56); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: g = t[0](*t[1:]); g                    # not tested (too long)
        Graph on 1024 vertices
        sage: g.is_strongly_regular(parameters=True) # not tested (too long)
        (1024, 231, 38, 56)
        sage: t = is_polhill(1024, 264,  56,  72); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 297,  76,  90); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 330,  98, 110); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 462, 206, 210); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]

    REFERENCE:

    .. [Polhill09] \J. Polhill,
       Negative Latin square type partial difference sets and
       amorphic association schemes with Galois rings,
       Journal of Combinatorial Designs 17, no. 3 (2009): 266-282.
       http://onlinelibrary.wiley.com/doi/10.1002/jcd.20206/abstract
    """
    if (v,k,l,mu) not in [(1024, 231,  38,  56),
                          (1024, 264,  56,  72),
                          (1024, 297,  76,  90),
                          (1024, 330,  98, 110),
                          (1024, 462, 206, 210)]:
        return

    from itertools import product
    from sage.categories.cartesian_product import cartesian_product
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from copy import copy

    def additive_cayley(vertices):
        g = Graph()
        g.add_vertices(vertices[0].parent())
        edges = [(x,x+vv)
                 for vv in set(vertices)
                 for x in g]
        g.add_edges(edges)
        g.relabel()
        return g

    # D is a Partial Difference Set of (Z4)^2, see section 2.
    G = cartesian_product([IntegerModRing(4),IntegerModRing(4)])
    D = [
        [(2,0),(0,1),(0,3),(1,1),(3,3)],
        [(1,0),(3,0),(0,2),(1,3),(3,1)],
        [(1,2),(3,2),(2,1),(2,3),(2,2)]
        ]
    D = [map(G,x) for x in D]

    # The K_i are hyperplanes partitionning the nonzero elements of
    # GF(2^s)^2. See section 6.
    s = 3
    G1 = GF(2**s,'x')
    Gp = cartesian_product([G1,G1])
    K = [Gp((x,1)) for x in G1]+[Gp((1,0))]
    K = [[x for x in Gp if x[0]*uu+x[1]*vv == 0] for (uu,vv) in K]

    # We now define the P_{i,j}. see section 6.

    P = {}
    P[0,1] = range((-1) + 1                  , 2**(s-2)+1)
    P[1,1] = range((-1) + 2**(s-2)+2         , 2**(s-1)+1)
    P[2,1] = range((-1) + 2**(s-1)+2         , 2**(s-1)+2**(s-2)+1)
    P[3,1] = range((-1) + 2**(s-1)+2**(s-2)+2, 2**(s)+1)

    P[0,2] = range((-1) + 2**(s-2)+2         , 2**(s-1)+2)
    P[1,2] = range((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+2)
    P[2,2] = range((-1) + 2**(s-1)+2**(s-2)+3, 2**(s)+1) + [0]
    P[3,2] = range((-1) + 2                  , 2**(s-2)+1)

    P[0,3] = range((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+3)
    P[1,3] = range((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1) + [0,1]
    P[2,3] = range((-1) + 3                  , 2**(s-2)+2)
    P[3,3] = range((-1) + 2**(s-2)+3         , 2**(s-1)+2)

    P[0,4] = range((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1)
    P[1,4] = range((-1) + 3                  , 2**(s-2)+1) + [2**(s-1)+1,2**(s-1)+2**(s-2)+2]
    P[2,4] = range((-1) + 2**(s-2)+3         , 2**(s-1)+1) + [2**(s-1)+2**(s-2)+1,1]
    P[3,4] = range((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+1) + [2**(s-2)+1,0]

    R = {x:copy(P[x]) for x in P}

    for x in P:
        P[x] = [K[i] for i in P[x]]
        P[x] = set(sum(P[x],[])).difference([Gp((0,0))])

    P[1,4].add(Gp((0,0)))
    P[2,4].add(Gp((0,0)))
    P[3,4].add(Gp((0,0)))

    # We now define the R_{i,j}. see *end* of section 6.

    R[0,3] = range((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+2)
    R[1,3] = range((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1) + [0,1,2**(s-1)+2**(s-2)+2]
    R[0,4] = range((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1) + [2**(s-1)+2**(s-2)+2]
    R[1,4] = range((-1) + 3                  , 2**(s-2)+1) + [2**(s-1)+1]

    for x in R:
        R[x] = [K[i] for i in R[x]]
        R[x] = set(sum(R[x],[])).difference([Gp((0,0))])

    R[1,3].add(Gp((0,0)))
    R[2,4].add(Gp((0,0)))
    R[3,4].add(Gp((0,0)))

    # Dabcd = Da, Db, Dc, Dd (cf. p273)
    # D1234 = D1, D2, D3, D4 (cf. p276)
    Dabcd = []
    D1234 = []

    Gprod = cartesian_product([G,Gp])
    for DD,PQ in [(Dabcd,P), (D1234,R)]:
        for i in range(1,5):
            Dtmp = [product([G.zero()],PQ[0,i]),
                    product(D[0],PQ[1,i]),
                    product(D[1],PQ[2,i]),
                    product(D[2],PQ[3,i])]
            Dtmp = map(set,Dtmp)
            Dtmp = map(Gprod,sum(map(list,Dtmp),[]))
            DD.append(Dtmp)

    # Now that we have the data, we can return the graphs.
    if k == 231:
        return [lambda :additive_cayley(Dabcd[0])]
    if k == 264:
        return [lambda :additive_cayley(D1234[2])]
    if k == 297:
        return [lambda :additive_cayley(D1234[0]+D1234[1]+D1234[2]).complement()]
    if k == 330:
        return [lambda :additive_cayley(Dabcd[0]+Dabcd[1]+Dabcd[2]).complement()]
    if k == 462:
        return [lambda :additive_cayley(Dabcd[0]+Dabcd[1])]

def is_RSHCD(int v,int k,int l,int mu):
    r"""
    Test whether some RSHCD graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see :func:`SRG_from_RSHCD`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_RSHCD
        sage: t = is_RSHCD(64,27,10,12); t
        [<built-in function SRG_from_RSHCD>, 64, 27, 10, 12]
        sage: g = t[0](*t[1:]); g
        Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 27, 10, 12)

    """
    if SRG_from_RSHCD(v,k,l,mu,existence=True):
        return [SRG_from_RSHCD,v,k,l,mu]

def SRG_from_RSHCD(v,k,l,mu, existence=False,check=True):
    r"""
    Return a `(v,k,l,mu)`-strongly regular graph from a RSHCD

    This construction appears in 8.D of [BvL84]_. For more information, see
    :func:`~sage.combinat.matrices.hadamard_matrix.regular_symmetric_hadamard_matrix_with_constant_diagonal`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    - ``existence`` (boolean) -- whether to return a graph or to test if Sage
      can build such a graph.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES:

    some graphs ::

        sage: from sage.graphs.strongly_regular_db import SRG_from_RSHCD
        sage: SRG_from_RSHCD(784, 0, 14, 38, existence=True)
        False
        sage: SRG_from_RSHCD(784, 377, 180, 182, existence=True)
        True
        sage: SRG_from_RSHCD(144, 65, 28, 30)
        Graph on 144 vertices

    an example with vertex-transitive automorphism group, found during the
    implementation of the case `v=324` ::

        sage: G=SRG_from_RSHCD(324,152,70,72)  # long time
        sage: a=G.automorphism_group()         # long time
        sage: a.order()                        # long time
        2592
        sage: len(a.orbits())                  # long time
        1

    TESTS::

        sage: SRG_from_RSHCD(784, 0, 14, 38)
        Traceback (most recent call last):
        ...
        ValueError: I do not know how to build a (784, 0, 14, 38)-SRG from a RSHCD

    """
    from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal
    sgn = lambda x: 1 if x>=0 else -1
    n = v
    a = (n-4*mu)//2
    e = 2*k - n + 1 + a
    t = abs(a//2)

    if (e**2 == 1              and
        k == (n-1-a+e)/2       and
        l == (n-2*a)/4 - (1-e) and
        mu== (n-2*a)/4         and
        regular_symmetric_hadamard_matrix_with_constant_diagonal(n,sgn(a)*e,existence=True)):
        if existence:
            return True
        from sage.matrix.constructor import identity_matrix as I
        from sage.matrix.constructor import ones_matrix     as J

        H = regular_symmetric_hadamard_matrix_with_constant_diagonal(n,sgn(a)*e)
        if list(H.column(0)[1:]).count(1) == k:
            H = -H
        G = Graph((J(n)-I(n)-H+H[0,0]*I(n))/2,loops=False,multiedges=False,format="adjacency_matrix")
        if check:
            assert G.is_strongly_regular(parameters=True) == (v,k,l,mu)
        return G

    if existence:
        return False
    raise ValueError("I do not know how to build a {}-SRG from a RSHCD".format((v,k,l,mu)))

@cached_function
def is_unitary_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Unitary Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_unitary_polar
        sage: t = is_unitary_polar(45, 12, 3, 3); t
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: g = t[0](*t[1:]); g
        Unitary Polar Graph U(4, 2); GQ(4, 2): Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)

        sage: t = is_unitary_polar(5,5,5,5); t

    TESTS:

    All the ``U(n,q)`` appear::

        sage: t = is_unitary_polar(45, 12, 3, 3); t
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: t = is_unitary_polar(165, 36, 3, 9); t
        (<function UnitaryPolarGraph at ...>, 5, 2)
        sage: t = is_unitary_polar(693, 180, 51, 45); t
        (<function UnitaryPolarGraph at ...>, 6, 2)
        sage: t = is_unitary_polar(1105, 80, 15, 5); t
        (<function UnitaryPolarGraph at ...>, 4, 4)
    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q = k/mu
    if q*mu != k or q < 2:
        return
    p,t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2 != 0:
        return
    # at this point we know that we should have U(n,q) for some n and q=p^t, t even
    if r > 0:
        q_pow_d_minus_one = r+1
    else:
        q_pow_d_minus_one = -s-1
    ppp,ttt = is_prime_power(q_pow_d_minus_one, get_data=True)
    d = ttt/t + 1
    if ppp != p or (d-1)*t != ttt:
        return
    t /= 2
    # U(2d+1,q); write q^(1/2) as p^t
    if (v == (q**d - 1)*((q**d)*p**t + 1)/(q - 1)               and
        k == q*(q**(d-1) - 1)*((q**d)/(p**t) + 1)/(q - 1)       and
        l == q*q*(q**(d-2)-1)*((q**(d-1))/(p**t) + 1)/(q - 1) + q - 1):
        from sage.graphs.generators.classical_geometries import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d+1, p**t)

    # U(2d,q);
    if (v == (q**d - 1)*((q**d)/(p**t) + 1)/(q - 1)             and
        k == q*(q**(d-1) - 1)*((q**(d-1))/(p**t) + 1)/(q - 1)   and
        l == q*q*(q**(d-2)-1)*((q**(d-2))/(p**t) + 1)/(q - 1) + q - 1):
        from sage.graphs.generators.classical_geometries import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d, p**t)

@cached_function
def is_unitary_dual_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Unitary Dual Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    This must be the U_5(q) on totally isotropic lines.
    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_unitary_dual_polar
        sage: t = is_unitary_dual_polar(297, 40, 7, 5); t
        (<function UnitaryDualPolarGraph at ...>, 5, 2)
        sage: g = t[0](*t[1:]); g
        Unitary Dual Polar Graph DU(5, 2); GQ(8, 4): Graph on 297 vertices
        sage: g.is_strongly_regular(parameters=True)
        (297, 40, 7, 5)
        sage: t = is_unitary_dual_polar(5,5,5,5); t

    TESTS::

        sage: is_unitary_dual_polar(6832, 270, 26, 10)
        (<function UnitaryDualPolarGraph at ...>, 5, 3)
    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q = mu - 1
    if q < 2:
        return
    p,t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2 != 0:
        return
    if (r < 0 and q != -r - 1) or (s < 0 and q != -s - 1):
       return
    t /= 2
    # we have correct mu, negative eigenvalue, and q=p^(2t)
    if (v == (q**2*p**t + 1)*(q*p**t + 1)  and
        k == q*p**t*(q + 1)                and
        l == k - 1 - q**2*p**t):
        from sage.graphs.generators.classical_geometries import UnitaryDualPolarGraph
        return (UnitaryDualPolarGraph, 5, p**t)

@cached_function
def is_GQqmqp(int v,int k,int l,int mu):
    r"""
    Test whether some `GQ(q-1,q+1)` or `GQ(q+1,q-1)`-graph is `(v,k,\lambda,\mu)`-srg.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_GQqmqp
        sage: t = is_GQqmqp(27,10,1,5); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 3, False)
        sage: g = t[0](*t[1:]); g
        AS(3); GQ(2, 4): Graph on 27 vertices
        sage: t = is_GQqmqp(45,12,3,3); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 3, True)
        sage: g = t[0](*t[1:]); g
        AS(3)*; GQ(4, 2): Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)
        sage: t = is_GQqmqp(16,6,2,2); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 2, True)
        sage: g = t[0](*t[1:]); g
        T2*(O,2)*; GQ(3, 1): Graph on 16 vertices
        sage: g.is_strongly_regular(parameters=True)
        (16, 6, 2, 2)
        sage: t = is_GQqmqp(64,18,2,6); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 4, False)
        sage: g = t[0](*t[1:]); g
        T2*(O,4); GQ(3, 5): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 18, 2, 6)

    TESTS::

        sage: (S,T)=(127,129)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 128, False)
        sage: (S,T)=(129,127)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 128, True)
        sage: (S,T)=(124,126)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 125, False)
        sage: (S,T)=(126,124)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 125, True)
        sage: t = is_GQqmqp(5,5,5,5); t
    """
    # do we have GQ(s,t)? we must have mu=t+1, s=l+1,
    # v=(s+1)(st+1), k=s(t+1)
    S=l+1
    T=mu-1
    q = (S+T)//2
    p, w = is_prime_power(q, get_data=True)
    if (v == (S+1)*(S*T+1)      and
        k == S*(T+1)            and
        q == p**w               and
        (S+T)/2 == q):
        if p % 2 == 0:
            from sage.graphs.generators.classical_geometries\
                    import T2starGeneralizedQuadrangleGraph as F
        else:
            from sage.graphs.generators.classical_geometries\
                    import AhrensSzekeresGeneralizedQuadrangleGraph as F
        if (S,T) == (q-1, q+1):
            return (F, q, False)
        elif (S,T) == (q+1, q-1):
            return (F, q, True)

@cached_function
def is_twograph_descendant_of_srg(int v, int k0, int l, int mu):
    r"""
    Test whether some descendant graph of a s.r.g. is `(v,k_0,\lambda,\mu)`-s.r.g.

    We check whether there can exist `(v+1,k,\lambda^*,\mu^*)`-s.r.g. `G` so
    that ``self`` is a descendant graph of the regular two-graph specified
    by `G`.
    Specifically, we must have that `v+1=2(2k-\lambda^*-\mu^*)`, and
    `k_0=2(k-\mu^*)`, `\lambda=k+\lambda^*-2\mu^*`, `\mu=k-\mu^*`, which give 2
    independent linear conditions, say `k-\mu^*=\mu` and
    `\lambda^*-\mu^*=\lambda-\mu`.  Further, there is a quadratic relation
    `2 k^2-(v+1+4 \mu) k+ 2 v \mu=0`.

    If we can contruct such `G` then we return a function to build a
    `(v,k_0,\lambda,\mu)`-s.r.g.  For more information,
    see 10.3 in http://www.win.tue.nl/~aeb/2WF02/spectra.pdf

    INPUT:

    - ``v,k0,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists and is known, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_twograph_descendant_of_srg
        sage: t = is_twograph_descendant_of_srg(27, 10, 1, 5); t
        (<cyfunction is_twograph_descendant_of_srg.<locals>.la at...
        sage: g = t[0](*t[1:]); g
        descendant of complement(Johnson graph with parameters 8,2) at {5, 7}: Graph on 27 vertices
        sage: g.is_strongly_regular(parameters=True)
        (27, 10, 1, 5)
        sage: t = is_twograph_descendant_of_srg(5,5,5,5); t

    TESTS::

        sage: graphs.strongly_regular_graph(279, 150, 85, 75, existence=True)
        True
        sage: graphs.strongly_regular_graph(279, 150, 85, 75).is_strongly_regular(parameters=True) # optional - gap_packages internet
        (279, 150, 85, 75)
    """
    cdef int b, k, s
    if k0 != 2*mu or v % 2 == 0:
        return
    b = v+1+4*mu
    D = sqrt(b**2-16*v*mu)
    if int(D)==D:
        for kf in [(-D+b)/4, (D+b)/4]:
            k = int(kf)
            if k == kf and \
                strongly_regular_graph(v+1, k, l - 2*mu + k , k - mu,  existence=True):
                def la(vv):
                    from sage.combinat.designs.twographs import twograph_descendant
                    g = strongly_regular_graph(vv, k, l - 2*mu + k)
                    return twograph_descendant(g, g.vertex_iterator().next(), name=True)
                return(la, v+1)
    return

@cached_function
def is_taylor_twograph_srg(int v,int k,int l,int mu):
    r"""
    Test whether some Taylor two-graph SRG is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see §7E of [BvL84]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph
    :func:`TaylorTwographSRG
    <sage.graphs.graph_generators.GraphGenerators.TaylorTwographSRG>` if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_taylor_twograph_srg
        sage: t = is_taylor_twograph_srg(28, 15, 6, 10); t
        (<function TaylorTwographSRG at ...>, 3)
        sage: g = t[0](*t[1:]); g
        Taylor two-graph SRG: Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)
        sage: t = is_taylor_twograph_srg(5,5,5,5); t

    TESTS::

        sage: is_taylor_twograph_srg(730, 369, 168, 205)
        (<function TaylorTwographSRG at ...>, 9)

    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    p,t = is_prime_power(v-1, get_data=True)
    if p**t+1 != v or t % 3 != 0 or p % 2 == 0:
        return
    q = p**(t//3)
    if (k, l, mu) == (q*(q**2+1)/2, (q**2+3)*(q-1)/4, (q**2+1)*(q+1)/4):
        from sage.graphs.generators.classical_geometries import TaylorTwographSRG
        return (TaylorTwographSRG, q)
    return

def is_switch_skewhad(int v, int k, int l, int mu):
    r"""
    Test whether some `switch skewhad^2+*` is `(v,k,\lambda,\mu)`-strongly regular.

    The `switch skewhad^2+*` graphs appear on `Andries Brouwer's database
    <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__ and are built by
    adding an isolated vertex to the complement of
    :func:`~sage.graphs.graph_generators.GraphGenerators.SquaredSkewHadamardMatrixGraph`,
    and a :meth:`Seidel switching <Graph.seidel_switching>` a set of disjoint
    `n`-cocliques.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: graphs.strongly_regular_graph(226, 105, 48, 49)
        switch skewhad^2+*_4: Graph on 226 vertices

    TESTS::

        sage: from sage.graphs.strongly_regular_db import is_switch_skewhad
        sage: t = is_switch_skewhad(5,5,5,5); t

    """
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    from sage.graphs.generators.families import SwitchedSquaredSkewHadamardMatrixGraph
    cdef int n
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    if r<s:
        r,s = s,r
    n = -s // 2
    if  int(r) == 2*n-1         and \
        v == (4*n-1)**2 + 1     and \
        k == (4*n-1)*(2*n-1)    and \
        skew_hadamard_matrix(4*n, existence=True):
            return (SwitchedSquaredSkewHadamardMatrixGraph, n)

def is_switch_OA_srg(int v, int k, int l, int mu):
    r"""
    Test whether some *switch* `OA(k,n)+*` is `(v,k,\lambda,\mu)`-strongly regular.

    The "switch* `OA(k,n)+*` graphs appear on `Andries Brouwer's database
    <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__ and are built by
    adding an isolated vertex to a
    :meth:`~sage.graphs.graph_generators.GraphGenerators.OrthogonalArrayBlockGraph`,
    and a :meth:`Seidel switching <Graph.seidel_switching>` a set of disjoint
    `n`-cocliques.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: graphs.strongly_regular_graph(170, 78, 35, 36) # indirect doctest
        Graph on 170 vertices

    TESTS::

        sage: from sage.graphs.strongly_regular_db import is_switch_OA_srg
        sage: t = is_switch_OA_srg(5,5,5,5); t
        sage: t = is_switch_OA_srg(170, 78, 35, 36);
        sage: t[0](*t[1:]).is_strongly_regular(parameters=True)
        (170, 78, 35, 36)
        sage: t = is_switch_OA_srg(290, 136,  63,  64);
        sage: t[0](*t[1:]).is_strongly_regular(parameters=True)
        (290, 136, 63, 64)
        sage: is_switch_OA_srg(626, 300, 143, 144)
        (<cyfunction is_switch_OA_srg.<locals>.switch_OA_srg at ..., 12, 25)
        sage: is_switch_OA_srg(842, 406, 195, 196)
        (<cyfunction is_switch_OA_srg.<locals>.switch_OA_srg at ..., 14, 29)

    """
    cdef int n_2_p_1 = v
    cdef int n = <int> floor(sqrt(n_2_p_1-1))

    if n*n != n_2_p_1-1: # is it a square?
        return None

    cdef int c = k/n
    if (k % n                            or
        l != c*c-1                       or
        k != 1+(c-1)*(c+1)+(n-c)*(n-c-1) or
        not orthogonal_array(c+1,n,existence=True,resolvable=True)):
        return None

    def switch_OA_srg(c,n):
        from itertools import izip
        OA = map(tuple,orthogonal_array(c+1,n,resolvable=True))
        g = Graph([OA,lambda x,y: any(xx==yy for xx,yy in izip(x,y))],loops=False)
        g.add_vertex(0)
        g.seidel_switching(OA[:c*n])
        return g

    return (switch_OA_srg,c,n)

cdef eigenvalues(int v,int k,int l,int mu):
    r"""
    Return the eigenvalues of a (v,k,l,mu)-strongly regular graph.

    If the set of parameters is not feasible, or if they correspond to a
    conference graph, the function returns ``(None,None)``. Otherwise
    it returns the pair [r,s] of eigenvalues, satisfying r>s.

    INPUT:

    - ``v,k,l,mu`` (integers)

    """
    # See 1.3.1 of [Distance-regular graphs]
    b = (mu-l)
    c = (mu-k)
    D = b**2-4*c
    if not is_square(D):
        return [None,None]
    return [(-b+sqrt(D))/2.0,
            (-b-sqrt(D))/2.0]

def eigenmatrix(int v,int k,int l,int mu):
    r"""
    Return the 1st eigenmatrix of a `(v,k,l,mu)`-strongly regular graph.

    The adjacency matrix `A` of an s.r.g. commutes with the adacency matrix
    `A'=J-A-I` of its complement (here `J` is all-1 matrix, and `I` the identity
    matrix).  Thus, they can be simultaneously diagonalized and so `A` and `A'`
    share eigenspaces.

    The eigenvalues of `J` are `v` with multiplicity 1, and 0 with multiplicity
    `v-1`. Thus the eigenvalue of `A'` corresponding to the 1-dimension
    `k`-eigenspace of `A` is `v-k-1`.  Respecively, the eigenvalues of `A'`
    corresponding to `t`-eigenspace of `A`, with `t` unequal to `k`, equals
    `-t-1`. The 1st eigenmatrix `P` of the C-algebra `C[A]` generated by `A`
    encodes this eigenvalue information in its three colums; the 2nd (resp. 3rd)
    column contains distinct eigenvalues of `A` (resp. of `A'`), and the 1st
    column contains the corresponding eigenvalues of `I`. The matrix `vP^{-1}`
    is called the 2nd eigenvalue matrix of `C[A]`.

    The most interesting feature of `vP^{-1}` is that it is the 1st eigenmatrix
    of the dual of `C[A]` if the dual is generated by the adjacency matrix of a
    strongly regular graph. See [BH12]_ and [BI84]_ for details.

    If the set of parameters is not feasible, or if they correspond to a
    conference graph, the function returns ``None``. Its output is stable, assuming
    that the eigenvalues r,s used satisfy r>s; this holds for the current
    implementation of eigenvalues().

    INPUT:

    - ``v,k,l,mu`` (integers)

    EXAMPLES:

    Petersen's graph's C-algebra does not have a dual coming from an s.r.g.::

        sage: from sage.graphs.strongly_regular_db import eigenmatrix
        sage: P=eigenmatrix(10,3,0,1); P
        [ 1  3  6]
        [ 1  1 -2]
        [ 1 -2  1]
        sage: 10*P^-1
        [   1    5    4]
        [   1  5/3 -8/3]
        [   1 -5/3  2/3]

    The line graph of `K_{3,3}` is self-dual::

        sage: P=eigenmatrix(9,4,1,2); P
        [ 1  4  4]
        [ 1  1 -2]
        [ 1 -2  1]
        sage: 9*P^-1
        [ 1  4  4]
        [ 1  1 -2]
        [ 1 -2  1]

    A strongly regular graph with a non-isomorphic dual coming from another
    strongly regular graph::

        sage: graphs.strongly_regular_graph(243,220,199,200, existence=True)
        True
        sage: graphs.strongly_regular_graph(243,110,37,60, existence=True)
        True
        sage: P=eigenmatrix(243,220,199,200); P
        [  1 220  22]
        [  1   4  -5]
        [  1  -5   4]
        sage: 243*P^-1
        [  1 110 132]
        [  1   2  -3]
        [  1 -25  24]
        sage: 243*P^-1==eigenmatrix(243,110,37,60)
        True

    TESTS::

        sage: eigenmatrix(5,5,5,-5)

    REFERENCE:

    .. [BI84] Eiichi Bannai, Tatsuro Ito,
      Algebraic Combinatorics I: Association Schemes,
      Benjamin/Cummings, 1984
    """
    from sage.rings.integer_ring import ZZ
    r,s = eigenvalues(v,k,l,mu)
    if r is not None:
        return Matrix(ZZ, [[1,k,v-k-1],[1,r,-r-1],[1,s,-s-1]])

cpdef latin_squares_graph_parameters(int v,int k, int l,int mu):
    r"""
    Check whether (v,k,l,mu)-strongly regular graph has parameters of an `L_g(n)` s.r.g.

    Also known as pseudo-OA(n,g) case, i.e. s.r.g. with parameters of an OA(n,g)-graph.
    Return g and n, if they exist. See Sect. 9.1 of [BH12]_ for details.

    INPUT:

    - ``v,k,l,mu`` -- (integers) parameters of the graph

    OUTPUT:

    - ``(g, n)`` -- parameters of an `L_g(n)` graph, or `None`

    TESTS::

        sage: from sage.graphs.strongly_regular_db import latin_squares_graph_parameters
        sage: latin_squares_graph_parameters(9,4,1,2)
        (2, 3)
        sage: latin_squares_graph_parameters(5,4,1,2)
    """
    cdef int g, n
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    if r < s:
        r, s = s, r
    g = -s
    n = r+g
    if v==n**2 and k==g*(n-1) and l==(g-1)*(g-2)+n-2 and mu==g*(g-1):
        return g, n
    return

def _H_3_cayley_graph(L):
    r"""
    return the `L`-Cayley graph of the group `H_3` from Prop. 12 in [JK03]_.

    INPUT:

    - the list of words for the generating set in the format ["abc",...,"xyz"] for
      a,b,...,z being integers between 0 and 4.

    TESTS::

        sage: from sage.graphs.strongly_regular_db import _H_3_cayley_graph
        sage: _H_3_cayley_graph(["100","110","130","140","200","230","240","300"])
        Graph on 100 vertices
    """
    from sage.groups.free_group import FreeGroup
    from sage.groups.finitely_presented import FinitelyPresentedGroup
    G = FreeGroup('x,y,z')
    x,y,z = G.gens()
    rels = (x**5,y**5,z**4,x*y*x**(-1)*y**(-1),z*x*z**(-1)*x**(-2),z*y*z**(-1)*y**(-2))
    G = FinitelyPresentedGroup(G,rels)
    x,y,z = G.gens()
    H = G.as_permutation_group()
    L = map(lambda x:map(int,x),L)
    x,y,z=(H.gen(0),H.gen(1),H.gen(2))
    L = [H(x**xx*y**yy*z**zz) for xx,yy,zz in L]
    return Graph(H.cayley_graph(generators=L, simple=True))

def SRG_100_44_18_20():
    r"""
    Return a `(100, 44, 18, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Delta_1`
    with group `H_3` presented in Table 8.1 of [JK03]_

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_100_44_18_20
        sage: G = SRG_100_44_18_20()                 # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (100, 44, 18, 20)

    REFERENCES:

    .. [JK03] \L. K. Jørgensen, M. Klin, M.,
      Switching of edges in strongly regular graphs.
      I. A family of partial difference sets on 100 vertices,
      Electronic Journal of Combinatorics 10(1), 2003.
    """
    return _H_3_cayley_graph(["100","110","130","140","200","230","240","300",
             "310","320","400","410","420","440","041","111","221","231","241",
             "321","331","401","421","441","002","042","112","122","142","212",
             "232","242","322","342","033","113","143","223","303","333","343",
             "413","433","443"])

def SRG_100_45_20_20():
    r"""
    Return a `(100, 45, 20, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Gamma_3`
    with group `H_3` presented in Table 8.1 of [JK03]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_100_45_20_20
        sage: G = SRG_100_45_20_20()              # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (100, 45, 20, 20)
    """
    return _H_3_cayley_graph(["120","140","200","210","201","401","411","321",
             "002","012","022","042","303","403","013","413","240","031","102",
             "323","300","231","132","133","310","141","142","233","340","241",
             "202","333","410","341","222","433","430","441","242","302","312",
             "322","332","442","143"])


def SRG_105_32_4_12():
    r"""
    Return a `(105, 32, 4, 12)`-strongly regular graph.

    The vertices are the flags of the projective plane of order 4. Two flags
    `(a,A)` and `(b,B)` are adjacent if the point `a` is on the line `B` or
    the point `b` is  on the line `A`, and `a \neq b`, `A \neq B`. See
    Theorem 2.7 in [GS70]_, and [Co06]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_105_32_4_12
        sage: G = SRG_105_32_4_12(); G
        Graph on 105 vertices
        sage: G.is_strongly_regular(parameters=True)
        (105, 32, 4, 12)

    REFERENCES:

    .. [GS70] \J.-M. Goethals and J. J. Seidel,
       Strongly regular graphs derived from combinatorial designs,
       Can. J. Math. 22 (1970) 597-614.
       http://dx.doi.org/10.4153/CJM-1970-067-9

    .. [Co06] \K. Coolsaet,
       The uniqueness of the strongly regular graph srg(105,32,4,12),
       Bull. Belg. Math. Soc. 12(2006), 707-718.
       http://projecteuclid.org/euclid.bbms/1136902608
    """
    from sage.combinat.designs.block_design import ProjectiveGeometryDesign
    P = ProjectiveGeometryDesign(2,1,GF(4,'a'))
    IG = P.incidence_graph().line_graph()
    a = IG.automorphism_group()
    h = a.stabilizer(a.domain()[0])
    o = filter(lambda x: len(x)==32, h.orbits())[0][0]
    e = a.orbit((a.domain()[0],o),action="OnSets")
    G = Graph()
    G.add_edges(e)
    return G

def SRG_120_77_52_44():
    r"""
    Return a `(120,77,52,44)`-strongly regular graph.

    To build this graph, we first build a `2-(21,7,12)` design, by removing two
    points from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_120_77_52_44
        sage: G = SRG_120_77_52_44()                 # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (120, 77, 52, 44)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x and 21 not in x])
    return H.intersection_graph(3)

def SRG_144_39_6_12():
    r"""
    Return a `(144,39,6,12)`-strongly regular graph.

    This graph is obtained as an orbit of length 2808 on sets of cardinality 2
    (among 2 such orbits) of the group `PGL_3(3)` acting on the (right) cosets of
    a subgroup of order 39.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_144_39_6_12
        sage: G = SRG_144_39_6_12()
        sage: G.is_strongly_regular(parameters=True)
        (144, 39, 6, 12)
    """

    from sage.libs.gap.libgap import libgap
    g=libgap.ProjectiveGeneralLinearGroup(3,3)
    ns=g.Normalizer(g.SylowSubgroup(13))
    G=g.Action(g.RightCosets(ns),libgap.OnRight)
    H=G.Stabilizer(1)
    for o in filter(lambda x: len(x)==39, H.Orbits()):
        h = Graph()
        h.add_edges(G.Orbit([1,o[0]],libgap.OnSets))
        if h.is_strongly_regular():
            h.relabel()
            return h

def SRG_176_49_12_14():
    r"""
    Return a `(176,49,12,14)`-strongly regular graph.

    This graph is built from the symmetric Higman-Sims design. In
    [BrouwerPolarities82]_, it is explained that there exists an involution
    `\sigma` exchanging the points and blocks of the Higman-Sims design, such
    that each point is mapped on a block that contains it (i.e. `\sigma` is a
    'polarity with all universal points'). The graph is then built by making two
    vertices `u,v` adjacent whenever `v\in \sigma(u)`.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_176_49_12_14
        sage: G = SRG_176_49_12_14()                 # optional - gap_packages # long time
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages # long time
        (176, 49, 12, 14)

    REFERENCE:

    .. [BrouwerPolarities82] \A. Brouwer,
       Polarities of G. Higman's symmetric design and a strongly regular graph on 176 vertices,
       Aequationes mathematicae 25, no. 1 (1982): 77-82.
    """
    from sage.combinat.designs.database import HigmanSimsDesign
    d = HigmanSimsDesign()
    g = d.incidence_graph(labels=True)
    ag=g.automorphism_group().conjugacy_classes_representatives()

    # Looking for an involution that maps a point of the design to one of the
    # blocks that contains it. It is called a polarity with only absolute
    # points in
    for aut in ag:
        try:
            0 in aut(0)
        except TypeError:
            continue
        if (aut.order() == 2 and
            all(i in aut(i) for i in d.ground_set())):
            g = Graph()
            g.add_edges((u,v) for u in d.ground_set() for v in aut(u))
            return g

def SRG_176_105_68_54():
    r"""
    Return a `(176, 105, 68, 54)`-strongly regular graph.

    To build this graph, we first build a `2-(22,7,16)` design, by removing one
    point from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3. Known as S.7 in [Hu75]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_176_105_68_54
        sage: G = SRG_176_105_68_54()                # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (176, 105, 68, 54)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x])
    return H.intersection_graph(3)

def SRG_210_99_48_45():
    r"""
    Return a strongly regular graph with parameters `(210, 99, 48, 45)`

    This graph is from Example 4.2 in [KPRWZ10]_. One considers the action of
    the symmetric group `S_7` on the 210 digraphs isomorphic to the
    disjoint union of `K_1` and the circulant 6-vertex digraph
    ``digraphs.Circulant(6,[1,4])``. It has 16 orbitals; the package [COCO]_
    found a megring of them, explicitly described in [KPRWZ10]_, resulting in
    this graph.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_210_99_48_45
        sage: g=SRG_210_99_48_45()
        sage: g.is_strongly_regular(parameters=True)
        (210, 99, 48, 45)

    REFERENCES:

    .. [KPRWZ10] \M. H. Klin, C. Pech, S. Reichard, A. Woldar, M. Zvi-Av,
       Examples of computer experimentation in algebraic combinatorics,
       ARS MATHEMATICA CONTEMPORANEA 3 (2010) 237–258
       http://amc-journal.eu/index.php/amc/article/viewFile/119/118

    .. [COCO] \I. A. Faradjev and M. H. Klin,
       Computer package for computations with coherent configurations,
       Proc. ISSAC-91, ACM Press, Bonn, 1991, pages 219–223;
       code, by I.A.Faradjev (with contributions by A.E.Brouwer, D.V.Pasechnik)
       https://github.com/dimpase/coco

    """
    from sage.libs.gap.libgap import libgap
    from sage.combinat.permutation import Permutation
    def ekg(g0): # return arcs of the Cayley digraph of <g> on {g,g^4}
        g = Permutation(g0)
        return libgap.Set(map(lambda x: (x,g(x)), range(1,8))\
                        + map(lambda x: (x,g(g(g(g(x))))), range(1,8)))

    kd=map(ekg,
        [(7, 1, 2, 3, 4, 5), (7, 1, 3, 4, 5, 6),
        (7, 3, 4, 5, 6, 2), (7, 1, 4, 3, 5, 6),
        (7, 3, 1, 4, 5, 6), (7, 2, 4, 3, 5, 6),
        (7, 3, 2, 4, 5, 1), (7, 2, 4, 3, 5, 1)])
    s=libgap.SymmetricGroup(7)
    O=s.Orbit(kd[0],libgap.OnSetsTuples)
    sa=s.Action(O,libgap.OnSetsTuples)
    G=Graph()
    for g in kd[1:]:
        G.add_edges(libgap.Orbit(sa,[libgap.Position(O,kd[0]),\
                                     libgap.Position(O,g)],libgap.OnSets))
    return G

def SRG_243_110_37_60():
    r"""
    Return a `(243, 110, 37, 60)`-strongly regular graph.

    Consider the orthogonal complement of the
    :func:`~sage.coding.code_constructions.TernaryGolayCode`, which has 243
    words. On them we define a graph, in which two words are adjacent
    whenever their Hamming distance is 9. This construction appears in [GS75]_.

    .. NOTE::

        A strongly regular graph with the same parameters is also obtained from
        the database of 2-weight codes.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_243_110_37_60
        sage: G = SRG_243_110_37_60()
        sage: G.is_strongly_regular(parameters=True)
        (243, 110, 37, 60)

    REFERENCE:

    .. [GS75] \J.M. Goethals, and J. J. Seidel,
       The regular two-graph on 276 vertices,
       Discrete Mathematics 12, no. 2 (1975): 143-158.
       http://dx.doi.org/10.1016/0012-365X(75)90029-1
    """
    from sage.coding.code_constructions import TernaryGolayCode
    M = TernaryGolayCode().generator_matrix()
    V = list(M.right_kernel())
    return Graph([range(len(V)), lambda x,y:(V[x]-V[y]).hamming_weight() == 9 ])

def SRG_253_140_87_65():
    r"""
    Return a `(253, 140, 87, 65)`-strongly regular graph.

    To build this graph, we first build the
    :func:`~sage.combinat.designs.block_design.WittDesign` on 23 points which is
    a `2-(23,7,21)` design. We then build the intersection graph of blocks with
    intersection size 3. Known as S.6 in [Hu75]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_253_140_87_65
        sage: G = SRG_253_140_87_65()                # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (253, 140, 87, 65)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    return W.intersection_graph(3)

def SRG_196_91_42_42():
    r"""
    Return a `(196,91,42,42)`-strongly regular graph.

    This strongly regular graph is built following the construction provided in
    Corollary 8.2.27 of [IS06]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_196_91_42_42
        sage: G = SRG_196_91_42_42()
        sage: G.is_strongly_regular(parameters=True)
        (196, 91, 42, 42)

    REFERENCE:

    .. [IS06] \Y.J. Ionin, S. Shrikhande,
      Combinatorics of symmetric designs.
      Cambridge University Press, 2006.
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.graphs.generators.intersection import IntersectionGraph
    k = 7
    G = IntegerModRing(91)
    A = map(G,{0, 10, 27, 28, 31, 43, 50})
    B = map(G,{0, 11, 20, 25, 49, 55, 57})
    H = map(G,[13*i for i in range(k)])
    U = map(frozenset,[[x+z for x in A] for z in G])
    V = map(frozenset,[[x+z for x in B] for z in G])
    W = map(frozenset,[[x+z for x in H] for z in G])
    G = IntersectionGraph(U+V+W)

    G.seidel_switching(U)

    G.add_edges((-1,x) for x in U)
    G.relabel()
    return G

def SRG_220_84_38_28():
    r"""
    Return a `(280, 135, 70, 60)`-strongly regular graph.

    This graph is obtained from the
    :meth:`~IncidenceStructure.intersection_graph` of a
    :func:`~sage.combinat.designs.database.BIBD_45_9_8`. This construction
    appears in VII.11.2 from [DesignHandbook]_

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_220_84_38_28
        sage: g=SRG_220_84_38_28()
        sage: g.is_strongly_regular(parameters=True)
        (220, 84, 38, 28)
    """
    from sage.combinat.designs.database import BIBD_45_9_8
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    G = IncidenceStructure(BIBD_45_9_8()).intersection_graph(3)
    G.relabel()
    return G

def SRG_276_140_58_84():
    r"""
    Return a `(276, 140, 58, 84)`-strongly regular graph.

    The graph is built from from
    :meth:`~sage.graphs.graph_generators.GraphGenerators.McLaughlinGraph`, with
    an added isolated vertex. We then perform a
    :meth:`~Graph.seidel_switching` on a set of 28 disjoint 5-cliques, which
    exist by cf. [HT96]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_276_140_58_84
        sage: g=SRG_276_140_58_84()                  # long time # optional - gap_packages
        sage: g.is_strongly_regular(parameters=True) # long time # optional - gap_packages
        (276, 140, 58, 84)

    REFERENCE:

    .. [HT96] \W. H. Haemers and V. D. Tonchev,
      Spreads in strongly regular graphs,
      Designs, Codes and Cryptography 8 (1996) 145-157.
    """
    from sage.graphs.generators.smallgraphs import McLaughlinGraph
    g = McLaughlinGraph()
    C = [[ 0,  72,  87, 131, 136], [ 1,  35,  61, 102, 168], [ 2,  32,  97, 125, 197], [ 3,  22,  96, 103, 202],
         [ 4,  46,  74, 158, 229], [ 5,  83,  93, 242, 261], [ 6,  26,  81, 147, 176], [ 7,  42,  63, 119, 263],
         [ 8,  49,  64, 165, 227], [ 9,  70,  85, 208, 273], [10,  73,  92, 230, 268], [11,  54,  95, 184, 269],
         [12,  55,  62, 185, 205], [13,  51,  65, 162, 254], [14,  78,  88, 231, 274], [15,  40,  59, 117, 252],
         [16,  24,  71, 137, 171], [17,  39,  43, 132, 163], [18,  57,  79, 175, 271], [19,  68,  80, 217, 244],
         [20,  75,  98, 239, 267], [21,  33,  56, 113, 240], [23, 127, 152, 164, 172], [25, 101, 128, 183, 264],
         [27, 129, 154, 160, 201], [28, 126, 144, 161, 228], [29, 100, 133, 204, 266], [30, 108, 146, 200, 219]]
    g.add_vertex(-1)
    g.seidel_switching(sum(C,[]))
    g.relabel()
    g.name('')
    return g

def SRG_280_135_70_60():
    r"""
    Return a strongly regular graph with parameters `(280, 135, 70, 60)`.

    This graph is built from the action of `J_2` on a `3.PGL(2,9)` subgroup it
    contains.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_280_135_70_60
        sage: g=SRG_280_135_70_60()                  # long time # optional - gap_packages
        sage: g.is_strongly_regular(parameters=True) # long time # optional - gap_packages
        (280, 135, 70, 60)
    """
    from sage.interfaces.gap import gap
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.graphs.graph import Graph

    gap.load_package("AtlasRep")

    # A representation of J2 acting on a 3.PGL(2,9) it contains.
    J2    = PermutationGroup(gap('AtlasGenerators("J2",2).generators'))
    edges = J2.orbit((1,2),"OnSets")
    g     = Graph()
    g.add_edges(edges)
    g.relabel()
    return g

def SRG_280_117_44_52():
    r"""
    Return a strongly regular graph with parameters `(280, 117, 44, 52)`.

    This graph is built according to a very pretty construction of Mathon and
    Rosa [MR85]_:

        The vertices of the graph `G` are all partitions of a set of 9 elements
        into `\{\{a,b,c\},\{d,e,f\},\{g,h,i\}\}`. The cross-intersection of two
        such partitions `P=\{P_1,P_2,P_3\}` and `P'=\{P'_1,P'_2,P'_3\}` being
        defined as `\{P_i \cap P'_j: 1\leq i,j\leq 3\}`, two vertices of `G` are
        set to be adjacent if the cross-intersection of their respective
        partitions does not contain exactly 7 nonempty sets.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_280_117_44_52
        sage: g=SRG_280_117_44_52()
        sage: g.is_strongly_regular(parameters=True)
        (280, 117, 44, 52)

    REFERENCE:

    .. [MR85] \R. Mathon and A. Rosa,
       A new strongly regular graph,
       Journal of Combinatorial Theory, Series A 38, no. 1 (1985): 84-86.
       http://dx.doi.org/10.1016/0097-3165(85)90025-1
    """
    from sage.graphs.hypergraph_generators import hypergraphs

    # V is the set of partions {{a,b,c},{d,e,f},{g,h,i}} of {0,...,8}
    H = hypergraphs.CompleteUniform(9,3)
    g = H.intersection_graph()
    V = g.complement().cliques_maximal()
    V = map(frozenset,V)

    # G is the graph defined on V in which two vertices are adjacent when they
    # corresponding partitions cross-intersect on 7 nonempty sets
    G = Graph([V, lambda x,y:
               sum(any(xxx in yy for xxx in xx) for xx in x for yy in y) != 7],
              loops=False)
    return G

def strongly_regular_from_two_weight_code(L):
    r"""
    Return a strongly regular graph from a two-weight code.

    A code is said to be a *two-weight* code the weight of its nonzero codewords
    (i.e. their number of nonzero coordinates) can only be one of two integer
    values `w_1,w_2`. It is said to be *projective* if the minimum weight of the
    dual code is `\geq 3`. A strongly regular graph can be built from a
    two-weight projective code with weights `w_1,w_2` (assuming `w_1<w_2`) by
    adding an edge between any two codewords whose difference has weight
    `w_1`. For more information, see [vLintSchrijver81]_ or [Delsarte72]_.

    INPUT:

    - ``L`` -- a two-weight linear code, or its generating matrix.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_from_two_weight_code
        sage: x=("100022021001111",
        ....:    "010011211122000",
        ....:    "001021112100011",
        ....:    "000110120222220")
        sage: M = Matrix(GF(3),[list(l) for l in x])
        sage: G = strongly_regular_from_two_weight_code(LinearCode(M))
        sage: G.is_strongly_regular(parameters=True)
        (81, 50, 31, 30)

    REFERENCES:

    .. [vLintSchrijver81] \J. H. van Lint, and A. Schrijver (1981),
      Construction of strongly regular graphs, two-weight codes and
      partial geometries by finite fields,
      Combinatorica, 1(1), 63-73.

    .. [Delsarte72] Ph. Delsarte,
      Weights of linear codes and strongly regular normed spaces,
      Discrete Mathematics (1972), Volume 3, Issue 1, Pages 47-64,
      http://dx.doi.org/10.1016/0012-365X(72)90024-6.

    """
    from sage.matrix.matrix import is_Matrix
    if is_Matrix(L):
        L = LinearCode(L)
    V = map(tuple,list(L))
    w1, w2 = sorted(set(sum(map(bool,x)) for x in V).difference([0]))
    G = Graph([V,lambda u,v: sum(uu!=vv for uu,vv in zip(u,v)) == w1])
    G.relabel()
    return G

def SRG_416_100_36_20():
    r"""
    Return a `(416,100,36,20)`-strongly regular graph.

    This graph is obtained as an orbit on sets of cardinality 2
    (among 2 that exists) of the group `G_2(4)`.
    This graph is isomorphic to the subgraph of the from :meth:`Suzuki Graph
    <sage.graphs.graph_generators.GraphGenerators.SuzukiGraph>` induced on
    the neighbors of a vertex. Known as S.14 in [Hu75]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_416_100_36_20
        sage: g = SRG_416_100_36_20()                # optional - gap_packages # long time
        sage: g.is_strongly_regular(parameters=True) # optional - gap_packages # long time
        (416, 100, 36, 20)
    """
    from sage.libs.gap.libgap import libgap
    libgap.LoadPackage("AtlasRep")
    g=libgap.AtlasGroup("G2(4)",libgap.NrMovedPoints,416)
    h = Graph()
    h.add_edges(g.Orbit([1,5],libgap.OnSets))
    h.relabel()
    return h

def SRG_560_208_72_80():
    r"""
    Return a `(560,208,72,80)`-strongly regular graph

    This graph is obtained as the union of 4 orbits of sets of cardinality 2
    (among the 13 that exists) of the group `Sz(8)`.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_560_208_72_80
        sage: g = SRG_560_208_72_80()                # optional - database_gap # not tested (~2s)
        sage: g.is_strongly_regular(parameters=True) # optional - database_gap # not tested (~2s)
        (560, 208, 72, 80)
    """
    from sage.libs.gap.libgap import libgap
    libgap.LoadPackage("AtlasRep")
    g=libgap.AtlasGroup("Sz8",libgap.NrMovedPoints,560)

    h = Graph()
    h.add_edges(g.Orbit([1,2],libgap.OnSets))
    h.add_edges(g.Orbit([1,4],libgap.OnSets))
    h.add_edges(g.Orbit([1,8],libgap.OnSets))
    h.add_edges(g.Orbit([1,27],libgap.OnSets))
    h.relabel()
    return h

def strongly_regular_from_two_intersection_set(M):
    r"""
    Return a strongly regular graph from a 2-intersection set.

    A set of points in the projective geometry `PG(k,q)` is said to be a
    2-intersection set if it intersects every hyperplane in either `h_1` or
    `h_2` points, where `h_1,h_2\in \\NN`.

    From a 2-intersection set `S` can be defined a strongly-regular graph in the
    following way:

    - Place the points of `S` on a hyperplane `H` in `PG(k+1,q)`

    - Define the graph `G` on all points of `PG(k+1,q)\backslash H`

    - Make two points of `V(G)=PG(k+1,q)\backslash H` adjacent if the line going
      through them intersects `S`

    For more information, see e.g. [CDB13]_ where this explanation has been
    taken from.

    INPUT:

    - `M` -- a `|S| \times k` matrix with entries in `F_q` representing the points of
      the 2-intersection set. We assume that the first non-zero entry of each row is
      equal to `1`, that is, they give points in homogeneous coordinates.

    The implementation does not check that `S` is actually a 2-intersection set.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_from_two_intersection_set
        sage: S=Matrix([(0,0,1),(0,1,0)] + map(lambda x: (1,x^2,x), GF(4,'b')))
        sage: g=strongly_regular_from_two_intersection_set(S)
        sage: g.is_strongly_regular(parameters=True)
        (64, 18, 2, 6)

    REFERENCES:

    .. [CDB13] \I. Cardinali and B. De Bruyn,
      Spin-embeddings, two-intersection sets and two-weight codes,
      Ars Comb. 109 (2013): 309-319.
      https://biblio.ugent.be/publication/4241842/file/4241845.pdf
    """
    from itertools import product, izip
    K = M.base_ring()
    k = M.ncols()
    g = Graph()

    M = [list(p) for p in M]

    # For every point in F_q^{k+1} not on the hyperplane of M
    for u in [tuple(x) for x in product(K,repeat=k)]:
        # For every v point of M
        for v in M:
            # u is adjacent with all vertices on a uv line.
            g.add_edges([[u,tuple([u[i]+qq*v[i] for i in range(k)])] \
                                            for qq in K if not qq==K.zero()])
    g.relabel()
    return g

def SRG_120_63_30_36():
    r"""
    Return a `(120,63,30,36)`-strongly regular graph

    It is the distance-2 graph of :meth:`JohnsonGraph(10,3)
    <sage.graphs.graph_generators.GraphGenerators.JohnsonGraph>`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_120_63_30_36
        sage: G =  SRG_120_63_30_36()
        sage: G.is_strongly_regular(parameters=True)
        (120, 63, 30, 36)
    """
    from sage.graphs.generators.families import JohnsonGraph
    return JohnsonGraph(10,3).distance_graph([2])

def SRG_126_25_8_4():
    r"""
    Return a `(126,25,8,4)`-strongly regular graph

    It is the distance-(1 or 4) graph of :meth:`JohnsonGraph(9,4)
    <sage.graphs.graph_generators.GraphGenerators.JohnsonGraph>`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_126_25_8_4
        sage: G =  SRG_126_25_8_4()
        sage: G.is_strongly_regular(parameters=True)
        (126, 25, 8, 4)
    """
    from sage.graphs.generators.families import JohnsonGraph
    return JohnsonGraph(9,4).distance_graph([1,4])

def SRG_175_72_20_36():
    r"""
    Return a `(175,72,20,36)`-strongly regular graph

    This graph is obtained from the line graph of
    :meth:`~sage.graphs.graph_generators.GraphGenerators.HoffmanSingletonGraph`. Setting
    two vertices to be adjacent if their distance in the line graph is exactly
    2 yields the graph. For more information, see 10.B.(iv) in [BvL84]_ and
    http://www.win.tue.nl/~aeb/graphs/McL.html.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_175_72_20_36
        sage: G = SRG_175_72_20_36()
        sage: G.is_strongly_regular(parameters=True)
        (175, 72, 20, 36)
    """
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    return HoffmanSingletonGraph().line_graph().distance_graph([2])

def SRG_176_90_38_54():
    r"""
    Return a `(176,90,38,54)`-strongly regular graph

    This graph is obtained from
    :func:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`
    by attaching a isolated vertex and doing Seidel switching
    with respect to disjoint union of 18 maximum cliques, following
    a construction by W.Haemers given in Sect.10.B.(vi) of [BvL84]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_176_90_38_54
        sage: G = SRG_176_90_38_54()
        sage: G.is_strongly_regular(parameters=True)
        (176, 90, 38, 54)
    """
    from sage.graphs.generators.basic import CompleteGraph
    from sage.misc.flatten import flatten
    g = SRG_175_72_20_36()
    g.relabel()
    # c=filter(lambda x: len(x)==5, g.cliques_maximal())
    # r=flatten(Hypergraph(c).packing()[:18]) # takes 3s, so we put the answer here
    r=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,28,29,32,\
       38,39,41,42,43,47,49,50,51,52,53,55,57,61,63,65,67,69,72,75,77,79,81,84,87,88,\
       89,92,95,96,97,99,101,102,104,105,107,112,114,117,118,123,125,129,132,139,140,\
       141,144,146,147,153,154,162,165,166,167,170,172,173,174]
    j=g.disjoint_union(CompleteGraph(1))
    j.relabel()
    j.seidel_switching(r)
    return j

def SRG_630_85_20_10():
    r"""
    Return a `(630,85,20,10)`-strongly regular graph

    This graph is the line graph of `pg(5,18,2)`; its point graph is
    :func:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`.
    One selects a subset of 630 maximum cliques in the latter following
    a construction by W.Haemers given in Sect.10.B.(v) of [BvL84]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_630_85_20_10
        sage: G = SRG_630_85_20_10()                    # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (630, 85, 20, 10)
    """
    from sage.graphs.generators.intersection import IntersectionGraph
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    hs = HoffmanSingletonGraph()
    P = range(5)+range(30,35)          # a Petersen in hs
    mc = [0, 1, 5, 6, 12, 13, 16, 17, 22, 23, 29, 33, 39, 42, 47]
    assert(hs.subgraph(mc).is_regular(k=0)) # a maximum coclique
    assert(hs.subgraph(P).is_regular(k=3))
    h = hs.automorphism_group().stabilizer(mc,action="OnSets")
    l = h.orbit(tuple(map(lambda x: (x[0],x[1]), hs.subgraph(P).matching())),"OnSetsSets")
    return IntersectionGraph(l)

def SRG_126_50_13_24():
    r"""
    Return a `(126,50,13,24)`-strongly regular graph

    This graph is a subgraph of
    :meth:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`.
    This construction, due to Goethals, is given in §10B.(vii) of [BvL84]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_126_50_13_24
        sage: G = SRG_126_50_13_24()
        sage: G.is_strongly_regular(parameters=True)
        (126, 50, 13, 24)
    """
    from sage.graphs.strongly_regular_db import SRG_175_72_20_36
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    hs = HoffmanSingletonGraph()
    s = set(hs.vertices()).difference(hs.neighbors(0)+[0])
    return SRG_175_72_20_36().subgraph(hs.edge_boundary(s,s))



def SRG_1288_792_476_504():
    r"""
    Return a `(1288, 792, 476, 504)`-strongly regular graph.

    This graph is built on the words of weight 12 in the
    :func:`~sage.coding.code_constructions.BinaryGolayCode`. Two of them are
    then made adjacent if their symmetric difference has weight 12 (cf
    [BvE92]_).

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_1288_792_476_504
        sage: G = SRG_1288_792_476_504()             # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (1288, 792, 476, 504)

    REFERENCE:

    .. [BvE92] \A. Brouwer and C. Van Eijl,
      On the p-Rank of the Adjacency Matrices of Strongly Regular Graphs
      Journal of Algebraic Combinatorics (1992), vol.1, n.4, pp329-346,
      http://dx.doi.org/10.1023/A%3A1022438616684
    """
    from sage.coding.code_constructions import BinaryGolayCode
    C = BinaryGolayCode()
    C = [[i for i,v in enumerate(c) if v]
         for c in C]
    C = [s for s in C if len(s) == 12]
    G = Graph([map(frozenset,C),
               lambda x,y:len(x.symmetric_difference(y))==12])
    G.relabel()
    return G

cdef bint seems_feasible(int v, int k, int l, int mu):
    r"""
    Tests is the set of parameters seems feasible

    INPUT:

    - ``v,k,l,mu`` (integers)
    """
    cdef int r,s,f,g
    cdef uint_fast32_t tmp[2]

    if (v<0 or k<=0 or l<0 or mu<0 or
        k>=v-1 or l>=k or mu>k or
        v-2*k+mu-2 < 0 or # lambda of complement graph >=0
        v-2*k+l    < 0 or # μ of complement graph >= 0
        mu*(v-k-1) != k*(k-l-1)):
        return False

    if mu == k: # complete multipartite graph
        r = v/(v-k) # number of parts (of size v-k each)
        return (l == (v-k)*(r-2) and v == r*(v-k))

    if mu == 0: # the complement of a complete multipartite graph
        r = v/(k+1) # number of parts (of size k+1 each)
        return (l == k-1 and v == r*(k+1))

    # Conference graphs. Only possible if 'v' is a sum of two squares (3.A of
    # [BvL84]
    if (v-1)*(mu-l)-2*k == 0:
        return two_squares_c(v,tmp)

    rr,ss = eigenvalues(v,k,l,mu)
    if rr is None:
        return False
    r,s = rr,ss

    # p.87 of [BvL84]_
    # "Integrality condition"
    if ((s+1)*(k-s)*k) % (mu*(s-r)) or ((r+1)*(k-r)*k) % (mu*(s-r)):
        return False

    # Theorem 21.3 of [WilsonACourse] or
    # 3.B of [BvL84]
    # (Krein conditions)
    if ((r+1)*(k+r+2*r*s) > (k+r)*(s+1)**2 or
        (s+1)*(k+s+2*r*s) > (k+s)*(r+1)**2):
        return False

    # multiplicity of eigenvalues 'r,s' (f=lambda_r, g=lambda_s)
    #
    # They are integers (checked by the 'integrality condition').
    f = -k*(s+1)*(k-s)/(mu*(r-s))
    g =  k*(r+1)*(k-r)/(mu*(r-s))
    if 1+f+g != v: # the only other eigenvalue, k, has multiplicity 1
        return False

    # 3.C of [BvL84]
    # (Absolute bound)
    if (2*v > f*(f+3) or
        2*v > g*(g+3)):
        return False

    # 3.D of [BvL84]
    # (Claw bound)
    if (mu != s**2    and
        mu != s*(s+1) and
        2*(r+1) > s*(s+1)*(mu+1)):
        return False

    # 3.E of [BvL84]
    # (the Case μ=1)
    if mu == 1:
        if (   k  % (l+1) or
            (v*k) % ((l+1)*(l+2))):
            return False

    # 3.F of [BvL84]
    # (the Case μ=2)
    if mu == 2 and 2*k < l*(l+3) and k%(l+1):
        return False

    return True

def strongly_regular_graph(int v,int k,int l,int mu=-1,bint existence=False,bint check=True):
    r"""
    Return a `(v,k,\lambda,\mu)`-strongly regular graph.

    This function relies partly on Andries Brouwer's `database of strongly
    regular graphs <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__. See
    the documentation of :mod:`sage.graphs.strongly_regular_db` for more
    information.

    INPUT:

    - ``v,k,l,mu`` (integers) -- note that ``mu``, if unspecified, is
      automatically determined from ``v,k,l``.

    - ``existence`` (boolean;``False``) -- instead of building the graph,
      return:

        - ``True`` -- meaning that a `(v,k,\lambda,\mu)`-strongly regular graph
          exists.

        - ``Unknown`` -- meaning that Sage does not know if such a strongly
          regular graph exists (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that no such strongly regular graph exists.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    EXAMPLES:

    Petersen's graph from its set of parameters::

        sage: graphs.strongly_regular_graph(10,3,0,1,existence=True)
        True
        sage: graphs.strongly_regular_graph(10,3,0,1)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices

    Now without specifying `\mu`::

        sage: graphs.strongly_regular_graph(10,3,0)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices

    An obviously infeasible set of parameters::

        sage: graphs.strongly_regular_graph(5,5,5,5,existence=True)
        False
        sage: graphs.strongly_regular_graph(5,5,5,5)
        Traceback (most recent call last):
        ...
        ValueError: There exists no (5, 5, 5, 5)-strongly regular graph

    An set of parameters proved in a paper to be infeasible::

        sage: graphs.strongly_regular_graph(324,57,0,12,existence=True)
        False
        sage: graphs.strongly_regular_graph(324,57,0,12)
        Traceback (most recent call last):
        ...
        EmptySetError: Andries Brouwer's database reports that no (324, 57, 0,
        12)-strongly regular graph exists. Comments: <a
        href="srgtabrefs.html#GavrilyukMakhnev05">Gavrilyuk & Makhnev</a> and <a
        href="srgtabrefs.html#KaskiOstergard07">Kaski & stergrd</a>

    A set of parameters unknown to be realizable in Andries Brouwer's database::

        sage: graphs.strongly_regular_graph(324,95,22,30,existence=True)
        Unknown
        sage: graphs.strongly_regular_graph(324,95,22,30)
        Traceback (most recent call last):
        ...
        RuntimeError: Andries Brouwer's database reports that no
        (324,95,22,30)-strongly regular graph is known to exist.
        Comments:

    A large unknown set of parameters (not in Andries Brouwer's database)::

        sage: graphs.strongly_regular_graph(1394,175,0,25,existence=True)
        Unknown
        sage: graphs.strongly_regular_graph(1394,175,0,25)
        Traceback (most recent call last):
        ...
        RuntimeError: Sage cannot figure out if a (1394,175,0,25)-strongly regular graph exists.

    Test the Claw bound (see 3.D of [BvL84]_)::

        sage: graphs.strongly_regular_graph(2058,242,91,20,existence=True)
        False

    TESTS:

    Check that all of our constructions are correct::

        sage: from sage.graphs.strongly_regular_db import apparently_feasible_parameters
        sage: for p in sorted(apparently_feasible_parameters(1300)):   # not tested
        ....:     if graphs.strongly_regular_graph(*p,existence=True): # not tested
        ....:         try:                                             # not tested
        ....:             _ = graphs.strongly_regular_graph(*p)        # not tested
        ....:             print p,"built successfully"                 # not tested
        ....:         except RuntimeError as e:                        # not tested
        ....:             if 'Brouwer' not in str(e):                  # not tested
        ....:                 raise                                    # not tested

    `\mu=0` behaves correctly (:trac:`19712`)::

        sage: graphs.strongly_regular_graph(10,2,1)
        Traceback (most recent call last):
        ...
        ValueError: There exists no (10, 2, 1, 0)-strongly regular graph
        sage: graphs.strongly_regular_graph(12,3,2)
        complement(Multipartite Graph with set sizes [4, 4, 4]): Graph on 12 vertices
        sage: graphs.strongly_regular_graph(6,3,0)
        Multipartite Graph with set sizes [3, 3]: Graph on 6 vertices
    """
    load_brouwer_database()
    if mu == -1:
        mu = k*(k-l-1)//(v-k-1)

    params = (v,k,l,mu)
    params_complement = (v,v-k-1,v-2*k+mu-2,v-2*k+l)

    if not seems_feasible(v,k,l,mu):
        if existence:
            return False
        raise ValueError("There exists no "+str(params)+"-strongly regular graph")

    def check_srg(G):
        if check and (v,k,l,mu) != G.is_strongly_regular(parameters=True):
            raise RuntimeError("Sage built an incorrect {}-SRG.".format((v,k,l,mu)))
        return G

    if _small_srg_database is None:
        _build_small_srg_database()

    if params in _small_srg_database:
        val = _small_srg_database[params]
        return True if existence else check_srg(val[0](*val[1:]))
    if params_complement in _small_srg_database:
        val = _small_srg_database[params_complement]
        return True if existence else check_srg(val[0](*val[1:]).complement())

    test_functions = [is_complete_multipartite, # must be 1st, to prevent 0-divisions
                      is_paley, is_johnson,
                      is_orthogonal_array_block_graph,
                      is_steiner, is_affine_polar,
                      is_goethals_seidel,
                      is_orthogonal_polar,
                      is_NOodd, is_NOperp_F5, is_NO_F2, is_NO_F3, is_NU,
                      is_unitary_polar, is_unitary_dual_polar, is_GQqmqp,
                      is_RSHCD,
                      is_twograph_descendant_of_srg,
                      is_taylor_twograph_srg,
                      is_switch_OA_srg,
                      is_polhill,
                      is_haemers,
                      is_cossidente_penttila,
                      is_mathon_PC_srg,
                      is_switch_skewhad]

    # Going through all test functions, for the set of parameters and its
    # complement.
    for f in test_functions:
        if f(*params):
            if existence:
                return True
            ans = f(*params)
            return check_srg(ans[0](*ans[1:]))
        if f(*params_complement):
            if existence:
                return True
            ans = f(*params_complement)
            return check_srg(ans[0](*ans[1:]).complement())

    # From now on, we have no idea how to build the graph.
    #
    # We try to return the most appropriate error message.

    global _brouwer_database
    brouwer_data = _brouwer_database.get(params,None)

    if brouwer_data is not None:
        if brouwer_data['status'] == 'impossible':
            if existence:
                return False
            raise EmptySetError("Andries Brouwer's database reports that no "+
                                str((v,k,l,mu))+"-strongly regular graph exists. "+
                                "Comments: "+brouwer_data['comments'].encode('ascii','ignore'))

        if brouwer_data['status'] == 'open':
            if existence:
                return Unknown
            raise RuntimeError(("Andries Brouwer's database reports that no "+
                                "({},{},{},{})-strongly regular graph is known "+
                                "to exist.\nComments: ").format(v,k,l,mu)
                               +brouwer_data['comments'].encode('ascii','ignore'))

        if brouwer_data['status'] == 'exists':
            if existence:
                return True
            raise RuntimeError(("Andries Brouwer's database claims that such a "+
                                "({},{},{},{})-strongly regular graph exists, but "+
                                "Sage does not know how to build it. If *you* do, "+
                                "please get in touch with us on sage-devel!\n"+
                                "Comments: ").format(v,k,l,mu)
                               +brouwer_data['comments'].encode('ascii','ignore'))
    if existence:
        return Unknown
    raise RuntimeError(("Sage cannot figure out if a ({},{},{},{})-strongly "+
                        "regular graph exists.").format(v,k,l,mu))

def apparently_feasible_parameters(int n):
    r"""
    Return a list of a priori feasible parameters `(v,k,\lambda,\mu)`, with `0<\mu<k`.

    Note that some of those that it returns may also be infeasible for more
    involved reasons. The condition `0<\mu<k` makes sure we skip trivial cases of
    complete multipartite graphs and their complements.

    INPUT:

    - ``n`` (integer) -- return all a-priori feasible tuples `(v,k,\lambda,\mu)`
      for `v<n`

    EXAMPLE:

    All sets of parameters with `v<20` which pass basic arithmetic tests are
    feasible::

        sage: from sage.graphs.strongly_regular_db import apparently_feasible_parameters
        sage: small_feasible = apparently_feasible_parameters(20); small_feasible
        {(5, 2, 0, 1),
         (9, 4, 1, 2),
         (10, 3, 0, 1),
         (10, 6, 3, 4),
         (13, 6, 2, 3),
         (15, 6, 1, 3),
         (15, 8, 4, 4),
         (16, 5, 0, 2),
         (16, 6, 2, 2),
         (16, 9, 4, 6),
         (16, 10, 6, 6),
         (17, 8, 3, 4)}
        sage: all(graphs.strongly_regular_graph(*x,existence=True) for x in small_feasible)
        True

    But that becomes wrong for `v<60` (because of the non-existence of a
    `(49,16,3,6)`-strongly regular graph)::

        sage: small_feasible = apparently_feasible_parameters(60)
        sage: all(graphs.strongly_regular_graph(*x,existence=True) for x in small_feasible)
        False

    """
    cdef int v,k,l,mu
    feasible = set()
    for v in range(n):
        for k in range(1,v-1):
            for l in range(k-1):
                mu = k*(k-l-1)//(v-k-1)
                if mu>0 and mu<k and seems_feasible(v,k,l,mu):
                    feasible.add((v,k,l,mu))
    return feasible

def _build_small_srg_database():
    r"""
    Build the database of small strongly regular graphs.

    This data is stored in the module-level variable ``_small_srg_database``.
    We use formulas from Cor.3.7 of [CK86]_ to compute parameters of the
    graph of the projective 2-intersection set associated with a 2-weight code `C`,
    and the usual theory of duality in association schemes to compute the
    parameters of the graph of words of `C`. Another relevant reference is
    Sect.9.8.3 of [BH12]_.

    EXAMPLE:

        sage: from sage.graphs.strongly_regular_db import _build_small_srg_database
        sage: _build_small_srg_database()

    TESTS:

    Make sure that all two-weight codes yield the strongly regular graphs we
    expect::

        sage: graphs.strongly_regular_graph(81, 50, 31, 30)
        Graph on 81 vertices
        sage: graphs.strongly_regular_graph(243, 220, 199, 200) # long time
        Graph on 243 vertices
        sage: graphs.strongly_regular_graph(256, 153, 92, 90)
        Graph on 256 vertices
        sage: graphs.strongly_regular_graph(256, 170, 114, 110)
        Graph on 256 vertices
        sage: graphs.strongly_regular_graph(256, 187, 138, 132)
        Graph on 256 vertices
        sage: graphs.strongly_regular_graph(512, 73, 12, 10)    # not tested (too long)
        Graph on 512 vertices
        sage: graphs.strongly_regular_graph(512, 219, 106, 84)  # not tested (too long)
        Graph on 512 vertices
        sage: graphs.strongly_regular_graph(512, 315, 202, 180) # not tested (too long)
        Graph on 512 vertices
        sage: graphs.strongly_regular_graph(625, 364, 213, 210) # long time
        Graph on 625 vertices
        sage: graphs.strongly_regular_graph(625, 416, 279, 272) # long time
        Graph on 625 vertices
        sage: graphs.strongly_regular_graph(625, 468, 353, 342) # not tested (too long)
        Graph on 625 vertices
        sage: graphs.strongly_regular_graph(729, 336, 153,156)  # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 420, 243, 240) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 448, 277, 272) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 476, 313, 306) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 532, 391, 380) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 560, 433, 420) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 616, 523, 506) # not tested (too long)
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(1024, 363, 122, 132)# not tested (too long)
        Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 396, 148, 156)# not tested (too long)
        Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 429, 176, 182)# not tested (too long)
        Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 825, 668, 650)# not tested (too long)
        Graph on 1024 vertices

    REFERENCES:

    .. [CK86] \R. Calderbank, W.M. Kantor,
      The geometry of two-weight codes,
      Bull. London Math. Soc. 18(1986) 97-122
    """

    from sage.graphs.generators.smallgraphs import McLaughlinGraph
    from sage.graphs.generators.smallgraphs import CameronGraph
    from sage.graphs.generators.smallgraphs import M22Graph
    from sage.graphs.generators.smallgraphs import SimsGewirtzGraph
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    from sage.graphs.generators.smallgraphs import SchlaefliGraph
    from sage.graphs.generators.smallgraphs import HigmanSimsGraph
    from sage.graphs.generators.smallgraphs import JankoKharaghaniGraph
    from sage.graphs.generators.smallgraphs import LocalMcLaughlinGraph
    from sage.graphs.generators.smallgraphs import SuzukiGraph
    from sage.graphs.generators.smallgraphs import MathonStronglyRegularGraph

    global _small_srg_database
    _small_srg_database = {
        ( 36,  14,  4,  6): [Graph,('c~rLDEOcKTPO`U`HOIj@MWFLQFAaRIT`HIWqPsQQJ'+
          'DXGLqYM@gRLAWLdkEW@RQYQIErcgesClhKefC_ygSGkZ`OyHETdK[?lWStCapVgKK')],
        ( 50,   7,  0,  1): [HoffmanSingletonGraph],
        ( 56,  10,  0,  2): [SimsGewirtzGraph],
        ( 77,  16,   0,  4): [M22Graph],
        (100,  22,   0,  6): [HigmanSimsGraph],
        (100,  44,  18, 20): [SRG_100_44_18_20],
        (100,  45,  20, 20): [SRG_100_45_20_20],
        (105,  32,   4, 12): [SRG_105_32_4_12],
        (120,  63,  30, 36): [SRG_120_63_30_36],
        (120,  77,  52, 44): [SRG_120_77_52_44],
        (126,  25,   8,  4): [SRG_126_25_8_4],
        (126,  50,  13, 24): [SRG_126_50_13_24],
        (144,  39,   6, 12): [SRG_144_39_6_12],
        (162,  56,  10, 24): [LocalMcLaughlinGraph],
        (175,  72,  20, 36): [SRG_175_72_20_36],
        (176,  49,  12, 14): [SRG_176_49_12_14],
        (176,  90,  38, 54): [SRG_176_90_38_54],
        (176, 105,  68, 54): [SRG_176_105_68_54],
        (196,  91,  42, 42): [SRG_196_91_42_42],
        (210,  99,  48, 45): [SRG_210_99_48_45],
        (220,  84,  38, 28): [SRG_220_84_38_28],
        (231,  30,   9,  3): [CameronGraph],
        (243, 110,  37, 60): [SRG_243_110_37_60],
        (253, 140,  87, 65): [SRG_253_140_87_65],
        (275, 112,  30, 56): [McLaughlinGraph],
        (276, 140,  58, 84): [SRG_276_140_58_84],
        (280, 117, 44,  52): [SRG_280_117_44_52],
        (280, 135,  70, 60): [SRG_280_135_70_60],
        (416, 100,  36, 20): [SRG_416_100_36_20],
        (560, 208,  72, 80): [SRG_560_208_72_80],
        (630,  85,  20, 10): [SRG_630_85_20_10],
        (784, 243,  82, 72): [MathonStronglyRegularGraph, 0],
        (784, 270, 98, 90):  [MathonStronglyRegularGraph, 1],
        (784, 297, 116, 110):[MathonStronglyRegularGraph, 2],
        (936, 375, 150,150): [JankoKharaghaniGraph, 936],
        (1288,792, 476,504): [SRG_1288_792_476_504],
        (1782,416, 100, 96): [SuzukiGraph],
        (1800,1029,588,588): [JankoKharaghaniGraph, 1800],
    }

    # Turns the known two-weight codes into SRG constructors
    #
    cdef int n,q,k,w1,w2,K,N,l,m,K_O,l_O,m_O
    import sage.coding.two_weight_db
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ
    cinv = matrix(ZZ, [[1,0,0],[0,0,1],[0,1,0]])
    for code in sage.coding.two_weight_db.data:
        n,q,k,w1,w2 = code['n'], code['K'].cardinality(), code['k'], code['w1'], code['w2']
        N = q**k
        K_O = n*(q-1)
        l_O = K_O**2+3*K_O-q*(w1+w2)-K_O*q*(w1+w2)+w1*w2*q**2
        m_O = (w1*w2*q**2)/N

        em = eigenmatrix(N,K_O,l_O,m_O) # 1st eigenmatrix
        assert((not em is None) and (em.det() != 0))
        emi = N*em.inverse()            # 2nd eigenmatrix
        # 1st and 2nd eigenmatrices equal up to renumbering graphs?
        selfdual = em==cinv*emi*cinv
        _small_srg_database[N,K_O,l_O,m_O] = \
            [lambda x: strongly_regular_from_two_intersection_set(x.transpose()), code['M']]
        if not selfdual: # we can build two graphs (not complements to each other!)
            K, s, r = emi[0,1], emi[1,1], emi[2,1] # by Thm 5.7 in [CK86]_.
            l = K+r*s+r+s
            m = K+r*s
            _small_srg_database[N,K,l,m] = [strongly_regular_from_two_weight_code, code['M']]

cdef load_brouwer_database():
    r"""
    Loads Andries Brouwer's database into _brouwer_database.
    """
    global _brouwer_database
    if _brouwer_database is not None:
        return
    import json

    from sage.env import SAGE_SHARE
    with open(SAGE_SHARE+"/graphs/brouwer_srg_database.json",'r') as datafile:
        _brouwer_database = {(v,k,l,mu):{'status':status,'comments':comments}
                             for (v,k,l,mu,status,comments) in json.load(datafile)}

def _check_database():
    r"""
    Checks the coherence of Andries Brouwer's database with Sage.

    The function also outputs some statistics on the database.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import _check_database
        sage: _check_database() # long time
        Sage cannot build a (196  60   14   20  ) that exists. Comment from Brouwer's database: pg(6,9,2)?
        ...
        In Andries Brouwer's database:
        - 452 impossible entries
        - 2936 undecided entries
        - 1150 realizable entries (Sage misses ... of them)

    """
    global _brouwer_database
    load_brouwer_database()

    # Check that all parameters detected as infeasible are actually infeasible
    # in Brouwer's database, for a test that was implemented.
    for params in set(_brouwer_database).difference(apparently_feasible_parameters(1301)):
        if _brouwer_database[params]['status'] != "impossible":
            raise RuntimeError("Brouwer's db does not seem to know that {} in unfeasible".format(params))
        comment = _brouwer_database[params]['comments']
        if ('Krein'    in comment or
            'Absolute' in comment or
            'Conf'     in comment or
            'mu=1'     in comment or
            '&mu;=2'   in comment):
            continue
        raise RuntimeError("We detected that {} was unfeasible, but maybe we should not have".format(params))

    # We empty the global database, to be sure that strongly_regular_graph does
    # not use its data to answer.
    _brouwer_database, saved_database = {}, _brouwer_database

    cdef int missed = 0
    for params,dic in sorted(saved_database.items()):
        sage_answer = strongly_regular_graph(*params,existence=True)
        if dic['status'] == 'open':
            if sage_answer:
                print "Sage can build a {}, Brouwer's database cannot".format(params)
            assert sage_answer is not False
        elif dic['status'] == 'exists':
            if sage_answer is not True:
                print (("Sage cannot build a ({:<4} {:<4} {:<4} {:<4}) that exists. "+
                       "Comment from Brouwer's database: ").format(*params)
                       +dic['comments'].encode('ascii','ignore'))
                missed += 1
            assert sage_answer is not False
        elif dic['status'] == 'impossible':
            assert sage_answer is not True
        else:
            assert False # must not happen

    status = [x['status'] for x in saved_database.values()]
    print "\nIn Andries Brouwer's database:"
    print "- {} impossible entries".format(status.count('impossible'))
    print "- {} undecided entries".format(status.count('open'))
    print "- {} realizable entries (Sage misses {} of them)".format(status.count('exists'),missed)

    # Reassign its value to the global database
    _brouwer_database = saved_database
