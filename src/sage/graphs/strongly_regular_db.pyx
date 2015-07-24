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

Functions
---------
"""
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown
from sage.rings.arith import is_square
from sage.rings.arith import is_prime_power
from sage.misc.cachefunc import cached_function
from sage.combinat.designs.orthogonal_arrays import orthogonal_array
from sage.combinat.designs.bibd import balanced_incomplete_block_design
from sage.graphs.generators.smallgraphs import McLaughlinGraph
from sage.graphs.generators.smallgraphs import CameronGraph
from sage.graphs.generators.smallgraphs import M22Graph
from sage.graphs.generators.smallgraphs import SimsGewirtzGraph
from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
from sage.graphs.generators.smallgraphs import SchlaefliGraph
from sage.graphs.graph import Graph
from libc.math cimport sqrt

cdef dict _brouwer_database = None

@cached_function
def is_paley(int v,int k,int l,int mu):
    r"""
    Test if a Paley graph is `(v,k,\lambda,\mu)`-strongly regular.

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
        return (lambda q : PaleyGraph(q),v)

@cached_function
def is_orthogonal_array_block_graph(int v,int k,int l,int mu):
    r"""
    Test if an Orthogonal Array graph is `(v,k,\lambda,\mu)`-strongly regular.

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

        sage: t = is_orthogonal_array_block_graph(5,5,5,5); t
    """
    # notations from
    # http://www.win.tue.nl/~aeb/graphs/OA.html
    if not is_square(v):
        return
    n = int(sqrt(v))
    if k % (n-1):
        return
    m = k//(n-1)
    if (l  != (m-1)*(m-2)+n-2 or
        mu != m*(m-1)):
        return
    if orthogonal_array(m,n,existence=True):
        from sage.graphs.generators.intersection import OrthogonalArrayBlockGraph
        return (lambda m,n : OrthogonalArrayBlockGraph(m, n), m,n)

@cached_function
def is_johnson(int v,int k,int l,int mu):
    r"""
    Test if a Johnson graph is `(v,k,\lambda,\mu)`-strongly regular.

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
    Test if a Steiner graph is `(v,k,\lambda,\mu)`-strongly regular.

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
    Test if an Affine Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

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
    from sage.rings.arith import divisors
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
            from sage.graphs.generators.families import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='+'),2*e,q)
        if (k == (q**(e-1) - 1)*(q**e+1) and
            l == q*(q**(e-2)- 1)*(q**(e-1)+1)+q-2 and
            mu== q**(e-1)*(q**(e-1) - 1)):
            from sage.graphs.generators.families import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='-'),2*e,q)

cdef eigenvalues(int v,int k,int l,int mu):
    r"""
    Return the eigenvalues of a (v,k,l,mu)-strongly regular graph.

    If the set of parameters is not feasible, ``(None,None)`` is returned
    instead.

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

def SRG_280_135_70_60():
    r"""
    Return a strongly regular graph with parameters (280, 135, 70, 60).

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

cdef bint seems_feasible(int v, int k, int l, int mu):
    r"""
    Tests is the set of parameters seems feasible

    INPUT:

    - ``v,k,l,mu`` (integers)
    """
    cdef int lambda_r, lambda_s,l1,l2,K2,D,F,e
    if (v<0 or k<=0 or l<0 or mu<0 or
        k>=v-1 or l>=k or mu>=k or
        v-2*k+mu-2 < 0 or # lambda of complement graph >=0
        v-2*k+l    < 0 or # mu of complement graph >=0
        mu*(v-k-1) != k*(k-l-1)):
        return False

    if (v-1)*(mu-l)-2*k == 0: # conference
        return True

    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return False

    # 1.3.1 of [Distance-regular graphs]
    if ((s+1)*(k-s)*k) % (mu*(s-r)):
        return False

    return True

def strongly_regular_graph(int v,int k,int l,int mu,bint existence=False):
    r"""
    Return a `(v,k,\lambda,\mu)`-strongly regular graph.

    This function relies partly on Andries Brouwer's `database of strongly
    regular graphs <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__. See
    the documentation of :mod:`sage.graphs.strongly_regular_db` for more
    information.

    INPUT:

    - ``v,k,l,mu`` (integers)

    - ``existence`` (boolean;``False``) -- instead of building the graph,
      return:

        - ``True`` -- meaning that a `(v,k,\lambda,\mu)`-strongly regular graph
          exists.

        - ``Unknown`` -- meaning that Sage does not know if such a strongly
          regular graph exists (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that no such strongly regular graph exists.


    EXAMPLES:

    Petersen's graph from its set of parameters::

        sage: graphs.strongly_regular_graph(10,3,0,1,existence=True)
        True
        sage: graphs.strongly_regular_graph(10,3,0,1)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices

    An obviously infeasible set of parameters::

        sage: graphs.strongly_regular_graph(5,5,5,5,existence=True)
        False
        sage: graphs.strongly_regular_graph(5,5,5,5)
        Traceback (most recent call last):
        ...
        ValueError: There exists no (5, 5, 5, 5)-strongly regular graph (basic arithmetic checks)

    An set of parameters proved in a paper to be infeasible::

        sage: graphs.strongly_regular_graph(324,57,0,12,existence=True)
        False
        sage: graphs.strongly_regular_graph(324,57,0,12)
        Traceback (most recent call last):
        ...
        EmptySetError: Andries Brouwer's database reports that no (324, 57, 0,
        12)-strongly regular graph exists.Comments: <a
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

    A realizable set of parameters that Sage cannot realize (help us!)::

        sage: graphs.strongly_regular_graph(279,128,52,64,existence=True)
        True
        sage: graphs.strongly_regular_graph(279,128,52,64)
        Traceback (most recent call last):
        ...
        RuntimeError: Andries Brouwer's database claims that such a
        (279,128,52,64)-strongly regular graph exists, but Sage does not know
        how to build it. If *you* do, please get in touch with us on sage-devel!
        Comments: pg(8,15,4)?; 2-graph*

    A large unknown set of parameters (not in Andries Brouwer's database)::

        sage: graphs.strongly_regular_graph(1394,175,0,25,existence=True)
        Unknown
        sage: graphs.strongly_regular_graph(1394,175,0,25)
        Traceback (most recent call last):
        ...
        RuntimeError: Sage cannot figure out if a (1394,175,0,25)-strongly regular graph exists.
    """
    load_brouwer_database()

    params = (v,k,l,mu)
    params_complement = (v,v-k-1,v-2*k+mu-2,v-2*k+l)

    if not seems_feasible(v,k,l,mu):
        if existence:
            return False
        raise ValueError("There exists no "+str(params)+"-strongly regular graph "+
                         "(basic arithmetic checks)")

    constructions = {
        ( 27,  16, 10,  8): [SchlaefliGraph],
        ( 36,  14,  4,  6): [Graph,('c~rLDEOcKTPO`U`HOIj@MWFLQFAaRIT`HIWqPsQQJ'+
          'DXGLqYM@gRLAWLdkEW@RQYQIErcgesClhKefC_ygSGkZ`OyHETdK[?lWStCapVgKK')],
        ( 40,  12,  2,  4): [Graph,('g}iS[A@_S@OA_BWQIGaPCQE@CcQGcQECXAgaOdS@a'+
          'CWEEAOIBH_HW?scb?f@GMBGGhIPGaQoh?q_bD_pGPq_WI`T_DBU?R_dECsSARGgogBO'+
          '{_IPBKZ?DI@Wgt_E?MPo{_?')],
        ( 45,  12,  3,  3): [Graph,('l~}CKMF_C?oB_FPCGaICQOaH@DQAHQ@Ch?aJHAQ@G'+
          'P_CQAIGcAJGO`IcGOY`@IGaGHGaKSCDI?gGDgGcE_@OQAg@PCSO_hOa`GIDADAD@XCI'+
          'ASDKB?oKOo@_SHCc?SGcGd@A`B?bOOHGQH?ROQOW`?XOPa@C_hcGo`CGJK')],
        ( 50,   7,  0,  1): [HoffmanSingletonGraph],
        ( 56,  10,  0,  2): [SimsGewirtzGraph],
        ( 64,  18,  2,  6): [Graph,('~?@?~aK[A@_[?O@_B_?O?K?B_?A??K??YQQPHGcQQ'+
          'CaPIOHAX?POhAPIC`GcgSAHDE?PCiC@BCcDADIG_QCocS@AST?OOceGG@QGcKcdCbCB'+
          'gIEHAScIDDOy?DAWaEg@IQO?maHPOhAW_dBCX?s@HOpKD@@GpOpHO?bCbHGOaGgpWQQ'+
          '?PDDDw@A_CSRIS_P?GeGpg`@?EOcaJGccbDC_dLAc_pHOe@`ocEGgo@sRo?WRAbAcPc'+
          '?iCiHEKBO_hOiOWpOSGSTBQCUAW_DDIWOqHBO?gghw_?`kOAXH?\\Ds@@@CpIDKOpc@'+
          'OCoeIS_YOgGATGaqAhKGA?cqDOwQKGc?')],
        ( 77,  16,  0,  4): [M22Graph],
        (231,  30,  9,  3): [CameronGraph],
        (275, 112, 30, 56): [McLaughlinGraph],
        (280, 135, 70, 60): [SRG_280_135_70_60],
    }

    if params in constructions:
        val = constructions[params]
        return True if existence else val[0](*val[1:])
    if params_complement in constructions:
        val = constructions[params_complement]
        return True if existence else val[0](*val[1:]).complement()

    test_functions = [is_paley, is_johnson,
                      is_orthogonal_array_block_graph,
                      is_steiner, is_affine_polar]

    # Going through all test functions, for the set of parameters and its
    # complement.
    for f in test_functions:
        if f(*params):
            if existence:
                return True
            ans = f(*params)
            return ans[0](*ans[1:])
        if f(*params_complement):
            if existence:
                return True
            ans = f(*params_complement)
            return ans[0](*ans[1:]).complement()

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
                                str((v,k,l,mu))+"-strongly regular graph exists."+
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
    Return a list of parameters `(v,k,\lambda,\mu)` which are a-priori feasible.

    Only basic arithmetic checks are performed on the parameters that this
    function return. Those that it does not return are infeasible for elementary
    reasons. Note that some of those that it returns may also be infeasible for
    more involved reasons.

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

    But that becomes wrong for `v<30`::

        sage: small_feasible = apparently_feasible_parameters(30)
        sage: all(graphs.strongly_regular_graph(*x,existence=True) for x in small_feasible)
        False

    """
    cdef int v,k,l,mu
    feasible = set()
    for v in range(n):
        for k in range(1,v-1):
            for l in range(k-1):
                mu = k*(k-l-1)//(v-k-1)
                if seems_feasible(v,k,l,mu):
                    feasible.add((v,k,l,mu))
    return feasible

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
        Sage cannot build a (45   22   10   11  ) that exists. Comment from Brouwer's database: <a href="srgtabrefs.html#Mathon78">Mathon</a>; 2-graph*
        ...
        In Andries Brouwer's database:
        - 448 impossible entries
        - 2950 undecided entries
        - 1140 realizable entries (Sage misses 342 of them)

    """
    global _brouwer_database
    load_brouwer_database()
    assert apparently_feasible_parameters(1301) == set(_brouwer_database)

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
