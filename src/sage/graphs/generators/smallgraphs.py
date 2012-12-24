r"""
Small graphs
==================

Add naury graphs to the list
Harary graph is not a small graph !
circle/line embeddings should not be here !
Virer GraphGenerators() from the code
"""


###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph
from sage.graphs import graph
from math import sin, cos, pi

####################
# Helper functions #
####################

def _circle_embedding(g, vertices, center=(0, 0), radius=1, shift=0):
    r"""
    Set some vertices on a circle in the embedding of a graph G.

    This method modifies the graph's embedding so that the vertices
    listed in ``vertices`` appear in this ordering on a circle of given
    radius and center. The ``shift`` parameter is actually a rotation of
    the circle. A value of ``shift=1`` will replace in the drawing the
    `i`-th element of the list by the `(i-1)`-th. Non-integer values are
    admissible, and a value of `\alpha` corresponds to a rotation of the
    circle by an angle of `\alpha 2\pi/n` (where `n` is the number of
    vertices set on the circle).

    EXAMPLE::

        sage: from sage.graphs.graph_generators import _circle_embedding
        sage: g = graphs.CycleGraph(5)
        sage: _circle_embedding(g, [0, 2, 4, 1, 3], radius=2, shift=.5)
        sage: g.show()
    """
    c_x, c_y = center
    n = len(vertices)
    d = g.get_pos()
    if d is None:
        d = {}

    for i,v in enumerate(vertices):
        i += shift
        v_x = c_x + radius * cos(2*i*pi / n)
        v_y = c_y + radius * sin(2*i*pi / n)
        d[v] = (v_x, v_y)

    g.set_pos(d)

def _line_embedding(g, vertices, first=(0, 0), last=(0, 1)):
    r"""
    Sets some vertices on a line in the embedding of a graph G.

    This method modifies the graph's embedding so that the vertices of
    ``vertices`` appear on a line, where the position of ``vertices[0]``
    is the pair ``first`` and the position of ``vertices[-1]`` is
    ``last``. The vertices are evenly spaced.

    EXAMPLE::

        sage: from sage.graphs.graph_generators import _line_embedding
        sage: g = graphs.PathGraph(5)
        sage: _line_embedding(g, [0, 2, 4, 1, 3], first=(-1, -1), last=(1, 1))
        sage: g.show()
    """
    n = len(vertices) - 1.

    fx, fy = first
    dx = (last[0] - first[0])/n
    dy = (last[1] - first[1])/n

    d = g.get_pos()
    if d is None:
        d = {}

    for v in vertices:
        d[v] = (fx, fy)
        fx += dx
        fy += dy


#######################################################################
#   Named Graphs
#######################################################################

def HararyGraph( self, k, n ):
    r"""
    Returns the Harary graph on `n` vertices and connectivity `k`, where
    `2 \leq k < n`.

    A `k`-connected graph `G` on `n` vertices requires the minimum degree
    `\delta(G)\geq k`, so the minimum number of edges `G` should have is
    `\lceil kn/2\rceil`. Harary graphs achieve this lower bound, that is,
    Harary graphs are minimal `k`-connected graphs on `n` vertices.

    The construction provided uses the method CirculantGraph.  For more
    details, see the book D. B. West, Introduction to Graph Theory, 2nd
    Edition, Prentice Hall, 2001, p. 150--151; or the `MathWorld article on
    Harary graphs <http://mathworld.wolfram.com/HararyGraph.html>`_.

    EXAMPLES:

    Harary graphs `H_{k,n}`::

        sage: h = graphs.HararyGraph(5,9); h
        Harary graph 5, 9: Graph on 9 vertices
        sage: h.order()
        9
        sage: h.size()
        23
        sage: h.vertex_connectivity()
        5

    TESTS:

    Connectivity of some Harary graphs::

        sage: n=10
        sage: for k in range(2,n):
        ...       g = graphs.HararyGraph(k,n)
        ...       if k != g.vertex_connectivity():
        ...          print "Connectivity of Harary graphs not satisfied."
    """
    if k < 2:
        raise ValueError("Connectivity parameter k should be at least 2.")
    if k >= n:
        raise ValueError("Number of vertices n should be greater than k.")

    if k%2 == 0:
        G = self.CirculantGraph( n, range(1,k/2+1) )
    else:
        if n%2 == 0:
            G = self.CirculantGraph( n, range(1,(k-1)/2+1) )
            for i in range(n):
                G.add_edge( i, (i+n/2)%n )
        else:
            G = self.HararyGraph( k-1, n )
            for i in range((n-1)/2+1):
                G.add_edge( i, (i+(n-1)/2)%n )
    G.name('Harary graph {0}, {1}'.format(k,n))
    return G

def HarriesGraph(self, embedding=1):
    r"""
    Returns the Harries Graph.

    The Harries graph is a Hamiltonian 3-regular graph on 70
    vertices. See the :wikipedia:`Wikipedia page on the Harries
    graph <Harries_graph>`.

    The default embedding here is to emphasize the graph's 4 orbits.
    This graph actually has a funny construction. The following
    procedure gives an idea of it, though not all the adjacencies
    are being properly defined.

    #. Take two disjoint copies of a :meth:`Petersen graph
       <PetersenGraph>`. Their vertices will form an orbit of the
       final graph.

    #. Subdivide all the edges once, to create 15+15=30 new
       vertices, which together form another orbit.

    #. Create 15 vertices, each of them linked to 2 corresponding
       vertices of the previous orbit, one in each of the two
       subdivided Petersen graphs. At the end of this step all
       vertices from the previous orbit have degree 3, and the only
       vertices of degree 2 in the graph are those that were just
       created.

    #. Create 5 vertices connected only to the ones from the
       previous orbit so that the graph becomes 3-regular.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.HarriesGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.HarriesGraph(embedding=2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.HarriesGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.

    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 70, [-29, -19, -13, 13, 21, -27, 27, 33, -13, 13,
                             19, -21, -33, 29], 5)
    g.name("Harries Graph")

    if embedding == 1:
        gpos = g.get_pos()
        ppos = PetersenGraph(self).get_pos()

        # The graph's four orbits
        o = [None]*4
        o[0] = [0, 2, 6, 8, 14, 16, 20, 22, 28, 30, 34, 36, 42, 44, 48, 50,
                56, 58, 62, 64]
        o[1] = [1, 3, 5, 7, 9, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35,
                37, 41, 43, 45, 47, 49, 51, 55, 57, 59, 61, 63, 65, 69]
        o[2] = [60, 10, 12, 4, 24, 26, 18, 38, 40, 32, 52, 54, 46, 66, 68]
        o[3] = [11, 25, 39, 53, 67]

        # Correspondence between the vertices of one of the two Petersen
        # graphs on o[0] and the vertices of a standard Petersen graph
        # object
        g_to_p = {0: 0, 2: 1, 42: 5, 44: 8, 14: 7, 16: 2, 56: 9, 58: 6,
                  28: 4, 30: 3}

        # Correspondence between the vertices of the other Petersen graph
        # on o[0] and the vertices of the first one
        g_to_g = {64: 44, 34: 0, 36: 28, 6: 2, 8: 58, 48: 16, 50: 30,
                  20: 14, 22: 56, 62: 42}

        # Position for the vertices from the first copy
        for v, i in g_to_p.iteritems():
            gpos[v] = ppos[i]

        # Position for the vertices in the second copy. Moves the first,
        # too.
        offset = 3.5
        for v, i in g_to_g.iteritems():
            x, y = gpos[i]
            gpos[v] = (x + offset*0.5, y)
            gpos[i] = (x - offset*0.5, y)

        # Vertices from o[1]. These are actually the "edges" of the
        # copies of Petersen.
        for v in o[1]:
            p1, p2 = [gpos[x] for x in g.neighbors(v) if x in o[0]]
            gpos[v] = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

        # 15 vertices from o[2]
        for i, v in enumerate(o[2]):
            gpos[v] = (-1.75 + i*.25, 2)

        # 5 vertices from o[3]
        for i, v in enumerate(o[3]):
            gpos[v] = (-1 + i*.5, 2.5)

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def HarriesWongGraph(self, embedding=1):
    r"""
    Returns the Harries-Wong Graph.

    See the :wikipedia:`Wikipedia page on the Harries-Wong graph
    <Harries-Wong_graph>`.

    *About the default embedding:*

    The default embedding is an attempt to emphasize the graph's
    8 (!!!) different orbits. In order to understand this better,
    one can picture the graph as being built in the following way:

        #. One first creates a 3-dimensional cube (8 vertices, 12
           edges), whose vertices define the first orbit of the
           final graph.

        #. The edges of this graph are subdivided once, to create 12
           new vertices which define a second orbit.

        #. The edges of the graph are subdivided once more, to
           create 24 new vertices giving a third orbit.

        #. 4 vertices are created and made adjacent to the vertices
           of the second orbit so that they have degree
           3. These 4 vertices also define a new orbit.

        #. In order to make the vertices from the third orbit
           3-regular (they all miss one edge), one creates a binary
           tree on 1 + 3 + 6 + 12 vertices. The leaves of this new
           tree are made adjacent to the 12 vertices of the third
           orbit, and the graph is now 3-regular. This binary tree
           contributes 4 new orbits to the Harries-Wong graph.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.HarriesWongGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: orbits = g.automorphism_group(orbits=True)[-1]
        sage: g.show(figsize=[15, 15], partition=orbits)   # long time

    Alternative embedding::

        sage: graphs.HarriesWongGraph(embedding=2).show()

    TESTS::

        sage: graphs.HarriesWongGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [9, 25, 31, -17, 17, 33, 9, -29, -15, -9, 9, 25, -25, 29, 17, -9,
         9, -27, 35, -9, 9, -17, 21, 27, -29, -9, -25, 13, 19, -9, -33,
         -17, 19, -31, 27, 11, -25, 29, -33, 13, -13, 21, -29, -21, 25,
         9, -11, -19, 29, 9, -27, -19, -13, -35, -9, 9, 17, 25, -9, 9, 27,
         -27, -21, 15, -9, 29, -29, 33, -9, -25]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 70, L, 1)
    g.name("Harries-Wong graph")

    if embedding == 1:
        d = g.get_pos()

        # Binary tree (left side)
        d[66] = (-9.5, 0)
        _line_embedding(g, [37, 65, 67], first=(-8, 2.25),
                last=(-8, -2.25))
        _line_embedding(g, [36, 38, 64, 24, 68, 30], first=(-7, 3),
                last=(-7, -3))
        _line_embedding(g, [35, 39, 63, 25, 59, 29, 11, 5, 55, 23, 69, 31],
                first=(-6, 3.5), last=(-6, -3.5))

        # Cube, corners: [9, 15, 21, 27, 45, 51, 57, 61]
        _circle_embedding(g, [61, 9], center=(0, -1.5), shift=.2,
                radius=4)
        _circle_embedding(g, [27, 15], center=(0, -1.5), shift=.7,
                radius=4*.707)
        _circle_embedding(g, [51, 21], center=(0, 2.5), shift=.2,
                radius=4)
        _circle_embedding(g, [45, 57], center=(0, 2.5), shift=.7,
                radius=4*.707)

        # Cube, subdivision
        _line_embedding(g, [21, 22, 43, 44, 45], first=d[21], last=d[45])
        _line_embedding(g, [21, 4, 3, 56, 57], first=d[21], last=d[57])
        _line_embedding(g, [57, 12, 13, 14, 15], first=d[57], last=d[15])
        _line_embedding(g, [15, 6, 7, 8, 9], first=d[15], last=d[9])
        _line_embedding(g, [9, 10, 19, 20, 21], first=d[9], last=d[21])
        _line_embedding(g, [45, 54, 53, 52, 51], first=d[45], last=d[51])
        _line_embedding(g, [51, 50, 49, 58, 57], first=d[51], last=d[57])
        _line_embedding(g, [51, 32, 33, 34, 61], first=d[51], last=d[61])
        _line_embedding(g, [61, 62, 41, 40, 27], first=d[61], last=d[27])
        _line_embedding(g, [9, 0, 1, 26, 27], first=d[9], last=d[27])
        _line_embedding(g, [27, 28, 47, 46, 45], first=d[27], last=d[45])
        _line_embedding(g, [15, 16, 17, 60, 61], first=d[15], last=d[61])

        # Top vertices
        _line_embedding(g, [2, 18, 42, 48], first=(-1, 7), last=(3, 7))

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def HallJankoGraph(self, from_string=True):
    r"""
    Returns the Hall-Janko graph.

    For more information on the Hall-Janko graph, see its
    :wikipedia:`Wikipedia page <Hall-Janko_graph>`.

    The construction used to generate this graph in Sage is by
    a 100-point permutation representation of the Janko group `J_2`,
    as described in version 3 of the ATLAS of Finite Group
    representations, in particular on the page `ATLAS: J2
    — Permutation representation on 100 points
    <http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0>`_.

    INPUT:

    - ``from_string`` (boolean) -- whether to build the graph from
      its sparse6 string or through GAP. The two methods return the
      same graph though doing it through GAP takes more time. It is
      set to ``True`` by default.

    EXAMPLES::

        sage: g = graphs.HallJankoGraph()
        sage: g.is_regular(36)
        True
        sage: g.is_vertex_transitive()
        True

    Is it really strongly regular with parameters 14, 12? ::

        sage: nu = set(g.neighbors(0))
        sage: for v in range(1, 100):
        ...     if v in nu:
        ...        expected = 14
        ...     else:
        ...        expected = 12
        ...     nv = set(g.neighbors(v))
        ...     nv.discard(0)
        ...     if len(nu & nv) != expected:
        ...        print "Something is wrong here!!!"
        ...        break

    Some other properties that we know how to check::

        sage: g.diameter()
        2
        sage: g.girth()
        3
        sage: factor(g.characteristic_polynomial())
        (x - 36) * (x - 6)^36 * (x + 4)^63

    TESTS::

        sage: gg = graphs.HallJankoGraph(from_string=False) # long time
        sage: g == gg # long time
        True
    """

    string = (":~?@c__E@?g?A?w?A@GCA_?CA`OWF`W?EAW?@?_OD@_[GAgcIaGGB@OcIA"
              "wCE@o_K_?GB@?WGAouC@OsN_?GB@O[GB`A@@_e?@OgLB_{Q_?GC@O[GAOs"
              "OCWGBA?kKBPA@?_[KB_{OCPKT`o_RD`]A?o[HBOwODW?DA?cIB?wRDP[X`"
              "ogKB_{QD@]B@o_KBPWXE`mC@o_JB?{PDPq@?oWGA_{OCPKTDp_YEwCA@_c"
              "IBOwOC`OX_OGB@?WPDPcYFg?C@_gKBp?SE@cYF`{_`?SGAOoOC`_\\FwCE"
              "A?gKBO{QD@k[FqI??_OFA_oQE@k\\Fq?`GgCB@pGRD@_XFP{a_?SE@ocIA"
              "ooNCPOUEqU@?oODA?cJB_{UEqYC@_kLC@CREPk]GAGbHgCA@?SMBpCSD`["
              "YFq?`Ga]BA?gPC`KSD`_\\Fa?cHWGB@?[IAooPD`[WF@s^HASeIg?@@OcP"
              "C`KYF@w^GQ[h`O[HAooMC@CQCpSVEPk\\GaSeIG?FA?kLB_{OC`OVE@cYG"
              "QUA@?WLBp?PC`KVEqKgJg?DA?sMBpCSDP[WEQKfIay@?_KD@_[GC`SUE@k"
              "[FaKdHa[k_?OLC@CRD@WVEpo^HAWfIAciIqoo_?CB@?kMCpOUE`o\\GAKg"
              "IQgq_?GD@_[GB?{OCpWVE@cYFACaHAWhJR?q_?CC@_kKBpC\\GACdHa[kJ"
              "a{o_?CA?oOFBpGRD@o\\GaKdIQonKrOt_?WHA`?PC`KTD`k]FqSeIaolJr"
              "CqLWCA@OkKCPGRDpcYGAKdIAgjJAsmJr?t__OE@ogJB_{XEps`HA[gIQwn"
              "KWKGAOoMBpGUE`k[Fa?aHqckJbSuLw?@?_SHA_kLC@OTFPw^GaOkLg?B@?"
              "[HA_{PDP_XFaCbHa[gIqooKRWx_?CFBpOTE@cZFPw^GACcHQgoKrSvMwWG"
              "BOwQCp_YFP{`HASfJAwnKRSx_OSSDP[WEq?aGqSfIQsoKR_zNWCE@o_HA_"
              "sREPg^GAGcHQWfIAciKbOxNg?A@__IAooMC`KTD`g\\GAKcIasoKrOtLb["
              "wMbyCA?cKBp?TD`[WE`s^GQGbHqcjJrK{NRw~_oODA?sNC@CQCpOZF@s]G"
              "QOfIaolJrGsLbk}_?OFA_sRD@SVE`k[HQcjJa{qLb[xMb|?_OOFA?cIAos"
              "RDP_ZFa?aGqOfIAsuMbk{Ns@@OsQAA_sPDPWXE`o\\FqKdIQkkJrCuLr_x"
              "Mro}NsDAPG?@@OWFApKUE@o`IQolKRKsLrc|NsQC@OWGAOgJCpOWE`o_GQ"
              "KiIqwnKr_~OcLCPS]A?oWHA_oMBpKSDP[\\FagjKBWxMbk{OSQ@@O_IAoo"
              "LBpCSD`g\\FaGbHQWgIQgmKRKwMRl?PgGC@OWHB@KSE@c[FqCaGqSeIAkk"
              "KBCqLBSuMBpGQWCA@?cKBOwRDPWVE@k^GqOfJr?pKbKtLrs}OSHDQwKIBO"
              "wPD@WWEQ?`HQWfIQglKBOtLbo}Ns@@OsTE_?kLCpWWHA[gIqomKBGwMRgz"
              "NBw~OSPDPc\\H_?CFAOoLCPSVE`o\\GAOeJAwpKbKtMrx?Qcq??OKFA?gJ"
              "B`?QDpcYEpo]FqKfIAgjJB?qKr_{NS@A__SE@o_HBO{PC`OTD`{_HaciIq"
              "{vMbt?OcPFQCeB@?SKBOwRD@SXE`k[FPw`HQ_lKRKxNRxBPC\\HQclK_?K"
              "EB?sOC`OTDa?`GqWgJRCrNBw~OSHFQStMRtDQ_?KC@OoQE`k_GaOdHa[gI"
              "q{tMBg|Nb|?OcPMSDDQSwCB@_cJB_{OCpOVFP{dHa[jJQwqKrk}NsHBQCd"
              "MRtMA?oSEA_wPDp_YEpo]GAOeIq{pLBk}NsLEQCtNTDU??OKEA_oLC@[[G"
              "aKnKBOtLbk~OCPFQStNSDLSTgGKC@GSD`[WEpw_GQGcIAciJAwpKb_xMbk"
              "~QShJRc|R`_wNCPcZF@s^GAGbHA_hJR?qKrOvMRg|NsDEPsxTTgCB@?gJB"
              "?sMC@CUDp_]FqCaHQcjJQwtLrhCPS\\IRCtQTw?B@?SHA_wPC`_aGqOiJa"
              "{oKRKvMRpFQChKRtXVUTi??ocNC@KUE@cYFaGdHa_mJrKsLb[yMro|OcXI"
              "RdPTTddZaOgJB@?UEPk[FQCfIaolJrSvMBczNR|AOsXFQCtOTtaB@?WGAP"
              "?TEPo\\GAGdHqgmKBCqLR[xMb|?PC`HQs|TTt`XUtu@?o[HB?sNCPGXF@{"
              "_GQKcIqolJb_yNCLDPs`MRtDRTTdYUwSEA?kLB`CWF@s]FqGgIqooLRgzN"
              "RxFQSlMSDDQTDXVUTi@?_KDAOoLBpKUEQOfIa{oLB_xMrt?Os\\HQcpMST"
              "HSTtl[VT}A@ocJBOwSD`_XEpo_Ha_mJrKtLbgzNSTGQspLRtDUUDp\\WG["
              "HB`CQCp[WFQGgIQgkJQ{rLbc{Nc@APsdLRt@PSt\\WUtt_Wn")

    if from_string:
        g = graph.Graph(string, loops = False, multiedges = False)
    else:

        # The following construction is due to version 3 of the ATLAS of
        # Finite Group Representations, specifically the page at
        # http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0 .

        from sage.interfaces.gap import gap
        gap.eval("g1 := (1,84)(2,20)(3,48)(4,56)(5,82)(6,67)(7,55)(8,41)"
                 "(9,35)(10,40)(11,78)(12,100)(13,49)(14,37)(15,94)(16,76)"
                 "(17,19)(18,44)(21,34)(22,85)(23,92)(24,57)(25,75)(26,28)"
                 "(27,64)(29,90)(30,97)(31,38)(32,68)(33,69)(36,53)(39,61)"
                 "(42,73)(43,91)(45,86)(46,81)(47,89)(50,93)(51,96)(52,72)"
                 "(54,74)(58,99)(59,95)(60,63)(62,83)(65,70)(66,88)(71,87)"
                 "(77,98)(79,80);")

        gap.eval("g2 := (1,80,22)(2,9,11)(3,53,87)(4,23,78)(5,51,18)"
                 "(6,37,24)(8,27,60)(10,62,47)(12,65,31)(13,64,19)"
                 "(14,61,52)(15,98,25)(16,73,32)(17,39,33)(20,97,58)"
                 "(21,96,67)(26,93,99)(28,57,35)(29,71,55)(30,69,45)"
                 "(34,86,82)(38,59,94)(40,43,91)(42,68,44)(46,85,89)"
                 "(48,76,90)(49,92,77)(50,66,88)(54,95,56)(63,74,72)"
                 "(70,81,75)(79,100,83);")

        gap.eval("G := Group([g1,g2]);")
        edges = gap('Orbit(G,[1,5],OnSets)').sage()
        g = graph.Graph([(int(u), int(v)) for u,v in edges])
        g.relabel()

    _circle_embedding(g, range(100))
    g.name("Hall-Janko graph")
    return g

def Balaban10Cage(self, embedding=1):
    r"""
    Returns the Balaban 10-cage.

    The Balaban 10-cage is a 3-regular graph with 70 vertices and
    105 edges. See its :wikipedia:`Wikipedia page
    <Balaban_10-cage>`.

    The default embedding gives a deeper understanding of the
    graph's automorphism group. It is divided into 4 layers (each
    layer being a set of points at equal distance from the drawing's
    center). From outside to inside:

    - L1: The outer layer (vertices which are the furthest from the
      origin) is actually the disjoint union of two cycles of length
      10.

    - L2: The second layer is an independent set of 20 vertices.

    - L3: The third layer is a matching on 10 vertices.

    - L4: The inner layer (vertices which are the closest from the
      origin) is also the disjoint union of two cycles of length 10.

    This graph is not vertex-transitive, and its vertices are
    partitioned into 3 orbits: L2, L3, and the union of L1 of L4
    whose elements are equivalent.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to be either 1 or 2.

    EXAMPLES::

        sage: g = graphs.Balaban10Cage()
        sage: g.girth()
        10
        sage: g.chromatic_number()
        2
        sage: g.diameter()
        6
        sage: g.is_hamiltonian()
        True
        sage: g.show(figsize=[10,10])   # long time

    TESTS::

        sage: graphs.Balaban10Cage(embedding='foo')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [-9, -25, -19, 29, 13, 35, -13, -29, 19, 25, 9, -29, 29, 17, 33,
          21, 9,-13, -31, -9, 25, 17, 9, -31, 27, -9, 17, -19, -29, 27,
          -17, -9, -29, 33, -25,25, -21, 17, -17, 29, 35, -29, 17, -17,
          21, -25, 25, -33, 29, 9, 17, -27, 29, 19, -17, 9, -27, 31, -9,
          -17, -25, 9, 31, 13, -9, -21, -33, -17, -29, 29]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 70, L, 1)
    g.name("Balaban 10-cage")

    if embedding == 2:
        return g
    elif embedding != 1:
        raise ValueError("The value of embedding must be 1 or 2.")

    L3 = [5, 24, 35, 46, 29, 40, 51, 34, 45, 56]
    _circle_embedding(g, L3, center=(0,0), radius = 4.3)

    L2  = [6, 4, 23, 25, 60, 36, 1, 47, 28, 30, 39, 41, 50, 52, 33, 9, 44,
            20, 55, 57]
    _circle_embedding(g, L2, center=(0,0), radius = 5, shift=-.5)


    L1a = [69, 68, 67, 66, 65, 64, 63, 62, 61, 0]
    L1b = [19, 18, 17, 16, 15, 14, 13, 12, 11, 10]
    _circle_embedding(g, L1a, center=(0,0), radius = 6, shift = 3.25)
    _circle_embedding(g, L1b, center=(0,0), radius = 6, shift = -1.25)

    L4a = [37, 2, 31, 38, 53, 32, 21, 54, 3, 22]
    _circle_embedding(g, L4a, center=(0,0), radius = 3, shift = 1.9)

    L4b = [26, 59, 48, 27, 42, 49, 8, 43, 58, 7]
    _circle_embedding(g, L4b, center=(0,0), radius = 3, shift = 1.1)

    return g

def Balaban11Cage(self, embedding = 1):
    r"""
    Returns the Balaban 11-cage.

    For more information, see this :wikipedia:`Wikipedia article on
    the Balaban 11-cage <Balaban_11-cage>`.

    INPUT:

    - ``embedding`` -- three embeddings are available, and can be
      selected by setting ``embedding`` to be 1, 2, or 3.

      - The first embedding is the one appearing on page 9 of the
        Fifth Annual Graph Drawing Contest report [FAGDC]_. It
        separates vertices based on their eccentricity (see
        :meth:`eccentricity()
        <sage.graphs.generic_graph.GenericGraph.eccentricity>`).

      - The second embedding has been produced just for Sage and is
        meant to emphasize the automorphism group's 6 orbits.

      - The last embedding is the default one produced by the
        :meth:`LCFGraph` constructor.

    .. NOTE::

        The vertex labeling changes according to the value of
        ``embedding=1``.

    EXAMPLES:

    Basic properties::

        sage: g = graphs.Balaban11Cage()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        11
        sage: g.diameter()
        8
        sage: g.automorphism_group().cardinality()
        64

    Our many embeddings::

        sage: g1 = graphs.Balaban11Cage(embedding=1)
        sage: g2 = graphs.Balaban11Cage(embedding=2)
        sage: g3 = graphs.Balaban11Cage(embedding=3)
        sage: g1.show(figsize=[10,10])   # long time
        sage: g2.show(figsize=[10,10])   # long time
        sage: g3.show(figsize=[10,10])   # long time

    Proof that the embeddings are the same graph::

        sage: g1.is_isomorphic(g2) # g2 and g3 are obviously isomorphic
        True

    TESTS::

        sage: graphs.Balaban11Cage(embedding='xyzzy')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1, 2, or 3.

    REFERENCES:

    .. [FAGDC] Fifth Annual Graph Drawing Contest
       P. Eaded, J. Marks, P.Mutzel, S. North
       http://www.merl.com/papers/docs/TR98-16.pdf
    """
    if embedding == 1:
        pos_dict = {}
        for j in range(8):
            for i in range(8):
                pos_dict[str(j) + str(i)]= [
                        0.8 * float(cos(2*((8*j + i)*pi/64 + pi/128))),
                        0.8 * float(sin(2*((8*j + i)*pi/64 + pi/128)))
                ]
            for i in range(4):
                pos_dict['1' + str(j) + str(i)] = [
                        1.1 * float(cos(2*((4*j + i)*pi/32 + pi/64))),
                        1.1 * float(sin(2*((4*j + i)*pi/32 + pi/64)))
                ]
            for i in range(2):
                pos_dict['1' + str(j) + str(i + 4)] = [
                        1.4 * float(cos(2*((2*j + i)*pi/16 + pi/32))),
                        1.4 * float(sin(2*((2*j + i)*pi/16 + pi/32)))
                ]

        edge_dict = {
            "00": ["11"], "01": ["10"],   "02": ["53"], "03": ["52"],
            "11": ["20"], "10": ["21"],   "53": ["22"], "52": ["23"],
            "20": ["31"], "21": ["30"],   "22": ["33"], "23": ["32"],
            "31": ["40"], "30": ["41"],   "33": ["43"], "32": ["42"],
            "40": ["50"], "41": ["51"],   "43": ["12"], "42": ["13"],
            "50": ["61"], "51": ["60"],   "12": ["63"], "13": ["62"],
            "61": ["70"], "60": ["71"],   "63": ["72"], "62": ["73"],
            "70": ["01"], "71": ["00"],   "72": ["03"], "73": ["02"],

            "04": ["35"], "05": ["34"],   "06": ["37"], "07": ["36"],
            "35": ["64"], "34": ["65"],   "37": ["66"], "36": ["67"],
            "64": ["55"], "65": ["54"],   "66": ["17"], "67": ["16"],
            "55": ["45"], "54": ["44"],   "17": ["46"], "16": ["47"],
            "45": ["74"], "44": ["75"],   "46": ["76"], "47": ["77"],
            "74": ["25"], "75": ["24"],   "76": ["27"], "77": ["26"],
            "25": ["14"], "24": ["15"],   "27": ["56"], "26": ["57"],
            "14": ["05"], "15": ["04"],   "56": ["07"], "57": ["06"],

            "100": ["03", "04"],   "110": ["10", "12"],
            "101": ["01", "06"],   "111": ["11", "13"],
            "102": ["00", "07"],   "112": ["14", "16"],
            "103": ["02", "05"],   "113": ["15", "17"],

            "120": ["22", "24"],   "130": ["33", "36"],
            "121": ["20", "26"],   "131": ["32", "37"],
            "122": ["21", "27"],   "132": ["31", "34"],
            "123": ["23", "25"],   "133": ["30", "35"],

            "140": ["43", "45"],   "150": ["50", "52"],
            "141": ["40", "46"],   "151": ["51", "53"],
            "142": ["41", "47"],   "152": ["54", "56"],
            "143": ["42", "44"],   "153": ["55", "57"],

            "160": ["60", "66"],   "170": ["73", "76"],
            "161": ["63", "65"],   "171": ["72", "77"],
            "162": ["62", "64"],   "172": ["71", "74"],
            "163": ["61", "67"],   "173": ["70", "75"],

            "104": ["100", "102", "105"],   "114": ["110", "111", "115"],
            "105": ["101", "103", "104"],   "115": ["112", "113", "114"],

            "124": ["120", "121", "125"],   "134": ["130", "131", "135"],
            "125": ["122", "123", "124"],   "135": ["132", "133", "134"],

            "144": ["140", "141", "145"],   "154": ["150", "151", "155"],
            "145": ["142", "143", "144"],   "155": ["152", "153", "154"],

            "164": ["160", "161", "165"],   "174": ["170", "171", "175"],
            "165": ["162", "163", "164"],   "175": ["172", "173", "174"]
        }

        return graph.Graph(edge_dict, pos=pos_dict, name="Balaban 11-cage")

    elif embedding == 2 or embedding == 3:
        L = [44, 26, -47, -15, 35, -39, 11, -27, 38, -37, 43, 14, 28, 51,
             -29, -16, 41, -11, -26, 15, 22, -51, -35, 36, 52, -14, -33,
             -26, -46, 52, 26, 16, 43, 33, -15, 17, -53, 23, -42, -35, -28,
             30, -22, 45, -44, 16, -38, -16, 50, -55, 20, 28, -17, -43,
             47, 34, -26, -41, 11, -36, -23, -16, 41, 17, -51, 26, -33,
             47, 17, -11, -20, -30, 21, 29, 36, -43, -52, 10, 39, -28, -17,
             -52, 51, 26, 37, -17, 10, -10, -45, -34, 17, -26, 27, -21,
             46, 53, -10, 29, -50, 35, 15, -47, -29, -41, 26, 33, 55, -17,
             42, -26, -36, 16]

        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(self, 112, L, 1)
        g.name("Balaban 11-cage")

        if embedding == 3:
            return g

        v1 = [34, 2, 54, 43, 66, 20, 89, 100, 72, 76, 6, 58, 16, 78, 74,
              70, 36, 94, 27, 25, 10, 8, 45, 60, 14, 64, 80, 82, 109, 107,
              49, 98]
        v2 = [88, 3, 19, 55, 67, 42, 101, 33, 77, 5, 17, 57, 69, 71, 73,
              75, 11, 61, 28, 9, 37, 26, 46, 95, 13, 63, 81, 83, 108, 106,
              48, 97]
        l1 = [35, 93, 1, 24, 53, 7, 44, 59, 15, 65, 79, 21, 110, 90, 50,
              99]
        l2 = [87, 4, 18, 56, 68, 41, 102, 32, 12, 62, 29, 84, 38, 105, 47,
              96]

        d = g.get_pos()
        for i,v in enumerate(v1):
            d[v] = (-2, 16.5-i)

        for i,v in enumerate(l1):
            d[v] = (-10, 8-i)

        for i,v in enumerate(l2):
            d[v] = (10, 8.5-i)

        for i,v in enumerate(v2):
            d[v] = (2, 16.5-i)

        for i,v in enumerate([0, 111, 92, 91, 52, 51, 23, 22]):
            d[v] = (-20, 14.5-4*i)

        for i,v in enumerate([104, 103, 86, 85, 40, 39, 31, 30]):
            d[v] = (20, 14.5-4*i)

        return g

    else:
        raise ValueError("The value of embedding must be 1, 2, or 3.")

def BidiakisCube(self):
    r"""
    Returns the Bidiakis cube.

    For more information, see this
    `Wikipedia article on the Bidiakis cube <http://en.wikipedia.org/wiki/Bidiakis_cube>`_.

    EXAMPLES:

    The Bidiakis cube is a 3-regular graph having 12 vertices and 18
    edges. This means that each vertex has a degree of 3. ::

        sage: g = graphs.BidiakisCube(); g
        Bidiakis cube: Graph on 12 vertices
        sage: g.show()  # long time
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.is_regular(3)
        True

    It is a Hamiltonian graph with diameter 3 and girth 4::

        sage: g.is_hamiltonian()
        True
        sage: g.diameter()
        3
        sage: g.girth()
        4

    It is a planar graph with characteristic polynomial
    `(x - 3) (x - 2) (x^4) (x + 1) (x + 2) (x^2 + x - 4)^2` and
    chromatic number 3::

        sage: g.is_planar()
        True
        sage: bool(g.characteristic_polynomial() == expand((x - 3) * (x - 2) * (x^4) * (x + 1) * (x + 2) * (x^2 + x - 4)^2))
        True
        sage: g.chromatic_number()
        3
    """
    edge_dict = {
        0:[1,6,11], 1:[2,5], 2:[3,10], 3:[4,9], 4:[5,8],
        5:[6], 6:[7], 7:[8,11], 8:[9], 9:[10], 10:[11]}
    pos_dict = {
        0: [0, 1],
        1: [0.5, 0.866025403784439],
        2: [0.866025403784439, 0.500000000000000],
        3: [1, 0],
        4: [0.866025403784439, -0.5],
        5: [0.5, -0.866025403784439],
        6: [0, -1],
        7: [-0.5, -0.866025403784439],
        8: [-0.866025403784439, -0.5],
        9: [-1, 0],
        10: [-0.866025403784439, 0.5],
        11: [-0.5, 0.866025403784439]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Bidiakis cube")

def BiggsSmithGraph(self, embedding=1):
    r"""
    Returns the Biggs-Smith graph.

    For more information, see this :wikipedia:`Wikipedia article on
    the Biggs-Smith graph <Biggs-Smith_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to be 1 or 2.

    EXAMPLES:

    Basic properties::

        sage: g = graphs.BiggsSmithGraph()
        sage: g.order()
        102
        sage: g.size()
        153
        sage: g.girth()
        9
        sage: g.diameter()
        7
        sage: g.automorphism_group().cardinality()
        2448
        sage: g.show(figsize=[10, 10])   # long time

    The other embedding::

        sage: graphs.BiggsSmithGraph(embedding=2).show()

    TESTS::

        sage: graphs.BiggsSmithGraph(embedding='xyzzy')
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.

    """
    L = [16, 24, -38, 17, 34, 48, -19, 41, -35, 47, -20, 34, -36,
         21, 14, 48, -16, -36, -43, 28, -17, 21, 29, -43, 46, -24,
         28, -38, -14, -50, -45, 21, 8, 27, -21, 20, -37, 39, -34,
         -44, -8, 38, -21, 25, 15, -34, 18, -28, -41, 36, 8, -29,
         -21, -48, -28, -20, -47, 14, -8, -15, -27, 38, 24, -48, -18,
         25, 38, 31, -25, 24, -46, -14, 28, 11, 21, 35, -39, 43, 36,
         -38, 14, 50, 43, 36, -11, -36, -24, 45, 8, 19, -25, 38, 20,
         -24, -14, -21, -8, 44, -31, -38, -28, 37]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 102, L, 1)
    g.name("Biggs-Smith graph")

    if embedding == 1:

        orbs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0],
                [17, 101, 25, 66, 20, 38, 53, 89, 48, 75, 56, 92, 45, 78,
                 34, 28, 63],
                [18, 36, 26, 65, 19, 37, 54, 90, 47, 76, 55, 91, 46, 77,
                 35, 27, 64],
                [21, 39, 52, 88, 49, 74, 57, 93, 44, 79, 33, 29, 62, 83,
                 100, 24, 67],
                [22, 97, 51, 96, 50, 95, 58, 94, 59, 80, 60, 81, 61, 82,
                 99, 23, 98],
                [30, 86, 84, 72, 70, 68, 42, 40, 31, 87, 85, 73, 71, 69,
                 43, 41, 32]]

        # central orbits
        _circle_embedding(g, orbs[1], center=(-.4, 0), radius=.2)
        _circle_embedding(g, orbs[3], center=(.4, 0), radius=.2, shift=4)

        # lower orbits
        _circle_embedding(g, orbs[0], center=(-.9, -.5), radius=.3,
                shift=2)
        _circle_embedding(g, orbs[2], center=(-.9, .5), radius=.3)

        # upper orbits
        _circle_embedding(g, orbs[4], center=(.9, -.5), radius=.3, shift=4)
        _circle_embedding(g, orbs[5], center=(.9, .5), radius=.3, shift=-2)

    elif embedding == 2:
        pass
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

    return g

def BrinkmannGraph(self):
    r"""
    Returns the Brinkmann graph.

    For more information, see the
    `Wikipedia article on the Brinkmann graph <http://en.wikipedia.org/wiki/Brinkmann_graph>`_.

    EXAMPLES:

    The Brinkmann graph is a 4-regular graph having 21 vertices and 42
    edges. This means that each vertex has degree 4. ::

        sage: G = graphs.BrinkmannGraph(); G
        Brinkmann graph: Graph on 21 vertices
        sage: G.show()  # long time
        sage: G.order()
        21
        sage: G.size()
        42
        sage: G.is_regular(4)
        True

    It is an Eulerian graph with radius 3, diameter 3, and girth 5. ::

        sage: G.is_eulerian()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        5

    The Brinkmann graph is also Hamiltonian with chromatic number 4::

        sage: G.is_hamiltonian()
        True
        sage: G.chromatic_number()
        4

    Its automorphism group is isomorphic to `D_7`::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(7))
        True
    """
    edge_dict = {
        0: [2,5,7,13],
        1: [3,6,7,8],
        2: [4,8,9],
        3: [5,9,10],
        4: [6,10,11],
        5: [11,12],
        6: [12,13],
        7: [15,20],
        8: [14,16],
        9: [15,17],
        10: [16,18],
        11: [17,19],
        12: [18,20],
        13: [14,19],
        14: [17,18],
        15: [18,19],
        16: [19,20],
        17: [20]}
    pos_dict = {
        0: [0, 4],
        1: [3.12732592987212, 2.49395920743493],
        2: [3.89971164872729, -0.890083735825258],
        3: [1.73553495647023, -3.60387547160968],
        4: [-1.73553495647023, -3.60387547160968],
        5: [-3.89971164872729, -0.890083735825258],
        6: [-3.12732592987212, 2.49395920743493],
        7: [0.867767478235116, 1.80193773580484],
        8: [1.94985582436365, 0.445041867912629],
        9: [1.56366296493606, -1.24697960371747],
        10: [0, -2],
        11: [-1.56366296493606, -1.24697960371747],
        12: [-1.94985582436365, 0.445041867912629],
        13: [-0.867767478235116, 1.80193773580484],
        14: [0.433883739117558, 0.900968867902419],
        15: [0.974927912181824, 0.222520933956314],
        16: [0.781831482468030, -0.623489801858733],
        17: [0, -1],
        18: [-0.781831482468030, -0.623489801858733],
        19: [-0.974927912181824, 0.222520933956315],
        20: [-0.433883739117558, 0.900968867902419]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Brinkmann graph")

def DoubleStarSnark(self):
    r"""
    Returns the double star snark.

    The double star snark is a 3-regular graph on 30 vertices. See
    the :wikipedia:`Wikipedia page on the double star snark
    <Double-star_snark>`.

    EXAMPLES::

        sage: g = graphs.DoubleStarSnark()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()
        False
        sage: g.automorphism_group().cardinality()
        80
        sage: g.show()
    """

    d = { 0: [1, 14, 15]
        , 1: [0, 2, 11]
        , 2: [1, 3, 7]
        , 3: [2, 4, 18]
        , 4: [3, 5, 14]
        , 5: [10, 4, 6]
        , 6: [5, 21, 7]
        , 7: [8, 2, 6]
        , 8: [9, 13, 7]
        , 9: [24, 8, 10]
        , 10: [9, 11, 5]
        , 11: [1, 10, 12]
        , 12: [11, 27, 13]
        , 13: [8, 12, 14]
        , 14: [0, 4, 13]
        , 15: [0, 16, 29]
        , 16: [15, 20, 23]
        , 17: [25, 18, 28]
        , 18: [3, 17, 19]
        , 19: [18, 26, 23]
        , 20: [16, 28, 21]
        , 21: [20, 6, 22]
        , 22: [26, 21, 29]
        , 23: [16, 24, 19]
        , 24: [25, 9, 23]
        , 25: [24, 17, 29]
        , 26: [27, 19, 22]
        , 27: [12, 26, 28]
        , 28: [17, 27, 20]
        , 29: [25, 22, 15]
        }

    g = graph.Graph(d, pos={}, name="Double star snark")
    _circle_embedding(g, range(15), radius=2)
    _circle_embedding(g, range(15, 30), radius=1.4)

    return g

def ChvatalGraph(self):
    r"""
    Returns the Chvatal graph.

    Chvatal graph is one of the few known graphs to satisfy Grunbaum's
    conjecture that for every m, n, there is an m-regular,
    m-chromatic graph of girth at least n. For more information, see this
    `Wikipedia article on the Chvatal graph <http://en.wikipedia.org/wiki/Chv%C3%A1tal_graph>`_.

    EXAMPLES:

    The Chvatal graph has 12 vertices and 24 edges. It is a 4-regular,
    4-chromatic graph with radius 2, diameter 2, and girth 4. ::

        sage: G = graphs.ChvatalGraph(); G
        Chvatal graph: Graph on 12 vertices
        sage: G.order(); G.size()
        12
        24
        sage: G.degree()
        [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        sage: G.chromatic_number()
        4
        sage: G.radius(); G.diameter(); G.girth()
        2
        2
        4
    """
    import networkx
    pos_dict = {}
    for i in range(5, 10):
        x = float(cos((pi / 2) + ((2 * pi) / 5) * i))
        y = float(sin((pi / 2) + ((2 * pi) / 5) * i))
        pos_dict[i] = (x, y)
    for i in range(5):
        x = float(2 * (cos((pi / 2) + ((2 * pi) / 5) * (i - 5))))
        y = float(2 * (sin((pi / 2) + ((2 * pi) / 5) * (i - 5))))
        pos_dict[i] = (x, y)
    pos_dict[10] = (0.5, 0)
    pos_dict[11] = (-0.5, 0)

    return graph.Graph(networkx.chvatal_graph(), pos=pos_dict, name="Chvatal graph")

def ClebschGraph(self):
    r"""
    Return the Clebsch graph.

    EXAMPLES::

        sage: g = graphs.ClebschGraph()
        sage: g.automorphism_group().cardinality()
        1920
        sage: g.girth()
        4
        sage: g.chromatic_number()
        4
        sage: g.diameter()
        2
        sage: g.show(figsize=[10, 10]) # long time
    """
    g = graph.Graph(pos={})
    x = 0
    for i in range(8):
        g.add_edge(x % 16, (x + 1) % 16)
        g.add_edge(x % 16, (x + 6) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1
        g.add_edge(x % 16, (x + 3) % 16)
        g.add_edge(x % 16, (x + 2) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1

    _circle_embedding(g, range(16), shift=.5)
    g.name("Clebsch graph")

    return g

def CoxeterGraph(self):
    r"""
    Return the Coxeter graph.

    See the :wikipedia:`Wikipedia page on the Coxeter graph
    <Coxeter_graph>`.

    EXAMPLES::

        sage: g = graphs.CoxeterGraph()
        sage: g.automorphism_group().cardinality()
        336
        sage: g.girth()
        7
        sage: g.chromatic_number()
        3
        sage: g.diameter()
        4
        sage: g.show(figsize=[10, 10]) # long time
    """
    g = graph.Graph({
            27: [6, 22, 14],
            24: [0, 7, 18],
            25: [8, 15, 2],
            26: [10, 16, 23],
            }, pos={})

    g.add_cycle(range(24))
    g.add_edges([(5, 11), (9, 20), (12, 1), (13, 19), (17, 4), (3, 21)])

    _circle_embedding(g, range(24))
    _circle_embedding(g, [24, 25, 26], radius=.5)
    g.get_pos()[27] = (0, 0)

    g.name("Coxeter Graph")

    return g

def DesarguesGraph(self):
    """
    Returns the Desargues graph.

    PLOTTING: The layout chosen is the same as on the cover of [1].

    EXAMPLE::

        sage: D = graphs.DesarguesGraph()
        sage: L = graphs.LCFGraph(20,[5,-5,9,-9],5)
        sage: D.is_isomorphic(L)
        True
        sage: D.show()  # long time

    REFERENCE:

    - [1] Harary, F. Graph Theory. Reading, MA: Addison-Wesley,
      1994.
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G = GeneralizedPetersenGraph(self,10,3)
    G.name("Desargues Graph")
    return G

def DurerGraph(self):
    r"""
    Returns the Dürer graph.

    For more information, see this
    `Wikipedia article on the Dürer graph <http://en.wikipedia.org/wiki/D%C3%BCrer_graph>`_.

    EXAMPLES:

    The Dürer graph is named after Albrecht Dürer. It is a planar graph
    with 12 vertices and 18 edges. ::

        sage: G = graphs.DurerGraph(); G
        Durer graph: Graph on 12 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Dürer graph has chromatic number 3, diameter 4, and girth 3. ::

        sage: G.chromatic_number()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3

    Its automorphism group is isomorphic to `D_6`. ::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,5,6],
        1: [2,7],
        2: [3,8],
        3: [4,9],
        4: [5,10],
        5: [11],
        6: [8,10],
        7: [9,11],
        8: [10],
        9: [11]}
    pos_dict = {
        0: [2, 0],
        1: [1, 1.73205080756888],
        2: [-1, 1.73205080756888],
        3: [-2, 0],
        4: [-1, -1.73205080756888],
        5: [1, -1.73205080756888],
        6: [1, 0],
        7: [0.5, 0.866025403784439],
        8: [-0.5, 0.866025403784439],
        9: [-1, 0],
        10: [-0.5, -0.866025403784439],
        11: [0.5, -0.866025403784439]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Durer graph")

def DyckGraph(self):
    """
    Returns the Dyck graph.

    For more information, see the `MathWorld article on the Dyck graph
    <http://mathworld.wolfram.com/DyckGraph.html>`_ or the `Wikipedia
    article on the Dyck graph <http://en.wikipedia.org/wiki/Dyck_graph>`_.

    EXAMPLES:

    The Dyck graph was defined by Walther von Dyck in 1881. It has `32`
    vertices and `48` edges, and is a cubic graph (regular of degree `3`)::

        sage: G = graphs.DyckGraph(); G
        Dyck graph: Graph on 32 vertices
        sage: G.order()
        32
        sage: G.size()
        48
        sage: G.is_regular()
        True
        sage: G.is_regular(3)
        True

    It is non-planar and Hamiltonian, as well as bipartite (making it a
    bicubic graph)::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()
        True
        sage: G.is_bipartite()
        True

    It has radius `5`, diameter `5`, and girth `6`::

        sage: G.radius()
        5
        sage: G.diameter()
        5
        sage: G.girth()
        6

    Its chromatic number is `2` and its automorphism group is of order
    `192`::

        sage: G.chromatic_number()
        2
        sage: G.automorphism_group().cardinality()
        192

    It is a non-integral graph as it has irrational eigenvalues::

        sage: G.characteristic_polynomial().factor()
        (x - 3) * (x + 3) * (x - 1)^9 * (x + 1)^9 * (x^2 - 5)^6

    It is a toroidal graph, and its embedding on a torus is dual to an
    embedding of the Shrikhande graph (:meth:`ShrikhandeGraph
    <GraphGenerators.ShrikhandeGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i]  = [0.75 * pos_dict[i][0],
                            0.75 * pos_dict[i][1]]
        pos_dict[16 + i] = [0.50 * pos_dict[i][0],
                            0.50 * pos_dict[i][1]]
        pos_dict[24 + i] = [0.25 * pos_dict[i][0],
                            0.25 * pos_dict[i][1]]

    edge_dict = {
        0O00: [0O07, 0O01,   0O10], 0O10: [0O00,   0O27, 0O21],
        0O01: [0O00, 0O02,   0O11], 0O11: [0O01,   0O20, 0O22],
        0O02: [0O01, 0O03,   0O12], 0O12: [0O02,   0O21, 0O23],
        0O03: [0O02, 0O04,   0O13], 0O13: [0O03,   0O22, 0O24],
        0O04: [0O03, 0O05,   0O14], 0O14: [0O04,   0O23, 0O25],
        0O05: [0O04, 0O06,   0O15], 0O15: [0O05,   0O24, 0O26],
        0O06: [0O05, 0O07,   0O16], 0O16: [0O06,   0O25, 0O27],
        0O07: [0O06, 0O00,   0O17], 0O17: [0O07,   0O26, 0O20],

        0O20: [0O17, 0O11,   0O30], 0O30: [0O20,   0O35, 0O33],
        0O21: [0O10, 0O12,   0O31], 0O31: [0O21,   0O36, 0O34],
        0O22: [0O11, 0O13,   0O32], 0O32: [0O22,   0O37, 0O35],
        0O23: [0O12, 0O14,   0O33], 0O33: [0O23,   0O30, 0O36],
        0O24: [0O13, 0O15,   0O34], 0O34: [0O24,   0O31, 0O37],
        0O25: [0O14, 0O16,   0O35], 0O35: [0O25,   0O32, 0O30],
        0O26: [0O15, 0O17,   0O36], 0O36: [0O26,   0O33, 0O31],
        0O27: [0O16, 0O10,   0O37], 0O37: [0O27,   0O34, 0O32],
    }

    return graph.Graph(edge_dict, pos=pos_dict, name="Dyck graph")

def EllinghamHorton54Graph(self):
    r"""
    Returns the Ellingham-Horton 54-graph.

    For more information, see the :wikipedia:`Wikipedia page on the
    Ellingham-Horton graphs <Ellingham-Horton_graph>`

    EXAMPLE:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton54Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity() # not tested - too long
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian() # not tested - too long
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10, 10]) # not tested - too long

    TESTS::

        sage: g.show() # long time
    """
    from sage.graphs.graph_generators import GraphGenerators
    up = GraphGenerators().CycleGraph(16)
    low = 2*GraphGenerators().CycleGraph(6)

    for v in range(6):
        low.add_edge(v, v + 12)
        low.add_edge(v + 6, v + 12)
    low.add_edge(12, 15)
    low.delete_edge(1, 2)
    low.delete_edge(8, 7)
    low.add_edge(1, 8)
    low.add_edge(7, 2)


    # The set of vertices on top is 0..15
    # Bottom left is 16..33
    # Bottom right is 34..52
    # The two other vertices are 53, 54
    g = up + 2*low
    g.name("Ellingham-Horton 54-graph")
    g.set_pos({})

    g.add_edges([(15, 4), (3, 8), (7, 12), (11, 0), (2, 13), (5, 10)])
    g.add_edges([(30, 6), (29, 9), (48, 14), (47, 1)])
    g.add_edge(32, 52)
    g.add_edge(50, 52)
    g.add_edge(33, 53)
    g.add_edge(51, 53)
    g.add_edge(52, 53)

    # Top
    _circle_embedding(g, range(16), center=(0, .5), shift=.5, radius=.5)

    # Bottom-left
    _circle_embedding(g, range(16, 22), center=(-1.5, -1))
    _circle_embedding(g, range(22, 28), center=(-1.5, -1), radius=.5)
    _circle_embedding(g, range(28, 34), center=(-1.5, -1), radius=.7)

    # Bottom right
    _circle_embedding(g, range(34, 40), center=(1.5, -1))
    _circle_embedding(g, range(40, 46), center=(1.5, -1), radius=.5)
    _circle_embedding(g, range(46, 52), center=(1.5, -1), radius=.7)

    d = g.get_pos()
    d[52] = (-.3, -2.5)
    d[53] = (.3, -2.5)
    d[31] = (-2.2, -.9)
    d[28] = (-.8, -.9)
    d[46] = (2.2, -.9)
    d[49] = (.8, -.9)


    return g

def EllinghamHorton78Graph(self):
    r"""
    Returns the Ellingham-Horton 78-graph.

    For more information, see the :wikipedia:`Wikipedia page on the
    Ellingham-Horton graphs
    <http://en.wikipedia.org/wiki/Ellingham%E2%80%93Horton_graph>`

    EXAMPLE:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton78Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity() # not tested - too long
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian() # not tested - too long
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10,10]) # not tested - too long

    TESTS::

        sage: g.show(figsize=[10, 10]) # not tested - too long
    """
    g = graph.Graph({
            0: [1, 5, 60], 1: [2, 12], 2: [3, 7], 3: [4, 14], 4: [5, 9],
            5: [6], 6: [7, 11], 7: [15], 8: [9, 13, 22], 9: [10],
            10: [11, 72], 11: [12], 12: [13], 13: [14], 14: [72],
            15: [16, 20], 16: [17, 27], 17: [18, 22], 18: [19, 29],
            19: [20, 24], 20: [21], 21: [22, 26], 23: [24, 28, 72],
            24: [25], 25: [26, 71], 26: [27], 27: [28], 28: [29],
            29: [69], 30: [31, 35, 52], 31: [32, 42], 32: [33, 37],
            33: [34, 43], 34: [35, 39], 35: [36], 36: [41, 63],
            37: [65, 66], 38: [39, 59, 74], 39: [40], 40: [41, 44],
            41: [42], 42: [74], 43: [44, 74], 44: [45], 45: [46, 50],
            46: [47, 57], 47: [48, 52], 48: [49, 75], 49: [50, 54],
            50: [51], 51: [52, 56], 53: [54, 58, 73], 54: [55],
            55: [56, 59], 56: [57], 57: [58], 58: [75], 59: [75],
            60: [61, 64], 61: [62, 71], 62: [63, 77], 63: [67],
            64: [65, 69], 65: [77], 66: [70, 73], 67: [68, 73],
            68: [69, 76], 70: [71, 76], 76: [77]}, pos={})

    _circle_embedding(g, range(15), center=(-2.5, 1.5))
    _circle_embedding(g, range(15, 30), center=(-2.5, -1.5))
    _circle_embedding(g, [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
        42, 74, 43, 44], center=(2.5, 1.5))
    _circle_embedding(g, [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
        57, 58, 75, 59], center=(2.5, -1.5))

    d = g.get_pos()

    d[76] = (-.2, -.1)
    d[77] = (.2, .1)
    d[38] = (2.2, .1)
    d[52] = (2.3, -.1)
    d[15] = (-2.1, -.1)
    d[72] = (-2.1, .1)

    _line_embedding(g, [60, 61, 62, 63], first=(-1, 2), last=(1, 2))
    _line_embedding(g, [64, 65, 37], first=(-.5, 1.5), last=(1.2, 1.5))
    _line_embedding(g, [66, 73, 67, 68, 69], first=(1.2, -2),
            last=(-.8, -2))
    _line_embedding(g, [66, 70, 71], first=(.7, -1.5), last=(-1, -1.5))

    g.name("Ellingham-Horton 78-graph")

    return g

def ErreraGraph(self):
    r"""
    Returns the Errera graph.

    For more information, see this
    `Wikipedia article on the Errera graph <http://en.wikipedia.org/wiki/Errera_graph>`_.

    EXAMPLES:

    The Errera graph is named after Alfred Errera. It is a planar graph
    on 17 vertices and having 45 edges. ::

        sage: G = graphs.ErreraGraph(); G
        Errera graph: Graph on 17 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        17
        sage: G.size()
        45

    The Errera graph is Hamiltonian with radius 3, diameter 4, girth 3,
    and chromatic number 4. ::

        sage: G.is_hamiltonian()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3
        sage: G.chromatic_number()
        4

    Each vertex degree is either 5 or 6. That is, if `f` counts the
    number of vertices of degree 5 and `s` counts the number of vertices
    of degree 6, then `f + s` is equal to the order of the Errera
    graph. ::

        sage: D = G.degree_sequence()
        sage: D.count(5) + D.count(6) == G.order()
        True

    The automorphism group of the Errera graph is isomorphic to the
    dihedral group of order 20. ::

        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(10))
        True
    """
    edge_dict = {
        0: [1,7,14,15,16],
        1: [2,9,14,15],
        2: [3,8,9,10,14],
        3: [4,9,10,11],
        4: [5,10,11,12],
        5: [6,11,12,13],
        6: [7,8,12,13,16],
        7: [13,15,16],
        8: [10,12,14,16],
        9: [11,13,15],
        10: [12],
        11: [13],
        13: [15],
        14: [16]}
    return graph.Graph(edge_dict, name="Errera graph")

def FlowerSnark(self):
    """
    Returns a Flower Snark.

    A flower snark has 20 vertices. It is part of the class of
    biconnected cubic graphs with edge chromatic number = 4, known as
    snarks. (i.e.: the Petersen graph). All snarks are not Hamiltonian,
    non-planar and have Petersen graph graph minors.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the nodes are
    drawn 0-14 on the outer circle, and 15-19 in an inner pentagon.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Flower Snark - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/FlowerSnark.html [2007, February 17]

    EXAMPLES: Inspect a flower snark::

        sage: F = graphs.FlowerSnark()
        sage: F
        Flower Snark: Graph on 20 vertices
        sage: F.graph6_string()
        'ShCGHC@?GGg@?@?Gp?K??C?CA?G?_G?Cc'

    Now show it::

        sage: F.show() # long time
    """
    pos_dict = {}
    for i in range(15):
        x = float(2.5*(cos((pi/2) + ((2*pi)/15)*i)))
        y = float(2.5*(sin((pi/2) + ((2*pi)/15)*i)))
        pos_dict[i] = (x,y)
    for i in range(15,20):
        x = float(cos((pi/2) + ((2*pi)/5)*i))
        y = float(sin((pi/2) + ((2*pi)/5)*i))
        pos_dict[i] = (x,y)
    return graph.Graph({0:[1,14,15],1:[2,11],2:[3,7],3:[2,4,16],4:[5,14], \
                        5:[6,10],6:[5,7,17],8:[7,9,13],9:[10,18],11:[10,12], \
                        12:[13,19],13:[14],15:[19],16:[15,17],18:[17,19]}, \
                        pos=pos_dict, name="Flower Snark")

def FosterGraph(self):
    """
    Returns the Foster graph.

    See the :wikipedia:`Wikipedia page on the Foster Graph
    <Foster_graph>`.

    EXAMPLE::

        sage: g = graphs.FosterGraph()
        sage: g.order()
        90
        sage: g.size()
        135
        sage: g.diameter()
        8
        sage: g.girth()
        10
        sage: g.automorphism_group().cardinality()
        4320
        sage: g.is_hamiltonian()
        True
    """
    from sage.graphs.generators.families import LCFGraph
    g= LCFGraph(self, 90, [17, -9, 37, -37, 9, -17], 15)
    g.name("Foster Graph")
    return g


def FranklinGraph(self):
    r"""
    Returns the Franklin graph.

    For more information, see this
    `Wikipedia article on the Franklin graph <http://en.wikipedia.org/wiki/Franklin_graph>`_.

    EXAMPLES:

    The Franklin graph is named after Philip Franklin. It is a
    3-regular graph on 12 vertices and having 18 edges. ::

        sage: G = graphs.FranklinGraph(); G
        Franklin graph: Graph on 12 vertices
        sage: G.is_regular(3)
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Franklin graph is a Hamiltonian, bipartite graph with radius 3,
    diameter 3, and girth 4. ::

        sage: G.is_hamiltonian()
        True
        sage: G.is_bipartite()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        4

    It is a perfect, triangle-free graph having chromatic number 2. ::

        sage: G.is_perfect()
        True
        sage: G.is_triangle_free()
        True
        sage: G.chromatic_number()
        2
    """
    edge_dict = {
        0: [1,5,6],
        1: [2,7],
        2: [3,8],
        3: [4,9],
        4: [5,10],
        5: [11],
        6: [7,9],
        7: [10],
        8: [9,11],
        10: [11]}
    pos_dict = {
        0: [2, 0],
        1: [1, 1.73205080756888],
        2: [-1, 1.73205080756888],
        3: [-2, 0],
        4: [-1, -1.73205080756888],
        5: [1, -1.73205080756888],
        6: [1, 0],
        7: [0.5, 0.866025403784439],
        8: [-0.5, 0.866025403784439],
        9: [-1, 0],
        10: [-0.5, -0.866025403784439],
        11: [0.5, -0.866025403784439]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Franklin graph")

def FruchtGraph(self):
    """
    Returns a Frucht Graph.

    A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic
    identity graph. It is planar and it is Hamiltonian.

    This constructor is dependent on NetworkX's numeric labeling.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the first
    seven nodes are on the outer circle, with the next four on an inner
    circle and the last in the center.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Frucht Graph - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/FruchtGraph.html [2007, February 17]

    EXAMPLES::

        sage: FRUCHT = graphs.FruchtGraph()
        sage: FRUCHT
        Frucht graph: Graph on 12 vertices
        sage: FRUCHT.graph6_string()
        'KhCKM?_EGK?L'
        sage: (graphs.FruchtGraph()).show() # long time
    """
    pos_dict = {}
    for i in range(7):
        x = float(2*(cos((pi/2) + ((2*pi)/7)*i)))
        y = float(2*(sin((pi/2) + ((2*pi)/7)*i)))
        pos_dict[i] = (x,y)
    pos_dict[7] = (0,1)
    pos_dict[8] = (-1,0)
    pos_dict[9] = (0,-1)
    pos_dict[10] = (1,0)
    pos_dict[11] = (0,0)
    import networkx
    G = networkx.frucht_graph()
    return graph.Graph(G, pos=pos_dict, name="Frucht graph")

def GoldnerHararyGraph(self):
    r"""
    Return the Goldner-Harary graph.

    For more information, see this
    `Wikipedia article on the Goldner-Harary graph <http://en.wikipedia.org/wiki/Goldner%E2%80%93Harary_graph>`_.

    EXAMPLES:

    The Goldner-Harary graph is named after A. Goldner and Frank Harary.
    It is a planar graph having 11 vertices and 27 edges. ::

        sage: G = graphs.GoldnerHararyGraph(); G
        Goldner-Harary graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        11
        sage: G.size()
        27

    The Goldner-Harary graph is chordal with radius 2, diameter 2, and
    girth 3. ::

        sage: G.is_chordal()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is 4 and its automorphism group is isomorphic to
    the dihedral group `D_6`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,3,4],
        1: [2,3,4,5,6,7,10],
        2: [3,7],
        3: [7,8,9,10],
        4: [3,5,9,10],
        5: [10],
        6: [7,10],
        7: [8,10],
        8: [10],
        9: [10]}

    pos = {
        0: (-2, 0),
        1: (0, 1.5),
        2: (2, 0),
        3: (0, -1.5),
        4: (-1.5, 0),
        5: (-0.5, 0.5),
        6: (0.5, 0.5),
        7: (1.5, 0),
        8: (0.5, -0.5),
        9: (-0.5, -0.5),
        10: (0, 0)}

    return graph.Graph(edge_dict, pos = pos, name="Goldner-Harary graph")

def GrayGraph(self, embedding=1):
    r"""
    Returns the Gray graph.

    See the :wikipedia:`Wikipedia page on the Gray Graph
    <Gray_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.GrayGraph()
        sage: g.order()
        54
        sage: g.size()
        81
        sage: g.girth()
        8
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.GrayGraph(embedding = 2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.GrayGraph(embedding = 3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1, 2, or 3.
    """

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 54, [-25,7,-7,13,-13,25], 9)
    g.name("Gray graph")

    if embedding == 1:
        o = g.automorphism_group(orbits=True)[-1]
        _circle_embedding(g, o[0], center=(0, 0), radius=1)
        _circle_embedding(g, o[1], center=(0, 0), radius=.6, shift=-.5)

    elif embedding != 2:
        raise ValueError("The value of embedding must be 1, 2, or 3.")

    return g

def GrotzschGraph(self):
    r"""
    Returns the Grötzsch graph.

    The Grötzsch graph is an example of a triangle-free graph with
    chromatic number equal to 4. For more information, see this
    `Wikipedia article on Grötzsch graph <http://en.wikipedia.org/wiki/Gr%C3%B6tzsch_graph>`_.

    REFERENCE:

    - [1] Weisstein, Eric W. "Grotzsch Graph."
      From MathWorld--A Wolfram Web Resource.
      http://mathworld.wolfram.com/GroetzschGraph.html

    EXAMPLES:

    The Grötzsch graph is named after Herbert Grötzsch. It is a
    Hamiltonian graph with 11 vertices and 20 edges. ::

        sage: G = graphs.GrotzschGraph(); G
        Grotzsch graph: Graph on 11 vertices
        sage: G.is_hamiltonian()
        True
        sage: G.order()
        11
        sage: G.size()
        20

    The Grötzsch graph is triangle-free and having radius 2, diameter 2,
    and girth 4. ::

        sage: G.is_triangle_free()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        4

    Its chromatic number is 4 and its automorphism group is isomorphic
    to the dihedral group `D_5`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(5))
        True
    """
    g = graph.Graph()
    g.add_vertices(range(11))

    edges = [];
    for u in range(1,6):
        edges.append( (0,u) )

    edges.append( (10,6) )

    for u in range(6,10):
        edges.append( (u,u+1) )
        edges.append( (u,u-4) )

    edges.append( (10,1) )

    for u in range(7,11):
        edges.append( (u,u-6) )

    edges.append((6,5))

    g.add_edges(edges)

    pos = {}
    pos[0] = (0,0)
    for u in range(1,6):
        theta = (u-1)*2*pi/5
        pos[u] = (float(5*sin(theta)),float(5*cos(theta)))
        pos[u+5] = (2*pos[u][0], 2*pos[u][1])

    g.set_pos(pos)
    g.name("Grotzsch graph")
    return g

def HeawoodGraph(self):
    """
    Returns a Heawood graph.

    The Heawood graph is a cage graph that has 14 nodes. It is a cubic
    symmetric graph. (See also the Moebius-Kantor graph). It is
    nonplanar and Hamiltonian. It has diameter = 3, radius = 3, girth =
    6, chromatic number = 2. It is 4-transitive but not 5-transitive.

    This constructor is dependent on NetworkX's numeric labeling.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the nodes are
    positioned in a circular layout with the first node appearing at
    the top, and then continuing counterclockwise.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Heawood Graph - from Wolfram
      MathWorld". [Online] Available:
      http://mathworld.wolfram.com/HeawoodGraph.html [2007, February 17]

    EXAMPLES::

        sage: H = graphs.HeawoodGraph()
        sage: H
        Heawood graph: Graph on 14 vertices
        sage: H.graph6_string()
        'MhEGHC@AI?_PC@_G_'
        sage: (graphs.HeawoodGraph()).show() # long time
    """
    pos_dict = {}
    for i in range(14):
        x = float(cos((pi/2) + (pi/7)*i))
        y = float(sin((pi/2) + (pi/7)*i))
        pos_dict[i] = (x,y)
    import networkx
    G = networkx.heawood_graph()
    return graph.Graph(G, pos=pos_dict, name="Heawood graph")

def HerschelGraph(self):
    r"""
    Returns the Herschel graph.

    For more information, see this
    `Wikipedia article on the Herschel graph <http://en.wikipedia.org/wiki/Herschel_graph>`_.

    EXAMPLES:

    The Herschel graph is named after Alexander Stewart Herschel. It is
    a planar, bipartite graph with 11 vertices and 18 edges. ::

        sage: G = graphs.HerschelGraph(); G
        Herschel graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.is_bipartite()
        True
        sage: G.order()
        11
        sage: G.size()
        18

    The Herschel graph is a perfect graph with radius 3, diameter 4, and
    girth 4. ::

        sage: G.is_perfect()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        4

    Its chromatic number is 2 and its automorphism group is
    isomorphic to the dihedral group `D_6`. ::

        sage: G.chromatic_number()
        2
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(6))
        True
    """
    edge_dict = {
        0: [1,3,4],
        1: [2,5,6],
        2: [3,7],
        3: [8,9],
        4: [5,9],
        5: [10],
        6: [7,10],
        7: [8],
        8: [10],
        9: [10]}
    pos_dict = {
        0: [2, 0],
        1: [0, 2],
        2: [-2, 0],
        3: [0, -2],
        4: [1, 0],
        5: [0.5, 0.866025403784439],
        6: [-0.5, 0.866025403784439],
        7: [-1, 0],
        8: [-0.5, -0.866025403784439],
        9: [0.5, -0.866025403784439],
        10: [0, 0]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Herschel graph")

def HigmanSimsGraph(self, relabel=True):
    r"""
    The Higman-Sims graph is a remarkable strongly regular
    graph of degree 22 on 100 vertices.  For example, it can
    be split into two sets of 50 vertices each, so that each
    half induces a subgraph isomorphic to the
    Hoffman-Singleton graph
    (:meth:`~HoffmanSingletonGraph`).
    This can be done in 352 ways (see [BROUWER-HS-2009]_).

    Its most famous property is that the automorphism
    group has an index 2 subgroup which is one of the
    26 sporadic groups. [HIGMAN1968]_

    The construction used here follows [HAFNER2004]_.

    INPUT:

    - ``relabel`` - default: ``True``.  If ``True`` the
      vertices will be labeled with consecutive integers.
      If ``False`` the labels are strings that are three
      digits long. "xyz" means the vertex is in group
      x (zero through three), pentagon or pentagram y
      (zero through four), and is vertex z (zero
      through four) of that pentagon or pentagram.
      See [HAFNER2004]_ for more.

    OUTPUT:

    The Higman-Sims graph.

    EXAMPLES:

    A split into the first 50 and last 50 vertices
    will induce two copies of the Hoffman-Singleton graph,
    and we illustrate another such split, which is obvious
    based on the construction used. ::

        sage: H = graphs.HigmanSimsGraph()
        sage: A = H.subgraph(range(0,50))
        sage: B = H.subgraph(range(50,100))
        sage: K = graphs.HoffmanSingletonGraph()
        sage: K.is_isomorphic(A) and K.is_isomorphic(B)
        True
        sage: C = H.subgraph(range(25,75))
        sage: D = H.subgraph(range(0,25)+range(75,100))
        sage: K.is_isomorphic(C) and K.is_isomorphic(D)
        True

    The automorphism group contains only one nontrivial
    proper normal subgroup, which is of index 2 and is
    simple.  It is known as the Higman-Sims group.  ::

        sage: H = graphs.HigmanSimsGraph()
        sage: G = H.automorphism_group()
        sage: g=G.order(); g
        88704000
        sage: K = G.normal_subgroups()[1]
        sage: K.is_simple()
        True
        sage: g//K.order()
        2

    REFERENCES:

        .. [BROUWER-HS-2009] `Higman-Sims graph
           <http://www.win.tue.nl/~aeb/graphs/Higman-Sims.html>`_.
           Andries E. Brouwer, accessed 24 October 2009.
        .. [HIGMAN1968] A simple group of order 44,352,000,
           Math.Z. 105 (1968) 110-113. D.G. Higman & C. Sims.
        .. [HAFNER2004] `On the graphs of Hoffman-Singleton and
           Higman-Sims
           <http://www.combinatorics.org/Volume_11/PDF/v11i1r77.pdf>`_.
           The Electronic Journal of Combinatorics 11 (2004), #R77,
           Paul R. Hafner, accessed 24 October 2009.

    AUTHOR:

        - Rob Beezer (2009-10-24)
    """
    HS = graph.Graph()
    HS.name('Higman-Sims graph')

    # Four groups of either five pentagons, or five pentagrams
    # 4 x 5 x 5 = 100 vertices
    # First digit is "group", second is "penta{gon|gram}", third is "vertex"
    vlist = ['%d%d%d'%(g,p,v)
                    for g in range(4) for p in range(5) for v in range(5)]
    for avertex in vlist:
        HS.add_vertex(avertex)

    # Edges: Within groups 0 and 2, joined as pentagons
    # Edges: Within groups 1 and 3, joined as pentagrams
    for g in range(4):
        shift = 1
        if g in [1,3]:
            shift += 1
        for p in range(5):
            for v in range(5):
                HS.add_edge(('%d%d%d'%(g,p,v), '%d%d%d'%(g,p,(v+shift)%5)))

    # Edges: group 0 to group 1
    for x in range(5):
        for m in range(5):
            for c in range(5):
                y = (m*x+c)%5
                HS.add_edge(('0%d%d'%(x,y), '1%d%d'%(m,c)))

    # Edges: group 1 to group 2
    for m in range(5):
        for A in range(5):
            for B in range(5):
                c = (2*(m-A)*(m-A)+B)%5
                HS.add_edge(('1%d%d'%(m,c), '2%d%d'%(A,B)))

    # Edges: group 2 to group 3
    for A in range(5):
        for a in range(5):
            for b in range(5):
                B = (2*A*A+3*a*A-a*a+b)%5
                HS.add_edge(('2%d%d'%(A,B), '3%d%d'%(a,b)))

    # Edges: group 3 to group 0
    for a in range(5):
        for b in range(5):
            for x in range(5):
                y = ((x-a)*(x-a)+b)%5
                HS.add_edge(('3%d%d'%(a,b), '0%d%d'%(x,y)))

    # Edges: group 0 to group 2
    for x in range(5):
        for A in range(5):
            for B in range(5):
                y = (3*x*x+A*x+B+1)%5
                HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))
                y = (3*x*x+A*x+B-1)%5
                HS.add_edge(('0%d%d'%(x,y), '2%d%d'%(A,B)))

    # Edges: group 1 to group 3
    for m in range(5):
        for a in range(5):
            for b in range(5):
                c = (m*(m-a)+b+2)%5
                HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))
                c = (m*(m-a)+b-2)%5
                HS.add_edge(('1%d%d'%(m,c), '3%d%d'%(a,b)))

    # Rename to integer vertex labels, creating dictionary
    # Or not, and create identity mapping
    if relabel:
        vmap = HS.relabel(return_map=True)
    else:
        vmap={}
        for v in vlist:
            vmap[v] = v
    # Layout vertices in a circle
    # In the order given in vlist
    # Using labels from vmap
    pos_dict = {}
    for i in range(100):
        x = float(cos((pi/2) + ((2*pi)/100)*i))
        y = float(sin((pi/2) + ((2*pi)/100)*i))
        pos_dict[vmap[vlist[i]]] = (x,y)
    HS.set_pos(pos_dict)
    return HS

def HoffmanSingletonGraph(self):
    r"""
    Returns the Hoffman-Singleton graph.

    The Hoffman-Singleton graph is the Moore graph of degree 7,
    diameter 2 and girth 5. The Hoffman-Singleton theorem states that
    any Moore graph with girth 5 must have degree 2, 3, 7 or 57. The
    first three respectively are the pentagon, the Petersen graph, and
    the Hoffman-Singleton graph. The existence of a Moore graph with
    girth 5 and degree 57 is still open.

    A Moore graph is a graph with diameter `d` and girth
    `2d + 1`. This implies that the graph is regular, and
    distance regular.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. A novel algorithm written by
    Tom Boothby gives a random layout which is pleasing to the eye.

    REFERENCES:

    - [1] Godsil, C. and Royle, G. Algebraic Graph Theory.
      Springer, 2001.

    EXAMPLES::

        sage: HS = graphs.HoffmanSingletonGraph()
        sage: Set(HS.degree())
        {7}
        sage: HS.girth()
        5
        sage: HS.diameter()
        2
        sage: HS.num_verts()
        50

    Note that you get a different layout each time you create the graph.
    ::

        sage: HS.layout()[1]
        (-0.844..., 0.535...)
        sage: graphs.HoffmanSingletonGraph().layout()[1]
        (-0.904..., 0.425...)

    """
    H = graph.Graph({ \
    'q00':['q01'], 'q01':['q02'], 'q02':['q03'], 'q03':['q04'], 'q04':['q00'], \
    'q10':['q11'], 'q11':['q12'], 'q12':['q13'], 'q13':['q14'], 'q14':['q10'], \
    'q20':['q21'], 'q21':['q22'], 'q22':['q23'], 'q23':['q24'], 'q24':['q20'], \
    'q30':['q31'], 'q31':['q32'], 'q32':['q33'], 'q33':['q34'], 'q34':['q30'], \
    'q40':['q41'], 'q41':['q42'], 'q42':['q43'], 'q43':['q44'], 'q44':['q40'], \
    'p00':['p02'], 'p02':['p04'], 'p04':['p01'], 'p01':['p03'], 'p03':['p00'], \
    'p10':['p12'], 'p12':['p14'], 'p14':['p11'], 'p11':['p13'], 'p13':['p10'], \
    'p20':['p22'], 'p22':['p24'], 'p24':['p21'], 'p21':['p23'], 'p23':['p20'], \
    'p30':['p32'], 'p32':['p34'], 'p34':['p31'], 'p31':['p33'], 'p33':['p30'], \
    'p40':['p42'], 'p42':['p44'], 'p44':['p41'], 'p41':['p43'], 'p43':['p40']})
    for j in range(5):
        for i in range(5):
            for k in range(5):
                con = (i+j*k)%5
                H.add_edge(('q%d%d'%(k,con),'p%d%d'%(j,i)))
    H.name('Hoffman-Singleton graph')
    from sage.combinat.combinat import permutations
    from sage.misc.prandom import randint
    P = permutations([1,2,3,4])
    qpp = [0]+P[randint(0,23)]
    ppp = [0]+P[randint(0,23)]
    qcycle = lambda i,s : ['q%s%s'%(i,(j+s)%5) for j in qpp]
    pcycle = lambda i,s : ['p%s%s'%(i,(j+s)%5) for j in ppp]
    l = 0
    s = 0
    D = []
    while l < 5:
        for q in qcycle(l,s):
            D.append(q)
        vv = 'p%s'%q[1]
        s = int([v[-1] for v in H.neighbors(q) if v[:2] == vv][0])
        for p in pcycle(l,s):
            D.append(p)
        vv = 'q%s'%(int(p[1])+1)
        v = [v[-1] for v in H.neighbors(p) if v[:2] == vv]
        if len(v):
            s = int(v[0])
        l+=1
    map = H.relabel(return_map=True)
    pos_dict = {}
    for i in range(50):
        x = float(cos((pi/2) + ((2*pi)/50)*i))
        y = float(sin((pi/2) + ((2*pi)/50)*i))
        pos_dict[map[D[i]]] = (x,y)
    H.set_pos(pos_dict)
    return H

def HoffmanGraph(self):
    r"""
    Returns the Hoffman Graph.

    See the :wikipedia:`Wikipedia page on the Hoffman graph
    <Hoffman_graph>`.

    EXAMPLES::

        sage: g = graphs.HoffmanGraph()
        sage: g.is_bipartite()
        True
        sage: g.is_hamiltonian() # long time
        True
        sage: g.radius()
        3
        sage: g.diameter()
        4
        sage: g.automorphism_group().cardinality()
        48
    """
    g = graph.Graph({
            0: [1, 7, 8, 13],
            1: [2, 9, 14],
            2: [3, 8, 10],
            3: [4, 9, 15],
            4: [5, 10, 11],
            5: [6, 12, 14],
            6: [7, 11, 13],
            7: [12, 15],
            8: [12, 14],
            9: [11, 13],
            10: [12, 15],
            11: [14],
            13: [15]})
    g.set_pos({})
    _circle_embedding(g, range(8))
    _circle_embedding(g, range(8, 14), radius=.7, shift=.5)
    _circle_embedding(g, [14, 15], radius=.1)

    g.name("Hoffman Graph")

    return g

def HoltGraph(self):
    r"""
    Returns the Holt graph (also called the Doyle graph)

    See the :wikipedia:`Wikipedia page on the Holt graph
    <Holt_graph>`.

    EXAMPLES::

        sage: g = graphs.HoltGraph();g
        Holt graph: Graph on 27 vertices
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian() # long time
        True
        sage: g.radius()
        3
        sage: g.diameter()
        3
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()
        54
    """
    g = graph.Graph(loops=False, name = "Holt graph", pos={})
    for x in range(9):
        for y in range(3):
            g.add_edge((x,y),((4*x+1)%9,(y-1)%3))
            g.add_edge((x,y),((4*x-1)%9,(y-1)%3))
            g.add_edge((x,y),((7*x+7)%9,(y+1)%3))
            g.add_edge((x,y),((7*x-7)%9,(y+1)%3))

    for j in range(0,6,2):
        _line_embedding(g, [(x,j/2) for x in range(9)],
                        first=(cos(2*j*pi/6),sin(2*j*pi/6)),
                        last=(cos(2*(j+1)*pi/6),sin(2*(j+1)*pi/6)))

    return g

def LjubljanaGraph(self, embedding=1):
    r"""
    Returns the Ljubljana Graph.

    The Ljubljana graph is a bipartite 3-regular graph on 112
    vertices and 168 edges. It is not vertex-transitive as it has
    two orbits which are also independent sets of size 56. See the
    :wikipedia:`Wikipedia page on the Ljubljana Graph
    <Ljubljana_graph>`.

    The default embedding is obtained from the Heawood graph.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.LjubljanaGraph()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        10
        sage: g.diameter()
        8
        sage: g.show(figsize=[10, 10])   # long time
        sage: graphs.LjubljanaGraph(embedding=2).show(figsize=[10, 10])   # long time

    TESTS::

        sage: graphs.LjubljanaGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
         -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
         -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
         -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 112, L, 2)
    g.name("Ljubljana graph")

    if embedding == 1:
        return g

    elif embedding == 2:
        dh = HeawoodGraph(self).get_pos()

        # Correspondence between the vertices of the Heawood Graph and
        # 8-sets of the Ljubljana Graph.

        d = {
            0: [1, 21, 39, 57, 51, 77, 95, 107],
            1: [2, 22, 38, 58, 50, 78, 94, 106],
            2: [3, 23, 37, 59, 49, 79, 93, 105],
            3: [4, 24, 36, 60, 48, 80, 92, 104],
            4: [5, 25, 35, 61, 15, 81, 91, 71],
            9: [6, 26, 44, 62, 16, 82, 100, 72],
            10: [7, 27, 45, 63, 17, 83, 101, 73],
            11: [8, 28, 46, 64, 18, 84, 102, 74],
            12: [9, 29, 47, 65, 19, 85, 103, 75],
            13: [10, 30, 0, 66, 20, 86, 56, 76],
            8: [11, 31, 111, 67, 99, 87, 55, 43],
            7: [12, 32, 110, 68, 98, 88, 54, 42],
            6: [13, 33, 109, 69, 97, 89, 53, 41],
            5: [14, 34, 108, 70, 96, 90, 52, 40]
            }

        # The vertices of each 8-set are plotted on a circle, and the
        # circles are slowly shifted to obtain a symmetric drawing.

        for i, (u, vertices) in enumerate(d.iteritems()):
            _circle_embedding(g, vertices, center=dh[u], radius=.1,
                    shift=8.*i/14)

        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def McGeeGraph(self, embedding=2):
    r"""
    Returns the McGee Graph.

    See the :wikipedia:`Wikipedia page on the McGee Graph
    <McGee_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.McGeeGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        7
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.McGeeGraph(embedding=1).show()

    TESTS::

        sage: graphs.McGeeGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
         -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
         -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
         -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 24, [12, 7, -7], 8)
    g.name('McGee graph')

    if embedding == 1:
        return g

    elif embedding == 2:

        o = [[7, 2, 13, 8, 19, 14, 1, 20],
             [5, 4, 11, 10, 17, 16, 23, 22],
             [3, 12, 9, 18, 15, 0, 21, 6]]

        _circle_embedding(g, o[0], radius=1.5)
        _circle_embedding(g, o[1], radius=3, shift=-.5)
        _circle_embedding(g, o[2], radius=2.25, shift=.5)

        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def MoebiusKantorGraph(self):
    """
    Returns a Moebius-Kantor Graph.

    A Moebius-Kantor graph is a cubic symmetric graph. (See also the
    Heawood graph). It has 16 nodes and 24 edges. It is nonplanar and
    Hamiltonian. It has diameter = 4, girth = 6, and chromatic number =
    2. It is identical to the Generalized Petersen graph, P[8,3].

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    REFERENCES:

    - [1] Weisstein, E. (1999). "Moebius-Kantor Graph - from
      Wolfram MathWorld". [Online] Available:
      http://mathworld.wolfram.com/Moebius-KantorGraph.html [2007,
      February 17]

    EXAMPLES::

        sage: MK = graphs.MoebiusKantorGraph()
        sage: MK
        Moebius-Kantor Graph: Graph on 16 vertices
        sage: MK.graph6_string()
        'OhCGKE?O@?ACAC@I?Q_AS'
        sage: (graphs.MoebiusKantorGraph()).show() # long time
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G=GeneralizedPetersenGraph(self,8,3)
    G.name("Moebius-Kantor Graph")
    return G

def MoserSpindle(self):
    r"""
    Returns the Moser spindle.

    For more information, see this
    `MathWorld article on the Moser spindle <http://mathworld.wolfram.com/MoserSpindle.html>`_.

    EXAMPLES:

    The Moser spindle is a planar graph having 7 vertices and 11 edges. ::

        sage: G = graphs.MoserSpindle(); G
        Moser spindle: Graph on 7 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        7
        sage: G.size()
        11

    It is a Hamiltonian graph with radius 2, diameter 2, and girth 3. ::

        sage: G.is_hamiltonian()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    The Moser spindle has chromatic number 4 and its automorphism
    group is isomorphic to the dihedral group `D_4`. ::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(4))
        True
    """
    edge_dict = {
        0: [1,4,5,6],
        1: [2,5],
        2: [3,5],
        3: [4,6],
        4: [6]}
    pos_dict = {
        0: [0, 2],
        1: [-1.90211303259031, 0.618033988749895],
        2: [-1.17557050458495, -1.61803398874989],
        3: [1.17557050458495, -1.61803398874989],
        4: [1.90211303259031, 0.618033988749895],
        5: [1, 0],
        6: [-1, 0]}
    return graph.Graph(edge_dict, pos=pos_dict, name="Moser spindle")


def NauruGraph(self, embedding=2):
    """
    Returns the Nauru Graph.

    See the :wikipedia:`Wikipedia page on the Nauru Graph
    <Nauru_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.NauruGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        6
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.NauruGraph(embedding=1).show()

    TESTS::

        sage: graphs.NauruGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
        sage: graphs.NauruGraph(embedding=1).is_isomorphic(g)
        True
    """

    if embedding == 1:
        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(self, 24, [5, -9, 7, -7, 9, -5], 4)
        g.name('Nauru Graph')
        return g
    elif embedding == 2:
        from sage.graphs.generators.families import GeneralizedPetersenGraph
        g = GeneralizedPetersenGraph(self, 12, 5)
        g.name("Nauru Graph")
        return g
    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def PappusGraph(self):
    """
    Returns the Pappus graph, a graph on 18 vertices.

    The Pappus graph is cubic, symmetric, and distance-regular.

    EXAMPLES::

        sage: G = graphs.PappusGraph()
        sage: G.show()  # long time
        sage: L = graphs.LCFGraph(18, [5,7,-7,7,-7,-5], 3)
        sage: L.show()  # long time
        sage: G.is_isomorphic(L)
        True
    """
    pos_dict = {}
    for i in range(6):
        pos_dict[i] = [float(cos(pi/2 + ((2*pi)/6)*i)),\
                       float(sin(pi/2 + ((2*pi)/6)*i))]
        pos_dict[6 + i] = [(2/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                           (2/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
        pos_dict[12 + i] = [(1/3.0)*float(cos(pi/2 + ((2*pi)/6)*i)),\
                            (1/3.0)*float(sin(pi/2 + ((2*pi)/6)*i))]
    return graph.Graph({0:[1,5,6],1:[2,7],2:[3,8],3:[4,9],4:[5,10],\
                        5:[11],6:[13,17],7:[12,14],8:[13,15],9:[14,16],\
                        10:[15,17],11:[12,16],12:[15],13:[16],14:[17]},\
                       pos=pos_dict, name="Pappus Graph")

def PetersenGraph(self):
    """
    The Petersen Graph is a named graph that consists of 10 vertices
    and 15 edges, usually drawn as a five-point star embedded in a
    pentagon.

    The Petersen Graph is a common counterexample. For example, it is
    not Hamiltonian.

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    EXAMPLES: We compare below the Petersen graph with the default
    spring-layout versus a planned position dictionary of [x,y]
    tuples::

        sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7], 3:[2,4,8], 4:[0,3,9], 5:[0,7,8], 6:[1,8,9], 7:[2,5,9], 8:[3,5,6], 9:[4,6,7]})
        sage: petersen_spring.show() # long time
        sage: petersen_database = graphs.PetersenGraph()
        sage: petersen_database.show() # long time
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    P=GeneralizedPetersenGraph(self, 5,2)
    P.name("Petersen graph")
    return P

def ShrikhandeGraph(self):
    """
    Returns the Shrikhande graph.

    For more information, see the `MathWorld article on the Shrikhande graph
    <http://mathworld.wolfram.com/ShrikhandeGraph.html>`_ or the `Wikipedia
    article on the Shrikhande graph
    <http://en.wikipedia.org/wiki/Shrikhande_graph>`_.

    EXAMPLES:

    The Shrikhande graph was defined by S. S. Shrikhande in 1959. It has
    `16` vertices and `48` edges, and is strongly regular of degree `6` with
    parameters `(2,2)`::

        sage: G = graphs.ShrikhandeGraph(); G
        Shrikhande graph: Graph on 16 vertices
        sage: G.order()
        16
        sage: G.size()
        48
        sage: G.is_regular(6)
        True
        sage: set([ len([x for x in G.neighbors(i) if x in G.neighbors(j)])
        ...     for i in range(G.order())
        ...     for j in range(i) ])
        set([2])

    It is non-planar, and both Hamiltonian and Eulerian::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()
        True
        sage: G.is_eulerian()
        True

    It has radius `2`, diameter `2`, and girth `3`::

        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is `4` and its automorphism group is of order
    `192`::

        sage: G.chromatic_number()
        4
        sage: G.automorphism_group().cardinality()
        192

    It is an integral graph since it has only integral eigenvalues::

        sage: G.characteristic_polynomial().factor()
        (x - 6) * (x - 2)^6 * (x + 2)^9

    It is a toroidal graph, and its embedding on a torus is dual to an
    embedding of the Dyck graph (:meth:`DyckGraph <GraphGenerators.DyckGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i] = [0.5 * pos_dict[i][0],
                           0.5 * pos_dict[i][1]]
    edge_dict = {
        0O00: [0O06, 0O07, 0O01, 0O02,   0O11, 0O17],
        0O01: [0O07, 0O00, 0O02, 0O03,   0O12, 0O10],
        0O02: [0O00, 0O01, 0O03, 0O04,   0O13, 0O11],
        0O03: [0O01, 0O02, 0O04, 0O05,   0O14, 0O12],
        0O04: [0O02, 0O03, 0O05, 0O06,   0O15, 0O13],
        0O05: [0O03, 0O04, 0O06, 0O07,   0O16, 0O14],
        0O06: [0O04, 0O05, 0O07, 0O00,   0O17, 0O15],
        0O07: [0O05, 0O06, 0O00, 0O01,   0O10, 0O16],

        0O10: [0O12, 0O13, 0O15, 0O16,   0O07, 0O01],
        0O11: [0O13, 0O14, 0O16, 0O17,   0O00, 0O02],
        0O12: [0O14, 0O15, 0O17, 0O10,   0O01, 0O03],
        0O13: [0O15, 0O16, 0O10, 0O11,   0O02, 0O04],
        0O14: [0O16, 0O17, 0O11, 0O12,   0O03, 0O05],
        0O15: [0O17, 0O10, 0O12, 0O13,   0O04, 0O06],
        0O16: [0O10, 0O11, 0O13, 0O14,   0O05, 0O07],
        0O17: [0O11, 0O12, 0O14, 0O15,   0O06, 0O00]
    }

    return graph.Graph(edge_dict, pos=pos_dict, name="Shrikhande graph")

def ThomsenGraph(self):
    """
    Returns the Thomsen Graph.

    The Thomsen Graph is actually a complete bipartite graph with (n1,
    n2) = (3, 3). It is also called the Utility graph.

    PLOTTING: See CompleteBipartiteGraph.

    EXAMPLES::

        sage: T = graphs.ThomsenGraph()
        sage: T
        Thomsen graph: Graph on 6 vertices
        sage: T.graph6_string()
        'EFz_'
        sage: (graphs.ThomsenGraph()).show() # long time
    """
    pos_dict = {0:(-1,1),1:(0,1),2:(1,1),3:(-1,0),4:(0,0),5:(1,0)}
    import networkx
    G = networkx.complete_bipartite_graph(3,3)
    return graph.Graph(G, pos=pos_dict, name="Thomsen graph")

def Tutte12Cage(self):
    r"""
    Returns Tutte's 12-Cage.

    See the :wikipedia:`Wikipedia page on the Tutte 12-Cage
    <Tutte_12-cage>`.

    EXAMPLES::

        sage: g = graphs.Tutte12Cage()
        sage: g.order()
        126
        sage: g.size()
        189
        sage: g.girth()
        12
        sage: g.diameter()
        6
        sage: g.show()
    """
    L = [17, 27, -13, -59, -35, 35, -11, 13, -53, 53, -27, 21, 57, 11,
         -21, -57, 59, -17]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 126, L, 7)
    g.name("Tutte 12-Cage")
    return g

def TutteCoxeterGraph(self, embedding=2):
    r"""
    Returns the Tutte-Coxeter graph.

    See the :wikipedia:`Wikipedia page on the Tutte-Coxeter Graph
    <Tutte-Coxeter_graph>`.

    INPUT:

    - ``embedding`` -- two embeddings are available, and can be
      selected by setting ``embedding`` to 1 or 2.

    EXAMPLES::

        sage: g = graphs.TutteCoxeterGraph()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.girth()
        8
        sage: g.diameter()
        4
        sage: g.show()
        sage: graphs.TutteCoxeterGraph(embedding=1).show()

    TESTS::

        sage: graphs.TutteCoxeterGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: The value of embedding must be 1 or 2.
    """

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 30, [-13, -9, 7, -7, 9, 13], 5)
    g.name("Tutte-Coxeter graph")

    if embedding == 1:
        d = {
            0: [1, 3, 5, 7, 29],
            1: [2, 4, 6, 28, 0],
            2: [8, 18, 26, 22, 12],
            3: [9, 13, 23, 27, 17],
            4: [11, 15, 21, 25, 19],
            5: [10, 14, 24, 20, 16]
            }

        _circle_embedding(g, d[0], center=(-1, 1), radius=.25)
        _circle_embedding(g, d[1], center=(1, 1), radius=.25)
        _circle_embedding(g, d[2], center=(-.8, 0), radius=.25, shift=2.5)
        _circle_embedding(g, d[3], center=(1.2, 0), radius=.25)
        _circle_embedding(g, d[4], center=(-1, -1), radius=.25, shift=2)
        _circle_embedding(g, d[5], center=(1, -1), radius=.25)

        return g

    elif embedding == 2:
        return g

    else:
        raise ValueError("The value of embedding must be 1 or 2.")

def WagnerGraph(self):
    """
    Returns the Wagner Graph.

    See the :wikipedia:`Wikipedia page on the Wagner Graph
    <Wagner_graph>`.

    EXAMPLES::

        sage: g = graphs.WagnerGraph()
        sage: g.order()
        8
        sage: g.size()
        12
        sage: g.girth()
        4
        sage: g.diameter()
        2
        sage: g.show()
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(self, 8, [4], 8)
    g.name("Wagner Graph")
    return g

