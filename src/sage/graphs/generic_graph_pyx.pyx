"""
GenericGraph Cython functions

AUTHORS:
    -- Robert L. Miller   (2007-02-13): initial version
    -- Robert W. Bradshaw (2007-03-31): fast spring layout algorithms
    -- Nathann Cohen                  : exhaustive search
"""

#*****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                         2007 Robert W. Bradshaw <robertwb@math.washington.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include "../ext/interrupt.pxi"
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'

# import from Python standard library
from sage.misc.prandom import random

# import from third-party library
from sage.graphs.base.dense_graph cimport DenseGraph

cdef extern from *:
    double sqrt(double)

cdef class GenericGraph_pyx(SageObject):
    pass

def spring_layout_fast_split(G, **options):
    """
    Graphs each component of G separately, placing them adjacent to
    each other. This is done because on a disconnected graph, the
    spring layout will push components further and further from each
    other without bound, resulting in very tight clumps for each
    component.

    NOTE:
        If the axis are scaled to fit the plot in a square, the
        horizontal distance may end up being "squished" due to
        the several adjacent components.

    EXAMPLES:

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast_split
        sage: spring_layout_fast_split(G)
        {0: [0.452..., 0.247...], ..., 502: [25.7..., 0.505...]}

    AUTHOR:
        Robert Bradshaw
    """
    Gs = G.connected_components_subgraphs()
    pos = {}
    left = 0
    buffer = 1/sqrt(len(G))
    for g in Gs:
        cur_pos = spring_layout_fast(g, **options)
        xmin = min([x[0] for x in cur_pos.values()])
        xmax = max([x[0] for x in cur_pos.values()])
        if len(g) > 1:
            buffer = (xmax - xmin)/sqrt(len(g))
        for v, loc in cur_pos.iteritems():
            loc[0] += left - xmin + buffer
            pos[v] = loc
        left += xmax - xmin + buffer
    return pos

def spring_layout_fast(G, iterations=50, int dim=2, vpos=None, bint rescale=True, bint height=False, by_component = False, **options):
    """
    Spring force model layout

    This function primarily acts as a wrapper around run_spring,
    converting to and from raw c types.

    This kind of speed cannot be achieved by naive Cythonification of the
    function alone, especially if we require a function call (let alone
    an object creation) every time we want to add a pair of doubles.

    INPUT:

     - ``by_component`` - a boolean

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: spring_layout_fast(G)
        {0: [-0.0733..., 0.157...], ..., 502: [-0.551..., 0.682...]}

    With ``split=True``, each component of G is layed out separately,
    placing them adjacent to each other. This is done because on a
    disconnected graph, the spring layout will push components further
    and further from each other without bound, resulting in very tight
    clumps for each component.

    NOTE:

        If the axis are scaled to fit the plot in a square, the
        horizontal distance may end up being "squished" due to
        the several adjacent components.

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: spring_layout_fast(G, by_component = True)
        {0: [2.12..., -0.321...], ..., 502: [26.0..., -0.812...]}
    """

    if by_component:
        return spring_layout_fast_split(G, iterations=iterations, dim = dim,
                                        vpos = vpos, rescale = rescale, height = height,
                                        **options)

    G = G.to_undirected()
    vlist = G.vertices() # this defines a consistent order

    cdef int i, j, x
    cdef int n = G.order()
    if n == 0:
        return {}

    cdef double* pos = <double*>sage_malloc(n * dim * sizeof(double))
    if pos is NULL:
            raise MemoryError, "error allocating scratch space for spring layout"

    # convert or create the starting positions as a flat list of doubles
    if vpos is None:  # set the initial positions randomly in 1x1 box
        for i from 0 <= i < n*dim:
            pos[i] = random()
    else:
        for i from 0 <= i < n:
            loc = vpos[vlist[i]]
            for x from 0 <= x <dim:
                pos[i*dim + x] = loc[x]


    # here we construct a lexicographically ordered list of all edges
    # where elist[2*i], elist[2*i+1] represents the i-th edge
    cdef int* elist = <int*>sage_malloc( (2 * len(G.edges()) + 2) * sizeof(int)  )
    if elist is NULL:
        sage_free(pos)
        raise MemoryError, "error allocating scratch space for spring layout"

    cdef int cur_edge = 0

    for i from 0 <= i < n:
        for j from i < j < n:
            if G.has_edge(vlist[i], vlist[j]):
                elist[cur_edge] = i
                elist[cur_edge+1] = j
                cur_edge += 2

    # finish the list with -1, -1 which never gets matched
    # but does get compared against when looking for the "next" edge
    elist[cur_edge] = -1
    elist[cur_edge+1] = -1

    run_spring(iterations, dim, pos, elist, n, height)

    # recenter
    cdef double* cen
    cdef double r, r2, max_r2 = 0
    if rescale:
        cen = <double *>sage_malloc(sizeof(double) * dim)
        if cen is NULL:
            sage_free(elist)
            sage_free(pos)
            raise MemoryError, "error allocating scratch space for spring layout"
        for x from 0 <= x < dim: cen[x] = 0
        for i from 0 <= i < n:
            for x from 0 <= x < dim:
                cen[x] += pos[i*dim + x]
        for x from 0 <= x < dim: cen[x] /= n
        for i from 0 <= i < n:
            r2 = 0
            for x from 0 <= x < dim:
                pos[i*dim + x] -= cen[x]
                r2 += pos[i*dim + x] * pos[i*dim + x]
            if r2 > max_r2:
                max_r2 = r2
        r = 1 if max_r2 == 0 else sqrt(max_r2)
        for i from 0 <= i < n:
            for x from 0 <= x < dim:
                pos[i*dim + x] /= r
        sage_free(cen)

    # put the data back into a position dictionary
    vpos = {}
    for i from 0 <= i < n:
        vpos[vlist[i]] = [pos[i*dim+x] for x from 0 <= x < dim]

    sage_free(pos)
    sage_free(elist)

    return vpos


cdef run_spring(int iterations, int dim, double* pos, int* edges, int n, bint height):
    """
    Find a locally optimal layout for this graph, according to the
    constraints that neighboring nodes want to be a fixed distance
    from each other, and non-neighboring nodes always repel.

    This is not a true physical model of mutually-repulsive particles
    with springs, rather it is more a model of such things traveling,
    without any inertia, through an (ever thickening) fluid.

    TODO: The inertial model could be incorporated (with F=ma)
    TODO: Are the hard-coded constants here optimal?

    INPUT:
        iterations -- number of steps to take
        dim        -- number of dimensions of freedom
        pos        -- already initialized initial positions
                      Each vertex is stored as [dim] consecutive doubles.
                      These doubles are then placed consecutively in the array.
                      For example, if dim=3, we would have
                      pos = [x_1, y_1, z_1, x_2, y_2, z_2, ... , x_n, y_n, z_n]
        edges      -- List of edges, sorted lexicographically by the first
                      (smallest) vertex, terminated by -1, -1.
                      The first two entries represent the first edge, and so on.
        n          -- number of vertices in the graph
        height     -- if True, do not update the last coordinate ever

    OUTPUT:
        Modifies contents of pos.

    AUTHOR:
        Robert Bradshaw
    """

    cdef int cur_iter, cur_edge
    cdef int i, j, x

    # k -- the equilibrium distance between two adjacent nodes
    cdef double t = 1, dt = t/(1e-20 + iterations), k = sqrt(1.0/n)
    cdef double square_dist, force, scale
    cdef double* disp_i
    cdef double* disp_j
    cdef double* delta

    cdef double* disp = <double*>sage_malloc((n+1) * dim * sizeof(double))
    if disp is NULL:
            raise MemoryError, "error allocating scratch space for spring layout"
    delta = &disp[n*dim]

    if height:
        update_dim = dim-1
    else:
        update_dim = dim

    _sig_on

    for cur_iter from 0 <= cur_iter < iterations:
      cur_edge = 1 # offset by one for fast checking against 2nd element first
      # zero out the disp vectors
      memset(disp, 0, n * dim * sizeof(double))
      for i from 0 <= i < n:
          disp_i = disp + (i*dim)
          for j from i < j < n:
              disp_j = disp + (j*dim)

              for x from 0 <= x < dim:
                  delta[x] = pos[i*dim+x] - pos[j*dim+x]

              square_dist = delta[0] * delta[0]
              for x from 1 <= x < dim:
                  square_dist += delta[x] * delta[x]

              if square_dist < 0.01:
                  square_dist = 0.01

              # they repel according to the (capped) inverse square law
              force = k*k/square_dist

              # and if they are neighbors, attract according Hooke's law
              if edges[cur_edge] == j and edges[cur_edge-1] == i:
                  force -= sqrt(square_dist)/k
                  cur_edge += 2

              # add this factor into each of the involved points
              for x from 0 <= x < dim:
                  disp_i[x] += delta[x] * force
                  disp_j[x] -= delta[x] * force

      # now update the positions
      for i from 0 <= i < n:
          disp_i = disp + (i*dim)

          square_dist = disp_i[0] * disp_i[0]
          for x from 1 <= x < dim:
              square_dist += disp_i[x] * disp_i[x]

          scale = t / (1 if square_dist < 0.01 else sqrt(square_dist))

          for x from 0 <= x < update_dim:
              pos[i*dim+x] += disp_i[x] * scale

      t -= dt

    _sig_off

    sage_free(disp)

def binary(n, length=None):
    """
    A quick python int to binary string conversion.

    EXAMPLE:
        sage: sage.graphs.generic_graph_pyx.binary(389)
        '110000101'
        sage: Integer(389).binary()
        '110000101'
        sage: sage.graphs.generic_graph_pyx.binary(2007)
        '11111010111'
    """
    cdef mpz_t i
    mpz_init(i)
    mpz_set_ui(i,n)
    cdef char* s=mpz_get_str(NULL, 2, i)
    t=str(s)
    free(s)
    mpz_clear(i)
    return t

def R(x):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import R
        sage: R('110111010110110010111000001100000001000000001')
        'vUqwK@?G'

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    # pad on the right to make a multiple of 6
    x += '0' * ( (6 - (len(x) % 6)) % 6)

    # split into groups of 6, and convert numbers to decimal, adding 63
    six_bits = ''
    cdef int i
    for i from 0 <= i < len(x)/6:
        six_bits += chr( int( x[6*i:6*(i+1)], 2) + 63 )
    return six_bits

def N(n):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import N
        sage: N(13)
        'L'
        sage: N(136)
        '~?AG'

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if n < 63:
        return chr(n + 63)
    else:
        # get 18-bit rep of n
        n = binary(n)
        n = '0'*(18-len(n)) + n
        return chr(126) + R(n)

def N_inverse(s):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import N_inverse
        sage: N_inverse('~??~?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?')
        (63, '?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?')
        sage: N_inverse('_???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????')
        (32, '???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????')

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if s[0] == chr(126): # first four bytes are N
        a = binary(ord(s[1]) - 63).zfill(6)
        b = binary(ord(s[2]) - 63).zfill(6)
        c = binary(ord(s[3]) - 63).zfill(6)
        n = int(a + b + c,2)
        s = s[4:]
    else: # only first byte is N
        o = ord(s[0])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        n = o - 63
        s = s[1:]
    return n, s

def R_inverse(s, n):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import R_inverse
        sage: R_inverse('?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?', 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: R_inverse('???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????', 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'

    """
    l = []
    cdef int i
    for i from 0 <= i < len(s):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        a = binary(o-63)
        l.append( '0'*(6-len(a)) + a )
    m = "".join(l)
    return m

def D_inverse(s, n):
    """
    A helper function for the dig6 format.

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import D_inverse
        sage: D_inverse('?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?', 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: D_inverse('???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????', 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'

    """
    l = []
    cdef int i
    for i from 0 <= i < len(s):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        a = binary(o-63)
        l.append( '0'*(6-len(a)) + a )
    m = "".join(l)
    return m[:n*n]


# Exhaustive search in graphs

cpdef list subgraph_search(G, H, bint induced=False):
    r"""
    Returns a set of vertices in ``G`` representing a copy of ``H``.

    ALGORITHM:

    This algorithm is a brute-force search.
    Let `V(H) = \{h_1,\dots,h_k\}`.  It first tries
    to find in `G` a possible representant of `h_1`, then a
    representant of `h_2` compatible with `h_1`, then
    a representant of `h_3` compatible with the first
    two, etc.

    This way, most of the time we need to test far less than
    `k! \binom{|V(G)|}{k}` subsets, and hope this brute-force
    technique can sometimes be useful.

    INPUT:

    - ``G``, ``H`` -- two graphs such that ``H`` is a subgraph of ``G``.

    - ``induced`` -- boolean (default: ``False``); whether to require that
      the subgraph is an induced subgraph.

    OUTPUT:

    A list of vertices inducing a copy of ``H`` in ``G``. If none is found,
    an empty list is returned.

    EXAMPLES:

    A Petersen graph contains an induced path graph `P_5`::

        sage: from sage.graphs.generic_graph_pyx import subgraph_search
        sage: g = graphs.PetersenGraph()
        sage: subgraph_search(g, graphs.PathGraph(5), induced=True)
        [0, 1, 2, 3, 8]

    It also contains a the claw `K_{1,3}`::

        sage: subgraph_search(g, graphs.ClawGraph())
        [0, 1, 4, 5]

    Though it contains no induced `P_6`::

        sage: subgraph_search(g, graphs.PathGraph(6), induced=True)
        []

    TESTS:

    Let `G` and `H` be graphs having orders `m` and `n`, respectively. If
    `m < n`, then there are no copies of `H` in `G`::

        sage: from sage.graphs.generic_graph_pyx import subgraph_search
        sage: m = randint(100, 200)
        sage: n = randint(m + 1, 300)
        sage: G = graphs.RandomGNP(m, random())
        sage: H = graphs.RandomGNP(n, random())
        sage: G.order() < H.order()
        True
        sage: subgraph_search(G, H)
        []
    """
    # TODO: This is a brute-force search and can be very inefficient. Write
    # a more efficient subgraph search implementation.
    cdef int ng = G.order()
    cdef int nh = H.order()
    if ng < nh:
        return []
    cdef int i, j, k
    cdef int *tmp_array
    cdef (bint) (*is_admissible) (int, int *, int *)
    if induced:
        is_admissible = vectors_equal
    else:
        is_admissible = vectors_inferior
    # static copies of the two graphs for more efficient operations
    cdef DenseGraph g = DenseGraph(ng)
    cdef DenseGraph h = DenseGraph(nh)
    # copying the adjacency relations in both G and H
    i = 0
    for row in G.adjacency_matrix():
        j = 0
        for k in row:
            if k:
                g.add_arc(i, j)
            j += 1
        i += 1
    i = 0
    for row in H.adjacency_matrix():
        j = 0
        for k in row:
            if k:
                h.add_arc(i, j)
            j += 1
        i += 1
    # A vertex is said to be busy if it is already part of the partial copy
    # of H in G.
    cdef int *busy = <int *>sage_malloc(ng * sizeof(int))
    memset(busy, 0, ng * sizeof(int))
    # 0 is the first vertex we use, so it is at first busy
    busy[0] = 1
    # stack -- list of the vertices which are part of the partial copy of H
    # in G.
    #
    # stack[i] -- the integer corresponding to the vertex of G representing
    # the i-th vertex of H.
    #
    # stack[i] = -1 means that i is not represented
    # ... yet!
    cdef int *stack = <int *>sage_malloc(nh * sizeof(int))
    stack[0] = 0
    stack[1] = -1
    # Number of representants we have already found. Set to 1 as vertex 0
    # is already part of the partial copy of H in G.
    cdef int active = 1
    # vertices is equal to range(nh), as an int *variable
    cdef int *vertices = <int *>sage_malloc(nh * sizeof(int))
    for 0 <= i < nh:
        vertices[i] = i
    # line_h[i] represents the adjacency sequence of vertex i
    # in h relative to vertices 0, 1, ..., i-1
    cdef int **line_h = <int **>sage_malloc(nh * sizeof(int *))
    for 0 <= i < nh:
        line_h[i] = <int *>h.adjacency_sequence(i, vertices, i)
    # the sequence of vertices to be returned
    cdef list value = []

    _sig_on

    # as long as there is a non-void partial copy of H in G
    while active:
        # If we are here and found nothing yet, we try the next possible
        # vertex as a representant of the active i-th vertex of H.
        i = stack[active] + 1
        # Looking for a vertex that is not busy and compatible with the
        # partial copy we have of H.
        while i < ng:
            if busy[i]:
                i += 1
            else:
                tmp_array = g.adjacency_sequence(active, stack, i)
                if is_admissible(active, tmp_array, line_h[active]):
                    sage_free(tmp_array)
                    break
                else:
                    sage_free(tmp_array)
                    i += 1
        # If we found none, it means that we cannot extend the current copy
        # of H so we update the status of stack[active] and prepare to change
        # the previous vertex.
        if i >= ng:
            if stack[active] != -1:
                busy[stack[active]] = 0
            stack[active] = -1
            active -= 1
        # If we have found a good representant of H's i-th vertex in G
        else:
            if stack[active] != -1:
                busy[stack[active]] = 0
            stack[active] = i
            busy[stack[active]] = 1
            active += 1
            # We have found our copy!!!
            if active == nh:
                g_vertices = G.vertices()
                value = [g_vertices[stack[i]] for i in xrange(nh)]
                break
            else:
                # we begin the search of the next vertex at 0
                stack[active] = -1

    _sig_off

    # Free the memory
    sage_free(busy)
    sage_free(stack)
    sage_free(vertices)
    for 0 <= i < nh:
        sage_free(line_h[i])
    sage_free(line_h)

    return value

cdef inline bint vectors_equal(int n, int *a, int *b):
    r"""
    Tests whether the two given vectors are equal. Two integer vectors
    `a = (a_1, a_2, \dots, a_n)` and `b = (b_1, b_2, \dots, b_n)` are equal
    iff `a_i = b_i` for all `i = 1, 2, \dots, n`. See the function
    ``_test_vectors_equal_inferior()`` for unit tests.

    INPUT:

    - ``n`` -- positive integer; length of the vectors.

    - ``a``, ``b`` -- two vectors of integers.

    OUTPUT:

    - ``True`` if ``a`` and ``b`` are the same vector; ``False`` otherwise.
    """
    cdef int i = 0
    for 0 <= i < n:
        if a[i] != b[i]:
            return False
    return True

cdef inline bint vectors_inferior(int n, int *a, int *b):
    r"""
    Tests whether the second vector of integers is inferior to the first. Let
    `u = (u_1, u_2, \dots, u_k)` and `v = (v_1, v_2, \dots, v_k)` be two
    integer vectors of equal length. Then `u` is said to be less than
    (or inferior to) `v` if `u_i \leq v_i` for all `i = 1, 2, \dots, k`. See
    the function ``_test_vectors_equal_inferior()`` for unit tests. Given two
    equal integer vectors `u` and `v`, `u` is inferior to `v` and vice versa.
    We could also define two vectors `a` and `b` to be equal if `a` is
    inferior to `b` and `b` is inferior to `a`.

    INPUT:

    - ``n`` -- positive integer; length of the vectors.

    - ``a``, ``b`` -- two vectors of integers.

    OUTPUT:

    - ``True`` if ``b`` is inferior to (or less than) ``a``; ``False``
      otherwise.
    """
    cdef int i = 0
    for 0 <= i < n:
        if a[i] < b[i]:
            return False
    return True

##############################
# Further tests. Unit tests for methods, functions, classes defined with cdef.
##############################

def _test_vectors_equal_inferior():
    """
    Unit testing the function ``vectors_equal()``. No output means that no
    errors were found in the random tests.

    TESTS::

        sage: from sage.graphs.generic_graph_pyx import _test_vectors_equal_inferior
        sage: _test_vectors_equal_inferior()
    """
    from sage.misc.prandom import randint
    n = randint(500, 10**3)
    cdef int *u = <int *>sage_malloc(n * sizeof(int))
    cdef int *v = <int *>sage_malloc(n * sizeof(int))
    cdef int i
    # equal vectors: u = v
    for 0 <= i < n:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    try:
        assert vectors_equal(n, u, v)
        assert vectors_equal(n, v, u)
        # Since u and v are equal vectors, then u is inferior to v and v is
        # inferior to u. One could also define u and v as being equal if
        # u is inferior to v and vice versa.
        assert vectors_inferior(n, u, v)
        assert vectors_inferior(n, v, u)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("Vectors u and v should be equal.")
    # Different vectors: u != v because we have u_j > v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j > v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} < v_{n-1}.
    cdef int j = randint(1, n//2)
    cdef int k
    for 0 <= i < j:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    u[j] = randint(-10**6, 10**6)
    v[j] = u[j] - randint(1, 10**6)
    for j < k < n:
        u[k] = randint(-10**6, 10**6)
        v[k] = randint(-10**6, 10**6)
    u[n - 1] = v[n - 1] - randint(1, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        # u is not inferior to v because at least u_j > v_j
        assert u[j] > v[j]
        assert not vectors_inferior(n, v, u)
        # v is not inferior to u because at least v_{n-1} > u_{n-1}
        assert v[n - 1] > u[n - 1]
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("".join([
                    "Vectors u and v should not be equal. ",
                    "u should not be inferior to v, and vice versa."]))
    # Different vectors: u != v because we have u_j < v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j < v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} > v_{n-1}.
    j = randint(1, n//2)
    for 0 <= i < j:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    u[j] = randint(-10**6, 10**6)
    v[j] = u[j] + randint(1, 10**6)
    for j < k < n:
        u[k] = randint(-10**6, 10**6)
        v[k] = randint(-10**6, 10**6)
    u[n - 1] = v[n - 1] + randint(1, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        # u is not inferior to v because at least u_{n-1} > v_{n-1}
        assert u[n - 1] > v[n - 1]
        assert not vectors_inferior(n, v, u)
        # v is not inferior to u because at least u_j < v_j
        assert u[j] < v[j]
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("".join([
                    "Vectors u and v should not be equal. ",
                    "u should not be inferior to v, and vice versa."]))
    # different vectors u != v
    # What's the probability of two random vectors being equal?
    for 0 <= i < n:
        u[i] = randint(-10**6, 10**6)
        v[i] = randint(-10**6, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("Vectors u and v should not be equal.")
    # u is inferior to v, but v is not inferior to u
    for 0 <= i < n:
        v[i] = randint(-10**6, 10**6)
        u[i] = randint(-10**6, 10**6)
        while u[i] > v[i]:
            u[i] = randint(-10**6, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        assert vectors_inferior(n, v, u)
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        raise AssertionError(
            "u should be inferior to v, but v is not inferior to u.")
    finally:
        sage_free(u)
        sage_free(v)
