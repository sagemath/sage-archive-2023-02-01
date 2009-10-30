"""
Graph Theory Cython functions

AUTHORS:
    -- Robert L. Miller   (2007-02-13): initial version
    -- Robert W. Bradshaw (2007-03-31): fast spring layout algorithms
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
from random import random

cdef extern from *:
    double sqrt(double)

def spring_layout_fast_split(G, iterations=50, dim=2, vpos=None, height=False):
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

    EXAMPLE:
        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.graph_fast import spring_layout_fast_split
        sage: spring_layout_fast_split(G)
        {0: [..., ...], ..., 502: [..., ...]}

    AUTHOR:
        Robert Bradshaw
    """
    Gs = G.connected_components_subgraphs()
    pos = {}
    left = 0
    buffer = 1/sqrt(len(G))
    for g in Gs:
        cur_pos = spring_layout_fast(g, iterations, dim, vpos, height)
        xmin = min([x[0] for x in cur_pos.values()])
        xmax = max([x[0] for x in cur_pos.values()])
        if len(g) > 1:
            buffer = (xmax - xmin)/sqrt(len(g))
        for v, loc in cur_pos.iteritems():
            loc[0] += left - xmin + buffer
            pos[v] = loc
        left += xmax - xmin + buffer
    return pos

def spring_layout_fast(G, iterations=50, int dim=2, vpos=None, bint rescale=True, bint height=False):
    """
    Spring force model layout

    This function primarily acts as a wrapper around run_spring,
    converting to and from raw c types.

    This kind of speed cannot be achieved by naive pyrexification of the
    function alone, especially if we require a function call (let alone
    an object creation) every time we want to add a pair of doubles.

    EXAMPLE:
        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.graph_fast import spring_layout_fast
        sage: spring_layout_fast(G)
        {0: [..., ...], ..., 502: [..., ...]}

    """
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
        sage: sage.graphs.graph_fast.binary(389)
        '110000101'
        sage: Integer(389).binary()
        '110000101'
        sage: sage.graphs.graph_fast.binary(2007)
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
        sage: from sage.graphs.graph_fast import R
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
        sage: from sage.graphs.graph_fast import N
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
        sage: from sage.graphs.graph_fast import N_inverse
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
        sage: from sage.graphs.graph_fast import R_inverse
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
        sage: from sage.graphs.graph_fast import D_inverse
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










