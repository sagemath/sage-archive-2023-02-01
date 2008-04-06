"""
Chromatic Polynomial Routine

AUTHORS:
    -- Gordon Royle - original C implementation
    -- Robert Miller - transplant

REFERENCE:
    Ronald C Read, An improved method for computing the chromatic polynomials of
    sparse graphs.
"""

#*****************************************************************************
#           Copyright (C) 2008 Robert Miller and Gordon Royle
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.misc.misc import prod
include '../ext/stdsage.pxi'

def chromatic_polynomial(G, return_tree_basis = False):
    """
    Computes the chromatic polynomial of the graph G.

    The algorithm used is a recursive one, based on the following observations
    of Read:
        - The chromatic polynomial of a tree on n vertices is x(x-1)^(n-1).
        - If e is an edge of G, G' is the result of deleting the edge e, and G''
        is the result of contracting e, then the chromatic polynomial of G is
        equal to that of G' minus that of G''.

    EXAMPLES:
        sage: graphs.CycleGraph(4).chromatic_polynomial()
        x^4 - 4*x^3 + 6*x^2 - 3*x
        sage: graphs.CycleGraph(3).chromatic_polynomial()
        x^3 - 3*x^2 + 2*x
        sage: graphs.CubeGraph(3).chromatic_polynomial()
        x^8 - 12*x^7 + 66*x^6 - 214*x^5 + 441*x^4 - 572*x^3 + 423*x^2 - 133*x
        sage: graphs.PetersenGraph().chromatic_polynomial()
        x^10 - 15*x^9 + 105*x^8 - 455*x^7 + 1353*x^6 - 2861*x^5 + 4275*x^4 - 4305*x^3 + 2606*x^2 - 704*x
        sage: graphs.CompleteBipartiteGraph(3,3).chromatic_polynomial()
        x^6 - 9*x^5 + 36*x^4 - 75*x^3 + 78*x^2 - 31*x
        sage: for i in range(2,7):
        ...     graphs.CompleteGraph(i).chromatic_polynomial().factor()
        (x - 1) * x
        (x - 2) * (x - 1) * x
        (x - 3) * (x - 2) * (x - 1) * x
        (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
        (x - 5) * (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
        sage: graphs.CycleGraph(5).chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^2 - 2*x + 2)
        sage: graphs.OctahedralGraph().chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^3 - 9*x^2 + 29*x - 32)
        sage: graphs.WheelGraph(5).chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^2 - 5*x + 7)
        sage: graphs.WheelGraph(6).chromatic_polynomial().factor()
        (x - 3) * (x - 2) * (x - 1) * x * (x^2 - 4*x + 5)

    """
    if not G.is_connected():
        return prod([chromatic_polynomial(g) for g in G.connected_components_subgraphs()])
    R = ZZ['x']
    x = R.gen()
    if G.is_tree():
        return x*(x-1)**(G.num_verts()-1)

    cdef int nverts, nedges, i, j, u, v, top, bot, num_chords, next_v, m, coeff
    cdef int *queue, *chords1, *chords2, *bfs_reorder, *parent, *tot, *coeffs
    G = G.relabel(inplace=False)
    G.remove_multiple_edges()
    G.remove_loops()
    nverts = G.num_verts()
    nedges = G.num_edges()
    queue = <int *> sage_malloc(nverts * sizeof(int))
    chords1 = <int *> sage_malloc((nedges - nverts + 1) * sizeof(int))
    chords2 = <int *> sage_malloc((nedges - nverts + 1) * sizeof(int))
    parent = <int *> sage_malloc(nverts * sizeof(int))
    bfs_reorder = <int *> sage_malloc(nverts * sizeof(int))
    tot = <int *> sage_malloc((nverts+1) * sizeof(int))
    if queue is NULL or \
       chords1 is NULL or \
       chords2 is NULL or \
       parent is NULL or \
       bfs_reorder is NULL or \
       tot is NULL:
        if queue is not NULL:
            sage_free(queue)
        if chords1 is not NULL:
            sage_free(chords1)
        if chords2 is not NULL:
            sage_free(chords2)
        if parent is not NULL:
            sage_free(parent)
        if bfs_reorder is not NULL:
            sage_free(bfs_reorder)
        if tot is not NULL:
            sage_free(tot)
        raise RuntimeError("Error allocating memory for chrompoly.")
    num_chords = 0

    for i from 0 <= i < nverts:
        queue[i] = 1000000000
        parent[i] = 1000000000
        bfs_reorder[i] = 1000000000
        tot[i] = 1000000000
    tot[nverts] = 1000000000
    for i from 0 <= i < nedges - nverts + 1:
        chords1[i] = 1000000000
        chords2[i] = 1000000000

    # Breadth first search from 0:
    bfs_reorder[0] = 0
    tot[0] = 0
    for i from 0 < i < nverts:
        bfs_reorder[i] = -1
        tot[i] = 0
    tot[nverts] = 0
    queue[0] = 0
    top = 1
    bot = 0
    next_v = 1
    while top > bot:
        v = queue[bot]
        bot += 1
        for u in G.neighbor_iterator(v):
            if bfs_reorder[u] == -1: # if u is not yet in tree
                bfs_reorder[u] = next_v
                next_v += 1
                queue[top] = u
                top += 1
                parent[bfs_reorder[u]] = bfs_reorder[v]
            else:
                if bfs_reorder[u] > bfs_reorder[v]:
                    chords1[num_chords] = bfs_reorder[u]
                    chords2[num_chords] = bfs_reorder[v]
                else:
                    continue
                i = num_chords
                num_chords += 1
                # bubble sort the chords
                while i > 0:
                    if chords1[i-1] > chords1[i]:
                        break
                    if chords1[i-1] == chords1[i] and chords2[i-1] > chords2[i]:
                        break
                    j = chords1[i-1]
                    chords1[i-1] = chords1[i]
                    chords1[i] = j
                    j = chords2[i-1]
                    chords2[i-1] = chords2[i]
                    chords2[i] = j
                    i -= 1
    try:
        contract_and_count(chords1, chords2, num_chords, nverts, tot, parent)
    except RuntimeError:
        sage_free(queue)
        sage_free(chords1)
        sage_free(chords2)
        sage_free(parent)
        sage_free(bfs_reorder)
        sage_free(tot)
        raise RuntimeError("Error allocating memory for chrompoly.")
    coeffs = <int *> sage_malloc((nverts+1) * sizeof(int))
    if coeffs is NULL:
        sage_free(queue)
        sage_free(chords1)
        sage_free(chords2)
        sage_free(parent)
        sage_free(bfs_reorder)
        sage_free(tot)
        raise RuntimeError("Error allocating memory for chrompoly.")
    for i from 0 <= i <= nverts:
        coeffs[i] = 0
    m = -1
    # start with the zero polynomial: f(x) = 0
    for i from nverts >= i > 0:
        if tot[i] == 0:
            continue
        m = -m

        # do this:
        # f += tot[i]*m*x*(x-1)**(i-1)
        coeff = 1
        coeffs[i] += tot[i]*m
        for j from 1 <= j < i:
            # an iterative method for binomial coefficients...
            coeff *= (j-i)
            coeff /= j
            coeffs[i-j] += tot[i]*m*coeff
    f = R([coeffs[i] for i from 0 <= i <= nverts])
    sage_free(queue)
    sage_free(chords1)
    sage_free(chords2)
    sage_free(parent)
    sage_free(bfs_reorder)
    sage_free(tot)
    sage_free(coeffs)
    return f

cdef int contract_and_count(int *chords1, int *chords2, int num_chords, int nverts, \
                         int *tot, int *parent):
    if num_chords == 0:
        tot[nverts] += 1
        return 0
    cdef int *new_chords1 = <int *> sage_malloc(num_chords * sizeof(int))
    cdef int *new_chords2 = <int *> sage_malloc(num_chords * sizeof(int))
    cdef int *ins_list1 = <int *> sage_malloc(num_chords * sizeof(int))
    cdef int *ins_list2 = <int *> sage_malloc(num_chords * sizeof(int))
    if new_chords1 is NULL or new_chords2 is NULL or ins_list1 is NULL or ins_list2 is NULL:
        if new_chords1 is not NULL:
            sage_free(new_chords1)
        if new_chords2 is not NULL:
            sage_free(new_chords2)
        if ins_list1 is not NULL:
            sage_free(ins_list1)
        if ins_list2 is not NULL:
            sage_free(ins_list2)
        raise RuntimeError("Error allocating memory for chrompoly.")
    cdef int i, j, k, x1, xj, z, num, insnum, parent_checked
    for i from 0 <= i < num_chords:
        # contract chord i, and recurse
        z = chords1[i]
        x1 = chords2[i]
        j = i + 1
        insnum = 0
        parent_checked = 0
        while j < num_chords and chords1[j] == z:
            xj = chords2[j]
            if parent[z] > xj:
                parent_checked = 1
                # now try adding {x1, parent[z]} to the list
                if not parent[x1] == parent[z]:
                    if x1 > parent[z]:
                        ins_list1[insnum] = x1
                        ins_list2[insnum] = parent[z]
                    else:
                        ins_list1[insnum] = parent[z]
                        ins_list2[insnum] = x1
                    insnum += 1
            if not parent[x1] == xj: # then {x1, xj} isn't already a tree edge
                ins_list1[insnum] = x1
                ins_list2[insnum] = xj
                insnum += 1
            j += 1
        if not parent_checked:
            if not parent[x1] == parent[z]:
                if x1 > parent[z]:
                    ins_list1[insnum] = x1
                    ins_list2[insnum] = parent[z]
                else:
                    ins_list1[insnum] = parent[z]
                    ins_list2[insnum] = x1
                insnum += 1

        # now merge new_chords and ins_list
        num = 0
        k = 0
        while k < insnum and j < num_chords:
            if chords1[j] > ins_list1[k] or \
              (chords1[j] == ins_list1[k] and chords2[j] > ins_list2[k]):
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
            elif chords1[j] < ins_list1[k] or \
              (chords1[j] == ins_list1[k] and chords2[j] < ins_list2[k]):
                new_chords1[num] = ins_list1[k]
                new_chords2[num] = ins_list2[k]
                num += 1
                k += 1
            else:
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
                k += 1
        if j == num_chords:
            while k < insnum:
                new_chords1[num] = ins_list1[k]
                new_chords2[num] = ins_list2[k]
                num += 1
                k += 1
        elif k == insnum:
            while j < num_chords:
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
        try:
            contract_and_count(new_chords1, new_chords2, num, nverts - 1, tot, parent)
        except RuntimeError:
            sage_free(new_chords1)
            sage_free(new_chords2)
            sage_free(ins_list1)
            sage_free(ins_list2)
            raise RuntimeError("Error allocating memory for chrompoly.")
    tot[nverts] += 1
    sage_free(new_chords1)
    sage_free(new_chords2)
    sage_free(ins_list1)
    sage_free(ins_list2)


