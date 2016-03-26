r"""
Permutation in Sage works by default on `{1, 2, ..., n}` but it might be much
more convenient to work on `{0, 1, ..., n-1}`. This module provide simple
functions for the latter representation.
"""

from sage.rings.integer import Integer
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement 
from sage.categories.groups import Groups


def cycles_to_list(t):
    r"""
    Returns a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    EXAMPLES::

        sage: from sage.misc.permutation import cycles_to_list
        sage: cycles_to_list([[1,3,5],[0,2,4],[6]])
        [2, 3, 4, 5, 0, 1, 6]
    """
    res = map(Integer, range(max(map(max, t))+1))

    for c in t:
        for j in xrange(len(c) - 1):
            res[c[j]] = c[j + 1]
        res[c[-1]] = c[0]

    return res


def str_to_cycles(s):
    """
    Returns a list of cycles from a string

    EXAMPLES::

        sage: from sage.misc.permutation import str_to_cycles
        sage: str_to_cycles('(0,1)')
        [[0, 1]]
        sage: str_to_cycles('(0,1)(3,2)')
        [[0, 1], [3, 2]]
    """
    return [map(Integer, c_str.replace(' ', '').split(','))
            for c_str in s[1:-1].split(')(')]


def init_perm(data):
    """
    Return a permutation from different kinds of data.

    INPUT:

    a list of permutations, each one as:

    - an element of a SymmetricGroup

    - or a list of values

    - or a tuple of cycles

    - or a string of cycles

    OUTPUT:

    a list of permutations in a SymmetricGroup

    EXAMPLES::

        sage: from sage.misc.permutation import init_perm
        sage: init_perm([[3,2,0,1]])
        [(0,3,1,2)]

        sage: init_perm([([2,1],[3,4,0])])
        [(0,3,4)(1,2)]

        sage: init_perm(['(0,1)(3,2)'])
        [(0,1)(2,3)]

        sage: S = SymmetricGroup(range(4))
        sage: init_perm([S([3,1,2,0])])
        [(0,3)]
    """
    if isinstance(data[0], list):  # already given as a list of values
        pass
    elif isinstance(data[0], tuple):  # given as a tuple of cycles
        data = map(cycles_to_list, data)
    elif isinstance(data[0], str):  # given as a string of cycles
        data = map(str_to_cycles, data)
        data = map(cycles_to_list, data)
    elif data[0].parent() in Groups:
        return data

    n = max(map(len, data))
    g = SymmetricGroup(range(n))
    for p in data:
        p.extend(xrange(len(p), n))
    # return data
    return [g(u) for u in data]


def perm_orbit(p, i):
    r"""
    Returns the orbit of an integer `i` under the permutation `p`

    EXAMPLES::

        sage: from sage.misc.permutation import perm_orbit
        sage: perm_orbit([0,3,1,2],2)
        [2, 1, 3]
    """
    res = [i]
    j = p[i]
    while j != i:
        res.append(j)
        j = p[j]
    return res


def perm_switch(p1, p2, i, j):
    """
    Exchanges the values at positions `i` and `j` in two permutations
    `p_1` and `p_2`

    EXAMPLES::

        sage: from sage.misc.permutation import perm_switch
        sage: a = [0,1,2]
        sage: b = [1,2,0]
        sage: perm_switch(a,b,0,1)
        sage: a;b
        [1, 0, 2]
        [2, 1, 0]
    """
    i1 = p1[i]
    j1 = p1[j]
    p1[i] = j1
    p1[j] = i1

    i2 = p2[i1]
    j2 = p2[j1]
    p2[i1] = j2
    p2[j1] = i2


def perms_are_connected(g, n):
    """
    Checks that the action of the generated group is transitive

    INPUT:

    - a list of permutations of `[0, n-1]` (in a SymmetricGroup)

    - an integer `n`

    EXAMPLES::

        sage: from sage.misc.permutation import perms_are_connected
        sage: S = SymmetricGroup(range(3))
        sage: perms_are_connected([S([0,1,2]),S([0,2,1])],3)
        False
        sage: perms_are_connected([S([0,1,2]),S([1,2,0])],3)
        True
    """
    from sage.graphs.graph import Graph
    G = Graph()
    if g:
        G.add_vertices(g[0].domain())
    for p in g:
        G.add_edges(enumerate(p.tuple()))
    return G.num_verts() == n and G.is_connected()


def perms_relabel(p, m):
    """
    Relabel the permutations in `p` according to `m`.

    INPUT:

    - `p` is a list of permutations (in a SymmetricGroup)

    - `m` is a permutation (in a SymmetricGroup)

    OUTPUT:

    a list of permutations (in a SymmetricGroup)

    EXAMPLES::

        sage: from sage.misc.permutation import perms_relabel
        sage: S = SymmetricGroup(range(4))
        sage: p0 = [S([0,3,1,2]),S([3,0,2,1]),S([1,2,3,0])]
        sage: p1 = perms_relabel(p0,S([2,1,3,0])); p1
        [(0,3,1), (0,1,2), (0,2,1,3)]

    Applying the inverse permutation we get back our original permutation::

        sage: p2 = perms_relabel(p1,S([3,1,0,2]))
        sage: p2 == p0
        True
    """
    S = p[0].parent()
    q = [list(k.tuple()) for k in p]
    for i in xrange(len(m.domain())):
        for j in xrange(len(p)):
            q[j][m(i)] = m(p[j](i))
    return [S(u) for u in q]


def perms_canonical_labels_from(x, y, j0, verbose=False):
    r"""
    Return canonical labels for ``x``, ``y`` that starts at ``j0``

    .. WARNING:

        The group generated by ``x`` and the elements of ``y`` should be
        transitive.

    INPUT:

    - ``x`` -- list - a permutation of `[0, ..., n]` as a list

    - ``y`` -- list of permutations of `[0, ..., n]` as a list of lists

    - ``j0`` -- an index in [0, ..., n]

    OUTPUT:

    mapping: a permutation that specify the new labels

    EXAMPLES::

        sage: from sage.misc.permutation import perms_canonical_labels_from
        sage: S = SymmetricGroup(range(3))
        sage: perms_canonical_labels_from(S([0,1,2]),[S([1,2,0])],0)
        ()
        sage: perms_canonical_labels_from(S([1,0,2]),[S([2,0,1])],0)
        ()
        sage: perms_canonical_labels_from(S([1,0,2]),[S([2,0,1])],1)
        (0,1)
        sage: perms_canonical_labels_from(S([1,0,2]),[S([2,0,1])],2)
        (0,2)
    """
    n = len(x.domain())
    S = x.parent()

    k = 0
    mapping = [None] * n
    waiting = [[] for i in xrange(len(y))]

    while k < n:
        if verbose:
            print "complete from", j0
        # initialize at j0
        mapping[j0] = k
        waiting[0].append(j0)
        k += 1
        # complete x cycle from j0
        j = x(j0)
        while j != j0:
            mapping[j] = k
            waiting[0].append(j)
            k += 1
            j = x(j)
        if verbose:
            print "completed cycle mapping=", mapping

        # find another guy
        if verbose:
            print "try to find somebody in", waiting
        l = 0
        while l < len(waiting):
            i = 0
            while i < len(waiting[l]):
                j1 = waiting[l][i]
                if mapping[y[l](j1)] is None:
                    break
                i += 1

            if i == len(waiting[l]):  # not found: go further in waiting
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l])
                waiting[l] = []
                l += 1
                i = 0

            else:  # found: complete cycle from new guy
                j0 = y[l](j1)
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l][:i+1])
                del waiting[l][:i+1]
                break

    return S(mapping)


def perms_canonical_labels(p, e=None):
    """
    Relabel a list with a common conjugation such that two conjugated lists are
    relabeled the same way.

    INPUT:

    - ``p`` is a list of at least 2 permutations

    - ``e`` is None or a list of integer in the domain of the permutations. If
      provided, then the renumbering algorithm is only performed from the
      elements of ``e``.

    OUTPUT:

    - a pair made of a list of permutations (as a list of lists) and a list that
      corresponds to the conjugacy used.

    EXAMPLES::

        sage: from sage.misc.permutation import perms_canonical_labels
        sage: S = SymmetricGroup(range(4))
        sage: l0 = [S([2,0,3,1]),S([3,1,2,0]),S([0,2,1,3])]
        sage: l,m = perms_canonical_labels(l0); l
        [(0,1,2,3), (1,3), (0,2)]

        sage: from sage.misc.permutation import perms_relabel
        sage: perms_relabel(l0,m) == l
        True

        sage: from sage.misc.permutation import perms_canonical_labels
        sage: perms_canonical_labels([])
        Traceback (most recent call last):
        ...
        ValueError: input must have length >= 2
    """
    if not len(p) > 1:
        raise ValueError('input must have length >= 2')
    n = len(p[0].domain())

    c_win = None
    m_win = range(n)

    x = p[0]
    y = p[1:]

    if e is None:
        e = range(n)

    # get canonical label from i in to_test and compare
    while e:
        i = e.pop()
        m_test = perms_canonical_labels_from(x, y, i)
        c_test = perms_relabel(p, m_test)
        if c_win is None or c_test < c_win:
            c_win = c_test
            m_win = m_test

    return c_win, m_win
