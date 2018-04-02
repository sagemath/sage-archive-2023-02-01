# -*- coding: utf-8 -*-
r"""
The Hillman-Grassl correspondence

Fix a partition `\lambda`.
We draw all partitions and tableaux in English notation.

A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
entries are nonnegative integers. (No conditions on the order of
these entries are made. Note that `0` is allowed.)

A *reverse plane partition of shape* `\lambda` (short:
`\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
increase along each row and weakly increase along each column.
(Some authors -- like Stanley in [EnumComb2]_ Section 7.22 --
call this a weak reverse plane partition.)

The Hillman-Grassl correspondence `H`
(:meth:`hillman_grassl`) is the map that sends a
`\lambda`-array `M` to a `\lambda`-rpp `H(M)` defined recursively
as follows:

* If all entries of `M` are `0`, then `H(M) = M`.

* Otherwise, let `s` be the index of the leftmost column of `M`
  containing a nonzero entry.
  Let `r` be the index of the bottommost nonzero entry in the
  `s`-th column of `M`.
  Let `M'` be the `\lambda`-array obtained from `M` by removing
  `1` from the `(r, s)`-th entry of `M`.
  Let `Q = (q_{i, j})` be the image `H(M')` (which is already
  defined by recursion).

* Define a sequence `((i_1, j_1), (i_2, j_2), \ldots,
  (i_n, j_n))` of boxes in the diagram of `\lambda` (actually a
  lattice path made of southward and westward steps) as follows:
  Set `(i_1, j_1) = (r, \lambda_r)` (the rightmost box in the
  `r`-th row of `\lambda`).
  If `(i_k, j_k)` is defined for some `k \geq 1`, then
  `(i_{k+1}, j_{k+1})` is constructed as follows:
  If `q_{i_k + 1, j_k}` is well-defined and equals `q_{i_k, j_k}`,
  then we set `(i_{k+1}, j_{k+1}) = (i_k + 1, j_k)`.
  Otherwise, if `j_k = s`, then the sequence ends here.
  Otherwise, we set `(i_{k+1}, j_{k+1}) = (i_k, j_k - 1)`.

* Let `H(M)` be the array obtained from `Q` by adding `1` to
  the `(i_k, j_k)`-th entry of `Q` for each
  `k \in \{1, 2, \ldots, n\}`.

See [Gans1981]_ (Section 3) for this construction.

The inverse `H^{-1}` of the Hillman-Grassl correspondence
(:meth:`hillman_grassl`_inverse) sends
a `\lambda`-rpp `\pi` to a `\lambda`-array `H^{-1}(\pi)` constructed
recursively as follows:

* If all entries of `\pi` are `0`, then `H^{-1}(\pi) = \pi`.

* Otherwise, let `s` be the index of the leftmost column of `\pi`
  containing a nonzero entry.
  Write the `\lambda`-array `M` as `(m_{i, j})`.

* Define a sequence `((i_1, j_1), (i_2, j_2), \ldots,
  (i_n, j_n))` of boxes in the diagram of `\lambda` (actually a
  lattice path made of northward and eastward steps) as follows:
  Let `(i_1, j_1)` be the bottommost box in the `s`-th column
  of `\pi`.
  If `(i_k, j_k)` is defined for some `k \geq 1`, then
  `(i_{k+1}, j_{k+1})` is constructed as follows:
  If `q_{i_k - 1, j_k}` is well-defined and equals `q_{i_k, j_k}`,
  then we set `(i_{k+1}, j_{k+1}) = (i_k - 1, j_k)`.
  Otherwise, we set `(i_{k+1}, j_{k+1}) = (i_k, j_k + 1)` if
  this is still a box of `\lambda`.
  Otherwise, the sequence ends here.

* Let `\pi'` be the `\lambda`-rpp obtained from `\pi` by
  subtracting `1` from the `(i_k, j_k)`-th entry of `\pi` for each
  `k \in \{1, 2, \ldots, n\}`.

* Let `N'` be the image `H^{-1}(\pi')` (which is already
  constructed by recursion).
  Then, `H^{-1}(\pi)` is obtained from `N'` by adding `1` to the
  `(i_n, s)`-th entry of `N'`.

This construction appears in [HilGra1976]_ Section 6 (where
`\lambda`-arrays are re-encoded as sequences of "hook number
multiplicities") and [EnumComb2]_ Section 7.22.

REFERENCES:

.. [Gans1981] Emden R. Gansner,
  *The Hillman-Grassl Correspondence and the
  Enumeration of Reverse Plane Partitions*,
  Journal of Combinatorial Theory, Series A
  30 (1981), pp. 71--89.
  :doi:`10.1016/0097-3165(81)90041-8`

.. [HilGra1976] \A. P. Hillman, R. M. Grassl,
  *Reverse plane partitions and tableau hook numbers*,
  Journal of Combinatorial Theory, Series A 21 (1976),
  pp. 216--221.
  :doi:`10.1016/0097-3165(76)90065-0`

- [EnumComb2]_

.. [Sulzgr2017] Robin Sulzgruber,
  *Inserting rim-hooks into reverse plane partitions*,
  :arxiv:`1710.09695v1`.

.. [Pak2002] Igor Pak,
  *Hook length formula and geometric combinatorics*,
  Seminaire Lotharingien de Combinatoire, 46 (2001),
  B46f,
  https://eudml.org/doc/121696

.. [Hopkins2017] Sam Hopkins,
  *RSK via local transformations*,
  http://web.mit.edu/~shopkins/docs/rsk.pdf

AUTHORS:

- Darij Grinberg and Tom Roby (2018): Initial implementation

.. TODO::

    * Connect this up with the tableaux and plane partitions.
      At the moment, everything is arrays. Probably tuples
      should be preferred, but also this should be a proper
      class? two proper classes? Wait until Tableaux have
      been refactored?

    * Properly document Pak and Sulzgruber correspondences?
"""
#*****************************************************************************
#       Copyright (C) 2018 Darij Grinberg <darijgrinberg@gmail.com>,
#                     2018 Tom Roby <tomrobyuconn@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def transpose(M):
    r"""
    Return the transpose of a `\lambda`-array.
    
    The transpose of a `\lambda`-array `(m_{i, j})` is the
    `\lambda^t`-array `(m_{j, i})`.

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import transpose
        sage: transpose([[1, 2, 3], [4, 5]])
        [[1, 4], [2, 5], [3]]
        sage: transpose([[5, 0, 3], [4, 1, 0], [7]])
        [[5, 4, 7], [0, 1], [3, 0]]

    TESTS::

        sage: transpose(((2, 1), (3,)))
        [[2, 3], [1]]
        sage: transpose([])
        []
    """
    if len(M) == 0:
        return []
    l = len(M[0])
    res = []
    for j in range(l):
        lj = []
        for r in M:
            if len(r) <= j:
                break
            lj.append(r[j])
        res.append(lj)
    return res


def hillman_grassl(M):
    r"""
    Return the image of the `\lambda`-array ``M``
    under the Hillman-Grassl correspondence.

    The Hillman-Grassl correspondence is a bijection
    between the tableaux with nonnegative entries
    (otherwise arbitrary) and the reverse plane
    partitions with nonnegative entries.
    This bijection preserves the shape of the
    tableau. See :mod:`~sage.combinat.hillman_grassl`.

    .. SEEALSO::

        :meth:`hillman_grassl_inverse`.

    EXAMPLES::

        sage: hillman_grassl([[2, 1, 1], [0, 2, 0], [1, 1]])
        [[2, 2, 4], [2, 3, 4], [3, 5]]
        sage: hillman_grassl([[1, 2, 0], [1, 0, 1], [1]])
        [[0, 1, 3], [2, 4, 4], [3]]
        sage: hillman_grassl([])
        []
        sage: hillman_grassl([[3, 1, 2]])
        [[3, 4, 6]]
        sage: hillman_grassl([[2, 2, 0], [1, 1, 1], [1]])
        [[1, 2, 4], [3, 5, 5], [4]]
        sage: hillman_grassl([[1, 1, 1, 1]]*3)
        [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

    TESTS::

        sage: hillman_grassl(((2, 2, 0), (1, 1, 1), (1,)))
        [[1, 2, 4], [3, 5, 5], [4]]
    """
    lam = [len(row) for row in M]
    l = len(lam)
    Mt = transpose(M)
    # We won't touch M from now on; we'll only modify Mt
    # in place.
    hook_mults = []
    for j, col_j in enumerate(Mt):
        col_j_hook_mults = []
        for r, entry in enumerate(col_j):
            if entry != 0:
                col_j_hook_mults += [(r, j)] * entry
        hook_mults += reversed(col_j_hook_mults)
    res = [[0] * rowlen for rowlen in lam]
    for (r, s) in reversed(hook_mults):
        i = r
        j = lam[r] - 1
        while True:
            old = res[i][j]
            res[i][j] += 1
            if i + 1 != l and j < lam[i + 1] and old == res[i + 1][j]:
                i += 1
            else:
                if j == s:
                    break
                j -= 1
    return res

def hillman_grassl_inverse(M):
    r"""
    Return the image of the `\lambda`-rpp ``M`` under the
    inverse of the Hillman-Grassl correspondence.

    See :mod:`~sage.combinat.hillman_grassl`.

    .. SEEALSO::

        :meth:`hillman_grassl`.

    EXAMPLES::

        sage: hillman_grassl_inverse([[2, 2, 4], [2, 3, 4], [3, 5]])
        [[2, 1, 1], [0, 2, 0], [1, 1]]
        sage: hillman_grassl_inverse([[0, 1, 3], [2, 4, 4], [3]])
        [[1, 2, 0], [1, 0, 1], [1]]

    Applying the inverse of the Hillman-Grassl correspondence
    to the transpose of a `\lambda`-rpp `M` yields the same
    result as applying it to `M` and then transposing the
    result ([Gans1981]_ Corollary 3.4)::

        sage: hillman_grassl_inverse([[1,3,5],[2,4]])
        [[1, 2, 2], [1, 1]]
        sage: hillman_grassl_inverse([[1,2],[3,4],[5]])
        [[1, 1], [2, 1], [2]]
        sage: hillman_grassl_inverse([[1, 2, 3], [1, 2, 3], [2, 4, 4]])
        [[1, 2, 0], [0, 1, 1], [1, 0, 1]]
        sage: hillman_grassl_inverse([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]

    TESTS::

        sage: a = [[3], [4], [6]]
        sage: hillman_grassl_inverse(a)
        [[3], [1], [2]]
        sage: a
        [[3], [4], [6]]

        sage: hillman_grassl_inverse(((1,2),(3,4),(5,)))
        [[1, 1], [2, 1], [2]]
    """
    lam = [len(row) for row in M]
    res = [[0] * rowlen for rowlen in lam]
    Mt = transpose(M)
    # We won't touch M from now on; we'll only modify Mt
    # in place.
    while True:
        for j, col_j in enumerate(Mt):
            if all(entry == 0 for entry in col_j):
                continue
            else:
                break
        else: # all entries of Mt are 0.
            break
        # Now, j is the index of the leftmost nonzero column of
        # the array.
        s = j
        i = len(col_j) - 1 # We already have j = s.
        while True:
            old = col_j[i]
            col_j[i] -= 1
            if i > 0 and old == col_j[i - 1]:
                i -= 1
            else:
                j += 1
                if j == lam[i]:
                    break
                col_j = Mt[j]
        res[i][s] += 1
    return res

def pak_correspondence(M):
    r"""
    Return the image of a `\lambda`-array ``M``
    under the Pak correspondence.

    The Pak correspondence is the map `\xi_\lambda`
    from [Sulzgr2017]_ Section 7, and is the map
    `\xi_\lambda^{-1}` from [Pak2002]_ Section 5.
    It is denoted by `\mathcal{RSK}` in [Hopkins2017]_.
    It is the inverse of the Sulzgruber correspondence
    (:meth:`sulzgruber_correspondence`).

    EXAMPLES::

        sage: pak_correspondence([[1, 0, 2], [0, 2, 0], [1, 1, 0]])
        [[1, 2, 3], [1, 2, 3], [2, 4, 4]]
        sage: pak_correspondence([[1, 1, 2], [0, 1, 0], [3, 0, 0]])
        [[1, 1, 4], [2, 3, 4], [4, 4, 4]]
        sage: pak_correspondence([[1, 0, 2], [0, 2, 0], [1, 1]])
        [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: pak_correspondence([[0, 2, 2], [1, 1], [2]])
        [[1, 2, 4], [1, 3], [3]]
        sage: pak_correspondence([[1, 1, 1, 1]]*3)
        [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

    TESTS::

        sage: pak_correspondence([])
        []
        sage: pak_correspondence(((0, 2, 2), (1, 1), (2,)))
        [[1, 2, 4], [1, 3], [3]]
    """
    lam = [len(row) for row in M]
    l = len(lam)
    if l == 0:
        return []
    # Finding a corner of lam...
    lam_0 = lam[0]
    for i, lam_i in enumerate(lam):
        if lam_i != lam_0:
            i -= 1
            break
    j = lam_0 - 1
    # Now, i is the index of the last row of ``M`` that
    # has the same length as the first row; hence, (i, j)
    # is a corner of lam.
    x = M[i][j]
    N = [list(row) for row in M]

    # remove the corner (i, j):
    N[i].pop()
    if not N[i]:
        N.pop()

    N = pak_correspondence(N)

    # toggling and inserting the new entry:
    for k in range(min(i, j) + 1):
        u = i - k
        v = j - k
        if u > 0 and v > 0:
            lower_bound = max(N[u-1][v], N[u][v-1])
        elif u > 0:
            lower_bound = N[u-1][v]
        elif v > 0:
            lower_bound = N[u][v-1]
        else:
            lower_bound = 0
        if k > 0:
            val = N[u][v]
            upper_bound = min(N[u+1][v], N[u][v+1])
            N[u][v] = lower_bound + upper_bound - val
        else:
            if len(N) <= u:
                N.append([])
            N[u].append(lower_bound + x)

    return N

def sulzgruber_correspondence(M, copy=True):
    r"""
    Return the image of a `\lambda`-array ``M``
    under the Sulzgruber correspondence.

    The Sulzgruber correspondence is the map `\Phi_\lambda`
    from [Sulzgr2017]_ Section 7, and is the map
    `\xi_\lambda` from [Pak2002]_ Section 4.
    It is the inverse of the Pak correspondence
    (:meth:`pak_correspondence`).

    INPUT:

    - ``copy`` (default: ``True``) -- boolean;
      if set to ``False``, the algorithm will clobber the
      input (but be more efficient)

    EXAMPLES::

        sage: sulzgruber_correspondence([[1, 2, 3], [1, 2, 3], [2, 4, 4]])
        [[1, 0, 2], [0, 2, 0], [1, 1, 0]]
        sage: sulzgruber_correspondence([[1, 1, 4], [2, 3, 4], [4, 4, 4]])
        [[1, 1, 2], [0, 1, 0], [3, 0, 0]]
        sage: sulzgruber_correspondence([[0, 2, 3], [1, 3, 3], [2, 4]])
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: sulzgruber_correspondence([[1, 2, 4], [1, 3], [3]])
        [[0, 2, 2], [1, 1], [2]]
        sage: sulzgruber_correspondence([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]

    TESTS::

        sage: sulzgruber_correspondence([])
        []
        sage: sulzgruber_correspondence(((1, 2, 4), (1, 3), (3,)))
        [[0, 2, 2], [1, 1], [2]]

        sage: a = [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: sulzgruber_correspondence(a)
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: a
        [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: sulzgruber_correspondence(a, copy=False)
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: a
        []
    """
    lam = [len(row) for row in M]
    l = len(lam)
    if l == 0:
        return []
    # Finding a corner of lam...
    lam_0 = lam[0]
    for i, lam_i in enumerate(lam):
        if lam_i != lam_0:
            i -= 1
            break
    j = lam_0 - 1
    # Now, i is the index of the last row of ``M`` that
    # has the same length as the first row; hence, (i, j)
    # is a corner of lam.
    x = M[i][j]

    if copy:
        N = [list(row) for row in M] # make a deep copy of M to avoid vandalizing M
    else:
        N = M

    # toggling and obtaining x (the corner entry of the output):
    for k in range(min(i, j) + 1):
        u = i - k
        v = j - k
        if u > 0 and v > 0:
            lower_bound = max(N[u-1][v], N[u][v-1])
        elif u > 0:
            lower_bound = N[u-1][v]
        elif v > 0:
            lower_bound = N[u][v-1]
        else:
            lower_bound = 0
        if k > 0:
            val = N[u][v]
            upper_bound = min(N[u+1][v], N[u][v+1])
            N[u][v] = lower_bound + upper_bound - val
        else:
            x -= lower_bound

    # remove the corner (i, j):
    N[i].pop()
    if not N[i]:
        N.pop()

    N = sulzgruber_correspondence(N, copy=False)
    if len(N) <= i:
        N.append([])
    N[i].append(x)
    return N


