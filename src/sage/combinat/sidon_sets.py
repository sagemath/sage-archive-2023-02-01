r"""
Implementations for Sidon sets and their generalizations, Sidon `g`-sets.

AUTHORS:

- Martin Raum (07-25-2011)
"""

from sage.misc.all import cached_function
from sage.rings.all import Integer

@cached_function
def sidon_sets(N, g = 1) :
    r"""
    Compute all Sidon-`g` sets, that have elements less or equal than `N`.

    INPUT:

    - `N` -- A positive integer.
    - `g` -- A positive integer (default: `1`).

    OUTPUT:

    - A frozen set of frozen sets of integers.

    EXAMPLES::

        sage: sidon_sets(3, 2)
        frozenset([frozenset([3]), frozenset([1, 2]), frozenset([]), frozenset([2, 3]), frozenset([1]), frozenset([1, 3]), frozenset([1, 2, 3]), frozenset([2])])

    TESTS::

        sage: sidon_sets(1,1)
        frozenset([frozenset([]), frozenset([1])])
        sage: frozenset([1,2,4,8,13]) in sidon_sets(13)
        True
        sage: from itertools import groupby
        sage: all(all(l <= 3 for l in map(lambda s: len(list(s[1])), groupby(sorted(a + ap for a in sid for ap in sid if a >= ap), lambda s: s))) for sid in sidon_sets(10, 3))
        True
        sage: all(all(l <= 5 for l in map(lambda s: len(list(s[1])), groupby(sorted(a + ap for a in sid for ap in sid if a >= ap), lambda s: s))) for sid in sidon_sets(10, 5))
        True
        sage: sidon_sets(-1,3)
        Traceback (most recent call last):
        ...
        ValueError: N must be a positive integer
        sage: sidon_sets(1, -3)
        Traceback (most recent call last):
        ...
        ValueError: g must be a positive integer
    """
    if not isinstance(N, (int, Integer)) or N < 1 :
        raise ValueError( "N must be a positive integer" )
    elif not isinstance(g, (int, Integer)) or g < 1 :
        raise ValueError( "g must be a positive integer" )

    if N == 1 :
        return frozenset([frozenset([]), frozenset([1])])

    pre_sidons = sidon_sets(N - 1, g)
    sidons = set(pre_sidons)
    for psid in pre_sidons :
        psid_shift = set(n - 1 for n in psid if n != 1)
        psid_shift.add(N - 1)
        if not psid_shift in pre_sidons :
            continue

        if not 1 in psid :
            add_sid = True
        else :
            add_sid = True
            Np1_count = 0
            for n in psid :
                if N + 1 - n in psid and 2 * n <= N + 1:
                    Np1_count += 1
                    if Np1_count >= g :
                        add_sid = False
                        break

        if add_sid :
            sid = set(psid)
            sid.add(N)
            sidons.add(frozenset(sid))

    return frozenset(sidons)
