r"""
Gray codes.

REFERENCES:

.. [Knuth-TAOCP2A] D. Knuth "The art of computer programming", fascicules 2A,
   "generating all n-tuples"

.. [Knuth-TAOCP3A] D. Knuth "The art of computer programming", fascicule 3A
   "generating all combinations"
"""

def integer_vector_gray_code_switch(m):
    r"""
    Return the switch to perform for the enumeration of tuples of non negative
    integers up to a specified upper bound.

    The iterator return at each step a pair `(position, increment)` which
    corresponds to a sequence of switches to perform in position ``position`` by
    the increment ``increment``. By construction, the increment is either ``+1`` or
    ``-1``.

    This is algorithm H in [Knuth-TAOCP2A]_: loopless reflected mixed-radix Gray
    generation.

    INPUT:

    - ``m`` -- a list or tuple of upper bound limits.

    EXAMPLES::

        sage: from sage.combinat.gray_codes import integer_vector_gray_code_switch
        sage: l = [0,0,0]
        sage: for i,j in integer_vector_gray_code_switch([3,3,3]):
        ....:     l[i] += j
        ....:     print l
        [1, 0, 0]
        [2, 0, 0]
        [2, 1, 0]
        [1, 1, 0]
        [0, 1, 0]
        [0, 2, 0]
        [1, 2, 0]
        [2, 2, 0]
        [2, 2, 1]
        [1, 2, 1]
        [0, 2, 1]
        [0, 1, 1]
        [1, 1, 1]
        [2, 1, 1]
        [2, 0, 1]
        [1, 0, 1]
        [0, 0, 1]
        [0, 0, 2]
        [1, 0, 2]
        [2, 0, 2]
        [2, 1, 2]
        [1, 1, 2]
        [0, 1, 2]
        [0, 2, 2]
        [1, 2, 2]
        [2, 2, 2]
        sage: l = [0,0]
        sage: for i,j in integer_vector_gray_code_switch([2,1]):
        ....:     l[i] += j
        ....:     print l
        [1, 0]

    TESTS::

        sage: for t in [[2,2,2],[2,1,2],[3,2,1],[2,1,3]]:
        ....:     assert sum(1 for _ in integer_vector_gray_code_switch(t)) == prod(t)-1
    """
    # n is the length of the vector (where we ignore upper bounds equal to 1)
    n = k = 0

    new_m = []
    mm = []
    for i in m:
        if i <= 0:
            raise ValueError("accept only positive integers")
        if i == 1:
            k += 1
        else:
            new_m.append(i-1)
            mm.append(k)
            n += 1
    m = new_m
    f = range(n+1)  # focus pointer
    o = [1] * n     # switch +1 or -1
    a = [0] * n     # current vector

    j = f[0]
    while j != n:
        f[0] = 0            
        oo = o[j]            
        a[j] += oo
        if a[j] == 0 or a[j] == m[j]:
            f[j] = f[j+1]
            f[j+1] = j+1
            o[j] = -oo

        yield (j+mm[j], oo)

        j = f[0]

def combination_gray_code_switch(n,t):
    r"""
    Iterator through the switches of the revolving door algorithm.

    The revolving door algorithm is a way to generate all combinations of a set
    (i.e. the subset of given cardinality) in such way that two consecutive
    subsets differ by one element. At each step, the iterator output a pair
    `(i,j)` where the item `i` has to be removed and `j` has to be added.

    The ground set is always `\{0, 1, ..., n-1\}`.

    See [Knuth-TAOCP3A]_.

    INPUT:

    - ``n`` -- size of the set

    - ``t`` -- size of subsets

    EXAMPLES::

        sage: from sage.combinat.gray_codes import combination_gray_code_switch
        sage: b = [1, 1, 1, 0, 0]
        sage: for i,j in combination_gray_code_switch(5,3):
        ....:     b[i] = 0; b[j] = 1
        ....:     print b
        [1, 0, 1, 1, 0]
        [0, 1, 1, 1, 0]
        [1, 1, 0, 1, 0]
        [1, 0, 0, 1, 1]
        [0, 1, 0, 1, 1]
        [0, 0, 1, 1, 1]
        [1, 0, 1, 0, 1]
        [0, 1, 1, 0, 1]
        [1, 1, 0, 0, 1]

        sage: b = [1,0,0,0]
        sage: for i,j in combination_gray_code_switch(4,1):
        ....:     b[i] = 0; b[j] = 1
        ....:     print b
        [0, 1, 0, 0]
        [0, 0, 1, 0]
        [0, 0, 0, 1]
    """
    t = int(t)
    n = int(n)
    assert 0 <= t and t <= n, "Parameters not admissible"
    if t == 0 or t == n:
        return iter([])
    if t % 2:
        return _revolving_door_switch_iterator_odd(n,t)
    else:
        return _revolving_door_switch_iterator_even(n,t)

def _revolving_door_switch_iterator_odd(n,t):
    r"""
    Revolving door switch for odd `t`.

    TESTS::

        sage: from sage.combinat.gray_codes import _revolving_door_switch_iterator_odd
        sage: sum(1 for _ in _revolving_door_switch_iterator_odd(13,3)) == binomial(13,3) - 1
        True
        sage: sum(1 for _ in _revolving_door_switch_iterator_odd(10,5)) == binomial(10,5) - 1
        True
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R3 : easy case
        if c[0] + 1 < c[1]:
            yield c[0], c[0]+1
            c[0] += 1
            continue

        j = 1
        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                yield c[j-1], c[j]+1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break

def _revolving_door_switch_iterator_even(n,t):
    r"""
    Revolving door algorithm for even `t`.

    TESTS::

        sage: from sage.combinat.gray_codes import _revolving_door_switch_iterator_even
        sage: sum(1 for _ in _revolving_door_switch_iterator_even(13,4)) == binomial(13,4) - 1
        True
        sage: sum(1 for _ in _revolving_door_switch_iterator_even(12,6)) == binomial(12,6) - 1
        True
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R3 : easy case
        if c[0] > 0:
            yield c[0], c[0]-1
            c[0] -= 1
            continue

        j = 1
        # R5 : try to increase c[j]
        # at this point c[j-1] = j-1
        if c[j] + 1 < c[j+1]:
            yield c[j-1], c[j]+1
            c[j-1] = c[j]
            c[j] += 1
            continue
        j += 1

        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                yield c[j-1], c[j] + 1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break

