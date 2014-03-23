r"""
Iterator over a conjugacy classe of the symmetric group
"""

from partition import Partition
from set_partition import SetPartitions
import itertools

def conjugacy_class_iterator(part, S=None):
    r"""
    Return an iterator over the conjugacy class associated to the partition
    ``part``.
    
    The elements are given as a list of tuples, each tuple being a cycle.

    INPUT:

    - ``part`` -- partition for the size of the subsets

    - ``S`` -- an optional set. If not specified, `{1,...,n}` is used.

    EXAMPLES::

        sage: from sage.combinat.permutation_conjugacy_class_iterator import conjugacy_class_iterator
        sage: for p in conjugacy_class_iterator([2,2]): print p
        [(1, 2), (3, 4)]
        [(1, 3), (2, 4)]
        [(1, 4), (2, 3)]

    In order to get permutations, one can use ``imap`` from the Python module
    ``itertools``::

        sage: from itertools import imap
        sage: S = SymmetricGroup(5)
        sage: for p in imap(S,conjugacy_class_iterator([3,2])): print p
        (1,2)(3,4,5)
        (1,2)(3,5,4)
        (1,3)(2,4,5)
        (1,3)(2,5,4)
        ...
        (1,4,2)(3,5)
        (1,2,3)(4,5)
        (1,3,2)(4,5)

    Check that the number of elements correspond to the number of elements in
    the conjugacy class

        sage: for p in Partitions(5):
        ....:     print p
        ....:     print sum(1 for _ in conjugacy_class_iterator(p))
        ....:     print p.conjugacy_class_size()
        [5]
        24
        24
        [4, 1]
        30
        30
        [3, 2]
        20
        20
        [3, 1, 1]
        20
        20
        [2, 2, 1]
        15
        15
        [2, 1, 1, 1]
        10
        10
        [1, 1, 1, 1, 1]
        1
        1

    It is also possible to specify any underlying set::

        sage: it = conjugacy_class_iterator([2,2,2],'abcdef')
        sage: it.next()
        [('a', 'c'), ('b', 'e'), ('d', 'f')]
        sage: it.next()
        [('a', 'c'), ('b', 'd'), ('e', 'f')]
    """
    part = Partition(part)
    n = sum(part)
    if S is None:
        S = range(1,n+1)
    else:
        S = list(S)
    if n != len(S):
        raise ValueError("the sum of the partition %s does not match the size of %s"%(part,S))

    m = len(part)
    for s in SetPartitions(S, part):
        firsts = [t[0] for t in s]
        rests = [t[1:] for t in s]
        iterator = tuple(itertools.permutations(r) for r in rests)
        for r in itertools.product(*iterator):
            yield [(firsts[i],)+r[i] for i in xrange(m)]

