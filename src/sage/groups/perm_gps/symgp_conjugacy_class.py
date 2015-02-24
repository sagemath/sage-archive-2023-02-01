r"""
Conjugacy Classes Of The Symmetric Group

AUTHORS:

- Vincent Delacroix, Travis Scrimshaw (2014-11-23)
"""

from sage.groups.conjugacy_classes import ConjugacyClassGAP
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.partition import Partitions_n
from sage.combinat.set_partition import SetPartitions
from sage.sets.set import Set
import itertools

class SymmetricGroupConjugacyClass(ConjugacyClassGAP):
    """
    A conjugacy class of the symmetric group.

    INPUT:

    - ``group`` -- the symmetric group
    - ``part`` -- a partition or an element of ``group``
    """
    def __init__(self, group, part):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(g)
            sage: TestSuite(C).run()
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: TestSuite(C).run()
        """
        P = Partitions_n( len(group.domain()) )
        if isinstance(part, PermutationGroupElement) and part.parent() is group:
            elt = part
            part = P( sorted(map(len, part.cycle_tuples(True)), reverse=True) )
        else:
            part = P(part)
            elt = default_representative(part, group)
        self._part = part
        self._set = None
        ConjugacyClassGAP.__init__(self, group, elt)

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: G.conjugacy_class(Partition([4]))
            Conjugacy class of cycle type [4] in
             Symmetric group of order 4! as a permutation group
        """
        return "Conjugacy class of cycle type %s in %s"%(self._part, self._parent)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: C = G.conjugacy_class(Partition([3,1]))
            sage: for x in C: x
            (2,3,4)
            (2,4,3)
            (1,3,4)
            (1,4,3)
            (1,2,4)
            (1,4,2)
            (1,2,3)
            (1,3,2)
        """
        if self._set:
            for x in self._set:
                yield x
            return

        for x in conjugacy_class_iterator(self._part, self._parent.domain()):
            yield PermutationGroupElement(x, self._parent, check=False)

    def __cmp__(self, other):
        r"""
        Comparison of conjugacy classes is done by comparing the
        defining cycle types.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(Partition([3,2]))
            sage: Cp = G.conjugacy_class(g)
            sage: C == Cp
            True
        """
        c = cmp(type(self), type(other))
        if c:
             return c
        return cmp(self._part, other._part)

    def partition(self):
        """
        Return the partition of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: g = G([(1,2), (3,4,5)])
            sage: C = G.conjugacy_class(g)
        """
        return self._part

    def set(self):
        r"""
        The set of all elements in the conjugacy class ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: g = G((1,2))
            sage: C = G.conjugacy_class(g)
            sage: S = [(2,3), (1,2), (1,3)]
            sage: C.set() == Set(G(x) for x in S)
            True
        """
        if self._set is None:
            to_elt = lambda x: PermutationGroupElement(x, self._parent, check=False)
            self._set = Set(to_elt(x) for x
                            in conjugacy_class_iterator(self._part, self._parent.domain()) )
        return self._set

#####################################################################
## Helper functions

def default_representative(part, G):
    r"""
    Construct the default representative for the conjugacy class of
    cycle type ``part`` of a symmetric group ``G``.

    Let `\lambda` be a partition of `n`. We pick a representative by

    .. MATH::

        (1, 2, \ldots, \lambda_1)
        (\lambda_1 + 1, \ldots, \lambda_1 + \lambda_2)
        (\lambda_1 + \lambda_2 + \cdots + \lambda_{\ell-1}, \ldots, n),

    where `\ell` is the length (or number of parts) of `\lambda`.

    INPUT:

    - ``part`` -- partition for the size of the subsets

    - ``G`` -- a symmetric group

    EXAMPLES::

        sage: from sage.groups.perm_gps.symgp_conjugacy_class import default_representative
        sage: S = SymmetricGroup(4)
        sage: for p in Partitions(4):
        ....:     print default_representative(p, S)
        (1,2,3,4)
        (1,2,3)
        (1,2)(3,4)
        (1,2)
        ()
    """
    D = G.domain()
    total = 0
    cycles = []
    for p in part:
        cycles.append(tuple(D[total:total+p]))
        total += p
    # TODO: Change this to G.element_class(cycles, check=False)
    #   once SymmetricGroup is a proper parent.
    return PermutationGroupElement(cycles, G, check=False)

def conjugacy_class_iterator(part, S=None):
    r"""
    Return an iterator over the conjugacy class associated to
    the partition ``part``.
    
    The elements are given as a list of tuples, each tuple being a cycle.

    INPUT:

    - ``part`` -- partition for the size of the subsets

    - ``S`` -- (optional) a set, if not specified `{1, ..., n}` is used

    EXAMPLES::

        sage: from sage.groups.perm_gps.symgp_conjugacy_class import conjugacy_class_iterator
        sage: for p in conjugacy_class_iterator([2,2]): print p
        [(1, 2), (3, 4)]
        [(1, 3), (2, 4)]
        [(1, 4), (2, 3)]

    In order to get permutations, one can use ``imap`` from the Python
    module ``itertools``::

        sage: from itertools import imap
        sage: S = SymmetricGroup(5)
        sage: for p in imap(S, conjugacy_class_iterator([3,2])): print p
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

        sage: s = lambda p: sum(1 for _ in conjugacy_class_iterator(p))
        sage: all(s(p) == p.conjugacy_class_size() for p in Partitions(5))
        True

    It is also possible to specify any underlying set::

        sage: it = conjugacy_class_iterator([2,2,2], 'abcdef')
        sage: next(it)
        [('a', 'c'), ('b', 'e'), ('d', 'f')]
        sage: next(it)
        [('a', 'c'), ('b', 'd'), ('e', 'f')]
    """
    n = sum(part)
    part = Partitions_n(n)(part)
    if S is None:
        S = range(1, n+1)
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

