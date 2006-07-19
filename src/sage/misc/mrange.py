"""
Multidimensional enumeration

AUTHORS:
    -- William Stein (2006-07-19)
    -- Jon Hanke
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import misc

def _xmrange(sizes, typ=list):
    n = len(sizes)
    if n == 0:
        return
    for i in sizes:
        if i <= 0:
            return
    v = [0] * n    # make a list of n 0's.
    v[-1] = -1
    ptr_max = n - 1
    ptr = ptr_max
    while True:
        while True:
            if ptr != -1 and v[ptr] + 1 < sizes[ptr]:
                v[ptr] += 1
                ptr = ptr_max
                break
            elif ptr != -1:
                v[ptr] = 0
                ptr -= 1
            else:
                return
        yield typ(v)   # make a copy of v!

def mrange(sizes, typ=list):
    """
    Return the multirange list with given sizes and type.

    This is the list version of xmrange.  Use xmrange for the
    iterator.

    More precisely, return the iterator over all objects of type typ
    of n-tuples of Python ints with entries between 0 and the integers
    in the sizes list.  The iterator is empty if sizes is empty or
    contains any non-positive integer.

    INPUT:
        sizes -- a list of nonnegative integers
        typ -- (default: list) a type or class; more generally,
               something that can be called with a list as input.

    OUTPUT:
        a list

    EXAMPLES:
        sage: mrange([3,2])
        [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
        sage: mrange([3,2], tuple)
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
        sage: mrange([3,2], sum)
        [0, 1, 1, 2, 2, 3]

    Examples that illustrate empty multi-ranges.
        sage: mrange([])
        []
        sage: mrange([5,3,-2])
        []
        sage: mrange([5,3,0])
        []

    AUTHORS:
        -- Jon Hanke
        -- William Stein
    """
    return list(_xmrange(sizes, typ))




class xmrange:
    """
    Return the multirange iterate with given sizes and type.

    More precisely, return the iterator over all objects of type typ
    of n-tuples of Python ints with entries between 0 and the integers
    in the sizes list.  The iterator is empty if sizes is empty or
    contains any non-positive integer.

    Use mrange for the non-iterator form.

    INPUT:
        sizes -- a list of nonnegative integers
        typ -- (default: list) a type or class; more generally,
               something that can be called with a list as input.

    OUTPUT:
        a generator

    EXAMPLES:
    We create multi-range iterators, print them and also iterate
    through a tuple version.
        sage: z = xmrange([3,2]);z
        xmrange([3, 2])
        sage: z = xmrange([3,2], tuple);z
        xmrange([3, 2], <type 'tuple'>)
        sage: for a in z:
        ...    print a
        (0, 0)
        (0, 1)
        (1, 0)
        (1, 1)
        (2, 0)
        (2, 1)

    We illustrate a few more iterations.
        sage: list(xmrange([3,2]))
        [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
        sage: list(xmrange([3,2], tuple))
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    Here we compute the sum of each element of the multi-range iterator:
        sage: list(xmrange([3,2], sum))
        [0, 1, 1, 2, 2, 3]

    Next we compute the product:
        sage: list(xmrange([3,2], prod))
        [0, 0, 0, 1, 0, 2]

    Examples that illustrate empty multi-ranges.
        sage: list(xmrange([]))
        []
        sage: list(xmrange([5,3,-2]))
        []
        sage: list(xmrange([5,3,0]))
        []

    We use a multi-range iterator to iterate through the cartesian
    product of sets.
        sage: X = ['red', 'apple', 389]
        sage: Y = ['orange', 'horse']
        sage: for i,j in xmrange([len(X), len(Y)]):
        ...    print (X[i], Y[j])
        ('red', 'orange')
        ('red', 'horse')
        ('apple', 'orange')
        ('apple', 'horse')
        (389, 'orange')
        (389, 'horse')

    AUTHORS:
        -- Jon Hanke
        -- William Stein
    """
    def __init__(self, sizes, typ=list):
        self.sizes = [int(x) for x in sizes]
        self.typ = typ

    def __repr__(self):
        if self.typ == list:
            return 'xmrange(%s)'%self.sizes
        else:
            return 'xmrange(%s, %s)'%(self.sizes, self.typ)

    def __len__(self):
        sizes = self.sizes
        n = len(sizes)
        if n == 0:
            return 0
        for i in sizes:
            if i <= 0:
                return 0
        return misc.prod(sizes, 1)

    def __iter__(self):
        return _xmrange(self.sizes, self.typ)

