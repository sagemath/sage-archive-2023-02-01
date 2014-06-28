"""
Rankers
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#  Ported from MuPAD-Combinat (combinat::rankers)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from collections import Iterable, Sequence
from sage.misc.cachefunc import cached_function
from sage.structure.parent import Parent
from sage.categories.enumerated_sets import EnumeratedSets

def from_list(l):
    """
    Returns a ranker from the list l.

    INPUT:

    -  ``l`` - a list

    OUTPUT:

    - ``[rank, unrank]`` - functions

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = [1,2,3]
        sage: r,u = ranker.from_list(l)
        sage: r(1)
        0
        sage: r(3)
        2
        sage: u(2)
        3
        sage: u(0)
        1
    """
    return [rank_from_list(l), unrank_from_list(l)]


def rank_from_list(l):
    """
    Returns a rank function given a list l.

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = [1,2,3]
        sage: r = ranker.rank_from_list(l)
        sage: r(1)
        0
        sage: r(3)
        2
    """
    rank = lambda obj: l.index(obj)

    return rank


def unrank_from_list(l):
    """
    Returns an unrank function from a list.

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = [1,2,3]
        sage: u = ranker.unrank_from_list(l)
        sage: u(2)
        3
        sage: u(0)
        1
    """
    unrank = lambda j: l[j]
    return unrank

def on_fly():
    """
    Returns a pair of enumeration functions rank / unrank.

    rank assigns on the fly an integer, starting from 0, to any object
    passed as argument. The object should be hashable. unrank is the
    inverse function; it returns None for indices that have not yet
    been assigned.

    EXAMPLES::

        sage: [rank, unrank] = sage.combinat.ranker.on_fly()
        sage: rank('a')
        0
        sage: rank('b')
        1
        sage: rank('c')
        2
        sage: rank('a')
        0
        sage: unrank(2)
        'c'
        sage: unrank(3)
        sage: rank('d')
        3
        sage: unrank(3)
        'd'

    .. todo:: add tests as in combinat::rankers
    """
    def count():
        i = 0
        while True:
            yield i
            i+=1

    counter = count()

    @cached_function
    def rank(x):
        i = counter.next()
        unrank.set_cache(x, i)
        return i

    @cached_function
    def unrank(i):
        return None

    return [rank, unrank]

def unrank(L, i):
    r"""
    Return the `i`-th element of `L`.

    INPUT:

    - ``L`` -- a list, tuple, finite enumerated set, ...
    - ``i`` -- an int or :class:`Integer`

    The purpose of this utility is to give a uniform idiom to recover
    the `i`-th element of an object ``L``, whether ``L`` is a list,
    tuple (or more generally a :class:`collections.Sequence`), an
    enumerated set, some old parent of Sage still implementing
    unranking in the method ``__getitem__``, or an iterable (see
    :class:`collections.Iterable`). See :trac:`15919`.

    EXAMPLES:

    Lists, tuples, and other :class:`sequences <collections.Sequence>`::

        sage: from sage.combinat.ranker import unrank
        sage: unrank(['a','b','c'], 2)
        'c'
        sage: unrank(('a','b','c'), 1)
        'b'
        sage: unrank(xrange(3,13,2), 1)
        5

    Enumerated sets::

        sage: unrank(GF(7), 2)
        2
        sage: unrank(IntegerModRing(29), 10)
        10

    An old parent with unranking implemented in ``__getitem__``::

        sage: M = MatrixSpace(GF(3), 2, 2)
        sage: hasattr(M, "unrank")
        False
        sage: M[42]
        [1 0]
        [2 1]
        sage: unrank(M, 42)
        [1 0]
        [2 1]

    An iterable::

        sage: unrank(NN,4)
        4

    An iterator::

        sage: unrank(('a{}'.format(i) for i in range(20)), 0)
        'a0'
        sage: unrank(('a{}'.format(i) for i in range(20)), 2)
        'a2'

    .. WARNING::

        When unranking an iterator, it returns the ``i``-th element
        beyond where it is currently at::

            sage: from sage.combinat.ranker import unrank
            sage: it = iter(range(20))
            sage: unrank(it, 2)
            2
            sage: unrank(it, 2)
            5

    TESTS::

        sage: from sage.combinat.ranker import unrank
        sage: unrank(range(3), 10)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range

        sage: unrank(('a{}'.format(i) for i in range(20)), 22)
        Traceback (most recent call last):
        ...
        IndexError: index out of range

        sage: M[100]
        Traceback (most recent call last):
        ...
        IndexError: list index out of range
    """
    if L in EnumeratedSets:
        return L.unrank(i)
    if isinstance(L, Sequence):
        return L[i]
    if isinstance(L, Parent):
        # handle parents still implementing unranking in __getitem__
        try:
            return L[i]
        except (AttributeError, TypeError, ValueError):
            pass
    if isinstance(L, Iterable):
        try:
            it = iter(L)
            for _ in range(i):
                it.next()
            return it.next()
        except StopIteration as e:
            raise IndexError("index out of range")
    raise ValueError("Don't know how to unrank on {}".format(L))

