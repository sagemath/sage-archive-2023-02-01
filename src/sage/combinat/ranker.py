#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                          Nicolas M. Thiéry <nthiery at users.sf.net>
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

from sage.misc.all import cached_function

def from_list(l):
    """
    Returns a ranker from the list l.

    INPUT:
        l -- a list

    OUTPUT
        [rank, unrank] -- functions

    EXAMPLES:
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

    EXAMPLES:
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

    EXAMPLES:
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
    passed as argument. The object should be hashable.  unrank is the
    inverse function; it returns None for indices that have not yet
    been assigned.

    EXAMPLES:
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

    TODO: add tests as in combinat::rankers
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
