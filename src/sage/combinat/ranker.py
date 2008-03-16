#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
