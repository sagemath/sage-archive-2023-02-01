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
        sage: l = [1,2,3]
        sage: r,u = sage.combinat.ranker.from_list(l)
        sage: r(1)
        0
        sage: r(3)
        2
        sage: u(2)
        3
        sage: u(0)
        1
    """

    n = len(l)
    def unrank(j):
        if j < 0 or j >= n:
            raise ValueError, "Argument j ( = %s ) must be between 0 and %d"%(str(j), n)

        return l[j]

    def rank(obj):
        return l.index(obj)

    return [rank, unrank]


def rank_from_list(l):
    """
    Returns a rank function given a list l.

    EXAMPLES:
        sage: l = [1,2,3]
        sage: r = sage.combinat.ranker.rank_from_list(l)
        sage: r(1)
        0
        sage: r(3)
        2
    """
    def rank(obj):
        return l.index(obj)

    return rank


def unrank_from_list(l):
    """
    Returns an unrank function from a list.

    EXAMPLES:
        sage: l = [1,2,3]
        sage: u = sage.combinat.ranker.unrank_from_list(l)
        sage: u(2)
        3
        sage: u(0)
        1
    """
    n = len(l)
    def unrank(j):
        if j < 0 or j >= n:
            raise ValueError, "Argument j ( = %s ) must be between 0 and %d"%(str(j), n)

        return l[j]

    return unrank
