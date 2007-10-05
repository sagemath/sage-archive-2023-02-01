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

from sage.rings.arith import binomial
import random as rnd

def count(n,k):
    """
    Returns the number of multichoices of k things from a list
    of n things.

    EXAMPLES:
        sage: multichoose_nk.count(3,2)
        6
    """
    return binomial(n+k-1,k)

def iterator(n,k):
    """
    An iterator for all multichoies of k thinkgs from range(n).

    EXAMPLES:
        sage: [c for c in multichoose_nk.iterator(3,2)]
        [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]

    """
    dif = 0
    if k == 0:
        yield []
        return

    if n < 1+(k-1)*dif:
        return
    else:
        subword = [ i*dif for i in range(k) ]

    yield subword[:]
    finished = False

    while not finished:
        #Find the biggest element that can be increased
        if subword[-1] < n-1:
            subword[-1] += 1
            yield subword[:]
            continue

        finished = True
        for i in reversed(range(k-1)):
            if subword[i]+dif < subword[i+1]:
                subword[i] += 1
                #Reset the bigger elements
                for j in range(1,k-i):
                    subword[i+j] = subword[i]+j*dif
                yield subword[:]
                finished = False
                break

    return

def list(n,k,repitition=False):
    """
    Returns a list of all the multichoices of k things from range(n).

    EXAMPLES:
        sage: multichoose_nk.list(3,2)
        [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]

    """

    return [c for c in iterator(n,k)]

def random(n,k):
    """
    Returns a random multichoice of k things from range(n).

    EXAMPLES:
        sage: multichoose_nk.random(5,2) #random
        [0,3]
        sage: multichoose_nk.random(5,2) #random
        [2,2]
    """
    rng = range(n)
    r = []
    for i in range(k):
        r.append( rnd.choice(rng))

    r.sort()
    return r
