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

import sage.combinat.generator as generator
from sage.rings.arith import binomial
from sage.rings.integer import Integer
from sage.misc.mrange import xmrange
import sage.combinat.permutation
import itertools
import __builtin__
from combinat import CombinatorialClass, CombinatorialObject


def Words(alphabet, k):
    """
    Returns the combinatorial class of words of length k
    from alphabet.

    EXAMPLES:
        sage: w = Words([1,2,3], 2); w
        Words from [1, 2, 3] of length 2
        sage: w.count()
        9
        sage: w.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    """
    if isinstance(alphabet, (int, Integer)):
        alphabet = range(1,alphabet+1)

    return Words_alphabetk(alphabet, k)

class Words_alphabetk(CombinatorialClass):
    def __init__(self, alphabet, k):
        """
        TESTS:
            sage: import sage.combinat.word as word
            sage: w = word.Words_alphabetk([1,2,3], 2); w
            Words from [1, 2, 3] of length 2
            sage: w.count()
            9
        """
        self.alphabet = alphabet
        self.k = k

    def __repr__(self):
        """
        TESTS:
            sage: repr(Words([1,2,3], 2))
            'Words from [1, 2, 3] of length 2'
        """
        return "Words from %s of length %s"%(self.alphabet, self.k)


    def count(self):
        r"""
        Returns the number of words of length n from
        alphabet.

        EXAMPLES:
            sage: Words(['a','b','c'], 4).count()
            81
            sage: Words(3, 4).count()
            81
            sage: Words(0,0).count()
            1
            sage: Words(5,0).count()
            1
            sage: Words(['a','b','c'],0).count()
            1
            sage: Words(0,1).count()
            0
            sage: Words(5,1).count()
            5
            sage: Words(['a','b','c'],1).count()
            3
            sage: Words(7,13).count()
            96889010407L
            sage: Words(['a','b','c','d','e','f','g'],13).count()
            96889010407L

        """
        n = len(self.alphabet)
        return n**self.k

    def iterator(self):
        """
        Returns an iterator for all of the words
        of length k from alphabet. The iterator
        outputs the words in lexicographic order.

        TESTS:
            sage: [w for w in Words(['a', 'b'], 2)]
            [['a', 'a'], ['a', 'b'], ['b', 'a'], ['b', 'b']]
            sage: [w for w in Words(['a', 'b'], 0)]
            [[]]
            sage: [w for w in Words([], 3)]
            []
        """

        n = len(self.alphabet)

        if self.k == 0:
            yield []
            return

        for w in xmrange([n]*self.k):
            yield map(lambda x: self.alphabet[x], w)


def alphabet_order(alphabet):
    """
    Returns the ordering function of the alphabet A.
    The function takes two elements x and y of A and
    returns True if x occurs before y in A.

    EXAMPLES:
        sage: f = word.alphabet_order(['a','b','c'])
        sage: f('a', 'b')
        True
        sage: f('c', 'a')
        False
        sage: f('b', 'b')
        False

    """

    return lambda x,y: alphabet.index(x) < alphabet.index(y)



def standard(w, alphabet = None, ordering = None):
    """
    Returns the standard permutation of the word
    w on the ordered alphabet.  It is defined as
    the permutation with exactly the same number of
    inversions as w.

    By default, the letters are ordered using <.  A
    custom ordering can be specified by passing an
    ordering function into ordering.

    EXAMPLES:
        sage: word.standard([1,2,3,2,2,1,2,1])
        [1, 4, 8, 5, 6, 2, 7, 3]
    """

    if alphabet == None:
        alphabet = range(1,max(w)+1)

    if ordering == None:
        ordering = lambda x,y: x < y

    def ordering_cmp(x,y):
        """
        Returns -1 if x < y under ordering otherwise it returns 1.
        """
        if ordering(x,y):
            return -1
        else:
            return 1


    d = evaluation_dict(w)
    keys = d.keys()
    keys.sort(cmp=ordering_cmp)

    offset = 0
    temp = 0
    for k in keys:
        temp = d[k]
        d[k] = offset
        offset += temp

    result = []
    for l in w:
        d[l] += 1
        result.append(d[l])

    return sage.combinat.permutation.Permutation(result)


def evaluation_dict(w):
    """
    Returns a dictionary that represents the evalution
    of the word w.

    EXAMPLES:
        sage: word.evaluation_dict([2,1,4,2,3,4,2])
        {1: 1, 2: 3, 3: 1, 4: 2}
        sage: d = word.evaluation_dict(['b','a','d','b','c','d','b'])
        sage: map(lambda i: d[i], ['a','b','c','d'])
        [1, 3, 1, 2]
        sage: word.evaluation_dict([])
        {}

    """

    d = {}
    for l in w:
        if l not in d:
            d[l] = 1
        else:
            d[l] += 1

    return d

def evaluation_sparse(w):
    """

    EXAMPLES:
        sage: word.evaluation_sparse([4,4,2,5,2,1,4,1])
        [[1, 2], [2, 2], [4, 3], [5, 1]]
    """
    ed = evaluation_dict(w)
    keys = ed.keys()
    keys.sort()
    return [[x, ed[x]] for x in keys]

def evaluation_partition(w):
    """
    Returns the evaluation of the word w as a partition.

    EXAMPLE:
        sage: word.evaluation_partition([2,1,4,2,3,4,2])
        [3, 2, 1, 1]
    """

    d = evaluation_dict(w)
    p = d.values()
    p.sort(reverse=True)
    return p

def evaluation(w, alphabet=None):
    r"""
    Returns the evalutation of the word as a list.
    Either the word must be a word of integers or you
    must provide the alphabet as a list.

    EXAMPLES:
        sage: word.evaluation(['b','a','d','b','c','d','b'],['a','b','c','d','e'])
        [1, 3, 1, 2, 0]
        sage: word.evaluation([1,2,2,1,3])
        [2, 2, 1]
    """

    if alphabet == None:
        if len(w) == 0:
            alpha = 0
        else:
            alpha = max(w)

        e = [0] * alpha

        for letter in w:
            e[letter-1] += 1
        return e

    else:
        d = evaluation_dict(w)
        e = [0] * len(alphabet)
        indices = {}
        for i in range(len(alphabet)):
            indices[alphabet[i]] = i
        for l in w:
            if l in d:
                e[indices[l]] += 1
        return e

def from_standard_and_evaluation(sp, e, alphabet=None):
    """
    Returns the word given by the standard permutation
    sp and the evaluation e.

    EXAMPLES:
        sage: a = [1,2,3,2,2,1,2,1]
        sage: e = word.evaluation(a)
        sage: sp = word.standard(a)
        sage: b = word.from_standard_and_evaluation(sp, e)
        sage: b
        [1, 2, 3, 2, 2, 1, 2, 1]
        sage: a == b
        True

    """

    if sum(e) != len(sp):
        return

    if alphabet == None:
        alphabet = range(1, len(e)+1)
    elif isinstance(alphabet,Integer):
        alphabet = range(1, alphabet+1)

    part_sum = 0
    subs = []
    word = sp[:]
    for i in range(len(e)):
        next_part_sum = part_sum + e[i]
        for k in range(part_sum, next_part_sum):
            subs.append([k+1,alphabet[i]])
        part_sum = next_part_sum

    for sub in subs:
        word[word.index(sub[0])] = sub[1]

    return word

def swap(w,i,j=None):
   """
   Returns the word w with entries at positions i and
   j swapped.  By default, j = i+1.

   EXAMPLES:
       sage: word.swap([1,2,3],0,2)
       [3, 2, 1]
       sage: word.swap([1,2,3],1)
       [1, 3, 2]
   """
   if j == None:
       j = i+1
   new = w[:]
   (new[i], new[j]) = (new[j], new[i])
   return new

def swap_increase(w, i):
    """
    Returns the word w with positions i and i+1 exchanged
    if w[i] > w[i+1].  Otherwise, it returns w.

    EXAMPLES:
        sage: word.swap_increase([1,3,2],0)
        [1, 3, 2]
        sage: word.swap_increase([1,3,2],1)
        [1, 2, 3]

    """

    if w[i] > w[i+1]:
        return swap(w,i)
    else:
        return w

def swap_decrease(w, i):
    """
    Returns the word w with positions i and i+1 exchanged
    if w[i] < w[i+1].  Otherwise, it returns w.

    EXAMPLES:
        sage: word.swap_decrease([1,3,2],0)
        [3, 1, 2]
        sage: word.swap_decrease([1,3,2],1)
        [1, 3, 2]

    """

    if w[i] < w[i+1]:
        return swap(w,i)
    else:
        return w

def lex_less(w1,w2):
    """
    Returns true if the word w1 is lexicographically
    less than w2.  It is restricted to words whose
    letters can be compared by <.

    EXAMPLES:
        sage: word.lex_less([1,2,3],[1,3,2])
        True
        sage: word.lex_less([3,2,1],[1,2,3])
        False
    """

    i = 0
    while (i <= len(w1) ) and (i <= len(w2)) and (w1[i] == w2[i]):
        i += 1

    if i > len(w2):
        return False
    if i > len(w1):
        return True

    return w1[i] < w2[i]

def lex_cmp(w1,w2):
    """
    Returns -1 if the word w1 is lexicographically
    less than w2; otherwise, it returns 1.
    It is restricted to words whose
    letters can be compared by <.

    Useful to pass into Python's sort()

    EXAMPLES:
        sage: word.lex_cmp([1,2,3],[1,3,2])
        -1
        sage: word.lex_cmp([3,2,1],[1,2,3])
        1
    """
    if lex_less(w1, w2):
        return -1
    else:
        return 1

def deg_lex_less(w1,w2):
    """
    Returns true if the word w1 is degree
    lexicographically less than w2.  It is
    restricted to words whose letters can be
    compared by < as well as summed.

    EXAMPLES:
        sage: word.deg_lex_less([1,2,3],[1,3,2])
        True
        sage: word.deg_lex_less([1,2,4],[1,3,2])
        False
        sage: word.deg_lex_less([3,2,1],[1,2,3])
        False
    """

    d1 = sum(w1)
    d2 = sum(w2)

    if d1 != d2:
        return d1 < d2

    return lex_less(w1,w2)

def inv_lex_less(w1, w2):
    """
    Returns true if the word w1 is inverse
    lexicographically less than w2.  It is
    restricted to words whose letters can be
    compared by <.

    EXAMPLES:
        sage: word.inv_lex_less([1,2,4],[1,3,2])
        False
        sage: word.inv_lex_less([3,2,1],[1,2,3])
        True
    """
    if len(w1) != len(w2):
        return len(w1) < len(w2)

    for i in reversed(range(0,len(w1))):
        if w1[i] != w2[i]:
            return w1[i] < w2[i]

    return False



def deg_inv_lex_less(w1,w2):
    """
    Returns true if the word w1 is degree inverse
    lexicographically less than w2.  It is
    restricted to words whose letters can be
    compared by < as well as summed.

    EXAMPLES:
        sage: word.deg_inv_lex_less([1,2,4],[1,3,2])
        False
        sage: word.deg_inv_lex_less([3,2,1],[1,2,3])
        True
    """

    d1 = sum(w1)
    d2 = sum(w2)

    if d1 != d2:
        return d1 < d2

    return inv_lex_less(w1,w2)


def rev_lex_less(w1,w2):
    """
    Returns true if the word w1 is reverse
    lexicographically less than w2.  It is
    restricted to words whose letters can be
    compared by <.

    EXAMPLES:
        sage: word.rev_lex_less([1,2,4],[1,3,2])
        True
        sage: word.rev_lex_less([3,2,1],[1,2,3])
        False
    """

    if len(w1) != len(w2):
        return len(w1) > len(w2)

    for i in reversed(range(0,len(w1))):
        if w1[i] != w2[i]:
            return w1[i] > w2[i]

    return False


def deg_rev_lex_less(w1, w2):
    """
    Returns true if the word w1 is degree reverse
    lexicographically less than w2.  It is
    restricted to words whose letters can be
    compared by < as well as summed.

    EXAMPLES:
        sage: word.deg_rev_lex_less([1,2,4],[1,3,2])
        False
        sage: word.deg_rev_lex_less([3,2,1],[1,2,3])
        False
    """

    d1 = sum(w1)
    d2 = sum(w2)

    if d1 != d2:
        return d1 < d2

    return rev_lex_less(w1,w2)


def min_lex(l):
    """Returns the minimal element in lexicographic
    order of a l of words.

    EXAMPLES:
        sage: word.min_lex([[1,2,3],[1,3,2],[3,2,1]])
        [1, 2, 3]
    """

    if len(l) == 0:
        return None
    new_l = l[:]
    new_l.sort(cmp=lex_cmp)
    return new_l[0]


###################
# Shuffle Product #
###################


## def shuffle(w1,w2):
##     """
##     Returns the list of words in the shuffle product
##     of w1 and w2.

##     """
##     return shuffle_list(w1,w2)

## def shuffle_shifted(w1, w2):
##     """
##     Returns the shifted shuffle of w1 and w2.

##     EXAMPLES:

##     """

##     return shuffle(w1, [x+len(w1) for x in v])

## def shuffle_count(w1,w2):
##     """
##     Returns the number of words in the shuffle product
##     of w1 and w2.

##     It is given by binomial(len(w1)+len(w2), len(w1)).

##     EXAMPLES:
##         sage: word.shuffle_count([2,3],[3,4])
##         6
##     """

##     return binomial(len(w1)+len(w2), len(w1))


## def shuffle_list(w1,w2):
##     """
##     Returns the list of words in the shuffle product
##     of w1 and w2.

##     """

##     return [w for w in shuffle_iterator(w1,w2)]

## def shuffle_iterator(w1,w2):
##     """
##     Returns an iterator for the words in the
##     shuffle product of w1 and w2.

##     EXAMPLES:

##     """

##     n1 = len(w1)
##     n2 = len(w2)

##     def proc(vect):
##         i1 = -1
##         i2 = -1
##         res = []
##         for v in vect:
##             if v == 1:
##                 i1 += 1
##                 res.append(w1[i1])
##             else:
##                 i2 += 1
##                 res.append(w2[i2])

##     return itertools.imap(proc, integer_vector.iterator(n1, n1+n2, max_parts=1))

