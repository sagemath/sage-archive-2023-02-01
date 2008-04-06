r"""
Words
"""
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
from sage.misc.mrange import xmrange
import sage.combinat.permutation
import itertools
from combinat import CombinatorialClass
from sage.rings.all import binomial, Integer, infinity
from sage.combinat.integer_vector import IntegerVectors
import copy
import partition
from subset import Subsets

def Words(*args):
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
    if len(args) == 0:
        return Words_all()
    elif len(args) == 1:
        if isinstance(args[0], (int, Integer)):
            return Words_n(args[0])
        elif isinstance(args[0], list):
            return Words_alphabet(args[0])
    elif len(args) == 2:
        if isinstance(args[0], (int, Integer)):
            alphabet = range(1,args[0]+1)
        else:
            alphabet = args[0]
        k = args[1]
        return Words_alphabetk(alphabet, k)

    raise ValueError, "do not know how to make a combinatorial class of words from %s"%args

class Words_all(CombinatorialClass):
    def __repr__(self):
        """
        EXAMPLES:
            sage: Words().__repr__()
            'Words'
        """
        return "Words"

    def count(self):
        """
        EXAMPLES:
            sage: Words().count()
            +Infinity
        """
        return infinity

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: 2 in Words()
            False
            sage: [1,2] in Words()
            True
        """
        return isinstance(x, list)

    def iterator(self):
        """
        EXAMPLES:
            sage: Words().list() #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class Words_n(CombinatorialClass):
    def __init__(self, n):
        """
        EXAMPLES:
            sage: w = Words(3)
            sage: w == loads(dumps(w))
            True
        """
        self.n = n

    def __repr__(self):
        """
        EXAMPLES:
            sage: Words(3).__repr__()
            'Words of length 3'
        """
        return "Words of length %s"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Words(3).count()
            +Infinity
        """
        return infinity

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: 2 in Words(3)
            False
            sage: [1,'a',3] in Words(3)
            True
            sage: [1,2] in Words(3)
            False
        """
        return isinstance(x, list) and len(x) == self.n

    def iterator(self):
        """
        EXAMPLES:
            sage: Words(3).list() #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class Words_alphabet(CombinatorialClass):
    def __init__(self, alphabet):
        """
        EXAMPLES:
            sage: w = Words([1,2,3])
            sage: w == loads(dumps(w))
            True
        """
        self.alphabet = alphabet

    def __repr__(self):
        """
        EXAMPLES:
            sage: Words([1,2,3]).__repr__()
            'Words from the alphabet [1, 2, 3]'
        """
        return "Words from the alphabet %s"%self.alphabet

    def count(self):
        """
        EXAMPLES:
            sage: Words([1,2,3]).count()
            +Infinity
        """
        return infinity

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: 2 in Words([1,2,3])
            False
            sage: [2] in Words([1,2,3])
            True
            sage: [1, 'a'] in Words([1,2,3])
            False
        """
        return isinstance(x, list) and all(map(lambda i: i in self.alphabet, x))

    def iterator(self):
        """
        EXAMPLES:
            sage: Words([1,2,3]).list() #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class Words_alphabetk(CombinatorialClass):
    def __init__(self, alphabet, k):
        """
        TESTS:
            sage: w = Words([1,2,3], 2); w
            Words from [1, 2, 3] of length 2
            sage: w.count()
            9
        """
        self.alphabet = alphabet
        self.k = k

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: w = Words([1,2,3], 2)
            sage: [1,2,3] in w
            False
            sage: [1,2] in w
            True
            sage: [3,4] in w
            False
        """
        return len(x) == self.k and all(i in self.alphabet for i in x)

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
            96889010407L               # 32-bit
            96889010407                # 64-bit
            sage: Words(['a','b','c','d','e','f','g'],13).count()
            96889010407L               # 32-bit
            96889010407                # 64-bit
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
        sage: word.standard([1,2,3,2,2,1,2,1])
        [1, 4, 8, 5, 6, 2, 7, 3]
    """

    if alphabet == None:
        alphabet = range(1,max(w)+1)

    if ordering == None:
        ordering = lambda x,y: x < y

    #Returns -1 if x < y under ordering otherwise it returns 1.
    ordering_cmp = lambda x,y: -1 if ordering(x,y) else 0

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


def charge(word, check=True):
    """

    EXAMPLES:
        sage: from sage.combinat.word import charge

        sage: charge([1,1,2,2,3])
        0
        sage: charge([3,1,1,2,2])
        1
        sage: charge([2,1,1,2,3])
        1
        sage: charge([2,1,1,3,2])
        2
        sage: charge([3,2,1,1,2])
        2
        sage: charge([2,2,1,1,3])
        3
        sage: charge([3,2,2,1,1])
        4

        sage: charge([3,3,2,1,1])
        Traceback (most recent call last):
        ...
        ValueError: the evaluation of w must be a partition
    """
    if check:
        if evaluation(word) not in partition.Partitions():
            raise ValueError, "the evaluation of w must be a partition"
    w = word[:]
    res = 0
    while len(w) != 0:
        i = 0
        l = 1
        index = 0
        while len(w) != 0 and l <= max(w):
            while w[i] != l:
                i += 1
                if i >= len(w):
                    i = 0
                    index += 1
            res += index
            l += 1
            w.pop(i)
            if i >= len(w):
                i = 0
                index += 1
    return res

def evaluation_dict(w):
    """
    Returns a dictionary that represents the evalution
    of the word w.

    EXAMPLES:
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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
        sage: import sage.combinat.word as word
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

def ShuffleProduct(w1, w2, shifted=False, overlapping=False):
    """
    Returns the combinatorial class representing the shuffle product between
    w1 and w2.  This consists of all words of length len(w1)+len(w2) which
    have both w1 and w2 as subwords.

    EXAMPLES:
        sage: sp = ShuffleProduct([1,2],[3,4]); sp
        Shuffle product of [1, 2] and [3, 4]
        sage: sp.count()
        6
        sage: sp.list()
        [[1, 2, 3, 4],
         [1, 3, 2, 4],
         [1, 3, 4, 2],
         [3, 1, 2, 4],
         [3, 1, 4, 2],
         [3, 4, 1, 2]]
    """
    if shifted is True:
        return ShuffleProduct_shifted(w1, w2)

    if overlapping is False:
        pass
    elif overlapping is True:
        return ShuffleProduct_overlapping(w1, w2)
    elif isinstance(overlapping, (int, Integer)) and overlapping:
        return ShuffleProduct_overlapping_r(w1, w2, overlapping)
    else:
        raise ValueError, 'overlapping must be True or an integer'

    return ShuffleProduct_w1w2(w1, w2)

class ShuffleProduct_w1w2(CombinatorialClass):
    def __init__(self, w1, w2):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4])
            sage: s == loads(dumps(s))
            True

        """
        self.w1 = w1
        self.w2 = w2

    def __repr__(self):
        """
        EXAMPLES:
            sage: repr(ShuffleProduct([1,2],[3,4]))
            'Shuffle product of [1, 2] and [3, 4]'

        """
        return "Shuffle product of %s and %s"%(self.w1, self.w2)

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4])
            sage: all(i in s for i in s)
            True
            sage: [2, 1, 3, 4] in s
            False
        """
        w1 = self.w1[:]
        w2 = self.w2[:]

        try:
            x = list(x)
        except TypeError:
            return False

        for _ in range(len(x)):
            try:
                letter = x.pop(0)
            except IndexError:
                return False
            if len(w1) > 0 and letter == w1[0]:
                w1.pop(0)
            elif len(w2) > 0 and letter == w2[0]:
                w2.pop(0)
            else:
                return False

        return len(x) == 0

    def count(self):
        """
        Returns the number of words in the shuffle product
        of w1 and w2.

        It is given by binomial(len(w1)+len(w2), len(w1)).

        EXAMPLES:
            sage: ShuffleProduct([2,3],[3,4]).count()
            6
         """
        return binomial(len(self.w1)+len(self.w2), len(self.w1))

    def _proc(self, vect):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4])
            sage: s._proc([0,1,0,1])
            [3, 1, 4, 2]
            sage: s._proc([1,1,0,0])
            [1, 2, 3, 4]
        """
        i1 = -1
        i2 = -1
        res = []
        for v in vect:
            if v == 1:
                i1 += 1
                res.append(self.w1[i1])
            else:
                i2 += 1
                res.append(self.w2[i2])
        return res


    def iterator(self):
        """
        Returns an iterator for the words in the
        shuffle product of w1 and w2.

        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4]).list() #indirect test
            [[1, 2, 3, 4],
             [1, 3, 2, 4],
             [1, 3, 4, 2],
             [3, 1, 2, 4],
             [3, 1, 4, 2],
             [3, 4, 1, 2]]
        """

        n1 = len(self.w1)
        n2 = len(self.w2)
        for iv in IntegerVectors(n1, n1+n2, max_part=1):
            yield self._proc(iv)



class ShuffleProduct_shifted(ShuffleProduct_w1w2):
    def __init__(self, w1, w2):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4], shifted=True)
            sage: s == loads(dumps(s))
            True
        """
        ShuffleProduct_w1w2.__init__(self, w1, [x+len(w1) for x in w2])

    def __repr__(self):
        """
        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4], shifted=True).__repr__()
            'Shifted shuffle product of [1, 2] and [5, 6]'
        """
        return "Shifted shuffle product of %s and %s"%(self.w1, self.w2)


class ShuffleProduct_overlapping_r(CombinatorialClass):
    def __init__(self, w1, w2, r):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4], overlapping=1)
            sage: s == loads(dumps(s))
            True
        """
        self.w1 = w1
        self.w2 = w2
        self.r  = r

    def __repr__(self):
        """
        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4], overlapping=1).__repr__()
            'Overlapping shuffle product of [1, 2] and [3, 4] with 1 overlaps'
        """
        return "Overlapping shuffle product of %s and %s with %s overlaps"%(self.w1, self.w2, self.r)

    def iterator(self):
        """
        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4], overlapping=1).list() #indirect doctest
            [[4, 2, 4], [1, 5, 4], [4, 4, 2], [1, 3, 6], [3, 5, 2], [3, 1, 6]]
        """
        m = len(self.w1)
        n = len(self.w2)
        r = self.r

        blank = [0]*(m+n-r)
        for iv in IntegerVectors(m, m+n-r, max_part=1):
            w = copy.copy(blank)
            filled_places = []
            unfilled_places = []
            #Fill in w1 into the iv slots
            i = 0
            for j in range(len(iv)):
                if iv[j] == 1:
                    w[j] = self.w1[i]
                    i += 1
                    filled_places.append(j)
                else:
                    unfilled_places.append(j)

            #Choose r of these filled places
            for subset in Subsets(filled_places, r):
                places_to_fill = unfilled_places + list(subset)
                places_to_fill.sort()

                #Fill in w2 into the places
                i = 0
                res = copy.copy(w)
                for j in places_to_fill:
                    res[j] += self.w2[i]
                    i += 1

                yield res

class ShuffleProduct_overlapping(CombinatorialClass):
    def __init__(self, w1, w2):
        """
        EXAMPLES:
            sage: s = ShuffleProduct([1,2],[3,4],overlapping=True)
            sage: s == loads(dumps(s))
            True
        """
        self.w1 = w1
        self.w2 = w2

    def __repr__(self):
        """
        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4],overlapping=True).__repr__()
            'Overlapping shuffle product of [1, 2] and [3, 4]'
        """
        return "Overlapping shuffle product of %s and %s"%(self.w1, self.w2)

    def iterator(self):
        """
        EXAMPLES:
            sage: ShuffleProduct([1,2],[3,4],overlapping=True).list() #indirect doctest
            [[1, 2, 3, 4],
             [1, 3, 2, 4],
             [1, 3, 4, 2],
             [3, 1, 2, 4],
             [3, 1, 4, 2],
             [3, 4, 1, 2],
             [4, 2, 4],
             [1, 5, 4],
             [4, 4, 2],
             [1, 3, 6],
             [3, 5, 2],
             [3, 1, 6],
             [4, 6]]
        """
        m = len(self.w1)
        n = len(self.w2)

        for r in range(min(m,n)+1):
            for w in ShuffleProduct_overlapping_r(self.w1, self.w2, r):
                yield w




##########################
# Symmetric group action #
##########################
def unmatched_places(w, open, close):
    """
    EXAMPLES:
        sage: from sage.combinat.word import unmatched_places
        sage: unmatched_places([2,2,2,1,1,1],2,1)
        ([], [])
        sage: unmatched_places([1,1,1,2,2,2],2,1)
        ([0, 1, 2], [3, 4, 5])
        sage: unmatched_places([], 2, 1)
        ([], [])
        sage: unmatched_places([1,2,4,6,2,1,5,3],2,1)
        ([0], [1])
        sage: unmatched_places([2,2,1,2,4,6,2,1,5,3], 2, 1)
        ([], [0, 3])
        sage: unmatched_places([3,1,1,1,2,1,2], 2, 1)
        ([1, 2, 3], [6])
    """
    lw = len(w)
    places_open = []
    places_close = []
    for i in range(lw):
        letter = w[i]
        if letter == open:
            places_open.append(i)
        elif letter == close:
            if places_open == []:
                places_close.append(i)
            else:
                places_open.pop()
    return places_close, places_open


def symmetric_group_action_on_values(w, perm):
    """
    EXAMPLES:
        sage: from sage.combinat.word import symmetric_group_action_on_values
        sage: symmetric_group_action_on_values([1,1,1],[1,3,2])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([1,1,1],[2,1,3])
        [2, 2, 2]
        sage: symmetric_group_action_on_values([1,2,1],[2,1,3])
        [2, 2, 1]
        sage: symmetric_group_action_on_values([2,2,2],[2,1,3])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([2,1,2],[2,1,3])
        [2, 1, 1]
        sage: symmetric_group_action_on_values([2,2,3,1,1,2,2,3],[1,3,2])
        [2, 3, 3, 1, 1, 2, 3, 3]
        sage: symmetric_group_action_on_values([2,1,1],[2,1])
        [2, 1, 2]
        sage: symmetric_group_action_on_values([2,2,1],[2,1])
        [1, 2, 1]
        sage: symmetric_group_action_on_values([1,2,1],[2,1])
        [2, 2, 1]
    """
    ts = sage.combinat.permutation.Permutation(perm).reduced_word()
    for j in reversed(range(len(ts))):
        r = ts[j]
        l = r + 1
        places_r, places_l = unmatched_places(w, l, r)

        #Now change the number of l's and r's in the new word
        nbl = len(places_l)
        nbr = len(places_r)
        ma = max(nbl, nbr)
        dif = ma - min(nbl, nbr)
        if ma == nbl:
            for i in range(dif):
                w[places_l[i]] = r
        else:
            for i in range(nbr-dif,ma):
                w[places_r[i]] = l
    return w

