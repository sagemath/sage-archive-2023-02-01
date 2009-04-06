# coding=utf-8
"""
Words
"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import count
from weakref import WeakValueDictionary
from sage.structure.sage_object import SageObject
from sage.misc.latex import latex
from sage.rings.integer import Integer
from sage.combinat.words.alphabet import OrderedAlphabet, OrderedAlphabet_class
from sage.combinat.words.word_content import is_WordContent, BuildWordContent
from sage.combinat.words.utils import *
from sage.combinat.combinat import InfiniteAbstractCombinatorialClass
from sage.rings.all import Infinity
from sage.misc.mrange import xmrange

def is_Words(obj):
    r"""
    Returns True if obj is a word set and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.words import is_Words
        sage: is_Words(33)
        False
        sage: is_Words(Words('ab'))
        True
    """
    return isinstance(obj, Words_all)

def Words(alphabet=None, length=None, finite=True, infinite=True):
    """
    Returns the combinatorial class of words of length k over an
    ordered alphabet.

    EXAMPLES::

        sage: Words()
        Words
        sage: Words(length=7)
        Finite Words of length 7
        sage: Words(5)
        Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words(5, 3)
        Finite Words over Ordered Alphabet [1, 2, 3, 4, 5] of length 3
        sage: Words(5, infinite=False)
        Finite Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words(5, finite=False)
        Infinite Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words('ab')
        Words over Ordered Alphabet ['a', 'b']
        sage: Words('ab', 2)
        Finite Words over Ordered Alphabet ['a', 'b'] of length 2
        sage: Words('ab', infinite=False)
        Finite Words over Ordered Alphabet ['a', 'b']
        sage: Words('ab', finite=False)
        Infinite Words over Ordered Alphabet ['a', 'b']
        sage: Words('positive integers', finite=False)
        Infinite Words over Ordered Alphabet of Positive Integers
        sage: Words('natural numbers')
        Words over Ordered Alphabet of Natural Numbers
    """
    if alphabet is None:
        if length is None:
            if finite and infinite:
                return Words_all()
            elif finite:
                raise NotImplementedError
            else:
                raise NotImplementedError
        elif isinstance(length, (int,Integer)) and finite:
            return Words_n(length)
    else:
        if not isinstance(alphabet, OrderedAlphabet_class):
            if isinstance(alphabet, (int,Integer)):
                alphabet = OrderedAlphabet(range(1,alphabet+1))
            else:
                if alphabet == "positive integers" or alphabet == "natural numbers":
                    alphabet = OrderedAlphabet(name=alphabet)
                else:
                    alphabet = OrderedAlphabet(alphabet)
        if length is None:
            if finite and infinite:
                return Words_over_OrderedAlphabet(alphabet)
            elif finite:
                return FiniteWords_over_OrderedAlphabet(alphabet)
            else:
                return InfiniteWords_over_OrderedAlphabet(alphabet)
        elif isinstance(length, (int,Integer)):
                return FiniteWords_length_k_over_OrderedAlphabet(alphabet, length)
    raise ValueError, "do not know how to make a combinatorial class of words from your input"

class Words_all(InfiniteAbstractCombinatorialClass):
    """
    TESTS::

        sage: from sage.combinat.words.words import Words_all
        sage: list(Words_all())
        Traceback (most recent call last):
        ...
        NotImplementedError
        sage: Words_all().list()
        Traceback (most recent call last):
        ...
        NotImplementedError: infinite list
        sage: Words_all().cardinality()
        +Infinity
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: Words_all().__repr__()
            'Words'
        """
        return "Words"

    def __contains__(self, x):
        """
        Return True if x is contained in self.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: 2 in Words_all()
            False
            sage: [1,2] in Words_all()
            False
            sage: Words('ab')('abba') in Words_all()
            True
        """
        from sage.combinat.words.word import is_Word
        return is_Word(x)



class Words_n(Words_all):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: w = Words_n(3)
            sage: w == loads(dumps(w))
            True
        """
        self._n = n

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: Words_n(3).__repr__()
            'Finite Words of length 3'
        """
        return "Finite Words of length %s"%self._n

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: 2 in Words_n(3)
            False
            sage: [1,'a',3] in Words_n(3)
            False
            sage: [1,2] in Words_n(3)
            False
            sage: "abc" in Words_n(3)
            False
            sage: Words("abc")("ababc") in Words_n(3)
            False
            sage: Words([0,1])([1,0,1]) in Words_n(3)
            True
        """
        from sage.combinat.words.word import is_FiniteWord
        return is_FiniteWord(x) and len(x) == self._n

class Words_over_Alphabet(Words_all):
    def __init__(self, alphabet):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: w = Words_over_Alphabet([1,2,3])
            sage: w == loads(dumps(w))
            True
        """
        self._alphabet = alphabet
        from sage.combinat.words.word import FiniteWord_over_Alphabet, InfiniteWord_over_Alphabet
        self._finite_word_class = FiniteWord_over_Alphabet
        self._infinite_word_class = InfiniteWord_over_Alphabet

    def alphabet(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet([1,2,3])
            sage: W.alphabet()
            [1, 2, 3]
            sage: from sage.combinat.words.words import OrderedAlphabet
            sage: W = Words_over_Alphabet(OrderedAlphabet('ab'))
            sage: W.alphabet()
            Ordered Alphabet ['a', 'b']
        """
        return self._alphabet

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: Words_over_Alphabet([1,2,3]).__repr__()
            'Words over [1, 2, 3]'
        """
        return "Words over %s"%self._alphabet

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: from sage.combinat.words.words import OrderedAlphabet
            sage: A = OrderedAlphabet('ab')
            sage: Words(A)('abba') in Words_over_Alphabet(A)
            True
            sage: Words(A)('aa') in Words_over_Alphabet(A)
            True
            sage: Words('a')('aa') in Words_over_Alphabet(A)
            False
            sage: 2 in Words_over_Alphabet([1,2,3])
            False
            sage: [2] in Words_over_Alphabet([1,2,3])
            False
            sage: [1, 'a'] in Words_over_Alphabet([1,2,3])
            False
        """
        from sage.combinat.words.word import is_Word
        return is_Word(x) and hasattr(x, "alphabet") and x.alphabet() == self.alphabet()

    def __call__(self, obj=None, part=slice(None), format=None):
        r"""
        TESTS::

            sage: from itertools import repeat, count, imap
            sage: W = Words('ab')
            sage: W(format='empty')
            word:
            sage: W(['a', 'b'], format='empty')
            Traceback (most recent call last):
              ...
            TypeError: trying to build an empty word with something other than None
            sage: W(['a', 'b', 'b'])
            word: abb
            sage: W(['a', 'b', 'b'], format='list')
            word: abb
            sage: W(10, format='list')
            Traceback (most recent call last):
              ...
            TypeError: trying to build a word backed by a list with a sequence not providing the required operations
            sage: Words([0, 1])(lambda x: x%2)
            Infinite word over [0, 1]
            sage: W = Words([0, 1, 2])
            sage: W(lambda x: x%3, slice(0,10))
            word: 0120120120
            sage: W(repeat(1))
            Infinite word over [0, 1, 2]
            sage: Words(xrange(10))(count(), slice(10))
            word: 0123456789
            sage: Words(xrange(10))(count(), slice(10, 0, -2))
            word: 97531
            sage: W(1)
            word: 1
            sage: W('a', format='letter')
            Traceback (most recent call last):
            ...
            ValueError: letter ('a') not in alphabet
            sage: W([0, 1])
            word: 01
        """
        from sage.combinat.words.word import is_Word
        if format is None:
            if is_WordContent(obj):
                format = 'content'
            elif obj in self.alphabet():
                format = 'letter'

        if format == 'letter':
            if obj not in self.alphabet():
                raise ValueError, "letter (%r) not in alphabet" % obj
            obj = (obj,)
            format = 'list'

        if not slice_ok(part):
            raise TypeError, "part is not a slice or has wrong types for elements"

        if format == 'content':
            if not is_WordContent(obj):
                raise TypeError, "trying to build a word based on raw content with a non-content"
            content = obj
        elif is_Word(obj) and obj.parent() is self:
            return obj[part]
        else:
            content = BuildWordContent(obj, self.alphabet().rank, format=format, part=part)
        if haslen(content):
            return self._finite_word_class(self, content)
        else:
            return self._infinite_word_class(self, content)

    def size_of_alphabet(self):
        r"""
        Returns the size of the alphabet.

        EXAMPLES::

            sage: Words('abcdef').size_of_alphabet()
            6
            sage: Words('').size_of_alphabet()
            0
        """
        return self.alphabet().cardinality()

    def __lt__(self, other):
        r"""
        Returns True if self is a proper subset of other and False
        otherwise.

        TESTS::

            sage: Words('ab') < Words('ab')
            False
            sage: Words('ab') < Words('abc')
            True
            sage: Words('abc') < Words('ab')
            False
        """
        if not is_Words(other):
            return NotImplemented
        return self <= other and self != other

    def __gt__(self, other):
        r"""
        Returns True if self is a proper superset of other and False
        otherwise.

        TESTS::

            sage: Words('ab') > Words('ab')
            False
            sage: Words('ab') > Words('abc')
            False
            sage: Words('abc') > Words('ab')
            True
        """
        if not is_Words(other):
            return NotImplemented
        return self >= other and self != other

    def __le__(self, other):
        r"""
        Returns True if self is a subset of other and False otherwise.

        TESTS::

            sage: Words('ab') <= Words('ab')
            True
            sage: Words('ab') <= Words('abc')
            True
            sage: Words('abc') <= Words('ab')
            False
        """
        if not is_Words(other):
            return NotImplemented
        if isinstance(other, type(self)):
            return self.alphabet() <= other.alphabet()
        else:
            return False

    def __ge__(self, other):
        r"""
        Returns True if self is a superset of other and False otherwise.

        TESTS::

            sage: Words('ab') >= Words('ab')
            True
            sage: Words('ab') >= Words('abc')
            False
            sage: Words('abc') >= Words('ab')
            True
        """
        if not is_Words(other):
            return NotImplemented
        if isinstance(self, type(other)):
            return self.alphabet() >= other.alphabet()
        else:
            return False

class Words_over_OrderedAlphabet(Words_over_Alphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import OrderedAlphabet
            sage: A = OrderedAlphabet("abc")
            sage: W = Words_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True
        """
        super(Words_over_OrderedAlphabet, self).__init__(alphabet)
        from sage.combinat.words.word import FiniteWord_over_OrderedAlphabet, InfiniteWord_over_OrderedAlphabet
        self._finite_word_class = FiniteWord_over_OrderedAlphabet
        self._infinite_word_class = InfiniteWord_over_OrderedAlphabet

    def iterate_by_length(self, l=1):
        r"""
        Returns an iterator over all the words of self of length l.

        INPUT:


        -  ``l`` - integer (default: 1) the length of the
           desired words


        EXAMPLES::

            sage: W = Words('ab')
            sage: list(W.iterate_by_length(1))
            [word: a, word: b]
            sage: list(W.iterate_by_length(2))
            [word: aa, word: ab, word: ba, word: bb]
            sage: list(W.iterate_by_length(3))
            [word: aaa,
             word: aab,
             word: aba,
             word: abb,
             word: baa,
             word: bab,
             word: bba,
             word: bbb]
            sage: list(W.iterate_by_length('a'))
            Traceback (most recent call last):
            ...
            TypeError: the parameter l (='a') must be an integer
        """
        if not isint(l):
            raise TypeError, "the parameter l (=%r) must be an integer"%l
        if l == Integer(0):
            yield self()
        for w in xmrange([self.size_of_alphabet()]*l):
            yield self(map(lambda x: self.alphabet().unrank(x), w))

    def __iter__(self):
        r"""
        Returns an iterator over all the words of self. The iterator
        outputs the words in lexicographic order, based on the order of the
        letters in the alphabet.

        EXAMPLES::

            sage: W = Words([4,5])
            sage: for w in W:
            ...     if len(w)>3:
            ...         break
            ...     else:
            ...         print w
            ...
            word:
            word: 4
            word: 5
            word: 44
            word: 45
            word: 54
            word: 55
            word: 444
            word: 445
            word: 454
            word: 455
            word: 544
            word: 545
            word: 554
            word: 555
            sage: W = Words([5,4])
            sage: for w in W:
            ...     if len(w)>3:
            ...         break
            ...     else:
            ...         print w
            ...
            word:
            word: 5
            word: 4
            word: 55
            word: 54
            word: 45
            word: 44
            word: 555
            word: 554
            word: 545
            word: 544
            word: 455
            word: 454
            word: 445
            word: 444
        """
        for l in count():
            for w in self.iterate_by_length(l):
                yield w

    def iter_morphisms(self, l):
        r"""
        Iterate over all endomorphisms `\varphi` of self satisfying
        `|\varphi(a[i])|=l[i]`, where a[i] is the i-th letter of
        self.

        INPUT:


        -  ``l`` - list of integers such that len(l) ==
           self.size_of_alphabet()


        OUTPUT:


        -  ``iterator outputting WordMorphism`` - outputs the
           endomorphisms


        EXAMPLES::

            sage: W = Words('ab')
            sage: map(str, W.iter_morphisms([2, 1]))
            ['WordMorphism: a->aa, b->a',
             'WordMorphism: a->aa, b->b',
             'WordMorphism: a->ab, b->a',
             'WordMorphism: a->ab, b->b',
             'WordMorphism: a->ba, b->a',
             'WordMorphism: a->ba, b->b',
             'WordMorphism: a->bb, b->a',
             'WordMorphism: a->bb, b->b']
            sage: map(str, W.iter_morphisms([2, 2]))
            ['WordMorphism: a->aa, b->aa',
             'WordMorphism: a->aa, b->ab',
             'WordMorphism: a->aa, b->ba',
             'WordMorphism: a->aa, b->bb',
             'WordMorphism: a->ab, b->aa',
             'WordMorphism: a->ab, b->ab',
             'WordMorphism: a->ab, b->ba',
             'WordMorphism: a->ab, b->bb',
             'WordMorphism: a->ba, b->aa',
             'WordMorphism: a->ba, b->ab',
             'WordMorphism: a->ba, b->ba',
             'WordMorphism: a->ba, b->bb',
             'WordMorphism: a->bb, b->aa',
             'WordMorphism: a->bb, b->ab',
             'WordMorphism: a->bb, b->ba',
             'WordMorphism: a->bb, b->bb']
            sage: map(str, W.iter_morphisms([0, 0]))
            ['WordMorphism: a->, b->']
            sage: map(str, W.iter_morphisms([0, 1]))
            ['WordMorphism: a->, b->a', 'WordMorphism: a->, b->b']
            sage: list(W.iter_morphisms([1,0]))
            [Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['a', 'b'], Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['a', 'b']]

        TESTS::

            sage: list(W.iter_morphisms([0, 1, 2]))
            Traceback (most recent call last):
            ...
            TypeError: l (=[0, 1, 2]) must be a list of 2 integers
            sage: list(W.iter_morphisms([0, 'a']))
            Traceback (most recent call last):
            ...
            TypeError: l (=[0, 'a']) must be a list of 2 integers
        """
        if not isinstance(l, list) or not len(l) == self.size_of_alphabet() \
            or not all(map(isint,l)):
            raise TypeError, "l (=%s) must be a list of %s integers" \
                             %(l, self.size_of_alphabet())

        from sage.combinat.words.morphism import WordMorphism

        cuts = [0] + l
        for i in range(1,len(cuts)):
            cuts[i] += cuts[i-1]

        s = cuts[-1] # same but better than s = sum(l)
        for big_word in self.iterate_by_length(s):
            d = {}
            i = 0
            for a in self.alphabet():
                d[a] = big_word[cuts[i]:cuts[i+1]]
                i += 1
            yield WordMorphism(d, codomain=self)

class InfiniteWords_over_OrderedAlphabet(Words_over_OrderedAlphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import InfiniteWords_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import OrderedAlphabet
            sage: A = OrderedAlphabet("abc")
            sage: W = InfiniteWords_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True
        """
        super(InfiniteWords_over_OrderedAlphabet, self).__init__(alphabet)
        from sage.combinat.words.word import InfiniteWord_over_OrderedAlphabet
        self._infinite_word_class = InfiniteWord_over_OrderedAlphabet
        self._finite_word_class = None

    def __repr__(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Words('ab', infinite=False).__repr__()
            "Finite Words over Ordered Alphabet ['a', 'b']"
        """
        return "Infinite Words over %s" % self.alphabet()

class FiniteWords_over_OrderedAlphabet(Words_over_OrderedAlphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import FiniteWords_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import OrderedAlphabet
            sage: A = OrderedAlphabet("abc")
            sage: W = FiniteWords_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True
        """
        super(FiniteWords_over_OrderedAlphabet, self).__init__(alphabet)
        from sage.combinat.words.word import FiniteWord_over_OrderedAlphabet
        self._infinite_word_class = None
        self._finite_word_class = FiniteWord_over_OrderedAlphabet

    def __repr__(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Words('ab', infinite=False).__repr__()
            "Finite Words over Ordered Alphabet ['a', 'b']"
        """
        return "Finite Words over %s" % self.alphabet()

    def __call__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import FiniteWords_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import OrderedAlphabet
            sage: A = OrderedAlphabet("abc")
            sage: W = FiniteWords_over_OrderedAlphabet(A)
            sage: W('abba')
            word: abba
        """
        from sage.combinat.words.word import is_FiniteWord
        w = super(FiniteWords_over_OrderedAlphabet, self).__call__(*args, **kwds)
        if not is_FiniteWord(w):
            raise TypeError, "infinite words are unsupported for this set"
        return w

class FiniteWords_length_k_over_OrderedAlphabet(FiniteWords_over_OrderedAlphabet):
    def __init__(self, alphabet, length):
        """
        TESTS::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.OrderedAlphabet([0,1])
            sage: W = FiniteWords_length_k_over_OrderedAlphabet(A, 3)
            sage: W == loads(dumps(W))
            True
        """
        super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).__init__(alphabet)
        self._length = length

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.OrderedAlphabet([0,1])
            sage: W = FiniteWords_length_k_over_OrderedAlphabet(A, 3)
            sage: [1,2,3] in W
            False
            sage: [1,2] in W
            False
            sage: Word([1,0,1]) in W
            True
            sage: Words([1,0])([1,0,1]) in W
            False
            sage: W([1,0,1]) in W
            True
            sage: Word([2,0]) in W
            False
        """
        if super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).__contains__(x) and len(x) == self._length:
            return True
        else:
            return False

    def __repr__(self):
        """
        TESTS::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.OrderedAlphabet([1,0])
            sage: FiniteWords_length_k_over_OrderedAlphabet(A,3).__repr__()
            'Finite Words over Ordered Alphabet [1, 0] of length 3'
        """
        return "Finite Words over %s of length %s"%(self.alphabet(), self._length)

    def cardinality(self):
        r"""
        Returns the number of words of length n from alphabet.

        EXAMPLES::

            sage: Words(['a','b','c'], 4).cardinality()
            81
            sage: Words(3, 4).cardinality()
            81
            sage: Words(0,0).cardinality()
            1
            sage: Words(5,0).cardinality()
            1
            sage: Words(['a','b','c'],0).cardinality()
            1
            sage: Words(0,1).cardinality()
            0
            sage: Words(5,1).cardinality()
            5
            sage: Words(['a','b','c'],1).cardinality()
            3
            sage: Words(7,13).cardinality()
            96889010407
            sage: Words(['a','b','c','d','e','f','g'],13).cardinality()
            96889010407
        """
        n = self.size_of_alphabet()
        return n**self._length

    def __iter__(self):
        """
        Returns an iterator for all of the words of length k from
        self.alphabet(). The iterator outputs the words in lexicographic
        order, with respect to the ordering of the alphabet.

        TESTS::

            sage: [w for w in Words(['a', 'b'], 2)]
            [word: aa, word: ab, word: ba, word: bb]
            sage: [w for w in Words(['b', 'a'], 2)]
            [word: bb, word: ba, word: ab, word: aa]
            sage: [w for w in Words(['a', 'b'], 0)]
            [word: ]
            sage: [w for w in Words([], 3)]
            []
        """
        return super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).iterate_by_length(self._length)

    def iterate_by_length(self, length):
        r"""
        All words in this class are of the same length, so use iterator
        instead.

        TESTS::

            sage: W = Words(['a', 'b'], 2)
            sage: list(W.iterate_by_length(2))
            [word: aa, word: ab, word: ba, word: bb]
            sage: list(W.iterate_by_length(1))
            []
        """
        if length == self._length:
            return iter(self)
        else:
            return iter([])
