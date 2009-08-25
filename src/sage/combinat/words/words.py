# coding=utf-8
r"""
Combinatorial classes of words.

To define a new class of words, please refer to the documentation file:
sage/combinat/words/notes/word_inheritance_howto.txt

AUTHORS:

    - Franco Saliola (2008-12-17): merged into sage
    - Sebastien Labbe (2008-12-17): merged into sage
    - Arnaud Bergeron (2008-12-17): merged into sage
    - Sebastien Labbe (2009-07-21): Improved morphism iterator (#6571).

EXAMPLES::

    sage: Words()
    Words
    sage: Words(5)
    Words over Ordered Alphabet [1, 2, 3, 4, 5]
    sage: Words('ab')
    Words over Ordered Alphabet ['a', 'b']
    sage: Words('natural numbers')
    Words over Ordered Alphabet of Natural Numbers
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
from sage.combinat.combinat import InfiniteAbstractCombinatorialClass
from sage.combinat.words.alphabet import OrderedAlphabet
from sage.misc.mrange import xmrange
from sage.rings.all import Infinity
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
import itertools

def Words(alphabet=None, length=None, finite=True, infinite=True):
    """
    Returns the combinatorial class of words of length k over an ordered
    alphabet.

    EXAMPLES::

        sage: Words()
        Words
        sage: Words(length=7)
        Finite Words of length 7
        sage: Words(5)
        Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words(5, 3)
        Finite Words of length 3 over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words(5, infinite=False)
        Finite Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words(5, finite=False)
        Infinite Words over Ordered Alphabet [1, 2, 3, 4, 5]
        sage: Words('ab')
        Words over Ordered Alphabet ['a', 'b']
        sage: Words('ab', 2)
        Finite Words of length 2 over Ordered Alphabet ['a', 'b']
        sage: Words('ab', infinite=False)
        Finite Words over Ordered Alphabet ['a', 'b']
        sage: Words('ab', finite=False)
        Infinite Words over Ordered Alphabet ['a', 'b']
        sage: Words('positive integers', finite=False)
        Infinite Words over Ordered Alphabet of Positive Integers
        sage: Words('natural numbers')
        Words over Ordered Alphabet of Natural Numbers
    """
    if isinstance(alphabet, Words_all):
        return alphabet
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
        if isinstance(alphabet, (int,Integer)):
            alphabet = OrderedAlphabet(range(1,alphabet+1))
        elif alphabet == "integers" \
                or alphabet == "positive integers" \
                or alphabet == "natural numbers":
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
    r"""
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
    def __call__(self, data=None, length=None, datatype=None, **kwds):
        r"""
        Construct a new word object with parent self.

        NOTE:

            We only check that the first 40 letters of the word are
            actually in the alphabet. This is a quick check implemented to
            test for small programming errors. Since we also support
            infinite words, we cannot really implement a more accurate
            check.

        EXAMPLES::

            sage: from itertools import count
            sage: Words()(count())
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: Words(range(10))(count())
            Traceback (most recent call last):
            ...
            ValueError: 10 not in alphabet!
            sage: Words()("abba")
            word: abba
            sage: Words("ab")("abba")
            word: abba
            sage: Words("ab")("abca")
            Traceback (most recent call last):
            ...
            ValueError: c not in alphabet!
        """
        from sage.combinat.words.word import Word
        kwds['data'] = data
        kwds['length'] = length
        kwds['datatype'] = datatype
        kwds['alphabet'] = self

        ## BACKWARD COMPATIBILITY / DEPRECATION WARNINGS.
        # In earlier versions, the first argument was ``obj``. We change
        # this to ``data`` for consistency with the Word interface.
        if kwds.has_key('obj'):
            from sage.misc.misc import deprecation
            deprecation("obj argument is deprecated, use data instead")
            kwds['data'] = kwds['obj']
            del kwds['obj']
        # In earlier versions, the third argument was ``format``, which
        # could have been: 'empty', 'letter', 'list', 'function',
        # 'iterator', 'content'.
        if kwds.has_key('format'):
            from sage.misc.misc import deprecation
            deprecation("format argument is deprecated, use datatype instead")
            datatype = kwds['format']
            del kwds['format']
        if datatype == 'empty':
            if kwds['data'] is not None:
                raise TypeError, "trying to build an empty word with something other than None"
            kwds['data'] = []
            datatype = 'list'
        elif datatype == 'letter':
            kwds['data'] = [kwds['data'],]
            datatype = 'list'
        elif datatype == 'function':
            datatype = 'callable'
        elif datatype == 'iterator':
            datatype = 'iter'
        elif datatype == 'content':
            datatype = None
        kwds['datatype'] = datatype
        # In earlier versions, the second argument was ``part``, which was
        # supposed to be a slice object.
        if isinstance(length, slice):
            kwds['part'] = length
            del kwds['length']
        if kwds.has_key('part'):
            from sage.misc.misc import deprecation
            deprecation("part argument is deprecated, use length instead")
            part = kwds['part']
            del kwds['part']
            alphabet = kwds['alphabet']
            del kwds['alphabet']
            w = Word(data=Word(**kwds)[part],alphabet=alphabet)
            return w
        ## END OF BACKWARD COMPATIBILITY / DEPRECATION WARNINGS.

        # The function Word handles the construction of the words.
        w = Word(**kwds)
        self._check(w)
        return w

    def _check(self, w, length=40):
        r"""
        Check that the first length elements are actually in the alphabet.

        NOTE:

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet(['a','b','c'])
            sage: W._check('abcabc') is None
            True
            sage: W._check('abcabcd')
            Traceback (most recent call last):
            ...
            ValueError: d not in alphabet!
            sage: W._check('abcabc'*10+'z') is None
            True
            sage: W._check('abcabc'*10+'z', length=80)
            Traceback (most recent call last):
            ...
            ValueError: z not in alphabet!
        """
        for a in itertools.islice(w, length):
            if a not in self._alphabet:
                raise ValueError, "%s not in alphabet!" % a

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: Words_all().__repr__()
            'Words'
        """
        return 'Words'

    def __contains__(self, x):
        """
        Returns True if x is contained in self.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: 2 in Words_all()
            False
            sage: [1,2] in Words_all()
            False
            sage: Words('ab')('abba') in Words_all()
            True
        """
        from sage.combinat.words.word import Word_class
        return isinstance(x, Word_class)

    class _python_object_alphabet(object):
        r"""
        The "ordered" alphabet of all Python objects.
        """
        def __contains__(self, x):
            r"""
            Returns always True.

            EXAMPLES::

                sage: W = Words()
                sage: A = W._python_object_alphabet()
                sage: 'a' in A
                True
                sage: 5 in A
                True
                sage: 'abcdef' in A
                True
                sage: [4,5] in A
                True
            """
            return True
        def __repr__(self):
            r"""
            EXAMPLES::

                sage: W = Words()
                sage: A = W._python_object_alphabet()
                sage: A
                Python objects
            """
            return "Python objects"

    _alphabet = _python_object_alphabet()

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

    def size_of_alphabet(self):
        r"""
        Returns the size of the alphabet.

        EXAMPLES::

            sage: Words().size_of_alphabet()
            +Infinity
        """
        return Infinity

    cmp_letters = cmp

    def has_letter(self, letter):
        r"""
        Returns True if the alphabet of self contains the given letter.

        INPUT:

        -  ``letter`` - a letter

        EXAMPLES::

            sage: W = Words()
            sage: W.has_letter('a')
            True
            sage: W.has_letter(1)
            True
            sage: W.has_letter({})
            True
            sage: W.has_letter([])
            True
            sage: W.has_letter(range(5))
            True
            sage: W.has_letter(Permutation([]))
            True

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet(['a','b','c'])
            sage: W.has_letter('a')
            True
            sage: W.has_letter('d')
            False
            sage: W.has_letter(8)
            False
        """
        return letter in self._alphabet

class Words_over_Alphabet(Words_all):
    def __init__(self, alphabet):
        """
        Words over Alphabet.

        INPUT:

        -  ``alphabet`` - assumed to be an instance of Alphabet, but no
           type checking is done here.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet([1,2,3])
            sage: W == loads(dumps(W))
            True

        The input alphabet must be an instance of Alphabet::

            sage: W = Words_over_Alphabet(Alphabet([1,2,3]))
            sage: W([1,2,2,3])
            word: 1223
        """
        self._alphabet = alphabet

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
        Tests whether self contains x.

        OUTPUT:
            This method returns True if x is a word of the appropriate
            length and the alphabets of the parents match. Returns False
            otherwise.

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
        from sage.combinat.words.word import Word_class
        return isinstance(x, Word_class) and x.parent().alphabet() == self.alphabet()

    def __lt__(self, other):
        r"""
        Returns whether self is a proper subset of other.

        TESTS::

            sage: Words('ab') < Words('ab')
            False
            sage: Words('ab') < Words('abc')
            True
            sage: Words('abc') < Words('ab')
            False
        """
        if not isinstance(other,Words_all):
            return NotImplemented
        return self <= other and self != other

    def __gt__(self, other):
        r"""
        Returns True if self is a proper superset of other and False otherwise.

        TESTS::

            sage: Words('ab') > Words('ab')
            False
            sage: Words('ab') > Words('abc')
            False
            sage: Words('abc') > Words('ab')
            True
        """
        if not isinstance(other,Words_all):
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
        if not isinstance(other,Words_all):
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
        if not isinstance(other,Words_all):
            return NotImplemented
        if isinstance(self, type(other)):
            return self.alphabet() >= other.alphabet()
        else:
            return False

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
        from sage.combinat.words.word import FiniteWord_class
        return isinstance(x, FiniteWord_class) and x.length() == self._n

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

    def iterate_by_length(self, l=1):
        r"""
        Returns an iterator over all the words of self of length l.

        INPUT:

        - ``l`` - integer (default: 1), the length of the desired words

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
        if not isinstance(l, (int,Integer)):
            raise TypeError, "the parameter l (=%r) must be an integer"%l
        #if l == Integer(0):
        #    yield self()
        for w in xmrange([self.size_of_alphabet()]*l):
            yield self(map(lambda x: self.alphabet().unrank(x), w))

    def __iter__(self):
        r"""
        Returns an iterator over all the words of self.

        The iterator outputs the words in lexicographic order,
        based on the order of the letters in the alphabet.

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
        for l in itertools.count():
            for w in self.iterate_by_length(l):
                yield w

    def iter_morphisms(self, l=None, codomain=None, min_length=1):
        r"""
        Iterate over all morphisms with domain ``self`` and the given
        codmain.

        INPUT:

        - ``l`` -- list of nonnegative integers (default: None). The length
          of the list must be the number of letters in the alphabet, and
          the `i`-th integer of ``l`` determines the length of the word
          mapped to by `i`-th letter of the (ordered) alphabet. If ``l`` is
          ``None``, then the method iterates through all morphisms.

        - ``codomain`` -- (default: None) a combinatorial class of words.
          By default, ``codomain`` is ``self``.

        - ``min_length`` -- (default: 1) nonnegative integer. If ``l`` is
          not specified, then iterate through all the morphisms where the
          length of the images of each letter in the alphabet is at least
          ``min_length``. This is ignored if ``l`` is not ``None``.

        OUTPUT:

            iterator

        EXAMPLES:

        Iterator over all non-erasing morphisms::

            sage: W = Words('ab')
            sage: it = W.iter_morphisms()
            sage: for _ in range(7): print it.next()
            WordMorphism: a->a, b->a
            WordMorphism: a->a, b->b
            WordMorphism: a->b, b->a
            WordMorphism: a->b, b->b
            WordMorphism: a->aa, b->a
            WordMorphism: a->aa, b->b
            WordMorphism: a->ab, b->a

        Iterator over all morphisms including erasing morphisms::

            sage: W = Words('ab')
            sage: it = W.iter_morphisms(min_length=0)
            sage: for _ in range(7): print it.next()
            WordMorphism: a->, b->
            WordMorphism: a->a, b->
            WordMorphism: a->b, b->
            WordMorphism: a->, b->a
            WordMorphism: a->, b->b
            WordMorphism: a->aa, b->
            WordMorphism: a->ab, b->

        Iterator over morphisms with specific image lengths::

            sage: for m in W.iter_morphisms([0, 0]): print m
            WordMorphism: a->, b->
            sage: for m in W.iter_morphisms([0, 1]): print m
            WordMorphism: a->, b->a
            WordMorphism: a->, b->b
            sage: for m in W.iter_morphisms([2, 1]): print m
            WordMorphism: a->aa, b->a
            WordMorphism: a->aa, b->b
            WordMorphism: a->ab, b->a
            WordMorphism: a->ab, b->b
            WordMorphism: a->ba, b->a
            WordMorphism: a->ba, b->b
            WordMorphism: a->bb, b->a
            WordMorphism: a->bb, b->b
            sage: for m in W.iter_morphisms([2, 2]): print m
            WordMorphism: a->aa, b->aa
            WordMorphism: a->aa, b->ab
            WordMorphism: a->aa, b->ba
            WordMorphism: a->aa, b->bb
            WordMorphism: a->ab, b->aa
            WordMorphism: a->ab, b->ab
            WordMorphism: a->ab, b->ba
            WordMorphism: a->ab, b->bb
            WordMorphism: a->ba, b->aa
            WordMorphism: a->ba, b->ab
            WordMorphism: a->ba, b->ba
            WordMorphism: a->ba, b->bb
            WordMorphism: a->bb, b->aa
            WordMorphism: a->bb, b->ab
            WordMorphism: a->bb, b->ba
            WordMorphism: a->bb, b->bb

        The codomain may be specified as well::

            sage: Y = Words('xyz')
            sage: for m in W.iter_morphisms([0, 2], codomain=Y): print m
            WordMorphism: a->, b->xx
            WordMorphism: a->, b->xy
            WordMorphism: a->, b->xz
            WordMorphism: a->, b->yx
            WordMorphism: a->, b->yy
            WordMorphism: a->, b->yz
            WordMorphism: a->, b->zx
            WordMorphism: a->, b->zy
            WordMorphism: a->, b->zz
            sage: for m in Y.iter_morphisms([0,2,1], codomain=W): print m
            WordMorphism: x->, y->aa, z->a
            WordMorphism: x->, y->aa, z->b
            WordMorphism: x->, y->ab, z->a
            WordMorphism: x->, y->ab, z->b
            WordMorphism: x->, y->ba, z->a
            WordMorphism: x->, y->ba, z->b
            WordMorphism: x->, y->bb, z->a
            WordMorphism: x->, y->bb, z->b
            sage: it = W.iter_morphisms(codomain=Y)
            sage: for _ in range(10): print it.next()
            WordMorphism: a->x, b->x
            WordMorphism: a->x, b->y
            WordMorphism: a->x, b->z
            WordMorphism: a->y, b->x
            WordMorphism: a->y, b->y
            WordMorphism: a->y, b->z
            WordMorphism: a->z, b->x
            WordMorphism: a->z, b->y
            WordMorphism: a->z, b->z
            WordMorphism: a->xx, b->x

        TESTS::

            sage: list(W.iter_morphisms([1,0]))
            [Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['a', 'b'], Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['a', 'b']]
            sage: list(W.iter_morphisms([0,0], codomain=Y))
            [Morphism from Words over Ordered Alphabet ['a', 'b'] to Words over Ordered Alphabet ['x', 'y', 'z']]
            sage: list(W.iter_morphisms([0, 1, 2]))
            Traceback (most recent call last):
            ...
            TypeError: l (=[0, 1, 2]) must be an iterable of 2 integers
            sage: list(W.iter_morphisms([0, 'a']))
            Traceback (most recent call last):
            ...
            TypeError: l (=[0, 'a']) must be an iterable of 2 integers
            sage: list(W.iter_morphisms([0, 1], codomain='a'))
            Traceback (most recent call last):
            ...
            TypeError: codomain (=a) must be an instance of Words_over_OrderedAlphabet
        """
        n = self.size_of_alphabet()
        # create an iterable of compositions (all "compositions" if l is
        # None, or [l] otherwise)
        if l is None:
            from sage.combinat.integer_list import IntegerListsLex
            compositions = IntegerListsLex(itertools.count(), \
                    length=n, min_part = max(0,min_length))
        else:
            l = list(l)
            if not len(l) == n or not \
                    all(isinstance(a, (int,Integer)) for a in l):
                raise TypeError, \
                    "l (=%s) must be an iterable of %s integers" %(l, n)
            compositions = [l]

        # set the codomain
        if codomain is None:
            codomain = self
        elif not isinstance(codomain, Words_over_OrderedAlphabet):
            raise TypeError, "codomain (=%s) must be an instance of Words_over_OrderedAlphabet"%codomain

        # iterate through the morphisms
        from sage.combinat.words.morphism import WordMorphism
        for composition in compositions:
            cuts = [0] + composition
            for i in range(1,len(cuts)):
                cuts[i] += cuts[i-1]
            s = cuts[-1] # same but better than s = sum(composition)
            for big_word in codomain.iterate_by_length(s):
                d = {}
                i = 0
                for a in self.alphabet():
                    d[a] = big_word[cuts[i]:cuts[i+1]]
                    i += 1
                yield WordMorphism(d, codomain=codomain)

    def cmp_letters(self, letter1, letter2):
        r"""
        Returns a negative number, zero or a positive number if
        ``letter1`` < ``letter2``, ``letter1`` == ``letter2`` or
        ``letter1`` > ``letter2`` respectively.

        INPUT:

        - ``letter1`` - a letter in the alphabet
        - ``letter2`` - a letter in the alphabet

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_OrderedAlphabet
            sage: from sage.combinat.words.words import OrderedAlphabet
            sage: A = OrderedAlphabet('woa')
            sage: W = Words_over_OrderedAlphabet(A)
            sage: W.cmp_letters('w','a')
            -2
            sage: W.cmp_letters('w','o')
            -1
            sage: W.cmp_letters('w','w')
            0
        """
        return self._alphabet.rank(letter1) - self._alphabet.rank(letter2)

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

    def __repr__(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Words('ab', infinite=False).__repr__()
            "Finite Words over Ordered Alphabet ['a', 'b']"
        """
        return "Finite Words over %s" % self.alphabet()

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
            sage: Words([0,1])([1,0,1]) in W
            True
            sage: Words([1,0])([1,0,1]) in W
            False
            sage: W([1,0,1]) in W
            True
            sage: Word([2,0]) in W
            False
        """
        if super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).__contains__(x) and x.length() == self._length:
            return True
        else:
            return False

    def __repr__(self):
        """
        TESTS::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.OrderedAlphabet([1,0])
            sage: FiniteWords_length_k_over_OrderedAlphabet(A,3).__repr__()
            'Finite Words of length 3 over Ordered Alphabet [1, 0]'
        """
        from sage.combinat.words.word_options import word_options
        if word_options['old_repr']:
            return "Finite Words over %s of length %s"%(self.alphabet(), self._length)
        return "Finite Words of length %s over %s"%(self._length, self.alphabet())

    def cardinality(self):
        r"""
        Returns the number of words of length `n` from alphabet.

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

    def list(self):
        r"""
        Returns a list of all the words contained in self.

        EXAMPLES::

            sage: Words(0,0).list()
            [word: ]
            sage: Words(5,0).list()
            [word: ]
            sage: Words(['a','b','c'],0).list()
            [word: ]
            sage: Words(5,1).list()
            [word: 1, word: 2, word: 3, word: 4, word: 5]
            sage: Words(['a','b','c'],2).list()
            [word: aa, word: ab, word: ac, word: ba, word: bb, word: bc, word: ca, word: cb, word: cc]
        """
        return list(self)

    def __iter__(self):
        """
        Returns an iterator for all of the words of length k from
        ``self.alphabet()``. The iterator outputs the words in lexicographic
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

###########################################################################
##### DEPRECATION WARNINGS ################################################
##### Added July 2009 #####################################################
###########################################################################

def is_Words(obj):
    r"""
    Returns True if obj is a word set and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.words import is_Words
        sage: is_Words(33)
        doctest:1: DeprecationWarning: is_Words is deprecated, use isinstance(your_object, Words_all) instead!
        False
        sage: is_Words(Words('ab'))
        True
    """
    from sage.misc.misc import deprecation
    deprecation("is_Words is deprecated, use isinstance(your_object, Words_all) instead!")
    return isinstance(obj, Words_all)

