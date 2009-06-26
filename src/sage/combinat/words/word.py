# coding=utf-8
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Amy Glen <amy.glen@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Words of all kinds

AUTHORS:

- Arnaud Bergeron

- Amy Glen

- Sébastien Labbé

- Franco Saliola

EXAMPLES:

We define the set of words over 'a' and 'b'.

::

    sage: W = Words('ab'); W
    Words over Ordered Alphabet ['a', 'b']

Then we can build some words in the set::

    sage: W('abba')
    word: abba

If you try to use letters that are not the alphabet of the set you get an error::

    sage: W([1, 2])
    Traceback (most recent call last):
    ...
    IndexError: letter not in alphabet: 1

You can also build infinite words backed by a function or an iterator!

::

    sage: def f(n):
    ...     if n % 2 == 1:
    ...         return 'a'
    ...     else:
    ...         return 'b'
    ...
    sage: W(f)
    Infinite word over ['a', 'b']
"""
from itertools import tee, islice, ifilter, ifilterfalse, imap, izip, \
                      starmap, count, repeat, dropwhile, chain, cycle
from sage.structure.sage_object import SageObject
from sage.sets.set import Set, is_Set
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.misc.latex import latex
from sage.combinat.words.alphabet import OrderedAlphabet_class
from sage.combinat.words.words import Words, is_Words
from sage.combinat.words.word_content import BuildWordContent, is_WordContent
from sage.combinat.words.utils import *
from sage.combinat.partition import Partition, Partitions
from sage.combinat.combinat import CombinatorialClass
from sage.combinat.permutation import Permutation, Permutation_class
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
import copy

word_options = {'identifier':'word: ', 'display':'string', 'truncate':True, 'truncate_length':40,
        'letter_separator':','}

def WordOptions(**kwargs):
    """
    Sets the global options for elements of the word class. The
    defaults are for words to be displayed in list notation.

    INPUT:


    -  ``display`` - 'string' (default), or 'list', words
       are displayed in string or list notation.

    -  ``truncate`` - boolean (default: True), whether to
       truncate the string output of long words (see truncate_length
       below).

    -  ``truncate_length`` - integer (default: 40), if the
       length of the word is greater than this integer, then the word is
       truncated.

    -  ``letter_separator`` - (string, default: ",") if
       the string representation of letters have length greater than 1,
       then the letters are separated by this string in the string
       representation of the word.


    If no parameters are set, then the function returns a copy of the
    options dictionary.

    EXAMPLES::

        sage: w = Word([2,1,3,12])
        sage: u = Word("abba")
        sage: WordOptions(display='list')
        sage: w
        word: [2, 1, 3, 12]
        sage: u
        word: ['a', 'b', 'b', 'a']
        sage: WordOptions(display='string')
        sage: w
        word: 2,1,3,12
        sage: u
        word: abba
    """
    global word_options
    if kwargs == {}:
        return copy.copy(word_options)

    if 'display' in kwargs:
        if kwargs['display'] not in ['list', 'string']:
            raise ValueError, "display must be either 'list' or 'string'"
        else:
            word_options['display'] = kwargs['display']
    elif 'truncate' in kwargs:
        if not isinstance(kwargs['truncate'], bool):
            raise ValueError, "truncate must be True or False"
        else:
            word_options['truncate'] = kwargs['truncate']
    elif 'truncate_length' in kwargs:
        if not isinstance(kwargs['truncate_length'], (int,Integer)) or kwargs['truncate_length'] <= 0:
            raise ValueError, "truncate_length must be a positive integer"
        else:
            word_options['truncate_length'] = kwargs['truncate_length']
    elif 'letter_separator' in kwargs:
        if not isinstance(kwargs['letter_separator'], str):
            raise ValueError, "letter_separator must be a string"
        else:
            word_options['letter_separator'] = kwargs['letter_separator']
    elif 'identifier' in kwargs:
        if not isinstance(kwargs['identifier'], str):
            raise ValueError, "identifier must be a string"
        else:
            word_options['identifier'] = kwargs['identifier']

def Word(data=None, alphabet=None):
    r"""
    Returns a word over the given alphabet.

    INPUT:


    -  ``data`` - An iterable (list, string, iterator) or a
       function defined on nonnegative integers that yields the letters of
       the word.

    -  ``alphabet`` - An iterable yielding letters in their
       comparative order. If alphabet is None and data is a string or
       list, then the letters appearing in data, ordered by Python's
       builtin sorted function, is used as alphabet.


    OUTPUT: A FiniteWord_over_OrderedAlphabet or
    InfiniteWord_over_OrderedAlphabet object.

    NOTE: Essentially, this function just returns
    Words(alphabet)(data).

    EXAMPLES: 0. The empty word.

    ::

        sage: Word()
        word:

    1. Finite words from strings.

    ::

        sage: Word("abbabaab")
        word: abbabaab

    ::

        sage: Word("abbabaab", alphabet="abc")
        word: abbabaab

    2. Finite words from lists.

    ::

        sage: Word(["a", "b", "b", "a", "b", "a", "a", "b"])
        word: abbabaab

    ::

        sage: Word([0,1,1,0,1,0])
        word: 011010

    3. Finite words from functions.

    ::

        sage: f = lambda n : n % 2
        sage: Word(f, [0, 1])[:17]
        word: 01010101010101010

    4. Finite words from iterators.

    ::

        sage: def tmword():
        ...    thuemorse = WordMorphism('a->ab,b->ba')
        ...    w = thuemorse('a')[:]
        ...    i = 0
        ...    while w:
        ...        for x in thuemorse(w[i]):
        ...            yield x
        ...        else:
        ...            w *= thuemorse(w[i+1])
        ...            i += 1
        sage: Word(tmword(), alphabet="ab")[:32]
        word: abbabaabbaababbabaababbaabbabaab

    5. Infinite words from functions.

    ::

        sage: f = lambda n : n % 2
        sage: Word(f, alphabet=[0, 1])
        Infinite word over [0, 1]

    6. Infinite words from iterators.

    ::

        sage: def tmword():
        ...    thuemorse = WordMorphism('a->ab,b->ba')
        ...    w = thuemorse('a')[:]
        ...    i = 0
        ...    while w:
        ...        for x in thuemorse(w[i]):
        ...            yield x
        ...        else:
        ...            w *= thuemorse(w[i+1])
        ...            i += 1
        sage: Word(tmword(), alphabet="ab")
        Infinite word over ['a', 'b']

    TESTS::

        sage: f = lambda n : n % 2
        sage: Word(f)
        Traceback (most recent call last):
          ...
        TypeError: alphabet is required for words not defined by lists or strings
        sage: Word(f, alphabet=[1])
        Infinite word over [1]
    """
    # If the alphabet is not specified, then try to build one from data.
    if alphabet is None:
        if data is None:
            alphabet = []
        elif isinstance(data, (str,list)):
            alphabet = sorted(set(data))
        else:
            raise TypeError, "alphabet is required for words not defined by lists or strings"
    # Return the word object
    return Words(alphabet)(data)

def is_Word(obj):
    r"""
    Returns True if obj is a word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_Word
        sage: is_Word(33)
        False
        sage: is_Word(Word('abba'))
        True
    """
    return isinstance(obj, AbstractWord)

def is_FiniteWord(obj):
    r"""
    Returns True if obj is a finite word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_FiniteWord
        sage: is_FiniteWord(33)
        False
        sage: is_FiniteWord(Word('baab'))
        True
    """
    if isinstance(obj, AbstractWord) and haslen(obj):
        return True
    else:
        return False

#######################################################################
#                                                                     #
#                    Abstract word classes                            #
#                                                                     #
#######################################################################

class AbstractWord(SageObject):
    def __init__(self, parent, word, mapping=None, format=None, part=slice(None)):
        r"""
        Initialize the AbstractWord object with the parent relationship.

        INPUT:


        -  ``word`` - a function or an iterable, possibly with
           a defined length

        -  ``mapping`` - function, a map sending elements
           output by word to nonnegative integers. The default is
           parent.alphabet().rank; that is, the index of the symbol in the
           alphabet.

        -  ``format`` - string (default None), the explicit
           type of word. Can be either 'empty', 'list', 'string', 'function'
           or 'iterator'. If set to None (the default), an attempt will be
           made to infer the type from properties of word.

        -  ``part`` - slice (default slice(None)), the portion
           of word to use.


        TESTS::

            sage: from sage.combinat.words.word import AbstractWord
            sage: from sage.combinat.words.word_content import BuildWordContent
            sage: W = Words('abc')
            sage: content = BuildWordContent('')
            sage: AbstractWord(W, content)
            Word
            sage: AbstractWord(W, AbstractWord(W, content))
            Word
            sage: AbstractWord(W, '', format='list')
            Word
            sage: AbstractWord(W, lambda n:n, format='function')
            Word
            sage: AbstractWord(W, '')
            Word
            sage: AbstractWord(33, '')
            Traceback (most recent call last):
            ...
            TypeError: the parent must be an instance of Words
            sage: g = AbstractWord(W, '')
            sage: loads(dumps(g))
            Word
        """
        if not is_Words(parent):
            raise TypeError, "the parent must be an instance of Words"
        self._parent = parent

        if isinstance(word, AbstractWord):
            # TODO: Is coercion is necessary here?
            self._word_content = word._word_content
        elif is_WordContent(word):
            # TODO: Is coercion is necessary here?
            self._word_content = word
        else:
            if mapping is None:
                if hasattr(parent, "alphabet"):
                    mapping = parent.alphabet().rank
                else:
                    mapping = id_f
            self._word_content = BuildWordContent(word,
                mapping=mapping, format=format, part=part)

    def size_of_alphabet(self):
        """
        Returns the size of the alphabet of the word.

        TESTS::

            sage: from sage.combinat.words.word import AbstractWord
            sage: AbstractWord(Words('ab'), '').size_of_alphabet()
            2
            sage: AbstractWord(Words('abc'), '').size_of_alphabet()
            3
        """
        return self._parent.size_of_alphabet()

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: from sage.combinat.words.word import AbstractWord
            sage: AbstractWord(Words('ab'), '')._repr_()
            'Word'
        """
        return "Word"

    def parent(self):
        r"""
        Returns the parent object from which the word was created.

        EXAMPLES::

            sage: from sage.combinat.words.word import AbstractWord
            sage: W = Words('ab')
            sage: w = AbstractWord(W, '')
            sage: w.parent()
            Words over Ordered Alphabet ['a', 'b']
            sage: w.parent() is W
            True
        """
        return self._parent

    def coerce(self, other):
        r"""
        Returns a pair of words with a common parent or raises an
        exception.

        This function begins by checking if both words have the same
        parent. If this is the case, then no work is done and both words
        are returned as-is.

        Otherwise it will attempt to convert other to the domain of self.
        If that fails, it will attempt to convert self to the domain of
        other. If both attempts fail, it raises a TypeError to signal
        failure.

        EXAMPLES::

            sage: W1 = Words('abc'); W2 = Words('ab')
            sage: w1 = W1('abc'); w2 = W2('abba'); w3 = W1('baab')
            sage: w1.parent() is w2.parent()
            False
            sage: a, b = w1.coerce(w2)
            sage: a.parent() is b.parent()
            True
        """
        if self._parent != other._parent:
            try:
                other = self.parent()(other)
            except:
                try:
                    self = other.parent()(self)
                except:
                    raise TypeError, "no coercion rule between %r and %r" % (self.alphabet(), other.alphabet())
        return self, other

    def __len__(self):
        r"""
        Returns the length of self.

        TESTS::

            sage: from sage.combinat.words.word import AbstractWord
            sage: w = AbstractWord(Words('ab'), 'abba', mapping='ab'.index)
            sage: len(w)
            4
            sage: w = AbstractWord(Words([0,1]), [0,1,1,0,1,0,0,1])
            sage: len(w)
            8
            sage: w = AbstractWord(Words([0,1]), lambda n:n%2)
            sage: len(w)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        return len(self._word_content)

class AbstractFiniteWord(AbstractWord):
    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: sage.combinat.words.word.AbstractFiniteWord(Words('ab'), sage.combinat.words.word_content.BuildWordContent(''))._repr_()
            'Finite Word'
        """
        return "Finite Word"

class AbstractInfiniteWord(AbstractWord):
    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: sage.combinat.words.word.AbstractInfiniteWord(Words('ab'), sage.combinat.words.word_content.BuildWordContent(''))._repr_()
            'Infinite Word'
        """
        return "Infinite Word"

class Word_over_Alphabet(AbstractWord):
    def __init__(self, parent, *args, **kwds):
        r"""
        Create an Word_over_Alphabet object.

        INPUT:


        -  ``parent`` - a parent object inheriting from
           Words_all that has the alphabet attribute defined

        -  ``*args, **kwds`` - arguments accepted by
           AbstractWord


        EXAMPLES::

            sage: from sage.combinat.words.word import Word_over_Alphabet
            sage: W = Words('abc')
            sage: Word_over_Alphabet(W, "abba")
            Word over Ordered Alphabet ['a', 'b', 'c']

        ::

            sage: from sage.combinat.words.words import Words_all
            sage: W = Words_all()
            sage: W.alphabet = lambda : Alphabet("ab")
            sage: Word_over_Alphabet(W, "abba")
            Word over Ordered Alphabet ['a', 'b']

        TESTS::

            sage: from sage.combinat.words.word import Word_over_Alphabet
            sage: W = Words('abc')
            sage: g = Word_over_Alphabet(W, "abba")
            sage: loads(dumps(g))
            Word over Ordered Alphabet ['a', 'b', 'c']
        """
        if not hasattr(parent, "alphabet"):
            raise TypeError, "parent object has no alphabet attribute"
        super(Word_over_Alphabet, self).__init__(parent, *args, **kwds)

    def alphabet(self):
        r"""
        Returns the alphabet of the parent.

        EXAMPLES::

            sage: Words('abc')('abbabaab').alphabet()
            Ordered Alphabet ['a', 'b', 'c']
        """
        return self._parent.alphabet()

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: sage.combinat.words.word.Word_over_Alphabet(Words([0,1]), sage.combinat.words.word_content.BuildWordContent([]))._repr_()
            'Word over Ordered Alphabet [0, 1]'
        """
        return "Word over %s" % self.alphabet()

    def __iter__(self):
        r"""
        Returns an iterator over the letters of self.

        EXAMPLES::

            sage: W = Words('ab')
            sage: w = W('abba')
            sage: list(w)           # indirect test
            ['a', 'b', 'b', 'a']
        """
        return imap(self.alphabet().unrank, self._word_content)

    def __getitem__(self, key):
        r"""
        TESTS::

            sage: w = Words('012345')('012345')
            sage: w[:]
            word: 012345
            sage: w[3]
            '3'
        """
        res = self._word_content[key]
        if is_WordContent(res):
            return self.parent()(res, format='content')
        else:
            return self.alphabet().unrank(res)

class Word_over_OrderedAlphabet(Word_over_Alphabet):
    def __init__(self, parent, *args, **kwds):
        r"""
        Create an Word_over_OrderedAlphabet object.

        INPUT:


        -  ``parent`` - a parent object inheriting from
           Words_all that has the alphabet attribute defined which returns an
           instance of an OrderedAlphabet_class

        -  ``*args, **kwds`` - arguments accepted by
           AbstractWord


        EXAMPLES::

            sage: from sage.combinat.words.word import Word_over_OrderedAlphabet
            sage: W = Words('abc')
            sage: Word_over_OrderedAlphabet(W, "abba")
            Word over Ordered Alphabet ['a', 'b', 'c']

        ::

            sage: from sage.combinat.words.words import Words_all
            sage: W = Words_all()
            sage: W.alphabet = lambda : Alphabet("ab")
            sage: Word_over_OrderedAlphabet(W, "abba")
            Word over Ordered Alphabet ['a', 'b']

        TESTS::

            sage: from sage.combinat.words.word import Word_over_OrderedAlphabet
            sage: W = Words('abc')
            sage: g = Word_over_OrderedAlphabet(W, "abba")
            sage: loads(dumps(g))
            Word over Ordered Alphabet ['a', 'b', 'c']
        """
        if not hasattr(parent, "alphabet"):
            raise TypeError, "parent object has no alphabet attribute"
        if not isinstance(parent.alphabet(), OrderedAlphabet_class):
            raise TypeError, "underlying alphabet must be an OrderedAlphabet_class"
        super(Word_over_OrderedAlphabet, self).__init__(parent, *args, **kwds)

#######################################################################
#                                                                     #
#                          Infinite words                             #
#                                                                     #
#######################################################################

class InfiniteWord_over_Alphabet(AbstractInfiniteWord, Word_over_Alphabet):
    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: Words([0, 1])(lambda n: n%2)._repr_()
            'Infinite word over [0, 1]'
        """
        return "Infinite word over %s" % self.alphabet().string_rep()

class InfiniteWord_over_OrderedAlphabet(InfiniteWord_over_Alphabet, Word_over_OrderedAlphabet):
    pass

#######################################################################
#                                                                     #
#                           Finite words                              #
#                                                                     #
#######################################################################

class FiniteWord_over_Alphabet(AbstractFiniteWord, Word_over_Alphabet):
    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: from sage.combinat.words.word import FiniteWord_over_Alphabet
            sage: FiniteWord_over_Alphabet(Words("ab"), "abba")._repr_()
            "Finite word of length 4 over Ordered Alphabet ['a', 'b']"
        """
        return "Finite word of length %s over %s" % (len(self), self.alphabet())

class FiniteWord_over_OrderedAlphabet(FiniteWord_over_Alphabet, Word_over_OrderedAlphabet):
    def __init__(self, parent, *args, **kwds):
        r"""
        TESTS::

            sage: sage.combinat.words.word.FiniteWord_over_OrderedAlphabet(Words('abcd'), sage.combinat.words.word_content.BuildWordContent('0123', int))
            word: abcd
        """
        super(FiniteWord_over_OrderedAlphabet, self).__init__(parent, *args, **kwds)
        self._hash = None

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: from sage.combinat.words.word import FiniteWord_over_OrderedAlphabet
            sage: mapping = lambda x : "abc".index(x)
            sage: w = FiniteWord_over_OrderedAlphabet(Words('abc'), 'cabc', mapping=mapping)
            sage: w._repr_()
            'word: cabc'
            sage: Word([0, 1, 0, 1, 1] * 10)._repr_()
            'Finite word of length 50 over [0, 1]'
        """
        global word_options
        if word_options['truncate'] and len(self) > word_options['truncate_length']:
            return "Finite word of length %s over %s" % (len(self), self.alphabet().string_rep())
        else:
            return str(self)

    def __mul__(self, other):
        r"""
        Concatenates two words.

        EXAMPLES::

            sage: W = Words([0, 1])
            sage: W() * W(lambda x:x%2, slice(3))
            word: 010
            sage: Words([0, 1, 2])(lambda x:x%3, slice(3)) * Words('01234')('01234')
            Traceback (most recent call last):
            ...
            TypeError: no coercion rule between Ordered Alphabet [0, 1, 2] and Ordered Alphabet ['0', '1', '2', '3', '4']
        """
        self, other = self.coerce(other)
        c = self._word_content.concatenate(other._word_content)
        return self.parent()(c, format='content')

    __add__ = __mul__

    def __pow__(self, exp):
        r"""
        Return the `exp`-th power of self.

        If `exp` is `\infty`, returns the infinite periodic
        word of base self. Otherwise, `|w|\cdot exp` must be an
        non-negative integer.

        INPUT:


        -  ``exp`` - an integer, a rational, a float number or
           plus infinity.


        OUTPUT:


        -  ``word`` - the exp-th power of self.


        EXAMPLES:

        You can take non-negative integer powers::

            sage: w = Word(range(6)); w
            word: 012345
            sage: w^2
            word: 012345012345
            sage: w^1
            word: 012345
            sage: w^0
            word:
            sage: w^(-1)
            Traceback (most recent call last):
            ...
            ValueError: The exponent must be non-negative.

        You can take non-negative rational powers::

            sage: w = Word(range(6)); w
            word: 012345
            sage: w^(.5)
            word: 012
            sage: w^(1/3)
            word: 01
            sage: (w*w)^(1/2) == w
            True
            sage: w^(5/2)
            word: 012345012345012

        ...but the length of the word times the exponent must be an
        integer::

            sage: w = Word(range(6))
            sage: w^(1/4)
            Traceback (most recent call last):
            ...
            ValueError: Power of the word is not defined on the exponent 1/4: the length of the word (6) times the exponent (1/4) must be a positive integer

        You can take infinite power::

            sage: w = Word(range(6)); w
            word: 012345
            sage: w^oo
            Infinite word over [0, 1, 2, 3, 4, 5]
            sage: (w^oo)[10000000:20000000]
            Finite word of length 10000000 over [0, 1, 2, 3, 4, 5]
            sage: (w^oo)[10000000:10000020]
            word: 45012345012345012345
            sage: Word()^oo
            word:
        """
        if exp is infinity:
            if self.is_empty():
                return self
            else:
                return self.parent()(lambda n: self[n % len(self)])

        nbCaracteres = len(self)*exp

        #If exp*|self| is not an integer
        if int(nbCaracteres) != nbCaracteres :
            raise ValueError, "Power of the word is not defined on the \
exponent %s: the length of the word (%s) times the exponent \
(%s) must be a positive integer"  % (exp, len(self), exp)

        #If exp is negative
        elif exp < 0 :
            raise ValueError, "The exponent must be non-negative."

        return self.parent()(lambda n: self[n % len(self)], part=slice(int(nbCaracteres)) )

    def __str__(self):
        r"""
        Returns the full string representation of the word.

        TESTS::

            sage: Word('abc').__str__()
            'word: abc'
            sage: Words([0, 1])([0, 1, 0, 0, 1] * 10).__str__()
            'word: 01001010010100101001010010100101001010010100101001'
            sage: Words(range(1000))([0,1,10,101]).__str__()
            'word: 0,1,10,101'
        """
        global word_options
        return word_options['identifier'] + self.string_rep()

    def string_rep(self):
        r"""
        Returns the raw sequence of letters as a string.

        EXAMPLES::

            sage: Words('ab')('abbabaab').string_rep()
            'abbabaab'
            sage: Words(range(2))([0, 1, 0, 0, 1]).string_rep()
            '01001'
            sage: Words(range(1000))([0,1,10,101]).string_rep()
            '0,1,10,101'
            sage: WordOptions(letter_separator='-')
            sage: Words(range(1000))([0,1,10,101]).string_rep()
            '0-1-10-101'
        """
        global word_options
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = map(str, list(self))
            if all(len(a)==1 for a in letters):
                return ''.join(letters)
            else:
                return ls.join(letters)
        elif word_options['display'] == 'list':
            return str(list(self))

    def __cmp__(self, other):
        """
        Compares two words lexicographically according to the order defined
        in the common alphabet. Provides for all normal comparison
        operators.

        EXAMPLES::

            sage: W1 = Words('12'); W2 = Words('123')
            sage: W2('123').__cmp__(W1('1211')) > 0
            True
            sage: W1('2111').__cmp__(W2('12')) > 0
            True
            sage: W2('123').__cmp__(W2('123')) == 0
            True
            sage: W1('121').__cmp__(W2('121')) == 0
            True
            sage: W2('123').__cmp__(W1('22')) < 0
            True
            sage: W1('122').__cmp__(W2('32')) < 0
            True
            sage: W1('122').__cmp__(Words('ab')('abba'))  # todo: not implemented
            Traceback (most recent call last):
            ...
            TypeError: no coercion rule between Alphabet: ['1', '2'] and Alphabet: ['a', 'b']
        """
        if not is_Word(other):
            return NotImplemented
        try:
            self, other = self.coerce(other)
        except TypeError:
            return NotImplemented
        for (c1, c2) in izip(self._word_content, other._word_content):
            r = c1 - c2
            if r != 0: return r
        return len(self) - len(other)

    def __hash__(self):
        r"""
        Returns the hash for this word.

        TESTS::

            sage: h = hash(Word('abc'))    # indirect test
            sage: Word('abc').__hash__() == Word('abc').__hash__()
            True
        """
        if self._hash is None:
            res = 5381
            for s in self._word_content:
                res = ((res << 5) + res) + s
            self._hash = res
        return self._hash

    def is_empty(self):
        """
        Returns True if the length of self is zero, and False otherwise.

        EXAMPLES::

            sage: W=Words('ab')
            sage: W().is_empty()
            True
            sage: W('a').is_empty()
            False
        """
        return len(self)==0

    def is_prefix_of(self, other):
        """
        Returns True if self is a prefix of other, and False otherwise.

        EXAMPLES::

            sage: V = Words('0123456789')
            sage: w = V('0123456789')
            sage: y = V('012345')
            sage: y.is_prefix_of(w)
            True
            sage: w.is_prefix_of(y)
            False
            sage: W = Words('ab')
            sage: w.is_prefix_of(W())
            False
            sage: W().is_prefix_of(w)
            True
            sage: W().is_prefix_of(W())
            True
        """
        return self == other[:len(self)]

    def is_proper_prefix_of(self, other):
        """
        Returns True if self is a proper prefix of other, and False
        otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('12').is_proper_prefix_of(W('123'))
            True
            sage: W('12').is_proper_prefix_of(W('12'))
            False
            sage: W().is_proper_prefix_of(W('123'))
            True
            sage: W('123').is_proper_prefix_of(W('12'))
            False
            sage: W().is_proper_prefix_of(W())
            False
        """
        return self.is_prefix_of(other) and len(self) < len(other)

    def is_suffix_of(self, other):
        """
        Returns True if w is a suffix of other, and False otherwise.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: w = W('0123456789')
            sage: y = W('56789')
            sage: y.is_suffix_of(w)
            True
            sage: w.is_suffix_of(y)
            False
            sage: W('579').is_suffix_of(w)
            False
            sage: W().is_suffix_of(y)
            True
            sage: w.is_suffix_of(W())
            False
            sage: W().is_suffix_of(W())
            True
        """
        return self.is_empty() or self == other[-len(self):]

    def is_proper_suffix_of(self, other):
        """
        Returns True if self is a proper suffix of other, and False
        otherwise.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('23').is_proper_suffix_of(W('123'))
            True
            sage: W('12').is_proper_suffix_of(W('12'))
            False
            sage: W().is_proper_suffix_of(W('123'))
            True
            sage: W('123').is_proper_suffix_of(W('12'))
            False
        """
        return self.is_suffix_of(other) and len(self) < len(other)

    def reversal(self):
        """
        Returns the reversal of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('124563').reversal()
            word: 365421
        """
        return self[::-1]

    def is_palindrome(self, f=None):
        r"""
        Returns True if self is a palindrome (or a `f` -palindrome),
        and False otherwise.

        - In French: Soit `f : \Sigma \rightarrow \Sigma` une
          involution qui s'étend évidemment à un morphisme sur
          `\Sigma^*`. On dit que `w\in\Sigma^*` est un
          *`f`-pseudo-palindrome* [1], ou plus simplement un
          *`f`-palindrome*, si `w=f(\tilde{w})` (extrait
          de [2]).

        - In English Let `f : \Sigma \rightarrow \Sigma` be an
          involution that extends to a morphism on
          `\Sigma^*`. We say that `w\in\Sigma^*` is a
          *`f`-palindrome* if `w=f(\tilde{w})` [2].

        Also called *`f` -pseudo-palindrome* [1].

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``boolean`` - if f is None, whether self is a
           palindrome; otherwise, whether self is a f-palindrome.


        EXAMPLES: Some palindromes...

        ::

            sage: W=Words('abcdefghijklmnopqrstuvwxyz I')
            sage: W('esope reste ici et se repose').is_palindrome()
            False
            sage: W('esoperesteicietserepose').is_palindrome()
            True
            sage: W('I saw I was I').is_palindrome()
            True
            sage: Word('abbcbba').is_palindrome()
            True
            sage: Word('abcbdba').is_palindrome()
            False

        Some `f`-palindromes... You can use a str::

            sage: Word('aababb').is_palindrome('a->b,b->a')
            True
            sage: Word('abacbacbab').is_palindrome('a->b,b->a,c->c')
            True

        You can use WordMorphism::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('aababb').is_palindrome(f)
            True

        You can use a dictionary::

            sage: f = {'a':'b','b':'a'}
            sage: Word('aababb').is_palindrome(f)
            True
            sage: w = words.ThueMorseWord()[:8]; w
            word: 01101001
            sage: w.is_palindrome(f={0:[1],1:[0]})
            True

        self must be in the domain of the involution::

            sage: f = WordMorphism('a->a')
            sage: Word('aababb').is_palindrome(f)
            Traceback (most recent call last):
            ...
            ValueError: self must be in the domain of the given involution

        The given involution must be an involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('abab').is_palindrome(f)
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution

        TESTS::

            sage: Y = Words('ab')
            sage: Y().is_palindrome()
            True
            sage: Y('a').is_palindrome()
            True
            sage: Y('ab').is_palindrome()
            False
            sage: Y('aba').is_palindrome()
            True
            sage: Y('aa').is_palindrome()
            True
            sage: E = 'a->b,b->a'
            sage: Y().is_palindrome(E)
            True
            sage: Y('a').is_palindrome(E)
            False
            sage: Y('ab').is_palindrome(E)
            True
            sage: Y('aa').is_palindrome(E)
            False
            sage: Y('aba').is_palindrome(E)
            False
            sage: Y('abab').is_palindrome(E)
            True

        REFERENCES:

        - [1] V. Anne, L.Q. Zamboni, I. Zorca, Palindromes and
          Pseudo- Palindromes in Episturmian and Pseudo-Palindromic Infinite
          Words, in : S. Brlek, C. Reutenauer (Eds.), Words 2005,
          Publications du LaCIM, Vol. 36 (2005) 91-100.

        - [2] S. Labbé, Propriétés combinatoires des `f`-palindromes, Mémoire de
          maîtrise en Mathématiques, Montréal, UQAM, 2008, 109 pages.
        """
        if f is None:
            l = len(self)
            return self[:l/2] == self[l/2 + l%2:].reversal()

        from sage.combinat.words.morphism import WordMorphism
        f = WordMorphism(f)

        if self not in f.domain():
            raise ValueError, "self must be in the domain of "\
                                 +"the given involution"

        if not f.is_involution():
            raise ValueError, "f must be an involution"

        l = len(self)
        return self[:l/2 + l%2] == f(self[l/2:].reversal())

    def lps(self, f=None):
        r"""
        Returns the longest palindromic (or `f`-palindromic) suffix
        of self.

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``word`` - If f is None, the longest palindromic
           suffix of self; otherwise, the longest f-palindromic suffix of
           self.


        EXAMPLES::

            sage: Word('0111').lps()
            word: 111
            sage: Word('011101').lps()
            word: 101
            sage: Word('6667').lps()
            word: 7
            sage: Word('abbabaab').lps()
            word: baab
            sage: Word().lps()
            word:
            sage: Word('abbabaab').lps('a->b,b->a')
            word: abbabaab
        """
        for i in range(len(self)+1):
            fact = self[i:]
            if fact.is_palindrome(f=f):
                return fact

    def _lps(self, l=None, f=None):
        r"""
        Returns the longest palindromic (or `f`-palindromic) suffix
        of self.

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).

        -  ``l`` - integer (default: None) the length of the
           longest palindrome suffix of self[:-1]


        OUTPUT:


        -  ``word`` - If f is None, the longest palindromic
           suffix of self; otherwise, the longest f-palindromic suffix of
           self.


        EXAMPLES::

            sage: W = Words('1234')
            sage: w = W('33412321')
            sage: w._lps(3)
            word: 12321
            sage: Y = Words('01')
            sage: w = Y('01101001')
            sage: w._lps(l=2)
            word: 1001
            sage: w._lps()
            word: 1001
            sage: w._lps(None)
            word: 1001
            sage: Y()._lps(2)
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
            sage: v=Words('ab')('abbabaab')
            sage: pal=v[:0]
            sage: for i in range(1,len(v)+1):
            ...     pal=v[:i]._lps(len(pal))
            ...     print pal
            ...
            word: a
            word: b
            word: bb
            word: abba
            word: bab
            word: aba
            word: aa
            word: baab
            sage: v=Words('ab')('abbabaab')
            sage: pal=v[:0]
            sage: for i in range(1,len(v)+1):
            ...     pal=v[:i]._lps(len(pal),'a->b,b->a')
            ...     print pal
            ...
            word:
            word: ab
            word:
            word: ba
            word: ab
            word: baba
            word: bbabaa
            word: abbabaab
        """
        if l == None:
            return self.lps(f=f)

        #If l == len(w[:-1]), there is no shortcut
        if len(self) == l + 1:
            return self.lps(f=f)

        #Obtain the letter to the left (g) and to the right (d) of the
        #precedent lps of self
        g = self[-l-2]
        d = self[-1]

        #If the word g*d is a $f$-palindrome, the result follows
        if self.parent()([g, d]).is_palindrome(f=f):
            return self[-l-2:]

        #Otherwise, the length of the lps of self is smallest than l+2
        else:
            return self[-l-1:].lps(f=f)

    def palindromic_lacunas_study(self, f=None):
        r"""
        Returns interesting statistics about longest
        (`f`-)palindromic suffixes and lacunas of self (see [1] and
        [2]).

        Note that a word `w` has at most `|w| + 1`
        different palindromic factors (see [4]).

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understand (dict, str, ...).


        OUTPUT:


        -  ``list`` - list of the length of the longest
           palindromic suffix (lps) for each non-empty prefix of self;

        -  ``list`` - list of all the lacunas, i.e. positions
           where there is no unioccurrent lps;

        -  ``set`` - set of palindromic factors of self.


        EXAMPLES::

            sage: W=Words('ab')
            sage: a,b,c = W('abbabaabbaab').palindromic_lacunas_study()
            sage: a
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]
            sage: b
            [8, 9]
            sage: c          # random order
            set([word: , word: b, word: bab, word: abba, word: bb, word: aa, word: baabbaab, word: baab, word: aba, word: aabbaa, word: a])
            sage: a,b,c = W('abbabaab').palindromic_lacunas_study(f='a->b,b->a')
            sage: a
            [0, 2, 0, 2, 2, 4, 6, 8]
            sage: b
            [0, 2, 4]
            sage: c           # random order
            set([word: , word: ba, word: baba, word: ab, word: bbabaa, word: abbabaab])
            sage: c == set([W(), W('ba'), W('baba'), W('ab'), W('bbabaa'), W('abbabaab')])
            True

        REFERENCES:

        - [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic
          lacunas of the Thue-Morse word, Proc. GASCOM 2008 (June
          16-20 2008, Bibbiena, Arezzo-Italia), 53-67.

        - [2] A. Blondin-Massé, S. Brlek, A.  Frosini, S. Labbé,
          S. Rinaldi, Reconstructing words from a fixed palindromic
          length sequence, Proc. TCS 2008, 5th IFIP International
          Conference on Theoretical Computer Science (September 8-10
          2008, Milano, Italia), accepted.

        - [3] S. Labbé, Propriétés combinatoires des
          `f`-palindromes, Mémoire de maitrise en Mathématiques,
          Montréal, UQAM, 2008, 109 pages.

        - [4] X. Droubay, J. Justin, G.  Pirillo, Episturmian words
          and some constructions of de Luca and Rauzy,
          Theoret. Comput. Sci. 255 (2001) 539-553.
        """
        #Initialize the results of computations
        palindromes = set()
        lengths_lps = [None] * len(self)
        lacunas = []

        #Initialize the first lps
        pal = self.parent()()
        palindromes.add(pal)

        #For all the non-empty prefixes of self,
        for i in xrange(len(self)):

            #Compute its longest $f$-palindromic suffix using the preceding lps (pal)
            pal = self[:i+1]._lps(l=len(pal),f=f)

            lengths_lps[i] = len(pal)

            if pal in palindromes:
                lacunas.append(i)
            else :
                palindromes.add(pal)

        return lengths_lps, lacunas, palindromes

    def lengths_lps(self, f=None):
        r"""
        Returns the list of the length of the longest palindromic suffix
        (lps) for each non-empty prefix of self.

        It corresponds to the function `G_w` defined in [2].

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understand (dict, str, ...).


        OUTPUT:


        -  ``list`` - list of the length of the longest
           palindromic suffix (lps) for each non-empty prefix of self.


        EXAMPLES::

            sage: Word().lengths_lps()
            []
            sage: Word('a').lengths_lps()
            [1]
            sage: Word('aaa').lengths_lps()
            [1, 2, 3]
            sage: Word('abbabaabbaab').lengths_lps()
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]
            sage: Word('abbabaabbaab').lengths_lps(f='a->b,b->a')
            [0, 2, 0, 2, 2, 4, 6, 8, 4, 6, 4, 6]
            sage: Word([5,8,5,5,8,8,5,5,8,8,5,8,5]).lengths_lps(f={5:[8],8:[5]})
            [0, 2, 2, 0, 2, 4, 6, 4, 6, 8, 10, 12, 4]

        REFERENCES:

        - [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic
          lacunas of the Thue-Morse word, Proc. GASCOM 2008 (June
          16-20 2008, Bibbiena, Arezzo-Italia), 53-67.

        - [2] A. Blondin-Massé, S. Brlek, A.  Frosini, S. Labbé,
          S. Rinaldi, Reconstructing words from a fixed palindromic
          length sequence, Proc. TCS 2008, 5th IFIP International
          Conference on Theoretical Computer Science (September 8-10
          2008, Milano, Italia), accepted.
        """
        return self.palindromic_lacunas_study(f=f)[0]

    def lacunas(self, f=None):
        r"""
        Returns the list of all the lacunas of self.

        A *lacuna* is a position in a word where the longest palindromic
        suffix is not unioccurrent (see [1]).

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``list`` - list of all the lacunas of self.


        EXAMPLES::

            sage: words.ThueMorseWord()[:100].lacunas()
            [8, 9, 24, 25, 32, 33, 34, 35, 36, 37, 38, 39, 96, 97, 98, 99]
            sage: words.ThueMorseWord()[:50].lacunas(f={0:[1],1:[0]})
            [0, 2, 4, 12, 16, 17, 18, 19, 48, 49]

        REFERENCES:

        - [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic
          lacunas of the Thue-Morse word, Proc. GASCOM 2008 (June
          16-20 2008, Bibbiena, Arezzo-Italia), 53-67.
        """
        return self.palindromic_lacunas_study(f=f)[1]

    def lengths_unioccurrent_lps(self, f=None):
        r"""
        Returns the list of the lengths of the unioccurrent longest
        palindromic suffixes (lps) for each non-empty prefix of self. No
        unioccurrent lps are indicated by None.

        It corresponds to the function `H_w` defined in [1] and
        [2].

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understand (dict, str, ...).


        OUTPUT:


        -  ``list`` - list of the length of the unioccurrent
           longest palindromic suffix (lps) for each non-empty prefix of
           self. No unioccurrent lps are indicated by None.


        EXAMPLES::

            sage: f = words.FibonacciWord()[:20]
            sage: f.lengths_unioccurrent_lps() == f.lengths_lps()
            True
            sage: words.ThueMorseWord()[:20].lengths_unioccurrent_lps()
            [1, 1, 2, 4, 3, 3, 2, 4, None, None, 6, 8, 10, 12, 14, 16, 6, 8, 10, 12]
            sage: words.ThueMorseWord()[:15].lengths_unioccurrent_lps(f={1:[0],0:[1]})
            [None, 2, None, 2, None, 4, 6, 8, 4, 6, 4, 6, None, 4, 6]

        REFERENCES:

        - [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic
          lacunas of the Thue-Morse word, Proc. GASCOM 2008 (June
          16-20 2008, Bibbiena, Arezzo-Italia), 53-67.

        - [2] A. Blondin-Massé, S. Brlek, A.  Frosini, S. Labbé,
          S. Rinaldi, Reconstructing words from a fixed palindromic
          length sequence, Proc. TCS 2008, 5th IFIP International
          Conference on Theoretical Computer Science (September 8-10
          2008, Milano, Italia), accepted.
        """
        l = self.lengths_lps(f=f)
        for i in self.lacunas(f=f):
            l[i] = None
        return l

    def palindromes(self, f=None):
        r"""
        Returns the set of all palindromic (or `f`-palindromic)
        factors of self.

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``set`` - If f is None, the set of all palindromic
           factors of self; otherwise, the set of all f-palindromic factors of
           self.


        EXAMPLES::

            sage: W = Words('01')
            sage: sorted(W('01101001').palindromes())
            [word: , word: 0, word: 00, word: 010, word: 0110, word: 1, word: 1001, word: 101, word: 11]
            sage: sorted(W('00000').palindromes())
            [word: , word: 0, word: 00, word: 000, word: 0000, word: 00000]
            sage: sorted(W('0').palindromes())
            [word: , word: 0]
            sage: sorted(W('').palindromes())
            [word: ]
            sage: sorted(W().palindromes())
            [word: ]
            sage: sorted(Word('abbabaab').palindromes('a->b,b->a'))
            [word: , word: ab, word: abbabaab, word: ba, word: baba, word: bbabaa]
        """
        return self.palindromic_lacunas_study(f=f)[2]

    def defect(self, f=None):
        r"""
        Returns the defect of self.

        The *defect* of a finite word `w` is given by
        `D(w)=|w|+1-|PAL(w)|`, where `PAL(w)` denotes the
        set of palindromic factors of `w` (including the empty
        word). See [1].

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``integer`` - If f is None, the palindromic defect
           of self; otherwise, the f-palindromic defect of self.


        EXAMPLES::

            sage: words.ThueMorseWord()[:100].defect()
            16
            sage: words.FibonacciWord()[:100].defect()
            0
            sage: W = Words('01')
            sage: W('000000000000').defect()
            0
            sage: W('011010011001').defect()
            2
            sage: W('0101001010001').defect()
            0
            sage: W().defect()
            0
            sage: Word('abbabaabbaababba').defect()
            2
            sage: Word('abbabaabbaababba').defect('a->b,b->a')
            4

        REFERENCES:

        - [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the
          Palindromic Complexity of Infinite Words, in J. Berstel, J.
          Karhumaki, D. Perrin, Eds, Combinatorics on Words with
          Applications, International Journal of Foundation of
          Computer Science, Vol. 15, No. 2 (2004) 293-306.
        """
        return len(self)+1-len(self.palindromes(f=f))

    def is_full(self, f=None):
        r"""
        Returns True if self has defect 0, and False otherwise.

        A word is *full* if its defect is zero (see [1]).

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``boolean`` - If f is None, whether self is full;
           otherwise, whether self is full of `f`-palindromes.


        EXAMPLES::

            sage: words.ThueMorseWord()[:100].is_full()
            False
            sage: words.FibonacciWord()[:100].is_full()
            True
            sage: Words('0')('000000000000000').is_full()
            True
            sage: Words('01')('011010011001').is_full()
            False
            sage: Words('123456789')('2194').is_full()
            True
            sage: Word().is_full()
            True
            sage: Word().is_full('a->b,b->a')
            True
            sage: w = Word('ab')
            sage: w.is_full()
            True
            sage: w.is_full('a->b,b->a')
            False

        REFERENCES:

        - [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the
          Palindromic Complexity of Infinite Words, in J. Berstel, J.
          Karhumaki, D. Perrin, Eds, Combinatorics on Words with
          Applications, International Journal of Foundation of
          Computer Science, Vol. 15, No. 2 (2004) 293-306.
        """
        return self.defect(f=f) == 0

    def palindromic_closure(self, side='right', f=None):
        r"""
        Returns the shortest palindrome having self as a prefix (or as a
        suffix if side=='left').

        Retourne le plus petit palindrome ayant self comme prefixe (ou
        comme suffixe si side=='left').

        INPUT:


        -  ``side`` - 'right' or 'left' (default: 'right') the
           direction of the closure

        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``word`` - If f is None, the right palindromic
           closure of self; otherwise, the right f-palindromic closure of
           self. If side is 'left', the left palindromic closure.


        EXAMPLES::

            sage: W = Words('1234567890')
            sage: W('1233').palindromic_closure()
            word: 123321
            sage: W('12332').palindromic_closure()
            word: 123321
            sage: W('0110343').palindromic_closure()
            word: 01103430110
            sage: W('0110343').palindromic_closure(side='left')
            word: 3430110343
            sage: W('01105678').palindromic_closure(side='left')
            word: 876501105678
            sage: w = Word('abbaba')
            sage: w.palindromic_closure()
            word: abbababba
            sage: w.palindromic_closure(f='a->b,b->a')
            word: abbabaab
            sage: w.palindromic_closure(f='a->b,b->a',side='left')
            word: babaabbaba
            sage: w.palindromic_closure(f='a->b,b->b',side='left')
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution
            sage: w.palindromic_closure(f='a->c,c->a',side='left')
            Traceback (most recent call last):
            ...
            ValueError: self must be in the domain of the given involution

        REFERENCES:

        - [1] A. de Luca, A. De Luca, Pseudopalindrome closure
          operators in free monoids, Theoret. Comput. Sci. 362 (2006)
          282-300.
        """
        if f is None:
            if side == 'right':
                l = len(self.lps())
                return self * self[-(l+1)::-1]
            elif side == 'left':
                l = len(self.reversal().lps())
                return self[:l-1:-1] * self
            else:
                raise ValueError, "side must be either 'left' or 'right' (not %s) " % side
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)

            if self not in f.domain():
                raise ValueError, "self must be in the domain of "\
                                     +"the given involution"

            if not f.is_involution():
                raise ValueError, "f must be an involution"

            if side == 'right':
                l = len(self.lps(f=f))
                return self * f(self[-(l+1)::-1])
            elif side == 'left':
                l = len(self.reversal().lps(f=f))
                return f(self[:l-1:-1]) * self
            else:
                raise ValueError, "side must be either 'left' or 'right' (not %s) " % side

    def iterated_palindromic_closure(self, side='right', f=None):
        r"""
        Returns the iterated (`f`-)palindromic closure of self.

        INPUT:


        -  ``side`` - 'right' or 'left' (default: 'right') the
           direction of the closure

        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``word`` - If f is None, the right iterated
           palindromic closure of self; otherwise, the right iterated
           f-palindromic closure of self. If side is 'left', the left
           palindromic closure.


        EXAMPLES::

            sage: W = Words('123')
            sage: W('123').iterated_palindromic_closure()
            word: 1213121
            sage: W('123').iterated_palindromic_closure(side='left')
            word: 3231323
            sage: W('1').iterated_palindromic_closure()
            word: 1
            sage: W().iterated_palindromic_closure()
            word:
            sage: W=Words('ab')
            sage: W('ab').iterated_palindromic_closure(f='a->b,b->a')
            word: abbaab
            sage: W('ab').iterated_palindromic_closure(f='a->b,b->a',side='left')
            word: abbaab
            sage: W('aab').iterated_palindromic_closure(f='a->b,b->a')
            word: ababbaabab
            sage: W('aab').iterated_palindromic_closure(f='a->b,b->a',side='left')
            word: abbaabbaab

        TESTS::

            sage: W('aab').iterated_palindromic_closure(f='a->b,b->a',side='leftt')
            Traceback (most recent call last):
            ...
            ValueError: side must be either 'left' or 'right' (not leftt)
            sage: W('aab').iterated_palindromic_closure(f='a->b,b->b',side='left')
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution

        REFERENCES:

        - [1] A. de Luca, A. De Luca, Pseudopalindrome closure
          operators in free monoids, Theoret. Comput. Sci. 362 (2006)
          282-300.
        """
        if side == 'right':
            par = self.parent()
            return reduce(lambda r, s: (r*par([s])).palindromic_closure(f=f),
                            self, par())
        elif side == 'left':
            if f is None:
                return self.reversal().iterated_palindromic_closure(f=f)
            else:
                from sage.combinat.words.morphism import WordMorphism
                f = WordMorphism(f)

                if self not in f.domain():
                    raise ValueError, "self must be in the domain of "\
                                         +"the given involution"

                if not f.is_involution():
                    raise ValueError, "f must be an involution"

                return f(self).reversal().iterated_palindromic_closure(f=f)

        else:
            raise ValueError, "side must be either 'left' or 'right' (not %s) " % side

    def is_symmetric(self, f=None):
        r"""
        Returns True if self is symmetric (or `f`-symmetric), and
        False otherwise.

        A word is *symmetric* (resp. *`f`-symmetric*) if it is the
        product of two palindromes (resp. `f`-palindromes). See [1]
        and [2].

        INPUT:


        -  ``f`` - involution (default: None) on the alphabet
           of self. It must be something that WordMorphism's constructor
           understands (dict, str, ...).


        OUTPUT:


        -  ``boolean`` - If f is None, whether self is
           symmetric; otherwise, whether self is `f`-symmetric.


        EXAMPLES::

            sage: W = Words('ab')
            sage: W('abbabab').is_symmetric()
            True
            sage: W('ababa').is_symmetric()
            True
            sage: W('aababaabba').is_symmetric()
            False
            sage: W('aabbbaababba').is_symmetric()
            False
            sage: W('aabbbaababba').is_symmetric('a->b,b->a')
            True

        REFERENCES:

        - [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the
          Palindromic Complexity of Infinite Words, in J. Berstel, J.
          Karhumaki, D. Perrin, Eds, Combinatorics on Words with
          Applications, International Journal of Foundation of
          Computer Science, Vol. 15, No. 2 (2004) 293-306.

        - [2] A. de Luca, A. De Luca, Pseudopalindrome closure
          operators in free monoids, Theoret.  Comput. Sci. 362 (2006)
          282-300.
        """
        for i in range(len(self)):
            if self[:i].is_palindrome(f=f) and self[i:].is_palindrome(f=f):
                return True
        return False

    def is_square(self):
        """
        Returns True if self is a square, and False otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('1212').is_square()
            True
            sage: W('1213').is_square()
            False
            sage: W('12123').is_square()
            False
            sage: W().is_square()
            True
        """
        if len(self) % 2 != 0:
            return False
        else:
            l = len(self) / 2
            return self[:l] == self[l:]

    def is_square_free(self):
        """
        Returns True if self does not contain squares, and False
        otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('12312').is_square_free()
            True
            sage: W('31212').is_square_free()
            False
            sage: W().is_square_free()
            True
        """
        l = len(self)
        if l < 2:
            return True
        suff = self
        for i in xrange(0, l - 2):
            for ll in xrange(2, l-i+1, 2):
                if suff[:ll].is_square():
                    return False
            suff = suff[1:]
        return True

    def is_cube(self):
        """
        Returns True if self is a cube, and False otherwise.

        EXAMPLES::

            sage: W = Words('012')
            sage: W('012012012').is_cube()
            True
            sage: W('01010101').is_cube()
            False
            sage: W().is_cube()
            True
            sage: W('012012').is_cube()
            False
        """
        if len(self) % 3 != 0:
            return False
        l = len(self) / 3
        return self[:l] == self[l:2*l] == self[2*l:]

    def is_cube_free(self):
        """
        Returns True if self does not contain cubes, and False otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('12312').is_cube_free()
            True
            sage: W('32221').is_cube_free()
            False
            sage: W().is_cube_free()
            True
        """
        l = len(self)
        if l < 3:
            return True
        suff = self
        for i in xrange(0, l - 3):
            for ll in xrange(3, l-i+1, 3):
                if suff[:ll].is_cube():
                    return False
            suff = suff[1:]
        return True

    def is_quasiperiodic(self):
        """
        Returns True if self is quasiperiodic, and False otherwise.

        A finite or infinite word `w` is *quasiperiodic* if it can
        be constructed by concatenations and superpositions of one of its
        proper factors `u`, which is called a *quasiperiod* of
        `w`. See for instance [1], [2], and [3].

        EXAMPLES::

            sage: W = Words('abc')
            sage: W('abaababaabaababaaba').is_quasiperiodic()
            True
            sage: W('abacaba').is_quasiperiodic()
            False
            sage: W('a').is_quasiperiodic()
            False
            sage: W().is_quasiperiodic()
            False
            sage: W('abaaba').is_quasiperiodic()
            True

        REFERENCES:

        - [1] A. Apostolico, A. Ehrenfeucht, Efficient detection of
          quasiperiodicities in strings, Theoret. Comput. Sci. 119
          (1993) 247-265.

        - [2] S. Marcus, Quasiperiodic infinite words, Bull. Eur.
          Assoc. Theor. Comput. Sci. 82 (2004) 170-174. [3] A. Glen,
          F. Levé, G. Richomme, Quasiperiodic and Lyndon episturmian
          words, Preprint, 2008, arXiv:0805.0730.
        """
        l = len(self)
        if l <= 1:
           return False
        for i in range(1, l - 1):
            return_lengths = [len(x) for x in self.return_words(self[:i])]
            if return_lengths != []:
               if (max(return_lengths) <= i and self[l-i:l] == self[:i]):
                  return True
        return False

    def _quasiperiods_list(self):
        """
        Returns the quasiperiods of self as a list ordered from shortest to
        longest.

        EXAMPLES::

            sage: W = Words('abc')
            sage: W('abaababaabaababaaba')._quasiperiods_list()
            [word: aba, word: abaaba, word: abaababaaba]
            sage: W('abaaba')._quasiperiods_list()
            [word: aba]
            sage: W('abacaba')._quasiperiods_list()
            []
        """
        l = len(self)
        if l <= 1:
           return []
        Q = []
        for i in range(1, l - 1):
            return_lengths = [len(x) for x in self.return_words(self[:i])]
            if return_lengths != []:
               if (max(return_lengths) <= i and self[l-i:l] == self[:i]):
                  Q.append(self[:i])
        return Q

    def quasiperiods(self):
        """
        Returns the set of quasiperiods of self.

        Let `w` be a finite or infinite word. A *quasiperiod* of
        `w` is a proper factor `u` of `w` such that
        the occurrences of `u` in `w` entirely cover
        `w`, i.e., every position of `w` falls within some
        occurrence of `u` in `w`. See for instance [1],
        [2], and [3].

        TESTS::

            sage: W = Words('abc')
            sage: W('abaababaabaababaaba').quasiperiods() == set([W('aba'), W('abaaba'), W('abaababaaba')])
            True
            sage: W('abaaba').quasiperiods() == set([W('aba')])
            True
            sage: W('abacaba').quasiperiods() == set([])
            True

        REFERENCES:

        - [1] A. Apostolico, A. Ehrenfeucht, Efficient detection of
          quasiperiodicities in strings, Theoret. Comput. Sci. 119
          (1993) 247-265.

        - [2] S. Marcus, Quasiperiodic infinite words, Bull. Eur.
          Assoc. Theor. Comput. Sci. 82 (2004) 170-174. [3] A. Glen,
          F. Levé, G. Richomme, Quasiperiodic and Lyndon episturmian
          words, Preprint, 2008, arXiv:0805.0730.
        """
        return set(self._quasiperiods_list())

    def commutes_with(self, other):
        """
        Returns True if self commutes with other, and False otherwise.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12').commutes_with(W('12'))
            True
            sage: W('12').commutes_with(W('11'))
            False
            sage: W().commutes_with(W('21'))
            True
        """
        return (self * other) == (other * self)

    def conjugate(self, pos):
        """
        Returns the conjugate at pos of self.

        pos can be any integer, the distance used is the modulo by the
        length of self.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12112').conjugate(1)
            word: 21121
            sage: W().conjugate(2)
            word:
            sage: W('12112').conjugate(8)
            word: 12121
            sage: W('12112').conjugate(-1)
            word: 21211
        """
        if self.is_empty():
            return self
        pos_mod = pos % len(self)
        return self[pos_mod:] * self[:pos_mod]

    def _conjugates_list(self):
        r"""
        Returns the list of conjugates of self, ordered from the 0-th to
        the (L-1)-st conjugate, where L is the length of self.

        TESTS::

            sage: W = Words('abc')
            sage: W('cbbca')._conjugates_list()
            [word: cbbca, word: bbcac, word: bcacb, word: cacbb, word: acbbc]
            sage: W('abcabc')._conjugates_list()
            [word: abcabc,
             word: bcabca,
             word: cabcab,
             word: abcabc,
             word: bcabca,
             word: cabcab]
            sage: W()._conjugates_list()
            [word: ]
            sage: W('a')._conjugates_list()
            [word: a]
        """
        S = [self]
        for i in range(1,len(self)):
            S.append(self.conjugate(i))
        return S

    def conjugates(self):
        """
        Returns the set of conjugates of self.

        TESTS::

            sage: W = Words('abc')
            sage: W('cbbca').conjugates() == set([W('cacbb'),W('bbcac'),W('acbbc'),W('cbbca'),W('bcacb')])
            True
            sage: W('abcabc').conjugates() == set([W('abcabc'),W('bcabca'),W('cabcab')])
            True
            sage: W().conjugates() == set([W()])
            True
            sage: W('a').conjugates() == set([W('a')])
            True
        """
        S = set([self])
        for i in range(1,self.primitive_length()):
            S.add(self.conjugate(i))
        return S

    def conjugate_position(self, other):
        """
        Returns the position where self is conjugate with other. Returns
        None if there is no such position.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('12113').conjugate_position(W('31211'))
            1
            sage: W('12131').conjugate_position(W('12113')) is None
            True
            sage: W().conjugate_position(W('123')) is None
            True
        """
        if len(self) != len(other):
            return None
        self, conj = self.coerce(other)
        if self == conj:
            return 0
        for l in xrange(1, len(conj) - 1):
            conj = conj.conjugate(1)
            if self == conj:
                return l
        return None

    def is_conjugate_with(self, other):
        """
        Returns True if self is a conjugate of other, and False otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('11213').is_conjugate_with(W('31121'))
            True
            sage: W().is_conjugate_with(W('123'))
            False
            sage: W('112131').is_conjugate_with(W('11213'))
            False
            sage: W('12131').is_conjugate_with(W('11213'))
            True
        """
        return self.conjugate_position(other) is not None

    def is_cadence(self, seq):
        """
        Returns True if seq is a cadence of self, and False otherwise.

        A *cadence* is an increasing sequence of indexes that all map to
        the same letter.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('121132123').is_cadence([0, 2, 6])
            True
            sage: W('121132123').is_cadence([0, 1, 2])
            False
            sage: W('121132123').is_cadence([])
            True
        """
        if len(seq) == 0:
            return True
        try:
            it = iter(self)
            s = islice(it, seq[0], None).next()
            for i in xrange(1, len(seq)):
                steps = seq[i] - seq[i-1]
                for n in xrange(steps-1): it.next()
                if it.next() != s:
                    return False
        except StopIteration:
            return False
        return True

    def longest_common_prefix(self, other):
        """
        Returns the longest common prefix of self and other.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: w=W('12345')
            sage: y=W('1236777')
            sage: w.longest_common_prefix(y)
            word: 123
            sage: w.longest_common_prefix(w)
            word: 12345
            sage: y.longest_common_prefix(w)
            word: 123
            sage: y.longest_common_prefix(w)==w.longest_common_prefix(y)
            True
            sage: y.longest_common_prefix(y)==y
            True
            sage: W().longest_common_prefix(w)
            word:
            sage: w.longest_common_prefix(W()) == w
            False
            sage: w.longest_common_prefix(W()) == W()
            True
            sage: w.longest_common_prefix(w[:3]) == w[:3]
            True
            sage: w.longest_common_prefix(w[:3]) == w
            False
            sage: Words('12')("11").longest_common_prefix(Words('1')("1"))
            word: 1
        """
        self, other = self.coerce(other)
        self_iter, other_iter = iter(self._word_content), iter(other._word_content)
        i = 0
        try:
            while True:
                if self_iter.next() != other_iter.next():
                    return self[:i]
                i += 1
        except StopIteration:  # one iterator finished
            if len(self) < len(other):
                return self
            else:
                return other

    def longest_common_suffix(self, other):
        """
        Returns the longest common suffix of self and other.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: y = W('549332345')
            sage: w = W('203945')
            sage: w.longest_common_suffix(y)
            word: 45
            sage: w.longest_common_suffix(y.reversal())
            word: 3945
            sage: w.longest_common_suffix(W())
            word:
            sage: W().longest_common_suffix(w)
            word:
            sage: W().longest_common_suffix(W())
            word:
            sage: w.longest_common_suffix(w[3:]) == w[3:]
            True
            sage: w.longest_common_suffix(w[3:]) == w
            False
        """
        self, other = self.coerce(other)
        i = 0
        while(i < min(len(self), len(other)) and self._word_content[-(i+1)] == other._word_content[-(i+1)]):
            i += 1
        if i == 0:
            return self.parent()()
        return self[-i:]

    def prefix_function_table(self):
        """
        Returns a vector containing the length of the proper
        prefix-suffixes for all the non-empty prefixes of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('121321').prefix_function_table()
            [0, 0, 1, 0, 0, 1]
            sage: W('1241245').prefix_function_table()
            [0, 0, 0, 1, 2, 3, 0]
            sage: W().prefix_function_table()
            []
        """
        k = 0
        res = list(repeat(0, len(self)))
        for q in xrange(1, len(self)):
            while k > 0 and self._word_content[k] != self._word_content[q]:
                k = res[k-1]
            if self._word_content[k] == self._word_content[q]:
                k += 1
            res[q] = k
        return res

    def length_border(self):
        """
        Returns the length of the border of self.

        The *border* of a word is the longest word that is both a proper
        prefix and a proper suffix of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('121').length_border()
            1
            sage: W('1').length_border()
            0
            sage: W('1212').length_border()
            2
            sage: W('111').length_border()
            2
            sage: W().length_border() is None
            True
        """
        if len(self) == 0:
            return None
        return self.prefix_function_table()[-1]

    def border(self):
        """
        Returns the longest word that is both a proper prefix and a proper
        suffix of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('121212').border()
            word: 1212
            sage: W('12321').border()
            word: 1
            sage: W().border() is None
            True
        """
        if len(self) == 0:
            return None
        return self[:self.length_border()]

    def minimal_period(self):
        r"""
        Returns the period of self.

        Let `A` be an alphabet. An integer `p\geq 1` is a
        *period* of a word `w=a_1a_2\cdots a_n` where
        `a_i\in A` if `a_i=a_{i+p}` for
        `i=1,\ldots,n-p`. The smallest period of `w` is
        called *the* period of `w`. See Chapter 1 of [1].

        EXAMPLES::

            sage: Word('aba').minimal_period()
            2
            sage: Word('abab').minimal_period()
            2
            sage: Word('ababa').minimal_period()
            2
            sage: Word('ababaa').minimal_period()
            5
            sage: Word('ababac').minimal_period()
            6
            sage: Word('aaaaaa').minimal_period()
            1
            sage: Word('a').minimal_period()
            1
            sage: Word().minimal_period()
            1

        REFERENCES:

        - [1] M. Lothaire, Algebraic Combinatorics On Words, vol.  90
          of Encyclopedia of Mathematics and its Applications,
          Cambridge University Press, U.K., 2002.
        """
        if self.is_empty():
            return 1
        return len(self)-self.length_border()

    def order(self):
        r"""
        Returns the order of self.

        Let `p(w)` be the period of a word `w`. The
        positive rational number `|w|/p(w)` is the *order* of
        `w`. See Chapter 8 of [1].

        OUTPUT:


        -  ``rational`` - the order


        EXAMPLES::

            sage: Word('abaaba').order()
            2
            sage: Word('ababaaba').order()
            8/5
            sage: Word('a').order()
            1
            sage: Word('aa').order()
            2
            sage: Word().order()
            0

        REFERENCES:

        - [1] M. Lothaire, Algebraic Combinatorics On Words, vol.  90
          of Encyclopedia of Mathematics and its Applications,
          Cambridge University Press, U.K., 2002.
        """
        from sage.rings.rational import Rational
        return Rational((len(self),self.minimal_period()))

    def critical_exponent(self):
        r"""
        Returns the critical exponent of self.

        The *critical exponent* of a word is the supremum of the order of
        all its (finite) factors. See [1].

        EXAMPLES::

            sage: Word('aaba').critical_exponent()
            2
            sage: Word('aabaa').critical_exponent()
            2
            sage: Word('aabaaba').critical_exponent()
            7/3
            sage: Word('ab').critical_exponent()
            1
            sage: Word('aba').critical_exponent()
            3/2
            sage: words.ThueMorseWord()[:20].critical_exponent()
            2

        REFERENCES:

        - [1] F. Dejean. Sur un théorème de Thue. J.  Combinatorial
          Theory Ser. A 13:90–99, 1972.
        """
        return max(map(FiniteWord_over_OrderedAlphabet.order, self.factor_iterator()))

    def is_overlap(self):
        """
        Returns True if self is an overlap, and False otherwise.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('12121').is_overlap()
            True
            sage: W('123').is_overlap()
            False
            sage: W('1231').is_overlap()
            False
            sage: W('123123').is_overlap()
            False
            sage: W('1231231').is_overlap()
            True
            sage: W().is_overlap()
            False
        """
        if len(self) == 0:
            return False
        return self.length_border() > len(self)/2

    def primitive_length(self):
        """
        Returns the length of the primitive of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('1231').primitive_length()
            4
            sage: W('121212').primitive_length()
            2
        """
        l = lu = len(self)
        if l == 0:
            return 0
        p = self.prefix_function_table()
        while l > 0:
            l = p[l-1]
            if lu % (lu - l) == 0:
                return lu - l

    def is_primitive(self):
        """
        Returns True if self is primitive, and False otherwise.

        A finite word `w` is *primitive* if it is not a positive
        integer power of a shorter word.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('1231').is_primitive()
            True
            sage: W('111').is_primitive()
            False
        """
        return len(self) == self.primitive_length()

    def primitive(self):
        """
        Returns the primitive of self.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('12312').primitive()
            word: 12312
            sage: W('121212').primitive()
            word: 12
        """
        return self[:self.primitive_length()]

    def exponent(self):
        """
        Returns the exponent of self.

        OUTPUT:


        -  ``integer`` - the exponent


        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('1231').exponent()
            1
            sage: W('121212').exponent()
            3
            sage: W().exponent()
            0
        """
        if len(self) == 0:
            return 0
        return len(self) / self.primitive_length()

    def colored_vector(self, x=0, y=0, width='default', height=1, cmap='hsv', thickness=1):
        r"""
        Returns a vector (Graphics object) illustrating self. Each letter
        is represented by a colored rectangle. There is a unique color for
        each letter of the alphabet.

        INPUT:


        -  ``x`` - (default: 0) bottom left x-coordinate of
           the vector

        -  ``y`` - (default: 0) bottom left y-coordinate of
           the vector

        -  ``width`` - (default: 'default') width of the
           vector. By default, the width is the length of self.

        -  ``height`` - (default: 1) height of the vector

        -  ``thickness`` - (default: 1) thickness of the
           contour

        -  ``cmap`` - (default: 'hsv') color map; type:
           ``import matplotlib.cm; matplotlib.cm.datad.keys()``
           for available colormap names.


        OUTPUT: Graphics

        EXAMPLES::

            sage: Word(range(20)).colored_vector()

        ::

            sage: Word(range(100)).colored_vector(0,0,10,1)

        ::

            sage: Words(range(100))(range(10)).colored_vector()

        ::

            sage: w = Word('abbabaab')
            sage: w.colored_vector()

        ::

            sage: w.colored_vector(cmap='autumn')

        TESTS::

            sage: Word(range(100)).colored_vector(cmap='jolies')
            Traceback (most recent call last):
            ...
            RuntimeError: Color map jolies not known
            sage: Word(range(100)).colored_vector(cmap='__doc__')
            Traceback (most recent call last):
            ...
            RuntimeError: Color map __doc__ not known
        """
        #Recognize the color map
        import matplotlib.cm as cm
        from matplotlib.colors import LinearSegmentedColormap as C
        key_error = False
        try:
            mpl_cmap = cm.__dict__[cmap]
        except KeyError:
            key_error = True

        if key_error or not isinstance(mpl_cmap, C):
            possibilities = ', '.join([str(x) for x in cm.__dict__.keys() if \
                                       isinstance(cm.__dict__[x], C)])
            import sage.misc.misc
            sage.misc.misc.verbose("The possible color maps include: %s"%possibilities, level=0)
            raise RuntimeError, "Color map %s not known"%cmap

        #Drawing the colored vector...
        from sage.plot.plot import polygon,line

        #The default width of the vector
        if width == 'default':
            width = len(self)

        #The black frame of the vector
        ymax = y + height
        L = [(x,y), (x+width,y), (x+width,ymax), (x,ymax), (x,y)]
        rep = line(L, rgbcolor=(0,0,0), thickness=thickness)

        #base : the width of each rectangle
        base = width / float(len(self))

        #A colored rectangle for each letter
        xp = x
        dim = float(self.parent().size_of_alphabet())
        for i in self._word_content:
            xq = xp + base
            L = [(xp,y), (xq,y), (xq,ymax), (xp,ymax) ]
            rgbcolor = mpl_cmap( i / dim ) [:3]
            rep += polygon(L, rgbcolor = rgbcolor)
            xp = xq
        rep.axes(False)
        return rep

    def last_position_table(self):
        """
        Returns a table (of size 256) that contains the last position of
        each letter in self. The letters not present in the word will have
        a position of None.

        EXAMPLES::

            sage: Words('01234')('1231232').last_position_table()
            [-1, 3, 6, 5, -1]
        """
        res = [-1]*self.size_of_alphabet()
        for (i, s) in izip(count(), self._word_content):
            res[s] = i
        return res

    def good_suffix_table(self):
        """
        Returns a table of the maximum skip you can do in order not to miss
        a possible occurrence of self in a word.

        This is a part of the Boyer-Moore algorithm to find factors. See
        [1].

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('121321').good_suffix_table()
            [5, 5, 5, 5, 3, 3, 1]
            sage: W('12412').good_suffix_table()
            [3, 3, 3, 3, 3, 1]

        REFERENCES:

        - [1] R.S. Boyer, J.S. Moore, A fast string searching
          algorithm, Communications of the ACM 20 (1977) 762-772.
        """
        l = len(self)
        p = self.reversal().prefix_function_table()
        res = [l - p[-1]]*(l+1)
        for i in xrange(1, l+1):
            j = l - p[i - 1]
            if res[j] > (i - p[i-1]):
                res[j] = i - p[i-1]
        return res

    def _pos_in(self, other, p):
        r"""
        Returns the position of the first occurrence of self starting at
        position p in other.

        EXAMPLES::

            sage: W = Words('123')
            sage: W('12')._pos_in(W('131231'), 2)
            2
            sage: W('12')._pos_in(W('131231'), 3) is None
            True
            sage: W('32')._pos_in(W('131231'), 0) is None
            True
        """
        lf = len(self)
        lm = len(other)
        if lf == 0 or lm == 0:
            return None
        occ = self.last_position_table()
        suff = self.good_suffix_table()
        s = p
        while s <= lm - lf:
            for j in xrange(lf-1, -1, -1):
                if self._word_content[j] != other._word_content[s+j]:
                    s += max(suff[j + 1], j - occ[other._word_content[s+j]])
                    break
            else:
                return s
        return None

    def first_pos_in(self, other):
        r"""
        Returns the position of the first occurrence of self in other, or
        None if self is not a factor of other.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('12').first_pos_in(W('131231'))
            2
            sage: W('32').first_pos_in(W('131231')) is None
            True
        """
        self, other = self.coerce(other)
        return self._pos_in(other, 0)

    def is_factor_of(self, other):
        r"""
        Returns True if self is a factor of other, and False otherwise.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('2113').is_factor_of(W('123121332131233121132123'))
            True
            sage: W('321').is_factor_of(W('1231241231312312312'))
            False
        """
        return self.first_pos_in(other) is not None

    def factor_occurrences_in(self, other):
        r"""
        Returns an iterator over all occurrences (including overlapping
        ones) of self in other in their order of appearance.

        EXAMPLES::

            sage: W = Words('123')
            sage: list(W('121').factor_occurrences_in(W('121213211213')))
            [0, 2, 8]
        """
        self, other = self.coerce(other)
        if len(self) == 0:
            raise NotImplementedError, "undefined value"
        p = self._pos_in(other, 0)
        while p is not None:
            yield p
            p = self._pos_in(other, p+1)

    def nb_factor_occurrences_in(self, other):
        r"""
        Returns the number of times self appears as a factor in other.

        EXAMPLES::

            sage: W = Words('123')
            sage: W().nb_factor_occurrences_in(W('123'))
            Traceback (most recent call last):
            ...
            NotImplementedError: undefined value
            sage: W('123').nb_factor_occurrences_in(W('112332312313112332121123'))
            4
            sage: W('321').nb_factor_occurrences_in(W('11233231231311233221123'))
            0
        """
        return len_it(self.factor_occurrences_in(other))

    def is_subword_of(self, other):
        """
        Returns True is self is a subword of other, and False otherwise.

        EXAMPLES::

            sage: W = Words('123')
            sage: W().is_subword_of(W('123'))
            True
            sage: W('123').is_subword_of(W('3211333213233321'))
            True
            sage: W('321').is_subword_of(W('11122212112122133111222332'))
            False
        """
        self, other = self.coerce(other)
        its = iter(self._word_content)
        try:
            s = its.next()
            for e in other._word_content:
                if s == e:
                    s = its.next()
            else:
                return False
        except StopIteration:
            return True

    def nb_subword_occurrences_in(self, other):
        r"""
        Returns the number of times self appears in other as a subword.

        EXAMPLES::

            sage: W = Words('1234')
            sage: W().nb_subword_occurrences_in(W('123'))
            Traceback (most recent call last):
              ...
            NotImplementedError: undefined value
            sage: W('123').nb_subword_occurrences_in(W('1133432311132311112'))
            11
            sage: W('4321').nb_subword_occurrences_in(W('1132231112233212342231112'))
            0
            sage: W('3').nb_subword_occurrences_in(W('122332112321213'))
            4
        """
        self, other = self.coerce(other)
        ls = len(self)
        if ls == 0:
            raise NotImplementedError, "undefined value"
        elif ls == 1:
            return self.nb_factor_occurrences_in(other)
        elif len(other) < ls:
            return 0
        symb = self[:1]
        suffword = other
        suffsm = self[1:]
        n = 0
        cpt = 0
        i = symb.first_pos_in(suffword)
        while i is not None:
            suffword = suffword[i+1:]
            m = suffsm.nb_subword_occurrences_in(suffword)
            if m == 0: break
            n += m
            i = symb.first_pos_in(suffword)
        return n

    def _return_words_list(self, fact):
        r"""
        Returns the return words as a list in the order they appear in the
        word.

        TESTS::

            sage: W = Words('abc')
            sage: W('baccabccbacbca')._return_words_list(W('b'))
            [word: bacca, word: bcc, word: bac]
        """
        self, fact = self.coerce(fact)
        i = fact.first_pos_in(self)
        if i is None:
            return []
        w = self[i+1:]
        j = fact.first_pos_in(w)
        res = []
        while j is not None:
            res.append(self[i:i+j+1])
            w = w[j+1:]
            i += j+1
            j = fact.first_pos_in(w)
        return res

    def return_words(self, fact):
        r"""
        Returns the set of return words of fact in self.

        This is the set of all factors starting by the given factor and
        ending just before the next occurrence of this factor. See [1] and
        [2].

        EXAMPLES::

            sage: W = Words('123')
            sage: W('21331233213231').return_words(W('2')) == set([W('21331'), W('233'), W('213')])
            True
            sage: W().return_words(W('213')) == set()
            True
            sage: W('121212').return_words(W('1212')) == set([W('12')])
            True

        REFERENCES:

        - [1] F. Durand, A characterization of substitutive sequences
          using return words, Discrete Math. 179 (1998) 89-101.

        - [2] C. Holton, L.Q. Zamboni, Descendants of primitive
          substitutions, Theory Comput. Syst. 32 (1999) 133-157.
        """
        return set(self._return_words_list(fact))

    def complete_return_words(self, fact):
        """
        Returns the set of complete return words of fact in self.

        This is the set of all factors starting by the given factor and
        ending just after the next occurrence of this factor. See for
        instance [1].

        EXAMPLES::

            sage: W = Words('123')
            sage: W('21331233213231').complete_return_words(W('2')) == set([W('213312'), W('2332'), W('2132')])
            True
            sage: W('').complete_return_words(W('213')) == set()
            True
            sage: W('121212').complete_return_words(W('1212')) == set([W('121212')])
            True

        REFERENCES:

        - [1] J. Justin, L. Vuillon, Return words in Sturmian and
          episturmian words, Theor. Inform. Appl. 34 (2000) 343-356.
        """
        i = fact.first_pos_in(self)
        if i is None:
            return set()
        w = self[i+1:]
        j = fact.first_pos_in(w)
        res = set()
        while j is not None:
            res.add(self[i:i+j+len(fact)+1])
            w = w[j+1:]
            i += j+1
            j = fact.first_pos_in(w)
        return res

    def return_words_derivate(self, fact, W=None):
        r"""
        Returns the word generated by mapping a letter to each occurrence
        of the return words for the given factor dropping any dangling
        prefix and suffix.

        The optional set of words parameter must be over an alphabet that
        contains as much letters as there are different return words in
        self, otherwise there will be some breakage in the function. The
        default value for this parameter always respects this property. The
        letters are attributed to the words in the order they are
        discovered.

        EXAMPLES::

            sage: Words('123')('12131221312313122').return_words_derivate(Words('1')('1'))
            word: 123242

        REFERENCES:

        - [1] F. Durand, A characterization of substitutive
          sequences using return words, Discrete Math. 179 (1998) 89-101.
        """
        self, fact = self.coerce(fact)
        idx = 0
        tab = {}
        ret = map(lambda w: tab.setdefault(w, len(tab)) + 1, self._return_words_list(fact))
        if W is None:
            W = Words(xrange(1, len(tab)+1))
        return W(ret)

    def delta(self):
        """
        Returns the delta equivalent of self.

        This is the word composed of the length of consecutive runs of the
        same letter in a given word.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('22112122').delta()
            word: 22112
            sage: W('555008').delta()
            word: 321
            sage: W().delta()
            word:
        """
        if len(self) == 0:
            return Words([])()
        ss = self[0]
        c = 0
        v = list()
        max_c = 0
        for s in self:
            if s == ss:
                c += 1
                if c > max_c:
                    max_c = c
            else:
                v.append(c)
                ss = s
                c = 1
        v.append(c)
        return Words(xrange(1, max_c+1))(v)

    def delta_inv(self, W=None, s=None):
        """
        Returns the inverse of the delta operator applied to self.

        The letters in the returned word will start at the specified letter
        or the first one if None is specified (the default). The default
        alphabet is [1, 2].

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 1, 1]).delta_inv()
            word: 112212
            sage: W([1, 1, 1, 1]).delta_inv(Words('123'))
            word: 1231
            sage: W([2, 2, 1, 1, 2]).delta_inv(s=2)
            word: 22112122
        """
        if not all(imap(isint, self.alphabet())):
            raise ValueError, "delta_inv() can only be applied to word composed of integer letters"
        if W is None:
            W = Words([1, 2])
        if len(self) == 0:
            return W()
        if not is_Words(W):
            raise TypeError, "W must be an instance of Words"
        if s is None:
            s = W.alphabet().unrank(0)
        if s not in W.alphabet():
            raise ValueError, "starting letter not in alphabet"
        v = []

        p = W.alphabet().rank(s)
        al = cycle(imap(W.alphabet().unrank, chain(xrange(p, W.size_of_alphabet()), xrange(p))))
        al.next()
        for e in self:
            v += ([s] * e)
            s = al.next()
        return W(v)

    def delta_derivate(self, W=None):
        """
        Returns the derivative under delta for self.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate()
            word: 22
            sage: W('1').delta_derivate(Words([1]))
            word: 1
            sage: W('2112').delta_derivate()
            word: 2
            sage: W('2211').delta_derivate()
            word: 22
            sage: W('112').delta_derivate()
            word: 2
            sage: W('11222').delta_derivate(Words([1, 2, 3]))
            word: 3
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[0] != W.alphabet().last():
            d = d[1:]
        if d[-1] != W.alphabet().last():
            d = d[:-1]
        return d

    def delta_derivate_left(self, W=None):
        """
        Returns the derivative under delta for self.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate_left()
            word: 22
            sage: W('1').delta_derivate_left(Words([1]))
            word: 1
            sage: W('2112').delta_derivate_left()
            word: 21
            sage: W('2211').delta_derivate_left()
            word: 22
            sage: W('112').delta_derivate_left()
            word: 21
            sage: W('11222').delta_derivate_left(Words([1, 2, 3]))
            word: 3
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[0] != W.alphabet().last():
            d = d[1:]
        return d

    def delta_derivate_right(self, W=None):
        """
        Returns the right derivative under delta for self.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate_right()
            word: 122
            sage: W('1').delta_derivate_right(Words([1]))
            word: 1
            sage: W('2112').delta_derivate_right()
            word: 12
            sage: W('2211').delta_derivate_right()
            word: 22
            sage: W('112').delta_derivate_right()
            word: 2
            sage: W('11222').delta_derivate_right(Words([1, 2, 3]))
            word: 23
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[-1] != W.alphabet().last():
            d = d[:-1]
        return d

    def phi(self):
        r"""
        Applies the phi function to self and returns the result.

        See for instance [1] and [2].

        INPUT:


        -  ``self`` - must be a word over integers


        OUTPUT:


        -  ``word`` - the result of the phi function


        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2]).phi()
            word: 222222
            sage: W().phi()
            word:
            sage: W([2, 1, 2, 2, 1, 2, 2, 1, 2, 1]).phi()
            word: 212113

        REFERENCES:

        - [1] S. Brlek, A. Ladouceur, A note on differentiable
          palindromes, Theoret. Comput. Sci. 302 (2003) 167-178.

        - [2] S. Brlek, S. Dulucq, A. Ladouceur, L. Vuillon,
          Combinatorial properties of smooth infinite words,
          Theoret. Comput. Sci. 352 (2006) 306-317.
        """
        if self.is_empty():
            return self
        m = self
        v = []
        s_max = 0
        while len(m) > 1:
            v.append(m[0])
            if m[0] > s_max:
                s_max = m[0]
            m = m.delta()
        v.append(m[0])
        if m[0] > s_max:
            s_max = m[0]
        return Words(xrange(1, s_max+1))(v)

    def phi_inv(self, W=None):
        r"""
        Applied the inverse of the phi function and returns the result.

        INPUT:


        -  ``self`` - must be a word over the integers

        -  ``W`` - the set of words of the result (must also be
           over the integers)


        OUTPUT:


        -  ``word`` - the inverse of the phi function


        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 2, 2, 1, 2]).phi_inv()
            word: 22112122
            sage: W([2, 2, 2]).phi_inv(Words([2, 3]))
            word: 2233
        """
        if W is None:
            W = self.parent()
        if self.is_empty():
            return W()
        v = self.parent()((self[-1],))
        for i in xrange(len(self) - 2, -1, -1):
            v = v.delta_inv(W, self[i])
        return v

    def _phi_inv_tab(self, tab):
        r"""
        Specialized version of phi_inv() for long or incremental words.

        TESTS::

            sage: Words([1, 2])([1, 1, 2, 2])._phi_inv_tab([2])
            word: 12211
        """
        res = self.delta_inv(s=tab[0])
        res = res[1:]
        for i in xrange(1, len(tab)):
            res = res.delta_inv(s=tab[i])
        return res

    def is_smooth_prefix(self):
        r"""
        Returns True if self is the prefix of a smooth word, and False
        otherwise.

        Let `A_k = \{1, \ldots ,k\}`, `k \geq 2`. An
        infinite word `w` in `A_k^\omega` is said to be
        *smooth* if and only if for all positive integers `m`,
        `\Delta^m(w)` is in `A_k^\omega`, where
        `\Delta(w)` is the word obtained from `w` by
        composing the length of consecutive runs of the same letter in
        `w`. See for instance [1] and [2].

        INPUT:


        -  ``self`` - must be a word over the integers to get
           something other than False


        OUTPUT:


        -  ``boolean`` - whether self is a smooth prefix or
           not


        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([1, 1, 2, 2, 1, 2, 1, 1]).is_smooth_prefix()
            True
            sage: W([1, 2, 1, 2, 1, 2]).is_smooth_prefix()
            False

        REFERENCES:

        - [1] S. Brlek, A. Ladouceur, A note on differentiable
          palindromes, Theoret. Comput. Sci. 302 (2003) 167-178.

        - [2] S.  Brlek, S. Dulucq, A. Ladouceur, L. Vuillon,
          Combinatorial properties of smooth infinite words,
          Theoret. Comput. Sci. 352 (2006) 306-317.
        """
        m = self
        while len(m) > 1:
            m = m.delta_derivate_right()
            if m not in self.parent():
                return False
        return True

    def is_lyndon(self):
        """
        Returns True if self is a Lyndon word, and False otherwise.

        A *Lyndon word* is a non-empty word that is lexicographically
        smaller than all of its proper suffixes for the given order on its
        alphabet. That is, `w` is a Lyndon word if `w` is
        non-empty and for each factorization `w = uv` (with
        `u`, `v` both non-empty), we have `w < v`.

        Equivalently, `w` is a Lyndon word iff `w` is a
        non-empty word that is lexicographically smaller than all of its
        proper conjugates for the given order on its alphabet.

        See for instance [1].

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('123132133').is_lyndon()
            True
            sage: W().is_lyndon()
            True
            sage: W('122112').is_lyndon()
            False

        REFERENCES:

        - [1] M. Lothaire, Combinatorics On Words, vol. 17 of
          Encyclopedia of Mathematics and its Applications,
          Addison-Wesley, Reading, Massachusetts, 1983.
        """
        if len(self) == 0:
            return True
        s = self._word_content[0]
        for (e, i) in izip(self[1:]._word_content, count(1)):
            if s < e:
                continue
            if not self < self[i:]:
                return False
        return True

    def BWT(self):
        """
        Returns the Burrows-Wheeler Transform (BWT) of self.

        The *Burrows-Wheeler transform* of a finite word `w` is
        obtained from `w` by first listing the conjugates of
        `w` in lexicographic order and then concatenating the final
        letters of the conjugates in this order. See [1].

        EXAMPLES::

            sage: W = Words('abc')
            sage: W('abaccaaba').BWT()
            word: cbaabaaca
            sage: W('abaab').BWT()
            word: bbaaa
            sage: W('bbabbaca').BWT()
            word: cbbbbaaa
            sage: W('aabaab').BWT()
            word: bbaaaa
            sage: Word().BWT()
            word:
            sage: W('a').BWT()
            word: a

        REFERENCES:

        - [1] M. Burrows, D.J. Wheeler, "A block-sorting lossless
          data compression algorithm", HP Lab Technical Report, 1994,
          available at
          http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.html
        """
        if self.is_empty():
           return self
        conjugates = self._conjugates_list()
        conjugates.sort()
        return self.parent()([x[len(x)-1] for x in conjugates])

    def _duval_algorithm(self):
        """
        TESTS::

            sage: Words('01')('010010010001000')._duval_algorithm()
            (01.001.001.0001)
            sage: Words('123')('122113112212')._duval_algorithm()
            (122.113.112212)

        REFERENCES:

        - [1] J.-P. Duval, Factorizing words over an ordered
          alphabet, J. Algorithms 4 (1983), no. 4, 363-381.
        """
        t = Factorization()
        cm = self
        c = iter(cm._word_content)
        cc = iter(cm._word_content)
        cc.next()
        i = 0
        d = 1
        j = k = 1
        l = len(cm)

        while k < l:
            c, c_s = peek_it(c)
            cc, cc_s = peek_it(cc)
            if c_s < cc_s:
                cc.next()
                k += 1
                j = k
                d = k - i
                c = iter(cm._word_content)
            elif c_s == cc_s:
                cc.next()
                k += 1
                if (k - j) == d:
                    j = k
                    c = iter(cm._word_content)
                else:
                    c.next()
            else:
                i += d
                while i <= j:
                    i += d
                    t.append(cm[:d])
                    cm = cm[d:]
                c = iter(cm._word_content)
                cc = iter(cm._word_content)
                cc.next()
                i = j
                j += 1
                k = j
                d = 1
        i += d
        while i <= j:
            i += d
            t.append(cm[:d])
            cm = cm[d:]
        return t

    def lyndon_factorization(self):
        r"""
        Returns the Lyndon factorization of self.

        The *Lyndon factorization* of a finite word `w` is the
        unique factorization of `w` as a non-increasing product of
        Lyndon words, i.e., `w = l_1\cdots l_n` where each
        `l_i` is a Lyndon word and
        `l_1\geq \cdots \geq l_n`. See for instance [1].

        OUTPUT:


        -  ``list`` - the list of factors obtained


        EXAMPLES::

            sage: Words('01')('010010010001000').lyndon_factorization()
            (01.001.001.0001.0.0.0)
            sage: Words('10')('010010010001000').lyndon_factorization()
            (0.10010010001000)
            sage: Words('ab')('abbababbaababba').lyndon_factorization()
            (abb.ababb.aababb.a)
            sage: Words('ba')('abbababbaababba').lyndon_factorization()
            (a.bbababbaaba.bba)

        TESTS::

            sage: Words('01')('01').lyndon_factorization()
            (01)
            sage: Words('10')('01').lyndon_factorization()
            (0.1)
            sage: lynfac = Words('ab')('abbababbaababba').lyndon_factorization()
            sage: map(lambda x:x.is_lyndon(), lynfac)
            [True, True, True, True]
            sage: lynfac = Words('ba')('abbababbaababba').lyndon_factorization()
            sage: map(lambda x:x.is_lyndon(), lynfac)
            [True, True, True]

        REFERENCES:

        - [1] J.-P. Duval, Factorizing words over an ordered
          alphabet, J. Algorithms 4 (1983) 363-381.
        """
        tab = self._duval_algorithm()
        l = sum(imap(len, tab))
        if l < len(self):
            tab += self[l:]._duval_algorithm()
        return tab

    def standard_factorization(self):
        r"""
        Returns the standard factorization of self.

        The *standard factorization* of a word `w` is the unique
        factorization: `w = uv` where `v` is the longest
        proper suffix of `w` that qualifies as a Lyndon word.

        Note that if `w` is a Lyndon word with standard
        factorization `w = uv`, then `u` and `v`
        are also Lyndon words and `u < v`.

        See for instance [1] and [2].

        OUTPUT:


        -  ``list`` - the list of factors


        EXAMPLES::

            sage: Words('01')('0010110011').standard_factorization()
            (001011.0011)
            sage: Words('123')('1223312').standard_factorization()
            (12233.12)

        REFERENCES:

        - [1] K.-T. Chen, R.H. Fox, R.C. Lyndon, Free differential
          calculus, IV. The quotient groups of the lower central
          series, Ann. of Math. 68 (1958) 81-95.

        - [2] J.-P. Duval, Factorizing words over an ordered alphabet,
          J. Algorithms 4 (1983) 363-381.
        """
        suff = self[1:]
        for l in xrange(1, len(self)):
            pref = self[:l]
            if pref.is_lyndon() and suff.is_lyndon():
                return Factorization([pref, suff])
            suff = suff[1:]
        return Factorization([self, self.parent()()])

    def standard_factorization_of_lyndon_factorization(self):
        r"""
        Returns the standard factorization of the Lyndon factorization of
        self.

        OUTPUT:


        -  ``list of lists`` - the factorization


        EXAMPLES::

            sage: Words('123')('1221131122').standard_factorization_of_lyndon_factorization()
            [(12.2), (1.13), (1.122)]
        """
        return map(FiniteWord_over_OrderedAlphabet.standard_factorization, self.lyndon_factorization())

    def crochemore_factorization(self):
        r"""
        Returns the Crochemore factorization of self as an ordered list of
        factors.

        The *Crochemore factorization* of a finite word `w` is the
        unique factorization: `(x_1, x_2, \ldots, x_n)` of
        `w` with each `x_i` satisfying either: C1.
        `x_i` is a letter that does not appear in
        `u = x_1\ldots x_{i-1}`; C2. `x_i` is the
        longest prefix of `v = x_i\ldots x_n` that also has an
        occurrence beginning within `u = x_1\ldots x_{i-1}`. See
        [1].

        This is not a very good implementation, and should be improved.

        EXAMPLES::

            sage: x = Words('ab')('abababb')
            sage: x.crochemore_factorization()
            (a.b.abab.b)
            sage: mul(x.crochemore_factorization()) == x
            True
            sage: y = Words('abc')('abaababacabba')
            sage: y.crochemore_factorization()
            (a.b.a.aba.ba.c.ab.ba)
            sage: mul(y.crochemore_factorization()) == y
            True
            sage: x = Words([0, 1])([0,1,0,1,0,1,1])
            sage: x.crochemore_factorization()
            (0.1.0101.1)
            sage: mul(x.crochemore_factorization()) == x
            True

        REFERENCES:

        - [1] M. Crochemore, Recherche linéaire d'un carré dans un
          mot, C. R. Acad. Sci. Paris Sér. I Math. 296 (1983) 14
          781-784.
        """
        c = Factorization([self[:1]])
        u = self[:sum(map(len,c))] # = x_1 ... x_{i-1}
        v = self[sum(map(len,c)):] # = x_i ... x_n
        while v:
            # C1. x_i is a letter that does not appear in u = x_1...x_{i-1}
            if v[0] not in u:
                c.append(v[:1])
            else:
            # C2. x_i is the longest prefix of v = x_i...x_n that also has an
            #     occurrence beginning within u = x_1...x_{i-1}.
                xi = v
                while True:
                    if xi.first_pos_in(self) < len(u):
                        c.append(xi)
                        break
                    else:
                        xi = xi[:-1]
            u = self[:sum(map(len,c))] # = x_1 ... x_{i-1}
            v = self[sum(map(len,c)):] # = x_i ... x_n
        return c

    def is_balanced(self, q=1):
        r"""
        Returns True if self is `q`-balanced, and False otherwise.

        A finite or infinite word `w` is said to be
        *`q`-balanced* if for any two factors `u`,
        `v` of `w` of the same length, the difference
        between the number of `x`'s in each of `u` and
        `v` is at most `q` for all letters `x` in
        the alphabet of `w`. A `1`-balanced word is simply
        said to be balanced. See for instance [1] and Chapter 2 of [2].

        INPUT:


        -  ``q`` - integer (default 1), the balance level


        OUTPUT:


        -  ``boolean`` - the result


        EXAMPLES::

            sage: Words('123')('1213121').is_balanced()
            True
            sage: Words('12')('1122').is_balanced()
            False
            sage: Words('123')('121333121').is_balanced()
            False
            sage: Words('123')('121333121').is_balanced(2)
            False
            sage: Words('123')('121333121').is_balanced(3)
            True
            sage: Words('12')('121122121').is_balanced()
            False
            sage: Words('12')('121122121').is_balanced(2)
            True
            sage: Words('12')('121122121').is_balanced(-1)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Words('12')('121122121').is_balanced(0)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Words('12')('121122121').is_balanced('a')
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer

        REFERENCES:

        - [1] J. Cassaigne, S. Ferenczi, L.Q. Zamboni, Imbalances
          in Arnoux-Rauzy sequences, Ann. Inst. Fourier (Grenoble) 50 (2000)
          1265-1276.

        - [2] M. Lothaire, Algebraic Combinatorics On Words, vol.  90
          of Encyclopedia of Mathematics and its Applications,
          Cambridge University Press, U.K., 2002.
        """
        if not isint(q) or q <= 0:
           raise TypeError, "the balance level must be a positive integer"
        for i in xrange(2, len(self)):
            tab = [None] * self.size_of_alphabet()
            for j in xrange(len(tab)):
                tab[j] = set()
            for fact in self.factor_iterator(i):
                for (n, sym) in izip(count(), self.alphabet()):
                    tab[n].add(self.parent()(sym).nb_factor_occurrences_in(fact))
            for t in tab:
                if len(t) > q+1:
                    return False
        return True

    def freq(self):
        r"""
        Returns a table of the frequencies of the letters in self.

        OUTPUT:


        -  ``dict`` - letters associated to their frequency


        EXAMPLES::

            sage: Words('123')('1213121').freq()    # keys appear in random order
            {'1': 4, '2': 2, '3': 1}

        TESTS::

            sage: f = Words('123')('1213121').freq()
            sage: f['1'] == 4
            True
            sage: f['2'] == 2
            True
            sage: f['3'] == 1
            True
        """
        res = dict()
        for sym in self.alphabet():
            res[sym] = 0
        for s in self:
            res[s] += 1
        return res

    def parikh_vector(self):
        """
        Returns the Parikh vector of self, i.e., the vector containing the
        number of occurrences of each letter, given in the order of the
        alphabet.

        See also evaluation, which returns a list truncated after the last
        nonzero entry.

        EXAMPLES::

            sage: Word('aabaa').parikh_vector()
            [4, 1]
            sage: Word('aabaacaccab').parikh_vector()
            [6, 2, 3]
            sage: Words('abc')('aabaa').parikh_vector()
            [4, 1, 0]
            sage: Word().parikh_vector()
            []
            sage: Word('a').parikh_vector()
            [1]
            sage: Words('abc')('a').parikh_vector()
            [1, 0, 0]
            sage: Words('ab')().parikh_vector()
            [0, 0]
            sage: Words('abc')().parikh_vector()
            [0, 0, 0]
            sage: Words('abcd')().parikh_vector()
            [0, 0, 0, 0]

        TESTS::

            sage: P = Alphabet(name="positive integers")
            sage: w = Words(P)(range(1,10))
            sage: w.parikh_vector()
            Traceback (most recent call last):
            ...
            TypeError: the alphabet is infinite; use evaluation() instead
        """
        n = self.parent().size_of_alphabet()
        if not isinstance(n, (int,Integer)):
            raise TypeError, "the alphabet is infinite; use evaluation() instead"
        v = [0]*n
        for c in self._word_content:
            v[c] += 1
        return v

    def apply_morphism(self,morphism):
        """
        Returns the word obtained by applying the morphism to self.

        INPUT:


        -  ``morphism`` - Can be an instance of WordMorphism,
           or anything that can be used to construct one.


        EXAMPLES::

            sage: w = Word("ab")
            sage: d = {'a':'ab', 'b':'ba'}
            sage: w.apply_morphism(d)
            word: abba
            sage: w.apply_morphism(WordMorphism(d))
            word: abba
        """
        from sage.combinat.words.morphism import WordMorphism
        if not isinstance(morphism, WordMorphism):
            morphism = WordMorphism(morphism)
        return morphism(self)

    def count(self, letter):
        """
        Counts the number of occurrences of letter in self.

        EXAMPLES::

            sage: Words('ab')('abbabaab').count('a')
            4
        """
        return self.parent()([letter]).nb_factor_occurrences_in(self)

    def suffix_trie(self):
        r"""
        Returns the suffix trie of self.

        The *suffix trie* of a finite word `w` is a data structure
        representing the factors of `w`. It is a tree whose edges
        are labelled with letters of `w`, and whose leafs
        correspond to suffixes of `w`.

        See sage.combinat.words.suffix_trees.SuffixTrie? for more
        information.

        EXAMPLES::

            sage: w = Words("cao")("cacao")
            sage: w.suffix_trie()
            Suffix Trie of the word: cacao

        ::

            sage: w = Words([0,1])([0,1,0,1,1])
            sage: w.suffix_trie()
            Suffix Trie of the word: 01011
        """
        from sage.combinat.words.suffix_trees import SuffixTrie
        return SuffixTrie(self)

    def implicit_suffix_tree(self):
        r"""
        Returns the implicit suffix tree of self.

        The *suffix tree* of a word `w` is a compactification of
        the suffix trie for `w`. The compactification removes all
        nodes that have exactly one incoming edge and exactly one outgoing
        edge. It consists of two components: a tree and a word. Thus,
        instead of labelling the edges by factors of `w`, we can
        labelled them by indices of the occurrence of the factors in
        `w`.

        See sage.combinat.words.suffix_trees.ImplicitSuffixTree? for more
        information.

        EXAMPLES::

            sage: w = Words("cao")("cacao")
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: cacao

        ::

            sage: w = Words([0,1])([0,1,0,1,1])
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: 01011
        """
        from sage.combinat.words.suffix_trees import ImplicitSuffixTree
        return ImplicitSuffixTree(self)

    def suffix_tree(self):
        r"""
        Alias for implicit_suffix_tree().

        EXAMPLES::

            sage: Words('ab')('abbabaab').suffix_tree()
            Implicit Suffix Tree of the word: abbabaab
        """
        return self.implicit_suffix_tree()

    def number_of_factors(self,n=None):
        r"""
        Counts the number of distinct factors of self.

        INPUT:


        -  ``n`` - an integer, or None.


        OUTPUT: If n is an integer, returns the number of distinct factors
        of length n. If n is None, returns the total number of distinct
        factors.

        EXAMPLES::

            sage: w = Words([1,2,3])([1,2,1,2,3])
            sage: w.number_of_factors()
            13
            sage: map(w.number_of_factors, range(6))
            [1, 3, 3, 3, 2, 1]

        ::

            sage: Words('123')('1213121').number_of_factors()
            22
            sage: Words('123')('1213121').number_of_factors(1)
            3

        ::

            sage: Words('a')('a'*100).number_of_factors()
            101
            sage: Words('a')('a'*100).number_of_factors(77)
            1

        ::

            sage: Words([])().number_of_factors()
            1
            sage: Words([])().number_of_factors(17)
            0

        ::

            sage: blueberry = Word("blueberry")
            sage: blueberry.number_of_factors()
            43
            sage: map(blueberry.number_of_factors, range(10))
            [1, 6, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        return self.suffix_tree().number_of_factors(n)

    def factor_iterator(self,n=None):
        r"""
        Generates distinct factors of self.

        INPUT:


        -  ``n`` - an integer, or None.


        OUTPUT: If n is an integer, returns an iterator over all distinct
        factors of length n. If n is None, returns an iterator generating
        all distinct factors.

        EXAMPLES::

            sage: w = Words('123')('1213121')
            sage: sorted( w.factor_iterator(0) )
            [word: ]
            sage: sorted( w.factor_iterator(10) )
            []
            sage: sorted( w.factor_iterator(1) )
            [word: 1, word: 2, word: 3]
            sage: sorted( w.factor_iterator(4) )
            [word: 1213, word: 1312, word: 2131, word: 3121]
            sage: sorted( w.factor_iterator() )
            [word: , word: 1, word: 12, word: 121, word: 1213, word: 12131, word: 121312, word: 1213121, word: 13, word: 131, word: 1312, word: 13121, word: 2, word: 21, word: 213, word: 2131, word: 21312, word: 213121, word: 3, word: 31, word: 312, word: 3121]

        ::

            sage: u = Words([1,2,3])([1,2,1,2,3])
            sage: sorted( u.factor_iterator(0) )
            [word: ]
            sage: sorted( u.factor_iterator(10) )
            []
            sage: sorted( u.factor_iterator(1) )
            [word: 1, word: 2, word: 3]
            sage: sorted( u.factor_iterator(5) )
            [word: 12123]
            sage: sorted( u.factor_iterator() )
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]

        ::

            sage: xxx = Word("xxx")
            sage: sorted( xxx.factor_iterator(0) )
            [word: ]
            sage: sorted( xxx.factor_iterator(4) )
            []
            sage: sorted( xxx.factor_iterator(2) )
            [word: xx]
            sage: sorted( xxx.factor_iterator() )
            [word: , word: x, word: xx, word: xxx]

        ::

            sage: e = Word()
            sage: sorted( e.factor_iterator(0) )
            [word: ]
            sage: sorted( e.factor_iterator(17) )
            []
            sage: sorted( e.factor_iterator() )
            [word: ]

        TESTS::

            sage: type( Words('cao')('cacao').factor_iterator() )
            <type 'generator'>
        """
        return self.suffix_tree().factor_iterator(n)

    def factor_set(self):
        r"""
        Returns the set of factors of self.

        EXAMPLES::

            sage: Words('123')('1213121').factor_set()   # random
            Set of elements of <generator object at 0xa8fde6c>
            sage: sorted(  Words([1,2,3])([1,2,1,2,3]).factor_set()  )
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]
            sage: sorted(  Words("x")("xx").factor_set()  )
            [word: , word: x, word: xx]
            sage: set( Words([])().factor_set() )
            set([word: ])
        """
        return Set(set(self.factor_iterator()))

    def overlap_partition(self, other, delay=0, p=None):
        """
        Returns the partition of the alphabet induced by the equivalence
        relation defined below.

        Let `u = u_0 u_1 \cdots u_{n-1}`,
        `v = v_0v_1\cdots v_{m-1}` be two words on the alphabet
        `A` where `u_i`, `v_j` are letters and
        let `d` be an integer. We define a relation
        `R_{u,v,d}\subset A \times A` by
        `R_{u,v,d} = \{ (u_k, v_{k-d}) : 0 \leq k < n, 0\leq k-d < m \}`.
        The equivalence relation returned is the symmetric, reflexive and
        transitive closure of `R_{self,other,delay}\cup p`
        (inspired from [1]).

        EXAMPLE ILLUSTRATING THE PRECEDENT DEFINITION: Let
        `A = \{a, b, c, d, e, f, h, l, v \}`, `u=cheval`,
        `v=abcdef` and `d=3`. Then `0 \leq k < 6`
        and `0\leq k-3 < 6` implies that `3\leq k \leq 5`.
        Then,

        .. math::

           R_{u,v,d} = \{ (u_3, v_0), (u_4, v_1), (u_5, v_2) \}
                     = \{ (v, a), (a, b), (l, c) \}

        These three couples correspond to the pairs of letters one above
        the other in the following overlap::

           cheval
              abcdef

        The symmetric, reflexive and transitive closure of `R_{u,v,d}`
        defines the following partition of the alphabet `A`:

        .. math::

           \{\{a, b, v\}, \{c, l\}, \{d\}, \{e\}, \{f\}, \{h\}\}.

        INPUT:


        -  ``other`` - word on the same alphabet as self

        -  ``delay`` - integer

        -  ``p`` - Set (default: None), a partition of the
           alphabet


        OUTPUT:


        -  ``p`` - Set, a set partition of the alphabet of self
           and other.


        EXAMPLES::

            sage: W = Words('abcdefhlv')
            sage: cheval = W('cheval')
            sage: abcdef = W('abcdef')
            sage: p = cheval.overlap_partition(abcdef,3); p
            {{'f'}, {'e'}, {'d'}, {'a', 'b', 'v'}, {'c', 'l'}, {'h'}}
            sage: cheval.overlap_partition(abcdef,2,p)
            {{'f'}, {'a', 'c', 'b', 'e', 'd', 'v', 'l'}, {'h'}}
            sage: W = Words('abcdef')
            sage: w = W('abc')
            sage: y = W('def')
            sage: w.overlap_partition(y, -3)
            {{'f'}, {'e'}, {'d'}, {'b'}, {'a'}, {'c'}}
            sage: w.overlap_partition(y, -2)
            {{'a', 'f'}, {'e'}, {'d'}, {'c'}, {'b'}}
            sage: w.overlap_partition(y, -1)
            {{'a', 'e'}, {'d'}, {'c'}, {'b', 'f'}}
            sage: w.overlap_partition(y, 0)
            {{'b', 'e'}, {'a', 'd'}, {'c', 'f'}}
            sage: w.overlap_partition(y, 1)
            {{'c', 'e'}, {'f'}, {'b', 'd'}, {'a'}}
            sage: w.overlap_partition(y, 2)
            {{'f'}, {'e'}, {'b'}, {'a'}, {'c', 'd'}}
            sage: w.overlap_partition(y, 3)
            {{'f'}, {'e'}, {'d'}, {'b'}, {'a'}, {'c'}}
            sage: w.overlap_partition(y, 4)
            {{'f'}, {'e'}, {'d'}, {'b'}, {'a'}, {'c'}}
            sage: W = Words(range(2))
            sage: w = W([0,1,0,1,0,1]); w
            word: 010101
            sage: w.overlap_partition(w, 0)
            {{0}, {1}}
            sage: w.overlap_partition(w, 1)
            {{0, 1}}

        TESTS::

            sage: Word().overlap_partition(Word(),'yo')
            Traceback (most recent call last):
            ...
            TypeError: delay (type given: <type 'str'>) must be an integer
            sage: Word().overlap_partition(Word(),2,'yo')
            Traceback (most recent call last):
            ...
            TypeError: p(=yo) is not a Set
            sage: Word('a').overlap_partition(Word('b'),0)
            Traceback (most recent call last):
            ...
            TypeError: no coercion rule between Ordered Alphabet ['a'] and Ordered Alphabet ['b']

        REFERENCES:

        - [1] S. Labbé, Propriétés combinatoires des
          `f`-palindromes, Mémoire de maîtrise en Mathématiques,
          Montréal, UQAM, 2008, 109 pages.
        """
        if not isint(delay):
            raise TypeError, \
                  "delay (type given: %s) must be an integer"%type(delay)
        elif delay < 0:
            return other.overlap_partition(self, -delay, p)

        self, other = self.coerce(other)

        if p is None:
            R = [[x] for x in self.alphabet()]
        else:
            if not is_Set(p):
                raise TypeError, "p(=%s) is not a Set" % p
            if Set(self.alphabet().list()) != reduce(lambda x,y:x.union(y), p):
                raise TypeError, "p(=%s) is not a partition of the alphabet" % p
            R = map(list,p)

        d, n, m = delay, len(self), len(other)
        S = [(self[k], other[k-d]) for k in range(max(0,d), min(n,m+d))]

        for (a,b) in S:
            index_a, index_b = -1, -1
            for (i, r) in enumerate(R):
                if index_a < 0 and a in r: index_a = i
                if index_b < 0 and b in r: index_b = i
                if index_a >= 0 and index_b >= 0: break
            if index_a != index_b:
                set_a = R.pop(max(index_a, index_b))
                set_b = R.pop(min(index_a, index_b))
                R += [set_a + set_b]
        return Set(map(Set, R))

    def evaluation(self):
        r"""
        Returns a list a where a[i] is the number of occurrences in self of
        the i-th letter in self.alphabet().

        NOTE: This is slightly different from self.parikh_vector() in that
        it truncates the list after the last nonzero output.

        EXAMPLES::

            sage: Words('ab')('aabaa').evaluation()
            [4, 1]
            sage: Words('bac')('aabaa').evaluation()
            [1, 4]
            sage: Words('bca')('aabaa').evaluation()
            [1, 0, 4]
            sage: Words('abc')('aabaa').evaluation()
            [4, 1]
            sage: Word('aabaacaccab').evaluation()
            [6, 2, 3]
            sage: Word().evaluation()
            []
            sage: Words('abcde')('badbcdb').evaluation()
            [1, 3, 1, 2]
            sage: Word([1,2,2,1,3]).evaluation()
            [2, 2, 1]
            sage: P = Alphabet(name="positive integers")
            sage: Words(P)([]).evaluation()
            []
        """
        if list(self) == []:
            return []
        d = self.evaluation_dict()
        m = max(self, key=self.alphabet().rank)
        ev = []
        for a in self.alphabet():
            ev.append(d.get(a,0))
            if a == m:
                break
        return ev

    def evaluation_dict(self):
        r"""
        Returns a dictionary keyed by the letters occurring in self with
        values the number of occurrences of the letter.

        EXAMPLES::

            sage: Word([2,1,4,2,3,4,2]).evaluation_dict()
            {1: 1, 2: 3, 3: 1, 4: 2}
            sage: Word('badbcdb').evaluation_dict()
            {'a': 1, 'c': 1, 'b': 3, 'd': 2}
            sage: Word().evaluation_dict()
            {}
        """
        d = {}
        alphabet = self.alphabet()
        for c in self._word_content:
            a = alphabet.unrank(c)
            d[a] = d.setdefault(a,0) + 1
        return d

    def evaluation_sparse(self):
        r"""
        Returns a list representing the evaluation of self. The entries of
        the list are two-element lists [a, n], where a is a letter
        occurring in self and n is the number of occurrences of a in self.

        EXAMPLES::

            sage: Word([4,4,2,5,2,1,4,1]).evaluation_sparse()
            [(1, 2), (2, 2), (4, 3), (5, 1)]
            sage: Words("acdb")("abcaccab").evaluation_sparse()
            [('a', 3), ('c', 3), ('b', 2)]
        """
        alphabet_rank = self.alphabet().rank
        sort_fcn = lambda a,b: cmp(alphabet_rank(a[0]), alphabet_rank(b[0]))
        return sorted(self.evaluation_dict().items(), sort_fcn)

    def evaluation_partition(self):
        """
        Returns the evaluation of the word w as a partition.

        EXAMPLES::

            sage: Word("acdabda").evaluation_partition()
            [3, 2, 1, 1]
            sage: Word([2,1,4,2,3,4,2]).evaluation_partition()
            [3, 2, 1, 1]
        """
        p = sorted(self.evaluation(), reverse=True)
        if 0 in p:
            return Partition(p[:p.index(0)])
        else:
            return Partition(p)

    def inversions(self):
        r"""
        Returns a list of the inversions of self. An inversion is a pair
        (i,j) of non-negative integers i j such that self[i] self[j].

        EXAMPLES::

            sage: Words([1,2,3])([1,2,3,2,2,1]).inversions()
            [[1, 5], [2, 3], [2, 4], [2, 5], [3, 5], [4, 5]]
            sage: Words([3,2,1])([1,2,3,2,2,1]).inversions()
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2]]
            sage: Words('ab')('abbaba').inversions()
            [[1, 3], [1, 5], [2, 3], [2, 5], [4, 5]]
            sage: Words('ba')('abbaba').inversions()
            [[0, 1], [0, 2], [0, 4], [3, 4]]
        """
        return Permutation(list(self._word_content)).inversions()

    def standard_permutation(self):
        r"""
        Returns the standard permutation of the word self on the ordered
        alphabet. It is defined as the permutation with exactly the same
        number of inversions as w. Equivalently, it is the permutation of
        minimal length whose inverse sorts self.

        EXAMPLES::

            sage: w = Word([1,2,3,2,2,1]); w
            word: 123221
            sage: p = w.standard_permutation(); p
            [1, 3, 6, 4, 5, 2]
            sage: v = Word(p.inverse().action(w)); v
            word: 112223
            sage: Permutations(len(w)).filter( \
            ...     lambda q: q.length() <= p.length() and \
            ...               q.inverse().action(w) == list(v) ).list()
            [[1, 3, 6, 4, 5, 2]]

        ::

            sage: w = Words([1,2,3])([1,2,3,2,2,1,2,1]); w
            word: 12322121
            sage: p = w.standard_permutation(); p
            [1, 4, 8, 5, 6, 2, 7, 3]
            sage: Word(p.inverse().action(w))
            word: 11122223

        ::

            sage: w = Words([3,2,1])([1,2,3,2,2,1,2,1]); w
            word: 12322121
            sage: p = w.standard_permutation(); p
            [6, 2, 1, 3, 4, 7, 5, 8]
            sage: Word(p.inverse().action(w))
            word: 32222111

        ::

            sage: w = Words('ab')('abbaba'); w
            word: abbaba
            sage: p = w.standard_permutation(); p
            [1, 4, 5, 2, 6, 3]
            sage: Word(p.inverse().action(w))
            word: aaabbb

        ::

            sage: w = Words('ba')('abbaba'); w
            word: abbaba
            sage: p = w.standard_permutation(); p
            [4, 1, 2, 5, 3, 6]
            sage: Word(p.inverse().action(w))
            word: bbbaaa
        """
        ev = self.evaluation()
        offset = 0
        temp = 0
        for k in range(len(ev)):
            temp = ev[k]
            ev[k] = offset
            offset += temp
        result = []
        for l in self._word_content:
            ev[l] += 1
            result.append(ev[l])
        return Permutation(result)

    def charge(self, check=True):
        """
        EXAMPLES::

            sage: Word([1,1,2,2,3]).charge()
            0
            sage: Word([3,1,1,2,2]).charge()
            1
            sage: Word([2,1,1,2,3]).charge()
            1
            sage: Word([2,1,1,3,2]).charge()
            2
            sage: Word([3,2,1,1,2]).charge()
            2
            sage: Word([2,2,1,1,3]).charge()
            3
            sage: Word([3,2,2,1,1]).charge()
            4

        TESTS::

            sage: Word([3,3,2,1,1]).charge()
            Traceback (most recent call last):
            ...
            ValueError: the evaluation of the word must be a partition
        """
        if check:
            if self.evaluation() not in Partitions():
                raise ValueError, "the evaluation of the word must be a partition"
        w = list(self._word_content)
        res = 0
        while len(w) != 0:
            i = 0
            l = min(w)
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

    def swap(self, i, j=None):
        """
        Returns the word w with entries at positions i and j swapped. By
        default, j = i+1.

        EXAMPLES::

            sage: Word([1,2,3]).swap(0,2)
            word: 321
            sage: Word([1,2,3]).swap(1)
            word: 132
            sage: Words("ba")("abba").swap(1,-1)
            word: aabb
        """
        if j == None:
            j = i+1
        new = list(self)
        (new[i], new[j]) = (new[j], new[i])
        return self.parent()(new)

    def swap_increase(self, i):
        """
        Returns the word with positions i and i+1 exchanged if self[i] >
        self[i+1]. Otherwise, it returns self.

        EXAMPLES::

            sage: w = Word([1,3,2])
            sage: w.swap_increase(1)
            word: 123
            sage: w.swap_increase(0)
            word: 132
            sage: w.swap_increase(0) is w
            True
            sage: Words("ab")("abba").swap_increase(0)
            word: abba
            sage: Words("ba")("abba").swap_increase(0)
            word: baba
        """
        if self.alphabet().rank(self[i]) > self.alphabet().rank(self[i+1]):
            return self.swap(i)
        else:
            return self

    def swap_decrease(self, i):
        """
        Returns the word with positions i and i+1 exchanged if self[i] <
        self[i+1]. Otherwise, it returns self.

        EXAMPLES::

            sage: w = Word([1,3,2])
            sage: w.swap_decrease(0)
            word: 312
            sage: w.swap_decrease(1)
            word: 132
            sage: w.swap_decrease(1) is w
            True
            sage: Words("ab")("abba").swap_decrease(0)
            word: baba
            sage: Words("ba")("abba").swap_decrease(0)
            word: abba
        """
        if self.alphabet().rank(self[i]) < self.alphabet().rank(self[i+1]):
            return self.swap(i)
        else:
            return self

    def lex_less(self, other):
        """
        Returns True if self is lexicographically less than other.

        EXAMPLES::

            sage: w = Word([1,2,3])
            sage: u = Word([1,3,2])
            sage: v = Word([3,2,1])
            sage: w.lex_less(u)
            True
            sage: v.lex_less(w)
            False
            sage: a = Word("abba")
            sage: b = Word("abbb")
            sage: a.lex_less(b)
            True
            sage: b.lex_less(a)
            False
        """
        return self < other

    def lex_greater(self, other):
        """
        Returns True if self is lexicographically greater than other.

        EXAMPLES::

            sage: w = Word([1,2,3])
            sage: u = Word([1,3,2])
            sage: v = Word([3,2,1])
            sage: w.lex_greater(u)
            False
            sage: v.lex_greater(w)
            True
            sage: a = Word("abba")
            sage: b = Word("abbb")
            sage: a.lex_greater(b)
            False
            sage: b.lex_greater(a)
            True
        """
        return self > other

    def degree(self, weights=None):
        """
        Returns the weighted degree of self, where the weighted degree of
        each letter in the ordered alphabet is given by weights, which
        defaults to [1, 2, 3, ...].

        INPUTS: weights - a list or tuple, or a dictionary keyed by the
        letters occurring in self.

        EXAMPLES::

            sage: Word([1,2,3]).degree()
            6
            sage: Word([3,2,1]).degree()
            6
            sage: Word("abba").degree()
            6
            sage: Word("abba").degree([0,2])
            4
            sage: Word("abba").degree([-1,-1])
            -4
            sage: Word("aabba").degree([1,1])
            5
            sage: Words([1,2,4])([1,2,4]).degree()
            6
            sage: Words([1,2,3,4])([1,2,4]).degree()
            7
            sage: Word("aabba").degree({'a':1,'b':2})
            7
            sage: Word([0,1,0]).degree({0:17,1:0})
            34
        """
        deg = 0
        if weights is None:
            for a in self._word_content:
                deg += a+1
        elif isinstance(weights, (list,tuple)):
            for a in self._word_content:
                deg += weights[a]
        elif isinstance(weights, dict):
            for a in self:
                deg += weights[a]
        else:
            raise TypeError, "incorrect type for weights"
        return deg

    def deg_lex_less(self, other, weights=None):
        """
        Returns True if self is degree lexicographically less than other,
        and False otherwise. The weight of each letter in the ordered
        alphabet is given by weights, which defaults to [1, 2, 3, ...].

        EXAMPLES::

            sage: Word([1,2,3]).deg_lex_less(Word([1,3,2]))
            True
            sage: Word([3,2,1]).deg_lex_less(Word([1,2,3]))
            False
            sage: W = Words(range(5))
            sage: W([1,2,4]).deg_lex_less(W([1,3,2]))
            False
            sage: Word("abba").deg_lex_less(Word("abbb"))
            True
            sage: Word("abba").deg_lex_less(Word("baba"))
            True
            sage: Word("abba").deg_lex_less(Word("aaba"))
            False
            sage: Word("abba").deg_lex_less(Word("aaba"),[1,0])
            True
        """
        deg_self = self.degree(weights)
        deg_other = other.degree(weights)
        if deg_self != deg_other:
            return deg_self < deg_other
        return self.lex_less(other)

    def inv_lex_less(self, other):
        """
        Returns True if self is inverse lexicographically less than other.

        EXAMPLES::

            sage: W = Words([1,2,3,4])
            sage: W([1,2,4]).inv_lex_less(W([1,3,2]))
            False
            sage: W([3,2,1]).inv_lex_less(W([1,2,3]))
            True
        """
        if len(self) != len(other):
            return len(self) < len(other)
        return self.reversal() < other.reversal()

    def deg_inv_lex_less(self,other,weights=None):
        """
        Returns True if the word self is degree inverse lexicographically
        less than other.

        EXAMPLES::

            sage: W = Words([1,2,3,4])
            sage: W([1,2,4]).deg_inv_lex_less(W([1,3,2]))
            False
            sage: W([3,2,1]).deg_inv_lex_less(W([1,2,3]))
            True
        """
        d1 = self.degree(weights)
        d2 = other.degree(weights)
        if d1 != d2:
            return d1 < d2
        return self.inv_lex_less(other)

    def rev_lex_less(self,other):
        """
        Returns True if the word self is reverse lexicographically less
        than other.

        EXAMPLES::

            sage: W = Words([1,2,3,4])
            sage: W([1,2,4]).rev_lex_less(W([1,3,2]))
            True
            sage: W([3,2,1]).rev_lex_less(W([1,2,3]))
            False
        """
        if len(self) != len(other):
            return len(self) > len(other)
        return self.reversal() > other.reversal()

    def deg_rev_lex_less(self, other, weights=None):
        """
        Returns True if self is degree reverse lexicographically less than
        other.

        EXAMPLES::

            sage: Word([3,2,1]).deg_rev_lex_less(Word([1,2,3]))
            False
            sage: W = Words([1,2,3,4])
            sage: W([1,2,4]).deg_rev_lex_less(W([1,3,2]))
            False
            sage: W([1,2,3]).deg_rev_lex_less(W([1,2,4]))
            True
        """
        d1 = self.degree(weights)
        d2 = other.degree(weights)
        if d1 != d2:
            return d1 < d2
        return self.rev_lex_less(other)

    def shuffle(self, other, overlap=0):
        """
        Returns the combinatorial class representing the shuffle product
        between words self and other. This consists of all words of length
        len(self)+len(other) that have both self and other as subwords.

        If overlap is non-zero, then the combinatorial class representing
        the shuffle product with overlaps is returned. The calculation of
        the shift in each overlap is done relative to the order of the
        alphabet. For example, "a" shifted by "a" is "b" in the alphabet
        [a, b, c] and 0 shifted by 1 in [0, 1, 2, 3] is 2.

        EXAMPLES::

            sage: W = Words("abcd")
            sage: ab = W("ab")
            sage: cd = W("cd")
            sage: sp = ab.shuffle(cd); sp
            Shuffle product of word: ab and word: cd
            sage: sp.cardinality()
            6
            sage: sp.list()
            [word: abcd, word: acbd, word: acdb, word: cabd, word: cadb, word: cdab]
            sage: W = Words(range(10))
            sage: w = W([0,1])
            sage: u = W([2,3])
            sage: w.shuffle(w)
            Shuffle product of word: 01 and word: 01
            sage: u.shuffle(u)
            Shuffle product of word: 23 and word: 23
            sage: w.shuffle(u)
            Shuffle product of word: 01 and word: 23
            sage: w.shuffle(u,2)
            Overlapping shuffle product of word: 01 and word: 23 with 2 overlaps
        """
        if overlap is True:
            from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping
            return ShuffleProduct_overlapping(self, other)
        elif overlap == 0:
            from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            return ShuffleProduct_w1w2(self, other)
        elif isinstance(overlap, (int,Integer)):
            from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping_r
            return ShuffleProduct_overlapping_r(self, other, overlap)
        raise ValueError, 'overlapping must be True or an integer'

    def shifted_shuffle(self, other):
        """
        Returns the combinatorial class representing the shifted shuffle
        product between words self and other. This is the same as the
        shuffle product of self with the word obtained from other by
        incrementing the values by len(self).

        EXAMPLES::

            sage: W = Words("abcd")
            sage: ab = W("ab")
            sage: cd = W("cd")
            sage: sp = ab.shifted_shuffle(ab); sp
            Shuffle product of word: ab and word: cd
            sage: sp.cardinality()
            6
            sage: sp.list()
            [word: abcd, word: acbd, word: acdb, word: cabd, word: cadb, word: cdab]
        """
        # TODO: type checking: what if the alphabet is too small? error?
        wc = BuildWordContent([x + len(self) for x in other._word_content])
        return self.shuffle(other.parent()(wc))

    def apply_permutation_to_positions(self, permutation):
        r"""
        Return the word obtained by permuting the positions of the letters
        in self.

        EXAMPLES::

            sage: w = Words('abcd')('abcd')
            sage: w.apply_permutation_to_positions([2,1,4,3])
            word: badc
            sage: u = Words('dabc')('abcd')
            sage: u.apply_permutation_to_positions([2,1,4,3])
            word: badc
            sage: w.apply_permutation_to_positions(Permutation([2,1,4,3]))
            word: badc
            sage: w.apply_permutation_to_positions(PermutationGroupElement([2,1,4,3]))
            word: badc
        """
        if not isinstance(permutation, Permutation_class):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.list())
            else:
                permutation = Permutation(permutation)
        return self.parent()(permutation.action(self))

    def apply_permutation_to_letters(self, permutation):
        r"""
        Return the word obtained by applying permutation to the letters of
        the alphabet of self.

        EXAMPLES::

            sage: w = Words('abcd')('abcd')
            sage: p = [2,1,4,3]
            sage: w.apply_permutation_to_letters(p)
            word: badc
            sage: u = Words('dabc')('abcd')
            sage: u.apply_permutation_to_letters(p)
            word: dcba
            sage: w.apply_permutation_to_letters(Permutation(p))
            word: badc
            sage: w.apply_permutation_to_letters(PermutationGroupElement(p))
            word: badc
        """
        if not isinstance(permutation, Permutation_class):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.list())
            else:
                permutation = Permutation(permutation)
        alphabet = self.alphabet().list()
        morphism = dict(zip(alphabet, permutation.action(alphabet)))
        return self.apply_morphism(morphism)
