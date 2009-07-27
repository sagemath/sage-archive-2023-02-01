# coding=utf-8
r"""
Words of all kinds and algorithms

AUTHORS:

- Arnaud Bergeron
- Amy Glen
- Sébastien Labbé
- Franco Saliola

EXAMPLES:

============
Finite words
============

Finite words from python strings, lists and tuples::

    sage: Word("abbabaab")
    word: abbabaab
    sage: Word([0, 1, 1, 0, 1, 0, 0, 1])
    word: 01101001
    sage: Word( ('a', 0, 5, 7, 'b', 9, 8) )
    word: a057b98

Finite words from functions::

    sage: f = lambda n : n%3
    sage: Word(f, length=13)
    word: 0120120120120

Finite words from iterators::

    sage: from itertools import count
    sage: Word(count(), length=10)
    word: 0123456789

::

    sage: Word( iter('abbccdef') )
    word: abbccdef

Finite words from words via concatenation::

    sage: u = Word("abcccabba")
    sage: v = Word([0, 4, 8, 8, 3])
    sage: u * v
    word: abcccabba04883
    sage: v * u
    word: 04883abcccabba
    sage: u + v
    word: abcccabba04883
    sage: u^3 * v^(8/5)
    word: abcccabbaabcccabbaabcccabba04883048

Finite words from infinite words::

    sage: vv = v^Infinity
    sage: vv[10000:10015]
    word: 048830488304883

Finite words in a specific combinatorial class::

    sage: W = Words("ab")
    sage: W
    Words over Ordered Alphabet ['a', 'b']
    sage: W("abbabaab")
    word: abbabaab
    sage: W(["a","b","b","a","b","a","a","b"])
    word: abbabaab
    sage: W( iter('ababab') )
    word: ababab

Finite word as the image under a morphism::

    sage: m = WordMorphism({0:[4,4,5,0],5:[0,5,5],4:[4,0,0,0]})
    sage: m(0)
    word: 4450
    sage: m(0, order=2)
    word: 400040000554450
    sage: m(0, order=3)
    word: 4000445044504450400044504450445044500550...

==============
Infinite words
==============

Periodic infinite words::

    sage: v = Word([0, 4, 8, 8, 3])
    sage: vv = v^Infinity
    sage: vv
    word: 0488304883048830488304883048830488304883...

Infinite words from a function `f:\mathbb{N}\rightarrow A`
over an alphabet `A`::

    sage: Word(lambda n: n%3)
    word: 0120120120120120120120120120120120120120...

::

    sage: def t(n):
    ...       return add(Integer(n).digits(base=2)) % 2
    sage: Word(t, alphabet = [0, 1])
    word: 0110100110010110100101100110100110010110...

or as a one-liner::

    sage: Word(lambda n : add(Integer(n).digits(base=2)) % 2, alphabet = [0, 1])
    word: 0110100110010110100101100110100110010110...

Infinite words from iterators::

    sage: from itertools import count,repeat
    sage: Word( repeat(4) )
    word: 4444444444444444444444444444444444444444...
    sage: Word( count() )
    word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...

Infinite words from morphism

For example, let `A=\{a,b\}` and `\mu : A^* \rightarrow A^*`
be the morphism defined by `a\mapsto ab, b\mapsto ba`::

    sage: mu = WordMorphism('a->ab,b->ba'); print mu
    WordMorphism: a->ab, b->ba
    sage: mu.fixed_point('a')
    Fixed point beginning with 'a' of the morphism WordMorphism: a->ab, b->ba

Infinite words in a specific combinatorial class::

    sage: W = Words("ab"); W
    Words over Ordered Alphabet ['a', 'b']
    sage: f = lambda n : 'a' if n % 2 == 1 else 'b'
    sage: W(f)
    word: babababababababababababababababababababa...

==============
Word functions
==============

::

    sage: w = Word('abaabbba'); w
    word: abaabbba
    sage: w.is_palindrome()
    False
    sage: w.is_lyndon()
    False
    sage: w.number_of_factors()
    28
    sage: w.critical_exponent()
    3

::

    sage: print w.lyndon_factorization()
    (ab, aabbb, a)
    sage: print w.crochemore_factorization()
    (a, b, a, ab, bb, a)

::

    sage: st = w.suffix_tree()
    sage: st
    Implicit Suffix Tree of the word: abaabbba
    sage: st.show(word_labels=True)

::

    sage: T = words.FibonacciWord('ab')
    sage: T.longest_common_prefix(Word('abaabababbbbbb'))
    word: abaababa

As matrix and many other sage objects, words have a parent::

    sage: u = Word('xyxxyxyyy')
    sage: u.parent()
    Words

::

    sage: v = Word('xyxxyxyyy', alphabet='xy')
    sage: v.parent()
    Words over Ordered Alphabet ['x', 'y']

"""
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
from itertools import tee, islice, ifilter, ifilterfalse, imap, izip, \
                      starmap, count, dropwhile, chain, cycle, groupby
from sage.structure.sage_object import SageObject
from sage.sets.set import Set, is_Set
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.misc.latex import latex
from sage.combinat.words.words import Words, Words_all
from sage.combinat.partition import Partition, Partitions
from sage.combinat.combinat import CombinatorialClass
from sage.combinat.permutation import Permutation, Permutation_class
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.words.word_options import word_options
import copy
from word_datatypes import (WordDatatype_str,
                            WordDatatype_list,
                            WordDatatype_tuple)

from word_infinite_datatypes import (
                            WordDatatype_iter_with_caching,
                            WordDatatype_iter,
                            WordDatatype_callable_with_caching,
                            WordDatatype_callable)
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.combinat import CombinatorialObject

# TODO. Word needs to be replaced by Word. Consider renameing
# Word_class to Word and imbedding Word as its __call__ method.

def Word(data=None, alphabet=None, length=None, datatype=None,
        caching=True):
    r"""
    Construct a word.

    INPUT:

    -  ``data`` - (default: None) list, string, tuple, iterator, None
       (shorthand for []), or a callable defined on [0,1,...,length].

    -  ``alphabet`` - any argument accepted by Words

    -  ``length`` - (default: None) This is dependent on the type of data.
       It is ignored for words defined by lists, strings, tuples,
       etc., because they have a naturally defined length.
       For callables, this defines the domain of definition,
       which is assumed to be [0, 1, 2, ..., length-1].
       For iterators: Infinity if you know the iterator will not
       terminate (default); "unknown" if you do not know whether the
       iterator terminates; "finite" if you know that the iterator
       terminates, but do know know the length.

    -  ``datatype`` - (default: None) None, "list", "str", "tuple", "iter",
       "callable". If None, then the function
       tries to guess this from the data.
    -  ``caching`` - (default: True) True or False. Whether to keep a cache
       of the letters computed by an iterator or callable.

    .. note::

       Be careful when defining words using callables and iterators. It
       appears that islice does not pickle correctly causing various errors
       when reloading. Also, most iterators do not support copying and
       should not support pickling by extension.

    EXAMPLES:

    Empty word::

        sage: Word()
        word:

    Word with string::

        sage: Word("abbabaab")
        word: abbabaab

    Word with string constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="str")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="str")
        word: 01101001

    Word with list::

        sage: Word([0,1,1,0,1,0,0,1])
        word: 01101001

    Word with list constructed from other types::

        sage: Word("01101001", datatype="list")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="list")
        word: 01101001

    Word with tuple::

        sage: Word((0,1,1,0,1,0,0,1))
        word: 01101001

    Word with tuple constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="tuple")
        word: 01101001
        sage: Word("01101001", datatype="str")
        word: 01101001

    Word with iterator::

        sage: from itertools import count
        sage: Word(count())
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: Word(iter("abbabaab")) # iterators default to infinite words
        word: abbabaab
        sage: Word(iter("abbabaab"), length="unknown")
        word: abbabaab
        sage: Word(iter("abbabaab"), length="finite")
        word: abbabaab

    Word with function (a 'callable')::

        sage: f = lambda n : add(Integer(n).digits(2)) % 2
        sage: Word(f)
        word: 0110100110010110100101100110100110010110...
        sage: Word(f, length=8)
        word: 01101001

    Word over a string with a parent::

        sage: w = Word("abbabaab", alphabet="abc"); w
        word: abbabaab
        sage: w.parent()
        Words over Ordered Alphabet ['a', 'b', 'c']

    The default parent is the combinatorial class of all words::

        sage: w = Word("abbabaab"); w
        word: abbabaab
        sage: w.parent()
        Words
    """
    # TODO: doctest this part!
    if isinstance(data, Word_class):
        from words import Words
        if data.parent() != Words(alphabet):
            import copy
            data = copy.copy(data)
            data._parent = Words(alphabet)
            data.parent()._check(data)
        return data

    if data is None:
        data = []

    # Guess the datatype if it is not given.
    if datatype is None:
        if isinstance(data, (list, CombinatorialObject)):
            datatype = "list"
        elif isinstance(data, (str)):
            datatype = "str"
        elif isinstance(data, tuple):
            datatype = "tuple"
        elif callable(data):
            datatype = "callable"
        elif hasattr(data,"__iter__"):
            datatype = "iter"
        elif isinstance(data, WordContent):
            # For backwards compatibility (picklejar)
            return _word_from_word_content(data=data, parent=alphabet)
        else:
            raise ValueError, "Cannot guess a datatype; please specify one"
    else:
        # type check the datatypes
        if datatype == "iter" and not hasattr(data, "__iter__"):
            raise ValueError, "Your data is not iterable"
        elif datatype == "callable" and not callable(data):
            raise ValueError, "Your data is not callable"
        elif datatype not in ("list", "tuple", "str",
                                "callable", "iter"):
            raise ValueError, "Unknown datatype"

    # Create the parent object
    from words import Words
    parent = Words() if alphabet is None else Words(alphabet)

    # Construct the word
    if datatype == 'list':
        w = FiniteWord_list(parent=parent,data=data)
    elif datatype == 'str':
        w = FiniteWord_str(parent=parent,data=data)
    elif datatype == 'tuple':
        w = FiniteWord_tuple(parent=parent,data=data)
    elif datatype == 'callable':
        if caching:
            if length is None or length is Infinity:
                cls = InfiniteWord_callable_with_caching
            else:
                cls = FiniteWord_callable_with_caching
        else:
            if length is None or length is Infinity:
                cls = InfiniteWord_callable
            else:
                cls = FiniteWord_callable
        w = cls(parent=parent,callable=data,length=length)
    elif datatype == 'iter':
        if caching:
            if length is None or length is Infinity:
                cls = InfiniteWord_iter_with_caching
            elif length == 'finite':
                cls = FiniteWord_iter_with_caching
            elif length == 'unknown':
                cls = Word_iter_with_caching
            elif length in ZZ and length >= 0:
                cls = FiniteWord_iter_with_caching
            else:
                raise ValueError, "not a correct value for length (%s)" % length
        else:
            if length is None or length is Infinity:
                cls = InfiniteWord_iter
            elif length == 'finite':
                cls = FiniteWord_iter
            elif length == 'unknown':
                cls = Word_iter
            elif length in ZZ and length >= 0:
                cls = FiniteWord_iter
            else:
                raise ValueError, "not a correct value for length (%s)" % length
        w = cls(parent=parent,iter=data,length=length)
    else:
        raise ValueError, "Not known datatype"

    # Do some minimal checking.
    w.parent()._check(w)
    return w

###########################################################################
##### DEPRECATION WARNINGS ################################################
##### Added July 2009 #####################################################
###########################################################################

def is_Word(obj):
    r"""
    Returns True if obj is a word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_Word
        sage: is_Word(33)
        doctest:1: DeprecationWarning: is_Word is deprecated, use isinstance(your_object, Word_all) instead!
        False
        sage: is_Word(Word('abba'))
        True
    """
    from sage.misc.misc import deprecation
    deprecation("is_Word is deprecated, use isinstance(your_object, Word_all) instead!")
    return isinstance(obj, Word_class)

def is_FiniteWord(obj):
    r"""
    Returns True if obj is a finite word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_FiniteWord
        sage: is_FiniteWord(33)
        doctest:1: DeprecationWarning: is_Word is deprecated, use isinstance(your_object, Word_all) instead!
        False
        sage: is_FiniteWord(Word('baab'))
        True
    """
    from sage.misc.misc import deprecation
    deprecation("is_Word is deprecated, use isinstance(your_object, Word_all) instead!")
    if isinstance(obj, Word_class):
        try:
            len(obj)
        except:
            return False
        return True
    else:
        return False

#######################################################################
#                                                                     #
#                    Abstract word classes                            #
#                                                                     #
#######################################################################

class Word_class(SageObject):
    def parent(self):
        r"""
        Returns the parent of self.

        TESTS::

            sage: Word(iter([1,2,3]), length="unknown").parent()
            Words
            sage: Word(range(12)).parent()
            Words
            sage: Word(range(4), alphabet=range(6)).parent()
            Words over Ordered Alphabet [0, 1, 2, 3, 4, 5]
            sage: Word(iter('abac'), alphabet='abc').parent()
            Words over Ordered Alphabet ['a', 'b', 'c']
        """
        return self._parent

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: Word(iter([1,2,3]), length="unknown")._repr_()
            'word: 123'
            sage: Word(xrange(100), length="unknown")._repr_()
            'word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...'
            sage: Word(lambda x:x%3)._repr_()
            'word: 0120120120120120120120120120120120120120...'
        """
        global word_options
        if word_options['old_repr']:
            rep = ""
            if isinstance(self, InfiniteWord_class):
                return "Infinite word over %s" % (str(self.parent().alphabet())[17:])
            elif isinstance(self, FiniteWord_class):
                if word_options['truncate'] and \
                        self.length() > word_options['truncate_length']:
                    return "Finite word of length %s over %s" % (self.length(), str(self.parent().alphabet())[17:])
                else:
                    return word_options['identifier'] + self.string_rep()
            else:
                return "Word over %s" % (str(self.parent().alphabet())[17:])
        return word_options['identifier'] + self.string_rep()

    def string_rep(self):
        r"""
        Returns the raw sequence of letters as a string.

        EXAMPLES::

            sage: Word('abbabaab').string_rep()
            'abbabaab'
            sage: Word([0, 1, 0, 0, 1]).string_rep()
            '01001'
            sage: Word([0,1,10,101]).string_rep()
            '0,1,10,101'
            sage: WordOptions(letter_separator='-')
            sage: Word([0,1,10,101]).string_rep()
            '0-1-10-101'
            sage: WordOptions(letter_separator=',')
        """
        global word_options
        l = word_options['truncate_length']
        letters = list(islice(self, l+1))
        if len(letters) == l+1:
            letters.pop()
            suffix = "..."
        else:
            suffix = ""
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = map(str, letters)
            if all(len(a)==1 for a in letters):
                return ''.join(letters) + suffix
            elif suffix == "...":
                return ls.join(letters) + ls + suffix
            else:
                return ls.join(letters)
        elif word_options['display'] == 'list':
            if suffix == "...":
                return "[%s, %s]" % (str(list(letters))[1:-1], suffix)
            else:
                return str(list(letters))

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word_class()
            sage: w.__iter__()
            Traceback (most recent call last):
            ...
            NotImplementedError: you need to define an iterator in __iter__
        """
        raise NotImplementedError, "you need to define an iterator in __iter__"

    def length(self):
        r"""
        Returns the length of self.

        TESTS::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word(iter('abba'*100), length="unknown")
            sage: w.length() is None
            True
            sage: w = Word(iter('abba'), length="finite")
            sage: w.length()
            4
            sage: w = Word(iter([0,1,1,0,1,0,0,1]*100), length="unknown")
            sage: w.length() is None
            True
            sage: w = Word(iter([0,1,1,0,1,0,0,1]), length="finite")
            sage: w.length()
            8
        """
        return self._len

    # TODO: deprecate the use of __len__ in a later version
    __len__ = length

    def __cmp__(self, other):
        r"""
        Compares two words lexicographically according to Python's built-in
        ordering. Provides for all normal comparison operators.

        NOTE:
            This function will not terminate if self and other are equal
            infinite words!

        EXAMPLES::

            sage: W = Word
            sage: from itertools import count
            sage: W(range(1,10)).__cmp__(W(range(10))) > 0
            True
            sage: W(range(10)).__cmp__(W(range(1,10))) < 0
            True
            sage: W(range(10)).__cmp__(W(range(10))) == 0
            True
            sage: W(range(10)).__cmp__(W(count())) < 0
            True
            sage: W(count()).__cmp__(W(range(10))) > 0
            True
        """
        if not isinstance(other, Word_class):
            return NotImplemented
        self_it, other_it = iter(self), iter(other)
        cmp_fcn = self._parent.cmp_letters
        for (c1, c2) in izip(self_it, other_it):
            r = cmp_fcn(c1,c2)
            if r != 0:
                return r
        else:
            # If self_it is not exhausted, then other_it must be,
            # so other is a proper prefix of self. So self > other;
            # return 1.
            try:
                self_it.next()
                return 1
            except StopIteration:
                # If self_it is exhausted, then we need to check other_it.
                # If other_it is exhausted also, then self == other. Return
                # 0. Otherwise, self is a proper prefix of other.
                # So self < other; return -1.
                try:
                    other_it.next()
                    return -1
                except StopIteration:
                    return 0

    def _longest_common_prefix_iterator(self, other):
        r"""
        Return an iterator of the longest common prefix of self and other.

        INPUT:

        -  ``other`` - word

        OUTPUT:

            iterator

        EXAMPLES::

            sage: f = words.FibonacciWord()
            sage: it = f._longest_common_prefix_iterator(f)
            sage: w = Word(it, length="unknown"); w
            word: 0100101001001010010100100101001001010010...
            sage: w[:6]
            word: 010010
            sage: it = w._longest_common_prefix_iterator(w[:10])
            sage: w = Word(it, length="finite"); w
            word: 0100101001
        """
        for (b,c) in izip(self, other):
            if b == c:
                yield b
            else:
                raise StopIteration
        else:
            raise StopIteration

    def longest_common_prefix(self, other):
        r"""
        Returns the longest common prefix of self and other.

        EXAMPLES::

            sage: f = lambda n : add(Integer(n).digits(2)) % 2
            sage: t = Word(f); t
            word: 0110100110010110100101100110100110010110...
            sage: u = t[:10]; u
            word: 0110100110
            sage: w = t.longest_common_prefix(u); w
            word: 0110100110

        The longest common prefix of two equal infinite words::

            sage: t1 = Word(f); t1
            word: 0110100110010110100101100110100110010110...
            sage: t2 = Word(f); t2
            word: 0110100110010110100101100110100110010110...
            sage: t1.longest_common_prefix(t2)
            word: 0110100110010110100101100110100110010110...
        """
        if isinstance(other, FiniteWord_class):
            return other.longest_common_prefix(self)
        return self._parent(self._longest_common_prefix_iterator(other), length="unknown")

    def is_empty(self):
        r"""
        Returns True if the length of self is zero, and False otherwise.

        EXAMPLES::

            sage: it = iter([])
            sage: Word(it).is_empty()
            True
            sage: it = iter([1,2,3])
            sage: Word(it).is_empty()
            False
            sage: from itertools import count
            sage: Word(count()).is_empty()
            False
        """
        try:
            iter(self).next()
            return False
        except StopIteration:
            return True

    def _to_integer_iterator(self):
        r"""
        Returns an iterator that iterators over the letters of an integer
        representation of self. The first letter occurring in self is
        mapped to zero, and every letter that hasn't yet occurred in the
        word is mapped to the next available integer.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word(count())
            sage: ir = w._to_integer_iterator()
            sage: [ir.next() for _ in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w = Word(iter("abbacabba"))
            sage: ir = w._to_integer_iterator()
            sage: list(ir)
            [0, 1, 1, 0, 2, 0, 1, 1, 0]
        """
        mapping = {}
        next_value = 0
        for letter in self:
            if not(letter in mapping):
                mapping[letter] = next_value
                next_value += 1
            yield mapping[letter]

    def to_integer_word(self):
        r"""
        Returns a word over the integers whose letters are those output by
        self._to_integer_iterator()

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word(count()); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: w.to_integer_word()
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: w = Word(iter("abbacabba"), length="finite"); w
            word: abbacabba
            sage: w.to_integer_word()
            word: 011020110
            sage: w = Word(iter("abbacabba"), length="unknown"); w
            word: abbacabba
            sage: w.to_integer_word()
            word: 011020110
        """
        length = "unknown" if self._len is None else self._len
        return Word(self._to_integer_iterator(), length=length)

    def lex_less(self, other):
        r"""
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

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t.lex_less(t[:10])
            False
            sage: t[:10].lex_less(t)
            True
        """
        return self < other

    def lex_greater(self, other):
        r"""
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

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t[:10].lex_greater(t)
            False
            sage: t.lex_greater(t[:10])
            True
        """
        return self > other

    def apply_morphism(self,morphism):
        r"""
        Returns the word obtained by applying the morphism to self.

        INPUT:

        -  ``morphism`` - Can be an instance of WordMorphism, or
           anything that can be used to construct one.

        EXAMPLES::

            sage: w = Word("ab")
            sage: d = {'a':'ab', 'b':'ba'}
            sage: w.apply_morphism(d)
            word: abba
            sage: w.apply_morphism(WordMorphism(d))
            word: abba

        ::

            sage: w = Word('ababa')
            sage: d = dict(a='ab', b='ba')
            sage: d
            {'a': 'ab', 'b': 'ba'}
            sage: w.apply_morphism(d)
            word: abbaabbaab

        For infinite words::

            sage: t = words.ThueMorseWord([0,1]); t
            Thue-Morse word over Ordered Alphabet [0, 1]
            sage: t.apply_morphism({0:8,1:9})
            word: 8998988998898998988989988998988998898998...
        """
        from sage.combinat.words.morphism import WordMorphism
        if not isinstance(morphism, WordMorphism):
            morphism = WordMorphism(morphism)
        return morphism(self)

    def _delta_iterator(self):
        r"""
        Returns an iterator of the image of self under the delta morphism.
        This is the word composed of the length of consecutive runs of the
        same letter in a given word.

        OUTPUT:

            generator object

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: it=W('22112122')._delta_iterator()
            sage: Word(it)
            word: 22112
            sage: Word(W('555008')._delta_iterator())
            word: 321
            sage: Word(W()._delta_iterator())
            word:

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: it = t._delta_iterator()
            sage: Word(it)
            word: 1211222112112112221122211222112112112221...
        """
        from itertools import groupby
        for letter, run in groupby(self):
            yield len(list(run))

    def delta(self):
        r"""
        Returns the image of self under the delta morphism. This is the
        word composed of the length of consecutive runs of the same letter
        in a given word.

        OUTPUT:

            Word over integers

        EXAMPLES:

        For finite words::

            sage: W = Words('0123456789')
            sage: W('22112122').delta()
            word: 22112
            sage: W('555008').delta()
            word: 321
            sage: W().delta()
            word:
            sage: Word('aabbabaa').delta()
            word: 22112

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t.delta()
            word: 1211222112112112221122211222112112112221...
        """
        return Word(self._delta_iterator())

    def _iterated_right_palindromic_closure_iterator(self, f=None):
        r"""
        Returns an iterator over the iterated (`f`-)palindromic closure of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            iterator -- the iterated (`f`-)palindromic closure of self

        EXAMPLES::

            sage: w = Word('abc')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: abacaba

        ::

            sage: w = Word('aaa')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: aaa

        ::

            sage: w = Word('abbab')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: ababaabababaababa

        An infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: it = t._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: ababaabababaababaabababaababaabababaabab...

        TESTS:

        The empty word::

            sage: w = Word()
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: it.next()
            Traceback (most recent call last):
            ...
            StopIteration

        REFERENCES:

        -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        par = self.parent()
        w = self[:0]
        for letter in self:
            length_before = w.length()
            w = (w*par([letter])).palindromic_closure(f=f)
            length_after = w.length()
            d = length_after - length_before
            for a in w[-d:]:
                yield a

    def iterated_right_palindromic_closure(self, f=None):
        r"""
        Returns the iterated (`f`-)palindromic closure of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            word -- the iterated (`f`-)palindromic closure of self

        EXAMPLES::

            sage: w = Word('abc')
            sage: w.iterated_right_palindromic_closure()
            word: abacaba

            sage: w = Word('aaa')
            sage: w.iterated_right_palindromic_closure()
            word: aaa

            sage: w = Word('abbab')
            sage: w.iterated_right_palindromic_closure()
            word: ababaabababaababa

        An right f-palindromic closure::

            sage: f = WordMorphism('a->b,b->a')
            sage: w = Word('abbab')
            sage: w.iterated_right_palindromic_closure(f=f)
            word: abbaabbaababbaabbaabbaababbaabbaab

        An infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: t.iterated_right_palindromic_closure()
            word: ababaabababaababaabababaababaabababaabab...

        TESTS:

        The empty word::

            sage: w = Word()
            sage: w.iterated_right_palindromic_closure()
            word:

        REFERENCES:

        -   A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        return Word(self._iterated_right_palindromic_closure_iterator(f=f), length='unknown')

    @lazy_attribute
    def _word_content(self):
        r"""
        EXAMPLES::

            sage: w = Word('abaccefa')
            sage: w._word_content
            doctest:...: DeprecationWarning: _word_content is deprecated! try to_integer_word instead! See the documentation for more information
            word: 01022340
        """
        from sage.misc.misc import deprecation
        deprecation("_word_content is deprecated! try to_integer_word instead! See the documentation for more information")
        return self.to_integer_word()

    def alphabet(self):
        r"""
        EXAMPLES::

            sage: w = Word('abaccefa')
            sage: w. alphabet()
            doctest:1: DeprecationWarning: alphabet() is deprecated, use parent().alphabet() instead
            Python objects
            sage: y = Words('456')('64654564')
            sage: y.alphabet()
            Ordered Alphabet ['4', '5', '6']

        """
        from sage.misc.misc import deprecation
        deprecation("alphabet() is deprecated, use parent().alphabet() instead")
        return self.parent().alphabet()

class FiniteWord_class(Word_class):
    def __str__(self):
        r"""
        Returns the full string representation of the word.

        TESTS::

            sage: Word('abc').__str__()
            'word: abc'
            sage: Word([0, 1, 0, 0, 1] * 10).__str__()
            'word: 01001010010100101001010010100101001010010100101001'
            sage: Word([0,1,10,101]).__str__()
            'word: 0,1,10,101'
        """
        global word_options
        letters = list(self)
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = map(str, letters)
            if all(len(a)==1 for a in letters):
                return word_options['identifier'] + ''.join(letters)
            else:
                return word_options['identifier'] + ls.join(letters)
        elif word_options['display'] == 'list':
            return word_options['identifier'] + str(list(letters))

    def coerce(self, other):
        r"""
        Tries to return a pair of words with a common parent; raises an
        exception if this is not possible.

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
            sage: w1.parent() is w2.parent()
            False
        """
        if self.parent() != other.parent():
            try:
                other = self.parent()(other)
                other.parent()._check(other, length=None)
            except:
                try:
                    self = other.parent()(self)
                    self.parent()._check(self, length=None)
                except:
                    raise TypeError, "no coercion rule between %r and %r" % (self.parent(), other.parent())
        return self, other

    def __cmp__(self, other):
        r"""
        Compares two finite words lexicographically according to Python's
        built-in ordering. Provides for all normal comparison operators.

        EXAMPLES::

            sage: W = Word
            sage: W('123').__cmp__(W('1211')) > 0
            True
            sage: W('2111').__cmp__(W('12')) > 0
            True
            sage: W('123').__cmp__(W('123')) == 0
            True
            sage: W('121').__cmp__(W('121')) == 0
            True
            sage: W('123').__cmp__(W('22')) < 0
            True
            sage: W('122').__cmp__(W('32')) < 0
            True
            sage: W([1,1,1]).__cmp__(W([1,1,1,1,1])) < 0
            True
            sage: W([1,1,1,1,1]).__cmp__(W([1,1,1])) > 0
            True
        """
        if isinstance(other, type(self)):
            try:
                self, other = self.coerce(other)
            except TypeError:
                return NotImplemented
            cmp_fcn = self._parent.cmp_letters
            for (c1, c2) in izip(self, other):
                r = cmp_fcn(c1,c2)
                if r != 0:
                    return r
            return self.length() - other.length()
        else:
            # for infinite words, use a super __cmp__.
            return super(FiniteWord_class, self).__cmp__(other)

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
            for s in self._to_integer_iterator():
                res = ((res << 5) + res) + s
            self._hash = res
        return self._hash

    def concatenate(self, other):
        r"""
        Returns the concatenation of self and other.

        INPUT:

        - ``other`` - a word over the same alphabet as self

        EXAMPLES:

        Concatenation may be made using ``+`` or ``*`` operations::

            sage: w = Word('abadafd')
            sage: y = Word([5,3,5,8,7])
            sage: w * y
            word: abadafd53587
            sage: w + y
            word: abadafd53587
            sage: w.concatenate(y)
            word: abadafd53587

        Both words must be defined over the same alphabet::

            sage: z = Word('12223', alphabet = '123')
            sage: z + y
            Traceback (most recent call last):
            ...
            ValueError: 5 not in alphabet!

        Eventually, it should work::

            sage: z = Word('12223', alphabet = '123')
            sage: z + y                   #todo: not implemented
            word: 1222353587
        """
        f = CallableFromListOfWords([self,other])
        length = self.length() + other.length()
        return self._parent(f, length=length, datatype='callable', caching=True)

    __mul__ = concatenate

    __add__ = concatenate

    # TODO: This function is using domain=range(n) for Word but
    # should be a domain=slice(n) # Seb : Feb 23th : I think this is fine now!!
    def __pow__(self, exp):
        r"""
        Return the `exp`-th power of self.

        If `exp` is `\infty`, returns the infinite periodic word of base self.
        Otherwise, `|w|\cdot exp` must be an non-negative integer.

        INPUT:

        -  ``exp``  - an integer, a rational, a float number or plus infinity.

        OUTPUT:

            word -- the exp-th power of self.

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
            ValueError: Power of the word is not defined on the exponent -1: the length of the word (6) times the exponent (-1) must be a positive integer


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

        ...but the length of the word times the exponent must be an integer::

            sage: w = Word(range(6))
            sage: w^(1/4)
            Traceback (most recent call last):
            ...
            ValueError: Power of the word is not defined on the exponent 1/4: the length of the word (6) times the exponent (1/4) must be a positive integer

        You can take infinite power::

            sage: w = Word(range(6)); w
            word: 012345
            sage: u = w^oo; u
            word: 0123450123450123450123450123450123450123...
            sage: u[10000000:20000000]
            word: 4501234501234501234501234501234501234501...
            sage: u[10000000:10000020]
            word: 45012345012345012345
            sage: Word()^oo
            word:
        """
        # powers of the empty word
        if self.is_empty():
            return self

        # infinite power of a non-empty word
        fcn = lambda n: self[n % self.length()]
        if exp is Infinity:
            return self._parent(fcn, length=Infinity)

        #If exp*|self| is not an integer
        length = exp* self.length()
        if length in ZZ and length >= 0:
            return self._parent(fcn, length=length)
        else:
            raise ValueError, "Power of the word is not defined on the \
exponent %s: the length of the word (%s) times the exponent \
(%s) must be a positive integer"  % (exp, self.length(), exp)

    def length(self):
        r"""
        Returns the length of self.

        TESTS::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word(iter('abba'*40), length="finite")
            sage: w._len is None
            True
            sage: w.length()
            160
            sage: w = Word(iter('abba'), length=4)
            sage: w._len
            4
            sage: w.length()
            4
            sage: def f(n):
            ...     return range(2,12,2)[n]
            sage: w = Word(f, length=5)
            sage: w.length()
            5
        """
        if self._len is None:
            self._len = Integer(sum(1 for _ in self))
        return self._len

    def is_empty(self):
        r"""
        Returns True if the length of self is zero, and False otherwise.

        EXAMPLES::

            sage: Word([]).is_empty()
            True
            sage: Word('a').is_empty()
            False
        """
        return self.length()==0

    def to_integer_word(self):
        r"""
        Returns a word defined over the integers [0,1,...,self.length()-1]
        whose letters are in the same relative order in the parent.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word('abbabaab')
            sage: w.to_integer_word()
            word: 01101001
            sage: w = Word(iter("cacao"), length="finite")
            sage: w.to_integer_word()
            word: 10102
            sage: w = Words([3,2,1])([2,3,3,1])
            sage: w.to_integer_word()
            word: 1002
        """
        return Word(self.to_integer_list())

    def to_integer_list(self):
        r"""
        Returns a list of integers from [0,1,...,self.length()-1] in the
        same relative order as the letters in self in the parent.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word('abbabaab')
            sage: w.to_integer_list()
            [0, 1, 1, 0, 1, 0, 0, 1]
            sage: w = Word(iter("cacao"), length="finite")
            sage: w.to_integer_list()
            [1, 0, 1, 0, 2]
            sage: w = Words([3,2,1])([2,3,3,1])
            sage: w.to_integer_list()
            [1, 0, 0, 2]
        """
        cmp_fcn = self._parent.cmp_letters
        ordered_alphabet = sorted(set(self), cmp=cmp_fcn)
        index = dict((b,a) for (a,b) in enumerate(ordered_alphabet))
        return [index[a] for a in self]

    def size_of_alphabet(self):
        r"""
        EXAMPLES::

            sage: w = Word('abaccefa')
            sage: w.size_of_alphabet()
            doctest:1: DeprecationWarning: size_of_alphabet() is deprecated, use parent().size_of_alphabet() instead!
            +Infinity
            sage: y = Words('456')('64654564')
            sage: y.size_of_alphabet()
            3
        """
        from sage.misc.misc import deprecation
        deprecation("size_of_alphabet() is deprecated, use parent().size_of_alphabet() instead!")
        return self.parent().size_of_alphabet()

    ###########################################################################
    ##### DEPRECATION WARNINGS (next 4 functions) #############################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def is_suffix_of(self, other):
        r"""
        Returns True if w is a suffix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: w = W('0123456789')
            sage: y = W('56789')
            sage: y.is_suffix_of(w)
            doctest:1: DeprecationWarning: is_suffix_of is deprecated, use is_suffix instead!
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
        from sage.misc.misc import deprecation
        deprecation("is_suffix_of is deprecated, use is_suffix instead!")
        return self.is_suffix(other)

    def is_proper_suffix_of(self, other):
        r"""
        Returns True if self is a proper suffix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W('23').is_proper_suffix_of(W('123'))
            doctest:1: DeprecationWarning: is_proper_suffix_of is deprecated, use is_proper_suffix instead!
            doctest:...: DeprecationWarning: is_suffix_of is deprecated, use is_suffix instead!
            True
            sage: W('12').is_proper_suffix_of(W('12'))
            False
            sage: W().is_proper_suffix_of(W('123'))
            True
            sage: W('123').is_proper_suffix_of(W('12'))
            False
        """
        from sage.misc.misc import deprecation
        deprecation("is_proper_suffix_of is deprecated, use is_proper_suffix instead!")
        return self.is_proper_suffix(other)

    def is_prefix_of(self, other):
        r"""
        Returns True if self is a prefix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: w = W('0123456789')
            sage: y = W('012345')
            sage: y.is_prefix_of(w)
            doctest:1: DeprecationWarning: is_prefix_of is deprecated, use is_prefix instead!
            True
            sage: w.is_prefix_of(y)
            False
            sage: w.is_prefix_of(W())
            False
            sage: W().is_prefix_of(w)
            True
            sage: W().is_prefix_of(W())
            True
        """
        from sage.misc.misc import deprecation
        deprecation("is_prefix_of is deprecated, use is_prefix instead!")
        return self.is_prefix(other)

    def is_proper_prefix_of(self, other):
        r"""
        Returns True if self is a proper prefix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W('12').is_proper_prefix_of(W('123'))
            doctest:1: DeprecationWarning: is_proper_prefix_of is deprecated, use is_proper_prefix instead!
            doctest:...: DeprecationWarning: is_prefix_of is deprecated, use is_prefix instead!
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
        from sage.misc.misc import deprecation
        deprecation("is_proper_prefix_of is deprecated, use is_proper_prefix instead!")
        return self.is_proper_prefix(other)

    # To fix : do not slice here ! (quite expensive in copy)
    def is_suffix(self, other):
        r"""
        Returns True if w is a suffix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: w = W('0123456789')
            sage: y = W('56789')
            sage: y.is_suffix(w)
            True
            sage: w.is_suffix(y)
            False
            sage: W('579').is_suffix(w)
            False
            sage: W().is_suffix(y)
            True
            sage: w.is_suffix(W())
            False
            sage: W().is_suffix(W())
            True
        """
        return self.is_empty() or self == other[-self.length():]

    def is_proper_suffix(self, other):
        r"""
        Returns True if self is a proper suffix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W('23').is_proper_suffix(W('123'))
            True
            sage: W('12').is_proper_suffix(W('12'))
            False
            sage: W().is_proper_suffix(W('123'))
            True
            sage: W('123').is_proper_suffix(W('12'))
            False
        """
        return self.is_suffix_of(other) and self.length() < other.length()

    def has_suffix(self, other):
        """
        Test whether ``self`` has ``other`` as a suffix.

        .. note::

           Some word datatype classes, like :class:`WordDatatype_str`,
           override this method.

        INPUT::

            - ``other`` - a word, or data describing a word

        OUTPUT::

            - boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("ababa")
            sage: w.has_suffix(u)
            True
            sage: u.has_suffix(w)
            False
            sage: u.has_suffix("ababa")
            True

        ::

            sage: w = Word([0,1,1,0,1,0,0,1,0,1,0,1,0])
            sage: u = Word([0,1,0,1,0])
            sage: w.has_suffix(u)
            True
            sage: u.has_suffix(w)
            False
            sage: u.has_suffix([0,1,0,1,0])
            True

        """
        w = Word(other)
        return w.is_suffix(self)

    def is_prefix(self, other):
        r"""
        Returns True if self is a prefix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: w = W('0123456789')
            sage: y = W('012345')
            sage: y.is_prefix(w)
            True
            sage: w.is_prefix(y)
            False
            sage: w.is_prefix(W())
            False
            sage: W().is_prefix(w)
            True
            sage: W().is_prefix(W())
            True
        """
        return self == other[:self.length()]

    def is_proper_prefix(self, other):
        r"""
        Returns True if self is a proper prefix of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W('12').is_proper_prefix(W('123'))
            True
            sage: W('12').is_proper_prefix(W('12'))
            False
            sage: W().is_proper_prefix(W('123'))
            True
            sage: W('123').is_proper_prefix(W('12'))
            False
            sage: W().is_proper_prefix(W())
            False
        """
        return self.is_prefix_of(other) and self.length() < other.length()

    def has_prefix(self, other):
        r"""
        Test whether ``self`` has ``other`` as a prefix.

        INPUT::

            - ``other`` - a word, or data describing a word

        OUTPUT::

            - boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("abbab")
            sage: w.has_prefix(u)
            True
            sage: u.has_prefix(w)
            False
            sage: u.has_prefix("abbab")
            True

        ::

            sage: w = Word([0,1,1,0,1,0,0,1,0,1,0,1,0])
            sage: u = Word([0,1,1,0,1])
            sage: w.has_prefix(u)
            True
            sage: u.has_prefix(w)
            False
            sage: u.has_prefix([0,1,1,0,1])
            True


        """
        w = Word(other)
        return w.is_prefix(self)

    def reversal(self):
        r"""
        Returns the reversal of self.

        EXAMPLES::

            sage: W = Word
            sage: W('124563').reversal()
            word: 365421
        """
        return self[::-1]

    def prefix_function_table(self):
        r"""
        Returns a vector containing the length of the proper prefix-suffixes
        for all the non-empty prefixes of self.

        EXAMPLES::

            sage: W = Word
            sage: W('121321').prefix_function_table()
            [0, 0, 1, 0, 0, 1]
            sage: W('1241245').prefix_function_table()
            [0, 0, 0, 1, 2, 3, 0]
            sage: W().prefix_function_table()
            []
        """
        k = 0
        res = [0]*self.length()
        for q in xrange(1, self.length()):
            while k > 0 and self[k] != self[q]:
                k = res[k-1]
            if self[k] == self[q]:
                k += 1
            res[q] = k
        return res

    def good_suffix_table(self):
        r"""
        Returns a table of the maximum skip you can do in order not to miss
        a possible occurrence of self in a word.

        This is a part of the Boyer-Moore algorithm to find factors. See [1].

        EXAMPLES::

            sage: W = Word
            sage: W('121321').good_suffix_table()
            [5, 5, 5, 5, 3, 3, 1]
            sage: W('12412').good_suffix_table()
            [3, 3, 3, 3, 3, 1]

        REFERENCES:

        -   [1] R.S. Boyer, J.S. Moore, A fast string searching algorithm,
            Communications of the ACM 20 (1977) 762--772.
        """
        l = self.length()
        p = self.reversal().prefix_function_table()
        res = [l - p[-1]]*(l+1)
        for i in xrange(1, l+1):
            j = l - p[i - 1]
            if res[j] > (i - p[i-1]):
                res[j] = i - p[i-1]
        return res

    def suffix_trie(self):
        r"""
        Returns the suffix trie of self.

        The *suffix trie* of a finite word `w` is a data structure
        representing the factors of `w`. It is a tree whose edges are
        labelled with letters of `w`, and whose leafs correspond to
        suffixes of `w`.

        See sage.combinat.words.suffix_trees.SuffixTrie? for more information.

        EXAMPLES::

            sage: w = Word("cacao")
            sage: w.suffix_trie()                   # not implemented
            Suffix Trie of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.suffix_trie()                   # not implemented
            Suffix Trie of the word: 01011
        """
        from sage.combinat.words.suffix_trees import SuffixTrie
        return SuffixTrie(self)

    def implicit_suffix_tree(self):
        r"""
        Returns the implicit suffix tree of self.

        The *suffix tree* of a word `w` is a compactification of the
        suffix trie for `w`. The compactification removes all nodes that have
        exactly one incoming edge and exactly one outgoing edge. It consists of
        two components: a tree and a word. Thus, instead of labelling the edges
        by factors of `w`, we can labelled them by indices of the occurrence of
        the factors in `w`.

        See sage.combinat.words.suffix_trees.ImplicitSuffixTree? for more information.

        EXAMPLES::

            sage: w = Word("cacao")
            sage: w.implicit_suffix_tree()          # not implemented
            Implicit Suffix Tree of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.implicit_suffix_tree()          # not implemented
            Implicit Suffix Tree of the word: 01011
        """
        from sage.combinat.words.suffix_trees import ImplicitSuffixTree
        return ImplicitSuffixTree(self)

    def suffix_tree(self):
        r"""
        Alias for implicit_suffix_tree().

        EXAMPLES::

            sage: Word('abbabaab').suffix_tree() # not implemented
            Implicit Suffix Tree of the word: abbabaab
        """
        return self.implicit_suffix_tree()

    def number_of_factors(self,n=None):
        r"""
        Counts the number of distinct factors of self.

        INPUT:

        -  ``n`` - an integer, or None.

        OUTPUT:

            If n is an integer, returns the number of distinct factors
            of length n. If n is None, returns the total number of
            distinct factors.

        EXAMPLES::

            sage: w = Word([1,2,1,2,3])
            sage: w.number_of_factors()
            13
            sage: map(w.number_of_factors, range(6))
            [1, 3, 3, 3, 2, 1]

        ::

            sage: Word('1213121').number_of_factors()
            22
            sage: Word('1213121').number_of_factors(1)
            3

        ::

            sage: Word('a'*100).number_of_factors()
            101
            sage: Word('a'*100).number_of_factors(77)
            1

        ::

            sage: Word().number_of_factors()
            1
            sage: Word().number_of_factors(17)
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
        Generates distinct factors of ``self``.

        INPUT:

        -  ``n`` - an integer, or ``None``.

        OUTPUT:

            If ``n`` is an integer, returns an iterator over all distinct
            factors of length ``n``. If ``n`` is ``None``, returns an iterator
            generating all distinct factors.

        EXAMPLES::

            sage: w = Word('1213121')
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

            sage: u = Word([1,2,1,2,3])
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

            sage: type( Word('cacao').factor_iterator() )
            <type 'generator'>
        """
        return self.suffix_tree().factor_iterator(n)

    def factor_set(self):
        r"""
        Returns the set of factors of self.

        EXAMPLES::

            sage: Word('1213121').factor_set()   # random
            Set of elements of <generator object at 0xa8fde6c>
            sage: sorted(  Word([1,2,1,2,3]).factor_set()  )
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]
            sage: sorted(  Word("xx").factor_set()  )
            [word: , word: x, word: xx]
            sage: set( Word().factor_set() )
            set([word: ])
        """
        return Set(set(self.factor_iterator()))

    def commutes_with(self, other):
        r"""
        Returns True if self commutes with other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W('12').commutes_with(W('12'))
            True
            sage: W('12').commutes_with(W('11'))
            False
            sage: W().commutes_with(W('21'))
            True
        """
        return (self * other) == (other * self)

    def conjugate(self, pos):
        r"""
        Returns the conjugate at pos of self.

        pos can be any integer, the distance used is the modulo by the length
        of self.

        EXAMPLES::

            sage: W = Word
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
        pos_mod = pos % self.length()
        return self[pos_mod:] * self[:pos_mod]

    def _conjugates_list(self):
        r"""
        Returns the list of conjugates of self, ordered from the 0-th to the
        (L-1)-st conjugate, where L is the length of self.

        TESTS::

            sage: W = Word
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
        for i in range(1,self.length()):
            S.append(self.conjugate(i))
        return S

    def conjugates(self):
        r"""
        Returns the set of conjugates of self.

        TESTS::

            sage: W = Word
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
        r"""
        Returns the position where self is conjugate with other.
        Returns None if there is no such position.

        EXAMPLES::

            sage: W = Word
            sage: W('12113').conjugate_position(W('31211'))
            1
            sage: W('12131').conjugate_position(W('12113')) is None
            True
            sage: W().conjugate_position(W('123')) is None
            True
        """
        if self.length() != other.length():
            return None
        if self == other:
            return 0
        for l in xrange(1, other.length() - 1):
            other = other.conjugate(1)
            if self == other:
                return l
        return None

    def is_conjugate_with(self, other):
        r"""
        Returns True if self is a conjugate of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
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
        r"""
        Returns True if seq is a cadence of self, and False otherwise.

        A *cadence* is an increasing sequence of indexes that all map to
        the same letter.

        EXAMPLES::

            sage: W = Word
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
        r"""
        Returns the longest common prefix of self and other.

        EXAMPLES::

            sage: W = Word
            sage: w = W('12345')
            sage: y = W('1236777')
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
            sage: Word("11").longest_common_prefix(Word("1"))
            word: 1

        With infinite words::

            sage: t = words.ThueMorseWord('ab')
            sage: u = t[:10]
            sage: u.longest_common_prefix(t)
            word: abbabaabba
            sage: u.longest_common_prefix(u)
            word: abbabaabba
        """
        i=0
        for i,(b,c) in enumerate(izip(self, other)):
            if b != c:
                return self[:i]
        else:
            return self[:i+1]

    def longest_common_suffix(self, other):
        r"""
        Returns the longest common suffix of self and other.

        EXAMPLES::

            sage: W = Word
            sage: w = W('112345678')
            sage: u = W('1115678')
            sage: w.longest_common_suffix(u)
            word: 5678
            sage: u.longest_common_suffix(u)
            word: 1115678
            sage: u.longest_common_suffix(w)
            word: 5678
            sage: w.longest_common_suffix(w)
            word: 112345678
            sage: y = W('549332345')
            sage: w.longest_common_suffix(y)
            word:

        TESTS:

        With the empty word::

            sage: w.longest_common_suffix(W())
            word:
            sage: W().longest_common_suffix(w)
            word:
            sage: W().longest_common_suffix(W())
            word:

        With an infinite word::

            sage: t=words.ThueMorseWord('ab')
            sage: w.longest_common_suffix(t)
            Traceback (most recent call last):
            ...
            TypeError: other must be a finite word
        """
        if not isinstance(other, FiniteWord_class):
            raise TypeError, "other must be a finite word"

        if self.is_empty():
            return self
        if other.is_empty():
            return other

        iter = enumerate(izip(reversed(self), reversed(other)))
        i,(b,c) = iter.next()
        if b != c:
            #In this case, return the empty word
            return self[:0]

        for i,(b,c) in iter:
            if b != c:
                return self[-i:]
        else:
            return self[-i-1:]

    def is_palindrome(self, f=None):
        r"""
        Returns True if self is a palindrome (or a `f`-palindrome), and
        False otherwise.

        Let `f : \Sigma \rightarrow \Sigma` be an involution that extends
        to a morphism on `\Sigma^*`. We say that `w\in\Sigma^*` is a
        *`f`-palindrome* if `w=f(\tilde{w})` [1]. Also called
        *`f`-pseudo-palindrome* [2].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        EXAMPLES::

            sage: W = Word
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

        Some `f`-palindromes::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('aababb').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: Word('abacbacbab').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism({'a':'b','b':'a'})
            sage: Word('aababb').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism({0:[1],1:[0]})
            sage: w = words.ThueMorseWord()[:8]; w
            word: 01101001
            sage: w.is_palindrome(f)
            True

        The word must be in the domain of the involution::

            sage: f = WordMorphism('a->a')
            sage: Word('aababb').is_palindrome(f)
            Traceback (most recent call last):
            ...
            ValueError: b not in alphabet!

        TESTS:

        If the given involution is not an involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('abab').is_palindrome(f)
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution

        ::

            sage: Y = Word
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
            sage: E = WordMorphism('a->b,b->a')
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

        -   [1] S. Labbé, Propriétés combinatoires des `f`-palindromes,
            Mémoire de maîtrise en Mathématiques, Montréal, UQAM, 2008,
            109 pages.
        -   [2] V. Anne, L.Q. Zamboni, I. Zorca, Palindromes and Pseudo-
            Palindromes in Episturmian and Pseudo-Palindromic Infinite Words,
            in : S. Brlek, C. Reutenauer (Eds.), Words 2005, Publications du
            LaCIM, Vol. 36 (2005) 91--100.
        """
        l = self.length()
        if f is None:
            return self[:l/2] == self[l/2 + l%2:].reversal()
        else:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError, "f must be an involution"
            return self[:l/2 + l%2] == f(self[l/2:].reversal())

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def _lps(self, l=None, f=None):
        r"""
        Returns the longest palindromic (or `f`-palindromic) suffix of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).
        -  ``l`` - integer (default: None) the length of the longest palindrome
           suffix of self[:-1], if known.

        OUTPUT:

            word -- If f is None, the longest palindromic suffix of self;
                    otherwise, the longest f-palindromic suffix of self.

        EXAMPLES::

            sage: Word('0111')._lps()
            doctest:1: DeprecationWarning: _lps is deprecated, use lps instead!
            word: 111
            sage: Word('011101')._lps()
            word: 101
            sage: Word('6667')._lps()
            word: 7
            sage: Word('abbabaab')._lps()
            word: baab
            sage: Word()._lps()
            word:
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaab')._lps(f=f)
            word: abbabaab
        """
        from sage.misc.misc import deprecation
        deprecation("_lps is deprecated, use lps instead!")
        return self.lps(l=l, f=f)

    def lps(self, f=None, l=None):
        r"""
        Returns the longest palindromic (or `f`-palindromic) suffix of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).
        -  ``l`` - integer (default: None) the length of the longest palindrome
           suffix of self[:-1], if known.

        OUTPUT:

            word -- If f is None, the longest palindromic suffix of self;
                    otherwise, the longest f-palindromic suffix of self.

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
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaab').lps(f=f)
            word: abbabaab
            sage: W = Word
            sage: w = W('33412321')
            sage: w.lps(l=3)
            word: 12321
            sage: Y = Word
            sage: w = Y('01101001')
            sage: w.lps(l=2)
            word: 1001
            sage: w.lps()
            word: 1001
            sage: w.lps(l=None)
            word: 1001
            sage: Y().lps(l=2)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: v = Word('abbabaab')
            sage: pal = v[:0]
            sage: for i in range(1, v.length()+1):
            ...     pal = v[:i].lps(l=pal.length())
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
            sage: f = WordMorphism('a->b,b->a')
            sage: v = Word('abbabaab')
            sage: pal = v[:0]
            sage: for i in range(1, v.length()+1):
            ...     pal = v[:i].lps(f=f, l=pal.length())
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
        #If the length of the lps of self[:-1] is not known:
        if l == None:
            for i in range(self.length()+1):
                fact = self[i:]
                if fact.is_palindrome(f=f):
                    return fact

        #If l == w[:-1].length(), there is no shortcut
        if self.length() == l + 1:
            return self.lps(f=f)

        #Obtain the letter to the left (g) and to the right (d) of the
        #precedent lps of self
        g = self[-l-2]
        d = self[-1]

        #If the word g*d is a `f`-palindrome, the result follows
        if f is None:
            if g == d:
                return self[-l-2:]
            else:
                #Otherwise, the length of the lps of self is smallest than l+2
                return self[-l-1:].lps()
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            if f(g)[0] == d:
                return self[-l-2:]
            else:
                return self[-l-1:].lps(f=f)

    def palindromic_lacunas_study(self, f=None):
        r"""
        Returns interesting statistics about longest (`f`-)palindromic suffixes
        and lacunas of self (see [1] and [2]).

        Note that a word `w` has at most `|w| + 1` different palindromic factors
        (see [3]).

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

        -  ``list`` - list of the length of the longest palindromic
           suffix (lps) for each non-empty prefix of self;
        -  ``list`` - list of all the lacunas, i.e. positions where there is no
           unioccurrent lps;
        -  ``set`` - set of palindromic factors of self.

        EXAMPLES::

            sage: W = Word
            sage: a,b,c = W('abbabaabbaab').palindromic_lacunas_study()
            sage: a
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]
            sage: b
            [8, 9]
            sage: c          # random order
            set([word: , word: b, word: bab, word: abba, word: bb, word: aa, word: baabbaab, word: baab, word: aba, word: aabbaa, word: a])

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: a,b,c = W('abbabaab').palindromic_lacunas_study(f=f)
            sage: a
            [0, 2, 0, 2, 2, 4, 6, 8]
            sage: b
            [0, 2, 4]
            sage: c           # random order
            set([word: , word: ba, word: baba, word: ab, word: bbabaa, word: abbabaab])
            sage: c == set([W(), W('ba'), W('baba'), W('ab'), W('bbabaa'), W('abbabaab')])
            True

        REFERENCES:

        -   [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic lacunas
            of the Thue-Morse word, Proc. GASCOM 2008 (June 16-20 2008,
            Bibbiena, Arezzo-Italia), 53--67.
        -   [2] A. Blondin-Massé, S. Brlek, A. Frosini, S. Labbé, S. Rinaldi,
            Reconstructing words from a fixed palindromic length sequence,
            Proc. TCS 2008, 5th IFIP International Conference on Theoretical
            Computer Science (September 8-10 2008, Milano, Italia), accepted.
        -   [3] X. Droubay, J. Justin, G. Pirillo, Episturmian words and
            some constructions of de Luca and Rauzy, Theoret. Comput. Sci.
            255 (2001) 539--553.
        """
        #Initialize the results of computations
        palindromes = set()
        lengths_lps = [None] * self.length()
        lacunas = []

        #Initialize the first lps
        pal = self[:0]
        palindromes.add(pal)

        #For all the non-empty prefixes of self,
        for i in xrange(self.length()):

            #Compute its longest `f`-palindromic suffix using the preceding lps (pal)
            pal = self[:i+1].lps(l=pal.length(),f=f)

            lengths_lps[i] = pal.length()

            if pal in palindromes:
                lacunas.append(i)
            else :
                palindromes.add(pal)

        return lengths_lps, lacunas, palindromes

    def lengths_lps(self, f=None):
        r"""
        Returns the list of the length of the longest palindromic
        suffix (lps) for each non-empty prefix of self.

        It corresponds to the function `G_w` defined in [1].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            list -- list of the length of the longest palindromic
                    suffix (lps) for each non-empty prefix of self.

        EXAMPLES::

            sage: Word().lengths_lps()
            []
            sage: Word('a').lengths_lps()
            [1]
            sage: Word('aaa').lengths_lps()
            [1, 2, 3]
            sage: Word('abbabaabbaab').lengths_lps()
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaabbaab').lengths_lps(f)
            [0, 2, 0, 2, 2, 4, 6, 8, 4, 6, 4, 6]

        ::

            sage: f = WordMorphism({5:[8],8:[5]})
            sage: Word([5,8,5,5,8,8,5,5,8,8,5,8,5]).lengths_lps(f)
            [0, 2, 2, 0, 2, 4, 6, 4, 6, 8, 10, 12, 4]

        REFERENCES:

        -   [1] A. Blondin-Massé, S. Brlek, A. Frosini, S. Labbé,
            S. Rinaldi, Reconstructing words from a fixed palindromic length
            sequence, Proc. TCS 2008, 5th IFIP International Conference on
            Theoretical Computer Science (September 8-10 2008, Milano,
            Italia), accepted.
        """
        return self.palindromic_lacunas_study(f=f)[0]

    def lacunas(self, f=None):
        r"""
        Returns the list of all the lacunas of self.

        A *lacuna* is a position in a word where the longest palindromic
        suffix is not unioccurrent (see [1]).

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            list -- list of all the lacunas of self.

        EXAMPLES::

            sage: w = Word([0,1,1,2,3,4,5,1,13,3])
            sage: w.lacunas()
            [7, 9]
            sage: words.ThueMorseWord()[:100].lacunas()
            [8, 9, 24, 25, 32, 33, 34, 35, 36, 37, 38, 39, 96, 97, 98, 99]
            sage: f = WordMorphism({0:[1],1:[0]})
            sage: words.ThueMorseWord()[:50].lacunas(f)
            [0, 2, 4, 12, 16, 17, 18, 19, 48, 49]

        REFERENCES:

        -   [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic lacunas
            of the Thue-Morse word, Proc. GASCOM 2008 (June 16-20 2008,
            Bibbiena, Arezzo-Italia), 53--67.
        """
        return self.palindromic_lacunas_study(f=f)[1]

    def lengths_unioccurrent_lps(self, f=None):
        r"""
        Returns the list of the lengths of the unioccurrent longest palindromic
        suffixes (lps) for each non-empty prefix of self. No unioccurrent lps
        are indicated by None.

        It corresponds to the function `H_w` defined in [1] and [2].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            list -- list of the length of the unioccurrent longest palindromic
                    suffix (lps) for each non-empty prefix of self.
                    No unioccurrent lps are indicated by None.

        EXAMPLES::

            sage: w = Word([0,1,1,2,3,4,5,1,13,3])
            sage: w.lengths_unioccurrent_lps()
            [1, 1, 2, 1, 1, 1, 1, None, 1, None]
            sage: f = words.FibonacciWord()[:20]
            sage: f.lengths_unioccurrent_lps() == f.lengths_lps()
            True
            sage: t = words.ThueMorseWord()
            sage: t[:20].lengths_unioccurrent_lps()
            [1, 1, 2, 4, 3, 3, 2, 4, None, None, 6, 8, 10, 12, 14, 16, 6, 8, 10, 12]
            sage: f = WordMorphism({1:[0],0:[1]})
            sage: t[:15].lengths_unioccurrent_lps(f)
            [None, 2, None, 2, None, 4, 6, 8, 4, 6, 4, 6, None, 4, 6]

        REFERENCES:

        -   [1] A. Blondin-Massé, S. Brlek, S. Labbé, Palindromic lacunas of
            the Thue-Morse word, Proc. GASCOM 2008 (June 16-20 2008, Bibbiena,
            Arezzo-Italia), 53--67.
        -   [2] A. Blondin-Massé, S. Brlek, A. Frosini, S. Labbé, S. Rinaldi,
            Reconstructing words from a fixed palindromic length sequence,
            Proc. TCS 2008, 5th IFIP International Conference on Theoretical
            Computer Science (September 8-10 2008, Milano, Italia), accepted.
        """
        l = self.lengths_lps(f=f)
        for i in self.lacunas(f=f):
            l[i] = None
        return l

    def palindromes(self, f=None):
        r"""
        Returns the set of all palindromic (or `f`-palindromic) factors of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            set -- If f is None, the set of all palindromic factors of self;
                   otherwise, the set of all f-palindromic factors of self.

        EXAMPLES::

            sage: W = Word
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
            sage: f = WordMorphism('a->b,b->a')
            sage: sorted(Word('abbabaab').palindromes(f))
            [word: , word: ab, word: abbabaab, word: ba, word: baba, word: bbabaa]
        """
        return self.palindromic_lacunas_study(f=f)[2]

    def defect(self, f=None):
        r"""
        Returns the defect of self.

        The *defect* of a finite word `w` is given by
        `D(w)=|w|+1-|PAL(w)|`, where `PAL(w)` denotes the set of palindromic
        factors of `w` (including the empty word). See [1].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            integer -- If f is None, the palindromic defect of self;
                       otherwise, the f-palindromic defect of self.

        EXAMPLES::

            sage: words.ThueMorseWord()[:100].defect()
            16
            sage: words.FibonacciWord()[:100].defect()
            0
            sage: W = Word
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

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaabbaababba').defect(f)
            4

        REFERENCES:

        -   [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the Palindromic
            Complexity of Infinite Words, in J. Berstel, J. Karhumaki,
            D. Perrin, Eds, Combinatorics on Words with Applications,
            International Journal of Foundation of Computer Science, Vol. 15,
            No. 2 (2004) 293--306.
        """
        return self.length()+1-len(self.palindromes(f=f))

    def is_full(self, f=None):
        r"""
        Returns True if self has defect 0, and False otherwise.

        A word is *full* if its defect is zero (see [1]).

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            boolean -- If f is None, whether self is full;
                       otherwise, whether self is full of `f`-palindromes.

        EXAMPLES::

            sage: words.ThueMorseWord()[:100].is_full()
            False
            sage: words.FibonacciWord()[:100].is_full()
            True
            sage: Word('000000000000000').is_full()
            True
            sage: Word('011010011001').is_full()
            False
            sage: Word('2194').is_full()
            True
            sage: Word().is_full()
            True

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word().is_full(f)
            True
            sage: w = Word('ab')
            sage: w.is_full()
            True
            sage: w.is_full(f)
            False

        REFERENCES:

        -   [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the Palindromic
            Complexity of Infinite Words, in J. Berstel, J. Karhumaki,
            D. Perrin, Eds, Combinatorics on Words with Applications,
            International Journal of Foundation of Computer Science, Vol. 15,
            No. 2 (2004) 293--306.
        """
        return self.defect(f=f) == 0

    def palindromic_closure(self, side='right', f=None):
        r"""
        Returns the shortest palindrome having self as a prefix
        (or as a suffix if side=='left').

        See [1].

        INPUT:

        -  ``side`` - 'right' or 'left' (default: 'right') the direction of the
           closure
        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            word -- If f is None, the right palindromic closure of self;
                    otherwise, the right f-palindromic closure of self.
                    If side is 'left', the left palindromic closure.

        EXAMPLES::

            sage: W = Word
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

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: w.palindromic_closure(f=f)
            word: abbabaab
            sage: w.palindromic_closure(f=f, side='left')
            word: babaabbaba

        TESTS::

            sage: f = WordMorphism('a->c,c->a')
            sage: w.palindromic_closure(f=f, side='left')
            Traceback (most recent call last):
            ...
            ValueError: b not in alphabet!

        REFERENCES:

        -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        if f is None:
            if side == 'right':
                l = self.lps().length()
                #return self * self[-(l+1)::-1]
                return self * self[:self.length()-l].reversal()
            elif side == 'left':
                l = self.reversal().lps().length()
                return self[:l-1:-1] * self
            else:
                raise ValueError, "side must be either 'left' or 'right' (not %s) " % side
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError, "f must be an involution"
            if side == 'right':
                l = self.lps(f=f).length()
                return self * f(self[-(l+1)::-1])
            elif side == 'left':
                l = self.reversal().lps(f=f).length()
                return f(self[:l-1:-1]) * self
            else:
                raise ValueError, "side must be either 'left' or 'right' (not %s) " % side

    def is_symmetric(self, f=None):
        r"""
        Returns True if self is symmetric (or `f`-symmetric), and
        False otherwise.

        A word is *symmetric* (resp. `f`-*symmetric*) if it is the
        product of two palindromes (resp. `f`-palindromes). See [1] and [2].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        EXAMPLES::

            sage: W = Word
            sage: W('abbabab').is_symmetric()
            True
            sage: W('ababa').is_symmetric()
            True
            sage: W('aababaabba').is_symmetric()
            False
            sage: W('aabbbaababba').is_symmetric()
            False
            sage: f = WordMorphism('a->b,b->a')
            sage: W('aabbbaababba').is_symmetric(f)
            True

        REFERENCES:

        -   [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the Palindromic
            Complexity of Infinite Words, in J. Berstel, J. Karhumaki,
            D. Perrin, Eds, Combinatorics on Words with Applications,
            International Journal of Foundation of Computer Science, Vol. 15,
            No. 2 (2004) 293--306.
        -   [2] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        for i in range(self.length()):
            if self[:i].is_palindrome(f=f) and self[i:].is_palindrome(f=f):
                return True
        return False

    def length_border(self):
        r"""
        Returns the length of the border of self.

        The *border* of a word is the longest word that is both a proper
        prefix and a proper suffix of self.

        EXAMPLES::

            sage: W = Word
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
        if self.is_empty():
            return None
        return self.prefix_function_table()[-1]

    def border(self):
        r"""
        Returns the longest word that is both a proper prefix and a proper
        suffix of self.

        EXAMPLES::

            sage: W = Word
            sage: W('121212').border()
            word: 1212
            sage: W('12321').border()
            word: 1
            sage: W().border() is None
            True
        """
        if self.is_empty():
            return None
        return self[:self.length_border()]

    def minimal_period(self):
        r"""
        Returns the period of self.

        Let `A` be an alphabet. An integer `p\geq 1` is a *period* of a
        word `w=a_1a_2\cdots a_n` where `a_i\in A` if `a_i=a_{i+p}` for
        `i=1,\ldots,n-p`. The smallest period of `w` is called *the*
        period of `w`. See Chapter 1 of [1].

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

        -   [1] M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
            Encyclopedia of Mathematics and its Applications, Cambridge
            University Press, U.K., 2002.
        """
        if self.is_empty():
            return 1
        return self.length()-self.length_border()

    def order(self):
        r"""
        Returns the order of self.

        Let `p(w)` be the period of a word `w`. The positive rational number
        `|w|/p(w)` is the *order* of `w`. See Chapter 8 of [1].

        OUTPUT:

            rational -- the order

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

        -   [1] M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
            Encyclopedia of Mathematics and its Applications, Cambridge
            University Press, U.K., 2002.
        """
        from sage.rings.rational import Rational
        return Rational((self.length(),self.minimal_period()))

    def critical_exponent(self):
        r"""
        Returns the critical exponent of self.

        The *critical exponent* of a word is the supremum of the order of
        all its (finite) factors. See [1].

        .. note::

            The implementation here uses the suffix tree to enumerate all the
            factors. It should be improved.

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

        -   [1] F. Dejean. Sur un théorème de Thue. J. Combinatorial Theory
            Ser. A 13:90–99, 1972.
        """
        return max(map(FiniteWord_class.order, self.factor_iterator()))

    def is_overlap(self):
        r"""
        Returns True if self is an overlap, and False otherwise.

        EXAMPLES::

            sage: W = Word
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
        if self.length() == 0:
            return False
        return self.length_border() > self.length()/2

    def primitive_length(self):
        r"""
        Returns the length of the primitive of self.

        EXAMPLES::

            sage: W = Word
            sage: W('1231').primitive_length()
            4
            sage: W('121212').primitive_length()
            2
        """
        l = lu = self.length()
        if l == 0:
            return 0
        p = self.prefix_function_table()
        while l > 0:
            l = p[l-1]
            if lu % (lu - l) == 0:
                return lu - l

    def is_primitive(self):
        r"""
        Returns True if self is primitive, and False otherwise.

        A finite word `w` is *primitive* if it is not a positive integer
        power of a shorter word.

        EXAMPLES::

            sage: W = Word
            sage: W('1231').is_primitive()
            True
            sage: W('111').is_primitive()
            False
        """
        return self.length() == self.primitive_length()

    def primitive(self):
        r"""
        Returns the primitive of self.

        EXAMPLES::

            sage: W = Word
            sage: W('12312').primitive()
            word: 12312
            sage: W('121212').primitive()
            word: 12
        """
        return self[:self.primitive_length()]

    def exponent(self):
        r"""
        Returns the exponent of self.

        OUTPUT:

            integer -- the exponent

        EXAMPLES::

            sage: W = Word
            sage: W('1231').exponent()
            1
            sage: W('121212').exponent()
            3
            sage: W().exponent()
            0
        """
        if self.length() == 0:
            return 0
        return self.length() / self.primitive_length()

    def is_subword_of(self, other):
        r"""
        Returns True is self is a subword of other, and False otherwise.

        EXAMPLES::

            sage: W = Word
            sage: W().is_subword_of(W('123'))
            True
            sage: W('123').is_subword_of(W('3211333213233321'))
            True
            sage: W('321').is_subword_of(W('11122212112122133111222332'))
            False
        """
        its = iter(self)
        try:
            s = its.next()
            for e in other:
                if s == e:
                    s = its.next()
            else:
                return False
        except StopIteration:
            return True

    def is_lyndon(self):
        r"""
        Returns True if self is a Lyndon word, and False otherwise.

        A *Lyndon word* is a non-empty word that is lexicographically
        smaller than all of its proper suffixes for the given order
        on its alphabet. That is, `w` is a Lyndon word if `w` is non-empty
        and for each factorization `w = uv` (with `u`, `v` both non-empty),
        we have `w < v`.

        Equivalently, `w` is a Lyndon word iff `w` is a non-empty word that is
        lexicographically smaller than all of its proper conjugates for the
        given order on its alphabet.

        See for instance [1].

        EXAMPLES::

            sage: W = Word
            sage: W('123132133').is_lyndon()
            True
            sage: W().is_lyndon()
            True
            sage: W('122112').is_lyndon()
            False

        TESTS:

        A sanity check: ``LyndonWords`` generators Lyndon words, so we
        filter all words of length `n<10` on the alphabet [1,2,3] for
        Lyndon words, and compare with the ``LyndonWords`` generator::

            sage: for n in range(1,10):
            ...       lw1 = [w for w in Words([1,2,3], n) if w.is_lyndon()]
            ...       lw2 = LyndonWords(3,n)
            ...       if set(lw1) != set(lw2): print False

        Filter all words of length 8 on the alphabet [c,a,b] for Lyndon
        words, and compare with the :class:`LyndonWords` generator after
        mapping [a,b,c] to [2,3,1]::

            sage: lw = [w for w in Words('cab', 8) if w.is_lyndon()]
            sage: phi = WordMorphism({'a':2,'b':3,'c':1})
            sage: set(map(phi, lw)) == set(LyndonWords(3,8))
            True

        REFERENCES:

        -   [1] M. Lothaire, Combinatorics On Words, vol. 17 of Encyclopedia
            of Mathematics and its Applications, Addison-Wesley, Reading,
            Massachusetts, 1983.

        """
        if self.is_empty():
            return True
        cmp = self.parent().cmp_letters
        n = self.length()
        i, j = 0, 1
        while j < n:
            c = cmp(self[i], self[j])
            if c == 0:
                # increment i and j
                i += 1
                j += 1
            elif c < 0:
                # reset i, increment j
                i = 0
                j += 1
            else:
                # we found the first word in the lyndon factorization;
                return False
        else:
            return i == 0

    def lyndon_factorization(self):
        r"""
        Returns the Lyndon factorization of self.

        The *Lyndon factorization* of a finite word `w` is the unique
        factorization of `w` as a non-increasing product of Lyndon words,
        i.e., `w = l_1\cdots l_n` where each `l_i` is a Lyndon word and
        `l_1\geq \cdots \geq l_n`. See for instance [1].

        OUTPUT:

            list -- the list of factors obtained

        EXAMPLES::

            sage: Word('010010010001000').lyndon_factorization()
            (01, 001, 001, 0001, 0, 0, 0)
            sage: Words('10')('010010010001000').lyndon_factorization()
            (0, 10010010001000)
            sage: Word('abbababbaababba').lyndon_factorization()
            (abb, ababb, aababb, a)
            sage: Words('ba')('abbababbaababba').lyndon_factorization()
            (a, bbababbaaba, bba)
            sage: Word([1,2,1,3,1,2,1]).lyndon_factorization()
            (1213, 12, 1)

        TESTS::

            sage: Words('01')('').lyndon_factorization()
            ()
            sage: Word('01').lyndon_factorization()
            (01)
            sage: Words('10')('01').lyndon_factorization()
            (0, 1)
            sage: lynfac = Word('abbababbaababba').lyndon_factorization()
            sage: [x.is_lyndon() for x in lynfac]
            [True, True, True, True]
            sage: lynfac = Words('ba')('abbababbaababba').lyndon_factorization()
            sage: [x.is_lyndon() for x in lynfac]
            [True, True, True]
            sage: w = words.ThueMorseWord()[:1000]
            sage: w == prod(w.lyndon_factorization())
            True

        REFERENCES:

        -   [1] J.-P. Duval, Factorizing words over an ordered alphabet,
            J. Algorithms 4 (1983) 363--381.

        -   [2] G. Melancon, Factorizing infinite words using Maple,
            MapleTech journal, vol. 4, no. 1, 1997, pp. 34-42.

        """
        cmp = self.parent().cmp_letters
        # We compute the indexes of the factorization.
        n = self.length()
        k = -1
        F = [0]
        while k < n-1:
            i = k+1
            j = k+2
            while j < n:
                c = cmp(self[i], self[j])
                if c < 0:
                    i = k+1
                    j += 1
                elif c == 0:
                    i += 1
                    j += 1
                else:
                    break
            while k < i:
                F.append(k + j - i + 1)
                k = k + j - i
        return Factorization([self[F[i]:F[i+1]] for i in range(len(F)-1)])

    def inversions(self):
        r"""
        Returns a list of the inversions of self. An inversion is a pair
        (i,j) of non-negative integers i < j such that self[i] > self[j].

        EXAMPLES::

            sage: Word([1,2,3,2,2,1]).inversions()
            [[1, 5], [2, 3], [2, 4], [2, 5], [3, 5], [4, 5]]
            sage: Words([3,2,1])([1,2,3,2,2,1]).inversions()
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2]]
            sage: Word('abbaba').inversions()
            [[1, 3], [1, 5], [2, 3], [2, 5], [4, 5]]
            sage: Words('ba')('abbaba').inversions()
            [[0, 1], [0, 2], [0, 4], [3, 4]]
        """
        inversion_list = []
        cmp_fcn = self._parent.cmp_letters
        for (i1, letter1) in enumerate(self):
            for (i2, letter2) in enumerate(self[i1+1:]):
                if cmp_fcn(letter1, letter2) > 0:
                    inversion_list.append([i1,i1+i2+1])
        return inversion_list

    # TODO: This function should be defined for words of integers, but it
    # naturally is defined over an alphabet with a rank function....
    def degree(self, weights=None):
        r"""
        Returns the weighted degree of self, where the weighted degree of
        each letter in the ordered alphabet is given by weights, which
        defaults to [1, 2, 3, ...].

        INPUTS:

        -  ``weights`` - a list or tuple, or a dictionary keyed by the
           letters occurring in self.

        EXAMPLES::

            sage: Word([1,2,3]).degree()
            6
            sage: Word([3,2,1]).degree()
            6
            sage: Words("ab")("abba").degree()
            6
            sage: Words("ab")("abba").degree([0,2])
            4
            sage: Words("ab")("abba").degree([-1,-1])
            -4
            sage: Words("ab")("aabba").degree([1,1])
            5
            sage: Words([1,2,4])([1,2,4]).degree()
            6
            sage: Word([1,2,4]).degree()
            7
            sage: Word("aabba").degree({'a':1,'b':2})
            7
            sage: Word([0,1,0]).degree({0:17,1:0})
            34
        """
        if isinstance(weights, dict):
            deg = 0
            for a in self:
                deg += weights[a]
            return deg

        if hasattr(self._parent._alphabet, "rank"):
            rank_fcn = self._parent._alphabet.rank
            deg = 0
            if weights is None:
                rank = {}
                for a in self:
                    if a not in rank:
                        rank[a] = rank_fcn(a)
                    deg += rank[a]+1
            elif isinstance(weights, (list,tuple)):
                rank = {}
                for a in self:
                    if a not in rank:
                        rank[a] = rank_fcn(a)
                    deg += weights[rank[a]]
            return deg

        if all(x in ZZ for x in self):
            return sum(self)

        raise TypeError, "degree is not defined for your word"

    def deg_lex_less(self, other, weights=None):
        r"""
        Returns True if self is degree lexicographically less than other,
        and False otherwise. The weight of each letter in the ordered
        alphabet is given by weights, which defaults to [1, 2, 3, ...].

        EXAMPLES::

            sage: WF = Word
            sage: WF([1,2,3]).deg_lex_less(WF([1,3,2]))
            True
            sage: WF([3,2,1]).deg_lex_less(WF([1,2,3]))
            False
            sage: W = Words(range(5))
            sage: W([1,2,4]).deg_lex_less(W([1,3,2]))
            False
            sage: WF("abba").deg_lex_less(WF("abbb"), dict(a=1,b=2))
            True
            sage: WF("abba").deg_lex_less(WF("baba"), dict(a=1,b=2))
            True
            sage: WF("abba").deg_lex_less(WF("aaba"), dict(a=1,b=2))
            False
            sage: WF("abba").deg_lex_less(WF("aaba"), dict(a=1,b=0))
            True
        """
        deg_self = self.degree(weights)
        deg_other = other.degree(weights)
        if deg_self != deg_other:
            return deg_self < deg_other
        return self.lex_less(other)

    def inv_lex_less(self, other):
        r"""
        Returns True if self is inverse lexicographically less than other.

        EXAMPLES::

            sage: W = Word
            sage: W([1,2,4]).inv_lex_less(W([1,3,2]))
            False
            sage: W([3,2,1]).inv_lex_less(W([1,2,3]))
            True
        """
        if self.length() != len(other):
            return self.length() < len(other)
        return self.reversal() < other.reversal()

    def deg_inv_lex_less(self,other,weights=None):
        r"""
        Returns True if the word self is degree inverse lexicographically
        less than other.

        EXAMPLES::

            sage: W = Word
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
        r"""
        Returns True if the word self is reverse
        lexicographically less than other.

        EXAMPLES::

            sage: W = Word
            sage: W([1,2,4]).rev_lex_less(W([1,3,2]))
            True
            sage: W([3,2,1]).rev_lex_less(W([1,2,3]))
            False
        """
        if self.length() != len(other):
            return self.length() > len(other)
        return self.reversal() > other.reversal()

    def deg_rev_lex_less(self, other, weights=None):
        r"""
        Returns True if self is degree reverse
        lexicographically less than other.

        EXAMPLES::

            sage: Word([3,2,1]).deg_rev_lex_less(Word([1,2,3]))
            False
            sage: W = Word
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

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def last_position_table(self):
        r"""
        Returns a dictionary that contains the last position of each letter
        in self.

        EXAMPLES::

            sage: Word('1231232').last_position_table()
            doctest:1: DeprecationWarning: last_position_table is deprecated, use last_position_dict instead!
            {'1': 3, '3': 5, '2': 6}
        """
        from sage.misc.misc import deprecation
        deprecation("last_position_table is deprecated, use last_position_dict instead!")
        lpd = self.last_position_dict()
        if self.parent().size_of_alphabet() in ZZ:
            return [lpd.get(a,-1) for a in self.parent().alphabet()]
        return lpd

    def last_position_dict(self):
        r"""
        Returns a dictionary that contains the last position of each letter
        in self.

        EXAMPLES::

            sage: Word('1231232').last_position_dict()
            {'1': 3, '3': 5, '2': 6}
        """
        d = {}
        for (i, letter) in enumerate(self):
            d[letter] = i
        return d

    def _pos_in(self, other, p):
        r"""
        Returns the position of the first occurrence of self starting at
        position p in other.

        EXAMPLES::

            sage: Word('12')._pos_in(Word('131231'), 2)
            2
            sage: Word('12')._pos_in(Word('131231'), 3) is None
            True
            sage: Word('32')._pos_in(Word('131231'), 0) is None
            True
        """
        lf = self.length()
        lm = len(other)
        if lf == 0 or lm == 0:
            return None
        occ = self.last_position_dict()
        suff = self.good_suffix_table()
        s = p
        while s <= lm - lf:
            for j in xrange(lf-1, -1, -1):
                a = other[s+j]
                if self[j] != a :
                    s += max(suff[j + 1], j - occ.get(a,-1))
                    break
            else:
                return s
        return None

    def first_pos_in(self, other):
        r"""
        Returns the position of the first occurrence of self in other,
        or None if self is not a factor of other.

        EXAMPLES::

            sage: Word('12').first_pos_in(Word('131231'))
            2
            sage: Word('32').first_pos_in(Word('131231')) is None
            True
        """
        return self._pos_in(other, 0)

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def is_factor_of(self, other):
        r"""
        Returns True if self is a factor of other, and False otherwise.

        EXAMPLES::

            sage: u = Word('2113')
            sage: w = Word('123121332131233121132123')
            sage: u.is_factor_of(w)
            doctest:1: DeprecationWarning: is_factor_of is deprecated, use is_factor instead!
            True
            sage: u = Word('321')
            sage: w = Word('1231241231312312312')
            sage: u.is_factor_of(w)
            False
        """
        from sage.misc.misc import deprecation
        deprecation("is_factor_of is deprecated, use is_factor instead!")
        return self.is_factor(other)

    def is_factor(self, other):
        r"""
        Returns True if self is a factor of other, and False otherwise.

        EXAMPLES::

            sage: u = Word('2113')
            sage: w = Word('123121332131233121132123')
            sage: u.is_factor(w)
            True
            sage: u = Word('321')
            sage: w = Word('1231241231312312312')
            sage: u.is_factor(w)
            False
        """
        return self.first_pos_in(other) is not None

    def factor_occurrences_in(self, other):
        r"""
        Returns an iterator over all occurrences (including overlapping ones)
        of self in other in their order of appearance.

        EXAMPLES::

            sage: u = Word('121')
            sage: w = Word('121213211213')
            sage: list(u.factor_occurrences_in(w))
            [0, 2, 8]
        """
        if self.length() == 0:
            raise NotImplementedError, "undefined value"
        p = self._pos_in(other, 0)
        while p is not None:
            yield p
            p = self._pos_in(other, p+1)

    def nb_factor_occurrences_in(self, other):
        r"""
        Returns the number of times self appears as a factor
        in other.

        EXAMPLES::

            sage: W = Word
            sage: W().nb_factor_occurrences_in(W('123'))
            Traceback (most recent call last):
            ...
            NotImplementedError: undefined value
            sage: W('123').nb_factor_occurrences_in(W('112332312313112332121123'))
            4
            sage: W('321').nb_factor_occurrences_in(W('11233231231311233221123'))
            0
        """
        n = 0
        for _ in self.factor_occurrences_in(other):
            n += 1
        return n

    def nb_subword_occurrences_in(self, other):
        r"""
        Returns the number of times self appears in other as a subword.

        EXAMPLES::

            sage: W = Word
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
        ls = self.length()
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
        Returns the return words as a list in the order they appear in the word.

        TESTS::

            sage: W = Word
            sage: W('baccabccbacbca')._return_words_list(W('b'))
            [word: bacca, word: bcc, word: bac]
        """
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

        This is the set of all factors starting by the given factor and ending
        just before the next occurrence of this factor. See [1] and [2].

        EXAMPLES::

            sage: W = Word
            sage: W('21331233213231').return_words(W('2'))
            set([word: 213, word: 21331, word: 233])
            sage: W().return_words(W('213'))
            set([])
            sage: W('121212').return_words(W('1212'))
            set([word: 12])

        REFERENCES:

        -   [1] F. Durand, A characterization of substitutive sequences using
            return words, Discrete Math. 179 (1998) 89-101.
        -   [2] C. Holton, L.Q. Zamboni, Descendants of primitive substitutions,
            Theory Comput. Syst. 32 (1999) 133-157.
        """
        return set(self._return_words_list(fact))

    def complete_return_words(self, fact):
        r"""
        Returns the set of complete return words of fact in self.

        This is the set of all factors starting by the given factor and ending
        just after the next occurrence of this factor. See for instance [1].

        EXAMPLES::

            sage: W = Word
            sage: s = W('21331233213231').complete_return_words(W('2'))
            sage: sorted(s)
            [word: 2132, word: 213312, word: 2332]
            sage: W('').complete_return_words(W('213'))
            set([])
            sage: W('121212').complete_return_words(W('1212'))
            set([word: 121212])

        REFERENCES:

        -   [1] J. Justin, L. Vuillon, Return words in Sturmian and
            episturmian words, Theor. Inform. Appl. 34 (2000) 343--356.
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

    def return_words_derivate(self, fact):
        r"""
        Returns the word generated by mapping a letter to each occurrence of
        the return words for the given factor dropping any dangling prefix and
        suffix. See for instance [1].

        EXAMPLES::

            sage: W = Word
            sage: W('12131221312313122').return_words_derivate(W('1'))
            word: 123242

        REFERENCES:

        -   [1] F. Durand, A characterization of substitutive sequences using
            return words, Discrete Math. 179 (1998) 89--101.
        """
        idx = 0
        tab = {}
        ret = map(lambda w: tab.setdefault(w, len(tab)) + 1, \
                                self._return_words_list(fact))
        return Word(ret)

    def is_quasiperiodic(self):
        r"""
        Returns True if self is quasiperiodic, and False otherwise.

        A finite or infinite word `w` is *quasiperiodic* if it can be
        constructed by concatenations and superpositions of one of its proper
        factors `u`, which is called a *quasiperiod* of `w`.
        See for instance [1], [2], and [3].

        EXAMPLES::

            sage: W = Word
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

        -   [1] A. Apostolico, A. Ehrenfeucht, Efficient detection of
            quasiperiodicities in strings, Theoret. Comput. Sci. 119 (1993)
            247--265.
        -   [2] S. Marcus, Quasiperiodic infinite words, Bull. Eur. Assoc.
            Theor. Comput. Sci. 82 (2004) 170-174.
        -   [3] A. Glen, F. Levé, G. Richomme, Quasiperiodic and Lyndon
            episturmian words, Preprint, 2008, arXiv:0805.0730.
        """
        l = self.length()
        if l <= 1:
           return False
        for i in range(1, l - 1):
            return_lengths = [x.length() for x in self.return_words(self[:i])]
            if return_lengths != []:
               if (max(return_lengths) <= i and self[l-i:l] == self[:i]):
                  return True
        return False

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def _quasiperiods_list(self):
        r"""
        Returns the quasiperiods of self as a list ordered from shortest to
        longest.

        EXAMPLES::

            sage: W = Word
            sage: l = W('abaababaabaababaaba')._quasiperiods_list()
            doctest:1: DeprecationWarning: _quasiperiods_list is deprecated, use quasiperiods instead!
            sage: l
            [word: aba, word: abaaba, word: abaababaaba]
        """
        from sage.misc.misc import deprecation
        deprecation("_quasiperiods_list is deprecated, use quasiperiods instead!")
        return self.quasiperiods()

    def quasiperiods(self):
        r"""
        Returns the quasiperiods of self as a list ordered from shortest to
        longest.

        Let `w` be a finite or infinite word. A *quasiperiod* of `w` is a
        proper factor `u` of `w` such that the occurrences of `u` in `w`
        entirely cover `w`, i.e., every position of `w` falls within some
        occurrence of `u` in `w`. See for instance [1], [2], and [3].

        EXAMPLES::

            sage: W = Word
            sage: W('abaababaabaababaaba').quasiperiods()
            [word: aba, word: abaaba, word: abaababaaba]
            sage: W('abaaba').quasiperiods()
            [word: aba]
            sage: W('abacaba').quasiperiods()
            []

        REFERENCES:

        -   [1] A. Apostolico, A. Ehrenfeucht, Efficient detection of
            quasiperiodicities in strings, Theoret. Comput. Sci. 119 (1993)
            247--265.
        -   [2] S. Marcus, Quasiperiodic infinite words, Bull. Eur. Assoc.
            Theor. Comput. Sci. 82 (2004) 170-174.
        -   [3] A. Glen, F. Levé, G. Richomme, Quasiperiodic and Lyndon
            episturmian words, Preprint, 2008, arXiv:0805.0730.
        """
        l = self.length()
        if l <= 1:
           return []
        Q = []
        for i in range(1, l - 1):
            return_lengths = [x.length() for x in self.return_words(self[:i])]
            if return_lengths != []:
               if (max(return_lengths) <= i and self[l-i:l] == self[:i]):
                  Q.append(self[:i])
        return Q

    def crochemore_factorization(self):
        r"""
        Returns the Crochemore factorization of self as an ordered list of
        factors.

        The *Crochemore factorization* of a finite word `w` is the unique
        factorization: `(x_1, x_2, \ldots, x_n)` of `w` with each `x_i`
        satisfying either:
        C1. `x_i` is a letter that does not appear in `u = x_1\ldots x_{i-1}`;
        C2. `x_i` is the longest prefix of `v = x_i\ldots x_n` that also
        has an occurrence beginning within `u = x_1\ldots x_{i-1}`. See [1].

        .. note::

            This is not a very good implementation, and should be improved.

        EXAMPLES::

            sage: x = Word('abababb')
            sage: x.crochemore_factorization()
            (a, b, abab, b)
            sage: mul(x.crochemore_factorization()) == x
            True
            sage: y = Word('abaababacabba')
            sage: y.crochemore_factorization()
            (a, b, a, aba, ba, c, ab, ba)
            sage: mul(y.crochemore_factorization()) == y
            True
            sage: x = Word([0,1,0,1,0,1,1])
            sage: x.crochemore_factorization()
            (0, 1, 0101, 1)
            sage: mul(x.crochemore_factorization()) == x
            True

        REFERENCES:

        -   [1] M. Crochemore, Recherche linéaire d'un carré dans un mot,
            C. R. Acad. Sci. Paris Sér. I Math. 296 (1983) 14 781--784.
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
                    if xi.first_pos_in(self) < u.length():
                        c.append(xi)
                        break
                    else:
                        xi = xi[:-1]
            u = self[:sum(map(len,c))] # = x_1 ... x_{i-1}
            v = self[sum(map(len,c)):] # = x_i ... x_n
        return c

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def freq(self):
        r"""
        Returns a table of the frequencies of the letters in self.

        OUTPUT:

            dict -- letters associated to their frequency

        EXAMPLES::

            sage: f = Word('1213121').freq()
            doctest:1: DeprecationWarning: freq is deprecated, use evaluation_dict instead!
            sage: f # keys appear in random order
            {'1': 4, '2': 2, '3': 1}

        TESTS::

            sage: f = Word('1213121').freq()
            sage: f['1'] == 4
            True
            sage: f['2'] == 2
            True
            sage: f['3'] == 1
            True
        """
        from sage.misc.misc import deprecation
        deprecation("freq is deprecated, use evaluation_dict instead!")
        return self.evaluation_dict()

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
        for a in self:
            d[a] = d.get(a,0) + 1
        return d

    def evaluation_sparse(self):
        r"""
        Returns a list representing the evaluation of self. The entries of
        the list are two-element lists [a, n], where a is a letter
        occurring in self and n is the number of occurrences of a in self.

        EXAMPLES::

            sage: Word([4,4,2,5,2,1,4,1]).evaluation_sparse()
            [(1, 2), (2, 2), (4, 3), (5, 1)]
            sage: Word("abcaccab").evaluation_sparse()
            [('a', 3), ('c', 3), ('b', 2)]
        """
        return self.evaluation_dict().items()

    def evaluation_partition(self):
        r"""
        Returns the evaluation of the word w as a partition.

        EXAMPLES::

            sage: Word("acdabda").evaluation_partition()
            [3, 2, 1, 1]
            sage: Word([2,1,4,2,3,4,2]).evaluation_partition()
            [3, 2, 1, 1]
        """
        p = sorted(self.evaluation_dict().values(), reverse=True)
        if 0 in p:
            return Partition(p[:p.index(0)])
        else:
            return Partition(p)

    def overlap_partition(self, other, delay=0, p=None):
        r"""
        Returns the partition of the alphabet induced by the equivalence
        relation obtained from the symmetric, reflexive and transitive
        closure of `R_{self,other,delay}\cup p` defined below.

        Let `u = u_0 u_1 \cdots u_{n-1}`, `v = v_0v_1\cdots v_{m-1}` be two
        words on the alphabet `A` where `u_i, v_j \in A` are letters and
        let `d` be an integer. We define a relation
        `R_{u,v,d}\subseteq A \times A` by
        `R_{u,v,d} = \{ (u_k, v_{k-d}) : 0 \leq k < n, 0\leq k-d < m \}`.
        The equivalence relation obtained from `R` is inspired from [1].

        EXAMPLE:

            Let `A = \{\tt{a}, \tt{b}, \tt{c}, \tt{d}, \tt{e}, \tt{f}, \tt{h},
            \tt{l}, \tt{v} \}`,
            `s=\tt{cheval}, t=\tt{abcdef} \in A^*` and `d=3`.
            Then `0 \leq k < 6` and `0\leq k-3 < 6` implies that
            `3\leq k \leq 5`. Then,
            `R_{s,t,d} = \{ (s_3, t_0), (s_4, t_1), (s_5, t_2) \} = \{ (\tt{v},
            \tt{a}), (\tt{a}, \tt{b}), (\tt{l}, \tt{c}) \}`.
            These three couples correspond to the pairs of letters one above
            the other in the following overlap :

                    `\tt{cheval}`
                       `\tt{abcdef}`

            The symmetric, reflexive and transitive closure of `R_{s,t,d}`
            defines the following partition of the alphabet `A`:
            `\{\{\tt{a}, \tt{b}, \tt{v}\}, \{\tt{c}, \tt{l}\}, \{\tt{d}\},
            \{\tt{e}\}, \{\tt{f}\}, \{\tt{h}\}\}`.

        INPUT:

        -  ``other`` - word on the same alphabet as self
        -  ``delay`` - integer
        -  ``p`` - Set (default: None), a partition of the alphabet

        OUTPUT:

        -  ``p`` - Set, a set partition of the alphabet of self and other.

        EXAMPLES:

        The above example::

            sage: W = Words('abcdefhlv')
            sage: cheval = W('cheval')
            sage: abcdef = W('abcdef')
            sage: p = cheval.overlap_partition(abcdef,3); p
            {{'f'}, {'e'}, {'d'}, {'a', 'b', 'v'}, {'c', 'l'}, {'h'}}

        The same example with delay 2 ::

            sage: cheval.overlap_partition(abcdef,2,p)
            {{'f'}, {'a', 'c', 'b', 'e', 'd', 'v', 'l'}, {'h'}}

        ::

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

        ::

            sage: W = Words(range(2))
            sage: w = W([0,1,0,1,0,1]); w
            word: 010101
            sage: w.overlap_partition(w, 0)
            {{0}, {1}}
            sage: w.overlap_partition(w, 1)
            {{0, 1}}

        TESTS::

            sage: empty = Word()
            sage: empty.overlap_partition(empty, 'yo')
            Traceback (most recent call last):
            ...
            TypeError: delay (type given: <type 'str'>) must be an integer
            sage: empty.overlap_partition(empty,2,'yo')
            Traceback (most recent call last):
            ...
            TypeError: p(=yo) is not a Set

        REFERENCES:

        -   [1] S. Labbé, Propriétés combinatoires des `f`-palindromes,
            Mémoire de maîtrise en Mathématiques, Montréal, UQAM, 2008,
            109 pages.
        """
        if not isinstance(delay, (int, Integer)):
            raise TypeError, \
                  "delay (type given: %s) must be an integer"%type(delay)
        elif delay < 0:
            return other.overlap_partition(self, -delay, p)

        alphabet = self.parent().alphabet()

        if p is None:
            if self.parent().size_of_alphabet() is Infinity:
                raise ValueError, 'cannot construct the default set partition of an infinite alphabet'
            else:
                R = [[x] for x in alphabet]
        else:
            if not is_Set(p):
                raise TypeError, "p(=%s) is not a Set" % p
            if Set(alphabet.list()) != reduce(lambda x,y:x.union(y), p):
                raise TypeError, "p(=%s) is not a partition of the alphabet" % p
            R = map(list,p)

        d, n, m = delay, self.length(), other.length()
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

    # TODO: requires a parent with a cmp_letters method
    def standard_permutation(self):
        r"""
        Returns the standard permutation of the word
        self on the ordered alphabet. It is defined as
        the permutation with exactly the same number of
        inversions as w. Equivalently, it is the permutation
        of minimal length whose inverse sorts self.

        EXAMPLES::

            sage: w = Word([1,2,3,2,2,1]); w
            word: 123221
            sage: p = w.standard_permutation(); p
            [1, 3, 6, 4, 5, 2]
            sage: v = Word(p.inverse().action(w)); v
            word: 112223
            sage: Permutations(w.length()).filter( \
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
        ev_dict = self.evaluation_dict()
        ordered_alphabet = sorted(ev_dict, cmp=self.parent().cmp_letters)
        offset = 0
        temp = 0
        for k in ordered_alphabet:
            temp = ev_dict[k]
            ev_dict[k] = offset
            offset += temp
        result = []
        for l in self:
            ev_dict[l] += 1
            result.append(ev_dict[l])
        return Permutation(result)

    def charge(self, check=True):
        r"""
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
            ev_dict = self.evaluation_dict()
            ordered_alphabet = sorted(ev_dict, cmp=self.parent().cmp_letters)
            evaluation = [ev_dict[a] for a in ordered_alphabet]
            if evaluation not in Partitions():
                raise ValueError, "the evaluation of the word must be a partition"
        res = 0
        w = self.to_integer_list()
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

    def BWT(self):
        r"""
        Returns the Burrows-Wheeler Transform (BWT) of self.

        The *Burrows-Wheeler transform* of a finite word `w` is obtained
        from `w` by first listing the conjugates of `w` in lexicographic order
        and then concatenating the final letters of the conjugates in this
        order. See [1].

        EXAMPLES::

            sage: W = Word
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

        -   [1] M. Burrows, D.J. Wheeler, "A block-sorting lossless data
            compression algorithm", HP Lab Technical Report, 1994, available
            at http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.html
        """
        if self.is_empty():
           return self
        conjugates = self._conjugates_list()
        conjugates.sort()
        return self.parent()([x[x.length()-1] for x in conjugates])

    ###########################################################################
    ##### DEPRECATION WARNINGS ################################################
    ##### Added July 2009 #####################################################
    ###########################################################################
    def iterated_palindromic_closure(self, side='right', f=None):
        r"""
        Returns the iterated (`f`-)palindromic closure of self.

        INPUT:

        -  ``side`` - 'right' or 'left' (default: 'right') the direction of the
           closure
        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            word -- If f is None, the right iterated palindromic closure of
            self; otherwise, the right iterated f-palindromic closure
            of self.  If side is 'left', the left palindromic closure.

        EXAMPLES::

            sage: W = Word
            sage: W('123').iterated_palindromic_closure()
            doctest:1: DeprecationWarning: iterated_palindromic_closure is deprecated, use iterated_left_palindromic_closure or iterated_right_palindromic_closure instead!
            word: 1213121
            sage: W('123').iterated_palindromic_closure(side='left')
            word: 3231323
            sage: W('1').iterated_palindromic_closure()
            word: 1
            sage: W().iterated_palindromic_closure()
            word:
            sage: W = Words('ab')
            sage: f = WordMorphism('a->b,b->a')
            sage: W('ab').iterated_palindromic_closure(f=f)
            word: abbaab
            sage: W('ab').iterated_palindromic_closure(f=f, side='left')
            word: abbaab
            sage: W('aab').iterated_palindromic_closure(f=f)
            word: ababbaabab
            sage: W('aab').iterated_palindromic_closure(f=f, side='left')
            word: abbaabbaab

        TESTS::

            sage: W('aab').iterated_palindromic_closure(f=f, side='leftt')
            Traceback (most recent call last):
            ...
            ValueError: side must be either 'left' or 'right' (not leftt)

        If f is not an involution:
            sage: f = WordMorphism('a->b,b->b')
            sage: W('aab').iterated_palindromic_closure(f=f, side='left')
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution

        REFERENCES:

        -   A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        from sage.misc.misc import deprecation
        deprecation("iterated_palindromic_closure is deprecated, "
                   +"use iterated_left_palindromic_closure or "
                   +"iterated_right_palindromic_closure instead!")

        if side == 'right':
            return self.iterated_right_palindromic_closure(f=f)
        elif side == 'left':
            return self.iterated_left_palindromic_closure(f=f)
        else:
            raise ValueError, "side must be either 'left' or 'right' (not %s) " % side

    def iterated_left_palindromic_closure(self, f=None):
        r"""
        Returns the iterated left (`f`-)palindromic closure of self.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

            word -- the left iterated f-palindromic closure of self.

        EXAMPLES::

            sage: W = Word
            sage: W('123').iterated_left_palindromic_closure()
            word: 3231323
            sage: f = WordMorphism('a->b,b->a')
            sage: W('ab').iterated_left_palindromic_closure(f=f)
            word: abbaab
            sage: W('aab').iterated_left_palindromic_closure(f=f)
            word: abbaabbaab

        TESTS:

        If f is not a involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: W('aab').iterated_left_palindromic_closure(f=f)
            Traceback (most recent call last):
            ...
            ValueError: f must be an involution

        REFERENCES:

        -   A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        if f is None:
            return self.reversal().iterated_right_palindromic_closure(f=f)
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            return f(self).reversal().iterated_right_palindromic_closure(f=f)

    def count(self, letter):
        r"""
        Counts the number of occurrences of letter in self.

        EXAMPLES::

            sage: Word('abbabaab').count('a')
            4
        """
        return Integer(sum(1 for a in self if a == letter))

    def is_balanced(self, q=1):
        r"""
        Returns True if self is `q`-balanced, and False otherwise.

        A finite or infinite word `w` is said to be *`q`-balanced* if for
        any two factors `u`, `v` of `w` of the same length, the difference
        between the number of `x`'s in each of `u` and `v` is at most `q`
        for all letters `x` in the alphabet of `w`. A `1`-balanced word is
        simply said to be balanced. See for instance [1] and Chapter 2 of
        [2].

        INPUT:

        -  ``q`` - integer (default 1), the balance level

        OUTPUT:

            boolean -- the result

        EXAMPLES::

            sage: Word('1213121').is_balanced()
            True
            sage: Word('1122').is_balanced()
            False
            sage: Word('121333121').is_balanced()
            False
            sage: Word('121333121').is_balanced(2)
            False
            sage: Word('121333121').is_balanced(3)
            True
            sage: Word('121122121').is_balanced()
            False
            sage: Word('121122121').is_balanced(2)
            True

        TESTS::

            sage: Word('121122121').is_balanced(-1)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Word('121122121').is_balanced(0)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Word('121122121').is_balanced('a')
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer

        REFERENCES:

        -   [1] J. Cassaigne, S. Ferenczi, L.Q. Zamboni, Imbalances in
            Arnoux-Rauzy sequences, Ann. Inst. Fourier (Grenoble) 50 (2000)
            1265--1276.
        -   [2] M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
            Encyclopedia of Mathematics and its Applications, Cambridge
            University Press, U.K., 2002.
        """
        if not isinstance(q, (int, Integer)) or q <= 0:
            raise TypeError, "the balance level must be a positive integer"
        alphabet = set(self)
        for i in xrange(2, self.length()):
            empty_sets = [set() for _ in range(len(alphabet))]
            tab = dict(zip(alphabet, empty_sets))
            for fact in self.factor_iterator(i):
                evaluation_dict = fact.evaluation_dict()
                for a in alphabet:
                    tab[a].add(evaluation_dict.get(a, 0))
            for t in tab.values():
                if len(t) > q+1:
                    return False
        return True

    # TODO.
    # 1. Those three swap functions should use the cmp of python.
    # 2. The actual code should then be copied as is in the Word_over_Alphabet
    # and continue to use the parent cmp
    # 3. Once Word can define Words over alphabet, the examples
    # should be updated appropriately.
    def swap(self, i, j=None):
        r"""
        Returns the word w with entries at positions i and
        j swapped. By default, j = i+1.

        EXAMPLES::

            sage: Word([1,2,3]).swap(0,2)
            word: 321
            sage: Word([1,2,3]).swap(1)
            word: 132
            sage: Word("abba").swap(1,-1)
            word: aabb
        """
        if j == None:
            j = i+1
        new = list(self)
        (new[i], new[j]) = (new[j], new[i])
        return Word(new)

    def swap_increase(self, i):
        r"""
        Returns the word with positions i and i+1 exchanged
        if self[i] > self[i+1]. Otherwise, it returns self.

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
        if self._parent.cmp_letters(self[i], self[i+1]) > 0:
            return self.swap(i)
        else:
            return self

    def swap_decrease(self, i):
        r"""
        Returns the word with positions i and i+1 exchanged
        if self[i] < self[i+1]. Otherwise, it returns self.

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
        if self._parent.cmp_letters(self[i], self[i+1]) < 0:
            return self.swap(i)
        else:
            return self

    def parikh_vector(self, alphabet=None):
        r"""
        Returns the Parikh vector of self, i.e., the vector containing the
        number of occurrences of each letter, given in the order of the
        alphabet.

        See also evaluation_dict.

        INPUT:

        -  ``alphabet`` - (default: None) finite ordered alphabet, if None it
           uses the set of letters in self with the ordering defined by the
           parent

        EXAMPLES::

            sage: Words('ab')().parikh_vector()
            [0, 0]
            sage: Word('aabaa').parikh_vector('abc')
            [4, 1, 0]
            sage: Word('a').parikh_vector('abc')
            [1, 0, 0]
            sage: Word('a').parikh_vector('cab')
            [0, 1, 0]
            sage: Word('a').parikh_vector('bca')
            [0, 0, 1]
            sage: Word().parikh_vector('ab')
            [0, 0]
            sage: Word().parikh_vector('abc')
            [0, 0, 0]
            sage: Word().parikh_vector('abcd')
            [0, 0, 0, 0]

        TESTS::

            sage: Word('aabaa').parikh_vector()
            Traceback (most recent call last):
            ...
            TypeError: the alphabet is infinite; specify a finite alphabet or use evaluation_dict() instead
        """
        if alphabet is None and self._parent.size_of_alphabet() is Infinity:
            raise TypeError, "the alphabet is infinite; specify a finite alphabet or use evaluation_dict() instead"
        if alphabet is None:
            alphabet = self._parent._alphabet
        ev_dict = self.evaluation_dict()
        return [ev_dict.get(a,0) for a in alphabet]

    evaluation = parikh_vector

    def shuffle(self, other, overlap=0):
        r"""
        Returns the combinatorial class representing the shuffle product
        between words self and other. This consists of all words of length
        self.length()+other.length() that have both self and other as
        subwords.

        If overlap is non-zero, then the combinatorial class representing
        the shuffle product with overlaps is returned. The calculation of
        the shift in each overlap is done relative to the order of the
        alphabet. For example, "a" shifted by "a" is "b" in the alphabet
        [a, b, c] and 0 shifted by 1 in [0, 1, 2, 3] is 2.

        INPUT:

        -  ``other`` - finite word
        -  ``overlap`` - (default: 0) integer or True

        OUTPUT:

            Combinatorial class of shuffle product of self and other

        EXAMPLES::

            sage: ab = Word("ab")
            sage: cd = Word("cd")
            sage: sp = ab.shuffle(cd); sp
            Shuffle product of word: ab and word: cd
            sage: sp.cardinality()
            6
            sage: sp.list()
            [word: abcd, word: acbd, word: acdb, word: cabd, word: cadb, word: cdab]
            sage: w = Word([0,1])
            sage: u = Word([2,3])
            sage: w.shuffle(w)
            Shuffle product of word: 01 and word: 01
            sage: u.shuffle(u)
            Shuffle product of word: 23 and word: 23
            sage: w.shuffle(u)
            Shuffle product of word: 01 and word: 23
            sage: w.shuffle(u,2)
            Overlapping shuffle product of word: 01 and word: 23 with 2 overlaps
        """
        if overlap == 0:
            from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            return ShuffleProduct_w1w2(self, other)
        else:
            if any(a not in ZZ for a in self) or any(a not in ZZ for a in other):
                raise ValueError, "for a nonzero overlap, words must contain integers as letters"
            if overlap is True:
                from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping
                return ShuffleProduct_overlapping(self, other)
            elif isinstance(overlap, (int,Integer)):
                from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping_r
                return ShuffleProduct_overlapping_r(self, other, overlap)
            raise ValueError, 'overlapping must be True or an integer'

    def shifted_shuffle(self, other, shift=None):
        r"""
        Returns the combinatorial class representing the shifted shuffle
        product between words self and other. This is the same as the
        shuffle product of self with the word obtained from other by
        incrementing its values (i.e. its letters) by the given shift.

        INPUT:

        -  ``other`` - finite word over the integers
        -  ``shift`` - integer or None (default: None) added to each letter of
           other. When shift is None, it is replaced by self.length()

        OUTPUT:

            Combinatorial class of shifted shuffle products of self and
            other.

        EXAMPLES::

            sage: w = Word([0,1,1])
            sage: sp = w.shifted_shuffle(w); sp
            Shuffle product of word: 011 and word: 344
            sage: sp = w.shifted_shuffle(w, 2); sp
            Shuffle product of word: 011 and word: 233
            sage: sp.cardinality()
            20
            sage: WordOptions(identifier='')
            sage: sp.list()
            [011233, 012133, 012313, 012331, 021133, 021313, 021331, 023113, 023131, 023311, 201133, 201313, 201331, 203113, 203131, 203311, 230113, 230131, 230311, 233011]
            sage: WordOptions(identifier='word: ')
            sage: y = Word('aba')
            sage: y.shifted_shuffle(w,2)
            Traceback (most recent call last):
            ...
            ValueError: for shifted shuffle, words must only contain integers as letters
        """
        if any(a not in ZZ for a in self) or any(a not in ZZ for a in other):
            raise ValueError, "for shifted shuffle, words must only contain integers as letters"
        if shift is None:
            from sage.combinat.words.shuffle_product import ShuffleProduct_shifted
            return ShuffleProduct_shifted(self, other)
        else:
            return self.shuffle(self._parent([x + shift for x in other]))

    def delta_inv(self, W=None, s=None):
        r"""
        Lifts self via the delta operator to obtain a word containing the
        letters in alphabet (default is [0, 1]). The letters used in the
        construction start with s (default is alphabet[0]) and cycle
        through alphabet.

        INPUT:

        -  ``alphabet`` - an iterable
        -  ``s`` - an object in the iterable

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 1, 1]).delta_inv()
            word: 112212
            sage: W([1, 1, 1, 1]).delta_inv(Words('123'))
            word: 1231
            sage: W([2, 2, 1, 1, 2]).delta_inv(s=2)
            word: 22112122
        """
        alphabet = [1, 2] if W is None else W.alphabet()
        cycle_alphabet = cycle(alphabet)
        if self.is_empty():
            return Words(alphabet)()
        if s is None:
            s = cycle_alphabet.next()
        else:
            if s not in alphabet:
                raise ValueError, "starting letter not in alphabet"
            t = cycle_alphabet.next()
            while t != s:
                t = cycle_alphabet.next()
        w = []
        for i in self:
            w.extend([s] * i)
            s = cycle_alphabet.next()
        return Words(alphabet)(w)

    def delta(self):
        r"""
        Returns the image of self under the delta morphism. This is the
        word composed of the length of consecutive runs of the same letter
        in a given word.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('22112122').delta()
            word: 22112
            sage: W('555008').delta()
            word: 321
            sage: W().delta()
            word:
            sage: Word('aabbabaa').delta()
            word: 22112

        """
        if self.is_empty():
            return Words()([])
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
        return Words(range(1,1+max_c))(v)

    # TODO. Decide whether delta_derivate* really need W.alphabet().last()....
    # RENAME: Should "derivate" be derivative?!

    def delta_derivate(self, W=None):
        r"""
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
        r"""
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
        r"""
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
        Applies the phi function to self and returns the result. This is
        the word obtained by taking the first letter of the words obtained
        by iterating delta on self.

        OUTPUT:

            word -- the result of the phi function

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2,2,1,1,2,1,2,2,1,2,2,1,1,2]).phi()
            word: 222222
            sage: W([2,1,2,2,1,2,2,1,2,1]).phi()
            word: 212113
            sage: W().phi()
            word:
            sage: Word([2,1,2,2,1,2,2,1,2,1]).phi()
            word: 212113
            sage: Word([2,3,1,1,2,1,2,3,1,2,2,3,1,2]).phi()
            word: 21215
            sage: Word("aabbabaabaabba").phi()
            word: a22222
            sage: w = Word([2,3,1,1,2,1,2,3,1,2,2,3,1,2])

        REFERENCES:

        -   S. Brlek, A. Ladouceur, A note on differentiable palindromes,
            Theoret. Comput. Sci. 302 (2003) 167--178.
        -   S. Brlek, S. Dulucq, A. Ladouceur, L. Vuillon, Combinatorial
            properties of smooth infinite words, Theoret. Comput. Sci. 352
            (2006) 306--317.
        """
        if self.is_empty():
            return self
        v = [self[0]]
        m = self.delta()
        while m.length() > 1:
            v.append(m[0])
            m = m.delta()
        v.append(m[0])
        return Words()(v)

    def phi_inv(self, W=None):
        r"""
        Apply the inverse of the phi function to self.

        INPUT:

        -  ``self`` - a word over the integers
        -  ``W`` - a parent object of words defined over integers

        OUTPUT:

            word -- the inverse of the phi function

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
        for i in xrange(self.length()-2, -1, -1):
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

        Let `A_k = \{1, \ldots ,k\}`, `k \geq 2`. An infinite word `w` in
        `A_k^\omega` is said to be *smooth* if and only if for all positive
        integers `m`, `\Delta^m(w)` is in `A_k^\omega`, where `\Delta(w)` is
        the word obtained from `w` by composing the length of consecutive
        runs of the same letter in `w`. See for instance [1] and [2].

        INPUT:

        -  ``self`` - must be a word over the integers to get something other
           than False

        OUTPUT:

            boolean -- whether self is a smooth prefix or not

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([1, 1, 2, 2, 1, 2, 1, 1]).is_smooth_prefix()
            True
            sage: W([1, 2, 1, 2, 1, 2]).is_smooth_prefix()
            False

        REFERENCES:

        -   [1] S. Brlek, A. Ladouceur, A note on differentiable palindromes,
            Theoret. Comput. Sci. 302 (2003) 167--178.
        -   [2] S. Brlek, S. Dulucq, A. Ladouceur, L. Vuillon, Combinatorial
            properties of smooth infinite words, Theoret. Comput. Sci. 352
            (2006) 306--317.
        """
        m = self
        W = self.parent()
        while m.length() > 1:
            m = m.delta_derivate_right()
            if not all(W.has_letter(a) for a in m.letters()):
                return False
        return True

    def letters(self):
        r"""
        Return a list of the letters that appear in self, listed in the
        order of first appearance.

        EXAMPLES::

            sage: Word([0,1,1,0,1,0,0,1]).letters()
            [0, 1]
            sage: Word("cacao").letters()
            ['c', 'a', 'o']
        """
        seen, res = {}, []
        for x in self:
            if not seen.has_key(x):
                res.append(x)
                seen[x] = True
        return res

    def standard_factorization(self):
        r"""
        Returns the standard factorization of self.

        The *standard factorization* of a word `w` is the unique
        factorization: `w = uv` where `v` is the longest proper suffix
        of `w` that is a Lyndon word.

        Note that if `w` is a Lyndon word with standard factorization
        `w = uv`, then `u` and `v` are also Lyndon words and `u < v`.

        See for instance [1] and [2].

        OUTPUT:

            list -- the list of factors

        EXAMPLES::

            sage: Words('01')('0010110011').standard_factorization()
            (001011, 0011)
            sage: Words('123')('1223312').standard_factorization()
            (12233, 12)
            sage: Word([3,2,1]).standard_factorization()
            (32, 1)
            sage: Words('123')('').standard_factorization()
            ()

        ::

            sage: w = Word('0010110011',alphabet='01')
            sage: w.standard_factorization()
            (001011, 0011)
            sage: w = Word('0010110011',alphabet='10')
            sage: w.standard_factorization()
            (001011001, 1)
            sage: w = Word('1223312',alphabet='123')
            sage: w.standard_factorization()
            (12233, 12)

        REFERENCES:

        -   [1] K.-T. Chen, R.H. Fox, R.C. Lyndon, Free differential calculus,
            IV. The quotient groups of the lower central series, Ann. of Math.
            68 (1958) 81--95.
        -   [2] J.-P. Duval, Factorizing words over an ordered alphabet,
            J. Algorithms 4 (1983) 363--381.
        """
        if self.is_empty():
            return Factorization([])
        for l in xrange(1, self.length()):
            suff = self[l:]
            if suff.is_lyndon():
                return Factorization([self[:l], suff])
        raise RuntimeError, 'Bug in standard factorization of words'

    def standard_factorization_of_lyndon_factorization(self):
        r"""
        Returns the standard factorization of the Lyndon factorization
        of self.

        OUTPUT:

            list of lists -- the factorization

        EXAMPLES::

            sage: Words('123')('1221131122').standard_factorization_of_lyndon_factorization()
            [(12, 2), (1, 13), (1, 122)]
        """
        return [x.standard_factorization() for x in self.lyndon_factorization()]

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
            sage: Word([1,2,3,4]).apply_permutation_to_positions([3,4,2,1])
            word: 3421
        """
        if not isinstance(permutation, Permutation_class):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.list())
            else:
                permutation = Permutation(permutation)
        return self.parent()(permutation.action(self))

    def apply_permutation_to_letters(self, permutation):
        r"""
        Return the word obtained by applying permutation to
        the letters of the alphabet of self.

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
        alphabet = self.parent().alphabet()
        morphism = dict(zip(alphabet, permutation.action(alphabet)))
        return self.apply_morphism(morphism)

    def colored_vector(self, x=0, y=0, width='default', height=1, cmap='hsv', thickness=1, label=None):
        r"""
        Returns a vector (Graphics object) illustrating self. Each letter
        is represented by a coloured rectangle.

        If the parent of self is a class of words over a finite alphabet,
        then each letter in the alphabet is assigned a unique colour, and
        this colour will be the same every time this method is called. This
        is especially useful when plotting and comparing words defined on
        the same alphabet.

        If the alphabet is infinite, then the letters appearing in the word
        are used as the alphabet.

        INPUT:

        -  ``x`` - (default: 0) bottom left x-coordinate of the vector
        -  ``y`` - (default: 0) bottom left y-coordinate of the vector
        -  ``width``  - (default: 'default') width of the vector. By default,
           the width is the length of self.
        -  ``height`` - (default: 1) height of the vector
        -  ``thickness`` - (default: 1) thickness of the contour
        -  ``cmap`` - (default: 'hsv') color map; for available color map names
            type: ``import matplotlib.cm; matplotlib.cm.datad.keys()``
        -  ``label`` - str (default: None) a label to add on the colored vector.

        OUTPUT:

            Graphics

        EXAMPLES::

            sage: Word(range(20)).colored_vector()
            sage: Word(range(100)).colored_vector(0,0,10,1)
            sage: Words(range(100))(range(10)).colored_vector()
            sage: w = Word('abbabaab')
            sage: w.colored_vector()
            sage: w.colored_vector(cmap='autumn')
            sage: Word(range(20)).colored_vector(label='Rainbow')

        When two words are defined under the same parent, same letters are
        mapped to same colors::

            sage: W = Words(range(20))
            sage: w = W(range(20))
            sage: y = W(range(10,20))
            sage: y.colored_vector(y=1, x=10) + w.colored_vector()

        TESTS:

        The empty word::

            sage: Word().colored_vector()
            sage: Word().colored_vector(label='empty')

        Unknown cmap::

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
        from sage.plot.plot import polygon,line,text

        #The default width of the vector
        if width == 'default':
            width = self.length()

        #The black frame of the vector
        ymax = y + height
        L = [(x,y), (x+width,y), (x+width,ymax), (x,ymax), (x,y)]
        rep = line(L, rgbcolor=(0,0,0), thickness=thickness)

        #The label
        if not label is None:
            hl = height/2.0 # height of the label rectangle
            ymax2 = ymax + hl
            rep += text(str(label), (x+width/2.0, ymax + hl/2.0), rgbcolor=(1,0,0))
            L = [(x,ymax), (x+width,ymax), (x+width,ymax2), (x,ymax2), (x,ymax)]
            rep += line(L, rgbcolor=(0,0,0), thickness=thickness)

        #base : the width of each rectangle
        base = width / float(self.length()) if not self.is_empty() else None

        #A colored rectangle for each letter
        dim = self.parent().size_of_alphabet()
        if dim is Infinity:
            ordered_alphabet = sorted(self.letters(), \
                    cmp=self.parent().cmp_letters)
            dim = float(len(ordered_alphabet))
        else:
            ordered_alphabet = self.parent().alphabet()
            dim = float(self.parent().size_of_alphabet())
        letter_to_integer_dict = dict((a,i) for (i,a) in
                enumerate(ordered_alphabet))
        xp = x
        for a in self:
            i = letter_to_integer_dict[a]
            xq = xp + base
            L = [(xp,y), (xq,y), (xq,ymax), (xp,ymax) ]
            rgbcolor = mpl_cmap( i / dim ) [:3]
            rep += polygon(L, rgbcolor = rgbcolor)
            xp = xq
        rep.axes(False)
        return rep

    def is_square(self):
        r"""
        Returns True if self is a square, and False otherwise.

        EXAMPLES::

            sage: Word([1,0,0,1]).is_square()
            False
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
        if self.length() % 2 != 0:
            return False
        else:
            l = self.length() / 2
            return self[:l] == self[l:]

    def is_square_free(self):
        r"""
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
        l = self.length()
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
        r"""
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
        if self.length() % 3 != 0:
            return False
        l = self.length() / 3
        return self[:l] == self[l:2*l] == self[2*l:]

    def is_cube_free(self):
        r"""
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
        l = self.length()
        if l < 3:
            return True
        suff = self
        for i in xrange(0, l - 3):
            for ll in xrange(3, l-i+1, 3):
                if suff[:ll].is_cube():
                    return False
            suff = suff[1:]
        return True

class InfiniteWord_class(Word_class):
    def length(self):
        r"""
        Returns the length of self.

        EXAMPLES::

            sage: f = lambda n : n % 6
            sage: w = Word(f); w
            word: 0123450123450123450123450123450123450123...
            sage: w.length()
            +Infinity
        """
        return Infinity

#######################################################################
#                                                                     #
#                    Concrete word classes                            #
#                                                                     #
#######################################################################

##### Finite Words #####

class FiniteWord_list(WordDatatype_list, FiniteWord_class):
    r"""
    TESTS::

        sage: w = Word([0,1,1,0])
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_str(WordDatatype_str, FiniteWord_class):
    r"""
    TESTS::

        sage: w = Word('abba')
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_tuple(WordDatatype_tuple, FiniteWord_class):
    r"""
    TESTS::

        sage: w = Word((0,1,1,0))
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_iter_with_caching(WordDatatype_iter_with_caching, FiniteWord_class):
    pass

class FiniteWord_iter(WordDatatype_iter, FiniteWord_class):
    pass

class FiniteWord_callable_with_caching(WordDatatype_callable_with_caching, FiniteWord_class):
    pass

class FiniteWord_callable(WordDatatype_callable, FiniteWord_class):
    pass

##### Infinite Words #####

class InfiniteWord_iter_with_caching(WordDatatype_iter_with_caching, InfiniteWord_class):
    pass

class InfiniteWord_iter(WordDatatype_iter, InfiniteWord_class):
    pass

class InfiniteWord_callable_with_caching(WordDatatype_callable_with_caching, InfiniteWord_class):
    pass

class InfiniteWord_callable(WordDatatype_callable, InfiniteWord_class):
    pass

##### Words of unknown length #####

class Word_iter_with_caching(WordDatatype_iter_with_caching, Word_class):
    pass

class Word_iter(WordDatatype_iter, Word_class):
    pass

#######################################################################

class Factorization(list):
    r"""
    A list subclass having a nicer representation for factorization of words.

    TESTS::

        sage: f = sage.combinat.words.word.Factorization()
        sage: f == loads(dumps(f))
        True
    """
    def __repr__(self):
        r"""
        Returns a string representation of the object.

        TESTS::

            sage: sage.combinat.words.word.Factorization()
            ()
            sage: sage.combinat.words.word.Factorization([Word('ab'), Word('ba')])
            (ab, ba)
        """
        return '(%s)' % ', '.join(w.string_rep() for w in self)

class CallableFromListOfWords(tuple):
    r"""
    A class to create a callable from a list of words. The concatenation of
    a list of words is obtained by creating a word from this callable.
    """
    def __new__(cls, words):
        r"""
        TESTS::

            sage: from sage.combinat.words.word import CallableFromListOfWords
            sage: w,u,x = Word([1,2,3]),Word([4,5]),Word([6,7,8])
            sage: f = CallableFromListOfWords([w,u,x]); f
            (word: 123, word: 45, word: 678)
            sage: f == loads(dumps(f))
            True
        """
        l = []
        for w in words:
            if isinstance(w, WordDatatype_callable) and \
                    isinstance(w._func, CallableFromListOfWords):
                l.extend(w._func)
            else:
                l.append(w)
        return tuple.__new__(cls, l)

    def __call__(self, i):
        r"""
        Returns the character at position i.

        TESTS::

            sage: from sage.combinat.words.word import CallableFromListOfWords
            sage: w,u,x = Word([1,2,3]),Word([4,5]),Word([6,7,8])
            sage: f = CallableFromListOfWords([w,u,x])
            sage: [f(i) for i in range(8)]
            [1, 2, 3, 4, 5, 6, 7, 8]
        """
        j = i
        for c in self:
            if (j - c.length() < 0):
                return c[j]
            j -= c.length()
        raise IndexError, "index (=%s) out of range" % i

#######################################################################
#                                                                     #
#                    Old word classes                                 #
#                                                                     #
#######################################################################

def _word_from_word_content(data, parent):
    r"""
    Create a word with the given parent from a :class:`WordContent` object.

    .. note::

       :class:`WordContent` is deprecated and will be deleted in a future
       version of Sage. This function exists to maintain backwards
       compatibility for unpickling objects saved with older versions of Sage.

    TESTS::

        sage: from sage.combinat.words.word import _word_from_word_content
        sage: import sage.combinat.words.word_content as word_content
        sage: cl = word_content.WordContentFromList([0, 1, 1, 2, 1])
        doctest:...: DeprecationWarning: WordContentFromList is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cl,Words([0,1,2]))
        word: 01121
        sage: _word_from_word_content(cl,Words([2,1,0]))
        word: 21101

    ::

        sage: cf = word_content.WordContentFromFunction(lambda x : x)[:10]
        doctest:...: DeprecationWarning: WordContentFromFunction is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: WordContentFromFunction is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cf, Words())
        word: 0123456789

    ::

        sage: from itertools import count
        sage: ci = word_content.WordContentFromIterator(count(3))[:5]
        doctest:...: DeprecationWarning: WordContentFromIterator is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: WordContentFromIterator is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(ci, Words())
        word: 34567

    ::

        sage: cc = word_content.ConcatenateContent((cl, cf, ci))
        doctest:...: DeprecationWarning: ConcatenateContent is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: ConcatenateContent is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cc, Words())
        word: 01121012345678934567

    """
    from sage.combinat.words import word_content
    # extract data from WordContent object
    length = None
    if isinstance(data, word_content.WordContentFromList):
        unranker = parent.alphabet().unrank
        data = map(unranker, data._list)
        datatype = 'list'
    elif isinstance(data, word_content.WordContentFromIterator):
        data = data._get_it()
        datatype = 'iter'
    elif isinstance(data, word_content.WordContentFromFunction):
        length = data._len
        data = data._func
        datatype = 'callable'
    elif isinstance(data, word_content.ConcatenateContent):
        l = map(_word_from_word_content, data)
        data = CallableFromListOfWords(l)
        datatype = 'callable'
    else:
        raise TypeError, 'data is not an instance of WordContent'
    return Word(data=data, alphabet=parent, datatype=datatype, length=length)

class DeprecatedWordClass(SageObject):
    def __setstate__(self, data):
        from sage.misc.misc import deprecation
        cls = self.__class__
        deprecation('Your word object is saved in an old file format '
            'since %s is deprecated and will be deleted in a future version '
            'of Sage (you can use %s instead). You can re-save your word by '
            'typing "word.save(filename)" to ensure that it will load in '
            'future versions of Sage.''' %
            (cls.__name__, cls.__bases__[1].__name__))
        parent = data['_parent']
        del data['_parent']
        wordcontent = data.get('_word_content', None)
        del data['_word_content']
        w = _word_from_word_content(wordcontent, parent)
        self.__class__ = type(w)
        self.__init__(parent, list(w))
        for key, item in data.iteritems():
            setattr(self, key, item)

class AbstractWord(DeprecatedWordClass, FiniteWord_list):
    pass

class AbstractFiniteWord(DeprecatedWordClass, FiniteWord_list):
    pass

class AbstractInfiniteWord(DeprecatedWordClass, InfiniteWord_class):
    pass

class Word_over_Alphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class Word_over_OrderedAlphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class FiniteWord_over_Alphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class FiniteWord_over_OrderedAlphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class InfiniteWord_over_Alphabet(DeprecatedWordClass, InfiniteWord_class):
    pass

class InfiniteWord_over_OrderedAlphabet(DeprecatedWordClass, InfiniteWord_class):
    pass

