# coding=utf-8
r"""
Finite word

AUTHORS:

- Arnaud Bergeron
- Amy Glen
- Sébastien Labbé
- Franco Saliola
- Julien Leroy (March 2010): reduced_rauzy_graph

EXAMPLES:

=========================
Creation of a finite word
=========================

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
    Words over {'a', 'b'}
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

========================
Functions and algorithms
========================

There are more than 100 functions defined on a finite word. Here are some
of them::

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
    Words over {'x', 'y'}

========================
Factors and Rauzy Graphs
========================

Enumeration of factors, the successive values returned by it.next()
can appear in a different order depending on hardware. Therefore we
mark the three first results of the test random. The important test
is that the iteration stops properly on the fourth call::

    sage: w = Word([4,5,6])^7
    sage: it = w.factor_iterator(4)
    sage: it.next() # random
    word: 6456
    sage: it.next() # random
    word: 5645
    sage: it.next() # random
    word: 4564
    sage: it.next()
    Traceback (most recent call last):
    ...
    StopIteration

The set of factors::

    sage: sorted(w.factor_set(3))
    [word: 456, word: 564, word: 645]
    sage: sorted(w.factor_set(4))
    [word: 4564, word: 5645, word: 6456]
    sage: w.factor_set().cardinality()
    61

Rauzy graphs::

    sage: f = words.FibonacciWord()[:30]
    sage: f.rauzy_graph(4)
    Looped digraph on 5 vertices
    sage: f.reduced_rauzy_graph(4)
    Looped multi-digraph on 2 vertices

Left-special and bispecial factors::

    sage: f.number_of_left_special_factors(7)
    1
    sage: f.bispecial_factors()
    [word: , word: 0, word: 010, word: 010010, word: 01001010010]
"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                     2008 Amy Glen <amy.glen@gmail.com>,
#                     2008-2012 Sébastien Labbé <slabqc@gmail.com>,
#                     2008-2010 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import islice, izip, cycle
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method
from sage.combinat.words.word_options import word_options
from sage.rings.all import Integer, Infinity, ZZ
from sage.sets.set import Set

class FiniteWord_class(Word_class):
    def __str__(self):
        r"""
        Returns the full (not truncated) string representation of the word
        without identifier.

        TESTS::

            sage: Word('abc').__str__()
            'abc'
            sage: Word([0, 1, 0, 0, 1] * 10).__str__()
            '01001010010100101001010010100101001010010100101001'
            sage: Word([0,1,10,101]).__str__()
            '0,1,10,101'

        Insertion in a str::

            sage: w = Word(range(5))
            sage: "Let's insert the word w = %s in this string." % w
            "Let's insert the word w = 01234 in this string."

        Using LatexExpr::

            sage: from sage.misc.latex import LatexExpr
            sage: LatexExpr(w)
            01234

        With the print statement::

            sage: print w
            01234

        No truncation is done for finite words::

            sage: w = Word([i % 5 for i in range(60)])
            sage: print w
            012340123401234012340123401234012340123401234012340123401234
        """
        global word_options
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = map(str, self)
            if all(len(a)==1 for a in letters):
                return ''.join(letters)
            else:
                return ls.join(letters)
        elif word_options['display'] == 'list':
            return str(list(self))

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TESTS::

            sage: Word(range(10))._repr_()
            'word: 0123456789'
            sage: Word(range(100))._repr_()
            'word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...'
        """
        global word_options
        if word_options['old_repr']:
            if word_options['truncate'] and \
                    self.length() > word_options['truncate_length']:
                return "Finite word of length %s over %s" % (self.length(), str(self.parent().alphabet())[17:])
        return word_options['identifier'] + self.string_rep()

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
            except Exception:
                try:
                    self = other.parent()(self)
                    self.parent()._check(self, length=None)
                except Exception:
                    raise TypeError("no coercion rule between %r and %r" % (self.parent(), other.parent()))
        return self, other

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

        TESTS:

        The empty word is not considered by concatenation::

            sage: type(Word([]) * Word('abcd'))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word())
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word([]))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word(()))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word([1,2,3]) * Word(''))
            <class 'sage.combinat.words.word.FiniteWord_list'>
        """
        if self.is_empty():
            return other
        if isinstance(other, Word_class) and other.is_empty():
            return self
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
            raise ValueError("Power of the word is not defined on the \
exponent %s: the length of the word (%s) times the exponent \
(%s) must be a positive integer"  % (exp, self.length(), exp))

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

    def schuetzenberger_involution(self, n = None):
        """
        Returns the Schuetzenberger involution of the word self, which is obtained
        by reverting the word and then complementing all letters within the
        underlying ordered alphabet. If `n` is specified, the underlying
        alphabet is assumed to be `[1,2,\ldots,n]`. If no alphabet is specified,
        `n` is the maximal letter appearing in self.

       INPUT:

        - ``self`` -- a word
        - ``n``    -- an integer specifying the maximal letter in the alphabet (optional)

        OUTPUT:

        - a word, the Schuetzenberger involution of self

        EXAMPLES::

            sage: w = Word([9,7,4,1,6,2,3])
            sage: v = w.schuetzenberger_involution(); v
            word: 7849631
            sage: v.parent()
            Words

            sage: w = Word([1,2,3],alphabet=[1,2,3,4,5])
            sage: v = w.schuetzenberger_involution();v
            word: 345
            sage: v.parent()
            Words over {1, 2, 3, 4, 5}

            sage: w = Word([1,2,3])
            sage: v = w.schuetzenberger_involution(n=5);v
            word: 345
            sage: v.parent()
            Words

            sage: w = Word([11,32,69,2,53,1,2,3,18,41])
            sage: w.schuetzenberger_involution()
            word: 29,52,67,68,69,17,68,1,38,59

            sage: w = Word([],alphabet=[1,2,3,4,5])
            sage: w.schuetzenberger_involution()
            word:

            sage: w = Word([])
            sage: w.schuetzenberger_involution()
            word:
        """
        if self.length() == 0:
            return self
        r = self.reversal()
        w = list(r)
        if n is None:
            alphsize = self.parent().size_of_alphabet()
            if not alphsize == +Infinity:
                n = max(self.parent().alphabet())
            elif r.length()>0:
                n = max(w)
        for k in range(r.length()):
            w[k] = n+1 - w[k]
        return self.parent()(w)

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

    def is_finite(self):
        r"""
        Returns True.

        EXAMPLES::

            sage: Word([]).is_finite()
            True
            sage: Word('a').is_finite()
            True
        """
        return True

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
        from sage.combinat.words.word import Word
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

    # To fix : do not slice here ! (quite expensive in copy)
    def is_suffix(self, other):
        r"""
        Returns True if w is a suffix of other, and False otherwise.

        EXAMPLES::

            sage: w = Word('0123456789')
            sage: y = Word('56789')
            sage: y.is_suffix(w)
            True
            sage: w.is_suffix(y)
            False
            sage: Word('579').is_suffix(w)
            False
            sage: Word().is_suffix(y)
            True
            sage: w.is_suffix(Word())
            False
            sage: Word().is_suffix(Word())
            True
        """
        return self.is_empty() or self == other[-self.length():]

    def is_proper_suffix(self, other):
        r"""
        Returns True if self is a proper suffix of other, and False otherwise.

        EXAMPLES::

            sage: Word('23').is_proper_suffix(Word('123'))
            True
            sage: Word('12').is_proper_suffix(Word('12'))
            False
            sage: Word().is_proper_suffix(Word('123'))
            True
            sage: Word('123').is_proper_suffix(Word('12'))
            False
        """
        return self.is_suffix(other) and self.length() < other.length()

    def has_suffix(self, other):
        """
        Test whether ``self`` has ``other`` as a suffix.

        .. note::

           Some word datatype classes, like :class:`WordDatatype_str`,
           override this method.

        INPUT:

            - ``other`` - a word, or data describing a word

        OUTPUT:

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
        from sage.combinat.words.word import Word
        w = Word(other)
        return w.is_suffix(self)

    def is_prefix(self, other):
        r"""
        Returns True if self is a prefix of other, and False otherwise.

        EXAMPLES::

            sage: w = Word('0123456789')
            sage: y = Word('012345')
            sage: y.is_prefix(w)
            True
            sage: w.is_prefix(y)
            False
            sage: w.is_prefix(Word())
            False
            sage: Word().is_prefix(w)
            True
            sage: Word().is_prefix(Word())
            True
        """
        return self == other[:self.length()]

    def is_proper_prefix(self, other):
        r"""
        Returns True if self is a proper prefix of other, and False otherwise.

        EXAMPLES::

            sage: Word('12').is_proper_prefix(Word('123'))
            True
            sage: Word('12').is_proper_prefix(Word('12'))
            False
            sage: Word().is_proper_prefix(Word('123'))
            True
            sage: Word('123').is_proper_prefix(Word('12'))
            False
            sage: Word().is_proper_prefix(Word())
            False
        """
        return self.is_prefix(other) and self.length() < other.length()

    def has_prefix(self, other):
        r"""
        Test whether ``self`` has ``other`` as a prefix.

        INPUT:

            - ``other`` - a word, or data describing a word

        OUTPUT:

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
        from sage.combinat.words.word import Word
        w = Word(other)
        return w.is_prefix(self)

    def reversal(self):
        r"""
        Returns the reversal of self.

        EXAMPLES::

            sage: Word('124563').reversal()
            word: 365421
        """
        return self[::-1]

    @cached_method
    def prefix_function_table(self):
        r"""
        Returns a vector containing the length of the proper prefix-suffixes
        for all the non-empty prefixes of self.

        EXAMPLES::

            sage: Word('121321').prefix_function_table()
            [0, 0, 1, 0, 0, 1]
            sage: Word('1241245').prefix_function_table()
            [0, 0, 0, 1, 2, 3, 0]
            sage: Word().prefix_function_table()
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

    @cached_method
    def good_suffix_table(self):
        r"""
        Returns a table of the maximum skip you can do in order not to miss
        a possible occurrence of self in a word.

        This is a part of the Boyer-Moore algorithm to find factors. See [1].

        EXAMPLES::

            sage: Word('121321').good_suffix_table()
            [5, 5, 5, 5, 3, 3, 1]
            sage: Word('12412').good_suffix_table()
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

    @cached_method
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
            sage: w.suffix_trie()
            Suffix Trie of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.suffix_trie()
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
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: 01011
        """
        from sage.combinat.words.suffix_trees import ImplicitSuffixTree
        return ImplicitSuffixTree(self)

    @cached_method
    def suffix_tree(self):
        r"""
        Alias for implicit_suffix_tree().

        EXAMPLES::

            sage: Word('abbabaab').suffix_tree()
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

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_factors(i) for i in range(10)]
            [1, 2, 4, 6, 10, 12, 16, 20, 22, 24]

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

    def factor_set(self, n=None):
        r"""
        Returns the set of factors (of length n) of self.

        INPUT:

        - ``n`` - an integer or ``None`` (default: None).

        OUTPUT:

            If ``n`` is an integer, returns the set of all distinct
            factors of length ``n``. If ``n`` is ``None``, returns the set
            of all distinct factors.

        EXAMPLES::

            sage: w = Word('121')
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: 1, word: 12, word: 121, word: 2, word: 21]

        ::

            sage: w = Word('1213121')
            sage: for i in range(w.length()): sorted(w.factor_set(i))
            [word: ]
            [word: 1, word: 2, word: 3]
            [word: 12, word: 13, word: 21, word: 31]
            [word: 121, word: 131, word: 213, word: 312]
            [word: 1213, word: 1312, word: 2131, word: 3121]
            [word: 12131, word: 13121, word: 21312]
            [word: 121312, word: 213121]

        ::

            sage: w = Word([1,2,1,2,3])
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]

        TESTS::

            sage: w = Word("xx")
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: x, word: xx]

        ::

            sage: Set(Word().factor_set())
            {word: }
        """
        return Set(set(self.factor_iterator(n)))

    def topological_entropy(self, n):
        r"""
        Return the topological entropy for the factors of length n.

        The topological entropy of a sequence `u` is defined as the
        exponential growth rate of the complexity of `u` as the length
        increases: `H_{top}(u)=\lim_{n\to\infty}\frac{\log_d(p_u(n))}{n}`
        where `d` denotes the cardinality of the alphabet and `p_u(n)` is
        the complexity function, i.e. the number of factors of length `n`
        in the sequence `u` [1].

        INPUT:

        - ``self`` - a word defined over a finite alphabet
        -  ``n`` - positive integer

        OUTPUT:

        real number (a symbolic expression)

        EXAMPLES::

            sage: W = Words([0, 1])
            sage: w = W([0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1])
            sage: t = w.topological_entropy(3); t
            1/3*log(7)/log(2)
            sage: n(t)
            0.935784974019201

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: topo = w.topological_entropy
            sage: for i in range(0, 41, 5): print i, n(topo(i), digits=5)
            0 1.0000
            5 0.71699
            10 0.48074
            15 0.36396
            20 0.28774
            25 0.23628
            30 0.20075
            35 0.17270
            40 0.14827

        If no alphabet is specified, an error is raised::

            sage: w = Word(range(20))
            sage: w.topological_entropy(3)
            Traceback (most recent call last):
            ...
            TypeError: The word must be defined over a finite alphabet

        The following is ok::

            sage: W = Words(range(20))
            sage: w = W(range(20))
            sage: w.topological_entropy(3)
            1/3*log(18)/log(20)

        REFERENCES:

           [1] N. Pytheas Fogg, Substitutions in Dynamics, Arithmetics,
           and Combinatorics, Lecture Notes in Mathematics 1794, Springer
           Verlag. V. Berthe, S. Ferenczi, C. Mauduit and A. Siegel, Eds.
           (2002).
        """
        d = self.parent().size_of_alphabet()
        if d is Infinity:
            raise TypeError("The word must be defined over a finite alphabet")
        if n == 0:
            return 1
        pn = self.number_of_factors(n)
        from sage.functions.all import log
        return log(pn, base=d)/n

    @cached_method
    def rauzy_graph(self, n):
        r"""
        Returns the Rauzy graph of the factors of length n of self.

        The vertices are the factors of length `n` and there is an edge from
        `u` to `v` if `ua = bv` is a factor of length `n+1` for some letters
        `a` and `b`.

        INPUT:

        - ``n`` - integer

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.rauzy_graph(3); g
            Looped digraph on 8 vertices
            sage: WordOptions(identifier='')
            sage: g.vertices()
            [012, 123, 234, 345, 456, 567, 678, 789]
            sage: g.edges()
            [(012, 123, 3),
             (123, 234, 4),
             (234, 345, 5),
             (345, 456, 6),
             (456, 567, 7),
             (567, 678, 8),
             (678, 789, 9)]
            sage: WordOptions(identifier='word: ')

        ::

            sage: f = words.FibonacciWord()[:100]
            sage: f.rauzy_graph(8)
            Looped digraph on 9 vertices

        ::

            sage: w = Word('1111111')
            sage: g = w.rauzy_graph(3)
            sage: g.edges()
            [(word: 111, word: 111, word: 1)]

        ::

            sage: w = Word('111')
            sage: for i in range(5) : w.rauzy_graph(i)
            Looped multi-digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 0 vertices

        Multi-edges are allowed for the empty word::

            sage: W = Words('abcde')
            sage: w = W('abc')
            sage: w.rauzy_graph(0)
            Looped multi-digraph on 1 vertex
            sage: _.edges()
            [(word: , word: , word: a),
             (word: , word: , word: b),
             (word: , word: , word: c)]
        """
        from sage.graphs.digraph import DiGraph
        multiedges = n == 0
        g = DiGraph(loops=True, multiedges=multiedges)
        if n == self.length():
            g.add_vertex(self)
        else:
            for w in self.factor_iterator(n+1):
                u = w[:-1]
                v = w[1:]
                a = w[-1:]
                g.add_edge(u,v,a)
        return g

    def reduced_rauzy_graph(self, n):
        r"""
        Returns the reduced Rauzy graph of order `n` of self.

        INPUT:

        - ``n`` - non negative integer. Every vertex of a reduced
          Rauzy graph of order `n` is a factor of length `n` of self.

        OUTPUT:

        Looped multi-digraph

        DEFINITION:

        For infinite periodic words (resp. for finite words of type `u^i
        u[0:j]`), the reduced Rauzy graph of order `n` (resp. for `n`
        smaller or equal to `(i-1)|u|+j`) is the directed graph whose
        unique vertex is the prefix `p` of length `n` of self and which has
        an only edge which is a loop on `p` labelled by `w[n+1:|w|] p`
        where `w` is the unique return word to `p`.

        In other cases, it is the directed graph defined as followed.  Let
        `G_n` be the Rauzy graph of order `n` of self. The vertices are the
        vertices of `G_n` that are either special or not prolongable to the
        right or to the left. For each couple (`u`, `v`) of such vertices
        and each directed path in `G_n` from `u` to `v` that contains no
        other vertices that are special, there is an edge from `u` to `v`
        in the reduced Rauzy graph of order `n` whose label is the label of
        the path in `G_n`.

        .. NOTE::

            In the case of infinite recurrent non periodic words, this
            definition correspond to the following one that can be found in
            [1] and [2]  where a simple path is a path that begins with a
            special factor, ends with a special factor and contains no
            other vertices that are special:

            The reduced Rauzy graph of factors of length `n` is obtained
            from `G_n` by replacing each simple path `P=v_1 v_2 ...
            v_{\ell}` with an edge `v_1 v_{\ell}` whose label is the
            concatenation of the labels of the edges of `P`.

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.reduced_rauzy_graph(3); g
            Looped multi-digraph on 2 vertices
            sage: g.vertices()
            [word: 012, word: 789]
            sage: g.edges()
            [(word: 012, word: 789, word: 3456789)]

        For the Fibonacci word::

            sage: f = words.FibonacciWord()[:100]
            sage: g = f.reduced_rauzy_graph(8);g
            Looped multi-digraph on 2 vertices
            sage: g.vertices()
            [word: 01001010, word: 01010010]
            sage: g.edges()
            [(word: 01001010, word: 01010010, word: 010), (word: 01010010, word: 01001010, word: 01010), (word: 01010010, word: 01001010, word: 10)]

        For periodic words::

            sage: from itertools import cycle
            sage: w = Word(cycle('abcd'))[:100]
            sage: g = w.reduced_rauzy_graph(3)
            sage: g.edges()
            [(word: abc, word: abc, word: dabc)]

        ::

            sage: w = Word('111')
            sage: for i in range(5) : w.reduced_rauzy_graph(i)
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped multi-digraph on 1 vertex
            Looped multi-digraph on 0 vertices

        For ultimately periodic words::

            sage: sigma = WordMorphism('a->abcd,b->cd,c->cd,d->cd')
            sage: w = sigma.fixed_point('a')[:100]; w
            word: abcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd...
            sage: g = w.reduced_rauzy_graph(5)
            sage: g.vertices()
            [word: abcdc, word: cdcdc]
            sage: g.edges()
            [(word: abcdc, word: cdcdc, word: dc), (word: cdcdc, word: cdcdc, word: dc)]

        AUTHOR:

        Julien Leroy (March 2010): initial version

        REFERENCES:

        - [1] M. Bucci et al.  A. De Luca, A. Glen, L. Q. Zamboni, A
          connection between palindromic and factor complexity using
          return words," Advances in Applied Mathematics 42 (2009) 60-74.

        - [2] L'ubomira Balkova, Edita Pelantova, and Wolfgang Steiner.
          Sequences with constant number of return words. Monatsh. Math,
          155 (2008) 251-263.
        """
        from sage.graphs.digraph import DiGraph
        from copy import copy
        g = copy(self.rauzy_graph(n))
        # Otherwise it changes the rauzy_graph function.
        l = [v for v in g if g.in_degree(v)==1 and g.out_degree(v)==1]
        if g.num_verts() != 0 and len(l) == g.num_verts():
            # In this case, the Rauzy graph is simply a cycle.
            g = DiGraph()
            g.allow_loops(True)
            g.add_vertex(self[:n])
            g.add_edge(self[:n],self[:n],self[n:n+len(l)])
        else:
            g.allow_loops(True)
            g.allow_multiple_edges(True)
            for v in l:
                [i] = g.neighbors_in(v)
                [o] = g.neighbors_out(v)
                g.add_edge(i,o,g.edge_label(i,v)[0]*g.edge_label(v,o)[0])
                g.delete_vertex(v)
        return g

    def left_special_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the left special factors (of length n).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all left special factors.

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.left_special_factors_iterator(3))
            [word: 000, word: 010]
            sage: sorted(w.left_special_factors_iterator(4))
            [word: 0000, word: 0101]
            sage: sorted(w.left_special_factors_iterator(5))
            [word: 00000, word: 01010]
        """
        if n is None:
            for i in range(self.length()):
                for w in self.left_special_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            in_d = g.in_degree
            for v in g:
                if in_d(v) > 1:
                    yield v

    def left_special_factors(self, n=None):
        r"""
        Returns the left special factors (of length n).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: ``None``). If ``None``, it
           returns all left special factors.

        OUTPUT:

        A list of words.

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: for i in range(5): print i, sorted(w.left_special_factors(i))
            0 [word: ]
            1 [word: 0]
            2 [word: 00, word: 01]
            3 [word: 000, word: 010]
            4 [word: 0000, word: 0101]
        """
        return list(self.left_special_factors_iterator(n))

    def right_special_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the right special factors (of length n).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all right special factors.

        EXAMPLES::

            sage: alpha, beta, x = 0.61, 0.54, 0.3
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.right_special_factors_iterator(3))
            [word: 010, word: 101]
            sage: sorted(w.right_special_factors_iterator(4))
            [word: 0101, word: 1010]
            sage: sorted(w.right_special_factors_iterator(5))
            [word: 00101, word: 11010]
        """
        if n is None:
            for i in range(self.length()):
                for w in self.right_special_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            out_d = g.out_degree
            for v in g:
                if out_d(v) > 1:
                    yield v

    def right_special_factors(self, n=None):
        r"""
        Returns the right special factors (of length n).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           all right special factors.

        OUTPUT:

        A list of words.

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(5): print i, sorted(w.right_special_factors(i))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 001, word: 010, word: 101, word: 110]
            4 [word: 0110, word: 1001]
        """
        return list(self.right_special_factors_iterator(n))

    def bispecial_factors_iterator(self, n=None):
        r"""
        Returns an iterator over the bispecial factors (of length n).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           an iterator over all bispecial factors.

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10):
            ...     for u in sorted(w.bispecial_factors_iterator(i)):
            ...         print i,u
            0
            1 0
            1 1
            2 01
            2 10
            3 010
            3 101
            4 0110
            4 1001
            6 011001
            6 100110
            8 10010110

        ::

            sage: key = lambda u : (len(u), u)
            sage: for u in sorted(w.bispecial_factors_iterator(), key=key): u
            word:
            word: 0
            word: 1
            word: 01
            word: 10
            word: 010
            word: 101
            word: 0110
            word: 1001
            word: 011001
            word: 100110
            word: 10010110
        """
        if n is None:
            for i in range(self.length()):
                for w in self.bispecial_factors_iterator(i):
                    yield w
        else:
            g = self.rauzy_graph(n)
            in_d = g.in_degree
            out_d = g.out_degree
            for v in g:
                if out_d(v) > 1 and in_d(v) > 1:
                    yield v

    def bispecial_factors(self, n=None):
        r"""
        Returns the bispecial factors (of length n).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        -  ``n`` - integer (optional, default: None). If None, it returns
           all bispecial factors.

        OUTPUT:

        A list of words.

        EXAMPLES::

            sage: w = words.FibonacciWord()[:30]
            sage: w.bispecial_factors()
            [word: , word: 0, word: 010, word: 010010, word: 01001010010]

        ::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10): print i, sorted(w.bispecial_factors(i))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 010, word: 101]
            4 [word: 0110, word: 1001]
            5 []
            6 [word: 011001, word: 100110]
            7 []
            8 [word: 10010110]
            9 []
        """
        return list(self.bispecial_factors_iterator(n))

    def number_of_left_special_factors(self, n):
        r"""
        Returns the number of left special factors of length n.

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer

        OUTPUT:

        Non negative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        L = self.rauzy_graph(n).in_degree()
        return sum(1 for i in L if i>1)

    def number_of_right_special_factors(self, n):
        r"""
        Returns the number of right special factors of length n.

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        -  ``n`` - integer

        OUTPUT:

        Non negative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        L = self.rauzy_graph(n).out_degree()
        return sum(1 for i in L if i>1)

    def commutes_with(self, other):
        r"""
        Returns True if self commutes with other, and False otherwise.

        EXAMPLES::

            sage: Word('12').commutes_with(Word('12'))
            True
            sage: Word('12').commutes_with(Word('11'))
            False
            sage: Word().commutes_with(Word('21'))
            True
        """
        return (self * other) == (other * self)

    def conjugate(self, pos):
        r"""
        Returns the conjugate at pos of self.

        pos can be any integer, the distance used is the modulo by the length
        of self.

        EXAMPLES::

            sage: Word('12112').conjugate(1)
            word: 21121
            sage: Word().conjugate(2)
            word:
            sage: Word('12112').conjugate(8)
            word: 12121
            sage: Word('12112').conjugate(-1)
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

            sage: Word('cbbca')._conjugates_list()
            [word: cbbca, word: bbcac, word: bcacb, word: cacbb, word: acbbc]
            sage: Word('abcabc')._conjugates_list()
            [word: abcabc,
             word: bcabca,
             word: cabcab,
             word: abcabc,
             word: bcabca,
             word: cabcab]
            sage: Word()._conjugates_list()
            [word: ]
            sage: Word('a')._conjugates_list()
            [word: a]
        """
        S = [self]
        for i in range(1,self.length()):
            S.append(self.conjugate(i))
        return S

    def conjugates_iterator(self):
        r"""
        Returns an iterator over the conjugates of self.

        EXAMPLES::

            sage: it = Word(range(4)).conjugates_iterator()
            sage: for w in it: w
            word: 0123
            word: 1230
            word: 2301
            word: 3012
        """
        yield self
        for i in range(1, self.primitive_length()):
            yield self.conjugate(i)

    def conjugates(self):
        r"""
        Returns the list of unique conjugates of self.

        EXAMPLES::

            sage: Word(range(6)).conjugates()
            [word: 012345,
             word: 123450,
             word: 234501,
             word: 345012,
             word: 450123,
             word: 501234]
            sage: Word('cbbca').conjugates()
            [word: cbbca, word: bbcac, word: bcacb, word: cacbb, word: acbbc]

        The result contains each conjugate only once::

            sage: Word('abcabc').conjugates()
            [word: abcabc, word: bcabca, word: cabcab]

        TESTS::

            sage: Word().conjugates()
            [word: ]
            sage: Word('a').conjugates()
            [word: a]
        """
        return list(self.conjugates_iterator())

    def conjugate_position(self, other):
        r"""
        Returns the position where self is conjugate with other.
        Returns None if there is no such position.

        EXAMPLES::

            sage: Word('12113').conjugate_position(Word('31211'))
            1
            sage: Word('12131').conjugate_position(Word('12113')) is None
            True
            sage: Word().conjugate_position(Word('123')) is None
            True

        TESTS:

        We check that trac #11128 is fixed::

            sage: w = Word([0,0,1,0,2,1])
            sage: [w.conjugate(i).conjugate_position(w) for i in range(w.length())]
            [0, 1, 2, 3, 4, 5]
        """
        if self.length() != other.length():
            return None
        other_square = other * other
        pos = other_square.find(self)
        return pos if pos != -1 else None

    def is_conjugate_with(self, other):
        r"""
        Returns True if self is a conjugate of other, and False otherwise.

        INPUT:

        - ``other`` - a finite word

        OUPUT

        bool

        EXAMPLES::

            sage: w = Word([0..20])
            sage: z = Word([7..20] + [0..6])
            sage: w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
            sage: z
            word: 7,8,9,10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5,6
            sage: w.is_conjugate_with(z)
            True
            sage: z.is_conjugate_with(w)
            True
            sage: u = Word([4]*21)
            sage: u.is_conjugate_with(w)
            False
            sage: u.is_conjugate_with(z)
            False

        Both words must be finite::

            sage: w = Word(iter([2]*100),length='unknown')
            sage: z = Word([2]*100)
            sage: z.is_conjugate_with(w) #TODO: Not implemented for word of unknown length
            True
            sage: wf = Word(iter([2]*100),length='finite')
            sage: z.is_conjugate_with(wf)
            True
            sage: wf.is_conjugate_with(z)
            True

        TESTS::

            sage: Word('11213').is_conjugate_with(Word('31121'))
            True
            sage: Word().is_conjugate_with(Word('123'))
            False
            sage: Word('112131').is_conjugate_with(Word('11213'))
            False
            sage: Word('12131').is_conjugate_with(Word('11213'))
            True

        We make sure that trac #11128 is fixed::

            sage: Word('abaa').is_conjugate_with(Word('aaba'))
            True
            sage: Word('aaba').is_conjugate_with(Word('abaa'))
            True
        """
        return self.length() == other.length() and self.is_factor(other * other)

    def is_cadence(self, seq):
        r"""
        Returns True if seq is a cadence of self, and False otherwise.

        A *cadence* is an increasing sequence of indexes that all map to
        the same letter.

        EXAMPLES::

            sage: Word('121132123').is_cadence([0, 2, 6])
            True
            sage: Word('121132123').is_cadence([0, 1, 2])
            False
            sage: Word('121132123').is_cadence([])
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

    def longest_common_suffix(self, other):
        r"""
        Returns the longest common suffix of self and other.

        EXAMPLES::

            sage: w = Word('112345678')
            sage: u = Word('1115678')
            sage: w.longest_common_suffix(u)
            word: 5678
            sage: u.longest_common_suffix(u)
            word: 1115678
            sage: u.longest_common_suffix(w)
            word: 5678
            sage: w.longest_common_suffix(w)
            word: 112345678
            sage: y = Word('549332345')
            sage: w.longest_common_suffix(y)
            word:

        TESTS:

        With the empty word::

            sage: w.longest_common_suffix(Word())
            word:
            sage: Word().longest_common_suffix(w)
            word:
            sage: Word().longest_common_suffix(Word())
            word:

        With an infinite word::

            sage: t=words.ThueMorseWord('ab')
            sage: w.longest_common_suffix(t)
            Traceback (most recent call last):
            ...
            TypeError: other must be a finite word
        """
        if not isinstance(other, FiniteWord_class):
            raise TypeError("other must be a finite word")

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

        -  ``f`` - involution (default: ``None``) on the alphabet of self. It
           must be callable on letters as well as words (e.g.
           :class:`~sage.combinat.words.morphism.WordMorphism`). The
           default value corresponds to usual palindromes, i.e., `f` equal to
           the identity.

        EXAMPLES::

            sage: Word('esope reste ici et se repose').is_palindrome()
            False
            sage: Word('esoperesteicietserepose').is_palindrome()
            True
            sage: Word('I saw I was I').is_palindrome()
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
            KeyError: 'b'

        TESTS:

        If the given involution is not an involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('abab').is_palindrome(f)
            Traceback (most recent call last):
            ...
            TypeError: self (=a->b, b->b) is not an endomorphism

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
            return self[:l//2] == self[l//2 + l%2:].reversal()
        else:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            return self[:l//2 + l%2] == f(self[l//2:].reversal())

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
            sage: w = Word('33412321')
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
            ...     pal
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
            ...     pal
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
        if l is None:
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

    @cached_method
    def palindromic_lacunas_study(self, f=None):
        r"""
        Returns interesting statistics about longest (`f`-)palindromic suffixes
        and lacunas of self (see [1] and [2]).

        Note that a word `w` has at most `|w| + 1` different palindromic factors
        (see [3]). For `f`-palindromes (or pseudopalidromes or theta-palindromes),
        the maximum number of `f`-palindromic factors is `|w|+1-g_f(w)`, where
        `g_f(w)` is the number of pairs `\{a, f(a)\}` such that `a` is a letter,
        `a` is not equal to `f(a)`, and `a` or `f(a)` occurs in `w`, see [4].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism). The
           default value corresponds to usual palindromes, i.e., `f` equal to
           the identity.

        OUTPUT:

        -  ``list`` - list of the length of the longest palindromic
           suffix (lps) for each non-empty prefix of self;
        -  ``list`` - list of all the lacunas, i.e. positions where there is no
           unioccurrent lps;
        -  ``set`` - set of palindromic factors of self.

        EXAMPLES::

            sage: a,b,c = Word('abbabaabbaab').palindromic_lacunas_study()
            sage: a
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]
            sage: b
            [8, 9]
            sage: c          # random order
            set([word: , word: b, word: bab, word: abba, word: bb, word: aa, word: baabbaab, word: baab, word: aba, word: aabbaa, word: a])

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: a,b,c = Word('abbabaab').palindromic_lacunas_study(f=f)
            sage: a
            [0, 2, 0, 2, 2, 4, 6, 8]
            sage: b
            [0, 2, 4]
            sage: c           # random order
            set([word: , word: ba, word: baba, word: ab, word: bbabaa, word: abbabaab])
            sage: c == set([Word(), Word('ba'), Word('baba'), Word('ab'), Word('bbabaa'), Word('abbabaab')])
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
        -   [4] Š. Starosta, On Theta-palindromic Richness, Theoret. Comp.
            Sci. 412 (2011) 1111--1121
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

        A *lacuna* is a position in a word where the longest (`f`-)palindromic
        suffix is not unioccurrent (see [1]).

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism). The
           default value corresponds to usual palindromes, i.e., `f` equal to
           the identity.

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
        Returns the list of the lengths of the unioccurrent longest
        (`f`)-palindromic suffixes (lps) for each non-empty prefix of self. No
        unioccurrent lps are indicated by None.

        It corresponds to the function `H_w` defined in [1] and [2].

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism). The
           default value corresponds to usual palindromes, i.e., `f` equal to
           the identity.

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

            sage: sorted(Word('01101001').palindromes())
            [word: , word: 0, word: 00, word: 010, word: 0110, word: 1, word: 1001, word: 101, word: 11]
            sage: sorted(Word('00000').palindromes())
            [word: , word: 0, word: 00, word: 000, word: 0000, word: 00000]
            sage: sorted(Word('0').palindromes())
            [word: , word: 0]
            sage: sorted(Word('').palindromes())
            [word: ]
            sage: sorted(Word().palindromes())
            [word: ]
            sage: f = WordMorphism('a->b,b->a')
            sage: sorted(Word('abbabaab').palindromes(f))
            [word: , word: ab, word: abbabaab, word: ba, word: baba, word: bbabaa]
        """
        return self.palindromic_lacunas_study(f=f)[2]

    def palindrome_prefixes(self):
        r"""
        Returns a list of all palindrome prefixes of self.

        OUTPUT:

            list -- A list of all palindrome prefixes of self.

        EXAMPLES::

            sage: w = Word('abaaba')
            sage: w.palindrome_prefixes()
            [word: , word: a, word: aba, word: abaaba]
            sage: w = Word('abbbbbbbbbb')
            sage: w.palindrome_prefixes()
            [word: , word: a]
        """
        return list(self.palindrome_prefixes_iterator())

    def defect(self, f=None):
        r"""
        Returns the defect of self.

        The *defect* of a finite word `w` is given by the difference between
        the maximum number of possible palindromic factors in a word of length
        `|w|` and the actual number of palindromic factors contained in `w`.
        It is well known that the maximum number of palindromic factors in `w`
        is `|w|+1` (see [DJP01]_).

        An optional involution on letters ``f`` can be given. In that case, the
        *f-palindromic defect* (or *pseudopalindromic defect*, or
        *theta-palindromic defect*) of `w` is returned. It is a
        generalization of defect to f-palindromes. More precisely, the defect is
        `D(w)=|w|+1-g_f(w)-|PAL_f(w)|`, where `PAL_f(w)` denotes the set of
        f-palindromic factors of `w` (including the empty word) and `g_f(w)` is
        the number of pairs `\{a, f(a)\}` such that `a` is a letter, `a` is not
        equal to `f(a)`, and `a` or `f(a)` occurs in `w`. In the case of usual
        palindromes (i.e., for ``f`` not given or equal to the identity),
        `g_f(w) = 0` for all `w`. See [BHNR04]_ for usual palindromes and [Sta11]_
        for f-palindromes.

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of self. It must
           be callable on letters as well as words (e.g. WordMorphism). The
           default value corresponds to usual palindromes, i.e., `f` equal to
           the identity.

        OUTPUT:

            integer -- If f is None, the palindromic defect of self;
                       otherwise, the f-palindromic defect of self.

        EXAMPLES::

            sage: Word('ara').defect()
            0
            sage: Word('abcacba').defect()
            1

        It is known that Sturmian words (see [DJP01]_) have zero defect::

            sage: words.FibonacciWord()[:100].defect()
            0

            sage: sa = WordMorphism('a->ab,b->b')
            sage: sb = WordMorphism('a->a,b->ba')
            sage: w = (sa*sb*sb*sa*sa*sa*sb).fixed_point('a')
            sage: w[:30].defect()
            0
            sage: w[110:140].defect()
            0

        It is even conjectured that the defect of an aperiodic word which is
        a fixed point of a primitive morphism is either `0` or infinite
        (see [BBGL08]_)::

            sage: w = words.ThueMorseWord()
            sage: w[:50].defect()
            12
            sage: w[:100].defect()
            16
            sage: w[:300].defect()
            52

        For generalized defect with an involution different from the identity,
        there is always a letter which is not a palindrome! This is the reason
        for the modification of the definition::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('a').defect(f)
            0
            sage: Word('ab').defect(f)
            0
            sage: Word('aa').defect(f)
            1
            sage: Word('abbabaabbaababba').defect(f)
            3

        ::

            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: Word('cabc').defect(f)
            0
            sage: Word('abcaab').defect(f)
            2

        Other examples::

            sage: Word('000000000000').defect()
            0
            sage: Word('011010011001').defect()
            2
            sage: Word('0101001010001').defect()
            0
            sage: Word().defect()
            0
            sage: Word('abbabaabbaababba').defect()
            2

        REFERENCES:

        .. [BBGL08] A. Blondin Massé, S. Brlek, A. Garon, and S. Labbé,
           Combinatorial properties of f -palindromes in the Thue-Morse
           sequence. Pure Math. Appl., 19(2-3):39--52, 2008.

        .. [BHNR04] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the
           Palindromic Complexity of Infinite Words, in J. Berstel, J.
           Karhumaki, D. Perrin, Eds, Combinatorics on Words with Applications,
           International Journal of Foundation of Computer Science, Vol. 15,
           No. 2 (2004) 293--306.

        .. [DJP01] X. Droubay, J. Justin, G. Pirillo, Episturmian words and some
           constructions of de Luca and Rauzy, Theoret. Comput. Sci. 255,
           (2001), no. 1--2, 539--553.

        .. [Sta11] Š. Starosta, On Theta-palindromic Richness, Theoret. Comp.
           Sci. 412 (2011) 1111--1121
        """
        g_w = 0
        if f is not None:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            D = f.domain()
            A = set(map(D,set(self)))
            while A:
                x = A.pop()
                if f(x) != x: # count only non f-palindromic letters
                    if f(x) in A:
                        A.remove(f(x))
                    g_w +=1

        return self.length()+1-g_w-len(self.palindromes(f=f))

    def is_full(self, f=None):
        r"""
        Returns True if self has defect 0, and False otherwise.

        A word is *full* (or *rich*) if its defect is zero (see [1]).
        If ``f`` is given, then the f-palindromic defect is used (see [2]).

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
            True

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abab').is_full(f)
            True
            sage: Word('abba').is_full(f)
            False

        A simple example of an infinite word full of f-palindromes::

            sage: p = WordMorphism({0:'abc',1:'ab'})
            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: p(words.FibonacciWord()[:50]).is_full(f)
            True
            sage: p(words.FibonacciWord()[:150]).is_full(f)
            True

        REFERENCES:

        -   [1] S. Brlek, S. Hamel, M. Nivat, C. Reutenauer, On the Palindromic
            Complexity of Infinite Words, in J. Berstel, J. Karhumaki,
            D. Perrin, Eds, Combinatorics on Words with Applications,
            International Journal of Foundation of Computer Science, Vol. 15,
            No. 2 (2004) 293--306.

        -   [2] E. Pelantová, Š. Starosta, Infinite words rich and almost rich
            in generalized palindromes, in: G. Mauri, A. Leporati (Eds.),
            Developments in Language Theory, volume 6795 of Lecture Notes
            in Computer Science, Springer-Verlag, Berlin, Heidelberg, 2011,
            pp. 406--416
        """
        return self.defect(f=f) == 0

    is_rich = is_full

    def palindromic_closure(self, side='right', f=None):
        r"""
        Return the shortest palindrome having ``self`` as a prefix
        (or as a suffix if ``side`` is ``'left'``).

        See [1].

        INPUT:

        -  ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``) the
           direction of the  closure

        -  ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
           It must be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

        word -- If f is ``None``, the right palindromic closure of ``self``;
                otherwise, the right ``f``-palindromic closure of ``self``.
                If side is ``'left'``, the left palindromic closure.

        EXAMPLES::

            sage: Word('1233').palindromic_closure()
            word: 123321
            sage: Word('12332').palindromic_closure()
            word: 123321
            sage: Word('0110343').palindromic_closure()
            word: 01103430110
            sage: Word('0110343').palindromic_closure(side='left')
            word: 3430110343
            sage: Word('01105678').palindromic_closure(side='left')
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
            KeyError: 'b'

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
                raise ValueError("side must be either 'left' or 'right' (not %s) " % side)
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            if side == 'right':
                l = self.lps(f=f).length()
                return self * f(self[-(l+1)::-1])
            elif side == 'left':
                l = self.reversal().lps(f=f).length()
                return f(self[:l-1:-1]) * self
            else:
                raise ValueError("side must be either 'left' or 'right' (not %s) " % side)

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

            sage: Word('abbabab').is_symmetric()
            True
            sage: Word('ababa').is_symmetric()
            True
            sage: Word('aababaabba').is_symmetric()
            False
            sage: Word('aabbbaababba').is_symmetric()
            False
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('aabbbaababba').is_symmetric(f)
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

            sage: Word('121').length_border()
            1
            sage: Word('1').length_border()
            0
            sage: Word('1212').length_border()
            2
            sage: Word('111').length_border()
            2
            sage: Word().length_border() is None
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

            sage: Word('121212').border()
            word: 1212
            sage: Word('12321').border()
            word: 1
            sage: Word().border() is None
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

            sage: Word('12121').is_overlap()
            True
            sage: Word('123').is_overlap()
            False
            sage: Word('1231').is_overlap()
            False
            sage: Word('123123').is_overlap()
            False
            sage: Word('1231231').is_overlap()
            True
            sage: Word().is_overlap()
            False
        """
        if self.length() == 0:
            return False
        return self.length_border() > self.length()//2

    def primitive_length(self):
        r"""
        Returns the length of the primitive of self.

        EXAMPLES::

            sage: Word('1231').primitive_length()
            4
            sage: Word('121212').primitive_length()
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

            sage: Word('1231').is_primitive()
            True
            sage: Word('111').is_primitive()
            False
        """
        return self.length() == self.primitive_length()

    def primitive(self):
        r"""
        Returns the primitive of self.

        EXAMPLES::

            sage: Word('12312').primitive()
            word: 12312
            sage: Word('121212').primitive()
            word: 12
        """
        return self[:self.primitive_length()]

    def exponent(self):
        r"""
        Returns the exponent of self.

        OUTPUT:

            integer -- the exponent

        EXAMPLES::

            sage: Word('1231').exponent()
            1
            sage: Word('121212').exponent()
            3
            sage: Word().exponent()
            0
        """
        if self.length() == 0:
            return 0
        return self.length() // self.primitive_length()

    def has_period(self, p):
        r"""
        Returns True if self has the period ``p``,
        False otherwise.

        .. NOTE::

            By convention, integers greater than the length
            of self are periods of self.

        INPUT:

        - ``p`` - an integer to check if it is a period
          of self.

        EXAMPLES::

            sage: w = Word('ababa')
            sage: w.has_period(2)
            True
            sage: w.has_period(3)
            False
            sage: w.has_period(4)
            True
            sage: w.has_period(-1)
            False
            sage: w.has_period(5)
            True
            sage: w.has_period(6)
            True
        """
        if p < 0:
            return False
        elif p >= len(self):
            return True
        else:
            for i in range(len(self) - p):
                if self[i] != self[i + p]:
                    return False
            return True

    def periods(self, divide_length=False):
        r"""
        Returns a list containing the periods of self
        between `1` and `n - 1`, where `n` is the length
        of self.

        INPUT:

        - ``divide_length`` - boolean (default: False).
          When set to True, then only periods that divide
          the length of self are considered.

        OUTPUT:

        List of positive integers

        EXAMPLES::

            sage: w = Word('ababab')
            sage: w.periods()
            [2, 4]
            sage: w.periods(divide_length=True)
            [2]
            sage: w = Word('ababa')
            sage: w.periods()
            [2, 4]
            sage: w.periods(divide_length=True)
            []
        """
        n = len(self)
        if divide_length:
            possible = (i for i in xrange(1,n) if n % i == 0)
        else:
            possible = xrange(1, n)
        return [x for x in possible if self.has_period(x)]

    def longest_common_subword(self,other):
        r"""
        Returns a longest subword of ``self`` and ``other``.

        A subword of a word is a subset of the word's letters, read in the
        order in which they appear in the word.

        For more information, see
        :wikipedia:`Longest_common_subsequence_problem`.

        INPUT:

        - ``other`` -- a word

        ALGORITHM:

        For any indices `i,j`, we compute the longest common subword ``lcs[i,j]`` of
        `self[:i]` and `other[:j]`. This can be easily obtained as the longest
        of

        - ``lcs[i-1,j]``

        - ``lcs[i,j-1]``

        - ``lcs[i-1,j-1]+self[i]`` if ``self[i]==other[j]``

        EXAMPLES::

            sage: v1 = Word("abc")
            sage: v2 = Word("ace")
            sage: v1.longest_common_subword(v2)
            word: ac

            sage: w1 = Word("1010101010101010101010101010101010101010")
            sage: w2 = Word("0011001100110011001100110011001100110011")
            sage: w1.longest_common_subword(w2)
            word: 00110011001100110011010101010

        TESTS::

            sage: Word().longest_common_subword(Word())
            word:

        .. SEEALSO::

            :meth:`is_subword_of`
        """
        from sage.combinat.words.word import Word
        if len(self) == 0 or len(other) == 0:
            return Word()

        w2 = list(other)

        # In order to avoid storing lcs[i,j] for each pair i,j of indices, we
        # only store the lcs[i,j] for two consecutive values of i. At any step
        # of the algorithm, lcs[i,j] is stored at lcs[0][j] and lcs[i-1,j] is
        # stored at lcs[1][j]

        # The weird +1 that follows exists to make sure that lcs[i,-1] returns
        # the empty word.
        lcs = [[[] for i in range(len(w2)+1)] for j in range(2)]

        for i,l1 in enumerate(self):
            for j,l2 in enumerate(other):
                lcs[0][j] = max(lcs[0][j-1], lcs[1][j],
                                lcs[1][j-1] + ([l1] if l1==l2 else []),key=len)

            # Maintaining the meaning of lcs for the next loop
            lcs.pop(1)
            lcs.insert(0,[[] for i in range(len(w2)+1)])

        return Word(lcs[1][-2])

    def is_subword_of(self, other):
        r"""
        Returns True is self is a subword of other, and False otherwise.

        EXAMPLES::

            sage: Word().is_subword_of(Word('123'))
            True
            sage: Word('123').is_subword_of(Word('3211333213233321'))
            True
            sage: Word('321').is_subword_of(Word('11122212112122133111222332'))
            False

        .. SEEALSO::

            :meth:`longest_common_subword`

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

            sage: Word('123132133').is_lyndon()
            True
            sage: Word().is_lyndon()
            True
            sage: Word('122112').is_lyndon()
            False

        TESTS:

        A sanity check: ``LyndonWords`` generates Lyndon words, so we
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

        raise TypeError("degree is not defined for your word")

    def deg_lex_less(self, other, weights=None):
        r"""
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
            sage: Word("abba").deg_lex_less(Word("abbb"), dict(a=1,b=2))
            True
            sage: Word("abba").deg_lex_less(Word("baba"), dict(a=1,b=2))
            True
            sage: Word("abba").deg_lex_less(Word("aaba"), dict(a=1,b=2))
            False
            sage: Word("abba").deg_lex_less(Word("aaba"), dict(a=1,b=0))
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

            sage: Word([1,2,4]).inv_lex_less(Word([1,3,2]))
            False
            sage: Word([3,2,1]).inv_lex_less(Word([1,2,3]))
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

            sage: Word([1,2,4]).deg_inv_lex_less(Word([1,3,2]))
            False
            sage: Word([3,2,1]).deg_inv_lex_less(Word([1,2,3]))
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

            sage: Word([1,2,4]).rev_lex_less(Word([1,3,2]))
            True
            sage: Word([3,2,1]).rev_lex_less(Word([1,2,3]))
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
            sage: Word([1,2,4]).deg_rev_lex_less(Word([1,3,2]))
            False
            sage: Word([1,2,3]).deg_rev_lex_less(Word([1,2,4]))
            True
        """
        d1 = self.degree(weights)
        d2 = other.degree(weights)
        if d1 != d2:
            return d1 < d2
        return self.rev_lex_less(other)

    @cached_method
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

        The empty word occurs in a word::

            sage: Word('')._pos_in(Word('123'), 0)
            0
            sage: Word('')._pos_in(Word(''), 0)
            0
        """
        lf = self.length()
        lm = other.length()
        if lf == 0:
            return p
        elif lm == 0:
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

    def find(self, sub, start=0, end=None):
        r"""
        Returns the index of the first occurrence of sub in self,
        such that sub is contained within self[start:end].
        Returns -1 on failure.

        INPUT:

        -  ``sub`` - string or word to search for.
        -  ``start`` - non negative integer (default: 0) specifying
           the position from which to start the search.
        -  ``end`` - non negative integer (default: None) specifying
           the position at which the search must stop. If None, then
           the search is performed up to the end of the string.

        OUTPUT:

            non negative integer or -1

        EXAMPLES::

            sage: w = Word([0,1,0,0,1])
            sage: w.find(Word([0,1]))
            0
            sage: w.find(Word([0,1]), start=1)
            3
            sage: w.find(Word([0,1]), start=1, end=5)
            3
            sage: w.find(Word([0,1]), start=1, end=4) == -1
            True
            sage: w.find(Word([1,1])) == -1
            True

        Instances of Word_str handle string inputs as well::

            sage: w = Word('abac')
            sage: w.find('a')
            0
            sage: w.find(Word('a'))
            0
        """
        w = self[start:end]
        if isinstance(sub, FiniteWord_class):
            p = sub.first_pos_in(w)
            if p is None:
                return -1
            else:
                return p + start
        else:
            L = len(sub)
            if start is None:
                i = len(self) - L
            else:
                i = start - L
            while i >= end:
                if self[i:i+L] == sub: return i
                i -= 1
            return -1

    def rfind(self, sub, start=0, end=None):
        r"""
        Returns the index of the last occurrence of sub in self,
        such that sub is contained within self[start:end].
        Returns -1 on failure.

        INPUT:

        -  ``sub`` - string or word to search for.
        -  ``start`` - non negative integer (default: 0) specifying
           the position at which the search must stop.
        -  ``end`` - non negative integer (default: None) specifying
           the position from which to start the search. If None, then
           the search is performed up to the end of the string.

        OUTPUT:

            non negative integer or -1

        EXAMPLES::

            sage: w = Word([0,1,0,0,1])
            sage: w.rfind(Word([0,1]))
            3
            sage: w.rfind(Word([0,1]), end=4)
            0
            sage: w.rfind(Word([0,1]), end=5)
            3
            sage: w.rfind(Word([0,0]), start=2, end=5)
            2
            sage: w.rfind(Word([0,0]), start=3, end=5) == -1
            True

        Instances of Word_str handle string inputs as well::

            sage: w = Word('abac')
            sage: w.rfind('a')
            2
            sage: w.rfind(Word('a'))
            2
        """
        L = len(sub)
        if end is None:
            i = len(self) - L
        else:
            i = end - L
        while i >= start:
            if self[i:i+L] == sub: return i
            i -= 1
        return -1

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

        The empty word is factor of another word::

            sage: Word().is_factor(Word())
            True
            sage: Word().is_factor(Word('a'))
            True
            sage: Word().is_factor(Word([1,2,3]))
            True
            sage: Word().is_factor(Word(lambda n:n, length=5))
            True
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
            raise NotImplementedError("The factor must be non empty")
        p = self._pos_in(other, 0)
        while p is not None:
            yield p
            p = self._pos_in(other, p+1)

    def nb_factor_occurrences_in(self, other):
        r"""
        Returns the number of times self appears as a factor
        in other.

        EXAMPLES::

            sage: Word().nb_factor_occurrences_in(Word('123'))
            Traceback (most recent call last):
            ...
            NotImplementedError: The factor must be non empty
            sage: Word('123').nb_factor_occurrences_in(Word('112332312313112332121123'))
            4
            sage: Word('321').nb_factor_occurrences_in(Word('11233231231311233221123'))
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

            sage: Word().nb_subword_occurrences_in(Word('123'))
            Traceback (most recent call last):
            ...
            NotImplementedError: undefined value
            sage: Word('123').nb_subword_occurrences_in(Word('1133432311132311112'))
            11
            sage: Word('4321').nb_subword_occurrences_in(Word('1132231112233212342231112'))
            0
            sage: Word('3').nb_subword_occurrences_in(Word('122332112321213'))
            4
        """
        ls = self.length()
        if ls == 0:
            raise NotImplementedError("undefined value")
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

        INPUT:

        - ``fact`` - a non empty finite word

        OUTPUT:

        Python list of finite words

        TESTS::

            sage: Word('baccabccbacbca')._return_words_list(Word('b'))
            [word: bacca, word: bcc, word: bac]
        """
        return list(self.return_words_iterator(fact))

    def return_words(self, fact):
        r"""
        Returns the set of return words of fact in self.

        This is the set of all factors starting by the given factor and ending
        just before the next occurrence of this factor. See [1] and [2].

        INPUT:

        - ``fact`` - a non empty finite word

        OUTPUT:

        Python set of finite words

        EXAMPLES::

            sage: Word('21331233213231').return_words(Word('2'))
            set([word: 213, word: 21331, word: 233])
            sage: Word().return_words(Word('213'))
            set([])
            sage: Word('121212').return_words(Word('1212'))
            set([word: 12])

        ::

            sage: TM = words.ThueMorseWord()[:10000]
            sage: TM.return_words(Word([0]))     # optional long time (1.34 s)
            set([word: 0, word: 01, word: 011])

        REFERENCES:

        -   [1] F. Durand, A characterization of substitutive sequences using
            return words, Discrete Math. 179 (1998) 89-101.
        -   [2] C. Holton, L.Q. Zamboni, Descendants of primitive substitutions,
            Theory Comput. Syst. 32 (1999) 133-157.
        """
        return set(self.return_words_iterator(fact))

    def complete_return_words(self, fact):
        r"""
        Returns the set of complete return words of fact in self.

        This is the set of all factors starting by the given factor and ending
        just after the next occurrence of this factor. See for instance [1].

        INPUT:

        - ``fact`` - a non empty finite word

        OUTPUT:

        Python set of finite words

        EXAMPLES::

            sage: s = Word('21331233213231').complete_return_words(Word('2'))
            sage: sorted(s)
            [word: 2132, word: 213312, word: 2332]
            sage: Word('').complete_return_words(Word('213'))
            set([])
            sage: Word('121212').complete_return_words(Word('1212'))
            set([word: 121212])

        REFERENCES:

        -   [1] J. Justin, L. Vuillon, Return words in Sturmian and
            episturmian words, Theor. Inform. Appl. 34 (2000) 343--356.
        """
        return set(self.complete_return_words_iterator(fact))

    def return_words_derivate(self, fact):
        r"""
        Returns the word generated by mapping a letter to each occurrence of
        the return words for the given factor dropping any dangling prefix and
        suffix. See for instance [1].

        EXAMPLES::

            sage: Word('12131221312313122').return_words_derivate(Word('1'))
            word: 123242

        REFERENCES:

        -   [1] F. Durand, A characterization of substitutive sequences using
            return words, Discrete Math. 179 (1998) 89--101.
        """
        idx = 0
        tab = {}
        ret = map(lambda w: tab.setdefault(w, len(tab)) + 1, \
                                self._return_words_list(fact))
        from sage.combinat.words.word import Word
        return Word(ret)

    def is_quasiperiodic(self):
        r"""
        Returns True if self is quasiperiodic, and False otherwise.

        A finite or infinite word `w` is *quasiperiodic* if it can be
        constructed by concatenations and superpositions of one of its proper
        factors `u`, which is called a *quasiperiod* of `w`.
        See for instance [1], [2], and [3].

        EXAMPLES::

            sage: Word('abaababaabaababaaba').is_quasiperiodic()
            True
            sage: Word('abacaba').is_quasiperiodic()
            False
            sage: Word('a').is_quasiperiodic()
            False
            sage: Word().is_quasiperiodic()
            False
            sage: Word('abaaba').is_quasiperiodic()
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

    def quasiperiods(self):
        r"""
        Returns the quasiperiods of self as a list ordered from shortest to
        longest.

        Let `w` be a finite or infinite word. A *quasiperiod* of `w` is a
        proper factor `u` of `w` such that the occurrences of `u` in `w`
        entirely cover `w`, i.e., every position of `w` falls within some
        occurrence of `u` in `w`. See for instance [1], [2], and [3].

        EXAMPLES::

            sage: Word('abaababaabaababaaba').quasiperiods()
            [word: aba, word: abaaba, word: abaababaaba]
            sage: Word('abaaba').quasiperiods()
            [word: aba]
            sage: Word('abacaba').quasiperiods()
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

        ::

            sage: f = Word('1213121').evaluation_dict() # keys appear in random order
            {'1': 4, '2': 2, '3': 1}

        TESTS::

            sage: f = Word('1213121').evaluation_dict()
            sage: f['1'] == 4
            True
            sage: f['2'] == 2
            True
            sage: f['3'] == 1
            True
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
        from sage.combinat.partition import Partition
        if 0 in p:
            return Partition(p[:p.index(0)])
        else:
            return Partition(p)

    def overlap_partition(self, other, delay=0, p=None, involution=None) :
        r"""
        Returns the partition of the alphabet induced by the overlap of
        self and other with the given delay.

        The partition of the alphabet is given by the equivalence
        relation obtained from the symmetric, reflexive and transitive
        closure of the set of pairs of letters
        `R_{u,v,d} = \{ (u_k, v_{k-d}) : 0 \leq k < n, 0\leq k-d < m \}`
        where `u = u_0 u_1 \cdots u_{n-1}`, `v = v_0v_1\cdots v_{m-1}` are
        two words on the alphabet `A` and `d` is an integer.

        The equivalence relation defined by `R` is inspired from [1].

        INPUT:

        -  ``other`` - word on the same alphabet as self
        -  ``delay`` - integer
        -  ``p`` - disjoint sets data structure (optional, default: None),
           a partition of the alphabet into disjoint sets to start with.
           If None, each letter start in distinct equivalence classes.
        -  ``involution`` - callable (optional, default: None), an
           involution on the alphabet. If involution is not None, the relation
           `R_{u,v,d} \cup R_{involution(u),involution(v),d}` is considered.

        OUTPUT:

        -  disjoint set data structure

        EXAMPLES::

            sage: W = Words(list('abc')+range(6))
            sage: u = W('abc')
            sage: v = W(range(5))
            sage: u.overlap_partition(v)
            {{0, 'a'}, {1, 'b'}, {2, 'c'}, {3}, {4}, {5}}
            sage: u.overlap_partition(v, 2)
            {{'a'}, {'b'}, {0, 'c'}, {1}, {2}, {3}, {4}, {5}}
            sage: u.overlap_partition(v, -1)
            {{0}, {1, 'a'}, {2, 'b'}, {3, 'c'}, {4}, {5}}

        You can re-use the same disjoint set and do more than one overlap::

            sage: p = u.overlap_partition(v, 2)
            sage: p
            {{'a'}, {'b'}, {0, 'c'}, {1}, {2}, {3}, {4}, {5}}
            sage: u.overlap_partition(v, 1, p)
            {{'a'}, {0, 1, 'b', 'c'}, {2}, {3}, {4}, {5}}

        The function  ``overlap_partition`` can be used to study equations
        on words. For example, if a word `w` overlaps itself with delay `d`, then
        `d` is a period of `w`::

            sage: W = Words(range(20))
            sage: w = W(range(14)); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
            sage: d = 5
            sage: p = w.overlap_partition(w, d)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 56789567895678
            sage: w2.minimal_period() == d
            True

        If a word is equal to its reversal, then it is a palindrome::

            sage: W = Words(range(20))
            sage: w = W(range(17)); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
            sage: p = w.overlap_partition(w.reversal(), 0)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 01234567876543210
            sage: w2.parent()
            Words over {0, 1, 2, 3, 4, 5, 6, 7, 8, 17, 18, 19}
            sage: w2.is_palindrome()
            True

        If the reversal of a word `w` is factor of its square `w^2`, then
        `w` is symmetric, i.e. the product of two palindromes::

            sage: W = Words(range(10))
            sage: w = W(range(10)); w
            word: 0123456789
            sage: p = (w*w).overlap_partition(w.reversal(), 4)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 0110456654
            sage: w2.is_symmetric()
            True

        If the image of the reversal of a word `w` under an involution `f`
        is factor of its square `w^2`, then `w` is `f`-symmetric::

            sage: W = Words([-11,-9,..,11])
            sage: w = W([1,3,..,11])
            sage: w
            word: 1,3,5,7,9,11
            sage: inv = lambda x:-x
            sage: f = WordMorphism(dict( (a, inv(a)) for a in W.alphabet()))
            sage: p = (w*w).overlap_partition(f(w).reversal(), 2, involution=f)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: m(w)
            word: 1,-1,5,7,-7,-5
            sage: m(w).is_symmetric(f)
            True

        TESTS::

            sage: W = Words('abcdef')
            sage: w = W('abc')
            sage: y = W('def')
            sage: w.overlap_partition(y, -3)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, -2)
            {{'a', 'f'}, {'b'}, {'c'}, {'d'}, {'e'}}
            sage: w.overlap_partition(y, -1)
            {{'a', 'e'}, {'b', 'f'}, {'c'}, {'d'}}
            sage: w.overlap_partition(y, 0)
            {{'a', 'd'}, {'b', 'e'}, {'c', 'f'}}
            sage: w.overlap_partition(y, 1)
            {{'a'}, {'b', 'd'}, {'c', 'e'}, {'f'}}
            sage: w.overlap_partition(y, 2)
            {{'a'}, {'b'}, {'c', 'd'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, 3)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, 4)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}

        ::

            sage: W = Words(range(2))
            sage: w = W([0,1,0,1,0,1]); w
            word: 010101
            sage: w.overlap_partition(w, 0)
            {{0}, {1}}
            sage: w.overlap_partition(w, 1)
            {{0, 1}}

        ::

            sage: empty = Word()
            sage: empty.overlap_partition(empty, 'yo')
            Traceback (most recent call last):
            ...
            TypeError: delay (=yo) must be an integer
            sage: empty.overlap_partition(empty,2,'yo')
            Traceback (most recent call last):
            ...
            TypeError: p(=yo) is not a DisjointSet

        The involution input can be any callable::

            sage: w = Words([-5,..,5])([-5..5])
            sage: inv = lambda x:-x
            sage: w.overlap_partition(w, 2, involution=inv)
            {{-4, -2, 0, 2, 4}, {-5, -3, -1, 1, 3, 5}}

        REFERENCES:

        -   [1] S. Labbé, Propriétés combinatoires des `f`-palindromes,
            Mémoire de maîtrise en Mathématiques, Montréal, UQAM, 2008,
            109 pages.
        """
        if not isinstance(delay, (int, Integer)):
            raise TypeError("delay (=%s) must be an integer"%delay)
        elif delay < 0:
            return other.overlap_partition(self, -delay, p)

        from sage.sets.disjoint_set import DisjointSet_class
        if p is None:
            if self.parent().size_of_alphabet() is Infinity:
                raise ValueError("The alphabet of the parent must be finite")
            from sage.sets.disjoint_set import DisjointSet
            p = DisjointSet(self.parent().alphabet())
        elif not isinstance(p, DisjointSet_class):
            raise TypeError("p(=%s) is not a DisjointSet" % p)

        #Join the classes of each pair of letters that are one above the other
        from itertools import islice, izip
        from sage.combinat.words.morphism import WordMorphism
        S = izip(islice(self, delay, None), other)
        if involution is None:
            for (a,b) in S:
                p.union(a, b)
        elif isinstance(involution, WordMorphism):
            for (a,b) in S:
                p.union(a, b)
                # take the first letter of the word
                p.union(involution(a)[0], involution(b)[0])
        elif callable(involution):
            for (a,b) in S:
                p.union(a, b)
                p.union(involution(a), involution(b))
        else:
            raise TypeError("involution (=%s) must be callable"%involution)
        return p

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
            sage: filter(lambda q: q.length() <= p.length() and \
            ....:       q.inverse().action(w) == list(v), \
            ....:       Permutations(w.length()) )
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
        from sage.combinat.permutation import Permutation
        return Permutation(result)

    def _s(self, i):
        r"""
        Implements Lascoux and Schutzenberger's `s_i` operator, swapping the
        number of `i` and `i+1`s in a word.

        EXAMPLES::

            sage: w = Word([1,1,2,1,2,2,1,3,1,2])
            sage: w._s(1)
            word: 1221221312

        TESTS::

            sage: w = Word([])
            sage: w._s(1)
            word:
            sage: w = Words(3)([2,1,2])
            sage: w._s(1).parent()
            Words over {1, 2, 3}
        """
        unpaired_i  = [] # positions of unpaired is
        unpaired_ip = [] # positions of unpaired i+1s
        for p,x in enumerate(self):
            if x == i:
                if unpaired_ip:
                    unpaired_ip.pop()
                else:
                    unpaired_i.append(p)
            elif x == i+1:
                unpaired_ip.append(p)

        unpaired = unpaired_i + unpaired_ip

        out = list(self)

        # replace the unpaired subword i^a (i+1)^b
        # with i^b (i+1)^a

        for j,p in enumerate(unpaired):
            if j < len(unpaired_ip):
                out[p] = i
            else:
                out[p] = i+1
        return self.parent()(out)

    def _to_partition_content(self):
        r"""
        Returns the conversion of self to a word with partition content using
        the `s_i` operators of Lascoux and Schutzenberger.

        EXAMPLES:
            sage: w = Word([1,3,2,1,2,3,4,6,4,2,3,2])
            sage: w._to_partition_content()
            word: 132112454132
            sage: Word([])._to_partition_content()
            word:
        """
        if self.length() == 0:
            return self

        from sage.combinat.words.word import Word
        n = max(self)
        ev = Word( Words(n)(self).evaluation() )
        sig = ev.reversal().standard_permutation().reduced_word()

        # sig is now the reverse complement of a reduced word for a minimal
        # length permutation which would sort the letters of ev into a
        # partition

        out = self
        for i in reversed(sig):
            out = out._s(n-i)
        return out

    def cocharge(self):
        r"""
        Returns the cocharge of self.  For a word `w`, this can be defined as
        `n_{ev} - ch(w)`, where `ch(w)` is the charge of `w` and `ev` is the
        evaluation of `w`, and `n_{ev}` is `\sum_{i<j} min(ev_i, ev_j)`.

        EXAMPLES::

            sage: Word([1,2,3]).cocharge()
            0
            sage: Word([3,2,1]).cocharge()
            3
            sage: Word([1,1,2]).cocharge()
            0
            sage: Word([2,1,2]).cocharge()
            1

        TESTS::

            sage: Word([]).cocharge()
            0
        """
        return self.evaluation_partition().weighted_size() - self.charge()

    def charge(self, check=True):
        r"""
        Returns the charge of self.  This is defined as follows.

        If `w` is a permutation of length `n`, (in other words, the evaluation
        of `w` is `(1, 1, \dots, 1)`), the statistic charge(`w`) is given by
        `\sum_{i=1}^n c_i(w)` where `c_1(w) = 0` and `c_i(w)` is defined
        recursively by setting `p_i` equal to `1` if `i` appears to the right
        of `i-1` in `w` and `0` otherwise.  Then we set `c_i(w) = c_{i-1}(w) +
        p_i`.

        EXAMPLES::

            sage: Word([1, 2, 3]).charge()
            3
            sage: Word([3, 5, 1, 4, 2]).charge() == 0 + 1 + 1 + 2 + 2
            True

        If `w` is not a permutation, but the evaluation of `w` is a partition,
        the charge of `w` is defined to be the sum of its charge subwords
        (each of which will be a permutation).  The first charge subword is
        found by starting at the end of `w` and moving left until the first
        `1` is found.  This is marked, and we continue to move to the left
        until the first `2` is found, wrapping around from the beginning of
        the word back to the end, if necessary.  We mark this `2`, and
        continue on until we have marked the largest letter in `w`.  The
        marked letters, with relative order preserved, form the first charge
        subword of `w`.  This subword is removed, and the next charge subword
        is found in the same manner from the remaining letters.  In the
        following example, `w1, w2, w3` are the charge subwords of `w`.

        EXAMPLE::

            sage: w = Word([5,2,3,4,4,1,1,1,2,2,3])
            sage: w1 = Word([5, 2, 4, 1, 3])
            sage: w2 = Word([3, 4, 1, 2])
            sage: w3 = Word([1, 2])
            sage: w.charge() == w1.charge() + w2.charge() + w3.charge()
            True

        Finally, if `w` does not have partition content, we apply the
        Lascoux-Schutzenberger standardization operators `s_i` in such a
        manner as to obtain a word with partition content. (The word we obtain
        is independent of the choice of operators.)  The charge is then
        defined to be the charge of this word::

            sage: Word([3,3,2,1,1]).charge()
            0
            sage: Word([1,2,3,1,2]).charge()
            2

        Note that this differs from the definition of charge given in
        Macdonald's book.  The difference amounts to a choice of
        reading a word from left-to-right or right-to-left.  The choice in
        Sage was made to agree with the definition of a reading word of a
        tableau in Sage, and seems to be the more common convention in the
        literature.

        REFERENCES:

        [1] Ian Macdonald, *Symmetric Functions and Hall Polynomials* second
        edition, 1995, Oxford University Press

        [2] A. Lascoux, L. Lapointe, and J. Morse.  *Tableau atoms and a new
        Macdonald positivity conjecture.* Duke Math Journal, **116 (1)**,
        2003.  Available at: [http://arxiv.org/abs/math/0008073]

        [3] A. Lascoux, B. Leclerc, and J.Y. Thibon.  *The Plactic Monoid*.
        Survey article available at
        [http://www-igm.univ-mlv.fr/~jyt/ARTICLES/plactic.ps]

        TESTS::

            sage: Word([1,1,2,2,3]).charge()
            4
            sage: Word([3,1,1,2,2]).charge()
            3
            sage: Word([2,1,1,2,3]).charge()
            2
            sage: Word([2,1,1,3,2]).charge()
            2
            sage: Word([3,2,1,1,2]).charge()
            1
            sage: Word([2,2,1,1,3]).charge()
            1
            sage: Word([3,2,2,1,1]).charge()
            0
            sage: Word([]).charge()
            0
        """
        if check:
            ev_dict = self.evaluation_dict()
            ordered_alphabet = sorted(ev_dict, cmp=self.parent().cmp_letters)
            evaluation = [ev_dict[a] for a in ordered_alphabet]
            from sage.combinat.partition import Partitions
            if evaluation not in Partitions():
                return self._to_partition_content().charge()
        res = 0
        w = self.to_integer_list()
        while len(w) != 0:
            i =len(w) - 1
            l = min(w)
            index = 0
            while len(w) != 0 and l <= max(w):
                while w[i] != l:
                    i -= 1
                    if i < 0:
                        i = len(w) - 1
                        index += 1
                res += index
                l += 1
                w.pop(i)
                i -= 1
                if i < 0:
                    i = len(w) - 1
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

            sage: Word('abaccaaba').BWT()
            word: cbaabaaca
            sage: Word('abaab').BWT()
            word: bbaaa
            sage: Word('bbabbaca').BWT()
            word: cbbbbaaa
            sage: Word('aabaab').BWT()
            word: bbaaaa
            sage: Word().BWT()
            word:
            sage: Word('a').BWT()
            word: a

        REFERENCES:

        -   [1] M. Burrows, D.J. Wheeler, "A block-sorting lossless data
            compression algorithm", HP Lab Technical Report, 1994, available
            at http://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.html
        """
        if self.is_empty():
           return self
        conjugates = sorted(self._conjugates_list())
        return self.parent()([x[x.length()-1] for x in conjugates])

    def iterated_left_palindromic_closure(self, f=None):
        r"""
        Returns the iterated left (`f`-)palindromic closure of self.

        INPUT:

        -  ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
           It must be callable on letters as well as words (e.g. WordMorphism).

        OUTPUT:

        word -- the left iterated `f`-palindromic closure of ``self``.

        EXAMPLES::

            sage: Word('123').iterated_left_palindromic_closure()
            word: 3231323
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('ab').iterated_left_palindromic_closure(f=f)
            word: abbaab
            sage: Word('aab').iterated_left_palindromic_closure(f=f)
            word: abbaabbaab

        TESTS:

        If ``f`` is not a involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('aab').iterated_left_palindromic_closure(f=f)
            Traceback (most recent call last):
            ...
            TypeError: self (=a->b, b->b) is not an endomorphism

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

    def balance(self):
        r"""
        Returns the balance of self.

        The balance of a word is the smallest number `q` such that self is
        `q`-balanced [1].

        A finite or infinite word `w` is said to be `q`-*balanced* if for
        any two factors `u`, `v` of `w` of the same length, the difference
        between the number of `x`'s in each of `u` and `v` is at most `q`
        for all letters `x` in the alphabet of `w`. A `1`-balanced word is
        simply said to be balanced. See Chapter 2 of [2].

        OUTPUT:

        integer

        EXAMPLES::

            sage: Word('1111111').balance()
            0
            sage: Word('001010101011').balance()
            2
            sage: Word('0101010101').balance()
            1

        ::

            sage: w = Word('11112222')
            sage: w.is_balanced(2)
            False
            sage: w.is_balanced(3)
            False
            sage: w.is_balanced(4)
            True
            sage: w.is_balanced(5)
            True
            sage: w.balance()
            4

        TESTS::

            sage: Word('1111122222').balance()
            5
            sage: Word('').balance()
            0
            sage: Word('1').balance()
            0
            sage: Word('12').balance()
            1
            sage: Word('1112').balance()
            1

        REFERENCES:

        -  [1] I. Fagnot, L. Vuillon, Generalized balances in Sturmian words,
           Discrete Applied Mathematics 121 (2002), 83--101.
        -  [2] M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
           Encyclopedia of Mathematics and its Applications, Cambridge
           University Press, U.K., 2002.
        """
        alphabet = set(self)
        best = 0
        for i in range(1, self.length()):
            start = iter(self)
            end = iter(self)
            abelian = dict(zip(alphabet, [0]*len(alphabet)))
            for _ in range(i):
                abelian[end.next()] += 1
            abel_max = abelian.copy()
            abel_min = abelian.copy()
            for _ in range(self.length() - i):
                lost = start.next()
                gain = end.next()
                abelian[gain] += 1
                abelian[lost] -= 1
                abel_max[gain] = max(abel_max[gain], abelian[gain])
                abel_min[lost] = min(abel_min[lost], abelian[lost])
            best = max(best, max(abel_max[a] - abel_min[a] for a in alphabet))
        return best

    def is_balanced(self, q=1):
        r"""
        Returns True if self is `q`-balanced, and False otherwise.

        A finite or infinite word `w` is said to be `q`-*balanced* if for
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
            raise TypeError("the balance level must be a positive integer")
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


    def sturmian_desubstitute_as_possible(self):
        r"""
        Sturmian desubstitutes the word ``self`` as much as possible.

        The finite word ``self`` must be defined on a two-letter
        alphabet or use at most two-letters.

        It can be Sturmian desubstituted if one letter appears
        isolated: the Sturmian desubstitution consists in removing one
        letter per run of the non-isolated letter. The accelerated
        Sturmian desubstitution consists in removing a run equal to
        the length of the shortest inner run from any run of the
        non-isolated letter (including possible leading and trailing
        runs even if they have shorter length). The (accelerated)
        Sturmian desubstitution is done as much as possible. A word is
        a factor of a Sturmian word if, and only if, the result is the
        empty word.

        OUTPUT:

        - A finite word defined on a two-letter alphabet.

        EXAMPLES::

            sage: u = Word('10111101101110111',alphabet='01') ; u
            word: 10111101101110111
            sage: v = u.sturmian_desubstitute_as_possible() ; v
            word: 01100101
            sage: v == v.sturmian_desubstitute_as_possible()
            True

            sage: Word('azaazaaazaaazaazaaaz', alphabet='az').sturmian_desubstitute_as_possible()
            word:


        TESTS::

            sage: w = Word('azazaza', alphabet='aze')
            sage: w.sturmian_desubstitute_as_possible()
            word:
            sage: Word('aze').sturmian_desubstitute_as_possible()
            Traceback (most recent call last):
            ...
            TypeError: your word must be defined on a binary alphabet or use at most two different letters
            sage: Word('azaaazaazaazaaazaaza', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azaaazaazaazaaazaaaza', alphabet='az').sturmian_desubstitute_as_possible()
            word: azzaa

        Boundary effects::

            sage: Word('', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azzzzz', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('zzzzza', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaazaaaaaaaaa', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaaaaaaaaaaa', alphabet='az').sturmian_desubstitute_as_possible()
            word:

        Boundary effects without alphabet::

            sage: Word('').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azzzzz').sturmian_desubstitute_as_possible()
            word:
            sage: Word('zzzzza').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaazaaaaaaaaa').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaaaaaaaaaaa').sturmian_desubstitute_as_possible()
            word:

        Idempotence::

            sage: r = words.RandomWord(randint(1,15)).sturmian_desubstitute_as_possible() ; r == r.sturmian_desubstitute_as_possible()
            True

        AUTHOR:

        -   Thierry Monteil
        """
        if self.is_empty():
            return self
        W = self.parent()
        if (W.size_of_alphabet() == 2):
            alphabet = W.alphabet()
        else:
            alphabet_as_set = set(self)
            if len(alphabet_as_set) > 2:
                raise TypeError('your word must be defined on a binary alphabet or use at most two different letters')
            elif len(alphabet_as_set) < 2:
                return W()
            else:
                alphabet = list(alphabet_as_set)
        word_from_letter = {l:W([l],datatype="list") for l in alphabet}
        is_prefix = True
        current_run_length = 0
        prefix_length = 0
        prefix_letter = self[0]
        is_isolated = {alphabet[0]:True,alphabet[1]:True}
        minimal_run = {alphabet[0]:Infinity,alphabet[1]:Infinity}
        maximal_run = {alphabet[0]:0,alphabet[1]:0}
        runs = {alphabet[0]:[],alphabet[1]:[]}
        for i in self:
            if is_prefix:
                if i == prefix_letter:
                    prefix_length += 1
                    if prefix_length > 1:
                        is_isolated[i] = False
                else:
                    is_prefix = False
                    current_run_length = 1
                    previous_letter = i
            else:
                if i == previous_letter:
                    current_run_length += 1
                    if current_run_length > 1:
                        is_isolated[i] = False
                else:
                    runs[previous_letter].append(current_run_length)
                    minimal_run[previous_letter] = min(minimal_run[previous_letter],current_run_length)
                    maximal_run[previous_letter] = max(maximal_run[previous_letter],current_run_length)
                    current_run_length = 1
                    previous_letter = i
        # at this point, previous_letter is the suffix letter and current_run_length is the suffix length
        if not (is_isolated[alphabet[0]] or is_isolated[alphabet[1]]):
            return self
        elif is_isolated[alphabet[0]] and is_isolated[alphabet[1]]:
            return W()
        else:
            if is_isolated[alphabet[0]]:
                l_isolated = alphabet[0] #the isolated letter
                l_running = alphabet[1] #the running letter (non-isolated)
            else:
                l_isolated = alphabet[1]
                l_running = alphabet[0]
            w_isolated = word_from_letter[l_isolated] #the word associated to the isolated letter
            w_running = word_from_letter[l_running] #the word associated to the running letter
            min_run = minimal_run[l_running]
            if (prefix_letter == l_isolated) or (prefix_length <= min_run):
                desubstitued_word = W()
            else:
                desubstitued_word = w_running ** (prefix_length - min_run)
            for i in runs[l_running]:
                desubstitued_word = desubstitued_word + w_isolated + w_running ** (i - min_run)
            if (current_run_length > 0):
                desubstitued_word = desubstitued_word + w_isolated
                if (previous_letter == l_running) and (current_run_length > min_run):
                    desubstitued_word = desubstitued_word + w_running ** (current_run_length - min_run)
            return desubstitued_word.sturmian_desubstitute_as_possible()


    def is_sturmian_factor(self):
        r"""
        Tells whether ``self`` is a factor of a Sturmian word.

        The finite word ``self`` must be defined on a two-letter alphabet.

        Equivalently, tells whether ``self`` is balanced. The
        advantage over the is_balanced method is that this one runs in
        linear time whereas is_balanced runs in quadratic time.

        OUTPUT:

        - boolean -- the result.

        EXAMPLES::

            sage: w = Word('0111011011011101101',alphabet='01')
            sage: w.is_sturmian_factor()
            True

        ::

            sage: words.LowerMechanicalWord(random(),alphabet='01')[:100].is_sturmian_factor()
            True
            sage: words.CharacteristicSturmianWord(random())[:100].is_sturmian_factor()
            True

        ::

            sage: w = Word('aabb',alphabet='ab')
            sage: w.is_sturmian_factor()
            False

            sage: s1 = WordMorphism('a->ab,b->b')
            sage: s2 = WordMorphism('a->ba,b->b')
            sage: s3 = WordMorphism('a->a,b->ba')
            sage: s4 = WordMorphism('a->a,b->ab')
            sage: W = Words('ab')
            sage: w = W('ab')
            sage: for i in xrange(8): w = choice([s1,s2,s3,s4])(w)
            sage: w
            word: abaaabaaabaabaaabaaabaabaaabaabaaabaaaba...
            sage: w.is_sturmian_factor()
            True

        Famous words::

            sage: words.FibonacciWord()[:100].is_sturmian_factor()
            True
            sage: words.ThueMorseWord()[:1000].is_sturmian_factor()
            False
            sage: words.KolakoskiWord()[:1000].is_sturmian_factor()
            False

        REFERENCES:

        .. [Arn2002] P. Arnoux, Sturmian sequences, in Substitutions in Dynamics,
           N. Pytheas Fogg (Ed.), Arithmetics, and Combinatorics (Lecture
           Notes in Mathematics, Vol. 1794), 2002.
        .. [Ser1985] C. Series. The geometry of Markoff numbers. The Mathematical
           Intelligencer, 7(3):20–29, 1985.
        .. [SU2009] J. Smillie and C. Ulcigrai. Symbolic coding for linear
           trajectories in the regular octagon, :arxiv:`0905.0871`, 2009.

        AUTHOR:

        -   Thierry Monteil
        """
        return self.sturmian_desubstitute_as_possible().is_empty()


    def is_tangent(self):
        r"""
        Tells whether ``self`` is a tangent word.

        The finite word ``self`` must be defined on a two-letter alphabet.

        A binary word is said to be *tangent* if it can appear in
        infintely many cutting sequences of a smooth curve, where each
        cutting sequence is observed on a progressively smaller grid.

        This class of words strictly contains the class of 1-balanced
        words, and is strictly contained in the class of 2-balanced words.

        This method runs in linear time.

        OUTPUT:

        - boolean -- the result.

        EXAMPLES::

            sage: w = Word('01110110110111011101',alphabet='01')
            sage: w.is_tangent()
            True

        Some tangent words may not be balanced::

            sage: Word('aabb',alphabet='ab').is_balanced()
            False
            sage: Word('aabb',alphabet='ab').is_tangent()
            True

        Some 2-balanced words may not be tangent::

            sage: Word('aaabb',alphabet='ab').is_tangent()
            False
            sage: Word('aaabb',alphabet='ab').is_balanced(2)
            True

        Famous words::

            sage: words.FibonacciWord()[:100].is_tangent()
            True
            sage: words.ThueMorseWord()[:1000].is_tangent()
            True
            sage: words.KolakoskiWord()[:1000].is_tangent()
            False

        REFERENCES:

        .. [Mon2010] T. Monteil, The asymptotic language of smooth curves, talk
           at LaCIM2010.

        AUTHOR:

        -   Thierry Monteil
        """
        if (self.parent().size_of_alphabet() != 2):
            raise TypeError('your word must be defined on a binary alphabet')
        [a,b] = self.parent().alphabet()
        mini = 0
        maxi = 0
        height = 0
        for i in self.sturmian_desubstitute_as_possible():
            if i == a:
                height = height + 1
                maxi = max(maxi , height)
            if i == b:
                height = height - 1
                mini = min(mini , height)
        return (maxi - mini <= 2)


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
        if j is None:
            j = i+1
        new = list(self)
        (new[i], new[j]) = (new[j], new[i])
        from sage.combinat.words.word import Word
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
            raise TypeError("the alphabet is infinite; specify a finite alphabet or use evaluation_dict() instead")
        if alphabet is None:
            alphabet = self._parent._alphabet
        ev_dict = self.evaluation_dict()
        return [ev_dict.get(a,0) for a in alphabet]

    evaluation = parikh_vector

    def robinson_schensted(self):
        """
        Return the semistandard tableau and standard tableau pair
        obtained by running the Robinson-Schensted algorithm on ``self``.

        This can also be done by running
        :func:`~sage.combinat.rsk.RSK` on ``self``.

        EXAMPLES::

            sage: Word([1,1,3,1,2,3,1]).robinson_schensted()
            [[[1, 1, 1, 1, 3], [2], [3]], [[1, 2, 3, 5, 6], [4], [7]]]
        """
        from sage.combinat.rsk import RSK
        return RSK(self)

    def _rsk_iter(self):
        r"""
        An iterator for :func:`~sage.combinat.rsk.RSK`.

        Yields pairs `(i, w_i)` for a word `w = w_1 w_2 \cdots w_k`.

        EXAMPLES::

            sage: for x in Word([1,1,3,1,2,3,1])._rsk_iter(): x
            ...
            (1, 1)
            (2, 1)
            (3, 3)
            (4, 1)
            (5, 2)
            (6, 3)
            (7, 1)
        """
        from itertools import izip
        return izip(xrange(1, len(self)+1), self)

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
                raise ValueError("for a nonzero overlap, words must contain integers as letters")
            if overlap is True:
                from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping
                return ShuffleProduct_overlapping(self, other)
            elif isinstance(overlap, (int,Integer)):
                from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping_r
                return ShuffleProduct_overlapping_r(self, other, overlap)
            raise ValueError('overlapping must be True or an integer')

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
            raise ValueError("for shifted shuffle, words must only contain integers as letters")
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
                raise ValueError("starting letter not in alphabet")
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
            if x not in seen:
                res.append(x)
                seen[x] = True
        return res

    def standard_factorization(self):
        r"""
        Returns the standard factorization of ``self``.

        The *standard factorization* of a word `w` of length greater than 1 is
        the unique factorization: `w = uv` where `v` is the longest proper
        suffix of `w` that is a Lyndon word.

        Note that if `w` is a Lyndon word of length greater than 1 with
        standard factorization `w = uv`, then `u` and `v` are also Lyndon words
        and `u < v`.

        See for instance [1], [2] and [3].

        INPUT:

        - ``self`` - finite word of length greater than 1

        OUTPUT:

            tuple -- tuple of two factors

        EXAMPLES::

            sage: Words('01')('0010110011').standard_factorization()
            (word: 001011, word: 0011)
            sage: Words('123')('1223312').standard_factorization()
            (word: 12233, word: 12)
            sage: Word([3,2,1]).standard_factorization()
            (word: 32, word: 1)

        ::

            sage: w = Word('0010110011',alphabet='01')
            sage: w.standard_factorization()
            (word: 001011, word: 0011)
            sage: w = Word('0010110011',alphabet='10')
            sage: w.standard_factorization()
            (word: 001011001, word: 1)
            sage: w = Word('1223312',alphabet='123')
            sage: w.standard_factorization()
            (word: 12233, word: 12)

        TESTS::

            sage: w = Word()
            sage: w.standard_factorization()
            Traceback (most recent call last):
            ...
            ValueError: Standard factorization not defined on words of
            length less than 2
            sage: w = Word('a')
            sage: w.standard_factorization()
            Traceback (most recent call last):
            ...
            ValueError: Standard factorization not defined on words of
            length less than 2

        REFERENCES:

        -   [1] K.-T. Chen, R.H. Fox, R.C. Lyndon, Free differential calculus,
            IV. The quotient groups of the lower central series, Ann. of Math.
            68 (1958) 81--95.
        -   [2] J.-P. Duval, Factorizing words over an ordered alphabet,
            J. Algorithms 4 (1983) 363--381.
        -   [3] M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
            Encyclopedia of Mathematics and its Applications, Cambridge
            University Press, U.K., 2002.
        """
        if self.length() < 2:
            raise ValueError("Standard factorization not defined on"
                             " words of length less than 2")
        for l in xrange(1, self.length()):
            suff = self[l:]
            if suff.is_lyndon():
                return self[:l], suff

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
        from sage.combinat.permutation import Permutation
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        if not isinstance(permutation, Permutation):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.domain())
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
        from sage.combinat.permutation import Permutation
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        if not isinstance(permutation, Permutation):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.domain())
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
            raise RuntimeError("Color map %s not known"%cmap)

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
            sage: Word('1212').is_square()
            True
            sage: Word('1213').is_square()
            False
            sage: Word('12123').is_square()
            False
            sage: Word().is_square()
            True
        """
        if self.length() % 2 != 0:
            return False
        else:
            l = self.length() // 2
            return self[:l] == self[l:]

    def is_square_free(self):
        r"""
        Returns True if self does not contain squares, and False
        otherwise.

        EXAMPLES::

            sage: Word('12312').is_square_free()
            True
            sage: Word('31212').is_square_free()
            False
            sage: Word().is_square_free()
            True

        TESTS:

        We make sure that #8490 is fixed::

            sage: Word('11').is_square_free()
            False
            sage: Word('211').is_square_free()
            False
            sage: Word('3211').is_square_free()
            False
        """
        L = self.length()
        if L < 2:
            return True
        for start in xrange(0, L-1):
            for end in xrange(start+2, L+1, 2):
                if self[start:end].is_square():
                    return False
        return True

    def is_cube(self):
        r"""
        Returns True if self is a cube, and False otherwise.

        EXAMPLES::

            sage: Word('012012012').is_cube()
            True
            sage: Word('01010101').is_cube()
            False
            sage: Word().is_cube()
            True
            sage: Word('012012').is_cube()
            False
        """
        if self.length() % 3 != 0:
            return False
        l = self.length() // 3
        return self[:l] == self[l:2*l] == self[2*l:]

    def is_cube_free(self):
        r"""
        Returns True if self does not contain cubes, and False otherwise.

        EXAMPLES::

            sage: Word('12312').is_cube_free()
            True
            sage: Word('32221').is_cube_free()
            False
            sage: Word().is_cube_free()
            True

        TESTS:

        We make sure that #8490 is fixed::

            sage: Word('111').is_cube_free()
            False
            sage: Word('2111').is_cube_free()
            False
            sage: Word('32111').is_cube_free()
            False
        """
        L = self.length()
        if L < 3:
            return True
        for start in xrange(0, L - 2):
            for end in xrange(start+3, L+1, 3):
                if self[start:end].is_cube():
                    return False
        return True

    def to_monoid_element(self):
        """
        Return ``self`` as an element the free monoid with the same alphabet
        as ``self``.

        EXAMPLES::

            sage: w = Word('aabb')
            sage: w.to_monoid_element()
            a^2*b^2
            sage: W = Words('abc')
            sage: w = W(w)
            sage: w.to_monoid_element()
            a^2*b^2

        TESTS:

        Check that ``w == w.to_monoid_element().to_word()``::

            sage: all(w.to_monoid_element().to_word() == w for i in range(6) for w in Words('abc', i))
            True
        """
        from sage.monoids.free_monoid import FreeMonoid
        try:
            l = list(self.parent().alphabet())
        except AttributeError:
            l = list(set(self))
        M = FreeMonoid(len(l), l)
        return M(self)

#######################################################################

class CallableFromListOfWords(tuple):
    r"""
    A class to create a callable from a list of words. The concatenation of
    a list of words is obtained by creating a word from this callable.
    """
    def __new__(cls, words):
        r"""
        TESTS::

            sage: from sage.combinat.words.finite_word import CallableFromListOfWords
            sage: w,u,x = Word([1,2,3]),Word([4,5]),Word([6,7,8])
            sage: f = CallableFromListOfWords([w,u,x]); f
            (word: 123, word: 45, word: 678)
            sage: f == loads(dumps(f))
            True
        """
        l = []
        for w in words:
            from word_infinite_datatypes import WordDatatype_callable
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

            sage: from sage.combinat.words.finite_word import CallableFromListOfWords
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
        raise IndexError("index (=%s) out of range" % i)

class Factorization(list):
    r"""
    A list subclass having a nicer representation for factorization of words.

    TESTS::

        sage: f = sage.combinat.words.finite_word.Factorization()
        sage: f == loads(dumps(f))
        True
    """
    def __repr__(self):
        r"""
        Returns a string representation of the object.

        TESTS::

            sage: sage.combinat.words.finite_word.Factorization()
            ()
            sage: sage.combinat.words.finite_word.Factorization([Word('ab'), Word('ba')])
            (ab, ba)
        """
        return '(%s)' % ', '.join(w.string_rep() for w in self)

