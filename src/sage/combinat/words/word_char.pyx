r"""
Fast word datatype using an array of unsigned char.
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
include "cysignals/memory.pxi"
include "sage/data_structures/bitset.pxi"

cimport cython
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.rational cimport Rational
from libc.string cimport memcpy, memcmp
from sage.combinat.words.word_datatypes cimport WordDatatype

from cpython.number cimport PyIndex_Check, PyNumber_Check
from cpython.sequence cimport PySequence_Check
from cpython.slice cimport PySlice_GetIndicesEx

import itertools

# the maximum value of a size_t
cdef size_t SIZE_T_MAX = -(<size_t> 1)

def reversed_word_iterator(WordDatatype_char w):
    r"""
    This function exists only because it is not possible to use yield in the
    special method ``__reversed__``.

    EXAMPLES::

        sage: W = Words([0,1,2])
        sage: w = W([0,1,0,0,1,2])
        sage: list(reversed(w)) # indirect doctest
        [2, 1, 0, 0, 1, 0]
    """
    cdef ssize_t i
    for i in range(w._length-1, -1, -1):
        yield w._data[i]

cdef class WordDatatype_char(WordDatatype):
    r"""
    A Fast class for words represented by an array ``unsigned char *``.

    Currently, only handles letters in [0,255].
    """
    cdef unsigned char * _data
    cdef size_t _length

    # _master is a just a reference to another Python object in case the finite
    # word is just a slice of another one. But because Cython takes care of
    # Python attributes *before* the call to __dealloc__ we need to duplicate
    # the information.
    cdef WordDatatype_char _master
    cdef int _is_slice

    def __cinit__(self):
        r"""
        Initialization of C attributes

        TESTS::

            sage: Words([0,1])([])
            word:
        """
        self._data = NULL
        self._length = 0
        self._is_slice = 0

    def __init__(self, parent, data):
        r"""
        Constructor

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: W([0,1,2,3])
            word: 0123
            sage: W(iter([0,1,2,3]))
            word: 0123
        """
        self._parent = parent

        if not PySequence_Check(data):
            data = list(data)
        if data:
            self._set_data(data)

    @cython.boundscheck(False) # assume that indexing will not cause any IndexErrors
    @cython.wraparound(False)  # not check not correctly handle negative indices
    cdef _set_data(self, data):
        r"""
        set the attribute ._data and ._length from the sequence data
        (usually data is a word, a tuple or a list)
        """
        cdef size_t i
        self._length = len(data)
        self._data = <unsigned char *> sig_malloc(self._length * sizeof(unsigned char))
        if self._data == NULL:
            raise MemoryError

        for i in range(self._length):
            self._data[i] = data[i]

    def __dealloc__(self):
        r"""
        Deallocate memory only if self uses it own memory.

        Note that ``sig_free`` will not deallocate memory if self is the
        master of another word.
        """
        # it is strictly forbidden here to access _master here! (it will be set
        # to None most of the time)
        if self._is_slice == 0:
            sig_free(self._data)

    def __nonzero__(self):
        r"""
        Test whether the word is not empty.

        EXAMPLES::

            sage: W = Words([0,3,5])
            sage: bool(W([0,3,3,5]))
            True
            sage: bool(W([]))
            False
        """
        return self._length != 0

    def is_empty(self):
        r"""
        Return whether the word is empty.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,2,2]).is_empty()
            False
            sage: W([]).is_empty()
            True
        """
        return not self

    def __len__(self):
        r"""
        Return the length of the word as a Python integer.

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: w = W([0,1,2,0,3,2,1])
            sage: len(w)
            7
            sage: type(len(w))
            <type 'int'>
        """
        return self._length

    def length(self):
        r"""
        Return the length of the word as a Sage integer.

        EXAMPLES::

            sage: W = Words([0,1,2,3,4])
            sage: w = W([0,1,2,0,3,2,1])
            sage: w.length()
            7
            sage: type(w.length())
            <type 'sage.rings.integer.Integer'>
            sage: type(len(w))
            <type 'int'>
        """
        return smallInteger(self._length)

    def letters(self):
        r"""
        Return the list of letters that appear in this word, listed in the
        order of first appearance.

        EXAMPLES::

            sage: W = Words(5)
            sage: W([1,3,1,2,2,3,1]).letters()
            [1, 3, 2]
        """
        cdef bitset_t seen
        bitset_init(seen, 256) # allocation + initialization to 0

        cdef size_t i
        cdef list res = []
        cdef unsigned char letter
        for i in range(self._length):
            letter = self._data[i]
            if not bitset_in(seen, letter):
                bitset_add(seen, letter)
                res.append(letter)
        bitset_free(seen)
        return res

    cdef _new_c(self, unsigned char * data, size_t length, WordDatatype_char master):
        r"""
        TO DISCUSS: in Integer (sage.rings.integer) this method is actually an
        external function. But we might want to have several possible inheritance.
        """
        cdef type t = type(self)
        cdef WordDatatype_char other = t.__new__(t)
        other._data = data
        other._master = master # can be None
        other._is_slice = 0 if master is None else 1
        other._length = length
        other._parent = self._parent

        return other

    def __hash__(self):
        r"""
        Return the hash value.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,0,1,0,0,0], datatype='list').__hash__()
            102060647
            sage: W([0,1,0,1,0,0,0], datatype='char').__hash__()
            102060647
        """
        cdef int res = 5381
        cdef size_t i
        if self._hash is None:
            for i in range(min(1024,self._length)):
                res = ((res << 5) + res) + self._data[i]
            self._hash = res
        return self._hash

    def __richcmp__(self, other, op):
        r"""
        INPUT:

        - ``other`` -- a word (WordDatatype_char)
        - ``op`` -- int, from 0 to 5

        TESTS::

            sage: W = Words(range(100))
            sage: w = W(range(10) * 2)
            sage: w == w
            True
            sage: w != w
            False
            sage: w[:-1] != w[1:]
            True
            sage: w < w[1:] and w[1:] > w
            True
            sage: w > w[1:] or w[1:] < w
            False
        """
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=
        if not isinstance(other, WordDatatype_char):
            return NotImplemented

        # word of different lengths are not equal!
        if (op == 2 or op == 3) and (<WordDatatype_char> self)._length != (<WordDatatype_char> other)._length:
            return op == 3

        cdef int test = (<WordDatatype_char> self)._lexico_cmp(other)
        if test < 0:
            return op < 2 or op == 3
        elif test > 0:
            return op > 3
        else:
            return op == 1 or op == 2 or op == 5

    def __cmp__(self, other):
        r"""
        INPUT:

        - ``other`` -- a word (WordDatatype_char)

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: cmp(W([0,1,0]), W([0,1,0]))
            0
            sage: cmp(W([0,1,0,0]), W([0,1,1]))
            -1
        """
        if not isinstance(other, WordDatatype_char):
            return NotImplemented

        cdef int test = self._lexico_cmp(other)
        if test:
            return test
        return (<Py_ssize_t> self._length) - (<Py_ssize_t> (<WordDatatype_char> other)._length)

    cdef int _lexico_cmp(self, WordDatatype_char other) except -2:
        r"""
        Lexicographic comparison of self and other up to
        the letter at position min(len(self),len(other))
        """
        cdef size_t l = min(self._length, other._length)

        sig_on()
        cdef int test = memcmp(<void *> (<WordDatatype_char> self)._data,
                      <void *> (<WordDatatype_char> other)._data,
                      l * sizeof(unsigned char))
        sig_off()

        if test == 0:
            return 0
        if test < 0:
            return -1
        else:
            return 1

    def __getitem__(self, key):
        r"""
        INPUT:

        - ``key`` -- index

        TESTS::

            sage: W = Words([0,1,2,3])
            sage: w = W([0,1,0,2,0,3,1,2,3])
            sage: w[0]
            0
            sage: w[2]
            0
            sage: w[1:]
            word: 10203123
            sage: w[5::-2]
            word: 321

            sage: w = W([randint(0,3) for _ in range(20)])
            sage: list(w) == [w[i] for i in range(len(w))]
            True

            sage: w['foo':'bar']
            Traceback (most recent call last):
            ...
            TypeError: slice indices must be integers or None or have an __index__ method

        Check a weird behavior of PySlice_GetIndicesEx (:trac:`17056`)::

            sage: w[1:0]
            word:
        """
        cdef Py_ssize_t i, start, stop, step, slicelength
        cdef unsigned char * data
        cdef size_t j,k
        if isinstance(key, slice):
            # here the key is a slice
            PySlice_GetIndicesEx(key,
                    self._length,
                    &start, &stop, &step,
                    &slicelength) 
            if slicelength == 0:
                return self._new_c(NULL, 0, None)
            if step == 1:
                return self._new_c(self._data+start, stop-start, self)
            data = <unsigned char *> sig_malloc(slicelength * sizeof(unsigned char))
            j = 0
            for k in range(start,stop,step):
                data[j] = self._data[k]
                j += 1
            return self._new_c(data, slicelength, None)

        elif PyIndex_Check(key):
            # here the key is an int
            i = key    # cast key into a size_t
            if i < 0:
                i += self._length;
            if i < 0 or i >= self._length:
                raise IndexError("word index out of range")
            return self._data[i]

        raise TypeError("word indices must be integers")

    def __iter__(self):
        r"""
        Iterator over the letter of self

        EXAMPLES::

            sage: W = Words([0,1,2,3])
            sage: for i in W([0,0,1,0]):  # indirect doctest
            ....:     print i,
            0 0 1 0
        """
        cdef size_t i
        for i in range(self._length):
            yield self._data[i]

    def __reversed__(self):
        r"""
        Reversed iterator over the letter of self

        EXAMPLES::

            sage: W = Words([0,1,2,3])
            sage: list(reversed(W([0,0,1,0]))) # indirect doctest
            [0, 1, 0, 0]

        TESTS::

            sage: list(reversed(W([])))
            []
            sage: list(reversed(W([1])))
            [1]
        """
        return reversed_word_iterator(self)

    cdef _concatenate(self, WordDatatype_char other):
        cdef unsigned char * data
        data = <unsigned char *> sig_malloc((self._length + other._length) * sizeof(unsigned char))
        if data == NULL:
            raise MemoryError

        sig_on()
        memcpy(data, self._data, self._length * sizeof(unsigned char))
        memcpy(data+self._length, other._data, other._length * sizeof(unsigned char))
        sig_off()

        return self._new_c(data, self._length + other._length, None)

    def __mul__(self, other):
        r"""
        Concatenation of ``self`` and ``other``.

        TESTS:

            sage: W = Words(IntegerRange(0,255))
            sage: W([0,1]) * W([2,0])
            word: 0120

        The result is automatically converted to a WordDatatype_char. Currently we can
        even do::

            sage: w = W([0,1,2,3])
            sage: w * [4,0,4,0]
            word: 01234040
        """
        cdef WordDatatype_char w

        if isinstance(other, WordDatatype_char):
            return (<WordDatatype_char> self)._concatenate(other)

        elif PySequence_Check(other):
            # we convert other to a WordDatatype_char and perform the concatenation
            w = (<WordDatatype_char> self)._new_c(NULL, 0, None)
            w._set_data(other)
            return (<WordDatatype_char> self)._concatenate(w)

        raise TypeError("not able to initialize a word from {}".format(other))

    def __pow__(self, exp, mod):
        r"""
        Power

        INPUT:

        -  ``exp``  - an integer, a rational, a float number or plus infinity.

        TESTS::

            sage: W = Words(range(20))
            sage: w = W([0,1,2,3])
            sage: w
            word: 0123
            sage: w ** (1/2)
            word: 01
            sage: w ** 2
            word: 01230123
            sage: w ** 3
            word: 012301230123
            sage: w ** (7/2)
            word: 01230123012301
            sage: len(((w ** 2) ** 3) ** 5) == len(w) * 2 * 3 * 5
            True

        Infinite exponents::

            sage: W([0,1]) ** Infinity
            word: 0101010101010101010101010101010101010101...
        """
        if not PyNumber_Check(exp):
            raise ValueError("the exponent must be a number or infinity")
        if mod is not None:
            raise ValueError("a word can not be taken modulo")

        if exp == float('inf'):
            from sage.rings.infinity import Infinity
            fcn = lambda n: self[n % self.length()]
            return self._parent.shift()(fcn, datatype='callable')

        if exp < 0:
            raise ValueError("can not take negative power of a word")

        cdef WordDatatype_char w = self
        cdef size_t i, rest

        if type(exp) is Rational:
            if w._length % exp.denominator():
                raise ValueError("undefined")
            i = exp.floor()
            rest = (exp - exp.floor()) * w._length
        else:
            i = exp
            rest = 0

        # first handle the cases (i*length + rest) <= length and return the
        # corresponding prefix of self
        if i == 1 and rest == 0:
            return self
        if w._length == 0:
            return w._new_c(NULL, 0, None)
        if i == 0:
            if rest == 0:
                return w._new_c(NULL, 0, None)
            else:
                return w._new_c(w._data, rest, self)

        # now consider non trivial powers
        if w._length > SIZE_T_MAX / (i+1):
            raise OverflowError("the length of the result is too large")
        cdef size_t new_length = w._length * i + rest
        cdef unsigned char * data = <unsigned char *> sig_malloc(new_length * sizeof(unsigned char))
        if data == NULL:
            raise MemoryError

        cdef Py_ssize_t j = w._length
        memcpy(data, w._data, j * sizeof(unsigned char))
        while 2*j < new_length:
            memcpy(data + j, data, j * sizeof(unsigned char))
            j *= 2
        memcpy(data + j, data, (new_length - j) * sizeof(unsigned char))

        return w._new_c(data, new_length, None)

    def has_prefix(self, other):
        r"""
        Test whether ``other`` is a prefix of ``self``.

        INPUT:

        - ``other`` -- a word or a sequence (e.g. tuple, list)

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,1,0,1,2,0])
            sage: w.has_prefix([0,1,1])
            True
            sage: w.has_prefix([0,1,2])
            False
            sage: w.has_prefix(w)
            True
            sage: w.has_prefix(w[:-1])
            True
            sage: w.has_prefix(w[1:])
            False

        TESTS:

        :trac:`19322`::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,0,2])
            sage: w.has_prefix(words.FibonacciWord())
            False

            sage: w.has_prefix([0,1,0,2,0])
            False
            sage: w.has_prefix([0,1,0,2])
            True
            sage: w.has_prefix([0,1,0])
            True
        """
        cdef size_t i
        cdef WordDatatype_char w

        if isinstance(other, WordDatatype_char):
            # C level
            w = <WordDatatype_char> other
            if w._length > self._length:
                return False
            return memcmp(self._data, w._data, w._length * sizeof(unsigned char)) == 0

        elif PySequence_Check(other):
            # python level
            from sage.combinat.words.infinite_word import InfiniteWord_class
            if isinstance(other, InfiniteWord_class) or len(other) > len(self):
                return False

            for i in range(len(other)):
                if other[i] != self[i]:
                    return False
            return True

        raise TypeError("not able to initialize a word from {}".format(other))

    def is_square(self):
        r"""
        Returns True if self is a square, and False otherwise.

        EXAMPLES::

            sage: w = Word([n % 4 for n in range(48)], alphabet=[0,1,2,3])
            sage: w.is_square()
            True

        ::

            sage: w = Word([n % 4 for n in range(49)], alphabet=[0,1,2,3])
            sage: w.is_square()
            False
            sage: (w*w).is_square()
            True

        TESTS:

        The above tests correspond to the present class (char)::

            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_char'>

        ::

            sage: Word([], alphabet=[0,1]).is_square()
            True
            sage: Word([0], alphabet=[0,1]).is_square()
            False
            sage: Word([0,0], alphabet=[0,1]).is_square()
            True
        """
        cdef size_t l
        if self._length % 2 != 0:
            return False
        else:
            l = self._length // 2
            return memcmp(self._data, 
                          self._data + l, 
                          l * sizeof(unsigned char)) == 0

    def longest_common_prefix(self, other):
        r"""
        Return the longest common prefix of this word and ``other``.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,0,2]).longest_common_prefix([0,1])
            word: 01
            sage: u = W([0,1,0,0,1])
            sage: v = W([0,1,0,2])
            sage: u.longest_common_prefix(v)
            word: 010
            sage: v.longest_common_prefix(u)
            word: 010

        Using infinite words is also possible (and the return type is also a
        of the same type as ``self``)::

            sage: W([0,1,0,0]).longest_common_prefix(words.FibonacciWord())
            word: 0100
            sage: type(_)
            <class 'sage.combinat.words.word.FiniteWord_char'>

        An example of an intensive usage::

            sage: W = Words([0,1])
            sage: w = words.FibonacciWord()
            sage: w = W(list(w[:5000]))
            sage: L = [[len(w[n:].longest_common_prefix(w[n+fibonacci(i):]))
            ....:      for i in range(5,15)] for n in range(1,1000)]
            sage: for n,l in enumerate(L):
            ....:     if l.count(0) > 4: print n+1,l
            375 [0, 13, 0, 34, 0, 89, 0, 233, 0, 233]
            376 [0, 12, 0, 33, 0, 88, 0, 232, 0, 232]
            608 [8, 0, 21, 0, 55, 0, 144, 0, 377, 0]
            609 [7, 0, 20, 0, 54, 0, 143, 0, 376, 0]
            985 [0, 13, 0, 34, 0, 89, 0, 233, 0, 610]
            986 [0, 12, 0, 33, 0, 88, 0, 232, 0, 609]

        TESTS::

            sage: W = Words([0,1,2])
            sage: w = W([0,2,1,0,0,1])
            sage: w.longest_common_prefix(0)
            Traceback (most recent call last):
            ...
            TypeError: unsupported input 0
        """
        cdef WordDatatype_char w
        cdef size_t i
        cdef size_t m

        if isinstance(other, WordDatatype_char):
            # C level
            # (this can be much faster if we allow to compare larger memory
            # zones)
            w = other
            m = min(self._length, w._length)
            for i in range(m):
                if self._data[i] != w._data[i]:
                    break
            else:
                if self._length <= w._length:
                    return self
                else:
                    return other

            return self._new_c(self._data, i, self)

        elif PySequence_Check(other):
            # Python level
            # we avoid to call len(other) since it might be an infinite word
            for i,a in enumerate(itertools.islice(other, self._length)):
                if self._data[i] != a:
                    break
            else:
                i += 1

            return self._new_c(self._data, i, self)

        raise TypeError("unsupported input {}".format(other))

    def longest_common_suffix(self, other):
        r"""
        Return the longest common suffix between this word and ``other``.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W([0,1,0,2]).longest_common_suffix([2,0,2])
            word: 02
            sage: u = W([0,1,0,0,1])
            sage: v = W([1,2,0,0,1])
            sage: u.longest_common_suffix(v)
            word: 001
            sage: v.longest_common_suffix(u)
            word: 001

        TESTS::

            sage: W = Words([0,1,2])
            sage: w = W([0,2,1,0,0,1])
            sage: w.longest_common_suffix(0)
            Traceback (most recent call last):
            ...
            TypeError: unsupported input 0
        """
        cdef WordDatatype_char w
        cdef size_t i
        cdef size_t m
        cdef size_t lo

        if isinstance(other, WordDatatype_char):
            # C level
            # (this can be much faster if we could compare larger memory
            # zones)
            w = other
            m = min(self._length, w._length)
            for i in range(m):
                if self._data[self._length-i-1] != w._data[w._length-i-1]:
                    break
            else:
                if self._length <= w._length:
                    return self
                else:
                    return other

            return self._new_c(self._data+self._length-i, i, self)

        elif PySequence_Check(other):
            # Python level
            lo = len(other)
            m = min(self._length, lo)
            for i in range(m):
                if self._data[self._length-i-1] != other[lo-i-1]:
                    break
            else:
                if self._length == m:
                    return self
                else:
                    i += 1

            return self._new_c(self._data+self._length-i, i, self)

        raise TypeError("unsupported input {}".format(other))

