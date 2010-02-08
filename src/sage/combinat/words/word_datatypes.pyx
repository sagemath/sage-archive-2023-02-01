r"""
Datatypes for finite words
"""
#*****************************************************************************
#       Copyright (C) 2009 Franco Saliola <saliola@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2
#  (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject

cdef class WordDatatype(SageObject):
    r"""
    The generic WordDatatype class.

    Any word datatype must contain two attributes (at least)::

      - _parent
      - _hash

    They are automatically defined here and it's not necessary (and forbidden)
    to define them anywhere else.

    TESTS::

        sage: w = Word([0,1,1,0,0,1])
        sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype)
        True

    """
    cdef public _parent
    cdef public _hash
    pass

    def __reduce__(self):
        r"""
        Default pickle support

        TESTS::

            sage: w = Word([0,1,1,0,0,1])
            sage: w.__reduce__()
            (<class 'sage.combinat.words.word.FiniteWord_list'>, (Words, [0, 1, 1, 0, 0, 1]))
        """
        return self.__class__, (self._parent, self._data)

cdef class WordDatatype_list(WordDatatype):
    r"""
    Datatype class for words defined by tuples.
    """
    cdef public list _data

    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with a given parent.

        .. note::

           It is slower than WordDatatype_str and WordDatatype_tuple.

        INPUT:

        - ``parent`` - an instance of :class:`Words_all`
        - ``data`` - an iterable

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_list)
            True

        """
        self._parent = parent
        if isinstance(data,list):
            self._data = data
        else:
            self._data = list(data)
        self._hash = None

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` - anything

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: 0 in w
            True
            sage: 3 in w
            False

        """
        return a in self._data

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of self.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: list(iter(w))
            [0, 1, 1, 0]

        """
        return iter(self._data)

    def __len__(self):
        r"""
        Return the length of the word.

        .. note::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: len(w)
            4

        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: w.length()
            4

        """
        return len(self._data)

    def __getitem__(self, key):
        r"""
        Implements :method:``__getitem__`` for words stored as lists.

        INPUT:

        - ``key`` - integer

        EXAMPLES::

            sage: w = Word(range(100))
            sage: w[4]
            4
            sage: w[-1]
            99
            sage: w[3:10:2]
            word: 3579

        """
        if isinstance(key, slice):
            return self._parent.__call__(self._data.__getitem__(key))
        else:
            return self._data[key]

    def count(self, a):
        r"""
        Returns the number of occurrences of the letter ``a`` in the word
        ``self``.

        INPUT:

        -  ``a`` - a letter

        OUTPUT:

        - integer

        EXAMPLES::

            sage: w = Word([0,1,1,0,1])
            sage: w.count(0)
            2
            sage: w.count(1)
            3
            sage: w.count(2)
            0

        """
        return self._data.count(a)

cdef class WordDatatype_str(WordDatatype):
    """
    Datatype for words defined by strings.
    """
    cdef public str _data

    # TODO : allow initialization from non string data
    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with parent ``parent`` from the string ``data``.

        INPUT:

        - ``parent`` - instance of :class:`Words_all`
        - ``data`` - string

        EXAMPLES::

            sage: w = Word("abba")
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_str)
            True

        """
        self._parent = parent
        if isinstance(data, str):
            self._data = data
        else:
            self._data = "".join(map(str,data))
        self._hash = None

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of ``self``.

        EXAMPLES::

            sage: w = Word('abba')
            sage: list(iter(w))
            ['a', 'b', 'b', 'a']

        """
        return iter(self._data)

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` - anything

        EXAMPLES::

            sage: w = Word('abba')
            sage: 'a' in w
            True
            sage: 'c' in w
            False

        """
        # we need to override the non standard comportement of
        # the comportment of the __contains__ of python str
        if not isinstance(a, str):
            return False
        if len(a) != 1:
            return False
        else:
            return a in self._data

    cpdef _has_factor_naive(self, w):
        r"""
        A naive test for testing whether the word contains ``w`` as a factor.

        .. note::

           This just wraps Python's builtin :method:`__contains__` for :class:`str`.

        INPUT:

        - ``w`` - a word, or something that behaves like one (list, tuple, str, ...)

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: w = Word('abba')
            sage: w._has_factor_naive('ba')
            True
            sage: w._has_factor_naive('bab')
            False

        """
        if isinstance(w, WordDatatype_str):
            return w._data in self._data
        elif isinstance(w, str):
            return w in self._data
        raise ValueError

    cpdef find(self, letter, start=None, stop=None):
        """
        Return the position of the first occurrence of ``letter`` in
        ``self``.

        INPUT:

        -  ``letter`` - the letter to search for
        - ``start`` - [default: 0] the position in which to start searching
        - ``stop`` - [default: None] the position in which to stop searching

        OUTPUT:

        - nonnegative integer if the ``letter`` occurs in the word
        - ``-1`` if ``letter`` does not occur in the word. (This is the same
          behaviour as Python's :class:`str`.)

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.find("a")
            0
            sage: w.find("a", 4)
            5
            sage: w.find("a", 4, 5)
            -1

        """
        if start is None:
            start = 0
        if stop is None:
            stop = len(self._data)
        return self._data.find(letter, start, stop)

    def rfind(self, letter, start=None, stop=None):
        """
        Similar to :method:`find`, but searches through the word in
        reverse.

        INPUT:

        - ``letter`` - the letter to search for
        - ``start`` - [default: None] the position in which to start searching
        - ``stop`` - [default: None] the position in which to stop searching

        OUTPUT:

        - nonnegative integer if the ``letter`` occurs in the word
        - ``-1`` if ``letter`` does not occur in the word. (This is the same
          behaviour as Python's :class:`str`.)

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.rfind("a")
            12
            sage: w.rfind("a", 4, 8)
            6
            sage: w.rfind("a", 4, 5)
            -1

        """
        if len(letter) != 1:
           return False
        if start is None:
            start = 0
        if stop is None:
            stop = len(self._data)
        return self._data.rfind(letter, start, stop)

    def __len__(self):
        r"""
        Return the length of the word.

        .. note::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: len(w)
            13

        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.length()
            13

        """
        return len(self._data)

    def __getitem__(self, key):
        r"""
        Implements the :method:`__getitem__`.

        TESTS::

            sage: alphabet = map(chr, range(97,123))
            sage: w = Word(alphabet)
            sage: w[4]
            'e'
            sage: w[-1]
            'z'
            sage: w[3:10:2]
            word: dfhj
            sage: all(chr(i+97) == w[i] for i in range(w.length()))
            True

        """
        if isinstance(key, slice):
            return self._parent(self._data[key])
        return self._data[key]

    def count(self, letter):
        r"""
        Count the number of occurrences of ``letter``.

        INPUT:

        - ``letter`` - a letter

        OUTPUT:

        - integer

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.count('a')
            7
            sage: w.count('b')
            6
            sage: w.count('c')
            0

        """
        return self._data.count(letter)

    def split(self, sep=None, maxsplit=None):
        r"""
        Returns a list of words, using sep as a delimiter string.
        If maxsplit is given, at most maxsplit splits are done.

        See also the partition method.

        .. note::

           This just wraps Python's builtin :meth:`str::split` for :class:`str`.
        INPUT:

        - ``sep`` - a string or a Word

        - ``maxsplit`` - None or a positive integer

        OUTPUT:

        - a list of words

        EXAMPLES:

        You can split along white space to find words in a sentence::

            sage: w = Word("My tailor is poor")
            sage: w.split(" ")
            [word: My, word: tailor, word: is, word: poor]

        You can split in two words letters to get the length of blocks in the
        other letter::

            sage: w = Word("ababbabaaba")
            sage: w.split('a')
            [word: , word: b, word: bb, word: b, word: , word: b, word: ]
            sage: w.split('b')
            [word: a, word: a, word: , word: a, word: aa, word: a]

        You can split along words::

            sage: w = Word("3230301030323212323032321")
            sage: w.split("32")
            [word: , word: 30301030, word: , word: 12, word: 30, word: , word: 1]

        If the separator is not a string a ValueError is raised::

            sage: w = Word("le papa du papa du papa etait un petit pioupiou")
            sage: w.split(Word(['p','a','p','a']))
            Traceback (most recent call last):
            ...
            ValueError: the separator must be a string.
        """
        if isinstance(sep, str):
            pass
        elif isinstance(sep, WordDatatype_str):
            sep = sep._data
        else:
            raise ValueError, "the separator must be a string."

        if maxsplit is None:
            return map(self._parent, self._data.split(sep))
        else:
            return map(self._parent, self._data.split(sep,maxsplit))

    def partition(self, sep):
        r"""
        Search for the separator sep in S, and return the part before it, the
        separator itself, and the part after it. The concatenation of the terms
        in the list gives back the initial word.

        See also the split method.

        .. note::

           This just wraps Python's builtin :meth:`str::partition` for :class:`str`.
        INPUT:

        - ``sep`` - a string or a Word

        EXAMPLES::

            sage: w = Word("MyTailorIsPoor")
            sage: w.partition("Tailor")
            [word: My, word: Tailor, word: IsPoor]

        ::

            sage: w = Word("3230301030323212323032321210121232121010")
            sage: l = w.partition("323")
            sage: print l
            [word: , word: 323, word: 0301030323212323032321210121232121010]
            sage: sum(l, Word('')) == w
            True

        If the separator is not a string an error is raised::

            sage: w = Word("le papa du papa du papa etait un petit pioupiou")
            sage: w.partition(Word(['p','a','p','a']))
            Traceback (most recent call last):
            ...
            ValueError: the separator must be a string.
        """
        if isinstance(sep, str):
            return map(self._parent, self._data.partition(sep))
        elif isinstance(sep, WordDatatype_str):
            return map(self._parent, self._data.partition(sep._data))
        else:
            raise ValueError, "the separator must be a string."

    def is_suffix(self, other):
        r"""
        Test whether ``self`` is a suffix of ``other``.

        INPUT:

        - ``other`` - a word (an instance of :class:`Word_class`) or a
          :class:`str`.

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("ababa")
            sage: w.is_suffix(u)
            False
            sage: u.is_suffix(w)
            True
            sage: u.is_suffix("abbabaabababa")
            True

        TESTS::

            sage: w = Word("abbabaabababa")
            sage: u = Word(['a','b','a','b','a'])
            sage: w.is_suffix(u)
            False
            sage: u.is_suffix(w)
            True

        """
        if isinstance(other, WordDatatype_str):
            return other._data.endswith(self._data)
        elif isinstance(other, str):
            return other.endswith(self._data)
        else:
            return super(WordDatatype_str, self).is_suffix(other)

    def has_suffix(self, other):
        """
        Test whether ``self`` has ``other`` as a suffix.

        INPUT:

        - ``other`` - a word (an instance of :class:`Word_class`) or a
          :class:`str`.

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

        """
        if isinstance(other, WordDatatype_str):
            return self._data.endswith(other._data)
        elif isinstance(other, str):
            return self._data.endswith(other)
        else:
            return super(WordDatatype_str, self).has_suffix(other)

    def is_prefix(self, other):
        r"""
        Test whether ``self`` is a prefix of ``other``.

        INPUT:

        - ``other`` - a word (an instance of :class:`Word_class`) or a
          :class:`str`.

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("abbab")
            sage: w.is_prefix(u)
            False
            sage: u.is_prefix(w)
            True
            sage: u.is_prefix("abbabaabababa")
            True

        TESTS::

            sage: ab = Word('ab')
            sage: abba = Word(['a','b','b','a'])
            sage: ab.is_prefix(abba)
            True
            sage: abba.is_prefix(ab)
            False

        """
        if isinstance(other, WordDatatype_str):
            return other._data.startswith(self._data)
        if isinstance(other ,str):
            return other.startswith(self._data)
        else:
            return super(WordDatatype_str, self).is_prefix(other)

    def has_prefix(self, other):
        r"""
        Test whether ``self`` has ``other`` as a prefix.

        INPUT:

        - ``other`` - a word (an instance of :class:`Word_class`) or a
          :class:`str`.

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

        TESTS::

            sage: ab = Word('ab')
            sage: abba = Word(['a','b','b','a'])
            sage: ab.has_prefix(abba)
            False
            sage: abba.has_prefix(ab)
            True

        """
        if isinstance(other, WordDatatype_str):
            return self._data.startswith(other._data)
        if isinstance(other, str):
            return self._data.startswith(other)
        else:
            return super(WordDatatype_str, self).has_prefix(other)

cdef class WordDatatype_tuple(WordDatatype):
    r"""
    Datatype class for words defined by tuples.
    """
    cdef public tuple _data

    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with parent ``parent`` from an iterable ``data``.

        INPUT:

        - ``parent`` - instance of :class:`Words_all`
        - ``data`` - iterable

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_tuple)
            True
            sage: u = Word([0,1,1,0], datatype='tuple')
            sage: isinstance(u, sage.combinat.words.word_datatypes.WordDatatype_tuple)
            True

        """
        self._parent = parent
        if isinstance(data,tuple):
            self._data = data
        else:
            self._data = tuple(data)
        self._hash = None

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of self.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: list(iter(w))
            [0, 1, 1, 0]

        """
        return iter(self._data)

    def __len__(self):
        r"""
        Return the length of the word.

        .. note::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: len(w)
            4

        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: w.length()
            4

        """
        return len(self._data)

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` - anything

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: 0 in w
            True
            sage: 3 in w
            False

        """
        return a in self._data

    def __getitem__(self, key):
        r"""
        Implements ``__getitem__`` for words stored as tuples.

        INPUT:

        - ``key`` - an integer

        OUTPUT:

        - can be anything (an object contained in the word)

        EXAMPLES::

            sage: w = Word(tuple(range(100)))
            sage: w[4]
            4
            sage: w[-1]
            99
            sage: w[3:10:2]
            word: 3579
            sage: all(w[i] == i for i in range(100))
            True

        """
        if isinstance(key, slice):
            return self._parent(self._data[key])
        return self._data[key]
