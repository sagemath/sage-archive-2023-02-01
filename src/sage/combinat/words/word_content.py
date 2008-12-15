# coding=utf-8
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import imap, izip, count
from sage.rings.infinity import Infinity
from sage.combinat.words.utils import *

def BuildWordContent(obj, mapping=id_f, format=None, part=slice(None)):
    r"""
    Builds the content for a word.

    INPUT:
        obj -- a function, an iterator or a list, potentially with a length defined
        mapping -- function, a map sending elements of the iterable to
                   nonnegative integers.
        format -- string (default None), the explicit type of obj.  Can be
                  either 'empty', 'list', 'function' or 'iterator'.  If set
                  to None (the default), the type will be guessed from the
                  properties of obj.
        part --  slice (default slice(None)), the portion of the object to use

    OUTPUT:
        word content -- an object respecting the word content protocol wrapping obj

    TESTS:
        sage: from sage.combinat.words.word_content import BuildWordContent
        sage: from itertools import count, imap, repeat
        sage: len(BuildWordContent(None))
        0
        sage: len(BuildWordContent(None, format='empty'))
        0
        sage: c = BuildWordContent(['a', 'b'], format='empty')
        Traceback (most recent call last):
        ...
        TypeError: trying to build an empty word with something other than None
        sage: len(BuildWordContent(['0', '1', '1'], int))
        3
        sage: len(BuildWordContent(['0', '1', '1'], int, format='list'))
        3
        sage: BuildWordContent(10, format='list')
        Traceback (most recent call last):
        ...
        TypeError: trying to build a word backed by a list with a sequence not providing the required operations
        sage: c = BuildWordContent(lambda x: x%2)
        sage: len(BuildWordContent(lambda x: x%3, part=slice(0,10)))
        10
        sage: c = BuildWordContent(repeat(1))
        sage: len(BuildWordContent(count(), part=slice(10)))
        10
        sage: len(BuildWordContent(count(), part=slice(10, 0, -2)))
        5
    """
    if format is None:
        if isinstance(obj, WordContent):
            return obj[part]
        elif obj is None:
            format = 'empty'
        elif is_iterable(obj):
            if haslen(obj):
                format = 'list'
            else:
                format = 'iterator'
        elif callable(obj):
            format = 'function'
        else:
            pass # Don't set format if we don't know

    if format == 'empty':
        if obj is not None:
            raise TypeError, "trying to build an empty word with something other than None"
        return WordContentFromList((), mapping)[part]
    elif format == 'list' or format == 'string':
        if not (is_iterable(obj) and haslen(obj)):
            raise TypeError, "trying to build a word backed by a list with a sequence not providing the required operations"
        return WordContentFromList(obj, mapping)[part]
    elif format == 'function':
        if not callable(obj):
            raise TypeError, "trying to build a word backed by a function with a non-callable"
        return WordContentFromFunction(obj, mapping)[part]
    elif format == 'iterator':
        if not is_iterable(obj):
            raise TypeError, "trying to build a word backed by an iterator with a non-iterable"
        return WordContentFromIterator(obj, mapping)[part]
    raise NotImplementedError

def is_WordContent(obj):
    r"""
    Returns True if obj is the content of a word.

    EXAMPLES:
        sage: from sage.combinat.words.word_content import BuildWordContent, is_WordContent
        sage: is_WordContent(33)
        False
        sage: is_WordContent(BuildWordContent([0, 1, 0], sage.combinat.words.utils.id_f))
        True
    """
    return isinstance(obj, WordContent)

class WordContent(object):
    def __cmp__(self, other):
        r"""
        Compares two contents according to the data they contain.

        The types of the contents doesn't matter.

        TESTS:
            sage: c1 = sage.combinat.words.word_content.WordContentFromList([1, 2, 1, 1])
            sage: c2 = sage.combinat.words.word_content.WordContentFromList([1, 2, 3])
            sage: c3 = sage.combinat.words.word_content.WordContentFromList([1, 4])
            sage: c2.__cmp__(c1) > 0
            True
            sage: c2.__cmp__(c2) == 0
            True
            sage: c2.__cmp__(c3) < 0
            True
        """
        if not is_WordContent(other):
            return NotImplemented
        if haslen(self) and haslen(other):
            for (i, j) in izip(self, other):
                if i != j:
                    return i - j
            return len(self) - len(other)
        return NotImplemented

    def concatenate(self, other):
        r"""
        Method to concatenate two contents together.

        TESTS:
            sage: from sage.combinat.words.word_content import BuildWordContent
            sage: list(BuildWordContent([1]).concatenate(BuildWordContent([2, 3, 4])))
            [1, 2, 3, 4]
        """
        if not (haslen(self) and haslen(other)):
            raise ValueError, "cannot concatenate infinite words"
        return ConcatenateContent(self._get_list() + other._get_list())

    def _get_list(self):
        r"""
        Returns self as a list of contents suitable for concatenate

        TESTS:
            sage: type(sage.combinat.words.word_content.BuildWordContent([1, 2, 3])._get_list())
            <type 'tuple'>
        """
        return (self,)

    def _check_getitem_args(self, key):
        r"""
        Does all the type and range checking for the key argument to
        __getitem__().  Also takes care of normalizing the ranges specified
        by slicing to values suitable for xrange() or similar in the finite
        case.  In the infinite case a stop value of None, will stay as None
        since no end can be defined.

        TESTS:
            # Generic Cases
            sage: w = sage.combinat.words.word_content.WordContent()
            sage: w._check_getitem_args(1)
            1
            sage: w._check_getitem_args("abc")
            Traceback (most recent call last):
            ...
            TypeError: word indices must be integers
            sage: w._check_getitem_args(slice(1, 2, 3))
            slice(1, 2, 3)
            sage: w._check_getitem_args(slice("a"))
            Traceback (most recent call last):
            ...
            TypeError: word indices must be integers
            sage: w._check_getitem_args(slice("a", 1))
            Traceback (most recent call last):
            ...
            TypeError: word indices must be integers
            sage: w._check_getitem_args(slice(1, 1, "a"))
            Traceback (most recent call last):
            ...
            TypeError: word indices must be integers

            # Finite Cases
            sage: f = sage.combinat.words.word_content.WordContentFromList(range(10))
            sage: f._check_getitem_args(slice(None, None, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(None, 0, None))
            slice(0, 0, 1)
            sage: f._check_getitem_args(slice(None, 1, None))
            slice(0, 1, 1)
            sage: f._check_getitem_args(slice(None, -1, None))
            slice(0, 9, 1)
            sage: f._check_getitem_args(slice(None, 11, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(None, 10, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(None, 9, None))
            slice(0, 9, 1)
            sage: f._check_getitem_args(slice(None, -9, None))
            slice(0, 1, 1)
            sage: f._check_getitem_args(slice(None, -10, None))
            slice(0, 0, 1)
            sage: f._check_getitem_args(slice(None, -11, None))
            slice(0, 0, 1)
            sage: f._check_getitem_args(slice(0, None, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(1, None, None))
            slice(1, 10, 1)
            sage: f._check_getitem_args(slice(-1, None, None))
            slice(9, 10, 1)
            sage: f._check_getitem_args(slice(11, None, None))
            slice(10, 10, 1)
            sage: f._check_getitem_args(slice(10, None, None))
            slice(10, 10, 1)
            sage: f._check_getitem_args(slice(9, None, None))
            slice(9, 10, 1)
            sage: f._check_getitem_args(slice(-9, None, None))
            slice(1, 10, 1)
            sage: f._check_getitem_args(slice(-10, None, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(-11, None, None))
            slice(0, 10, 1)
            sage: f._check_getitem_args(slice(None, None, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(None, 0, -1))
            slice(9, 0, -1)
            sage: f._check_getitem_args(slice(None, 1, -1))
            slice(9, 1, -1)
            sage: f._check_getitem_args(slice(None, -1, -1))
            slice(9, 9, -1)
            sage: f._check_getitem_args(slice(None, 11, -1))
            slice(9, 9, -1)
            sage: f._check_getitem_args(slice(None, 10, -1))
            slice(9, 9, -1)
            sage: f._check_getitem_args(slice(None, 9, -1))
            slice(9, 9, -1)
            sage: f._check_getitem_args(slice(None, -10, -1))
            slice(9, 0, -1)
            sage: f._check_getitem_args(slice(None, -11, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(None, -12, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(0, None, -1))
            slice(0, -1, -1)
            sage: f._check_getitem_args(slice(1, None, -1))
            slice(1, -1, -1)
            sage: f._check_getitem_args(slice(-1, None, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(10, None, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(9, None, -1))
            slice(9, -1, -1)
            sage: f._check_getitem_args(slice(8, None, -1))
            slice(8, -1, -1)
            sage: f._check_getitem_args(slice(-10, None, -1))
            slice(0, -1, -1)
            sage: f._check_getitem_args(slice(-11, None, -1))
            slice(-1, -1, -1)
            sage: f._check_getitem_args(slice(-12, None, -1))
            slice(-1, -1, -1)
            sage: f._check_getitem_args(slice(None, None, 0))
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: f._check_getitem_args(9)
            9
            sage: f._check_getitem_args(10)
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
            sage: f._check_getitem_args(0)
            0
            sage: f._check_getitem_args(-1)
            9
            sage: f._check_getitem_args(-10)
            0
            sage: f._check_getitem_args(-11)
            Traceback (most recent call last):
            ...
            IndexError: word index out of range

            # Infinite Cases
            sage: from itertools import count
            sage: i = sage.combinat.words.word_content.WordContentFromIterator(count())
            sage: i._check_getitem_args(slice(None, None, None))
            slice(0, None, 1)
            sage: i._check_getitem_args(slice(None, 0, None))
            slice(0, 0, 1)
            sage: i._check_getitem_args(slice(None, 1, None))
            slice(0, 1, 1)
            sage: i._check_getitem_args(slice(None, -1, None))
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
            sage: i._check_getitem_args(slice(0, None, None))
            slice(0, None, 1)
            sage: i._check_getitem_args(slice(1, None, None))
            slice(1, None, 1)
            sage: i._check_getitem_args(slice(-1, None, None))
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
            sage: i._check_getitem_args(slice(None, None, -1))
            Traceback (most recent call last):
            ...
            ValueError: negative step and no defined start
            sage: i._check_getitem_args(slice(None, 0, -1))
            Traceback (most recent call last):
            ...
            ValueError: negative step and no defined start
            sage: i._check_getitem_args(slice(None, -1, -1))
            Traceback (most recent call last):
            ...
            ValueError: negative step and no defined start
            sage: i._check_getitem_args(slice(0, None, -1))
            slice(0, -1, -1)
            sage: i._check_getitem_args(slice(-1, None, -1))
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
            sage: i._check_getitem_args(slice(None, None, 0))
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: i._check_getitem_args(1)
            1
            sage: i._check_getitem_args(-1)
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
        """
        if slice_ok(key):
            start, stop, step = key.start, key.stop, key.step
            if step is None:
                step = 1
            if step == 0:
                raise ValueError, "slice step cannot be zero"
            if haslen(self):
                if start is None:
                    if step < 0:
                        start = len(self) - 1
                    else:
                        start = 0
                else:
                    if start < 0:
                        start += len(self)
                    if step < 0:
                        start = clamp(start, -1, len(self) - 1)
                    else:
                        start = clamp(start, 0, len(self))
                if stop is None:
                    if step < 0:
                        stop = -1
                    else:
                        stop = len(self)
                else:
                    if stop < 0:
                        stop += len(self)
                    if step < 0:
                        stop = clamp(stop, -1, len(self) - 1)
                    else:
                        stop = clamp(stop, 0, len(self))
            else:
                if start is None:
                    if step < 0:
                        raise ValueError, "negative step and no defined start"
                    else:
                        start = 0
                else:
                    start = self._check_getitem_args(start)
                if stop is None:
                    if step < 0:
                       stop = -1
                else:
                    stop = self._check_getitem_args(stop)
            return slice(start, stop, step)
        elif isint(key):
            key = int(key)
            if haslen(self):
                if key < 0:
                    key += len(self)
                if key < 0 or key >= len(self):
                    raise IndexError, "word index out of range"
            else:
                if key < 0:
                    raise IndexError, "negative index on infinite word"
            return key
        else:
            raise TypeError, "word indices must be integers"

class WordContentFromList(WordContent):
    def __init__(self, l, trans=id_f):
        r"""
        INPUT:
            l -- list
            trans -- function (default identity), defines a mapping between the objects in l and the nonnegative integers

        TESTS:
            sage: c = sage.combinat.words.word_content.WordContentFromList([0, 1, 1, 2, 1])
            sage: c == loads(dumps(c))
            True
        """
        # FIXME: the line below
        #typecode = {1:'B', 2:'I', 4:'L'}[dropwhile(lambda x: x < len(parent.alphabet())/256, [1, 2, 4]).next()]
        #from array import array
        #self._list = array('B', imap(trans, l))
        self._list = map(trans, l)

    def __iter__(self):
        r"""
        TESTS:
            sage: list(sage.combinat.words.word_content.WordContentFromList('012345', int)) # indirect test
            [0, 1, 2, 3, 4, 5]
        """
        return iter(self._list)

    def __len__(self):
        r"""
        TESTS:
            sage: len(sage.combinat.words.word_content.WordContentFromList([0, 1, 0, 0, 1]))
            5
        """
        return len(self._list)

    def __getitem__(self, key):
        r"""
        TESTS:
            sage: from sage.combinat.words.word_content import WordContentFromList
            sage: w = WordContentFromList('012345', int)
            sage: e = WordContentFromList('', int)
            sage: w__2 = WordContentFromList('45', int)
            sage: w_2 = WordContentFromList('01', int)
            sage: w__1 = WordContentFromList('01234', int)
            sage: w_s2 = WordContentFromList('024', int)
            sage: w_s_2 = WordContentFromList('531', int)
            sage: w[:] == w
            True
            sage: w[0:] == w
            True
            sage: w[10:] == e
            True
            sage: w[-2:] == w__2
            True
            sage: w[-10:] == w
            True
            sage: w[:0] == e
            True
            sage: w[:2] == w_2
            True
            sage: w[:10] == w
            True
            sage: w[:-1] == w__1
            True
            sage: w[:-10] == e
            True
            sage: w[::2] == w_s2
            True
            sage: w[::-2] == w_s_2
            True
            sage: w[::0]
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: w[0]
            0
            sage: w[5]
            5
            sage: w[6]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
            sage: w[-1]
            5
            sage: w[-6]
            0
            sage: w[-7]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
        """
        key = self._check_getitem_args(key)
        if isinstance(key, slice):
            if key.stop == -1:
                key = slice(key.start, None, key.step)
            return WordContentFromList(self._list[key])
        else:
            return self._list[key]

class WordContentFromFunction(WordContent):
    def __init__(self, func, trans=id_f):
        r"""
        INPUT:
            func -- callable
            trans -- callable (default: identity), defines a mapping between the objects returned by func and the nonnegative integers.

        NOTE:
            Unnamed functions do not pickle.  And if you arbitrarily
            slice this type of container then it may use lambdas
            and break pickling.

        TESTS:
            sage: c = sage.combinat.words.word_content.WordContentFromFunction(sage.combinat.words.utils.id_f)[:10]
            sage: c == loads(dumps(c))
            True
        """
        self._func = func
        self._len = Infinity
        self._trans = trans

    def __iter__(self):
        r"""
        TESTS:
            sage: list(sage.combinat.words.word_content.WordContentFromFunction(sage.combinat.words.utils.id_f)[:6]) # indirect test
            [0, 1, 2, 3, 4, 5]
        """
        for i in xrange(self._len) if haslen(self) else count():
            yield self._trans(self._func(i))

    def __len__(self):
        r"""
        TESTS:
            sage: c = sage.combinat.words.word_content.WordContentFromFunction(sage.combinat.words.utils.id_f)
            sage: len(c)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
            sage: c = c[:10]
            sage: len(c)
            10
        """
        return self._len

    def __getitem__(self, key):
        r"""
        TESTS:
            sage: from sage.combinat.words.utils import id_f
            sage: from sage.combinat.words.word_content import WordContentFromFunction, WordContentFromList
            sage: w = WordContentFromFunction(id_f);
            sage: e = WordContentFromList('', int)
            sage: w__2 = WordContentFromList('45', int)
            sage: w_2 = WordContentFromList('01', int)
            sage: w__1 = WordContentFromList('01234', int)
            sage: w_s2 = WordContentFromList('024', int)
            sage: w_s_2 = WordContentFromList('531', int)
            sage: w[:]          # random results
            <sage.combinat.words.word_content.WordContentFromFunction object at 0xe499c30>
            sage: w[0:]         # random results
            <sage.combinat.words.word_content.WordContentFromFunction object at 0xe499c50>
            sage: w[1:]         # random results
            <sage.combinat.words.word_content.WordContentFromFunction object at 0xe499c30>
            sage: w[:0] == e
            True
            sage: w[:5] == w__1
            True
            sage: w[::2]        # random results
            <sage.combinat.words.word_content.WordContentFromFunction object at 0xe499bf0>
            sage: w[::0]
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: w[5]
            5
            sage: w[-1]
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
            sage: w._len = 6
            sage: w[:] == w
            True
            sage: w[0:] == w
            True
            sage: w[10:] == e
            True
            sage: w[-2:] == w__2
            True
            sage: w[-10:] == w
            True
            sage: w[:0] == e
            True
            sage: w[:2] == w_2
            True
            sage: w[:10] == w
            True
            sage: w[:-1] == w__1
            True
            sage: w[:-10] == e
            True
            sage: w[::2] == w_s2
            True
            sage: w[::-2] == w_s_2
            True
            sage: w[::0]
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: w[0]
            0
            sage: w[5]
            5
            sage: w[6]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
            sage: w[-1]
            5
            sage: w[-6]
            0
            sage: w[-7]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
        """
        key = self._check_getitem_args(key)
        if isinstance(key, slice):
            if key.start != 0 or key.step != 1:
                # TODO: maybe use an accessory class for the slicing
                fn = lambda x: self._func(key.start + x*key.step)
            else:
                fn = self._func
            res = WordContentFromFunction(fn, self._trans)
            if key.stop is not None:
                res._len = int((key.stop - key.start)/key.step)
            return res
        else:
            return self._trans(self._func(key))

class WordContentFromIterator(WordContent):
    def __init__(self, it, trans=id_f):
        r"""
        INPUT:
            it -- iterator
            trans -- function (default identity), defines a mapping between the objects returned by it and the nonnegative integers.

        NOTE:
            It appears that islice does not pickle correctly causing
            various errors when reloading.  Also, most iterators
            do not support copying and should not support pickling by
            extension.

        TESTS:
            sage: from itertools import count
            sage: c = sage.combinat.words.word_content.WordContentFromIterator(count())[:10]

            #sage: c == loads(dumps(c)) # broken because of islice
            #True
        """
        self._iter = iter(it)
        self._len = Infinity
        self._trans = trans

    def __iter__(self):
        r"""
        TESTS:
            sage: from itertools import count
            sage: list(sage.combinat.words.word_content.WordContentFromIterator(count())[:6]) # indirect test
            [0, 1, 2, 3, 4, 5]
        """
        return imap(self._trans, self._get_it())

    def _get_it(self):
        r"""
        TESTS:
            sage: from itertools import count
            sage: c = sage.combinat.words.word_content.WordContentFromIterator(count())
            sage: it1 = c._get_it()
            sage: it2 = c._get_it()
            sage: it1.next()
            0
            sage: it2.next()
            0
            sage: list(c[:4])
            [0, 1, 2, 3]
            sage: it1.next()
            1
            sage: it2.next()
            1
        """
        self._iter, new = copy_it(self._iter)
        return new

    def __len__(self):
        r"""
        TESTS:
            sage: from itertools import count
            sage: c = sage.combinat.words.word_content.WordContentFromIterator(count())
            sage: len(c)
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
            sage: c = c[:10]
            sage: len(c)
            10
        """
        return self._len

    def __getitem__(self, key):
        r"""
        TESTS:
            sage: from itertools import count
            sage: from sage.combinat.words.word_content import WordContentFromIterator, WordContentFromList
            sage: w = WordContentFromIterator(count())
            sage: e = WordContentFromList('', int)
            sage: w__2 = WordContentFromList('45', int)
            sage: w_2 = WordContentFromList('01', int)
            sage: w__1 = WordContentFromList('01234', int)
            sage: w_s2 = WordContentFromList('024', int)
            sage: w_s_2 = WordContentFromList('531', int)
            sage: w[:]          # random results
            <sage.combinat.words.word_content.WordContentFromIterator object at 0xe499c30>
            sage: w[0:]         # random results
            <sage.combinat.words.word_content.WordContentFromIterator object at 0xe499c50>
            sage: w[1:]         # random results
            <sage.combinat.words.word_content.WordContentFromIterator object at 0xe499c30>
            sage: w[:0] == e
            True
            sage: w[:5] == w__1
            True
            sage: w[::2]        # random results
            <sage.combinat.words.word_content.WordContentFromIterator object at 0xe499bf0>
            sage: w[::0]
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: w[5]
            5
            sage: w[-1]
            Traceback (most recent call last):
            ...
            IndexError: negative index on infinite word
            sage: w._len = 6
            sage: w[:] == w
            True
            sage: w[0:] == w
            True
            sage: w[10:] == e
            True
            sage: w[-2:] == w__2
            True
            sage: w[-10:] == w
            True
            sage: w[:0] == e
            True
            sage: w[:2] == w_2
            True
            sage: w[:10] == w
            True
            sage: w[:-1] == w__1
            True
            sage: w[:-10] == e
            True
            sage: w[::2] == w_s2
            True
            sage: w[::-2] == w_s_2
            True
            sage: w[::0]
            Traceback (most recent call last):
            ...
            ValueError: slice step cannot be zero
            sage: w[0]
            0
            sage: w[5]
            5
            sage: w[6]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
            sage: w[-1]
            5
            sage: w[-6]
            0
            sage: w[-7]
            Traceback (most recent call last):
            ...
            IndexError: word index out of range
        """
        key = self._check_getitem_args(key)
        if isinstance(key, slice):
            if key.step < 0:
                if key.stop == -1:
                    key = slice(key.start, None, key.step)
                if haslen(self):
                    return WordContentFromList(slice_it(self, len(self), key))
                else:
                    return WordContentFromList(slice_it(self, \
                                int(key.start - key.stop), key))
            else:
                res = WordContentFromIterator(islice(self._get_it(), key.start, \
                                              key.stop, key.step), self._trans)
                if key.stop is not None:
                    res._len = int((key.stop - key.start)/key.step)
                return res
        else:
            return islice(self, key, None).next()

class ConcatenateContent(WordContentFromFunction):
    r"""
    Class that acts as a function to concatenate contents.
    """
    def __init__(self, l):
        r"""
        TESTS:
            sage: from sage.combinat.words.word_content import ConcatenateContent, BuildWordContent
            sage: c = ConcatenateContent((BuildWordContent([1, 2]),))
            sage: len(c)
            2
            sage: c = ConcatenateContent((BuildWordContent('1221', int), BuildWordContent('2112', int)))
            sage: len(c)
            8
        """
        self._list = l
        super(ConcatenateContent, self).__init__(self)
        self._len = sum(imap(len, self._list))

    def __iter__(self):
        r"""
        TESTS:
            sage: from sage.combinat.words.word_content import ConcatenateContent, BuildWordContent
            sage: list(ConcatenateContent((BuildWordContent('012', int), BuildWordContent([3, 4, 5])))) # indirect test
            [0, 1, 2, 3, 4, 5]
        """
        for c in self._list:
            for s in c:
                yield s

    def _get_list(self):
        r"""
        Optimization: for ConcatenateContents return the internal list.

        TESTS:
            sage: from sage.combinat.words.word_content import ConcatenateContent, BuildWordContent
            sage: type(ConcatenateContent((BuildWordContent([1, 2]),))._get_list())
            <type 'tuple'>
        """
        return self._list

    def __call__(self, key):
        r"""
        Returns the character at position key in the word.

        EXAMPLES:
            sage: from sage.combinat.words.word_content import ConcatenateContent, BuildWordContent
            sage: c = ConcatenateContent((BuildWordContent('1221', int), BuildWordContent('2112', int)))
            sage: c(0)
            1
            sage: c(7)
            2
        """
        for c in self._list:
            if (key - len(c) < 0):
                return c[key]
            key -= len(c)
