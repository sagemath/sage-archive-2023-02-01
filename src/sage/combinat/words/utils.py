"""
Word utilities
"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import islice, tee, izip, starmap, imap, ifilter, ifilterfalse
from copy import copy

def copy_it(it):
    r"""
    Copy an iterator using its builtin __copy__ method if
    available, otherwise use itertools.tee(). Define __copy__ for
    your iterators. (See PEP 323)

    TESTS::

        sage: from sage.combinat.words.utils import copy_it
        sage: it = iter([1, 2, 3, 4, 5])
        sage: it, it2 = copy_it(it)
        sage: list(it)
        [1, 2, 3, 4, 5]
        sage: it2.next()
        1
        sage: it2, it3 = copy_it(it2)
        sage: list(it2)
        [2, 3, 4, 5]
        sage: it3.next()
        2
    """
    it = iter(it)
    if hasattr(it, '__copy__'):
        return (it, copy(it))
    else:
        return tee(it)

def slice_it(it, l, key):
    r"""
    Slice an iterator, supporting negative step sizes by expliciting
    the elements if needed.

    NOTE: The iterator returned depends on it. You must pass in a copy
    of your iterator if you intend to keep using the original
    iterator.

    TESTS::

        sage: from sage.combinat.words.utils import slice_it
        sage: list(slice_it(iter(range(5)), 5, slice(None)))
        [0, 1, 2, 3, 4]
        sage: list(slice_it(iter(range(5)), 5, slice(3)))
        [0, 1, 2]
        sage: list(slice_it(iter(range(5)), 5, slice(-1)))
        [0, 1, 2, 3]
        sage: list(slice_it(iter(range(5)), 5, slice(2, 4)))
        [2, 3]
        sage: list(slice_it(iter(range(5)), 5, slice(3, 1, -1)))
        [3, 2]
    """
    start, stop, step = key.start, key.stop, key.step
    if step is None:
        step = 1
    if step == 0:
        raise ValueError, "attempt to use 0 as step value"
    elif step > 0:
        return islice(it, *slice_indices(key, l))
    else:
        return tuple(islice(it, l))[key]

def peek_it(it):
    r"""
    Returns the first element of an iterator and returns an iterator at
    the same position as the original iterator

    EXAMPLES::

        sage: from sage.combinat.words.utils import peek_it
        sage: it = iter([1, 2, 3])
        sage: it, n = peek_it(it); n
        1
        sage: it.next()
        1
    """
    it, cp = copy_it(it)
    return it, cp.next()

def len_it(it):
    r"""
    Returns the number of elements in it.

    This function will modify the iterator, so if you want to access
    the elements later, make a copy of the iterator and pass it to this
    function.

    EXAMPLES::

        sage: from sage.combinat.words.utils import len_it
        sage: len_it(iter([1, 2, 3]))
        3
    """
    n = 0
    try:
        while True:
            it.next()
            n += 1
    except StopIteration:
        return n

def haslen(obj):
    r"""
    Returns true if obj has a properly defined length

    EXAMPLES::

        sage: from sage.combinat.words.utils import haslen
        sage: haslen([1, 2, 3])
        True
        sage: haslen(33)
        False

    TESTS::

        sage: class test(object):
        ...     def __len__(self):
        ...         return -1
        ...
        sage: haslen(test())
        False
    """
    try:
        len(obj)
    except:
        return False
    return True

def sliceable(obj):
    r"""
    Returns true if obj is completely sliceable, including negative
    step sizes

    EXAMPLES::

        sage: from sage.combinat.words.utils import sliceable
        sage: sliceable([1, 2, 3])
        True
        sage: sliceable(33)
        False

    TESTS::

        sage: class test(object):
        ...     def __getitem__(self, key):
        ...         return islice(iter([1]), *key)
        ...
        sage: sliceable(test())
        False
    """
    try:
        obj[0:0:-1]
    except:
        return False
    return True

def isint(obj):
    r"""
    Returns True if obj is an integer or a custom object representing
    an integer and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.utils import isint
        sage: isint(1)
        True
        sage: isint("2")
        False
        sage: isint(1.0)
        False
    """
    return hasattr(obj, '__index__')

def is_iterable(obj):
    r"""
    Returns true if the obj is iterable.

    EXAMPLES::

        sage: from sage.combinat.words.utils import is_iterable
        sage: is_iterable('123')
        True
        sage: is_iterable([1, 2, 3])
        True
        sage: is_iterable(xrange(1, 4))
        True
        sage: is_iterable(33)
        False
    """
    try:
        iter(obj)
        return True
    except:
        return False

def slice_ok(part):
    r"""
    Returns true if part is a slice and doesn't have funny values.

    EXAMPLES::

        sage: from sage.combinat.words.utils import slice_ok
        sage: slice_ok(slice(None))
        True
        sage: slice_ok(slice(2))
        True
        sage: slice_ok(slice(1, 2, 3))
        True
        sage: slice_ok(slice(None, 2, 3))
        True
        sage: slice_ok(slice(1, None, 3))
        True
        sage: slice_ok(slice(1, 2, None))
        True
        sage: slice_ok(slice("a"))
        False
        sage: slice_ok(slice("a", 1))
        False
        sage: slice_ok(slice(1, 2, "a"))
        False
        sage: slice_ok(1)
        False
    """
    return isinstance(part, slice) and \
           (part.start is None or isint(part.start)) and \
           (part.stop is None or isint(part.stop)) and \
           (part.step is None or isint(part.step))

def slice_indices(s, l):
    r"""
    Implement slice.indices whitout bugs.

    TESTS::

        sage: from sage.combinat.words.utils import slice_indices
        sage: slice_indices(slice(None), 8)
        (0, 8, 1)
        sage: slice_indices(slice(1), 8)
        (0, 1, 1)
        sage: slice_indices(slice(10), 8)
        (0, 8, 1)
        sage: slice_indices(slice(-4), 8)
        (0, 4, 1)
        sage: slice_indices(slice(-10), 8)
        (0, 0, 1)
        sage: slice_indices(slice(1, None), 8)
        (1, 8, 1)
        sage: slice_indices(slice(10, None), 8)
        (8, 8, 1)
        sage: slice_indices(slice(-4, None), 8)
        (4, 8, 1)
        sage: slice_indices(slice(-10, None), 8)
        (0, 8, 1)
        sage: slice_indices(slice(None, None, 2), 8)
        (0, 8, 2)
        sage: slice_indices(slice(None, None, -2), 8)
        (7, -1, -2)
    """
    start, stop, step = s.indices(l)
    if stop*step < start*step:
        stop = start
    return start, stop, step


def reverse_map(d):
    r"""
    Return a new dict with swapped keys and values

    EXAMPLES::

        sage: from sage.combinat.words.utils import reverse_map
        sage: reverse_map({'a': 1, 'b': 2}) == {1: 'a', 2: 'b'}
        True
    """
    return dict(izip(d.itervalues(), d))

def id_f(x):
    r"""
    Dummy identity function for when a function is required but none is
    desired.

    TESTS::

        sage: l = [1, 2, 3]
        sage: l is sage.combinat.words.utils.id_f(l)
        True
    """
    return x

def clamp(x, min, max):
    r"""
    Clamp a value between a maximum and a minimum.

    EXAMPLES::

        sage: from sage.combinat.words.utils import clamp
        sage: clamp(0, -1, 1)
        0
        sage: clamp(-2, 0, 4)
        0
        sage: clamp(10, 0, 4)
        4
    """
    if x < min:
        return min
    elif x > max:
        return max
    else:
        return x

class Factorization(list):
    r"""
    A list subclass having a nicer representation for factorization of
    words.

    TESTS::

        sage: f = sage.combinat.words.utils.Factorization()
        sage: f == loads(dumps(f))
        True
    """
    def __repr__(self):
        r"""
        Returns a string representation of the object.

        TESTS::

            sage: sage.combinat.words.utils.Factorization()
            ()
            sage: sage.combinat.words.utils.Factorization([Word('ab'), Word('ba')])
            (ab.ba)
        """
        from word import FiniteWord_over_OrderedAlphabet
        return '(%s)' % '.'.join(imap(FiniteWord_over_OrderedAlphabet.string_rep, self))
