"""
Families
"""
#*****************************************************************************
#       Copyright (C) 2008 Nicolas Thiery <nthiery at users.sf.net>,
#                          Mike Hansen <mhansen@gmail.com>,
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
from sage.combinat.combinat import CombinatorialClass
from sage.combinat.finite_class import FiniteCombinatorialClass
from sage.rings.integer import Integer

def Family(indices, function = None, name = None, hidden_keys = [], hidden_function = None):
    r"""
    A Family is an associative container which models a family
    `(f_i)_{i \in I}`. Then, f[i] returns the element of the family
    indexed by i. Whenever available, set and combinatorial class
    operations (counting, iteration, listing) on the family are induced
    from those of the index set.

    There are several available implementations (classes) for different
    usages; Family serves as a factory, and will create instances of
    the appropriate classes depending on its arguments.

    EXAMPLES:

    In its simplest form, a list l or a tuple by itself is considered as the
    family `(l[i]_{i \in I})` where `I` is the range `0\dots,len(l)`. So
    Family(l) returns the corresponding family.

    ::

        sage: f = Family([1,2,3])
        sage: f
        Family (1, 2, 3)
        sage: f = Family((1,2,3))
        sage: f
        Family (1, 2, 3)

    A family can also be constructed from a dictionary t. The resulting
    family is very close to t, except that the elements of the family
    are the values of t. Here, we define the family `(f_i)_{i \in \{3,4,7\}}`
    with f_3='a', f_4='b', and f_7='d'::

        sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
        sage: f
        Finite family {3: 'a', 4: 'b', 7: 'd'}
        sage: f[7]
        'd'
        sage: len(f)
        3
        sage: list(f)
        ['a', 'b', 'd']
        sage: [ x for x in f ]
        ['a', 'b', 'd']
        sage: f.keys()
        [3, 4, 7]
        sage: 'b' in f
        True
        sage: 'e' in f
        False

    A family can also be constructed by its index set `I` and
    a function `f`, as in `(f(i))_{i \in I}`::

        sage: f = Family([3,4,7], lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    By default, if the index set is a list or a tuple, all images are
    computed right away, and stored in an internal dictionary::

        sage: f = Family((3,4,7), lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}

    Note that this requires all the elements of the list to be
    hashable. One can ask instead for the images `f(i)` to be computed
    lazily, when needed::

        sage: f = LazyFamily([3,4,7], lambda i: 2r*i)
        sage: f
        Lazy family (f(i))_{i in [3, 4, 7]}
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]

    This allows in particular for modeling infinite families::

        sage: f = Family(ZZ, lambda i: 2r*i)
        sage: f
        Lazy family (f(i))_{i in Integer Ring}
        sage: f.keys()
        Integer Ring
        sage: f[1]
        2
        sage: f[-5]
        -10
        sage: i = iter(f)
        sage: i.next(), i.next(), i.next(), i.next(), i.next()
        (0, 2, -2, 4, -4)

    Beware that for those kind of families len(f) is not supposed to
    work. As a replacement, use the .cardinality() method::

       sage: f = LazyFamily(Permutations(3), attrcall("to_lehmer_code"))
       sage: list(f)
       [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [2, 0, 0], [2, 1, 0]]
       sage: f.cardinality()
       6

    Caveat: Only certain families with lazy behavior can be pickled. In
    particular, only functions that work with Sage's pickle_function
    and unpickle_function (in sage.misc.fpickle) will correctly
    unpickle. The following two work::

       sage: loads(dumps(LazyFamily(Permutations(3), lambda p: p.to_lehmer_code())))
       Lazy family (f(i))_{i in Standard permutations of 3}

       sage: loads(dumps(LazyFamily(Permutations(3), attrcall("to_lehmer_code"))))
       Lazy family (f(i))_{i in Standard permutations of 3}

    But this one dont::

       sage: def plus_n(n): return lambda x: x+n
       sage: loads(dumps(LazyFamily([1,2,3], plus_n(3))))
       Traceback (most recent call last):
       ...
       ValueError: Cannot pickle code objects from closures

    Finally, it can occasionally be useful to add some hidden elements
    in a family, which are accessible as f[i], but do not appear in the
    keys or the container operations.

    ::

        sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        4
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    The following example illustrates when the function is actually
    called::

        sage: def compute_value(i):
        ...       print('computing 2*'+str(i))
        ...       return 2*i
        sage: f = Family([3,4,7], compute_value, hidden_keys=[2])
        computing 2*3
        computing 2*4
        computing 2*7
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        computing 2*2
        4
        sage: f[2]
        4
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    Here is a close variant where the function for the hidden keys is
    different from that for the other keys::

        sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2], hidden_function = lambda i: 3*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        6
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    Family behaves the same way with FiniteCombinatorialClass instances
    and lists. This feature will eventually disappear when
    FiniteCombinatorialClass won't be needed anymore.

    ::

        sage: f = Family(FiniteCombinatorialClass([1,2,3]))
        sage: f
        Combinatorial class with elements in [1, 2, 3]

    ::

        sage: f = Family(FiniteCombinatorialClass([3,4,7]), lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    TESTS::

        sage: f = Family({1:'a', 2:'b', 3:'c'})
        sage: f
        Finite family {1: 'a', 2: 'b', 3: 'c'}
        sage: f[2]
        'b'
        sage: loads(dumps(f)) == f
        True

    ::

        sage: f = Family(range(1,27), lambda i: chr(i+96))
        sage: f
            Finite family {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h', 9: 'i', 10: 'j', 11: 'k', 12: 'l', 13: 'm', 14: 'n', 15: 'o', 16: 'p', 17: 'q', 18: 'r', 19: 's', 20: 't', 21: 'u', 22: 'v', 23: 'w', 24: 'x', 25: 'y', 26: 'z'}
        sage: f[2]
        'b'

    The factory ``Family`` is supposed to be idempotent. We test this feature here::

        sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = Family([3,4,7], lambda i: 2r*i, hidden_keys=[2])
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = LazyFamily([3,4,7], lambda i: 2r*i)
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = TrivialFamily([3,4,7])
        sage: g = Family(f)
        sage: f == g
        True
    """
    assert(type(hidden_keys) == list)
    if function is None and hidden_keys == []:
        if isinstance(indices, dict):
            return FiniteFamily(indices)
        if isinstance(indices, (list, tuple) ):
            return TrivialFamily(indices)
        if isinstance(indices, (FiniteCombinatorialClass,
                                FiniteFamily, LazyFamily, TrivialFamily) ):
            return indices
    else:
        if isinstance(indices, (list, tuple, FiniteCombinatorialClass) ):
            if not hidden_keys == []:
                if hidden_function is None:
                    hidden_function = function
                return FiniteFamilyWithHiddenKeys(dict([(i, function(i)) for i in indices]),
                                                  hidden_keys, hidden_function)
            else:
                return FiniteFamily(dict([(i, function(i)) for i in indices]), keys = indices)
        elif hidden_keys == [] and hidden_function is None:
            return LazyFamily(indices, function)
    raise NotImplementedError

class AbstractFamily(CombinatorialClass):
    """
    The abstract class for family

    Any family belongs to a class which inherits from ``AbstractFamily``.
    """
    def hidden_keys(self):
        """
        Returns the hidden keys of the family, if any.

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: f.hidden_keys()
            []
        """
        return []

    def zip(self, f, other, name = None):
        """
        Given two families with same index set `I` (and same hidden
        keys if relevant), returns the family
        `( f(self[i], other[i]) )_{i \in I}`

        TODO: generalize to any number of families and merge with map?

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: g = Family({3: '1', 4: '2', 7: '3'})
            sage: h = f.zip(lambda x,y: x+y, g)
            sage: list(h)
            ['a1', 'b2', 'd3']
        """
        assert(self.keys() == other.keys())
        assert(self.hidden_keys() == other.hidden_keys())
        return Family(self.keys(), lambda i: f(self[i],other[i]), hidden_keys = self.hidden_keys(), name = name)

    def map(self, f, name = None):
        """
        Returns the family `( f(\mathtt{self}[i]) )_{i \in I}`, where
        `I` is the index set of self.

        TODO: good name?

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: g = f.map(lambda x: x+'1')
            sage: list(g)
            ['a1', 'b1', 'd1']
        """
        return Family(self.keys(), lambda i: f(self[i]), hidden_keys = self.hidden_keys(), name = name)

class FiniteFamily(AbstractFamily):
    r"""
    A FiniteFamily is an associative container which models a finite
    family `(f_i)_{i \in I}`. Its elements `f_i` are therefore
    its values. Instances should be created via the Family factory,
    which see for further examples and tests.

    EXAMPLES: We define the family `(f_i)_{i \in \{3,4,7\}}` with f_3=a,
    f_4=b, and f_7=d

    ::

        sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})

    Individual elements are accessible as in a usual dictionary::

        sage: f[7]
        'd'

    And the other usual dictionary operations are also available::

        sage: len(f)
        3
        sage: f.keys()
        [3, 4, 7]

    However f behaves as a container for the `f_i`'s::

        sage: list(f)
        ['a', 'b', 'd']
        sage: [ x for x in f ]
        ['a', 'b', 'd']
    """

    def __init__(self, dictionary, keys = None):
        """
        TESTS::

            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: f == loads(dumps(f))
            True

        Check for bug #5538::

            sage: d = {1:"a", 3:"b", 4:"c"}
            sage: f = Family(d)
            sage: d[2] = 'DD'
            sage: f
            Finite family {1: 'a', 3: 'b', 4: 'c'}
            """
        # TODO: use keys to specify the order of the elements
        self.dictionary = dict(dictionary)
        self.keys = dictionary.keys
        self.values = dictionary.values

    def __repr__(self):
        """
        EXAMPLES::

            sage: FiniteFamily({3: 'a'})
            Finite family {3: 'a'}
        """
        return "Finite family %s"%self.dictionary

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: f = FiniteFamily({3: 'a'})
            sage: 'a' in f
            True
            sage: 'b' in f
            False
        """
        return x in self.values()

    def __len__(self):
        """
        Returns the number of elements in self.

        EXAMPLES::

            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: len(f)
            3
        """
        return len(self.dictionary)

    def cardinality(self):
        """
        Returns the number of elements in self.

        EXAMPLES::

            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: f.cardinality()
            3
        """
        return Integer(len(self.dictionary))

    def __iter__(self):
        """
        EXAMPLES::

            sage: f = FiniteFamily({3: 'a'})
            sage: i = iter(f)
            sage: i.next()
            'a'
        """
        return iter(self.values())

    def __getitem__(self, i):
        """
        Note that we can't just do self.__getitem__ =
        dictionary.__getitem__ in the __init__ method since Python
        queries the object's type/class for the special methods rather than
        querying the object itself.

        EXAMPLES::

            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: f[3]
            'a'
        """
        return self.dictionary.__getitem__(i)

    # For the pickle and copy modules
    def __getstate__(self):
        """
        TESTS::

            sage: f = FiniteFamily({3: 'a'})
            sage: f.__getstate__()
            {'dictionary': {3: 'a'}}
        """
        return {'dictionary': self.dictionary}

    def __setstate__(self, state):
        """
        TESTS::

            sage: f = FiniteFamily({3: 'a'})
            sage: f.__setstate__({'dictionary': {4:'b'}})
            sage: f
            Finite family {4: 'b'}
        """
        self.__init__(state['dictionary'])

class FiniteFamilyWithHiddenKeys(FiniteFamily):
    r"""
    A close variant of FiniteFamily where the family contains some
    hidden keys whose corresponding values are computed lazily (and
    remembered). Instances should be created via the Family factory,
    which see for examples and tests.

    Caveat: Only instances of this class whose functions are compatible
    with sage.misc.fpickle can be pickled.
    """
    def __init__(self, dictionary, hidden_keys, hidden_function):
        """
        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2r*i, hidden_keys=[2])
            sage: f == loads(dumps(f))
            True
        """
        FiniteFamily.__init__(self, dictionary)
        self._hidden_keys = hidden_keys
        self.hidden_function = hidden_function
        self.hidden_dictionary = {}

        # would be better to define as usual method
        # any better to unset the def of __getitem__ by FiniteFamily?
        #self.__getitem__ = lambda i: dictionary[i] if dictionary.has_key(i) else hidden_dictionary[i]

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: f[3]
            6
            sage: f[2]
            4
            sage: f[5]
            Traceback (most recent call last):
            ...
            KeyError
        """
        if i in self.dictionary:
            return self.dictionary[i]

        if i not in self.hidden_dictionary:
            if i not in self._hidden_keys:
                raise KeyError
            self.hidden_dictionary[i] = self.hidden_function(i)

        return self.hidden_dictionary[i]

    def hidden_keys(self):
        """
        Returns self's hidden keys.

        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: f.hidden_keys()
            [2]
        """
        return self._hidden_keys

    def __getstate__(self):
        """
        TESTS::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: d = f.__getstate__()
            sage: d['hidden_keys']
            [2]
        """
        from sage.misc.fpickle import pickle_function
        f = pickle_function(self.hidden_function)
        return {'dictionary': self.dictionary,
                'hidden_keys': self._hidden_keys,
                'hidden_dictionary': self.hidden_dictionary,
                'hidden_function': f}

    def __setstate__(self, d):
        """
        TESTS::

            sage: f = Family([3,4,7], lambda i: 2r*i, hidden_keys=[2])
            sage: d = f.__getstate__()
            sage: f = Family([4,5,6], lambda i: 2r*i, hidden_keys=[2])
            sage: f.__setstate__(d)
            sage: f.keys()
            [3, 4, 7]
            sage: f[3]
            6

            sage: f = LazyFamily(Permutations(3), lambda p: p.to_lehmer_code())
            sage: f == loads(dumps(f))
            True

            sage: f = LazyFamily(Permutations(3), attrcall("to_lehmer_code"))
            sage: f == loads(dumps(f))
            True
        """
        hidden_function = d['hidden_function']
        if isinstance(hidden_function, str):
        # Let's assume that hidden_function is an unpicled function.
            from sage.misc.fpickle import unpickle_function
            hidden_function = unpickle_function(hidden_function)
        self.__init__(d['dictionary'], d['hidden_keys'], hidden_function)
        self.hidden_dictionary = d['hidden_dictionary']


class LazyFamily(AbstractFamily):
    r"""
    A LazyFamily(I, f) is an associative container which models the
    (possibly infinite) family `(f(i))_{i \in I}`.

    Instances should be created via the Family factory, which see for
    examples and tests.
    """
    def __init__(self, set, function, name = "f"):
        """
        TESTS::

            sage: f = LazyFamily([3,4,7], lambda i: 2r*i); f
            Lazy family (f(i))_{i in [3, 4, 7]}
            sage: f == loads(dumps(f))
            True

        Check for bug #5538::

            sage: l = [3,4,7]
            sage: f = LazyFamily(l, lambda i: 2r*i);
            sage: l[1] = 18
            sage: f
            Lazy family (f(i))_{i in [3, 4, 7]}

        """
        from copy import copy
        self.set = copy(set)
        self.name = name
        self.function = function

    def __repr__(self):
        """
        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2*i); f
            Lazy family (f(i))_{i in [3, 4, 7]}
        """
        return "Lazy family (%s(i))_{i in %s}"%(self.name,self.set)

    def keys(self):
        """
        Returns self's keys.

        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f.keys()
            [3, 4, 7]
        """
        return self.set

    def cardinality(self):
        """
        Return the number of elements in self.

        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f.cardinality()
            3
        """
        try:
            return Integer(len(self.set))
        except AttributeError:
            return self.set.cardinality()

    def __iter__(self):
        """
        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: [i for i in f]
            [6, 8, 14]
        """
        for i in self.set:
            yield self[i]

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f[3]
            6

        TESTS::

            sage: f[5]
            10
        """
        return self.function(i)

    def __getstate__(self):
        """
        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2r*i)
            sage: d = f.__getstate__()
            sage: d['set']
            [3, 4, 7]
        """
        f = self.function
        # This should be done once for all by registering
        # sage.misc.fpickle.pickle_function to copy_reg
        if type(f) is type(Family): # TODO: where is the python `function` type?
            from sage.misc.fpickle import pickle_function
            f = pickle_function(f)

        return {'set': self.set,
                'name': self.name,
                'function': f}

    def __setstate__(self, d):
        """
        EXAMPLES::

            sage: f = LazyFamily([3,4,7], lambda i: 2r*i)
            sage: d = f.__getstate__()
            sage: f = LazyFamily([4,5,6], lambda i: 2r*i)
            sage: f.__setstate__(d)
            sage: f.keys()
            [3, 4, 7]
            sage: f[3]
            6
        """
        function = d['function']
        if isinstance(function, str):
        # Let's assume that function is an unpicled function.
            from sage.misc.fpickle import unpickle_function
            function = unpickle_function(function)

        self.__init__(d['set'], function, d['name'])


class TrivialFamily(AbstractFamily):
    r"""
    ``TrivialFamily(c)`` turn the container c into a family indexed by
    the set `{0,\dots, len(c)}`. The container `c` can be either a list or a
    tuple.

    Instances should be created via the Family factory, which see for
    examples and tests.
    """
    def __init__(self, set):
        """
        EXAMPLES::

            sage: f = TrivialFamily((3,4,7)); f
            Family (3, 4, 7)
            sage: f = TrivialFamily([3,4,7]); f
            Family (3, 4, 7)
            sage: f == loads(dumps(f))
            True
        """
        self.set = tuple(set)

    def __repr__(self):
        """
        EXAMPLES::

            sage: f = TrivialFamily([3,4,7]); f
            Family (3, 4, 7)
        """
        return "Family %s"%((self.set),)

    def keys(self):
        """
        Returns self's keys.

        EXAMPLES::

            sage: f = TrivialFamily([3,4,7])
            sage: f.keys()
            [0, 1, 2]
        """
        return range(len(self.set))

    def cardinality(self):
        """
        Return the number of elements in self.

        EXAMPLES::

            sage: f = TrivialFamily([3,4,7])
            sage: f.cardinality()
            3
        """
        return Integer(len(self.set))

    def __iter__(self):
        """
        EXAMPLES::

            sage: f = TrivialFamily([3,4,7])
            sage: [i for i in f]
            [3, 4, 7]
        """
        for i in self.set:
            yield i

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: f = TrivialFamily([3,4,7])
            sage: f[1]
            4
        """
        return self.set[i]

    def __getstate__(self):
        """
        TESTS::

            sage: f = TrivialFamily([3,4,7])
            sage: f.__getstate__()
            {'set': (3, 4, 7)}
        """
        return {'set': self.set}

    def __setstate__(self, state):
        """
        TESTS::

            sage: f = TrivialFamily([3,4,7])
            sage: f.__setstate__({'set': (2, 4, 8)})
            sage: f
            Family (2, 4, 8)
        """
        self.__init__(state['set'])
