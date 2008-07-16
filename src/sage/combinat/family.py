from sage.combinat.combinat import CombinatorialClass
from sage.combinat.finite_class import FiniteCombinatorialClass_l

def Family(indices, function = None, name = None, hidden_keys = [], hidden_function = None):
    r"""
    A Family is an associative container which models a family
    (f_i)_{i in I}. Then, f[i] returns the element of the family
    indexed by i. Whenever available, set and combinatorial class
    operations (counting, iteration, listing) on the family are
    induced from those of the index set.

    There are several available implementations (classes) for
    different usages; Family serves as a factory, and will create
    instances of the appropriate classes depending on its arguments.

    EXAMPLES:

        In its simplest form, a list l by itself is considered as the
        family $(l[i]_{i in I})$ where $I$ is the range $0\dots,len(l)$.
        So Family(l) just returns it.

            sage: f = Family([1,2,3])
            sage: f
            [1, 2, 3]

        A family can also be constructed from a dictionary t. The
        resulting family is very close to t, except that the elements
        of the family are the values of t. Here, we define the family
        (f_i)_{i in \{3,4,7\}} with f_3='a', f_4='b', and f_7='d':

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
            sage: 'b' in f       # todo: not implemented
            True
            sage: 'e' in f       # todo: not implemented
            False

        A familly can also be constructed by its index set $I$ and a
        function $f$, as in $(f(i))_{i in I}$:

            sage: f = Family([3,4,7], lambda i: 2*i)
            sage: f
            Finite family {3: 6, 4: 8, 7: 14}
            sage: f.keys()
            [3, 4, 7]
            sage: f[7]
            14
            sage: list(f)
            [6, 8, 14]
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        By default, if the index set is a list, all images are
        computed right away, and stored in an internal
        dictionary. Note that this requires all the elements of the
        list to be hashable. One can ask instead for the images $f(i)$
        to be computed lazily, when needed:

            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f
            Lazy family (f(i))_{i in [3, 4, 7]}
            sage: f[7]
            14
            sage: list(f)
            [6, 8, 14]
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        This allows in particular for modeling infinite families:
            sage: f = Family(ZZ, lambda i: 2*i)
            sage: f
            Lazy family (f(i))_{i in Integer Ring}
            sage: f.keys()
            Integer Ring
            sage: f[1]
            2
            sage: f[-5]
            -10
            sage: i = f.__iter__()
            sage: i.next(), i.next(), i.next(), i.next(), i.next()
            (0, 2, -2, 4, -4)

        Caveat: families with lazy behavior cannot be pickled


        Finally, it can occasionally be useful to add some hidden
        elements in a familly, which are accessible as f[i], but
        do not appear in the keys or the container operations.

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
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        The following example illustrates when the function is actually called:
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
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        Here is a close variant where the function for the hidden keys
        is different from that for the other keys:

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
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        Family behaves the same way with FiniteCombinatorialClass
        instances and lists.  This feature will eventually disapear
        when FiniteCombinatorialClass won't be needed anymore.

            sage: f = Family(FiniteCombinatorialClass([1,2,3]))
            sage: f
            Combinatorial class with elements in [1, 2, 3]

            sage: f = Family(FiniteCombinatorialClass([3,4,7]), lambda i: 2*i)
            sage: f
            Finite family {3: 6, 4: 8, 7: 14}
            sage: f.keys()
            [3, 4, 7]
            sage: f[7]
            14
            sage: list(f)
            [6, 8, 14]
            sage: [ x for x in f]
            [6, 8, 14]
            sage: len(f)
            3

        TESTS:
            sage: f = Family({1:'a', 2:'b', 3:'c'})
            sage: f
            Finite family {1: 'a', 2: 'b', 3: 'c'}
            sage: f[2]
            'b'
            sage: loads(dumps(f)) == f
            True

            sage: f = Family(range(1,27), lambda i: chr(i+96))
            sage: f
                Finite family {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h', 9: 'i', 10: 'j', 11: 'k', 12: 'l', 13: 'm', 14: 'n', 15: 'o', 16: 'p', 17: 'q', 18: 'r', 19: 's', 20: 't', 21: 'u', 22: 'v', 23: 'w', 24: 'x', 25: 'y', 26: 'z'}
            sage: f[2]
            'b'
    """
    assert(type(hidden_keys) == list)
    if function == None and hidden_keys == []:
        if type(indices) == list or isinstance(indices, FiniteCombinatorialClass_l) or isinstance(indices, FiniteFamily) or isinstance(indices, LazyFamily):
            return indices
        if type(indices) == dict:
            return FiniteFamily(indices)
    else:
        if type(indices) == list or isinstance(indices, FiniteCombinatorialClass_l):
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

    def hidden_keys(self):
        """
        Returns the hidden keys of the family, if any
        """
        return []

    def zip(self, f, other, name = None):
        """
        Given two families with same index set $I$ (and same hidden keys
        if relevant), returns the family $( f(self[i], other[i]) )_{i in I}$

        TODO: generalize to any number of families and merge with map?

        EXAMPLES:
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
        Returns the family $( f(self[i]) )_{i in I}$, where
        $I$ is the index set of self.

        TODO: good name?

        EXAMPLES:
            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: g = f.map(lambda x: x+'1')
            sage: list(g)
            ['a1', 'b1', 'd1']
        """
        return Family(self.keys(), lambda i: f(self[i]), hidden_keys = self.hidden_keys(), name = name)

class FiniteFamily(AbstractFamily):
    r"""
    A FiniteFamily is an associative container which models a finite
    family (f_i)_{i in I}. Its elements $f_i$ are therefore its
    values. Instances should be created via the Family factory, which
    see for further examples and tests.

    EXAMPLES:
    We define the family (f_i)_{i in \{3,4,7\}} with f_3=a, f_4=b, and f_7=d
        sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})

    Individual elements are accessible as in a usual dictionary:
        sage: f[7]
        'd'

    And the other usual dictionary operations are also available:
        sage: len(f)
        3
        sage: f.keys()
        [3, 4, 7]

    However f behaves as a container for the $f_i$'s:
        sage: list(f)
        ['a', 'b', 'd']
        sage: [ x for x in f ]
        ['a', 'b', 'd']
        sage: f == loads(dumps(f))
        True
    """

    def __init__(self, dictionary, keys = None):
        # TODO: use keys to specify the order of the elements
        self.dictionary = dictionary
        self.keys = dictionary.keys
        self.__iter__ = dictionary.itervalues
        self.list = dictionary.values # should not be required
        self.values = dictionary.values
        self.__getitem__ = dictionary.__getitem__

    def __repr__(self):
        return "Finite family %s"%self.dictionary

    def count(self):
        return len(self.dictionary)

    # Why isn't it sufficient to set self.__getitem__ as above?
    def __getitem__(self, i):
        return self.__getitem__(i)

    # For the pickle and copy modules
    def __getstate__(self):
        return {'dictionary': self.dictionary}

    def __setstate__(self, state):
        self.__init__(state['dictionary'])

class FiniteFamilyWithHiddenKeys(FiniteFamily):
    r"""
    A close variant of FiniteFamily where the family contains some
    hidden keys whose corresponding values are computed lazily (and
    remembered). Instances should be created via the Family factory,
    which see for examples and tests.

    Caveat: instances of this class cannot be pickled, because the
    hidden_function itself cannot.
    """
    def __init__(self, dictionary, hidden_keys, hidden_function):
        FiniteFamily.__init__(self, dictionary)
        self._hidden_keys = hidden_keys
        self.hidden_function = hidden_function
        self.hidden_dictionary = {}

        # would be better to define as usual method
        # any better to unset the def of __getitem__ by FiniteFamily?
        #self.__getitem__ = lambda i: dictionary[i] if dictionary.has_key(i) else hidden_dictionary[i]

    def __getitem__(self, i):
        if self.dictionary.has_key(i):
            return self.dictionary[i]
        if not self.hidden_dictionary.has_key(i):
            if not i in self._hidden_keys:
                raise KeyError
            self.hidden_dictionary[i] = self.hidden_function(i)
        return self.hidden_dictionary[i]

    def hidden_keys(self):
        return self._hidden_keys

    def __getstate__(self):
        raise NotImplementedError

class LazyFamily(AbstractFamily):
    r"""
    A LazyFamily(I, f) is an associative container which models the
    (possibly infinite) family (f(i))_{i in I}.

    Instances should be created via the Family factory, which see for
    examples and tests.
    """

    def __init__(self, set, function, name = "f"):
        self.set = set
        self.name = name
        self.__getitem__ = function

    def __repr__(self):
        return "Lazy family (%s(i))_{i in %s}"%(self.name,self.set)

    def keys(self):
        return self.set

    def __iter__(self):
        for i in self.set.__iter__():
            yield self[i]

    # Should disappear
    iterator = __iter__

    # Why isn't it sufficient to set self.__getitem__ as above?
    def __getitem__(self, i):
        return self.__getitem__(i)
