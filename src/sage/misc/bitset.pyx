"""
Bitsets

A Python interface to the fast bitsets in Sage.  Bitsets are fast
binary sets that store elements by toggling bits in an array of
numbers.  A bitset can store values between `0` and ``capacity-1``
(where ``capacity`` is finite, but arbitrary).  The storage cost is
linear in ``capacity``.

.. warning::

    This class is most likely to be useful as a way to store Cython
    bitsets in Python datastructures, acting on them using the Cython
    inline functions.  If you want to use these classes for a Python
    set type, the Python ``set`` or ``frozenset`` data types may be
    faster.
"""
#*****************************************************************************
#       Copyright (C) 2009 Jason Grout <jason-sage@creativetrax.com>
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

include 'bitset.pxi'

cdef class FrozenBitset:
    """
    A frozen bitset class which leverages inline Cython functions for creating
    and manipulating bitsets.

    A bitset can be thought of in two ways.  First, as a set of
    elements from the universe of the `n` natural numbers `0`, `1`,
    ..., `n-1` (where `n`, the capacity, can be specified), with
    typical set operations (intersection, union, symmetric difference,
    etc.).  Secondly, a bitset can be thought of as a binary vector
    with typical binary operations (i.e., ``and``, ``or``, ``xor``,
    etc.).  This class supports both interfaces.

    The interface in this class mirrors the interface in the
    ``frozenset`` datatype of Python.

    .. warning::

        This class is most likely to be useful as a way to store
        Cython bitsets in Python datastructures, acting on them using
        the Cython inline functions.  If you want to use this class
        for a Python set type, the Python ``frozenset`` data type may
        be faster.

    EXAMPLES::

        sage: a=FrozenBitset('1101')
        sage: loads(dumps(a))==a
        True

        sage: a=FrozenBitset('1101'*64)
        sage: loads(dumps(a))==a
        True
    """

    def __cinit__(self, iter=None, capacity=None):
        """
        Allocate the bitset.

        EXAMPLES::

            sage: FrozenBitset('1101')
            1101
            sage: FrozenBitset('1101'*32)
            11011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101110111011101
        """
        if capacity is None:
            bitset_init(self._bitset, 1)
        else:
            bitset_init(self._bitset, capacity)

    def __dealloc__(self):
        """
        Deallocate the C bitset data structure.

        EXAMPLES::

            sage: a=FrozenBitset('11010')
            sage: del a
            sage: b=FrozenBitset('11010'*64)
            sage: del b
        """
        bitset_free(self._bitset)

    def __init__(self, iter=None, capacity=None):
        """
        Initialize a bitset.

        INPUTS:

        - ``iter`` - initialization parameter.  If this is a Bitset, then
        it is copied.  If it is None, then the bitset is set to the
        empty set.  If it is a string, then the bitset is initialized
        by including an element if the index of the string is '1'.  If
        it is an iterable, then it is assumed to contain a list of
        integers, and those integers are placed in the set.

        - ``capacity`` - The maximum capacity of the bitset.  If this is not
        specified, then it is automatically calculated from the passed
        iterable.  It must be at least one.

        EXAMPLES::

            sage: FrozenBitset(capacity=3)
            000
            sage: FrozenBitset('11011')
            11011
            sage: FrozenBitset('110'*32)
            110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110
            sage: FrozenBitset([0,3,2])
            1011
            sage: FrozenBitset(set([0,3,2]))
            1011
            sage: FrozenBitset(FrozenBitset('11011'))
            11011
            sage: FrozenBitset([0,3,2], capacity=10)
            1011000000
            sage: FrozenBitset([i for i in range(100) if i%2==0])
            101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101

        TESTS:

        The capacity must be at least one::

            sage: FrozenBitset()
            0

        If capacity is specified, it must match up with the
        initialization items::

            sage: FrozenBitset('10', capacity=3)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed string
            sage: FrozenBitset([0,3,2], capacity=2)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not allow storing the passed iterable
            sage: FrozenBitset(FrozenBitset('1101'), capacity=2)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed bitset
            sage: FrozenBitset(FrozenBitset('1101'), capacity=5)
            Traceback (most recent call last):
            ...
            ValueError: bitset capacity does not match passed bitset
        """
        cdef unsigned long n
        cdef FrozenBitset b
        if iter is not None and not isinstance(iter, (str,tuple,list,FrozenBitset,Bitset)):
            iter = list(iter)

        if iter is None:
            bitset_clear(self._bitset)
        elif isinstance(iter, (FrozenBitset, Bitset)):
            b = iter
            if capacity is None:
                bitset_realloc(self._bitset, b._bitset.size)
            elif self._bitset.size != b._bitset.size:
                raise ValueError, "bitset capacity does not match passed bitset"
            bitset_copy(self._bitset, b._bitset)
        elif isinstance(iter, str):
            if capacity is None:
                bitset_realloc(self._bitset, len(iter))
            elif self._bitset.size != len(iter):
                raise ValueError, "bitset capacity does not match passed string"
            bitset_from_str(self._bitset, iter)
        else: # an iterable
            iter = list(iter)
            if len(iter)>0:
                need_capacity = max(iter)+1
            else:
                need_capacity=1
            if capacity is None:
                bitset_realloc(self._bitset, need_capacity)
            elif self._bitset.size < need_capacity:
                raise ValueError, "bitset capacity does not allow storing the passed iterable"
            bitset_clear(self._bitset)
            for n in iter:
                bitset_add(self._bitset, n)

    cdef FrozenBitset _new(self,long int capacity):
        """
        Return an object of the same type as self, initialized with a bitset of capacity ``capacity``.

        """
        cdef FrozenBitset b
        b = FrozenBitset.__new__(FrozenBitset,None, capacity)
        return b

    def __getstate__(self):
        """
        Return the current state of the object as a string.

        EXAMPLES::

            sage: FrozenBitset('1101').__getstate__()
            '1101'
            sage: FrozenBitset('110'*32).__getstate__()
            '110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110'

        """
        return str(self)

    def __setstate__(self,state):
        """
        Set the state of the object from the string state.

        EXAMPLES::

            sage: a=FrozenBitset()
            sage: a.__setstate__('1101')
            sage: a
            1101
            sage: a.__setstate__('110'*32)
            sage: a
            110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110
        """
        bitset_realloc(self._bitset, len(state))
        bitset_from_str(self._bitset, state)

    def __iter__(self):
        """
        Return an iterator over self.

        EXAMPLES::

            sage: list(FrozenBitset('11011'))
            [0, 1, 3, 4]
            sage: list(FrozenBitset('00001'*20))
            [4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, 89, 94, 99]
            sage: set(FrozenBitset('11011'))
            set([0, 1, 3, 4])
        """
        return iter(bitset_list(self._bitset))




    cpdef FrozenBitset _larger_capacity_(self, long capacity):
        """
        Return a copy of self where the bitset has the maximum of the
        current capacity and the capacity passed.  If no resizing is needed,
        return self.

        This function is mainly used to satisfy the assumption of the
        underlying bitset functions that all bitsets are of the same
        capacity.

        INPUTS:

        - ``capacity`` - the underlying bitset of the returned bitset
          will have this capacity if it is bigger than the current
          capacity.


        EXAMPLES::

            sage: a=FrozenBitset('11010')
            sage: a.capacity()
            5
            sage: a._larger_capacity_(4) is a
            True
            sage: a._larger_capacity_(5) is a
            True
            sage: b=a._larger_capacity_(6)
            sage: b
            110100
            sage: b is a
            False
            sage: b.capacity()
            6
            sage: c=a._larger_capacity_(98)
            sage: c
            11010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: c.capacity()
            98
        """
        cdef FrozenBitset temp
        if self._bitset.size >= capacity:
            return self
        else:
            temp = self._new(self._bitset.size)
            bitset_copy(temp._bitset, self._bitset)
            bitset_realloc(temp._bitset, capacity)
            return temp

    cpdef long capacity(self):
        """
        Return the size of the underlying bitset.  The maximum value
        that can be stored in the current underlying bitset is
        ``self.capacity()-1``.

        EXAMPLES::

            sage: FrozenBitset('11000').capacity()
            5
            sage: FrozenBitset('110'*32).capacity()
            96
            sage: FrozenBitset(range(20),capacity=450).capacity()
            450
        """
        return self._bitset.size

    def __hash__(self):
        """
        Return a hash value for a bitset

        EXAMPLES::

            sage: hash(FrozenBitset(capacity=5))
            0
            sage: hash(FrozenBitset('10110'))
            13
            sage: hash(FrozenBitset('10110'+'0'*120,capacity=125))
            13
        """
        cdef long hash
        hash = bitset_hash(self._bitset)
        if hash == -1:
            return 0
        else:
            return hash

    cpdef bint isempty(self):
        """
        Return True if the bitset is empty; False otherwise.

        EXAMPLES::

            sage: FrozenBitset().isempty()
            True
            sage: FrozenBitset([1]).isempty()
            False
            sage: FrozenBitset([],capacity=110).isempty()
            True
            sage: FrozenBitset(range(99)).isempty()
            False
        """
        return bitset_isempty(self._bitset)

    def __richcmp__(FrozenBitset self, FrozenBitset other not None, int op):
        """
        Implement comparisons, using the Cython richcmp convention.
        Comparison is done by inclusion (a set is less than another if
        it is a subset).

        EXAMPLES::

            sage: FrozenBitset('11') < FrozenBitset('101')
            False
            sage: FrozenBitset('11') <= FrozenBitset('110')
            True
            sage: FrozenBitset('11') != FrozenBitset('10')
            True
            sage: FrozenBitset('11') == FrozenBitset('10')
            False
            sage: FrozenBitset('11') == FrozenBitset('110')
            True
            sage: FrozenBitset('11') > FrozenBitset('10')
            True
            sage: FrozenBitset('11') >= FrozenBitset('10')
            True
            sage: FrozenBitset('11') < FrozenBitset('110'*32)
            True
        """
        cdef FrozenBitset left, right

        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        if op == 2: # ==
            return bitset_eq(left._bitset, right._bitset)
        elif op == 3: # !=
            return not bitset_eq(left._bitset, right._bitset)
        elif op == 0: # <
            return bitset_issubset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == 1: # <=
            return bitset_issubset(left._bitset, right._bitset)
        elif op == 4: # >
            return bitset_issuperset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == 5: # >=
            return bitset_issuperset(left._bitset, right._bitset)

    cpdef bint issubset(self, FrozenBitset other):
        """
        Test to see if the self is a subset of other.

        EXAMPLES::

            sage: FrozenBitset('11').issubset(FrozenBitset('01'))
            False
            sage: FrozenBitset('01').issubset(FrozenBitset('11'))
            True
            sage: FrozenBitset('01').issubset(FrozenBitset('01'*45))
            True
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset left, right
        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        return bitset_issubset(left._bitset, right._bitset)

    cpdef bint issuperset(self, FrozenBitset other):
        """
        Test to see if the self is a superset of other.

        EXAMPLES::

            sage: FrozenBitset('11').issuperset(FrozenBitset('01'))
            True
            sage: FrozenBitset('01').issuperset(FrozenBitset('11'))
            False
            sage: FrozenBitset('01').issuperset(FrozenBitset('10'*45))
            False
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset left, right
        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        return bitset_issuperset(left._bitset, right._bitset)

    cpdef bint isdisjoint(self, FrozenBitset other):
        """
        Test to see if the self is disjoint from other.

        EXAMPLES::

            sage: FrozenBitset('11').isdisjoint(FrozenBitset('01'))
            False
            sage: FrozenBitset('01').isdisjoint(FrozenBitset('001'))
            True
            sage: FrozenBitset('00101').isdisjoint(FrozenBitset('110'*35))
            False
        """
        cdef bint retval
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset smaller, larger
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_init(temp, self._bitset.size)
            bitset_intersection(temp, self._bitset, other._bitset)
            retval = bitset_isempty(temp)
            bitset_free(temp)
            return retval
        elif self._bitset.size < other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        bitset_init(temp, smaller._bitset.size)
        bitset_copy(temp, smaller._bitset)
        bitset_realloc(temp, larger._bitset.size)
        bitset_intersection(temp, temp, larger._bitset)
        retval = bitset_isempty(temp)
        bitset_free(temp)
        return retval

    def __contains__(self, unsigned long n):
        """
        Test to see if the `n` is in self.

        EXAMPLES::

            sage: 0 in FrozenBitset([0,1])
            True
            sage: 0 in FrozenBitset([1,2])
            False
            sage: 10 in FrozenBitset([0,1])
            False
            sage: 121 in FrozenBitset('110'*50)
            True
            sage: 122 in FrozenBitset('110'*50)
            False
        """
        if n < self._bitset.size:
            return bitset_in(self._bitset, n)
        else:
            return False

    def __len__(self):
        """
        Return the number of elements in the bitset (which may be
        different from the capacity of the bitset).

        EXAMPLES::

            sage: len(FrozenBitset([0,1],capacity=10))
            2
            sage: len(FrozenBitset(range(98)))
            98
        """
        return bitset_len(self._bitset)

    def __str__(self):
        """
        Return a string representing the bitset as a binary vector.

        EXAMPLES::

            sage: a=FrozenBitset('10110')
            sage: str(a)
            '10110'
            sage: str(FrozenBitset('110'*32))
            '110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110'
        """
        return bitset_string(self._bitset)

    def __repr__(self):
        """
        Return a string representing the bitset as a binary vector.

        EXAMPLES::

            sage: a=FrozenBitset('10110')
            sage: repr(a)
            '10110'
            sage: repr(FrozenBitset('110'*32))
            '110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110110'
        """
        return self.__str__()

    cpdef _union(self, FrozenBitset other):
        """
        Return the union of self and other.

        In order to get a Cython "union" function, we have to use the
        underscore since "union" is a C keyword.

        EXAMPLES::

            sage: FrozenBitset('10101')._union(FrozenBitset('11100'))
            11101
            sage: FrozenBitset('10101'*10)._union(FrozenBitset('01010'*10))
            11111111111111111111111111111111111111111111111111

        TESTS::

            sage: set(FrozenBitset('10101'*10)._union(FrozenBitset('01010'*10)))==set(FrozenBitset('10101'*10)).union(FrozenBitset('01010'*10))
            True
            sage: set(FrozenBitset('10101')._union(FrozenBitset('01010'*20)))==set(FrozenBitset('10101')).union(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20)._union(FrozenBitset('01010')))==set(FrozenBitset('10101'*20)).union(FrozenBitset('01010'))
            True
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_union(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def union(self, FrozenBitset other):
        """
        Return the union of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101').union(FrozenBitset('11100'))
            11101
            sage: FrozenBitset('10101'*10).union(FrozenBitset('01010'*10))
            11111111111111111111111111111111111111111111111111

        TESTS::

            sage: set(FrozenBitset('10101'*10).union(FrozenBitset('01010'*10)))==set(FrozenBitset('10101'*10)).union(FrozenBitset('01010'*10))
            True
            sage: set(FrozenBitset('10101').union(FrozenBitset('01010'*20)))==set(FrozenBitset('10101')).union(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20).union(FrozenBitset('01010')))==set(FrozenBitset('10101'*20)).union(FrozenBitset('01010'))
            True
        """
        return self._union(other)

    def __or__(self, FrozenBitset other not None):
        """
        Return the union of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101') | FrozenBitset('11100')
            11101
            sage: FrozenBitset('10101'*10) | FrozenBitset('01010'*10)
            11111111111111111111111111111111111111111111111111

       TESTS::

            sage: set(FrozenBitset('10101'*10) | FrozenBitset('01010'*10))==set(FrozenBitset('10101'*10)) | set(FrozenBitset('01010'*10))
            True
            sage: set(FrozenBitset('10101') | FrozenBitset('01010'*20))==set(FrozenBitset('10101')) | set(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20) | FrozenBitset('01010'))==set(FrozenBitset('10101'*20)) | set(FrozenBitset('01010'))
            True
        """
        return self._union(other)

    cpdef intersection(self, FrozenBitset other):
        """
        Return the intersection of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101').intersection(FrozenBitset('11100'))
            10100
            sage: FrozenBitset('11111'*10).intersection(FrozenBitset('010101'*10))
            010101010101010101010101010101010101010101010101010000000000

        TESTS::

            sage: set(FrozenBitset('11111'*10).intersection(FrozenBitset('010101'*10)))==set(FrozenBitset('11111'*10)).intersection(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5).intersection(FrozenBitset('01010'*20)))==set(FrozenBitset('1'*5)).intersection(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20).intersection(FrozenBitset('1'*5)))==set(FrozenBitset('10101'*20)).intersection(FrozenBitset('1'*5))
            True
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_intersection(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def __and__(self, FrozenBitset other not None):
        """
        Return the intersection of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101') & FrozenBitset('11100')
            10100
            sage: FrozenBitset('11111'*10) & FrozenBitset('010101'*10)
            010101010101010101010101010101010101010101010101010000000000

        TESTS::

            sage: set(FrozenBitset('11111'*10) & FrozenBitset('010101'*10))==set(FrozenBitset('11111'*10)) & set(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5) & FrozenBitset('01010'*20))==set(FrozenBitset('1'*5)) & set(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20) & FrozenBitset('1'*5))==set(FrozenBitset('10101'*20)) & set(FrozenBitset('1'*5))
            True
        """
        return self.intersection(other)

    cpdef difference(self, FrozenBitset other):
        """
        Return the difference of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101').difference(FrozenBitset('11100'))
            00001
            sage: FrozenBitset('11111'*10).difference(FrozenBitset('010101'*10))
            101010101010101010101010101010101010101010101010100000000000

        TESTS::

            sage: set(FrozenBitset('11111'*10).difference(FrozenBitset('010101'*10)))==set(FrozenBitset('11111'*10)).difference(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5).difference(FrozenBitset('01010'*20)))==set(FrozenBitset('1'*5)).difference(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20).difference(FrozenBitset('1'*5)))==set(FrozenBitset('10101'*20)).difference(FrozenBitset('1'*5))
            True
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_copy(temp._bitset, self._bitset)

        if temp._bitset.size == other._bitset.size:
            bitset_difference(temp._bitset, temp._bitset, other._bitset)
        elif temp._bitset.size < other._bitset.size:
            bitset_realloc(temp._bitset, other._bitset.size)
            bitset_difference(temp._bitset, temp._bitset, other._bitset)
        else:
            bitset_difference(temp._bitset, temp._bitset, other._larger_capacity_(temp._bitset.size)._bitset)

        return temp

    def __sub__(self, FrozenBitset other not None):
        """
        Return the difference of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101') - FrozenBitset('11100')
            00001
            sage: FrozenBitset('11111'*10)-FrozenBitset('010101'*10)
            101010101010101010101010101010101010101010101010100000000000

        TESTS::

            sage: set(FrozenBitset('11111'*10)-FrozenBitset('010101'*10))==set(FrozenBitset('11111'*10))-set(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5)-FrozenBitset('01010'*20))==set(FrozenBitset('1'*5))-set(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20)-FrozenBitset('1'*5))==set(FrozenBitset('10101'*20))-set(FrozenBitset('1'*5))
            True
        """
        return self.difference(other)

    cpdef symmetric_difference(self, FrozenBitset other):
        """
        Return the symmetric difference of self and other.

        EXAMPLES::

            sage: FrozenBitset('10101').symmetric_difference(FrozenBitset('11100'))
            01001
            sage: FrozenBitset('11111'*10).symmetric_difference(FrozenBitset('010101'*10))
            101010101010101010101010101010101010101010101010100101010101

        TESTS::

            sage: set(FrozenBitset('11111'*10).symmetric_difference(FrozenBitset('010101'*10)))==set(FrozenBitset('11111'*10)).symmetric_difference(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5).symmetric_difference(FrozenBitset('01010'*20)))==set(FrozenBitset('1'*5)).symmetric_difference(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20).symmetric_difference(FrozenBitset('1'*5)))==set(FrozenBitset('10101'*20)).symmetric_difference(FrozenBitset('1'*5))
            True
        """
        if other is None:
            raise ValueError, "other can not be None"
        cdef FrozenBitset temp, smaller, larger
        if self._bitset.size <= other._bitset.size:
            smaller = self
            larger = other
        else:
            smaller = other
            larger = self

        temp = self._new(smaller._bitset.size)
        bitset_copy(temp._bitset, smaller._bitset)
        bitset_realloc(temp._bitset, larger._bitset.size)
        bitset_symmetric_difference(temp._bitset, temp._bitset, larger._bitset)
        return temp

    def __xor__(self, FrozenBitset other not None):
        """
        Return the symmetric difference of self and other.

        Note that because of the Sage preprocessor, in Sage, ``^^`` is the exclusive or, rather than ``^``.

        EXAMPLES::

            sage: FrozenBitset('10101') ^^ FrozenBitset('11100')
            01001
            sage: FrozenBitset('11111'*10) ^^ FrozenBitset('010101'*10)
            101010101010101010101010101010101010101010101010100101010101

        TESTS::

            sage: set(FrozenBitset('11111'*10) ^^ FrozenBitset('010101'*10))==set(FrozenBitset('11111'*10)) ^^ set(FrozenBitset('010101'*10))
            True
            sage: set(FrozenBitset('1'*5) ^^ FrozenBitset('01010'*20))==set(FrozenBitset('1'*5)) ^^ set(FrozenBitset('01010'*20))
            True
            sage: set(FrozenBitset('10101'*20) ^^ FrozenBitset('1'*5))==set(FrozenBitset('10101'*20)) ^^ set(FrozenBitset('1'*5))
            True
        """
        return self.symmetric_difference(other)

    cpdef complement(self):
        """
        Return the complement of self.

        EXAMPLES::

            sage: ~FrozenBitset('10101')
            01010
            sage: ~FrozenBitset('11111'*10)
            00000000000000000000000000000000000000000000000000
            sage: x = FrozenBitset('10'*40)
            sage: x == ~x
            False
            sage: x == ~~x
            True
            sage: x|(~x) == FrozenBitset('11'*40)
            True
            sage: ~x == FrozenBitset('01'*40)
            True
        """
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_complement(temp._bitset, self._bitset)
        return temp

    def __invert__(self):
        """
        Return the complement of self.

        EXAMPLES::

            sage: ~FrozenBitset('10101')
            01010
            sage: ~FrozenBitset('11111'*10)
            00000000000000000000000000000000000000000000000000
            sage: x = FrozenBitset('10'*40)
            sage: x == ~x
            False
            sage: x == ~~x
            True
            sage: x|(~x) == FrozenBitset('11'*40)
            True
            sage: ~x == FrozenBitset('01'*40)
            True
        """
        return self.complement()

    cpdef  __copy__(self):
        """
        Return self (since self is immutable).

        EXAMPLES::

            sage: a=FrozenBitset('10101')
            sage: from copy import copy
            sage: b=copy(a)
            sage: b is a
            True
            sage: c=FrozenBitset('1010'*32)
            sage: d=copy(c)
            sage: d is c
            True
        """
        return self

cdef class Bitset(FrozenBitset):
    """
    A bitset class which leverages inline Cython functions for creating
    and manipulating bitsets.

    A bitset can be thought of in two ways.  First, as a set of
    elements from the universe of the `n` natural numbers `0`, `1`,
    ..., `n-1` (where `n`, the capacity, can be specified), with
    typical set operations (intersection, union, symmetric difference,
    etc.).  Secondly, a bitset can be thought of as a binary vector
    with typical binary operations (i.e., ``and``, ``or``, ``xor``,
    etc.).  This class supports both interfaces.

    The interface in this class mirrors the interface in the ``set``
    datatype of Python.

    .. warning::

        This class is most likely to be useful as a way to store
        Cython bitsets in Python datastructures, acting on them using
        the Cython inline functions.  If you want to use this class
        for a Python set type, the Python ``set`` data type may be
        faster.


    EXAMPLES::

        sage: a=Bitset('1101')
        sage: loads(dumps(a))==a
        True
        sage: a=Bitset('1101'*32)
        sage: loads(dumps(a))==a
        True
    """

    cpdef __copy__(self):
        """
        Return a copy of self.

        EXAMPLES::

            sage: a=Bitset('10101')
            sage: from copy import copy
            sage: b=copy(a)
            sage: b is a
            False
            sage: b==a
            True
            sage: c=Bitset('1010'*32)
            sage: d=copy(c)
            sage: d is c
            False
            sage: d==c
            True

        """
        cdef FrozenBitset temp = self._new(self._bitset.size)
        bitset_copy(temp._bitset, self._bitset)
        return temp

    def __hash__(self):
        """
        Raise an error, since mutable Bitsets are not hashable.

        EXAMPLE::

            sage: hash(Bitset('110'))
            Traceback (most recent call last):
            ...
            TypeError: Bitset objects are unhashable; use FrozenBitset

        """
        raise TypeError, "Bitset objects are unhashable; use FrozenBitset"


    def __richcmp__(FrozenBitset self, FrozenBitset other not None, int op):
        """
        Implement comparisons, using the Cython richcmp convention.
        Comparison is done by inclusion (a set is less than another if
        it is a subset).

        EXAMPLES::

            sage: Bitset('11') < Bitset('101')
            False
            sage: Bitset('11') <= Bitset('110')
            True
            sage: Bitset('11') != Bitset('10')
            True
            sage: Bitset('11') == Bitset('10')
            False
            sage: Bitset('11') == Bitset('110')
            True
            sage: Bitset('11') > Bitset('10')
            True
            sage: Bitset('11') >= Bitset('10')
            True
            sage: FrozenBitset('11') < FrozenBitset('110'*32)
            True
        """
        cdef FrozenBitset left, right

        if self._bitset.size == other._bitset.size:
            left = self
            right = other
        elif self._bitset.size < other._bitset.size:
            left = self._larger_capacity_(other._bitset.size)
            right = other
        else:
            left = self
            right = other._larger_capacity_(self._bitset.size)

        if op == 2: # ==
            return bitset_eq(left._bitset, right._bitset)
        elif op == 3: # !=
            return not bitset_eq(left._bitset, right._bitset)
        elif op == 0: # <
            return bitset_issubset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == 1: # <=
            return bitset_issubset(left._bitset, right._bitset)
        elif op == 4: # >
            return bitset_issuperset(left._bitset, right._bitset) and not bitset_eq(left._bitset, right._bitset)
        elif op == 5: # >=
            return bitset_issuperset(left._bitset, right._bitset)


    cdef FrozenBitset _new(self,long int capacity):
        """
        Return an object of the same type as self, initialized with a bitset of capacity ``capacity``.

        """
        cdef Bitset b
        b = Bitset.__new__(Bitset,None, capacity)
        return b


    cpdef update(self, FrozenBitset other):
        """
        Update the bitset to include items in other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.update(Bitset('0101'))
            sage: a
            1101
            sage: a_set = set(a)
            sage: a.update(Bitset('00011'*25))
            sage: a
            11011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011
            sage: a_set.update(Bitset('00011'*25))
            sage: set(a)==a_set
            True
        """
        if other is None:
            raise TypeError, "other can not be None"
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_union(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_union(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_union(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __ior__(self, FrozenBitset other not None):
        """
        Update the bitset to include items in other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a |= Bitset('0101')
            sage: a
            1101
            sage: a_set = set(a)
            sage: a |= Bitset('00011'*25)
            sage: a
            11011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011000110001100011
            sage: a_set |= set(Bitset('00011'*25))
            sage: set(a)==a_set
            True
        """
        self.update(other)
        return self

    cpdef intersection_update(self, FrozenBitset other):
        """
        Update the bitset to the intersection of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.intersection_update(Bitset('0101'))
            sage: a
            0100
            sage: a_set = set(a)
            sage: a.intersection_update(Bitset('0110'*25))
            sage: a
            0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set.intersection_update(Bitset('0110'*25))
            sage: set(a)==a_set
            True
        """
        if other is None:
            raise TypeError, "other can not be None"
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_intersection(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_intersection(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_intersection(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __iand__(self, FrozenBitset other not None):
        """
        Update the bitset to the intersection of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a &= Bitset('0101')
            sage: a
            0100
            sage: a_set = set(a)
            sage: a &= Bitset('0110'*25)
            sage: a
            0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set &= set(Bitset('0110'*25))
            sage: a_set==set(a)
            True
        """
        self.intersection_update(other)
        return self

    cpdef difference_update(self, FrozenBitset other):
        """
        Update the bitset to the difference of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.difference_update(Bitset('0101'))
            sage: a
            1000
            sage: a_set=set(a)
            sage: a.difference_update(FrozenBitset('010101'*10)); a
            100000000000000000000000000000000000000000000000000000000000
            sage: a_set.difference_update(FrozenBitset('010101'*10))
            sage: a_set==set(a)
            True
            sage: a.difference_update(FrozenBitset('110'))
            sage: a_set.difference_update(FrozenBitset('110'))
            sage: a_set==set(a)
            True
            sage: a.difference_update(FrozenBitset('01010'*20)); a
            0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set.difference_update(FrozenBitset('01010'*20))
            sage: a_set==set(a)
            True
            sage: b=Bitset('10101'*20)
            sage: b_set=set(b)
            sage: b.difference_update(FrozenBitset('1'*5)); b
            0000010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set.difference_update(FrozenBitset('1'*5))
            sage: b_set==set(b)
            True
        """
        if other is None:
            raise TypeError, "other can not be None"
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_difference(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_difference(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_difference(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __isub__(self, FrozenBitset other not None):
        """
        Update the bitset to the difference of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a -= Bitset('0101')
            sage: a
            1000
            sage: a_set=set(a)
            sage: a-=FrozenBitset('010101'*10); a
            100000000000000000000000000000000000000000000000000000000000
            sage: a_set-=set(FrozenBitset('010101'*10))
            sage: a_set==set(a)
            True
            sage: a=Bitset('110')
            sage: a_set=set(a)
            sage: a-=FrozenBitset('01010'*20); a
            1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a_set-=set(FrozenBitset('01010'*20))
            sage: a_set==set(a)
            True
            sage: b=FrozenBitset('10101'*20)
            sage: b_set = set(b)
            sage: b -= FrozenBitset('1'*5); b
            0000010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set -= FrozenBitset('1'*5)
            sage: b_set==set(b)
            True
        """
        self.difference_update(other)
        return self

    cpdef symmetric_difference_update(self, FrozenBitset other):
        """
        Update the bitset to the symmetric difference of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.symmetric_difference_update(Bitset('0101'))
            sage: a
            1001
            sage: a_set=set(a)
            sage: a.symmetric_difference_update(FrozenBitset('010101'*10)); a
            110001010101010101010101010101010101010101010101010101010101
            sage: a_set.symmetric_difference_update(FrozenBitset('010101'*10))
            sage: a_set==set(a)
            True
            sage: a.symmetric_difference_update(FrozenBitset('01010'*20)); a
            1001011111000001111100000111110000011111000001111100000111110101001010010100101001010010100101001010
            sage: a_set.symmetric_difference_update(FrozenBitset('01010'*20))
            sage: a_set==set(a)
            True
            sage: b=Bitset('10101'*20)
            sage: b_set=set(b)
            sage: b.symmetric_difference_update( FrozenBitset('1'*5)); b
            0101010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set.symmetric_difference_update( FrozenBitset('1'*5))
            sage: b_set==set(b)
            True
        """
        if other is None:
            raise TypeError, "other can not be None"
        cdef bitset_t temp
        if self._bitset.size == other._bitset.size:
            bitset_symmetric_difference(self._bitset, self._bitset, other._bitset)
        elif self._bitset.size < other._bitset.size:
            bitset_realloc(self._bitset, other._bitset.size)
            bitset_symmetric_difference(self._bitset, self._bitset, other._bitset)
        else:
            bitset_init(temp, other._bitset.size)
            bitset_copy(temp, other._bitset)
            bitset_realloc(temp, self._bitset.size)
            bitset_symmetric_difference(self._bitset, self._bitset, temp)
            bitset_free(temp)

    def __ixor__(self, FrozenBitset other not None):
        """
        Update the bitset to the symmetric difference of self and other.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a ^^=Bitset('0101')
            sage: a
            1001
            sage: a_set=set(a)
            sage: a ^^=FrozenBitset('010101'*10); a
            110001010101010101010101010101010101010101010101010101010101
            sage: a_set ^^=set(FrozenBitset('010101'*10))
            sage: a_set==set(a)
            True
            sage: a ^^=FrozenBitset('01010'*20); a
            1001011111000001111100000111110000011111000001111100000111110101001010010100101001010010100101001010
            sage: a_set ^^=set(FrozenBitset('01010'*20))
            sage: a_set==set(a)
            True
            sage: b=Bitset('10101'*20)
            sage: b_set=set(b)
            sage: b ^^=FrozenBitset('1'*5); b
            0101010101101011010110101101011010110101101011010110101101011010110101101011010110101101011010110101
            sage: b_set ^^=set(FrozenBitset('1'*5))
            sage: b_set==set(b)
            True
        """
        self.symmetric_difference_update(other)
        return self

    cpdef add(self, unsigned long n):
        """
        Update the bitset by adding `n`.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.add(5)
            sage: a
            110001
            sage: a.add(100)
            sage: sorted(list(a))
            [0, 1, 5, 100]
            sage: a.capacity()
            101
        """
        if n >= self._bitset.size:
            bitset_realloc(self._bitset, n+1)
        bitset_add(self._bitset, n)

    cpdef remove(self, unsigned long n):
        """
        Update the bitset by removing `n`.  Raises KeyError if `n` is
        not contained in the bitset.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.remove(1)
            sage: a
            100
            sage: a.remove(2)
            Traceback (most recent call last):
            ...
            KeyError: 2L
            sage: a.remove(4)
            Traceback (most recent call last):
            ...
            KeyError: 4L
            sage: a
            100
            sage: a=Bitset('000001'*15); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89]
            sage: a.remove(83); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]
        """
        if n >= self._bitset.size:
            raise KeyError, n
        else:
            bitset_remove(self._bitset, n)

    cpdef discard(self, unsigned long n):
        """
        Update the bitset by removing `n`.

        EXAMPLES::

            sage: a = Bitset('110')
            sage: a.discard(1)
            sage: a
            100
            sage: a.discard(2)
            sage: a.discard(4)
            sage: a
            100
            sage: a=Bitset('000001'*15); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89]
            sage: a.discard(83); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]
            sage: a.discard(82); sorted(list(a))
            [5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 89]
        """
        if n < self._bitset.size:
            bitset_discard(self._bitset, n)


    cpdef pop(self):
        """
        Remove and return an arbitrary element from the set. Raises
        KeyError if the set is empty.

        EXAMPLES::

            sage: a = Bitset('011')
            sage: a.pop()
            1
            sage: a
            001
            sage: a.pop()
            2
            sage: a
            000
            sage: a.pop()
            Traceback (most recent call last):
            ...
            KeyError: 'pop from an empty set'
            sage: a = Bitset('0001'*32)
            sage: a.pop()
            3
            sage: [a.pop() for _ in range(20)]
            [7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83]
        """
        return bitset_pop(self._bitset)


    cpdef clear(self):
        """
        Removes all elements from the set

        EXAMPLES::

            sage: a = Bitset('011')
            sage: a.clear()
            sage: a
            000
            sage: a = Bitset('011'*32)
            sage: a.clear()
            sage: set(a)
            set([])
        """
        bitset_clear(self._bitset)
