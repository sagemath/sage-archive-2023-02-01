# coding=utf-8
"""
Alphabets
"""
#*****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.combinat import CombinatorialClass
from sage.misc.misc import uniq
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
import itertools

def Alphabet(data=None, name=None):
    r"""
    Returns an object representing an ordered alphabet.

    EXAMPLES::

        sage: Alphabet("ab")
        Ordered Alphabet ['a', 'b']
        sage: Alphabet([0,1,2])
        Ordered Alphabet [0, 1, 2]
        sage: Alphabet(name="positive integers")
        Ordered Alphabet of Positive Integers
        sage: Alphabet(name="PP")
        Ordered Alphabet of Positive Integers
        sage: Alphabet(name="natural numbers")
        Ordered Alphabet of Natural Numbers
        sage: Alphabet(name="NN")
        Ordered Alphabet of Natural Numbers
    """
    if data is None and name is None:
        raise TypeError, "provide at least one argument"
    if data is None:
        if name == "positive integers" or name == "PP":
            return OrderedAlphabet_PositiveIntegers()
        elif name == "natural numbers" or name == "NN":
            return OrderedAlphabet_NaturalNumbers()
        else:
            raise TypeError, "name is not recognized"
    else:
        try:
            return OrderedAlphabet_Finite(data)
        except:
            raise TypeError, "cannot construct an alphabet from given data"

OrderedAlphabet = Alphabet

class OrderedAlphabet_class(CombinatorialClass):
    r"""
    Generic class for ordered alphabets.
    """
    pass

class OrderedAlphabet_Finite(OrderedAlphabet_class):
    def __init__(self, alphabet):
        """
        Builds an ordered alphabet from an iterable. The order is that
        given by the order the items appear in the iterable. There must be
        no duplicates.

        NOTE: The alphabet is expanded in memory and stored as a list.

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: A = OrderedAlphabet_Finite([0,1,2])
            sage: A == loads(dumps(A))
            True
            sage: A = OrderedAlphabet_Finite("abc")
            sage: A == loads(dumps(A))
            True

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('aba')
            Traceback (most recent call last):
            ...
            ValueError: duplicate elements in alphabet
            sage: OrderedAlphabet_Finite(33)
            Traceback (most recent call last):
            ...
            TypeError: cannot build an alphabet from 33
        """
        try:
            self._alphabet = list(alphabet)
        except TypeError:
            raise TypeError, "cannot build an alphabet from %s" % (alphabet)
        if len(uniq(self._alphabet)) != len(self._alphabet):
            raise ValueError, "duplicate elements in alphabet"

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite([0,1,2]).__repr__()
            'Ordered Alphabet [0, 1, 2]'
            sage: OrderedAlphabet_Finite("cba").__repr__()
            "Ordered Alphabet ['c', 'b', 'a']"
        """
        return "Ordered Alphabet " + str(self._alphabet)

    def iterator(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: type(OrderedAlphabet_Finite("abc").iterator())
            <type 'generator'>
            sage: list(OrderedAlphabet_Finite("abc").iterator())
            ['a', 'b', 'c']
            sage: list(OrderedAlphabet_Finite([10, 17, 3]).iterator())
            [10, 17, 3]
        """
        for a in self._alphabet:
            yield a

    def __contains__(self, a):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: A = OrderedAlphabet_Finite([0,1,2])
            sage: 2 in A
            True
            sage: 17 in A
            False
            sage: B = OrderedAlphabet_Finite("abc")
            sage: "b" in B
            True
            sage: "z" in B
            False
        """
        try:
            return a in self._alphabet
        except:
            return False

    def __le__(self, other):
        r"""
        Returns True if the elements of self appear among the elements of
        other in the same respective order, and False otherwise.

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('abc') <= OrderedAlphabet_Finite('daefbgc')
            True
            sage: OrderedAlphabet_Finite('abc') <= OrderedAlphabet_Finite('abc')
            True
            sage: OrderedAlphabet_Finite('abc') <= OrderedAlphabet_Finite('ac')
            False
            sage: OrderedAlphabet_Finite('abc') <= OrderedAlphabet_Finite([1, 2, 3])
            False
            sage: OrderedAlphabet_Finite('abc') <= OrderedAlphabet_Finite("bca")
            False
        """
        return self.list() == other.filter(lambda a: a in self).list()

    def __ge__(self, other):
        r"""
        Returns True if the elements of other appear among the elements of
        self in the same respective order, and False otherwise.

        The ordering of the alphabet is taken into consideration.

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('abc') >= OrderedAlphabet_Finite('daefbgc')
            False
            sage: OrderedAlphabet_Finite('abc') >= OrderedAlphabet_Finite('abc')
            True
            sage: OrderedAlphabet_Finite('abc') >= OrderedAlphabet_Finite('ac')
            True
            sage: OrderedAlphabet_Finite('abc') >= OrderedAlphabet_Finite('cba')
            False
            sage: OrderedAlphabet_Finite('abc') >= OrderedAlphabet_Finite([1, 2, 3])
            False
        """
        return other.list() == self.filter(lambda a: a in other).list()

    def string_rep(self):
        r"""
        Returns the string representation of the alphabet.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('cba').string_rep()
            "['c', 'b', 'a']"
            sage: OrderedAlphabet_Finite([1, 3, 2]).string_rep()
            '[1, 3, 2]'
        """
        return "[%s]" % ', '.join(itertools.imap(repr, self))

    def rank(self, letter):
        r"""
        Returns the index of letter in self.

        INPUT:


        -  ``letter`` - a letter contained in this alphabet


        OUTPUT:


        -  ``integer`` - the integer mapping for the letter


        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('abcd').rank('a')
            0
            sage: OrderedAlphabet_Finite('abcd').rank('d')
            3
            sage: OrderedAlphabet_Finite('abcd').rank('e')
            Traceback (most recent call last):
            ...
            IndexError: letter not in alphabet: 'e'
            sage: OrderedAlphabet_Finite('abcd').rank('')
            Traceback (most recent call last):
            ...
            IndexError: letter not in alphabet: ''
        """
        try:
            return self._alphabet.index(letter)
        except ValueError:
            raise IndexError, "letter not in alphabet: %s" % repr(letter)

    def unrank(self, n):
        r"""
        Returns the letter in position n of the alphabet self.

        INPUT:


        -  ``n`` - a nonnegative integer


        OUTPUT: the (n+1)-th object output by iter(self)

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Finite
            sage: OrderedAlphabet_Finite('abcd').unrank(0)
            'a'
            sage: OrderedAlphabet_Finite('abcd').unrank(3)
            'd'
            sage: OrderedAlphabet_Finite('abcd').unrank(5)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._alphabet[n]

class OrderedAlphabet_Infinite(OrderedAlphabet_class):
    def __le__(self, other):
        r"""
        Returns NotimplementedError since it is not clear how to define
        this for infinite ordered alphabets.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Infinite
            sage: A1 = OrderedAlphabet_Infinite()
            sage: A2 = OrderedAlphabet_Infinite()
            sage: A1.__le__(A2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __ge__(self, other):
        r"""
        Returns NotimplementedError since it is not clear how to define
        this for infinite ordered alphabets.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Infinite
            sage: A1 = OrderedAlphabet_Infinite()
            sage: A2 = OrderedAlphabet_Infinite()
            sage: A1.__le__(A2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def count(self):
        r"""
        Return the number of elements in self.

        OUTPUT: +Infinity

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Infinite
            sage: OrderedAlphabet_Infinite().count()
            +Infinity
        """
        return Infinity

    def list(self):
        r"""
        Returns NotImplementedError since we cannot list all the
        nonnegative integers.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_Infinite
            sage: OrderedAlphabet_Infinite().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class OrderedAlphabet_NaturalNumbers(OrderedAlphabet_Infinite):
    r"""
    The alphabet of nonnegative integers, ordered in the usual way.

    TESTS::

        sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
        sage: NN = OrderedAlphabet_NaturalNumbers()
        sage: NN == loads(dumps(NN))
        True
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: OrderedAlphabet_NaturalNumbers().__repr__()
            'Ordered Alphabet of Natural Numbers'
        """
        return "Ordered Alphabet of Natural Numbers"

    def iterator(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: it = OrderedAlphabet_NaturalNumbers().iterator()
            sage: type(it)
            <type 'generator'>
            sage: it.next()
            0
            sage: it.next()
            1
            sage: it.next()
            2
        """
        for i in itertools.count(0):
            yield Integer(i)

    def __contains__(self, a):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: A = OrderedAlphabet_NaturalNumbers()
            sage: 2 in A
            True
            sage: 17 in A
            True
            sage: 0 in A
            True
            sage: -1 in A
            False
            sage: "z" in A
            False
        """
        return isinstance(a, (int, Integer)) and a >= 0

    def string_rep(self):
        r"""
        Returns the string representation of the alphabet.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: OrderedAlphabet_NaturalNumbers().string_rep()
            'Natural Numbers'
        """
        return "Natural Numbers"

    def rank(self, letter):
        r"""
        Returns the index of letter in self.

        INPUT:


        -  ``letter`` - a letter contained in this alphabet


        OUTPUT:


        -  ``integer`` - the integer mapping for the letter


        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: NN = OrderedAlphabet_NaturalNumbers()
            sage: NN.rank(0)
            0
            sage: NN.rank(17)
            17

        TESTS::

            sage: NN.rank(-1)
            Traceback (most recent call last):
            ...
            ValueError: letter(=-1) not in the alphabet
            sage: NN.rank("a")
            Traceback (most recent call last):
            ...
            ValueError: letter(=a) not in the alphabet
        """
        if letter in self:
            return letter
        else:
            raise ValueError, "letter(=%s) not in the alphabet" % letter

    def unrank(self, n):
        r"""
        Returns the letter in position n in self, which in this case is n.

        INPUT:


        -  ``n`` - nonnegative integer


        OUTPUT:


        -  ``n`` - nonnegative integer


        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: NN = OrderedAlphabet_NaturalNumbers()
            sage: NN.unrank(0)
            0
            sage: NN.unrank(17)
            17

        TESTS::

            sage: NN.unrank(-1)
            Traceback (most recent call last):
            ...
            ValueError: -1 is not a nonnegative integer
            sage: NN.unrank("a")
            Traceback (most recent call last):
            ...
            ValueError: a is not a nonnegative integer
        """
        if isinstance(n, (int, Integer)) and n >= 0:
            return n
        else:
            raise ValueError, "%s is not a nonnegative integer" % n

    def next(self, n):
        r"""
        Returns the letter following n in the alphabet self.

        INPUT:


        -  ``n`` - nonnegative integer


        OUTPUT: n+1

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_NaturalNumbers
            sage: NN = OrderedAlphabet_NaturalNumbers()
            sage: NN.next(0)
            1
            sage: NN.next(117)
            118
            sage: NN.next(-1)
            Traceback (most recent call last):
            ...
            ValueError: letter(=-1) not in the alphabet
        """
        if n in self:
            return Integer(n+1)
        else:
            raise ValueError, "letter(=%s) not in the alphabet" % n

class OrderedAlphabet_PositiveIntegers(OrderedAlphabet_Infinite):
    r"""
    The alphabet of nonnegative integers, ordered in the usual way.

    TESTS::

        sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
        sage: PP = OrderedAlphabet_PositiveIntegers()
        sage: PP == loads(dumps(PP))
        True
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: OrderedAlphabet_PositiveIntegers().__repr__()
            'Ordered Alphabet of Positive Integers'
        """
        return "Ordered Alphabet of Positive Integers"

    def iterator(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: it = OrderedAlphabet_PositiveIntegers().iterator()
            sage: type(it)
            <type 'generator'>
            sage: it.next()
            1
            sage: it.next()
            2
            sage: it.next()
            3
        """
        for i in itertools.count(1):
            yield Integer(i)

    def __contains__(self, a):
        """
        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: A = OrderedAlphabet_PositiveIntegers()
            sage: 2 in A
            True
            sage: 17 in A
            True
            sage: 0 in A
            False
            sage: -1 in A
            False
            sage: "z" in A
            False
        """
        return isinstance(a, (int, Integer)) and a > 0

    def string_rep(self):
        r"""
        Returns the string representation of the alphabet.

        TESTS::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: OrderedAlphabet_PositiveIntegers().string_rep()
            'Positive Integers'
        """
        return "Positive Integers"

    def rank(self, letter):
        r"""
        Returns the index of letter in self.

        INPUT:


        -  ``letter`` - a positive integer


        OUTPUT: letter-1

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: OrderedAlphabet_PositiveIntegers().rank(1)
            0
            sage: OrderedAlphabet_PositiveIntegers().rank(8)
            7

        TESTS::

            sage: OrderedAlphabet_PositiveIntegers().rank(-1)
            Traceback (most recent call last):
            ...
            TypeError: -1 not in alphabet
        """
        if letter in self:
            return letter-1
        else:
            raise TypeError, "%s not in alphabet" % letter

    def unrank(self, i):
        r"""
        Returns the i-th letter in self, where the first letter is the 0-th
        letter.

        INPUT:


        -  ``i`` - an integer


        OUTPUT: i+1

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: OrderedAlphabet_PositiveIntegers().unrank(0)
            1
            sage: OrderedAlphabet_PositiveIntegers().unrank(7)
            8
        """
        if i == 0 or i in self:
            return i+1
        elif isinstance(i, slice):
            return range(1,1+i.stop)[i]
        else:
            raise IndexError

    def next(self, n):
        r"""
        Returns the letter following n in the alphabet self.

        INPUT:


        -  ``n`` - positive integer


        OUTPUT: n+1

        EXAMPLES::

            sage: from sage.combinat.words.alphabet import OrderedAlphabet_PositiveIntegers
            sage: PP = OrderedAlphabet_PositiveIntegers()
            sage: PP.next(1)
            2
            sage: PP.next(117)
            118
            sage: PP.next(0)
            Traceback (most recent call last):
            ...
            ValueError: letter(=0) not in the alphabet
        """
        if n in self:
            return Integer(n+1)
        else:
            raise ValueError, "letter(=%s) not in the alphabet" % n
