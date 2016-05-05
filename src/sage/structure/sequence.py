r"""
Finite Homogenous Sequences

A mutable sequence of elements with a common guaranteed category,
which can be set immutable.

Sequence derives from list, so has all the functionality of lists and
can be used wherever lists are used.  When a sequence is created
without explicitly given the common universe of the elements, the
constructor coerces the first and second element to some
*canonical* common parent, if possible, then the second and
third, etc.  If this is possible, it then coerces everything into the
canonical parent at the end.  (Note that canonical coercion is very
restrictive.)  The sequence then has a function ``universe()``
which returns either the common canonical parent (if the coercion
succeeded), or the category of all objects (Objects()).  So if you
have a list `v` and type::

    sage: v = [1, 2/3, 5]
    sage: w = Sequence(v)
    sage: w.universe()
    Rational Field

then since ``w.universe()`` is `\QQ`, you're guaranteed that all
elements of `w` are rationals::

    sage: v[0].parent()
    Integer Ring
    sage: w[0].parent()
    Rational Field

If you do assignment to `w` this property of being rationals is guaranteed
to be preserved::

    sage: w[0] = 2
    sage: w[0].parent()
    Rational Field
    sage: w[0] = 'hi'
    Traceback (most recent call last):
    ...
    TypeError: unable to convert 'hi' to a rational

However, if you do ``w = Sequence(v)`` and the resulting universe
is ``Objects()``, the elements are not guaranteed to have any
special parent.  This is what should happen, e.g., with finite field
elements of different characteristics::

    sage: v = Sequence([GF(3)(1), GF(7)(1)])
    sage: v.universe()
    Category of objects

You can make a list immutable with ``v.freeze()``.  Assignment is
never again allowed on an immutable list.

Creation of a sequence involves making a copy of the input list, and
substantial coercions.  It can be greatly sped up by explicitly
specifying the universe of the sequence::

    sage: v = Sequence(range(10000), universe=ZZ)

TESTS::

    sage: v = Sequence([1..5])
    sage: loads(dumps(v)) == v
    True

"""


##########################################################################
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################


from sage.misc.latex import list_function as list_latex_function
import sage.structure.sage_object

#from mutability import Mutability #we cannot inherit from Mutability and list at the same time

def Sequence(x, universe=None, check=True, immutable=False, cr=False, cr_str=None, use_sage_types=False):
    """
    A mutable list of elements with a common guaranteed universe,
    which can be set immutable.

    A universe is either an object that supports coercion (e.g., a
    parent), or a category.

    INPUT:

    - ``x`` - a list or tuple instance

    - ``universe`` - (default: None) the universe of elements; if None
      determined using canonical coercions and the entire list of
      elements.  If list is empty, is category Objects() of all
      objects.

    - ``check`` -- (default: True) whether to coerce the elements of x
      into the universe

    - ``immutable`` - (default: True) whether or not this sequence is
      immutable

    - ``cr`` - (default: False) if True, then print a carriage return
      after each comma when printing this sequence.

    - ``cr_str`` - (default: False) if True, then print a carriage return
      after each comma when calling ``str()`` on this sequence.

    - ``use_sage_types`` -- (default: False) if True, coerce the
       built-in Python numerical types int, long, float, complex to the
       corresponding Sage types (this makes functions like vector()
       more flexible)

    OUTPUT:

    - a sequence

    EXAMPLES::

        sage: v = Sequence(range(10))
        sage: v.universe()
        <type 'int'>
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We can request that the built-in Python numerical types be coerced
    to Sage objects::

        sage: v = Sequence(range(10), use_sage_types=True)
        sage: v.universe()
        Integer Ring
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    You can also use seq for "Sequence", which is identical to using
    Sequence::

        sage: v = seq([1,2,1/1]); v
        [1, 2, 1]
        sage: v.universe()
        Rational Field

    Note that assignment coerces if possible,::

        sage: v = Sequence(range(10), ZZ)
        sage: a = QQ(5)
        sage: v[3] = a
        sage: parent(v[3])
        Integer Ring
        sage: parent(a)
        Rational Field
        sage: v[3] = 2/3
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Sequences can be used absolutely anywhere lists or tuples can be used::

        sage: isinstance(v, list)
        True

    Sequence can be immutable, so entries can't be changed::

        sage: v = Sequence([1,2,3], immutable=True)
        sage: v.is_immutable()
        True
        sage: v[0] = 5
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    Only immutable sequences are hashable (unlike Python lists),
    though the hashing is potentially slow, since it first involves
    conversion of the sequence to a tuple, and returning the hash of
    that.::

        sage: v = Sequence(range(10), ZZ, immutable=True)
        sage: hash(v)
        1591723448             # 32-bit
        -4181190870548101704   # 64-bit


    If you really know what you are doing, you can circumvent the type
    checking (for an efficiency gain)::

        sage: list.__setitem__(v, int(1), 2/3)        # bad circumvention
        sage: v
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list.__setitem__(v, int(1), int(2))     # not so bad circumvention

    You can make a sequence with a new universe from an old sequence.::

        sage: w = Sequence(v, QQ)
        sage: w
        [0, 2, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: w.universe()
        Rational Field
        sage: w[1] = 2/3
        sage: w
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]

    The default universe for any sequence, if no compatible parent structure
    can be found, is the universe of all Sage objects.

    This example illustrates how every element of a list is taken into account
    when constructing a sequence.::

        sage: v = Sequence([1,7,6,GF(5)(3)]); v
        [1, 2, 1, 3]
        sage: v.universe()
        Finite Field of size 5
    """
    from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal


    if isinstance(x, Sequence_generic) and universe is None:
        universe = x.universe()
        x = list(x)

    if isinstance(x, MPolynomialIdeal) and universe is None:
        universe = x.ring()
        x = x.gens()

    if universe is None:
        if not isinstance(x, (list, tuple)):
            x = list(x)
            #raise TypeError("x must be a list or tuple")

        if len(x) == 0:
            import sage.categories.all
            universe = sage.categories.all.Objects()
        else:
            import sage.structure.element as coerce
            y = x
            x = list(x)   # make a copy, or we'd change the type of the elements of x, which would be bad.
            if use_sage_types:
                # convert any Python built-in numerical types to Sage objects
                from sage.rings.integer_ring import ZZ
                from sage.rings.real_double import RDF
                from sage.rings.complex_double import CDF
                for i in range(len(x)):
                    if isinstance(x[i], int) or isinstance(x[i], long):
                        x[i] = ZZ(x[i])
                    elif isinstance(x[i], float):
                        x[i] = RDF(x[i])
                    elif isinstance(x[i], complex):
                        x[i] = CDF(x[i])
            # start the pairwise coercion
            for i in range(len(x)-1):
                try:
                    x[i], x[i+1] = coerce.canonical_coercion(x[i],x[i+1])
                except TypeError:
                    import sage.categories.all
                    universe = sage.categories.all.Objects()
                    x = list(y)
                    check = False  # no point
                    break
            if universe is None:   # no type errors raised.
                universe = coerce.parent(x[len(x)-1])

    from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
    from sage.rings.polynomial.pbori import BooleanMonomialMonoid
    from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
    from sage.rings.quotient_ring import is_QuotientRing

    if is_MPolynomialRing(universe) or isinstance(universe, BooleanMonomialMonoid) or (is_QuotientRing(universe) and is_MPolynomialRing(universe.cover_ring())):
        return PolynomialSequence(x, universe, immutable=immutable, cr=cr, cr_str=cr_str)
    else:
        return Sequence_generic(x, universe, check, immutable, cr, cr_str, use_sage_types)


class Sequence_generic(sage.structure.sage_object.SageObject, list):
    """
    A mutable list of elements with a common guaranteed universe,
    which can be set immutable.

    A universe is either an object that supports coercion (e.g., a parent),
    or a category.

    INPUT:

    - ``x`` - a list or tuple instance

    - ``universe`` - (default: None) the universe of elements; if None
      determined using canonical coercions and the entire list of
      elements.  If list is empty, is category Objects() of all
      objects.

    - ``check`` -- (default: True) whether to coerce the elements of x
      into the universe

    - ``immutable`` - (default: True) whether or not this sequence is
      immutable

    - ``cr`` - (default: False) if True, then print a carriage return
      after each comma when printing this sequence.

    - ``use_sage_types`` -- (default: False) if True, coerce the
       built-in Python numerical types int, long, float, complex to the
       corresponding Sage types (this makes functions like vector()
       more flexible)

    OUTPUT:

    - a sequence

    EXAMPLES::

        sage: v = Sequence(range(10))
        sage: v.universe()
        <type 'int'>
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We can request that the built-in Python numerical types be coerced
    to Sage objects::

        sage: v = Sequence(range(10), use_sage_types=True)
        sage: v.universe()
        Integer Ring
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    You can also use seq for "Sequence", which is identical to using Sequence::

        sage: v = seq([1,2,1/1]); v
        [1, 2, 1]
        sage: v.universe()
        Rational Field

    Note that assignment coerces if possible,

    ::

        sage: v = Sequence(range(10), ZZ)
        sage: a = QQ(5)
        sage: v[3] = a
        sage: parent(v[3])
        Integer Ring
        sage: parent(a)
        Rational Field
        sage: v[3] = 2/3
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Sequences can be used absolutely anywhere lists or tuples can be used::

        sage: isinstance(v, list)
        True

    Sequence can be immutable, so entries can't be changed::

        sage: v = Sequence([1,2,3], immutable=True)
        sage: v.is_immutable()
        True
        sage: v[0] = 5
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    Only immutable sequences are hashable (unlike Python lists),
    though the hashing is potentially slow, since it first involves
    conversion of the sequence to a tuple, and returning the hash of
    that.

    ::

        sage: v = Sequence(range(10), ZZ, immutable=True)
        sage: hash(v)
        1591723448             # 32-bit
        -4181190870548101704   # 64-bit


    If you really know what you are doing, you can circumvent the type
    checking (for an efficiency gain)::

        sage: list.__setitem__(v, int(1), 2/3)        # bad circumvention
        sage: v
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list.__setitem__(v, int(1), int(2))     # not so bad circumvention

    You can make a sequence with a new universe from an old sequence.

    ::

        sage: w = Sequence(v, QQ)
        sage: w
        [0, 2, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: w.universe()
        Rational Field
        sage: w[1] = 2/3
        sage: w
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]

    The default universe for any sequence, if no compatible parent structure
    can be found, is the universe of all Sage objects.

    This example illustrates how every element of a list is taken into account
    when constructing a sequence.

    ::

        sage: v = Sequence([1,7,6,GF(5)(3)]); v
        [1, 2, 1, 3]
        sage: v.universe()
        Finite Field of size 5

    """
    def __init__(self, x, universe=None, check=True, immutable=False,
                 cr=False, cr_str=None, use_sage_types=False):
        """
        Create a sequence.

        EXAMPLES::

            sage: Sequence([1..5])
            [1, 2, 3, 4, 5]
            sage: a = Sequence([1..3], universe=QQ, check=False, immutable=True, cr=True, cr_str=False, use_sage_types=True)
            sage: a
            [
            1,
            2,
            3
            ]
            sage: a = Sequence([1..5], universe=QQ, check=False, immutable=True, cr_str=True, use_sage_types=True)
            sage: a
            [1, 2, 3, 4, 5]
            sage: a._Sequence_generic__cr_str
            True
            sage: a.__str__()
            '[\n1,\n2,\n3,\n4,\n5\n]'
        """
        self.__hash = None

        self.__cr = cr
        if cr_str is None:
            self.__cr_str = cr
        else:
            self.__cr_str = cr_str

        if isinstance(x, Sequence_generic):
            if universe is None or universe == x.__universe:
                list.__init__(self, x)
                self.__universe = x.__universe
                self._is_immutable = immutable
                return

        self.__universe = universe
        if check:
            x = [universe(t) for t in x]
        list.__init__(self, x)
        self._is_immutable = immutable

    def reverse(self):
        """
        Reverse the elements of self, in place.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.reverse(); B
            [3, 2, 1]
        """
        self._require_mutable()
        list.reverse(self)

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: a = Sequence([1..5])
            sage: a[2] = 19
            sage: a
            [1, 2, 19, 4, 5]
            sage: a[2] = 'hello'
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'hello' to an integer
            sage: a[2] = '5'
            sage: a
            [1, 2, 5, 4, 5]
            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: v[1:3] = [5,7]
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: v = Sequence([1,2,3,4])
            sage: v[1:3] = [5, 3/1]
            sage: v
            [1, 5, 3, 4]
            sage: type(v[2])
            <type 'sage.rings.integer.Integer'>
        """
        self._require_mutable()
        if isinstance(n, slice):
            y = [self.__universe(x) for x in value]
        else:
            y = self.__universe(value)
        list.__setitem__(self, n, y)
        self.__hash=None

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: w = v[2:]
            sage: w
            [3, 4]
            sage: type(w)
            <class 'sage.structure.sequence.Sequence_generic'>
            sage: w[0] = 5; w
            [5, 4]
            sage: v
            [1, 2, 3, 4]
        """
        if isinstance(n, slice):
            return Sequence(list.__getitem__(self, n),
                            universe = self.__universe,
                            check = False,
                            immutable = False,
                            cr = self.__cr)
        else:
            return list.__getitem__(self,n)

    # We have to define the *slice functions as long as Sage uses Python 2.*
    # otherwise the inherited *slice functions from list are called
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i,j))

    def __setslice__(self, i, j, value):
        return self.__setitem__(slice(i,j), value)

    def append(self, x):
        """
        EXAMPLES::

            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: v.append(34)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: v = Sequence([1/3,2,3,4])
            sage: v.append(4)
            sage: type(v[4])
            <type 'sage.rings.rational.Rational'>
        """
        self._require_mutable()
        y = self.__universe(x)
        list.append(self, y)

    def extend(self, iterable):
        """
        Extend list by appending elements from the iterable.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.extend(range(4))
            sage: B
            [1, 2, 3, 0, 1, 2, 3]
        """
        self._require_mutable()
        v = [self.__universe(x) for x in iterable]
        list.extend(self, v)

    def insert(self, index, object):
        """
        Insert object before index.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.insert(10, 5)
            sage: B
            [1, 2, 3, 5]
        """
        self._require_mutable()
        list.insert(self, index, self.__universe(object))

    def pop(self, index=-1):
        """
        Remove and return item at index (default last)

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.pop(1)
            2
            sage: B
            [1, 3]
        """
        self._require_mutable()
        return list.pop(self, index)

    def remove(self, value):
        """
        Remove first occurrence of value

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.remove(2)
            sage: B
            [1, 3]
        """
        self._require_mutable()
        list.remove(self, value)

    def sort(self, cmp=None, key=None, reverse=False):
        """
        Sort this list *IN PLACE*.

        INPUT:

        - ``cmp`` - see Python ``list sort``
        
        - ``key`` - see Python ``list sort``
        
        - ``reverse`` - see Python ``list sort``

        EXAMPLES::

            sage: B = Sequence([3,2,1/5])
            sage: B.sort()
            sage: B
            [1/5, 2, 3]
            sage: B.sort(reverse=True); B
            [3, 2, 1/5]
            sage: B.sort(cmp = lambda x,y: cmp(y,x)); B
            [3, 2, 1/5]
            sage: B.sort(cmp = lambda x,y: cmp(y,x), reverse=True); B
            [1/5, 2, 3]
        """
        self._require_mutable()
        list.sort(self, cmp=cmp, key=key, reverse=reverse)

    def __hash__(self):
        """
        EXAMPLES::

            sage: a = Sequence([1..5])
            sage: a.__hash__()
            Traceback (most recent call last):
            ...
            ValueError: mutable sequences are unhashable
            sage: a[0] = 10
            sage: a.set_immutable()
            sage: a.__hash__()
            -123014399  # 32-bit
            -5823618793256324351  # 64-bit
            sage: hash(a)
            -123014399  # 32-bit
            -5823618793256324351  # 64-bit
        """
        if not self._is_immutable:
            raise ValueError("mutable sequences are unhashable")
        if self.__hash is None:
            self.__hash = hash(tuple(self))
        return self.__hash

    def _repr_(self):
        """
        EXAMPLES::

            sage: Sequence([1,2/3,-2/5])._repr_()
            '[1, 2/3, -2/5]'
            sage: print Sequence([1,2/3,-2/5], cr=True)._repr_()
            [
            1,
            2/3,
            -2/5
            ]
        """
        if self.__cr:
            return '[\n' + ',\n'.join([repr(x) for x in self]) + '\n]'
        else:
            return list.__repr__(self)

    def _latex_(self):
        r"""
        TESTS::

            sage: t= Sequence([sqrt(x), exp(x), x^(x-1)], universe=SR); t
            [sqrt(x), e^x, x^(x - 1)]
            sage: t._latex_()
            '\\left[\\sqrt{x}, e^{x}, x^{x - 1}\\right]'
            sage: latex(t)
            \left[\sqrt{x}, e^{x}, x^{x - 1}\right]
        """
        return list_latex_function(self)

    def __str__(self):
        """
        EXAMPLES::

            sage: s = Sequence([1,2,3], cr=False)
            sage: str(s)
            '[1, 2, 3]'
            sage: repr(s)
            '[1, 2, 3]'
            sage: print s
            [1, 2, 3]
            sage: s = Sequence([1,2,3], cr=True)
            sage: str(s)
            '[\n1,\n2,\n3\n]'
        """
        if self.__cr_str:
            return '[\n' + ',\n'.join([str(x) for x in self]) + '\n]'
        else:
            return list.__str__(self)

    def universe(self):
        """
        EXAMPLES::

            sage: Sequence([1,2/3,-2/5]).universe()
            Rational Field
            sage: Sequence([1,2/3,'-2/5']).universe()
            Category of objects
        """
        return self.__universe

    def _require_mutable(self):
        """
        EXAMPLES::

            sage: a = Sequence([1,2/3,'-2/5'])
            sage: a._require_mutable()
            sage: a.set_immutable()
            sage: a._require_mutable()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead.")

    def set_immutable(self):
        """
        Make this object immutable, so it can never again be changed.

        EXAMPLES::

            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.set_immutable()
            sage: v[3] = 7
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        self._is_immutable = True

    def is_immutable(self):
        """
        Return True if this object is immutable (can not be changed)
        and False if it is not.

        To make this object immutable use :meth:`set_immutable`.

        EXAMPLE::

            sage: v = Sequence([1,2,3,4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.is_immutable()
            False
            sage: v.set_immutable()
            sage: v.is_immutable()
            True
        """
        try:
            return self._is_immutable
        except AttributeError:
            return False

    def is_mutable(self):
        """
        EXAMPLES::

            sage: a = Sequence([1,2/3,-2/5])
            sage: a.is_mutable()
            True
            sage: a[0] = 100
            sage: type(a[0])
            <type 'sage.rings.rational.Rational'>
            sage: a.set_immutable()
            sage: a[0] = 50
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: a.is_mutable()
            False
        """
        try:
            return not self._is_immutable
        except AttributeError:
            return True


    def __copy__(self):
        """
        Return a copy of this sequence

        EXAMPLES::

            sage: s = seq(range(10))
            sage: t = copy(s)
            sage: t == s
            True
            sage: t.is_immutable == s.is_immutable
            True
            sage: t.is_mutable == s.is_mutable
            True

        """
        return Sequence(self,universe = self.__universe,
                        check = False,
                        immutable = self._is_immutable,
                        cr = self.__cr_str)

    def __getattr__(self, name):
        """
        Strictly for unpickling old 'Sequences'

        INPUT:

        - ``name`` - some string

        TESTS::

            sage: S = Sequence([])
            sage: del S._Sequence_generic__universe
            sage: S.universe()
            Traceback (most recent call last):
            ...
            AttributeError: 'Sequence_generic' object has no attribute '_Sequence_generic__universe'
            sage: S._Sequence__universe = 'foobar'
            sage: S.universe()
            'foobar'

        We test that :trac:`13998` is fixed::

            sage: S = Sequence([])
            sage: S.set_immutable()
            sage: del S._Sequence_generic__hash
            sage: hash(S)
            Traceback (most recent call last):
            ...
            AttributeError: 'Sequence_generic' object has no attribute '_Sequence_generic__hash'
            sage: S._Sequence__hash = 34
            sage: hash(S)
            34
        """
        if name == "_Sequence_generic__cr" and hasattr(self,"_Sequence__cr"):
            self.__cr = self._Sequence__cr
            return self.__cr
        elif name == "_Sequence_generic__cr_str" and hasattr(self,"_Sequence__cr_str"):
            self.__cr_str = self._Sequence__cr_str
            return self.__cr_str
        elif name == "_Sequence_generic__immutable" and hasattr(self,"_Sequence__immutable"):
            self.__immutable = self._Sequence__immutable
            return self.__immutable
        elif name == "_Sequence_generic__universe" and hasattr(self,"_Sequence__universe"):
            self.__universe = self._Sequence__universe
            return self.__universe
        elif name == "_Sequence_generic__hash" and hasattr(self,"_Sequence__hash"):
            self.__hash = self._Sequence__hash
            return self.__hash
        else:
            raise AttributeError("'Sequence_generic' object has no attribute '%s'"%name)
seq = Sequence

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.structure.sequence', 'Sequence', Sequence_generic)
