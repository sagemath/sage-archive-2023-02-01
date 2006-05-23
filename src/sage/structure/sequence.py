r"""
Sequences

A mutable sequence of elements with a common guaranteed category,
which can be set immutable.

Sequence derives from list, so have all the functionality of lists and
can be used wherever lists are used.  When a sequence is created
without explicitly given the common universe of the elements, the
constructor coerces the first and second element to come
\emph{canonical} common parent, if possible, then the second and
third, etc..  If this is possible, it then coerces everything into the
canonical parent at the end.  (Note that canonical coercion is very
restrictive.)  The sequence then has a function \code{universe()}
which returns either the common canonical parent (if the coercion
succeeded), or the category of all objects (Objects()).  So if you
have a list $v$ and type

    sage: v = [1, 2/3, 5]
    sage: w = Sequence(v)
    sage: w.universe()
    Rational Field

then since \code{w.universe()} is $\Q$ then you're guaranteed that all
elements of $w$ are rationals:

    sage: v[0].parent()
    Integer Ring
    sage: w[0].parent()
    Rational Field

If you do assignment to $w$ this property of being rationals is guaranteed
to be preserved.

    sage: w[0] = 2
    sage: w[0].parent()
    Rational Field
    sage: w[0] = 'hi'
    Traceback (most recent call last):
    ...
    TypeError: unable to convert hi to a rational

However, if you do \code{w = Sequence(v)} and the resulting universe
is \code{Objects()}, the elements are not guaranteed to have any
special parent.  This is what should happen, e.g., with finite field
elements of different characteristics:

    sage: v = Sequence([GF(3)(1), GF(7)(1)])
    sage: v.universe()
    Category of objects

You can make a list immutable with \code{v.freeze()}.  Assignment is
never again allowed on an immutable list.

Creation of a sequence involves making a copy of the input list, and
substantial coercions.  It can be greatly sped up by explicitly
specifying the universe of the sequence:

    sage: v = Sequence(range(10000), universe=ZZ)
"""


##########################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################


import sage.misc.latex as latex
import sage.ext.sage_object
from mutability import Mutability

class Sequence(Mutability, sage.ext.sage_object.SageObject, list):
    """
    A mutable list of elements with a common guaranteed universe,
    which can be set immutable.

    A universe is either an object that supports coercion (e.g., a parent),
    or a category.

    INPUT:
        x -- a list or tuple instance
        universe -- (default: None) the universe of elements; if None determined from
                     first element; if list is empty, is category Objects() of all objects.
        check -- (default: True) whether to coerce the elements of x into the universe
        cr -- (default: False) if True, then print a carriage return after each comma
                               when printing this sequence.

    OUTPUT:
        a sequence

    EXAMPLES:
        sage: v = Sequence(range(10))
        sage: v.universe()
        <type 'int'>
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    You can also use seq for "Sequence", which is identical to using Sequence:
        sage: v = seq([1,2,1/1]); v
        [1, 2, 1]
        sage: v.universe()
        Rational Field
        sage: v.parent()
        Category of sequences in Rational Field
        sage: v.parent()([3,4/3])
        [3, 4/3]


    Note that assignment coerces if possible,
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
        TypeError: Unable to coerce rational (=2/3) to an Integer.

    Sequences can be used absolutely anywhere lists or tuples can be used:
        sage: isinstance(v, list)
        True

    Sequences are hashable (unlike Python lists), though the hashing
    is potentially slow, since it first involves conversion of the
    sequence to a tuple, and returning the hash of that.  The hash
    is cached, and is only recomputed if the sequence is changed
    (which has a small performance penalty for assignment).

        sage: hash(v)
        2083920238            # 32-bit
        -8049699692026128018  # 64-bit
        sage: v[0] = 10
        sage: hash(v)
        -377547984            # 32-bit
        -2271601447248391376  # 64-bit

    If you really know what you are doing, you can circumvent the type
    checking (for an efficiency gain):
        sage: list.__setitem__(v, int(1), 2/3)        # bad circumvention
        sage: v
        [10, 2/3, 2, 5, 4, 5, 6, 7, 8, 9]
        sage: list.__setitem__(v, int(1), int(2))     # not so bad circumvention

    You can make a sequence with a new universe from an old sequence.
        sage: w = Sequence(v, QQ)
        sage: w
        [10, 2, 2, 5, 4, 5, 6, 7, 8, 9]
        sage: w.universe()
        Rational Field
        sage: w[1] = 2/3
        sage: w
        [10, 2/3, 2, 5, 4, 5, 6, 7, 8, 9]

    Sequences themselves live in a category, the category of all sequences
    in the given universe.
        sage: w.category()
        Category of sequences in Rational Field

    This is also the parent of any sequence:
        sage: w.parent()
        Category of sequences in Rational Field

    The default universe for any sequence, if no compatible parent structure
    can be found, is the universe of all SAGE objects.

    This example illustrates how every element of a list is taken into account
    when constructing a sequence.
        sage: v = Sequence([1,7,6,GF(5)(3)]); v
        [1, 2, 1, 3]
        sage: v.universe()
        Finite Field of size 5
        sage: v.parent()
        Category of sequences in Finite Field of size 5
        sage: v.parent()([7,8,9])
        [2, 3, 4]
    """
    def __init__(self, x, universe=None, check=True,
                 immutable=False, cr=False):
        if not isinstance(x, (list, tuple)):
            x = list(x)
            #raise TypeError, "x must be a list or tuple"
        self.__hash = None
        self.__cr = cr
        if isinstance(x, Sequence):
            if universe is None or universe == x.__universe:
                list.__init__(self, x)
                self.__universe = x.__universe
                return
        if universe is None:
            if len(x) == 0:
                import sage.categories.all
                universe = sage.categories.all.Objects()
            else:
                import sage.ext.coerce
                x = list(x)   # make a copy, or we'd change the type of the elements of x, which would be bad.
                for i in range(len(x)-1):
                    try:
                        x[i], x[i+1] = sage.ext.coerce.canonical_coercion(x[i],x[i+1])
                    except TypeError:
                        import sage.categories.all
                        universe = sage.categories.all.Objects()
                        check = False  # no point
                        break
                if universe is None:   # no type errors raised.
                    universe = sage.ext.coerce.parent(x[len(x)-1])
                #universe = sage.ext.coerce.parent(x[0])
        self.__universe = universe
        if check:
            x = [universe(t) for t in x]
        list.__init__(self, x)
        Mutability.__init__(self, immutable)

    def __setitem__(self, n, x):
        self._require_mutable()
        y = self.__universe(x)
        list.__setitem__(self, n, y)
        self.__hash=None

    def __hash__(self):
        if self.__hash is None:
            self.__hash = hash(tuple(self))
        return self.__hash

    def _repr_(self):
        if self.__cr:
            return '[\n' + ',\n'.join([str(x) for x in self]) + '\n]'
        else:
            return list.__repr__(self)

    def category(self):
        import sage.categories.all
        return sage.categories.all.Sequences(self.universe())

    def parent(self):
        return self.category()

    def universe(self):
        return self.__universe


seq = Sequence
