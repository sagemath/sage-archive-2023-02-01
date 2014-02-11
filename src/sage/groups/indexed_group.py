"""
Indexed Free Groups

Free groups and free abelian groups implemented using an indexed set of
generators.

AUTHORS:

- Travis Scrimshaw (2013-10-16): Initial version
"""

##############################################################################
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from copy import copy
from sage.categories.groups import Groups
from sage.groups.group import Group, AbelianGroup
from sage.monoids.indexed_monoid import IndexedMonoidElement, IndexedFreeMonoid, \
        IndexedFreeAbelianMonoid
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.infinity import infinity

class IndexedFreeGroup(IndexedFreeMonoid, Group):
    """
    An indexed free group.

    EXAMPLES::

        sage: G = FreeGroup(index_set=ZZ)
        sage: G
        Free group indexed by Integer Ring
        sage: G = FreeGroup(index_set='abcde')
        sage: G
        Free group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = FreeGroup(index_set=ZZ)
            sage: TestSuite(G).run()
            sage: G = FreeGroup(index_set='abc')
            sage: TestSuite(G).run()
        """
        Group.__init__(self)
        category = Groups().or_subcategory(category)
        IndexedFreeMonoid.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        TESTS::

            sage: FreeGroup(index_set=ZZ)
            Free group indexed by Integer Ring
        """
        return 'Free group indexed by {}'.format(self._indices)

    def order(self):
        r"""
        Return the number of elements of ``self``, which is `\infty` unless
        this is the trivial group.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ)
            sage: G.order()
            +Infinity
            sage: G = FreeGroup(index_set='abc')
            sage: G.order()
            +Infinity
            sage: G = FreeGroup(index_set=[])
            sage: G.order()
            1
        """
        if self.is_finite():
            return Integer(1)
        return infinity

    def is_finite(self):
        """
        Return ``True`` if ``self`` is finite.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ)
            sage: G.is_finite()
            False
            sage: G = FreeGroup(index_set='abc')
            sage: G.is_finite()
            False
            sage: G = FreeGroup(index_set=[])
            sage: G.is_finite()
            True
        """
        return self.rank() == 0

    def rank(self):
        """
        Return the rank of ``self``, which is the number of
        generators of ``self``.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ)
            sage: G.rank()
            +Infinity
            sage: G = FreeGroup(index_set='abc')
            sage: G.rank()
            3
            sage: G = FreeGroup(index_set=[])
            sage: G.rank()
            0
        """
        return self.gens().cardinality()

    class Element(IndexedFreeMonoid.Element):
        def __lt__(self, y):
            """
            Check less than.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: a < b
                True
                sage: a^-1*b < b^-1*a
                True
                sage: a*b < a*a^-1
                False
                sage: a^-1*a < a^2
                True
                sage: a^2*b < a*b^-1*a*b
                True
            """
            if not isinstance(y, IndexedMonoidElement):
                return False
            return self.to_word_list() < y.to_word_list()

        def __len__(self):
            """
            Return the length of ``self``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: elt = a*c^-3*b^-2*a
                sage: elt.length()
                7
                sage: len(elt)
                7

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: elt = a*c^-3*b^-2*a
                sage: elt.length()
                7
                sage: len(elt)
                7
            """
            return sum(abs(exp) for gen,exp in self._sorted_items())

        length = __len__

        def _mul_(self, other):
            """
            Multiply ``self`` by ``other``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: a*b^2*e*d
                F[0]*F[1]^2*F[4]*F[3]
                sage: (a*b^2*d^2) * (d^-4*b*e)
                F[0]*F[1]^2*F[3]^-2*F[1]*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            if len(self._monomial) == 0:
                return other
            if len(other._monomial) == 0:
                return self

            ret = list(self._monomial)
            rhs = list(other._monomial)
            while len(ret) > 0 and len(rhs) > 0 and ret[-1][0] == rhs[0][0]:
                rhs[0] = (rhs[0][0], rhs[0][1] + ret.pop()[1])
                if rhs[0][1] == 0:
                    rhs.pop(0)
            ret += rhs
            return self.__class__(self.parent(), tuple(ret))

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; ~x
                F[3]^-1*F[4]*F[1]^-2*F[0]^-1
                sage: x * ~x
                1
            """
            return self.__class__(self.parent(), tuple((x[0], -x[1]) for x in reversed(self._monomial)))

        def __pow__(self, n):
            """
            Raise ``self`` to the power of ``n``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: x = a*b^2*e*a^-1; x
                F[0]*F[1]^2*F[4]*F[0]^-1
                sage: x^3
                F[0]*F[1]^2*F[4]*F[1]^2*F[4]*F[1]^2*F[4]*F[0]^-1
                sage: x^0
                1
                sage: x^-3
                F[0]*F[4]^-1*F[1]^-2*F[4]^-1*F[1]^-2*F[4]^-1*F[1]^-2*F[0]^-1
            """
            if not isinstance(n, (int, long, Integer)):
                raise TypeError("Argument n (= {}) must be an integer".format(n))
            if n == 0:
                return self.parent().one()
            if n == 1:
                return self
            if n == -1:
                return ~self
            if len(self._monomial) == 1:
                gen,exp = self._monomial[0]
                return self.__class__(self.parent(), ((gen, exp*n),))
            if n < 0:
                self = ~self
                n = -n
            ret = self
            for i in range(n-1):
                ret *= self
            return ret

        def to_word_list(self):
            """
            Return ``self`` as a word represented as a list whose entries
            are the pairs ``(i, s)`` where ``i`` is the index and ``s`` is
            the sign.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: x = a*b^2*e*a^-1
                sage: x.to_word_list()
                [(0, 1), (1, 1), (1, 1), (4, 1), (0, -1)]
            """
            sign = lambda x: 1 if x > 0 else -1 # It is never 0
            return [ (k, sign(e)) for k,e in self._sorted_items() for dummy in range(abs(e))]

class IndexedFreeAbelianGroup(IndexedFreeAbelianMonoid, AbelianGroup):
    """
    An indexed free abelian group.

    EXAMPLES::

        sage: G = FreeGroup(index_set=ZZ, abelian=True)
        sage: G
        Free abelian group indexed by Integer Ring
        sage: G = FreeGroup(index_set='abcde', abelian=True)
        sage: G
        Free abelian group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = FreeGroup(index_set=ZZ, abelian=True)
            sage: TestSuite(G).run()
            sage: G = FreeGroup(index_set='abc', abelian=True)
            sage: TestSuite(G).run()
        """
        AbelianGroup.__init__(self)
        category = Groups().or_subcategory(category)
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        TESTS::

            sage: FreeGroup(index_set=ZZ, abelian=True)
            Free abelian group indexed by Integer Ring
        """
        return 'Free abelian group indexed by {}'.format(self._indices)

    def order(self):
        r"""
        Return the number of elements of ``self``, which is `\infty` unless
        this is the trivial group.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ, abelian=True)
            sage: G.order()
            +Infinity
            sage: G = FreeGroup(index_set='abc', abelian=True)
            sage: G.order()
            +Infinity
            sage: G = FreeGroup(index_set=[], abelian=True)
            sage: G.order()
            1
        """
        if self.is_finite():
            return Integer(1)
        return infinity

    def is_finite(self):
        """
        Return ``True`` if ``self`` is finite.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ, abelian=True)
            sage: G.is_finite()
            False
            sage: G = FreeGroup(index_set='abc', abelian=True)
            sage: G.is_finite()
            False
            sage: G = FreeGroup(index_set=[], abelian=True)
            sage: G.is_finite()
            True
        """
        return self.rank() == 0

    def rank(self):
        """
        Return the rank of ``self``, which is the number of
        generators of ``self``.

        EXAMPLES::

            sage: G = FreeGroup(index_set=ZZ, abelian=True)
            sage: G.rank()
            +Infinity
            sage: G = FreeGroup(index_set='abc', abelian=True)
            sage: G.rank()
            3
            sage: G = FreeGroup(index_set=[], abelian=True)
            sage: G.rank()
            0
        """
        return self.gens().cardinality()

    class Element(IndexedFreeAbelianMonoid.Element, IndexedFreeGroup.Element):
        def _mul_(self, other):
            """
            Multiply ``self`` by ``other``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ, abelian=True)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: a*b^2*e^-1*d
                F[0]*F[1]^2*F[3]*F[4]^-1
                sage: (a*b^2*d^2) * (d^-4*b^-2*e)
                F[0]*F[3]^-2*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            ret = copy(self._monomial)
            for k,v in other._monomial.iteritems():
                ret[k] = ret.get(k, 0) + v
                if ret[k] == 0:
                    del ret[k]
            return self.__class__(self.parent(), ret)

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ, abelian=True)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; ~x
                F[0]^-1*F[1]^-2*F[3]^-1*F[4]
                sage: x * ~x
                1
            """
            return self.__pow__(-1)

        def __pow__(self, n):
            """
            Raise ``self`` to the power of ``n``.

            EXAMPLES::

                sage: F = FreeGroup(index_set=ZZ, abelian=True)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; x
                F[0]*F[1]^2*F[3]*F[4]^-1
                sage: x^3
                F[0]^3*F[1]^6*F[3]^3*F[4]^-3
                sage: x^0
                1
                sage: x^-3
                F[0]^-3*F[1]^-6*F[3]^-3*F[4]^3
            """
            if not isinstance(n, (int, long, Integer)):
                raise TypeError("Argument n (= {}) must be an integer".format(n))
            if n == 1:
                return self
            if n == 0:
                return self.parent().one()
            return self.__class__(self.parent(), {k:v*n for k,v in self._monomial.iteritems()})

