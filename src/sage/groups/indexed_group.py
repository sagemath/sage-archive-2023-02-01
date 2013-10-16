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
from sage.monoids.indexed_monoid import IndexedFreeMonoid, IndexedFreeAbelianMonoid
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.infinity import infinity

class IndexedFreeGroup(IndexedFreeMonoid, Group):
    """
    An indexed free group.

    EXAMPLES::

        sage: G = IndexedFreeGroup(ZZ)
        sage: G
        Free group indexed by Integer Ring
        sage: G = IndexedFreeGroup('abcde')
        sage: G
        Free group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = IndexedFreeGroup(ZZ)
            sage: TestSuite(G).run()
            sage: G = IndexedFreeGroup('abc')
            sage: TestSuite(G).run()
        """
        Group.__init__(self)
        category = Groups().or_subcategory(category)
        IndexedFreeMonoid.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        TESTS::

            sage: IndexedFreeGroup(ZZ)
            Free group indexed by Integer Ring
        """
        return 'Free group indexed by {}'.format(self._indices)

    def __len__(self):
        """
        Return the length of ``self``.

        EXAMPLES::

            sage: F = IndexedFreeGroup(ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*c^-3*b^-2*a
            sage: elt.length()
            7
            sage: len(elt)
            7

            sage: F = IndexedFreeAbelianGroup(ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*c^-3*b^-2*a
            sage: elt.length()
            7
            sage: len(elt)
            7
        """
        return sum(abs(exp) for gen,exp in self._sorted_items())

    length = __len__

    def order(self):
        r"""
        Return the number of elements of ``self``, which is `\infty` unless
        this is the trivial group.

        EXAMPLES::

            sage: G = IndexedFreeGroup(ZZ)
            sage: G.order()
            +Infinity
            sage: G = IndexedFreeGroup('abc')
            sage: G.order()
            +Infinity
            sage: G = IndexedFreeGroup([])
            sage: G.order()
            1
        """
        if self.ngens() == 0:
            return Integer(1)
        return infinity

    def is_finite(self):
        """
        Return ``True`` if ``self`` is finite.

        EXAMPLES::

            sage: G = IndexedFreeGroup(ZZ)
            sage: G.is_finite()
            False
            sage: G = IndexedFreeGroup('abc')
            sage: G.is_finite()
            False
            sage: G = IndexedFreeGroup([])
            sage: G.is_finite()
            True
        """
        return self.ngens() == 0

    rank = IndexedFreeMonoid.ngens

    class Element(IndexedFreeMonoid.Element):
        def _mul_(self, y):
            """
            Multiply ``self`` by ``y``.

            EXAMPLES::

                sage: F = IndexedFreeGroup(ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: a*b^2*e*d
                F[0]*F[1]^2*F[4]*F[3]
                sage: (a*b^2*d^2) * (d^-4*b*e)
                F[0]*F[1]^2*F[3]^-2*F[1]*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            if len(self._monomial) == 0:
                return y
            if len(y._monomial) == 0:
                return self

            ret = list(self._monomial)
            rhs = list(y._monomial)
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

                sage: F = IndexedFreeGroup(ZZ)
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

                sage: F = IndexedFreeGroup(ZZ)
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

class IndexedFreeAbelianGroup(IndexedFreeAbelianMonoid, AbelianGroup):
    """
    An indexed free abelian group.

    EXAMPLES::

        sage: G = IndexedFreeAbelianGroup(ZZ)
        sage: G
        Free abelian group indexed by Integer Ring
        sage: G = IndexedFreeAbelianGroup('abcde')
        sage: G
        Free abelian group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = IndexedFreeAbelianGroup(ZZ)
            sage: TestSuite(G).run()
            sage: G = IndexedFreeAbelianGroup('abc')
            sage: TestSuite(G).run()
        """
        AbelianGroup.__init__(self)
        category = Groups().or_subcategory(category)
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        TESTS::

            sage: IndexedFreeAbelianGroup(ZZ)
            Free abelian group indexed by Integer Ring
        """
        return 'Free abelian group indexed by {}'.format(self._indices)

    def order(self):
        r"""
        Return the number of elements of ``self``, which is `\infty` unless
        this is the trivial group.

        EXAMPLES::

            sage: G = IndexedFreeAbelianGroup(ZZ)
            sage: G.order()
            +Infinity
            sage: G = IndexedFreeAbelianGroup('abc')
            sage: G.order()
            +Infinity
            sage: G = IndexedFreeAbelianGroup([])
            sage: G.order()
            1
        """
        if self.ngens() == 0:
            return Integer(1)
        return infinity

    def is_finite(self):
        """
        Return ``True`` if ``self`` is finite.

        EXAMPLES::

            sage: G = IndexedFreeAbelianGroup(ZZ)
            sage: G.is_finite()
            False
            sage: G = IndexedFreeAbelianGroup('abc')
            sage: G.is_finite()
            False
            sage: G = IndexedFreeAbelianGroup([])
            sage: G.is_finite()
            True
        """
        return self.ngens() == 0

    rank = IndexedFreeAbelianMonoid.ngens

    class Element(IndexedFreeAbelianMonoid.Element):
        def _mul_(self, y):
            """
            Multiply ``self`` by ``y``.

            EXAMPLES::

                sage: F = IndexedFreeAbelianGroup(ZZ)
                sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
                sage: a*b^2*e^-1*d
                F[0]*F[1]^2*F[3]*F[4]^-1
                sage: (a*b^2*d^2) * (d^-4*b^-2*e)
                F[0]*F[3]^-2*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            ret = copy(self._monomial)
            for k,v in y._monomial.iteritems():
                ret[k] = ret.get(k, 0) + v
                if ret[k] == 0:
                    del ret[k]
            return self.__class__(self.parent(), ret)

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: F = IndexedFreeAbelianGroup(ZZ)
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

                sage: F = IndexedFreeAbelianGroup(ZZ)
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

