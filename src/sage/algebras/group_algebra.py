r"""
Class for group algebras of arbitrary groups (over a general commutative base
ring).

NOTE:
    -- It seems to be impossible to make this fit nicely with Sage's coercion
    model. The problem is that (for example) if G is the additive group (ZZ,+),
    and R = ZZ[G] is its group ring, then the integer 2 can be coerced into R
    in two ways -- via G, or via the base ring -- and *the answers are
    different*. In practice we get around this by preventing elements of G
    coercing automatically into ZZ[G], which is a shame, but makes more sense
    than preventing elements of the base ring doing so.

AUTHOR:
    -- David Loeffler (2008-08-24): initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008 David Loeffler <d.loeffler.01@cantab.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
from sage.rings.all import IntegerRing
from sage.groups.group import Group
from sage.structure.formal_sum import FormalSums_generic, FormalSums, FormalSum
from sage.sets.set import Set


class GroupAlgebra(Algebra):

    def __init__(self, group, base_ring = IntegerRing()):
        r""" Create the given group algebra.
        INPUT:
            -- (Group) group: a generic group.
            -- (Ring) base_ring: a commutative ring.
        OUTPUT:
            -- a GroupAlgebra instance.

        EXAMPLES:
            sage: GroupAlgebra(GL(3, GF(7)))
            Group algebra of group "General Linear Group of degree 3 over Finite
            Field of size 7" over base ring Integer Ring
            sage: GroupAlgebra(1)
            Traceback (most recent call last):
            ...
            TypeError: "1" is not a group
        """
        if not base_ring.is_commutative():
            raise NotImplementedError, "Base ring must be commutative"

        if not isinstance(group, Group):
            raise TypeError, '"%s" is not a group' % group

        Algebra.__init__(self, base_ring)
        self._formal_sum_module = FormalSums(base_ring)
        self._group = group

    def group(self):
        r""" Return the group of this group algebra.
        EXAMPLES:
            sage: GroupAlgebra(GL(3, GF(11))).group()
            General Linear Group of degree 3 over Finite Field of size 11
            sage: GroupAlgebra(SymmetricGroup(10)).group()
            Symmetric group of order 10! as a permutation group
        """
        return self._group

    def is_commutative(self):
        r""" Return True if self is a commutative ring. True if and only if
        self.group() is abelian.

        EXAMPLES:
            sage: GroupAlgebra(SymmetricGroup(2)).is_commutative()
            True
            sage: GroupAlgebra(SymmetricGroup(3)).is_commutative()
            False
        """
        return self.group().is_abelian()

    def is_field(self):
        r""" Return True if self is a field. This is always false unless
        self.group() is trivial and self.base_ring() is a field.
        EXAMPLES:
            sage: GroupAlgebra(SymmetricGroup(2)).is_field()
            False
            sage: GroupAlgebra(SymmetricGroup(1)).is_field()
            False
            sage: GroupAlgebra(SymmetricGroup(1), QQ).is_field()
            True
        """
        if not self.base_ring().is_field():
            return False
        return (self.group().order() == 1)

    def is_finite(self):
        r""" Return True if self is finite, which is true if and only if
        self.group() and self.base_ring() are both finite.

        EXAMPLES:
            sage: GroupAlgebra(SymmetricGroup(2), IntegerModRing(10)).is_finite()
            True
            sage: GroupAlgebra(SymmetricGroup(2)).is_finite()
            False
            sage: GroupAlgebra(AbelianGroup(1), IntegerModRing(10)).is_finite()
            False
        """
        return (self.base_ring().is_finite() and self.group().is_finite())

    def is_exact(self):
        r""" Return True if elements of self have exact representations,
        which is true of self if and only if it is true of self.group()
        and self.base_ring().

        EXAMPLES:
            sage: GroupAlgebra(GL(3, GF(7))).is_exact()
            True
            sage: GroupAlgebra(GL(3, GF(7)), RR).is_exact()
            False
            sage: GroupAlgebra(GL(3, pAdicRing(7))).is_exact() # not implemented correctly (not my fault)!
            False
        """
        return self.group().is_exact() and self.base_ring().is_exact()

    def is_integral_domain(self):
        r""" Return True if self is an integral domain.

        This is false unless
        self.base_ring() is an integral domain, and even then it is false unless
        self.group() has no nontrivial elements of finite order. I don't know if
        this condition suffices, but it obviously does if the group is abelian and
        finitely generated.

        EXAMPLES:
            sage: GroupAlgebra(SymmetricGroup(2)).is_integral_domain()
            False
            sage: GroupAlgebra(SymmetricGroup(1)).is_integral_domain()
            True
            sage: GroupAlgebra(SymmetricGroup(1), IntegerModRing(4)).is_integral_domain()
            False
            sage: GroupAlgebra(AbelianGroup(1)).is_integral_domain()
            True
            sage: GroupAlgebra(AbelianGroup(2, [0,2])).is_integral_domain()
            False
            sage: GroupAlgebra(GL(2, ZZ)).is_integral_domain() # not implemented
            False
        """
        if not self.base_ring().is_integral_domain():
            return False
        if self.group().is_finite():
            if self.group().order() > 1:
                return False
            else:
                return True
        if self.group().is_abelian():
            invs = self.group().invariants()
            if Set(invs) != Set([0]):
                return False
            else:
                return True
        if not self.group().is_abelian():
            raise NotImplementedError

    # I haven't written is_noetherian(), because I don't know when group
    # algebras are noetherian, and I haven't written is_prime_field(), because
    # I don't know if that means "is canonically isomorphic to a prime field"
    # or "is identical to a prime field".

    def _coerce_impl(self, x):
        return self(self.base_ring().coerce(x))

    def _an_element_impl(self):
        """
        Return an element of self.

        EXAMPLE:
            sage: GroupAlgebra(SU(2, 13), QQ).an_element() # random; hideous formatting!
            -1/95*[       9 2*a + 12]
            [       0        3] - 4*[      9 9*a + 2]
            [3*a + 5       1]
        """
        try:
            return self(self._formal_sum_module([
                (self.base_ring().random_element(), self.group().random_element()),
                (self.base_ring().random_element(), self.group().random_element()),
                ]))
        except: # base ring or group might not implement .random_element()
            return self(self._formal_sum_module([ (self.base_ring().an_element(), self.group().an_element()) ]))

    def __call__(self, x, check=True):
        r"""
        Create an element of this group algebra.

        INPUT:
            -- x: either a FormalSum element consisting of elements of
            self.group(), an element of self.base_ring(), or an element
            of self.group().
            -- check (boolean): whether or not to check that the given elements
            really do lie in self.group().

        OUTPUT:
            -- a GroupAlgebraElement instance whose parent is self.

        EXAMPLES:
            sage: G = AbelianGroup(1)
            sage: f = G.gen()
            sage: ZG = GroupAlgebra(G)
            sage: ZG(f)
            f
            sage: ZG(1) == ZG(G(1))
            True
            sage: ZG(FormalSum([(1,f), (2, f**2)]))
            2*f^2 + f
            sage: G = GL(2,7)
            sage: OG = GroupAlgebra(G, ZZ[sqrt(5)])
            sage: OG(2)
            2*[1 0]
            [0 1]
            sage: OG(G(2)) # conversion is not the obvious one
            [2 0]
            [0 2]
            sage: OG(FormalSum([ (1, G(2)), (2, RR(0.77)) ]) )
            Traceback (most recent call last):
            ...
            TypeError: 0.770000000000000 is not an element of group General Linear Group of degree 2 over Finite Field of size 7
            sage: OG(FormalSum([ (1, G(2)), (2, RR(0.77)) ]), check=False)
            [2 0]
            [0 2] + 2*0.770000000000000
            sage: OG(OG.base_ring().gens()[1])
            sqrt5*[1 0]
            [0 1]
            """
        return GroupAlgebraElement(self, x, check)

    def __eq__(self, other):
        r""" Test for equality.
        EXAMPLES:
            sage: GroupAlgebra(AbelianGroup(1)) == GroupAlgebra(AbelianGroup(1))
            True
            sage: GroupAlgebra(AbelianGroup(1), QQ) == GroupAlgebra(AbelianGroup(1), ZZ)
            False
            sage: GroupAlgebra(AbelianGroup(2)) == GroupAlgebra(AbelianGroup(1))
            False
            """
        if not isinstance(other, GroupAlgebra):
            return False
        else:
            return self.base_ring() == other.base_ring() and self.group() == other.group()

    def _repr_(self):
        r""" String representation of self. See GroupAlgebra.__init__ for a
        doctest."""
        return "Group algebra of group \"%s\" over base ring %s" % (self.group(), self.base_ring())

    def element_class(self):
        r"""
        The class of elements of self, which is GroupAlgebraElement.

        EXAMPLES:
            sage: GroupAlgebra(SU(2, GF(4,'a'))).element_class()
            <class 'sage.algebras.group_algebra.GroupAlgebraElement'>
        """
        return GroupAlgebraElement


    def category(self):
        r"""
        The category to which self belongs, which is the category of group algebras over self.base_ring().

        EXAMPLES:
            sage: GroupAlgebra(SU(2, GF(4, 'a')), IntegerModRing(12)).category()
            Category of group algebras over Ring of integers modulo 12
        """
        from sage.categories.category_types import GroupAlgebras
        return GroupAlgebras(self.base_ring())



class GroupAlgebraElement(AlgebraElement):

    def __init__(self, parent, x, check):
        r""" Create an element of the parent group algebra. Not intended to be
        called by the user; see GroupAlgebra.__call__ for examples and
        doctests."""
        AlgebraElement.__init__(self, parent)

        if not hasattr(x, 'parent'):
            x = IntegerRing()(x) # occasionally coercion framework tries to pass a Python int

        if isinstance(x, FormalSum):
            if check:
                for c,d in x._data:
                    if d.parent() != self.parent().group():
                        raise TypeError, "%s is not an element of group %s" % (d, self.parent().group())
                self._fs = x
            else:
                self._fs = x

        elif self.base_ring().has_coerce_map_from(x.parent()):
            self._fs = self.parent()._formal_sum_module([ (x, self.parent().group()(1)) ])
        elif self.parent().group().has_coerce_map_from(x.parent()):
            self._fs = self.parent()._formal_sum_module([ (1, self.parent().group()(x)) ])
        else:
            raise TypeError, "Don't know how to create an element of %s from %s" % (self.parent(), x)

    def _repr_(self):
        return self._fs._repr_()

    def _add_(self, other):
        r"""
        Add self to other.

        EXAMPLE:
            sage: G = GL(3, GF(7))
            sage: ZG = GroupAlgebra(G)
            sage: g1 = G([0,0,2,2,5,0,6,6,2])
            sage: s = ZG(g1)
            sage: s + s
            2*[0  0  2]
            [2  5  0]
            [6  6  2]
"""
        fs_sum = self._fs + other._fs
        return self.parent()(fs_sum, check=False)

    def _mul_(self, right):
        r""" Calculate self*right, where both self and right are GroupAlgebraElements.

        EXAMPLE:
            sage: G = GL(3, GF(7))
            sage: ZG = GroupAlgebra(G)
            sage: a, b = G.random_element(), G.random_element()
            sage: za, zb = ZG(a), ZG(b)
            sage: za*ZG(2) # random
            2*[4,5,0]
            [0,5,1]
            [2,5,1]
            sage: za*2 == za*ZG(2)
            True
            sage: (ZG(1) + za)*(ZG(2) + zb) == ZG(FormalSum([ (2,G(1)), (2,a), (1, b), (1, a*b)]))
            True
            sage: za*za == za^2
            True
        """
        d1 = self._fs._data
        d2 = right._fs._data
        new = []
        for (a1, g1) in d1:
            for a2,g2 in d2:
                if self.parent().group().is_multiplicative():
                    new.append( (a1*a2, g1*g2) )
                else:
                    new.append( (a1*a2, g1 + g2) )
        return self.parent()( self.parent()._formal_sum_module(new), check=False)

    def __eq__(self, other):
        r""" Test if self is equal to other.

        EXAMPLES:
            sage: G = AbelianGroup(1,[4])
            sage: a = GroupAlgebra(G)(1)
            sage: b = GroupAlgebra(G)(2)
            sage: a + a == b
            True
            sage: a == b
            False
            sage: a == GroupAlgebra(AbelianGroup(1, [5]))(1)
            False
        """
        if isinstance(other, GroupAlgebraElement) and self.parent() == other.parent():
            return self._fs == other._fs
        else:
            return False
