"""
Formal sums

AUTHORS:

 - David Harvey (2006-09-20): changed FormalSum not to derive from
   "list" anymore, because that breaks new Element interface

 - Nick Alexander (2006-12-06): added test cases.

 - William Stein (2006, 2009): wrote the first version in 2006, documented it in 2009.

FUNCTIONS:
    - ``FormalSums(ring)`` -- create the module of formal finite sums with
                          coefficients in the given ring.

    - ``FormalSum(list of pairs (coeff, number))`` -- create a formal sum

EXAMPLES::

    sage: A = FormalSum([(1, 2/3)]); A
    2/3
    sage: B = FormalSum([(3, 1/5)]); B
    3*1/5
    sage: -B
    -3*1/5
    sage: A + B
    3*1/5 + 2/3
    sage: A - B
    -3*1/5 + 2/3
    sage: B*3
    9*1/5
    sage: 2*A
    2*2/3
    sage: list(2*A + A)
    [(3, 2/3)]

TESTS::

    sage: R = FormalSums(QQ)
    sage: loads(dumps(R)) == R
    True
    sage: a = R(2/3) + R(-5/7); a
    -5/7 + 2/3
    sage: loads(dumps(a)) == a
    True
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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

import sage.misc.misc
import element
import sage.misc.latex

from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.parent_base import ParentWithBase

class FormalSums_generic(Module):
    def __init__(self, base=ZZ):
        """
        EXAMPLES::

            sage: FormalSums(ZZ)
            Abelian Group of all Formal Finite Sums over Integer Ring
            sage: FormalSums(GF(7))
            Abelian Group of all Formal Finite Sums over Finite Field of size 7
        """
        ParentWithBase.__init__(self, base)

    def _repr_(self):
        """
        EXAMPLES::
            sage: FormalSums(GF(7))
            Abelian Group of all Formal Finite Sums over Finite Field of size 7
            sage: FormalSums(GF(7))._repr_()
            'Abelian Group of all Formal Finite Sums over Finite Field of size 7'
        """
        return "Abelian Group of all Formal Finite Sums over %s"%self.base_ring()

    def __call__(self, x, check=True, reduce=True):
        """
        Make a formal sum in self from the element x.

        INPUT:
            - ``x`` -- formal sum, list or number
            - ``check`` -- bool (default: True)
            - ``reduce`` -- bool (default: True); whether to combine terms

        EXAMPLES::
            sage: P = FormalSum([(1,2/3)]).parent()
            sage: P([(1,2/3), (5,-2/9)])
            5*-2/9 + 2/3
        """
        if isinstance(x, FormalSum):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return FormalSum(x._data, check=False, reduce=False, parent=self)
            else:
                x = x._data
        if isinstance(x, list):
            return FormalSum(x, check=check,reduce=reduce,parent=self)
        if x == 0:
            return FormalSum([], check=False, reduce=False, parent=self)
        else:
            return FormalSum([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def _coerce_impl(self, x):
        """
        EXAMPLES:

            sage: P = FormalSum([(1,2/3)]).parent()
            sage: P._coerce_impl(2)
            2
        """
        if x == 0:
            return FormalSum([], check=False, reduce=False, parent=self)
        elif isinstance(x, FormalSum) and x.parent().has_coerce_map_from(self.base_ring()):
            return self(x)
        else:
            return FormalSum([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def base_extend(self, R):
        """
        EXAMPLES::

            sage: FormalSums(ZZ).base_extend(GF(7))
            Abelian Group of all Formal Finite Sums over Finite Field of size 7
        """
        if self.base_ring().has_coerce_map_from(R):
            return self
        elif R.has_coerce_map_from(self.base_ring()):
            return FormalSums(R)

    def __cmp__(left, right):
        """
        EXAMPLES::

            sage: FormalSums(ZZ) < FormalSums(QQ)
            True
            sage: FormalSums(QQ) == FormalSums(QQ)
            True
            sage: FormalSums(QQ) > FormalSums(ZZ)
            True
        """
        c = cmp(type(left), type(right))
        if c: return c
        return cmp(left.base_ring(), right.base_ring())

    def get_action_impl(self, other, op, self_is_left):
        """
        EXAMPLES::

            sage: A = FormalSums(RR).get_action(RR); A     # indirect doctest
            Right scalar multiplication by Real Field with 53 bits of precision on Abelian Group of all Formal Finite Sums over Real Field with 53 bits of precision

            sage: A = FormalSums(ZZ).get_action(QQ); A
            Right scalar multiplication by Rational Field on Abelian Group of all Formal Finite Sums over Rational Field
            with precomposition on left by Call morphism:
              From: Abelian Group of all Formal Finite Sums over Integer Ring
              To:   Abelian Group of all Formal Finite Sums over Rational Field
            sage: A = FormalSums(QQ).get_action(ZZ); A
            Right scalar multiplication by Integer Ring on Abelian Group of all Formal Finite Sums over Rational Field
        """
        import operator
        from sage.structure.coerce import LeftModuleAction, RightModuleAction
        from sage.categories.action import PrecomposedAction
        if op is operator.mul and isinstance(other, Parent):
            extended = self.base_extend(other)
            if self_is_left:
                action = RightModuleAction(other, extended)
                if extended is not self:
                    action = PrecomposedAction(action, extended.coerce_map_from(self), None)
            else:
                action = LeftModuleAction(other, extended)
                if extended is not self:
                    action = PrecomposedAction(action, None, extended.coerce_map_from(self))
            return action

    def _an_element_impl(self):
        """
        EXAMPLES::

            sage: FormalSums(ZZ)._an_element_impl()
            1
            sage: FormalSums(QQ)._an_element_impl()
            1/2*1
            sage: QQ.an_element()
            1/2
        """
        return FormalSum([(self.base_ring().an_element(), 1)], check=False, reduce=False, parent=self)

import weakref
cache = {}
def FormalSums(R=ZZ):
    """
    Return the R-module of finite formal sums with coefficients in R.

    INPUT:
        R -- a ring (default: ZZ)

    EXAMPLES::

        sage: FormalSums()
        Abelian Group of all Formal Finite Sums over Integer Ring
        sage: FormalSums(ZZ[sqrt(2)])
        Abelian Group of all Formal Finite Sums over Order in Number Field in sqrt2 with defining polynomial x^2 - 2
        sage: FormalSums(GF(9,'a'))
        Abelian Group of all Formal Finite Sums over Finite Field in a of size 3^2
    """
    try:
        F = cache[R]()
        if not F is None:
            return F
    except KeyError:
        pass
    F = FormalSums_generic(R)
    cache[R] = weakref.ref(F)
    return F


formal_sums = FormalSums_generic()

class FormalSum(ModuleElement):
    """
    A formal sum over a ring.
    """
    def __init__(self, x, parent=formal_sums, check=True, reduce=True):
        """
        INPUT:
            - ``x`` -- object
            - ``parent`` -- FormalSums(R) module (default: FormalSums(ZZ))
            - ``check`` -- bool (default: True) if False, might not coerce
                           coefficients into base ring, which can speed
                           up constructing a formal sum.
            - ``reduce`` -- reduce (default: True) if False, do not
                            combine common terms

        EXAMPLES::

            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)])
            4*2/3 - 5*7
            sage: a = FormalSum([(1,2/3), (3,2/3), (-5, 7)], reduce=False); a
            2/3 + 3*2/3 - 5*7
            sage: a.reduce()
            sage: a
            4*2/3 - 5*7
            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)], parent=FormalSums(GF(5)))
            4*2/3

        Notice below that the coefficient 5 doesn't get reduced modulo 5::

            sage: FormalSum([(1,2/3), (3,2/3), (-5, 7)], parent=FormalSums(GF(5)), check=False)
            4*2/3 - 5*7
        """
        if x == 0:
            x = []
        if check:
            k = parent.base_ring()
            try:
                x = [(k(t[0]), t[1]) for t in x]
            except (IndexError, KeyError), msg:
                raise TypeError, "%s\nInvalid formal sum"%msg
        self._data = x
        if parent is None:
            parent = formal_sums
        ModuleElement.__init__(self, parent)
        if reduce:
            self.reduce()

    def __iter__(self):
        """
        EXAMPLES::       indirect doctest

            sage: for z in FormalSum([(1,2), (5, 'a'), (-3, 7)]): print z
            (5, 'a')
            (1, 2)
            (-3, 7)
        """
        return iter(self._data)

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: v = FormalSum([(1,2), (5, 'a'), (-3, 7)]); v
            5*a + 2 - 3*7
            sage: v[0]         # indirect doctest
            (5, 'a')
            sage: v[1]
            (1, 2)
            sage: v[2]
            (-3, 7)
            sage: list(v)
            [(5, 'a'), (1, 2), (-3, 7)]
        """
        return self._data[n]

    def __len__(self):
        """
        EXAMPLES::

            sage: v = FormalSum([(1,2), (5, 'a'), (-3, 7)]); v
            5*a + 2 - 3*7
            sage: len(v)            # indirect test
            3
        """
        return len(self._data)

    def _repr_(self):
        """
        EXAMPLES::

            sage: a = FormalSum([(1,2/3), (-3,4/5), (7,Mod(2,3))])
            sage: a   # random, since comparing Mod(2,3) and rationals ill-defined
            sage: a._repr_()    # random
            '2/3 - 3*4/5 + 7*2'
        """
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.misc.repr_lincomb(symbols, coeffs)

    def _latex_(self):
        """
        EXAMPLES::

            sage: latex(FormalSum([(1,2), (5, 8/9), (-3, 7)]))  # indirect doctest
            5\cdot \frac{8}{9} + 2 - 3\cdot 7
        """
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.latex.repr_lincomb(symbols, coeffs)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: a = FormalSum([(1,3),(2,5)]); a
            3 + 2*5
            sage: b = FormalSum([(1,3),(2,7)]); b
            3 + 2*7
            sage: a < b                                      # indirect doctest
            True
            sage: b < a
            False
            sage: a == 1
            False
            sage: b == 0
            False
        """
        return cmp(self._data, right._data)

    def _neg_(self):
        """
        EXAMPLES::

            sage: -FormalSum([(1,3),(2,5)])                # indirect doctest
            -3 - 2*5
        """
        return self.__class__([(-c, s) for (c, s) in self._data], check=False, parent=self.parent())

    def _add_(self, other):
        """
        EXAMPLES::

            sage: FormalSum([(1,3/7),(2,5/8)]) + FormalSum([(1,3/7),(-2,5)])  # indirect doctest
            2*3/7 + 2*5/8 - 2*5
        """
        return self.__class__(self._data + other._data, check=False, parent=self.parent())

    def _lmul_(self, s):
        """
        EXAMPLES::

            sage: -3*FormalSum([(1,3/7),(-2,5)])      # indirect doctest
            -3*3/7 + 6*5
        """
        return self.__class__([(c*s, x) for (c, x) in self], check=False, parent=self.parent())

    def _rmul_(self, s):
        """
        EXAMPLES::

            sage: FormalSum([(1,3/7),(-2,5)])*(-3)      # indirect doctest
            -3*3/7 + 6*5
        """
        return self.__class__([(s*c, x) for (c, x) in self], check=False, parent=self.parent())

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: bool(FormalSum([(1,3/7),(-2,5)]))
            True
            sage: bool(FormalSums(QQ)(0))
            False
            sage: bool(FormalSums(QQ)(1))
            True
        """
        if len(self._data) == 0:
            return False
        for c, _ in self._data:
            if not c.is_zero():
                return True
        return False

    def reduce(self):
        """
        EXAMPLES::

            sage: a = FormalSum([(-2,3), (2,3)], reduce=False); a
            -2*3 + 2*3
            sage: a.reduce()
            sage: a
            0
        """
        if len(self) == 0:
            return
        v = [(x,c) for c, x in self if not c.is_zero()]
        if len(v) == 0:
            self._data = v
            return
        v.sort()
        w = []
        last = v[0][0]
        coeff = v[0][1]
        for x, c in v[1:]:
            if x == last:
                coeff += c
            else:
                if coeff != 0:
                    w.append((coeff, last))
                last = x
                coeff = c
        if not coeff.is_zero():
            w.append((coeff,last))
        self._data = w
