r"""
Factorizations

AUTHORS:
    -- William Stein (2006-01-22): added unit part as suggested by David Kohel.
    -- Nick Alexander (2008-01-19): added support for non-commuting factors.

EXAMPLES:

This example illustrates that the unit part is not discarded from
factorizations.

    sage: x = QQ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-5) * (x - 3) * (x - 2)
    sage: F.unit()
    -5
    sage: expand(F)
    -5*x^2 + 25*x - 30

The underlying list is the list of pairs $(p_i, e_i)$, where $p_i$
is prime and $e_i$ is an integer. The unit part is discarded by
the list.

    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]
    sage: len(F)
    2
    sage: F[1]
    (x - 2, 1)

In the ring $\Z[x]$, the integer $-5$ is not a unit, so the
factorization has three factors:

    sage: x = ZZ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-5) * (x - 3) * (x - 2)
    sage: F.unit()
    1
    sage: list(F)
    [(-5, 1), (x - 3, 1), (x - 2, 1)]
    sage: expand(F)
    -5*x^2 + 25*x - 30
    sage: len(F)
    3

On the other hand, -1 is a unit in $\Z$, so it is included in the unit.
    sage: x = ZZ['x'].0
    sage: f = -1*(x-2)*(x-3)
    sage: F = f.factor(); F
    (-1) * (x - 3) * (x - 2)
    sage: F.unit()
    -1
    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.misc.latex as latex
from sage.structure.sage_object import SageObject

class Factorization(SageObject, list):
    """
    EXAMPLES:
        sage: N = 2006
        sage: F = N.factor(); F
        2 * 17 * 59
        sage: F.unit()
        1
        sage: F = factor(-2006); F
        -1 * 2 * 17 * 59
        sage: F.unit()
        -1
        sage: loads(F.dumps()) == F
        True
    """
    def __init__(self, x, unit=None, cr=False, sort=True, simplify=True):
        if not isinstance(x, list):
            raise TypeError, "x must be a list"
        if isinstance(x, Factorization):
            if unit is None:
                unit = x.__unit
            else:
                unit = x.__unit * unit
        from sage.rings.integer import Integer
        for t in x:
            if not (isinstance(t, tuple) and len(t) == 2 and isinstance(t[1],(int, long, Integer))):
                raise TypeError, "x must be a list of tuples of length 2"
        list.__init__(self, x)
        if unit is None:
            if len(x) > 0:
                try:
                    unit = self.base_ring()(1)
                except AttributeError:
                    unit = Integer(1)
            else:
                unit = Integer(1)
        self.__unit = unit
        self.__cr = cr
        if self.is_commutative() and sort:
            self.sort()
        if simplify:
            self.simplify()

    def base_ring(self):
        """
        Return the parent structure of my factors.

        EXAMPLES:
            sage: F = factor(2006)
            sage: F.base_ring()
            Integer Ring

            sage: R.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: F = Factorization([(z, 2)], 3)
            sage: (F*F^-1).base_ring()
            Rational Field
        """
        if len(self) > 0:
            return self[0][0].parent()
        else:
            return self.unit().parent()

    def is_commutative(self):
        """
        Return True if my factors commute.

        EXAMPLES:
            sage: F = factor(2006)
            sage: F.is_commutative()
            True
            sage: K = QuadraticField(23, 'a')
            sage: F = K.factor_integer(13)
            sage: F.is_commutative()
            True
            sage: R.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: F = Factorization([(z, 2)], 3)
            sage: F.is_commutative()
            False
            sage: (F*F^-1).is_commutative()
            True
        """
        try:
            return self.base_ring().is_commutative()
        except:
            # This is not the mathematically correct default, but agrees with
            # history -- we've always assumed factored things commute
            return True

    def _set_cr(self, cr):
        self.__cr = cr

    def simplify(self):
        """
        Combine adjacent products that commute as much as possible.

        TESTS:
            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (y,2)], simplify=False); F
            x^3 * y^2 * y^2
            sage: F.simplify(); F
            x^3 * y^4
            sage: F * Factorization([(y, -2)], 2)
            (2) * x^3 * y^2
        """
        repeat = False
        simp = []
        import itertools
        for obj, agroup in itertools.groupby(list(self), lambda x: x[0]):
            xs = list(agroup)
            if len(xs) > 1:
                repeat = True
            n = sum([x[1] for x in xs])
            if n != 0:
                simp.append((obj, n))
        self[0:] = simp
        if repeat:
            self.simplify()

    def sort(self):
        """
        Try to order my factors: by .degree(), then .dimension(), then size.
        """
        if len(self) == 0:
            return
        try:
            a = self[0][0].degree()
            def _cmp(f,g):
                try:
                    return cmp((f[0].degree(),f[1],f[0]), (g[0].degree(),g[1],g[0]))
                except (AttributeError, NotImplementedError):
                    return cmp(f[0], g[0])
        except (AttributeError, NotImplementedError):
            try:
                a = self[0][0].dimension()
                def _cmp(f,g):
                    try:
                        return cmp(f[0].dimension()*f[1], g[0].dimension()*g[1])
                    except (AttributeError, NotImplementedError):
                        return cmp(f[0], g[0])
            except (AttributeError, NotImplementedError):
                def _cmp(f,g):
                    return cmp(f[0],g[0])
        list.sort(self, _cmp)

    def __reduce__(self):
        """
        Prepare me for pickling.
        """
        x = list(self)
        return Factorization, (x,)

    def unit(self):
        """
        Return the unit part of this factorization.

        EXAMPLES:
            sage: F = factor(-2006); F
            -1 * 2 * 17 * 59
            sage: F.unit()
            -1
        """
        return self.__unit

    def unit_part(self):
        r"""
        Same as \code{self.unit()}.
        """
        return self.__unit

    def _cr(self):
        try:
            return self.__cr
        except AttributeError:
            return False

    def _repr_(self):
        cr = self._cr()
        if len(self) == 0:
            return str(self.__unit)
        try:
            atomic = ((isinstance(self[0][0], (int, long)) or \
                       self.base_ring().is_atomic_repr()))
        except AttributeError:
            atomic = False
        s = ''
        mul =  ' * '
        if cr:
            mul += '\n'
        for i in range(len(self)):
            t = str(self[i][0])
            n = self[i][1]
            if (n>1 or len(self) > 1 or self.__unit != 1) and not atomic  and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)'%t
            if n != 1:
                t += '^%s'%n
            s += t
            if i < len(self)-1:
                s += mul
        if self.__unit != 1:
            if atomic:
                u = str(self.__unit)
            else:
                u = '(%s)'%self.__unit
            s =  u + mul + s
        return s

    def _latex_(self):
        if len(self) == 0:
            return latex.latex(self.__unit)
        try:
            atomic = ((isinstance(self[0][0], (int, long)) or \
                       self.base_ring().is_atomic_repr()))
        except AttributeError:
            atomic = False
        s = ''
        for i in range(len(self)):
            t = latex.latex(self[i][0])
            if not atomic and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)'%t
            n = self[i][1]
            if n != 1:
                t += '^{%s}'%latex.latex(n)
            s += t
            if i < len(self)-1:
                s += ' \\cdot '
        if self.__unit != 1:
            if atomic:
                u = latex.latex(self.__unit)
            else:
                u = '\\left(%s\\right)'%latex.latex(self.__unit)
            s =  u + ' \\cdot ' + s
        return s

    def __add__(self, other):
        """
        Return the sum of self and other.

        EXAMPLES:
            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
        """
        return self.value() + other


    def __sub__(self, other):
        """
        Return the sum of self and other.

        EXAMPLES:
            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
        """
        return self.value() - other

    def __mul__(self, other):
        """
        Return the product of two factorizations.

        EXAMPLES:
            sage: factor(-10) * factor(-16)
            2^5 * 5
            sage: factor(-10) * factor(16)
            -1 * 2^5 * 5

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F*F
            x^3 * y^2 * x^4 * y^2 * x
        """
        if not isinstance(other, Factorization):
            return self.value() * other
        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1.keys()).union(set(d2.keys())):
                s[a] = 0
                if d1.has_key(a):
                    s[a] += d1[a]
                if d2.has_key(a):
                    s[a] += d2[a]
            return Factorization(list(s.iteritems()), unit=self.unit()*other.unit())
        else:
            return Factorization(list(self) + list(other), unit=self.unit()*other.unit())

    def __pow__(self, n):
        """
        Return the n-th power of the factors in the factorization.

        EXAMPLES:
            sage: F = factor(2006); F
            2 * 17 * 59
            sage: F**2
            2^2 * 17^2 * 59^2

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F**2
            x^3 * y^2 * x^4 * y^2 * x
        """
        from sage.rings.arith import generic_power
        return generic_power(self, n, Factorization([]))

    def __invert__(self):
        """
        Return the formal inverse of the factors in the factorization.

        EXAMPLES:
            sage: F = factor(2006); F
            2 * 17 * 59
            sage: F^-1
            2^-1 * 17^-1 * 59^-1

            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)], 2); F
            (2) * x^3 * y^2 * x
            sage: F^-1
            (1/2) * x^-1 * y^-2 * x^-3
        """
        return Factorization([(p,-e) for p,e in reversed(self)],
            cr=self._cr(), unit=self.unit()**(-1))

    def value(self):
        """
        Return the product of the factors in the factorization, multiplied out.

        EXAMPLES:
            sage: F = factor(2006); F
            2 * 17 * 59
            sage: F.value()
            2006

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F.value()
            x^3*y^2*x
        """
        x = self.__unit
        for f, e in self:
            x *= (f**e)
        return x

    def expand(self):
        r"""
        Same as \code{self.value()}.
        """
        return self.value()

    def prod(self):
        r"""
        Same as \code{self.value()}.

        EXAMPLES:
            sage: F = factor(100)
            sage: F.prod()
            100
        """
        return self.value()


def Factorization_deduce_unit(x, mul):
    F = Factorization(x)
    z = F.value()
    u = mul/z
    F._Factorization__unit = u
    return F
