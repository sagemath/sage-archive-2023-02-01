r"""
Factorizations

The \code{Factorization} class derives from \code{list}, so it can
print nicely and be manipulated like a list of prime-exponent pairs or
easily turned into a list.  For example, we factor the integer $-45$:

    sage: F = factor(-45)

This returns an object of type \code{Factorization}:
    sage: type(F)
    <class 'sage.structure.factorization.Factorization'>

It prints in a nice factored form:
    sage: F
    -1 * 3^2 * 5

There is an underlying list representation, \emph{which ignores the
unit part} (!).
    sage: list(F)
    [(3, 2), (5, 1)]
    sage: isinstance(F, list)
    False

We can access the \code{Factorization} F itself as if it were a list:
    sage: F[0]
    (3, 2)
    sage: F[1]
    (5, 1)

To get at the unit part, use the \code{unit_part} function:
    sage: F.unit_part()
    -1

All factorizations are immutable.  Thus if you write a function that
returns a cached version of a factorization, you do not have to return
a copy.
    sage: F = factor(-12); F
    -1 * 2^2 * 3
    sage: F[0] = (5,4)
    Traceback (most recent call last):
    ...
    TypeError: 'Factorization' object does not support item assignment


EXAMPLES:

This more complicated example involving polynomials also illustrates
+that the unit part is not discarded from factorizations.

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

In the ring $\ZZ[x]$, the integer $-5$ is not a unit, so the
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

On the other hand, -1 is a unit in $\ZZ$, so it is included in the unit.
    sage: x = ZZ['x'].0
    sage: f = -1*(x-2)*(x-3)
    sage: F = f.factor(); F
    (-1) * (x - 3) * (x - 2)
    sage: F.unit()
    -1
    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]

Factorizations can involve fairly abstract mathematical objects:
    sage: F = ModularSymbols(11,4).factorization()
    sage: F
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field)
    sage: type(F)
    <class 'sage.structure.factorization.Factorization'>

TESTS:
    sage: F = factor(-20); F
    -1 * 2^2 * 5
    sage: G = loads(dumps(F)); G
    -1 * 2^2 * 5
    sage: G == F
    True
    sage: G is F
    False

AUTHORS:
    -- William Stein (2006-01-22): added unit part as suggested by David Kohel.
    -- William Stein (2008-01-17): wrote much of the documentation and fixed
                                   a couple of bugs.
    -- Nick Alexander (2008-01-19): added support for non-commuting factors.
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

class Factorization(SageObject):
    """
    A formal factorization of an object.

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
        sage: F = Factorization([(x,1/3)])
        Traceback (most recent call last):
        ...
        TypeError: powers of factors must be integers
    """
    def __init__(self, x, unit=None, cr=False, sort=True, simplify=True):
        """
        Create a \code{Factorization} object.

        INPUT:
            x    -- a list of (p, e) pairs with e an integer (or a TypeError
                    is raised).
            unit -- (default: 1) the unit part of the factorization
            cr   -- (default: False) if True, print the factorization with
                    carriage returns between factors
            sort -- (default: True) if True, sort the factors by calling
                    the sort function after creating the factorization.
                    See the documentation for self.sort for how this works.

        OUTPUT:
            a Factorization object

        EXAMPLES:
        We create a factorization with all the default options:
            sage: Factorization([(2,3), (5, 1)])
            2^3 * 5

        We create a factorization with a specified unit part:
            sage: Factorization([(2,3), (5, 1)], unit=-1)
            -1 * 2^3 * 5

        We try to create a factorization but with a string an exponent, which
        results in a TypeError:
            sage: Factorization([(2,3), (5, 'x')])
            Traceback (most recent call last):
            ...
            TypeError: powers of factors must be integers

        We create a factorization that puts newlines after each multiply sign when
        printing.  This is mainly useful when the primes are large.
            sage: Factorization([(2,3), (5, 2)], cr=True)
            2^3 *
            5^2

        Another factorization with newlines and nontrivial unit part (which appears
        on a line by itself):
            sage: Factorization([(2,3), (5, 2)], cr=True, unit=-2)
            -2 *
            2^3 *
            5^2

        A factorization, but where we do not sort the factors:
            sage: Factorization([(5,3), (2, 3)], sort=False)
            5^3 * 2^3

        By default factorizations are sorted by the prime base (for commutative bases):
            sage: Factorization([(2, 7), (5,2), (2, 5)])
            2^12 * 5^2
            sage: R.<a,b> = FreeAlgebra(QQ,2)
            sage: Factorization([(a,1),(b,1),(a,2)])
            a * b * a^2

        Autosorting (the default) swaps around the factors below:
            sage: F = Factorization([(ZZ^3, 2), (ZZ^2, 5)], cr=True); F
            (Ambient free module of rank 2 over the principal ideal domain Integer Ring)^5 *
            (Ambient free module of rank 3 over the principal ideal domain Integer Ring)^2
        """
        if not isinstance(x, list):
            raise TypeError, "x must be a list"
        if isinstance(x, Factorization):
            if unit is None:
                unit = x.__unit
            else:
                unit = x.__unit * unit
        from sage.rings.integer import Integer
        for i in xrange(len(x)):
            t=x[i]
            if not (isinstance(t, tuple) and len(t) == 2):
                raise TypeError, "x must be a list of tuples (p, e) of length 2 with e an integer"
            if not isinstance(t[1],(int, long, Integer)):
                try: # try coercing to an integer
                    x[i]= (t[0], Integer(t[1]))
                except TypeError:
                    raise TypeError, "powers of factors must be integers"

        self.__x = [ (t[0],int(t[1])) for t in x]
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

    def __getitem__(self, i):
        """
        Return i-th factor of self.

        EXAMPLES:
            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0]
            (3, 1)
            sage: a[1]
            (5, 2)
            sage: a[-1]
            (5, 2)
            sage: a[5]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.__x.__getitem__(i)

    def __setitem__(self, i, v):
        """
        Set the i-th factor of self.

        NOT ALLOWED -- Factorizations are immutable.

        EXAMPLES:
            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0] = (2,3)
            Traceback (most recent call last):
            ...
            TypeError: 'Factorization' object does not support item assignment
        """
        raise TypeError, "'Factorization' object does not support item assignment"
        #from sage.rings.integer import Integer
        #if len(v) != 2 or not isinstance(v[1],(int,long,Integer)):
        #    raise TypeError, "right hand side must be a pair (p,e) with e an integer."
        #return self.__x.__setitem__(i, v)

    def __len__(self):
        """
        Return the number of prime factors of self, not counting
        the unit part.

        EXAMPLES:
            sage: len(factor(15))
            2

        Note that the unit part is not included in the count.
            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: len(a)
            2
            sage: list(a)
            [(3, 1), (5, 2)]
            sage: len(list(a))
            2
        """
        return len(self.__x)

    def __cmp__(self, other):
        """
        Compare self and other.  This compares the underlying
        lists of self and other (ignoring the unit!)

        EXAMPLES:
        We compare two contrived formal factorizations:
            sage: a = Factorization([(2, 7), (5,2), (2, 5)])
            sage: b = Factorization([(2, 7), (5,10), (7, 3)])
            sage: a
            2^12 * 5^2
            sage: b
            2^7 * 5^10 * 7^3
            sage: a < b
            True
            sage: b < a
            False
            sage: a.expand()
            102400
            sage: b.expand()
            428750000000

        We compare factorizations of some polynomials:
            sage: x = polygen(QQ)
            sage: x^2 - 1 > x^2 - 4
            True
            sage: factor(x^2 - 1) > factor(x^2 - 4)
            True
        """
        if not isinstance(other, Factorization):
            return cmp(type(self), type(other))
        try:
            return cmp(self.expand(), other.expand())
        except:
            c = cmp(self.__unit, other.__unit)
            if c: return c
            return list.__cmp__(self, other)

    def __copy__(self):
        r"""
        Return a copy of self.

        This is of course not a deepcopy -- only references to the
        factors are returned, not copies of them.  Use
        \code{deepcopy(self)} if you need a deep copy of self.

        EXAMPLES:
        We create a factorization that has mutable primes:
            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        We make a copy of it:
            sage: G = copy(F); G
            ([1, 2])^5 * ([5, 6])^10
            sage: G is F
            False

        Note that if we change one of the mutable "primes" of F, this does
        change G.
            sage: F[1][0][0] = 'hello'
            sage: G
            ([1, 2])^5 * (['hello', 6])^10
        """
        # no need to sort, since the factorization is already sorted in whatever
        # order is desired.
        return Factorization(self.__x, unit=self.__unit, cr=self.__cr, sort=False, simplify=False)

    def __deepcopy__(self, memo):
        r"""
        Return a deep copy of self.

        This is of course not a deepcopy -- only references to the factors
        are returned, not copies of them.

        EXAMPLES:
        We make a factorization that has mutable entries:
            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        Now we make a copy of it and a deep copy.
            sage: K = copy(F)
            sage: G = deepcopy(F); G
            ([1, 2])^5 * ([5, 6])^10

        We change one of the mutable entries of F:
            sage: F[0][0][0] = 10

        This of course changes F:
            sage: F
            ([10, 2])^5 * ([5, 6])^10

        It also changes the copy K of F:
            sage: K
            ([10, 2])^5 * ([5, 6])^10

        It does \emph{not} change the deep copy G:
            sage: G
            ([1, 2])^5 * ([5, 6])^10
        """
        import copy
        return Factorization(copy.deepcopy(list(self), memo), cr=self.__cr, sort=False)

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
        """
        Change whether or not the factorization is printed with
        carriage returns after each factor.

        EXAMPLES:
            sage: x = polygen(QQ,'x')
            sage: F = factor(x^6 - 1); F
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
            sage: F._set_cr(True); F
            (x - 1) *
            (x + 1) *
            (x^2 - x + 1) *
            (x^2 + x + 1)
            sage: F._set_cr(False); F
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
        """
        self.__cr = bool(cr)

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
        self.__x[0:] = simp
        if repeat:
            self.simplify()

    def sort(self, _cmp=None):
        r"""
        Sort the factors in this factorization.

        INPUT:
            _cmp -- (default: None) comparison function

        OUTPUT:
            changes this factorization to be sorted

        If _cmp is None, we determine the comparison function as
        follows: If the prime in the first factor has a dimension
        method, then we sort based first on \emph{dimension} then on
        the exponent.  If there is no dimension method, we next
        attempt to sort based on a degree method, in which case, we
        sort based first on \emph{degree}, then exponent to break ties
        when two factors have the same degree, and if those match
        break ties based on the actual prime itself.  If there is no
        degree method, we sort based on dimension.

        EXAMPLES:
        We create a factored polynomial:
            sage: x = polygen(QQ,'x')
            sage: F = factor(x^3 + 1); F
            (x + 1) * (x^2 - x + 1)

        Then we sort it but using the negated version of the standard
        Python cmp function:
            sage: F.sort(_cmp = lambda x,y: -cmp(x,y))
            sage: F
            (x^2 - x + 1) * (x + 1)
        """
        if len(self) == 0:
            return
        if _cmp is None:
            try:
                a = self.__x[0][0].dimension()
                def _cmp(f,g):
                    """
                    This is used internally for comparing.  (indirect doctest)

                    EXAMPLES:
                        sage: factor(6)
                        2 * 3
                    """
                    try:
                        return cmp((f[0].dimension(), f[1]), (g[0].dimension(),g[1]))
                    except (AttributeError, NotImplementedError):
                        return cmp((f[0],f[1]), (g[0], g[1]))
            except (AttributeError, NotImplementedError):
                try:
                    a = self.__x[0][0].degree()
                    def _cmp(f,g):
                        """
                        This is used internally for comparing.  (indirect doctest)

                        EXAMPLES:
                            sage: factor(6)
                            2 * 3
                        """
                        try:
                            return cmp((f[0].degree(),f[1],f[0]), (g[0].degree(),g[1],g[0]))
                        except (AttributeError, NotImplementedError):
                            return cmp(f[0], g[0])
                except (AttributeError, NotImplementedError):
                    self.__x.sort()
                    return

        self.__x.sort(_cmp)

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

        EXAMPLES:
        We create a polynomial over the real double field and factor it:
            sage: x = polygen(RDF, 'x')
            sage: F = factor(-2*x^2 - 1); F
            (-2.0) * (1.0*x^2 + 0.5)

        Note that the unit part of the factorization is $-2.0$.
            sage: F.unit_part()
            -2.0
       """
        return self.__unit

    def _cr(self):
        """
        Return whether or not factorizations are printed with carriage returns
        between factors.

        EXAMPLES:
        Our first example involves factoring an integer:
            sage: F = factor(-93930); F
            -1 * 2 * 3 * 5 * 31 * 101
            sage: F._cr()
            False
            sage: F._set_cr(True)
            sage: F._cr()
            True

        This of course looks funny:
            sage: F
            -1 *
            2 *
            3 *
            5 *
            31 *
            101

        Next we factor a modular symbols space:
            sage: F = ModularSymbols(11).factor(); F
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...)
        """
        try:
            return self.__cr
        except AttributeError:
            self.__cr = False
            return False

    def _repr_(self):
        """
        Return the string representation of this factorization.

        EXAMPLES:
            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        Note that the default printing of a factorization can be overloaded
        using the rename method.
            sage: f.rename('factorization of -100')
            sage: f
            factorization of -100

        However _repr_ always prints normally.
            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        EXAMPLES:
           sage: x = polygen(QQ)
           sage: Factorization([(x-1,1), (x-2,2)])
            (x - 1) * (x - 2)^2
        """
        cr = self._cr()
        if len(self) == 0:
            return repr(self.__unit)
        try:
            atomic = ((isinstance(self.__x[0][0], (int, long)) or \
                       self.base_ring().is_atomic_repr()))
        except AttributeError:
            atomic = False
        s = ''
        mul =  ' * '
        if cr:
            mul += '\n'
        for i in range(len(self)):
            t = repr(self.__x[i][0])
            n = self.__x[i][1]
            if (n>1 or len(self) > 1 or self.__unit != 1) and not atomic  and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)'%t
            if n != 1:
                t += '^%s'%n
            s += t
            if i < len(self)-1:
                s += mul
        if self.__unit != 1:
            if atomic:
                u = repr(self.__unit)
            else:
                u = '(%s)'%self.__unit
            s =  u + mul + s
        return s

    def _latex_(self):
        r"""
        Return the \LaTeX{} representation of this factorization.

        EXAMPLES:
            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: latex(f)
            -1 \cdot 2^{2} \cdot 5^{2}
            sage: f._latex_()
            '-1 \\cdot 2^{2} \\cdot 5^{2}'
        """
        if len(self) == 0:
            return latex.latex(self.__unit)
        try:
            atomic = ((isinstance(self.__x[0][0], (int, long)) or \
                       self.base_ring().is_atomic_repr()))
        except AttributeError:
            atomic = False
        s = ''
        for i in range(len(self)):
            t = latex.latex(self.__x[i][0])
            if not atomic and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)'%t
            n = self.__x[i][1]
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
            sage: factor(100) + factor(19)
            119
        """
        if isinstance(other, Factorization):
            other = other.value()
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
        if isinstance(other, Factorization):
            other = other.value()
        return self.value() - other

    def __neg__(self):
        """
        Return negative of this factorization.

        EXAMPLES:
            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: -a
            3 * 5^2
            sage: (-a).unit()
            1
        """
        unit = -self.__unit
        return Factorization(list(self), unit, self.__cr, sort=False, simplify=False)

    def __rmul__(self, left):
        """
        Return the product left * self, where left is not a Factorization.

        EXAMPLES:
            sage: a = factor(15); a
            3 * 5
            sage: -2 * a
            -2 * 3 * 5
            sage: a * -2
            -2 * 3 * 5
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: f = Factorization([(x,2),(y,3)]); f
            x^2 * y^3
            sage: x * f
            x^3 * y^3
            sage: f * x
            x^2 * y^3 * x
        """
        return Factorization([(left, 1)]) * self

    def __mul__(self, other):
        r"""
        Return the product of two factorizations, which is obtained by
        combining together like factors.

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
            sage: -1 * F
            -1 * x^4 * y^2
        """
        if not isinstance(other, Factorization):
            return self * Factorization([(other, 1)])
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
        Return the $n$-th power of a factorization, which is got by
        combining together like factors.

        EXAMPLES:
            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f^3
            -1 * 2^6 * 5^6
            sage: f^4
            2^8 * 5^8

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
        from sage.groups.generic import power
        return power(self, n, Factorization([]))

    def __invert__(self):
        r"""
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
        Same as \code{self.value()}, so this returns the product of
        the factors, multiplied out.

            sage: x = polygen(QQ, 'x')
            sage: F = factor(-x^5 + 1); F
            (-1) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)
            sage: F.expand()
            -x^5 + 1
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


