"""
Manin symbols
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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


from sage.misc.search import search

import sage.matrix.all
import sage.modular.congroup as congroup
import sage.modular.cusps as cusps
import sage.modular.modsym.p1list as p1list
import sage.modular.modsym.g1list as g1list
import sage.modular.modsym.ghlist as ghlist
import sage.modular.modsym.modular_symbols
import sage.rings.arith as arith
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.all as rings

from sage.structure.sage_object import SageObject

R = polynomial_ring.PolynomialRing(rings.QQ, 'X')

X = R([0,1])

def is_ManinSymbol(x):
    return isinstance(x, ManinSymbol)

class ManinSymbolList(SageObject):
    """
    All Manin symbols for a given group, weight, and character.
    """
    def __init__(self, weight, list):
        self._weight = weight
        self._list = list

    def __cmp__(self, right):
        if not isinstance(right, ManinSymbolList):
            return -1
        return cmp([self._weight, self._list], [right._weight, right._list])

    def __getitem__(self, n):
        return self._list[n]

    def __len__(self):
        return len(self._list)

    def apply(self, j, X):
        raise NotImplementedError

    def apply_S(self, j):
        raise NotImplementedError

    def apply_I(self, j):
        raise NotImplementedError

    def apply_T(self, j):
        raise NotImplementedError

    def apply_TT(self, j):
        raise NotImplementedError

    def index(self, x):
        """
        Return the index into the list of Manin symbols of x, where x
        is a 3-tuple of ints.  If x is not in the list, then this
        function returns -1.

        INPUT:
            x -- 3-tuple of ints. Something equivalent to an element
                 of Manin symbols list, which need not be normalized.
        OUTPUT:
            int -- the index of the Manin symbol equivalent to (i,u,v).
        """
        t, i = search(self._list, x)
        if t: return i
        x = self.normalize(x)
        t, i = search(self._list, x)
        if t: return i
        return -1

    def manin_symbol_list(self):
        try:
            return self.__manin_symbol_list
        except AttributeError:
            self.__manin_symbol_list = [self.manin_symbol(i) \
                                        for i in xrange(len(self))]
        return self.__manin_symbol_list

    def manin_symbol(self, i):
        return ManinSymbol(self, self._list[int(i)])

    def normalize(self, x):
        raise NotImplementedError

    def weight(self):
        return self._weight

class ManinSymbolList_group(ManinSymbolList):
    """
    Base class for Manin symbol lists for a given group.

    ManinSymbolList_group(level, weight, syms):

    INPUT:
        level -- integer level
        weight -- integer weight >= 2
        syms -- something with a normalize and list method, e.g., P1List.
    """
    def __init__(self, level, weight, syms):
        self.__level = level
        self.__syms = syms  # syms is anything with a normalize and list method.

        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight - 2 + 1) \
                            for u, v in syms.list()]
        ManinSymbolList.__init__(self, weight, L)

    def apply_S(self, j):
        """
        """
        i, u, v = self._list[j]
        k = self.index((self._weight-2-i, v, -u))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def apply_I(self, j):
        """
        """
        i, u, v = self._list[j]
        k = self.index((i, -u, v))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def apply_T(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__syms.normalize(v,-u-v)
        if (k-2) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        for j in range(k-2-i +1):
            m = self.index((j, u, v))
            z.append((m,s*arith.binomial(k-2-i,j)))
            s *= -1
        return z

    def apply_TT(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__syms.normalize(-u-v,u)
        if (k-2-i) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        for j in range(i+1):
            m = self.index((k-2-i+j, u, v))
            z.append((m,s*arith.binomial(i,j)))
            s *= -1
        return z

    def apply(self, j, m):
        """
        Apply the matrix m=[a,b,c,d] to the j-th Manin symbol.

        INPUT:
            j -- integer
            m = [a, b, c, d] a list of 4 integers, which defines a 2x2 matrix.

        OUTPUT:
            list -- a list of pairs (j_i, alpha_i), where each alpha_i
            is a nonzero integer, j_i is an integer (the j_i-th Manin
            symbol), and the sum alpha_i*x_{j_i} is the image of the
            j-th Manin symbol under the right action of the matrix
            [a,b;c,d].  Here the right action of g=[a,b;c,d] on a
            Manin symbol [P(X,Y),(u,v)] is [P(aX+bY,cX+dY),(u,v)*g].
        """
        a, b, c, d = m[0], m[1], m[2], m[3]
        i, u, v = self[j]
        P = apply_to_monomial(i, self._weight-2, a, b, c, d)
        m = self.index((0, u*a+v*c, u*b+v*d))
        if m == -1:
            return []
        r = len(self.__syms)
        return [(m + r*k, P[k]) for k in range(self._weight-2+1)
                            if P[k] != 0]


    def level(self):
        return self.__level

    def normalize(self, x):
        u,v = self.__syms.normalize(x[1],x[2])
        return (x[0],u,v)


class ManinSymbolList_gamma0(ManinSymbolList_group):
    """
    EXAMPLE:
        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
        sage: m = ManinSymbolList_gamma0(5,2); m
        Manin Symbol List of weight 2 for Gamma0(5)
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (1,4)]
        sage: m = ManinSymbolList_gamma0(6,4); m
        Manin Symbol List of weight 4 for Gamma0(6)
        sage: len(m)
        36
    """
    def __init__(self, level, weight):
        ManinSymbolList_group.__init__(self, level, weight, p1list.P1List(level))

    def __repr__(self):
        return "Manin Symbol List of weight %s for Gamma0(%s)"%(
                    self.weight(), self.level())


class ManinSymbolList_gamma1(ManinSymbolList_group):
    def __init__(self, level, weight):
        ManinSymbolList_group.__init__(self, level, weight, g1list.G1list(level))

    def __repr__(self):
        return "Manin Symbol List of weight %s for Gamma1(%s)"%(
                    self.weight(), self.level())


class ManinSymbolList_gamma_h(ManinSymbolList_group):
    def __init__(self, group, weight):
        ManinSymbolList_group.__init__(self, group.level(), weight, ghlist.GHlist(group))

    def __repr__(self):
        return "Manin Symbol List of weight %s for %s"%(
                    self.weight(), self.group())


class ManinSymbolList_character(ManinSymbolList):
    """
    List of Manin Symbols with character.

    ManinSymbolList_character(character, weight):
    INPUT:
        character -- a dirichlet character
        weight -- integer weight >= 2
    EXAMPLE:
        sage: eps = DirichletGroup(4).gen(0)
        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
        sage: m = ManinSymbolList_character(eps,2); m
        Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (2,1)]
    """
    def __init__(self, character, weight):
        self.__level = character.modulus()
        self.__P1 = p1list.P1List(self.level())

        # We make a copy of the character *only* to program around what seems
        # to be a bug in the cPickle module in some obscure case.
        # If we don't due this, then this doctest fails.
        #    sage: M = ModularSymbols(DirichletGroup(5).0)
        #    sage: loads(dumps(M)) == M

        self.__character = character.__copy__()

        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight-2+1) \
                            for u, v in self.__P1.list()]
        self.__list = L
        ManinSymbolList.__init__(self, weight, L)

    def __repr__(self):
        return "Manin Symbol List of weight %s for Gamma1(%s) with character %s"%(
                    self.weight(), self.level(), self.character())

    def apply(self, j, m):
        """
        Apply the matrix m=[a,b,c,d] to the j-th Manin symbol.
        INPUT:
            j -- integer
            m = [a, b, c, d] a list of 4 integers.
        OUTPUT:
            a list of pairs (j, alpha_i), where each alpha_i is an integer,
            j is an integer (the j-th Manin symbol), and the sum alpha_i*x_i
            is the image of self under the right action of the matrix [a,b;c,d].
            Here the right action of g=[a,b;c,d] on a Manin symbol [P(X,Y),(u,v)]
            is [P(aX+bY,cX+dY),(u,v)*g].
        EXAMPLES:
            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4)
            sage: m[6]
            (1, 0, 1)
            sage: m.apply(4, [1,0,0,1])
            [(4, 1)]
            sage: m.apply(1, [-1,0,0,1])
            [(1, -1)]
        """
        a, b, c, d = m[0], m[1], m[2], m[3]
        i, u, v = self[int(j)]
        P = apply_to_monomial(i, self._weight-2, a, b, c, d)
        m, s = self.index((0, u*a+v*c, u*b+v*d))
        if m == -1 or s == 0:
            return []
        r = len(self.__P1)
        return [(m + r*k, s*P[k]) for k in range(self._weight-2+1)
                            if P[k] != 0]

    def apply_S(self, j):
        """
        """
        i, u, v = self._list[j]
        k, s = self.index((self._weight-2-i, v, -u))
        if i%2==0:
            return k, s
        else:
            return k, -s

    def apply_I(self, j):
        """
        """
        i, u, v = self._list[j]
        k, s = self.index((i, -u, v))
        if i%2==0:
            return k, s
        else:
            return k, -s

    def apply_T(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v, r = self.__P1.normalize_with_scalar(v,-u-v)
        r = self.__character(r)
        if (k-2) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        for j in range(k-2-i +1):
            m, r = self.index((j, u, v))
            z.append((m,s*r*arith.binomial(k-2-i,j)))
            s *= -1
        return z

    def apply_TT(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v, r = self.__P1.normalize_with_scalar(-u-v,u)
        r = self.__character(r)
        if (k-2-i) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        for j in range(i+1):
            m, r = self.index((k-2-i+j, u, v))
            z.append((m,s*r*arith.binomial(i,j)))
            s *= -1
        return z

    def character(self):
        return self.__character

    def index(self, x):
        """
        Compute the index into the list of standard Manin symbols of a
        symbol that is equivalent, modulo a scalar s, to x.  Returns
        the index and the scalar.

        If x is not in the list, then this function returns -1, 0.

        INPUT:
            x -- 3-tuple of ints. Something equivalent to an element
                 of Manin symbols list, which need not be normalized.
        OUTPUT:
            int -- the index of the Manin symbol equivalent to (i,u,v).
            scalar -- element of the base field or the int 0.
        """
        t, i = search(self._list, x)
        if t: return i, 1
        x, s = self.normalize(x)
        t, i = search(self._list, x)
        if t: return i, s
        return -1, 0

    def level(self):
        return self.__level

    def normalize(self, x):
        u,v,s = self.__P1.normalize_with_scalar(x[1],x[2])
        return (x[0],u,v), self.__character(s)


class x__ManinSymbolList_gamma1(ManinSymbolList):
    """
    List of Manin symbols for Gamma0(N).

    EXAMPLE:
        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
        sage: m = ManinSymbolList_gamma0(5,2); m
        Manin Symbol List of weight 2 for Gamma0(5)
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (1,4)]
        sage: m = ManinSymbolList_gamma0(6,4); m
        Manin Symbol List of weight 4 for Gamma0(6)
        sage: len(m)
        36
    """
    def __init__(self, level, weight):
        self.__level = level
        self.__G1 = g1list.G1list(self.level())
        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight-2+1) for \
                        u, v in self.__G1.list()]
        ManinSymbolList.__init__(self, weight, L)

    def __repr__(self):
        return "Manin Symbol List of weight %s for Gamma1(%s)"%(
                    self.weight(), self.level())

    def apply_S(self, j):
        """
        """
        i, u, v = self._list[j]
        k = self.index((self._weight-2-i, v, -u))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def apply_I(self, j):
        """
        """
        i, u, v = self._list[j]
        k = self.index((i, -u, v))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def apply_J(self, j):
        """
        Apply 2x2 matrix J = [-1,0,0,-1].
        """
        i, u, v = self._list[j]
        N = self.__level
        return self.index((i, -u, -v)), 1

    def apply_T(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__G1.normalize(v,-u-v)
        if (k-2) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        for j in range(k-2-i +1):
            m = self.index((j, u, v))
            z.append((m,s*arith.binomial(k-2-i,j)))
            s *= -1
        return z

    def apply_TT(self, j):
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__G1.normalize(-u-v,u)
        if (k-2-i) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        for j in range(i+1):
            m = self.index((k-2-i+j, u, v))
            z.append((m,s*arith.binomial(i,j)))
            s *= -1
        return z

    def apply(self, j, m):
        """
        Apply the matrix m=[a,b,c,d] to the j-th Manin symbol.
        INPUT:
            j -- integer
            m = [a, b, c, d] a list of 4 integers.
        OUTPUT:
            a list of pairs (j, alpha_i), where each alpha_i is an integer,
            j is an integer (the j-th Manin symbol), and the sum alpha_i*x_i
            is the image of self under the right action of the matrix [a,b;c,d].
            Here the right action of g=[a,b;c,d] on a Manin symbol [P(X,Y),(u,v)]
            is [P(aX+bY,cX+dY),(u,v)*g].

        """
        a, b, c, d = m[0], m[1], m[2], m[3]
        i, u, v = self[j]
        P = apply_to_monomial(i, self._weight-2, a, b, c, d)
        m = self.index((0, u*a+v*c, u*b+v*d))
        if m == -1:
            return []
        r = len(self.__G1)
        return [(m + r*k, P[k]) for k in range(self._weight-2+1)
                            if P[k] != 0]

    def level(self):
        return self.__level

    def normalize(self, x):
        u,v = self.__G1.normalize(x[1],x[2])
        return (x[0],u,v)


class ManinSymbol(SageObject):
    r"""
    A Manin symbol $[X^i\cdot Y^{k-2-i},(u,v)]$.
    """
    def __init__(self, parent, t):
        r"""
        Create a Manin symbol $[X^i*Y^{k-2-i},(u,v)]$, where $k$ is
        the weight.

        INPUT:
            parent -- ManinSymbolList
            t -- a 3-tuple (i,u,v) of int's.
        """
        if not isinstance(parent, ManinSymbolList):
            raise TypeError, "parent (=%s) must be of type ManinSymbolList."%(
                parent)
        self.__parent = parent
        if not isinstance(t, tuple):
            raise TypeError, "t (=%s) must be of type tuple."%t
        self.__t = t

    def tuple(self):
        return self.__t

    def __get_i(self):
        return self.__t[0]
    i = property(__get_i)
    def __get_u(self):
        return self.__t[1]
    u = property(__get_u)
    def __get_v(self):
        return self.__t[2]
    v = property(__get_v)

    def _repr_(self):
        if self.weight() > 2:
            polypart = _print_polypart(self.i, self.weight()-2-self.i)
            return "[%s,(%s,%s)]"%\
                   (polypart, self.u, self.v)
        return "(%s,%s)"%(self.u, self.v)

    def _latex_(self):
        return self._repr_()

    def __cmp__(self, other):
        if not isinstance(other, ManinSymbol):
            return -1
        if self.__t == other.__t:
            return 0
        return 1

    def __mul__(self, matrix):
        if self.weight > 2:
            raise NotImplementedError
        if sage.matrix.all.is_Matrix(matrix):
            assert matrix.nrows() == 2 and matrix.ncols()==2, "matrix must be 2x2"
            matrix = matrix.list()
        return ManinSymbol(self.weight, \
                           matrix[0]*self.u + matrix[2]*self.v,\
                           matrix[1]*self.u + matrix[3]*self.v)
        raise ArithmeticError, "Multiplication of %s by %s not defined."%(self, matrix)


    def apply(self, a,b,c,d):
        """
        Return the image of self under the matrix [a,b;c,d].

        INPUT:
            a, b, c, d -- integers

        OUTPUT:
            a list of pairs (alpha_i, x_i), where each alpha_i is an integer,
            x_i is a Manin symbol, and the sum alpha_i*x_i is the image of
            self under the right action of the matrix [a,b;c,d].
            Here the right action of g=[a,b;c,d] on a Manin symbol [P(X,Y),(u,v)]
            is [P(aX+bY,cX+dY),(u,v)*g].
        """
        raise NotImplementedError

    def copy(self):
        return ManinSymbol(self.parent, self.u, self.v, self.i)

    def lift_to_sl2z(self, N):
        """
        If this Manin symbol is (c,d) viewed modulo N, this function
        computes and returns a list [a,b, c',d'] that defines a 2x2
        matrix with determinant 1 and integer entries, such that
        c=c'(mod N) and d=d'(mod N).
        """
        if N == 1:
            return [1,0,0,1]
        c = self.u
        d = self.v
        g, z1, z2 = arith.XGCD(c,d)

        # We're lucky: z1*c + z2*d = 1.
        if g==1:
            return [z2, -z1, c, d]

        # Have to try harder.
        if c == 0:
            c += N;
        if d == 0:
            d += N;
        m = c;

        # compute prime-to-d part of m.
        while True:
            g = arith.GCD(m,d)
            if g == 1:
                break
            m //= g

        # compute prime-to-N part of m.
        while True:
            g = arith.GCD(m,N);
            if g == 1:
                break
            m //= g
        d += N*m
        g, z1, z2 = arith.XGCD(c,d)
        assert g==1
        return [z2, -z1, c, d]

    def endpoints(self, N=None):
        r"""
        Returns cusps alpha, beta such that this Manin symbol, viewed
        as a symbol for level N, is $X^i*Y^{k-2-i} \{alpha, beta\}$.
        """
        if N is None:
            N = self.parent().level()
        else:
            N=int(N)
            if N < 1:
                raise ArithmeticError, "N must be positive"
        a,b,c,d = self.lift_to_sl2z(N)
        return cusps.Cusp(b,d), cusps.Cusp(a,c)

    def parent(self):
        return self.__parent

    def weight(self):
        return self.__parent.weight()

    def modular_symbol_rep(self):
        """
        Returns a representation of self as a formal sum of modular symbols.
        (The result is not cached.)
        """
        # TODO: It would likely be much better to do this slightly more directly
        x = sage.modular.modsym.modular_symbols.ModularSymbol(\
              self.__parent, self.i, 0, rings.infinity)
        a,b,c,d = self.lift_to_sl2z(self.__parent.level())
        return x.apply([a,b,c,d])

def apply_to_monomial(i, j, a, b, c, d):
    """
    Returns a list of the coefficients of
    $$
             (aX + bY)^i (cX + dY)^{j-i},
    $$
    where $0 \leq i \leq j$, and $a,b,c,d$ are integers.

    One should think of $j$ as being $k-2$ for the application to
    modular symbols.

    INPUT:
        i, j, a, b, c, d -- all ints

    OUTPUT:
        list of ints, which are the coefficients
        of Y^j, Y^(j-1)*X, ..., X^j, respectively.

    EXAMPLE:
    We compute that $(X+Y)^2(X-Y) = X^3 + X^2Y - XY^2 - Y^3$.
        sage: from sage.modular.modsym.manin_symbols import apply_to_monomial
        sage: apply_to_monomial(2, 3, 1,1,1,-1)
        [-1, -1, 1, 1]
        sage: apply_to_monomial(5, 8, 1,2,3,4)
        [2048, 9728, 20096, 23584, 17200, 7984, 2304, 378, 27]
        sage: apply_to_monomial(6,12, 1,1,1,-1)
        [1, 0, -6, 0, 15, 0, -20, 0, 15, 0, -6, 0, 1]
    """
    f = R([b,a])
    g = R([d,c])
    if i < 0 or j-i < 0:
        raise ValueError, "i (=%s) and j-i (=%s) must both be nonnegative."%(i,j-i)
    h = (f**i)*(g**(j-i))
    return [int(h[k]) for k in xrange(j+1)]

def _print_polypart(i, j):
    if i > 1:
        xpart = "X^%s"%i
    elif i == 1:
        xpart = "X"
    else:
        xpart = ""
    if j > 1:
        ypart = "Y^%s"%j
    elif j == 1:
        ypart = "Y"
    else:
        ypart = ""
    if len(xpart) > 0 and len(ypart) > 0:
        times = "*"
    else:
        times = ""
    if len(xpart + ypart) > 0:
        polypart = "%s%s%s"%(xpart, times, ypart)
    else:
        polypart = ""
    return polypart
