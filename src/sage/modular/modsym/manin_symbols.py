# -*- coding: utf-8 -*-
"""
Manin symbols

This module defines the class ManinSymbol.  A Manin Symbol of
weight `k`, level `N` has the form `[P(X,Y),(u:v)]` where
`P(X,Y)\in\mathbb{Z}[X,Y]` is homogeneous of weight `k-2` and
`(u:v)\in\mathbb{P}^1(\mathbb{Z}/N\mathbb{Z}).`  The ManinSymbol class
holds a "monomial Manin Symbol" of the simpler form
`[X^iY^{k-2-i},(u:v)]`, which is stored as a triple `(i,u,v)`; the
weight and level are obtained from the parent structure, which is a
ManinSymbolList.

Integer matrices `[a,b;c,d]` act on Manin Symbols on the right,
sending `[P(X,Y),(u,v)]` to `[P(aX+bY,cX+dY),(u,v)g]`.  Diagonal
matrices (with `b=c=0`, such as `I=[-1,0;0,1]` and `J=[-1,0;0,-1]`)
and anti-diagonal matrices (with `a=d=0`, such as `S=[0,-1;1,0]`) map
monomial Manin Symbols to monomial Manin Symbols, up to a scalar
factor.  For general matrices (such as `T=[0,1,-1,-1]` and
`T^2=[-1,-1;0,1]`) the image of a monomial Manin Symbol is expressed
as a formal sum of monomial Manin Symbols, with integer coefficients.

There are various different classes holding lists of Manin symbols of
different types.  The hierarchy is as follows:

::

    class ManinSymbolList(SageObject)

    class ManinSymbolList_group(ManinSymbolList)
        class ManinSymbolList_gamma0(ManinSymbolList_group)
        class ManinSymbolList_gamma1(ManinSymbolList_group)
        class ManinSymbolList_gamma_h(ManinSymbolList_group)

    class ManinSymbolList_character(ManinSymbolList)

"""

#*****************************************************************************
#       Sage: System for Algebra and Geometry Experimentation
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

import sage.matrix.all
import sage.modular.cusps as cusps
import sage.modular.modsym.p1list as p1list
import sage.modular.modsym.g1list as g1list
import sage.modular.modsym.ghlist as ghlist
import sage.modular.modsym.modular_symbols
import sage.rings.arith as arith
import sage.rings.all as rings

from sage.structure.sage_object import SageObject

from apply import apply_to_monomial

def is_ManinSymbol(x):
    """
    Returns True if ``x`` is a ManinSymbol.

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0, is_ManinSymbol
        sage: m = ManinSymbolList_gamma0(6, 4)
        sage: is_ManinSymbol(m[3])
        False
        sage: s = ManinSymbol(m,m[3])
        sage: s
        [Y^2,(1,2)]
        sage: is_ManinSymbol(s)
        True

    """
    return isinstance(x, ManinSymbol)

class ManinSymbolList(SageObject):
    """
    Base class for lists of all Manin symbols for a given weight, group or character.
    """
    def __init__(self, weight, list):
        """
        Constructor for a ManinSymbolList.

        INPUT:

        -  ``weight``- the weight of the symbols.
        -  ``list``- the list of symbols.

        On construction, a ManinSymbolList constructs a dict for
        rapid determination of the index of any given symbol.

        This is a base class only; users will only directly construct
        objects in the derived classes ManinSymbolList_gamma0,
        ManinSymbolList_gamma1, ManinSymbolList_gamma_h,
        ManinSymbolList_gamma_character.  Many standard methods are
        only implemented in the derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: ManinSymbolList(6,P1List(11))
            <class 'sage.modular.modsym.manin_symbols.ManinSymbolList'>

        """
        self._weight = weight
        self._list = list
        self._index = dict([(list[i],i) for i in range(len(list))])

    def __cmp__(self, right):
        """
        Comparison function for ManinSymbolList objects.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m1 = ManinSymbolList(6,P1List(11))
            sage: m2 = ManinSymbolList(6,P1List(13))
            sage: m3 = ManinSymbolList(4,P1List(11))
            sage: m1 < m2
            True
            sage: m2 < m3
            False
            sage: m1 < m3
            False
        """
        if not isinstance(right, ManinSymbolList):
            return cmp(type(self), type(right))
        return cmp((self._weight, self._list), (right._weight, right._list))

    def __getitem__(self, n):
        """
        Returns the `n`th ManinSymbol in this ManinSymbolList

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m[4]
            (1, 3)
        """
        return self._list[n]

    def __len__(self):
        """
        Returns the length of this ManinSymbolList

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: len(m)
            12
        """
        return len(self._list)

    def apply(self, j, X):
        """
        Apply the matrix `X=[a,b;c,d]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply(10, [1,2,0,1])
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes

        """
        raise NotImplementedError, "Only implemented in derived classes"

    def _apply_S_only_0pm1(self):
        """
        Return True if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            True
            sage: eps = DirichletGroup(7).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            False
        """
        return False # derived classes could overload and put True

    def apply_S(self, j):
        """
        Apply the matrix `S=[0,-1;1,0]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_S(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError, "Only implemented in derived classes"

    def apply_I(self, j):
        """
        Apply the matrix `I=[-1,0;0,1]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_I(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError, "Only implemented in derived classes"

    def apply_T(self, j):
        """
        Apply the matrix `T=[0,1;-1,-1]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_T(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError, "Only implemented in derived classes"

    def apply_TT(self, j):
        """
        Apply the matrix `TT=T^2=[-1,-1;0,1]` to the `j`-th Manin symbol.

        Implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.apply_TT(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in derived classes
        """
        raise NotImplementedError, "Only implemented in derived classes"

    def index(self, x):
        """
        Return the index of ``x`` in the list of Manin symbols, where ``x`` is a
        3-tuple of ints. If ``x`` is not in the list, then return -1.

        INPUT:

        - ``x`` - 3-tuple of integers, `(i,u,v)` defining a valid Manin symbol, which need not be normalized.

        OUTPUT:

        (``int``) the index of the normalized Manin symbol equivalent to `(i,u,v)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.index(m[2])
            2
            sage: all([i == m.index(m[i]) for i in xrange(len(m))])
            True
        """
        if x in self._index:
            return self._index[x]
        x = self.normalize(x)
        try:
            return self._index[x]
        except KeyError:
            return -1

    def manin_symbol_list(self):
        """
        Returns all the ManinSymbols in this ManinSymbolList as a list

        Cached for subsequent calls.

        OUTPUT:

        a list of ``ManinSymbol`` objects, which is a copy of the complete list
        of Manin symbols.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.manin_symbol_list() # not implemented for the base class

        ::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: m.manin_symbol_list()
            [[Y^2,(0,1)],
            [Y^2,(1,0)],
            [Y^2,(1,1)],
            ...
            [X^2,(3,1)],
            [X^2,(3,2)]]

        """
        import copy
        try:
            return copy.copy(self.__manin_symbol_list)
        except AttributeError:
            self.__manin_symbol_list = [self.manin_symbol(i) \
                                        for i in xrange(len(self))]
        return copy.copy(self.__manin_symbol_list)

    def manin_symbol(self, i):
        """
        Returns the i'th ManinSymbol in this ManinSymbolList.

        INPUT:


        - ``i`` - integer, a valid index of a symbol in this list.


        OUTPUT:

        ``ManinSymbol`` - the `i`'th Manin symbol in the list.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.manin_symbol(3) # not implemented for base class

        ::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: s = m.manin_symbol(3); s
            [Y^2,(1,2)]
            sage: type(s)
            <class 'sage.modular.modsym.manin_symbols.ManinSymbol'>
        """
        return ManinSymbol(self, self._list[int(i)])

    def normalize(self, x):
        """
        Returns a normalized ManinSymbol from x.

        To be implemented in derived classes.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList
            sage: m = ManinSymbolList(6,P1List(11))
            sage: m.normalize((0,6,7)) # not implemented in base class

        """
        raise NotImplementedError, "Only implemented in derived classes"

    def weight(self):
        """
        Returns the weight of the ManinSymbols in this ManinSymbolList.


        OUTPUT:

        ``integer`` - the weight of the Manin symbols in the list.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(6, 4)
            sage: m.weight()
            4
        """
        return self._weight

class ManinSymbolList_group(ManinSymbolList):
    """
    Base class for Manin symbol lists for a given group.

    INPUT:


    -  ``level`` - integer level

    -  ``weight`` - integer weight

    -  ``syms`` - something with a normalize and list
       method, e.g., P1List.

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_group
        sage: ManinSymbolList_group(11, 2, P1List(11))
        <class 'sage.modular.modsym.manin_symbols.ManinSymbolList_group'>

    """
    def __init__(self, level, weight, syms):
        """
        Constructor for class ManinSymbolList_group.

        INPUT:

        -  ``level`` - integer level

        -  ``weight`` - integer weight

        -  ``syms`` - something with a normalize and list method, e.g., P1List.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_group
            sage: ManinSymbolList_group(11, 2, P1List(11))
            <class 'sage.modular.modsym.manin_symbols.ManinSymbolList_group'>

        """
        self.__level = level
        self.__syms = syms  # syms is anything with a normalize and list method.

        # The list returned from P1List is guaranteed to be sorted.
        # Thus each list constructed below is also sorted.  This is
        # important since the index function assumes the list is sorted.
        L = [(i, u, v) for i in range(weight - 2 + 1) \
                            for u, v in syms.list()]
        ManinSymbolList.__init__(self, weight, L)

    def level(self):
        """
        Return the level of this ManinSymbolList.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: ManinSymbolList_gamma0(5,2).level()
            5

        ::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma1
            sage: ManinSymbolList_gamma1(51,2).level()
            51

        ::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma_h
            sage: ManinSymbolList_gamma_h(GammaH(117, [4]),2).level()
            117
        """
        return self.__level

    def apply_S(self, j):
        """
        Apply the matrix `S=[0,-1,1,0]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (int) a symbol index

        OUTPUT:

        ``(k, s)`` where k is the index of the symbol obtained by acting on the
        `j`'th symbol with `S`, and `s` is the parity of of the `j`'th symbol
        (a Python ``int``, either 1 or -1).

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_S(4)
            (40, 1)
            sage: [m.apply_S(i) for i in xrange(len(m))]
            [(37, 1),
            (36, 1),
            (41, 1),
            (39, 1),
            (40, 1),
            (38, 1),
            (31, -1),
            (30, -1),
            (35, -1),
            (33, -1),
            (34, -1),
            (32, -1),
            ...
            (4, 1),
            (2, 1)]
        """
        i, u, v = self._list[j]
        k = self.index((self._weight-2-i, v, -u))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def _apply_S_only_0pm1(self):
        """
        Return True if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: ManinSymbolList_gamma0(5,8)._apply_S_only_0pm1()
            True
        """
        return True

    def apply_I(self, j):
        """
        Apply the matrix `I=[-1,0,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (int) a symbol index

        OUTPUT:

        ``(k, s)`` where k is the index of the symbol obtained by acting on the
        `j`'th symbol with `I`, and `s` is the parity of of the `j`'th symbol
        (a Python ``int``, either 1 or -1)

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_I(4)
            (3, 1)
            sage: [m.apply_I(i) for i in xrange(10)]
            [(0, 1),
            (1, 1),
            (5, 1),
            (4, 1),
            (3, 1),
            (2, 1),
            (6, -1),
            (7, -1),
            (11, -1),
            (10, -1)]
        """
        i, u, v = self._list[j]
        k = self.index((i, -u, v))
        if i%2==0:
            return k, 1
        else:
            return k, -1

    def apply_T(self, j):
        """
        Apply the matrix `T=[0,1,-1,-1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (int) a symbol index

        OUTPUT: see documentation for apply()

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_T(4)
            [(3, 1), (9, -6), (15, 15), (21, -20), (27, 15), (33, -6), (39, 1)]
            sage: [m.apply_T(i) for i in xrange(10)]
            [[(5, 1), (11, -6), (17, 15), (23, -20), (29, 15), (35, -6), (41, 1)],
            [(0, 1), (6, -6), (12, 15), (18, -20), (24, 15), (30, -6), (36, 1)],
            [(4, 1), (10, -6), (16, 15), (22, -20), (28, 15), (34, -6), (40, 1)],
            [(2, 1), (8, -6), (14, 15), (20, -20), (26, 15), (32, -6), (38, 1)],
            [(3, 1), (9, -6), (15, 15), (21, -20), (27, 15), (33, -6), (39, 1)],
            [(1, 1), (7, -6), (13, 15), (19, -20), (25, 15), (31, -6), (37, 1)],
            [(5, 1), (11, -5), (17, 10), (23, -10), (29, 5), (35, -1)],
            [(0, 1), (6, -5), (12, 10), (18, -10), (24, 5), (30, -1)],
            [(4, 1), (10, -5), (16, 10), (22, -10), (28, 5), (34, -1)],
            [(2, 1), (8, -5), (14, 10), (20, -10), (26, 5), (32, -1)]]
        """
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__syms.normalize(v,-u-v)
        if (k-2) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        a = rings.ZZ(k-2-i)
        for j in range(k-2-i+1):
            m = self.index((j, u, v))
            z.append((m, s * a.binomial(j)))
            s *= -1
        return z

    def apply_TT(self, j):
        """
        Apply the matrix `TT=[-1,-1,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (int) a symbol index

        OUTPUT: see documentation for apply()

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply_TT(4)
            [(38, 1)]
            sage: [m.apply_TT(i) for i in xrange(10)]
            [[(37, 1)],
            [(41, 1)],
            [(39, 1)],
            [(40, 1)],
            [(38, 1)],
            [(36, 1)],
            [(31, -1), (37, 1)],
            [(35, -1), (41, 1)],
            [(33, -1), (39, 1)],
            [(34, -1), (40, 1)]]
        """
        k = self._weight
        i, u, v = self._list[j]
        u, v = self.__syms.normalize(-u-v,u)
        if (k-2-i) % 2 == 0:
            s = 1
        else:
            s = -1
        z = []
        a = rings.ZZ(i)
        for j in range(i+1):
            m = self.index((k-2-i+j, u, v))
            z.append((m, s * a.binomial(j)))
            s *= -1
        return z

    def apply(self, j, m):
        r"""
        Apply the matrix `m=[a,b;c,d]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (int) a symbol index

        - ``m = [a, b, c, d]`` a list of 4 integers, which defines a 2x2 matrix.


        OUTPUT:

        a list of pairs `(j_i, \alpha_i)`, where each `\alpha_i` is a nonzero
        integer, `j_i` is an integer (index of the `j_i`-th Manin symbol), and
        `\sum_i \alpha_i\*x_{j_i}` is the image of the j-th Manin symbol under
        the right action of the matrix [a,b;c,d]. Here the right action of
        g=[a,b;c,d] on a Manin symbol `[P(X,Y),(u,v)]` is
        `[P(aX+bY,cX+dY),(u,v)\*g]`.

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: m.apply(40,[2,3,1,1])
            [(0, 729), (6, 2916), (12, 4860), (18, 4320), (24, 2160), (30, 576), (36, 64)]

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

    def normalize(self, x):
        """
        Returns the normalization of the ModSym ``x`` with respect to this
        list.

        INPUT:

        - ``x`` - (3-tuple of ints) a tuple defining a ManinSymbol.

        OUTPUT:

        ``(i,u,v)`` - (3-tuple of ints) another tuple defining the associated
        normalized ManinSymbol.

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: [m.normalize(s.tuple()) for s in m.manin_symbol_list()][:10]
            [(0, 0, 1),
            (0, 1, 0),
            (0, 1, 1),
            (0, 1, 2),
            (0, 1, 3),
            (0, 1, 4),
            (1, 0, 1),
            (1, 1, 0),
            (1, 1, 1),
            (1, 1, 2)]
        """
        u,v = self.__syms.normalize(x[1],x[2])
        return (x[0],u,v)


class ManinSymbolList_gamma0(ManinSymbolList_group):
    r"""
    Class for Manin Symbols for `\Gamma_0(N)`.

    INPUT:

    - ``level`` - (integer): the level.

    - ``weight`` - (integer): the weight.

    EXAMPLE::

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
        """
        Constructor for a ModularSymbolList for Gamma_0(N)

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: M11 = ManinSymbolList_gamma0(11,2)
            sage: M11
            Manin Symbol List of weight 2 for Gamma0(11)
            sage: M11 == loads(dumps(M11))
            True
        """
        ManinSymbolList_group.__init__(self, level, weight, p1list.P1List(level))

    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
            sage: M11 = ManinSymbolList_gamma0(11,2)
            sage: str(M11)
            'Manin Symbol List of weight 2 for Gamma0(11)'

        """
        return "Manin Symbol List of weight %s for Gamma0(%s)"%(
                    self.weight(), self.level())


class ManinSymbolList_gamma1(ManinSymbolList_group):
    r"""
    Class for Manin Symbols for `\Gamma_1(N)`.

    INPUT:

    - ``level`` - (integer): the level.

    - ``weight`` - (integer): the weight.

    EXAMPLE::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma1
        sage: m = ManinSymbolList_gamma1(5,2); m
        Manin Symbol List of weight 2 for Gamma1(5)
        sage: m.manin_symbol_list()
        [(0,1),
        (0,2),
        (0,3),
        ...
        (4,3),
        (4,4)]
        sage: m = ManinSymbolList_gamma1(6,4); m
        Manin Symbol List of weight 4 for Gamma1(6)
        sage: len(m)
        72
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, level, weight):
        """
        Constructor for a ModularSymbolList for `\Gamma_0(N)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma1
            sage: M11 = ManinSymbolList_gamma1(11,2)
            sage: M11
            Manin Symbol List of weight 2 for Gamma1(11)
        """
        ManinSymbolList_group.__init__(self, level, weight, g1list.G1list(level))

    def __repr__(self):
        """
        Return the string representation of this ManinSymbbol list.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma1
            sage: M11 = ManinSymbolList_gamma1(11,4)
            sage: str(M11)
            'Manin Symbol List of weight 4 for Gamma1(11)'
        """
        return "Manin Symbol List of weight %s for Gamma1(%s)"%(
                    self.weight(), self.level())


class ManinSymbolList_gamma_h(ManinSymbolList_group):
    r"""
    Class for Manin Symbols for `\Gamma_H(N)`.

    INPUT:

    - ``group`` - (integer): the congruence subgroup.

    - ``weight`` - (integer): the weight.

    EXAMPLE::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma_h
        sage: G = GammaH(117, [4])
        sage: m = ManinSymbolList_gamma_h(G,2); m
        Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(117) with H generated by [4]
        sage: m.manin_symbol_list()[100:110]
        [(1,88),
        (1,89),
        (1,90),
        (1,91),
        (1,92),
        (1,93),
        (1,94),
        (1,95),
        (1,96),
        (1,97)]
        sage: len(m.manin_symbol_list())
        2016
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, group, weight):
        r"""
        Constructor for Manin Symbols for `\Gamma_H(N)`.

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma_h
            sage: G = GammaH(117, [4])
            sage: m = ManinSymbolList_gamma_h(G,2); m
            Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(117) with H generated by [4]
            """
        self.__group = group
        ManinSymbolList_group.__init__(self, group.level(), weight, ghlist.GHlist(group))

    def group(self):
        """
        Return the group associated to self.

        EXAMPLES::

            sage: ModularSymbols(GammaH(12, [5]), 2).manin_symbols().group()
            Congruence Subgroup Gamma_H(12) with H generated by [5]
        """
        return self.__group

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: ModularSymbols(GammaH(12, [5]), 2).manin_symbols().__repr__()
            'Manin Symbol List of weight 2 for Congruence Subgroup Gamma_H(12) with H generated by [5]'
        """
        return "Manin Symbol List of weight %s for %s"%(
                    self.weight(), self.group())


class ManinSymbolList_character(ManinSymbolList):
    """
    List of Manin Symbols with character.

    INPUT:

    -  ``character`` - (DirichletCharacter) the Dirichlet character.

    -  ``weight`` - (integer) the weight.

    EXAMPLE::

        sage: eps = DirichletGroup(4).gen(0)
        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
        sage: m = ManinSymbolList_character(eps,2); m
        Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
        sage: m.manin_symbol_list()
        [(0,1), (1,0), (1,1), (1,2), (1,3), (2,1)]
        sage: m == loads(dumps(m))
        True
    """
    def __init__(self, character, weight):
        """
        Constructor for objects of class ManinSymbolList_character

        INPUT:


        -  ``character`` - (DirichletCharacter) the Dirichlet character.

        -  ``weight`` - (integer) the weight.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.manin_symbol_list()
            [(0,1), (1,0), (1,1), (1,2), (1,3), (2,1)]

        """
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
        """
        Standard function returning string representation.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: str(m) # indirect doctest
            'Manin Symbol List of weight 2 for Gamma1(4) with character [-1]'
        """
        return "Manin Symbol List of weight %s for Gamma1(%s) with character %s"%(
                    self.weight(), self.level(), self.character()._repr_short_())

    def level(self):
        """
        Return the level of this ManinSymbolList.

        OUTPUT:

        ``integer`` - the level of the symbols in this list.

        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,4).level()
            4
        """
        return self.__level

    def apply(self, j, m):
        """
        Apply the integer matrix `m=[a,b;c,d]` to the `j`-th Manin symbol.

        INPUT:


        - ``j`` (integer): the index of the symbol to act on.

        - ``m`` (list of ints):  `[a,b,c,d]` where `m = [a, b; c, d]` is the matrix to be applied.


        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of self under the right action
        of the matrix `[a,b;c,d]`. Here the right action of
        `g=[a,b;c,d]` on a Manin symbol `[P(X,Y),(u,v)]` is by
        definition `[P(aX+bY,cX+dY),(u,v)*g]`.

        EXAMPLES::

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
        Apply the matrix `S=[0,1;-1,0]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (integer) a symbol index.

        OUTPUT:

        ``(k, s)`` where `k` is the index of the symbol obtained by acting
        on the `j`'th symbol with `S`, and `s` is the parity of the
        `j`'th symbol.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_S(4)
            (2, -1)
            sage: [m.apply_S(i) for i in xrange(len(m))]
            [(1, 1), (0, -1), (4, 1), (5, -1), (2, -1), (3, 1)]
        """
        i, u, v = self._list[j]
        k, s = self.index((self._weight-2-i, v, -u))
        if i%2==0:
            return k, s
        else:
            return k, -s

    def _apply_S_only_0pm1(self):
        """
        Return True if the coefficient when applying the S relation is
        always 0, 1, or -1.  This is useful for optimizing code in
        relation_matrix.py.

        EXAMPLES::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: ManinSymbolList_character(eps,2)._apply_S_only_0pm1()
            True
            sage: ManinSymbolList_character(DirichletGroup(13).0,2)._apply_S_only_0pm1()
            False
        """
        return self.__character.order() <= 2

    def apply_I(self, j):
        """
        Apply the matrix `I=[-1,0,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (integer) a symbol index

        OUTPUT:

        ``(k, s)`` where `k` is the index of the symbol obtained by acting
        on the `j`'th symbol with `I`, and `s` is the parity of of the
        `j`'th symbol.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_I(4)
            (2, -1)
            sage: [m.apply_I(i) for i in xrange(len(m))]
            [(0, 1), (1, -1), (4, -1), (3, -1), (2, -1), (5, 1)]
        """
        i, u, v = self._list[j]
        k, s = self.index((i, -u, v))
        if i%2==0:
            return k, s
        else:
            return k, -s

    def apply_T(self, j):
        """
        Apply the matrix `T=[0,1,-1,-1]` to the j-th Manin symbol.

        INPUT:

        - ``j`` - (integer) a symbol index.

        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of self under the right action
        of the matrix `T`.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_T(4)
            [(1, -1)]
            sage: [m.apply_T(i) for i in xrange(len(m))]
            [[(4, 1)], [(0, -1)], [(3, 1)], [(5, 1)], [(1, -1)], [(2, 1)]]
        """
        k = self._weight
        i, u, v = self._list[j]
        u, v, r = self.__P1.normalize_with_scalar(v,-u-v)
        r = self.__character(r)
        if (k-2) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        a = rings.ZZ(k-2-i)
        for j in range(k-2-i+1):
            m, r = self.index((j, u, v))
            z.append((m, s * r * a.binomial(j)))
            s *= -1
        return z

    def apply_TT(self, j):
        """
        Apply the matrix `TT=[-1,-1,0,1]` to the `j`-th Manin symbol.

        INPUT:

        - ``j`` - (integer) a symbol index

        OUTPUT:

        A list of pairs `(j, c_i)`, where each `c_i` is an
        integer, `j` is an integer (the `j`-th Manin symbol), and the
        sum `c_i*x_i` is the image of self under the right action
        of the matrix `T^2`.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.apply_TT(4)
            [(0, 1)]
            sage: [m.apply_TT(i) for i in xrange(len(m))]
            [[(1, -1)], [(4, -1)], [(5, 1)], [(2, 1)], [(0, 1)], [(3, 1)]]
        """
        k = self._weight
        i, u, v = self._list[j]
        u, v, r = self.__P1.normalize_with_scalar(-u-v,u)
        r = self.__character(r)
        if (k-2-i) % 2 == 0:
            s = r
        else:
            s = -r
        z = []
        a = rings.ZZ(i)
        for j in range(i+1):
            m, r = self.index((k-2-i+j, u, v))
            z.append((m, s * r * a.binomial(j)))
            s *= -1
        return z

    def character(self):
        """
        Return the character of this ManinSymbolList_character object.

        OUTPUT:

        The Dirichlet character of this Manin symbol list.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,2); m
            Manin Symbol List of weight 2 for Gamma1(4) with character [-1]
            sage: m.character()
            Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1

        """
        return self.__character

    def index(self, x):
        """
        Returns the index in the list of standard Manin symbols of a
        symbol that is equivalent, modulo a scalar `s`, to
        ``x``. Returns the index and the scalar.

        If ``x`` is not in the list, return (-1, 0).

        INPUT:

        - ``x`` - 3-tuple of integers `(i,u,v)`, defining an element of this list of Manin symbols, which need not be normalized.

        OUTPUT:

        ``(i, s)`` where i (``int``) is the index of the Manin symbol
        equivalent to `(i,u,v)` (or -1) and ``s`` is the scalar (an element of
        the base field) or the int 0.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4); m
            Manin Symbol List of weight 4 for Gamma1(4) with character [-1]
            sage: [m.index(s.tuple()) for s in m.manin_symbol_list()]
            [(0, 1),
            (1, 1),
            (2, 1),
            (3, 1),
            ...
            (16, 1),
            (17, 1)]
        """
        if x in self._index:
            return self._index[x], 1
        x, s= self.normalize(x)
        try:
            return self._index[x], s
        except KeyError:
            return -1, 0

    def normalize(self, x):
        """
        Returns the normalization of the Manin Symbol ``x`` with respect to this
        list, together with the normalizing scalar.

        INPUT:

        - ``x`` - 3-tuple of integers ``(i,u,v)``, defining an element of this
          list of Manin symbols, which need not be normalized.

        OUTPUT:

        ``((i,u,v),s)``, where ``(i,u,v)`` is the normalized Manin symbol equivalent
        to ``x``, and ``s`` is the normalizing scalar.

        EXAMPLE::

            sage: eps = DirichletGroup(4).gen(0)
            sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_character
            sage: m = ManinSymbolList_character(eps,4); m
            Manin Symbol List of weight 4 for Gamma1(4) with character [-1]
            sage: [m.normalize(s.tuple()) for s in m.manin_symbol_list()]
            [((0, 0, 1), 1),
            ((0, 1, 0), 1),
            ((0, 1, 1), 1),
            ...
            ((2, 1, 3), 1),
            ((2, 2, 1), 1)]
        """
        u,v,s = self.__P1.normalize_with_scalar(x[1],x[2])
        return (x[0],u,v), self.__character(s)


# class x__ManinSymbolList_gamma1(ManinSymbolList):
#     r"""
#     List of Manin symbols for `\Gamma_1(N)`.

#     EXAMPLE::

#         sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
#         sage: m = ManinSymbolList_gamma0(5,2); m
#         Manin Symbol List of weight 2 for Gamma0(5)
#         sage: m.manin_symbol_list()
#         [(0,1), (1,0), (1,1), (1,2), (1,3), (1,4)]
#         sage: m = ManinSymbolList_gamma0(6,4); m
#         Manin Symbol List of weight 4 for Gamma0(6)
#         sage: len(m)
#         36
#     """
#     def __init__(self, level, weight):
#         r"""
#         Constructor for list of Manin symbols for `\Gamma_1(N)`.
#         """
#         self.__level = level
#         self.__G1 = g1list.G1list(self.level())
#         # The list returned from P1List is guaranteed to be sorted.
#         # Thus each list constructed below is also sorted.  This is
#         # important since the index function assumes the list is sorted.
#         L = [(i, u, v) for i in range(weight-2+1) for \
#                         u, v in self.__G1.list()]
#         ManinSymbolList.__init__(self, weight, L)

#     def __repr__(self):
#         """
#         Returns a string representation for this ManinSymbol list.
#         """
#         return "Manin Symbol List of weight %s for Gamma1(%s)"%(
#                     self.weight(), self.level())

#     def apply_S(self, j):
#         """
#         Apply the matrix `S=[0,1,-1,0]` to the j-th Manin symbol.

#         INPUT:

#         - `j` - (int) a symbol index

#         OUTPUT:

#         (k, s) where k is the index of the symbol obtained by acting
#         on the `j`'th symbol with `S`, and `s` is the parity of of the
#         `j`'th symbol.

#         """
#         i, u, v = self._list[j]
#         k = self.index((self._weight-2-i, v, -u))
#         if i%2==0:
#             return k, 1
#         else:
#             return k, -1

#     def apply_I(self, j):
#         """
#         Apply the matrix `I=[-1,0,0,1]` to the j-th Manin symbol.

#         INPUT:

#         - `j` - (int) a symbol index

#         OUTPUT:

#         (k, s) where k is the index of the symbol obtained by acting
#         on the `j`'th symbol with `I`, and `s` is the parity of of the
#         `j`'th symbol.

#         """
#         i, u, v = self._list[j]
#         k = self.index((i, -u, v))
#         if i%2==0:
#             return k, 1
#         else:
#             return k, -1

#     def apply_J(self, j):
#         """
#         Apply the matrix `J=[-1,0,0,-1]` to the j-th Manin symbol.

#         INPUT:

#         - `j` - (int) a symbol index

#         OUTPUT:

#         (k, s) where k is the index of the symbol obtained by acting
#         on the `j`'th symbol with `J`, and `s` is 1.

#         """
#         """
#         Apply 2x2 matrix J = [-1,0,0,-1].
#         """
#         i, u, v = self._list[j]
#         N = self.__level
#         return self.index((i, -u, -v)), 1

#     def apply_T(self, j):
#         """
#         Apply the matrix `T=[0,1,-1,-1]` to the j-th Manin symbol.

#         INPUT:

#         - `j` - (int) a symbol index

#         OUTPUT: see documentation for apply()

#         """
#         k = self._weight
#         i, u, v = self._list[j]
#         u, v = self.__G1.normalize(v,-u-v)
#         if (k-2) % 2 == 0:
#             s = 1
#         else:
#             s = -1
#         z = []
#         for j in range(k-2-i +1):
#             m = self.index((j, u, v))
#             z.append((m,s*arith.binomial(k-2-i,j)))
#             s *= -1
#         return z

#     def apply_TT(self, j):
#         """
#         Apply the matrix `TT=[-1,-1,0,1]` to the j-th Manin symbol.

#         INPUT:

#         - `j` - (int) a symbol index

#         OUTPUT: see documentation for apply()

#         """
#         k = self._weight
#         i, u, v = self._list[j]
#         u, v = self.__G1.normalize(-u-v,u)
#         if (k-2-i) % 2 == 0:
#             s = 1
#         else:
#             s = -1
#         z = []
#         for j in range(i+1):
#             m = self.index((k-2-i+j, u, v))
#             z.append((m,s*arith.binomial(i,j)))
#             s *= -1
#         return z

#     def apply(self, j, m):
#         """
#         Apply the integer matrix `m=[a,b;c,d]` to the `j`-th Manin symbol.

#         INPUT:

#         - ``j`` - (integer): the index of the symbol to act on.

#         - ``m`` -  `m = [a, b; c, d]`, a list of 4 integers.


#         OUTPUT: a list of pairs (j, c_i), where each c_i is an
#         integer, j is an integer (the j-th Manin symbol), and the sum
#         c_i\*x_i is the image of self under the right action of the
#         matrix [a,b;c,d]. Here the right action of g=[a,b;c,d] on a Manin
#         symbol [P(X,Y),(u,v)] is [P(aX+bY,cX+dY),(u,v)\*g].
#         """
#         a, b, c, d = m[0], m[1], m[2], m[3]
#         i, u, v = self[j]
#         P = apply_to_monomial(i, self._weight-2, a, b, c, d)
#         m = self.index((0, u*a+v*c, u*b+v*d))
#         if m == -1:
#             return []
#         r = len(self.__G1)
#         return [(m + r*k, P[k]) for k in range(self._weight-2+1)
#                             if P[k] != 0]

#     def normalize(self, x):
#         """
#         Returns the normalization of the ModSym x with respect to this
#         list.
#         """
#         u,v = self.__G1.normalize(x[1],x[2])
#         return (x[0],u,v)


class ManinSymbol(SageObject):
    r"""
    A Manin symbol `[X^i\cdot Y^{k-2-i},(u,v)]`.

    INPUT:

    -  ``parent`` - ManinSymbolList.

    -  ``t`` - a 3-tuple `(i,u,v)` of integers.

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
        sage: m = ManinSymbolList_gamma0(5,2)
        sage: s = ManinSymbol(m,(2,2,3)); s
        (2,3)
        sage: s == loads(dumps(s))
        True

        ::

        sage: m = ManinSymbolList_gamma0(5,8)
        sage: s = ManinSymbol(m,(2,2,3)); s
        [X^2*Y^4,(2,3)]

    """
    def __init__(self, parent, t):
        r"""
        Create a Manin symbol `[X^i*Y^{k-2-i},(u,v)]`, where
        `k` is the weight.

        INPUT:


        -  ``parent`` - ManinSymbolList

        -  ``t`` - a 3-tuple (i,u,v) of int's.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,2)
            sage: s = ManinSymbol(m,(2,2,3)); s
            (2,3)

        ::

            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3)); s
            [X^2*Y^4,(2,3)]

        """
        if not isinstance(parent, ManinSymbolList):
            raise TypeError, "parent (=%s) must be of type ManinSymbolList."%(
                parent)
        self.__parent = parent
        if not isinstance(t, tuple):
            raise TypeError, "t (=%s) must be of type tuple."%t
        self.__t = t

    def tuple(self):
        r"""
        Return the 3-tuple `(i,u,v)` of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.tuple()
            (2, 2, 3)

        """
        return self.__t

    def __get_i(self):
        """
        Return the `i` field of this ManinSymbol `(i,u,v)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s._ManinSymbol__get_i()
            2
            sage: s.i
            2
        """
        return self.__t[0]

    i = property(__get_i)
    def __get_u(self):
        """
        Return the `u` field of this ManinSymbol `(i,u,v)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.u # indirect doctest
            2
        """
        return self.__t[1]
    u = property(__get_u)

    def __get_v(self):
        """
        Return the `v` field of this ManinSymbol `(i,u,v)`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.v # indirect doctest
            3
        """
        return self.__t[2]
    v = property(__get_v)

    def _repr_(self):
        """
        Returns a string representation of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: str(s)  # indirect doctest
            '[X^2*Y^4,(2,3)]'
        """
        if self.weight() > 2:
            polypart = _print_polypart(self.i, self.weight()-2-self.i)
            return "[%s,(%s,%s)]"%\
                   (polypart, self.u, self.v)
        return "(%s,%s)"%(self.u, self.v)

    def _latex_(self):
        """
        Returns a LaTeX representation of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: latex(s) # indirect doctest
            [X^2*Y^4,(2,3)]
        """
        return self._repr_()

    def __cmp__(self, other):
        """
        Comparison function for ManinSymbols.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: slist = m.manin_symbol_list()
            sage: cmp(slist[10],slist[20])
            -1
            sage: cmp(slist[20],slist[10])
            1
            sage: cmp(slist[20],slist[20])
            0
        """
        if not isinstance(other, ManinSymbol):
            return -1
        if self.__t == other.__t:
            return 0
        return cmp(self.tuple(), other.tuple())

    def __mul__(self, matrix):
        """
        Returns the result of applying a matrix to this ManinSymbol.


        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,2)
            sage: s = ManinSymbol(m,(0,2,3))
            sage: s*[1,2,0,1]
            (2,7)

        ::

            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s*[1,2,0,1]
            Traceback (most recent call last):
            ...
            NotImplementedError: ModSym * Matrix only implemented in weight 2
        """
        if self.weight() > 2:
            raise NotImplementedError("ModSym * Matrix only implemented "
                                      "in weight 2")
        from sage.matrix.matrix import is_Matrix
        if is_Matrix(matrix):
            if (not matrix.nrows() == 2) or (not matrix.ncols() == 2):
                raise ValueError("matrix(=%s) must be 2x2" % matrix)
            matrix = matrix.list()
        return ManinSymbol(self.parent(), \
                           (self.i,
                           matrix[0]*self.u + matrix[2]*self.v,\
                           matrix[1]*self.u + matrix[3]*self.v))
        raise ArithmeticError("Multiplication of %s by %s not defined." % (self, matrix))


    def apply(self, a,b,c,d):
        """
        Return the image of self under the matrix `[a,b;c,d]`.

        Not implemented for raw ManinSymbol objects, only for members
        of ManinSymbolLists.

        EXAMPLE::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,2)
            sage: m.apply(10,[1,0,0,1]) # not implemented for base class
        """
        raise NotImplementedError

    def __copy__(self):
        """
        Return a copy of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s2 = copy(s)
            sage: s2
            [X^2*Y^4,(2,3)]

        """
        return ManinSymbol(self.parent(), (self.i, self.u, self.v))

    def lift_to_sl2z(self, N=None):
        r"""
        Returns a lift of this Manin Symbol to `SL_2(\mathbb{Z})`.

        If this Manin symbol is `(c,d)` and `N` is its level, this
        function returns a list `[a,b, c',d']` that defines a 2x2
        matrix with determinant 1 and integer entries, such that
        `c=c'` (mod `N`) and `d=d'` (mod `N`).

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s
            [X^2*Y^4,(2,3)]
            sage: s.lift_to_sl2z()
            [1, 1, 2, 3]

        """
        if N is None:
            N = self.level()
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
        Returns cusps `alpha`, `beta` such that this Manin symbol, viewed as a
        symbol for level `N`, is `X^i*Y^{k-2-i} \{alpha, beta\}`.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3)); s
            [X^2*Y^4,(2,3)]
            sage: s.endpoints()
            (1/3, 1/2)
        """
        if N is None:
            N = self.parent().level()
        else:
            N=int(N)
            if N < 1:
                raise ArithmeticError, "N must be positive"
        a,b,c,d = self.lift_to_sl2z()
        return cusps.Cusp(b,d), cusps.Cusp(a,c)

    def parent(self):
        """
        Return the parent of this ManinSymbol.


        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.parent()
            Manin Symbol List of weight 8 for Gamma0(5)
        """
        return self.__parent

    def weight(self):
        """
        Return the weight of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.weight()
            8

        """
        return self.__parent.weight()

    def level(self):
        """
        Return the level of this ManinSymbol.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.level()
            5

        """
        return self.__parent.level()

    def modular_symbol_rep(self):
        """
        Returns a representation of self as a formal sum of modular
        symbols.

        The result is not cached.

        EXAMPLES::

            sage: from sage.modular.modsym.manin_symbols import ManinSymbol, ManinSymbolList_gamma0
            sage: m = ManinSymbolList_gamma0(5,8)
            sage: s = ManinSymbol(m,(2,2,3))
            sage: s.modular_symbol_rep()
             144*X^6*{1/3, 1/2} - 384*X^5*Y*{1/3, 1/2} + 424*X^4*Y^2*{1/3, 1/2} - 248*X^3*Y^3*{1/3, 1/2} + 81*X^2*Y^4*{1/3, 1/2} - 14*X*Y^5*{1/3, 1/2} + Y^6*{1/3, 1/2}


        """
        # TODO: It would likely be much better to do this slightly more directly
        from sage.modular.modsym.modular_symbols import ModularSymbol
        x = ModularSymbol(self.__parent, self.i, 0, rings.infinity)
        a,b,c,d = self.lift_to_sl2z()
        return x.apply([a,b,c,d])

def _print_polypart(i, j):
    r"""
    Helper function for printing the polynomial part `X^iY^j` of a ManinSymbol.

    EXAMPLES::

    sage: from sage.modular.modsym.manin_symbols import _print_polypart
    sage: _print_polypart(2,3)
    'X^2*Y^3'
    sage: _print_polypart(2,0)
    'X^2'
    sage: _print_polypart(0,1)
    'Y'
    """
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
