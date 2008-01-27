"""
Permutation group elements

AUTHORS:
    - David Joyner (2006-02)
    - David Joyner (2006-03), word problem method and reorganization
    - Robert Bradshaw (2007-11), convert to Cython

EXAMPLES:
The Rubik's cube group:
    sage: f= [(17,19,24,22),(18,21,23,20),(6,25,43,16),(7,28,42,13),(8,30,41,11)]
    sage: b=[(33,35,40,38),(34,37,39,36),( 3, 9,46,32),( 2,12,47,29),( 1,14,48,27)]
    sage: l=[( 9,11,16,14),(10,13,15,12),( 1,17,41,40),( 4,20,44,37),( 6,22,46,35)]
    sage: r=[(25,27,32,30),(26,29,31,28),( 3,38,43,19),( 5,36,45,21),( 8,33,48,24)]
    sage: u=[( 1, 3, 8, 6),( 2, 5, 7, 4),( 9,33,25,17),(10,34,26,18),(11,35,27,19)]
    sage: d=[(41,43,48,46),(42,45,47,44),(14,22,30,38),(15,23,31,39),(16,24,32,40)]
    sage: cube = PermutationGroup([f,b,l,r,u,d])
    sage: F=cube.gens()[0]
    sage: B=cube.gens()[1]
    sage: L=cube.gens()[2]
    sage: R=cube.gens()[3]
    sage: U=cube.gens()[4]
    sage: D=cube.gens()[5]
    sage: cube.order()
    43252003274489856000
    sage: F.order()
    4

The interested user may wish to explore the following commands:
move = cube.random_element() and time word_problem([F,B,L,R,U,D], move, False).
This typically takes about 5 minutes (on a 2 Ghz machine) and outputs
a word ('solving' the cube in the position move) with about 60 terms
or so.

OTHER EXAMPLES:
We create element of a permutation group of large degree.
    sage: G = SymmetricGroup(30)
    sage: s = G(srange(30,0,-1)); s
    (1,30)(2,29)(3,28)(4,27)(5,26)(6,25)(7,24)(8,23)(9,22)(10,21)(11,20)(12,19)(13,18)(14,17)(15,16)
"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


import random

import sage.groups.group as group

include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/python_list.pxi"

from sage.rings.all      import ZZ, Integer, is_MPolynomial, MPolynomialRing, is_Polynomial
from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement

import sage.structure.coerce as coerce

import operator

from sage.rings.integer import Integer

from sage.ext.arith cimport arith_llong
cdef arith_llong arith = arith_llong()
cdef extern from *:
    long long LONG_LONG_MAX

#import permgroup_named

def make_permgroup_element(G, x):
    G._deg = len(x)
    return G(x, check=False)

def is_PermutationGroupElement(x):
    return isinstance(x, PermutationGroupElement)


def gap_format(x):
    """
    Put a permutation in Gap format, as a string.
    """
    if isinstance(x, list) and not isinstance(x[0], tuple):
        return gap.eval('PermList(%s)' % x)
    x = str(x).replace(' ','').replace('\n','')
    return x.replace('),(',')(').replace('[','').replace(']','')

cdef class PermutationGroupElement(MultiplicativeGroupElement):
    """
    An element of a permutation group.

    EXAMPLES:
        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5)]
        sage: g = G.random_element()
        sage: g in G
        True
        sage: g = G.gen(0); g
        (1,2,3)(4,5)
        sage: print g
        (1,2,3)(4,5)
        sage: g*g
        (1,3,2)
        sage: g**(-1)
        (1,3,2)(4,5)
        sage: g**2
        (1,3,2)
        sage: G = PermutationGroup([(1,2,3)])
        sage: g = G.gen(0); g
        (1,2,3)
        sage: g.order()
        3

    This example illustrates how permutations act on multivariate
    polynomials.

        sage: R = MPolynomialRing(RationalField(), 5, ["x","y","z","u","v"])
        sage: x, y, z, u, v = R.gens()
        sage: f = x**2 - y**2 + 3*z**2
        sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
        sage: sigma = G.gen(0)
        sage: f * sigma
        3*x^2 + y^2 - z^2

    """
    def __init__(self, g, parent = None, check = True):
        r"""
        Create element of a permutation group.

        There are several ways to define a permutation group element:
        \begin{itemize}
        \item  Define a permutation group $G$, then use
               \code{G.gens()} and multiplication * to construct
               elements.
        \item Define a permutation group $G$, then use e.g.,
               \code{G([(1,2),(3,4,5)])} to construct an element of
               the group.  You could also use \code{G('(1,2)(3,4,5)')}
        \item Use e.g., \code{PermutationGroupElement([(1,2),(3,4,5)])}
        or \code{PermutationGroupElement('(1,2)(3,4,5)')}
              to make a permutation group element with parent $S_5$.
        \end{itemize}

        INPUT:
            g -- defines element
            parent (optional) -- defines parent group (g must be in parent if
                                 specified, or a TypeError is raised).
            check -- bool (default: True), if False assumes g is a
                     gap element in parent (if specified).

        EXAMPLES:
        We illustrate construction of permutation using several
        different methods.

        First we construct elements by multiplying together generators
        for a group.

            sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'])
            sage: s = G.gens()
            sage: s[0]
            (1,2)(3,4)
            sage: s[1]
            (3,4,5,6)
            sage: s[0]*s[1]
            (1,2)(3,5,6)
            sage: (s[0]*s[1]).parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]

        Next we illustrate creation of a permutation using
        coercion into an already-created group.

            sage: g = G([(1,2),(3,5,6)])
            sage: g
            (1,2)(3,5,6)
            sage: g.parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]
            sage: g == s[0]*s[1]
            True

        We can also use a string instead of a list to specify
        the permutation.

            sage: h = G('(1,2)(3,5,6)')
            sage: g == h
            True

        We can also make a permutation group element directly
        using the \code{PermutationGroupElement} command.  Note
        that the parent is then the full symmetric group $S_n$,
        where $n$ is the largest integer that is moved by the
        permutation.

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k
            (1,2)(3,5,6)
            sage: k.parent()
            Symmetric group of order 6! as a permutation group

        Note the comparison of permutations doesn't require that the
        parent groups are the same.

            sage: k == g
            True

        Arithmetic with permutations having different parents is also defined:

            sage: k*g
            (3,6,5)
            sage: (k*g).parent()
            Symmetric group of order 6! as a permutation group

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: loads(dumps(G.0)) == G.0
            True

        EXAMPLES:
            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k._gap_()
            (1,2)(3,5,6)
            sage: k._gap_().parent()
            Gap
        """
        import sage.groups.perm_gps.permgroup_named
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.groups.perm_gps.permgroup import PermutationGroup_generic
        if check:
            if not (parent is None or isinstance(parent, PermutationGroup_generic)):
                raise TypeError, 'parent must be a permutation group'
            self.__gap = gap_format(g)
            if not parent is None:
                P = parent._gap_()
                if not P.parent()(self.__gap) in P:
                    raise TypeError, 'permutation %s not in %s'%(self.__gap, parent)
        else:
            self.__gap = str(g)

        if parent is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
            parent = sage.groups.perm_gps.permgroup_named.SymmetricGroup(G(self.__gap).LargestMovedPoint())

        Element.__init__(self, parent)

        self.n = parent.degree()
        if self.perm is NULL:
            self.perm = <int *>sage_malloc(sizeof(int) * self.n)
        else:
            self.perm = <int *>sage_realloc(self.perm, sizeof(int) * self.n)

        if PyList_CheckExact(g) and not PY_TYPE_CHECK(g[0], tuple):
            v = g
        else:
            v = eval(gap.eval('ListPerm(%s)' % self.__gap))

        cdef int i
        assert(len(v) <= self.n)
        for i from 0 <= i < len(v):
            self.perm[i] = v[i] - 1
        for i from len(v) <= i < self.n:
            self.perm[i] = i


    def __cinit__(self, g = None, parent = None, check = True):
        self.perm = NULL

    def __dealloc__(self):
        sage_free(self.perm)

    def __reduce__(self):
        return make_permgroup_element, (self._parent, self.list())

    cdef PermutationGroupElement _new_c(self):
        cdef PermutationGroupElement other = PY_NEW_SAME_TYPE(self)
        if HAS_DICTIONARY(self):
            other.__class__ = self.__class__
        other._parent = self._parent
        other.n = self.n
        other.perm = <int *>sage_malloc(sizeof(int) * other.n)
        return other

    def _gap_(self, G=None):
        if self._gap_element is None or \
            (G is not None and self._gap_element._parent is not G):
            if G is None:
                import sage.interfaces.gap
                G = sage.interfaces.gap.gap
            self._gap_element = G("PermList(%s)" % self.list())
        return self._gap_element

    def _repr_(self):
        """
        Return string representation of this permutation.

        EXAMPLES:
        We create the permutation $(1,2,3)(4,5)$ and print it.

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._repr_()
            '(1,2,3)(4,5)'

        Permutation group elements support renaming them so
        they print however you want, as illustred below:

            sage: g.rename('sigma')
            sage: g
            sigma
            sage: g.rename()
            sage: g
            (1,2,3)(4,5)
        """
        cycles = self.cycle_tuples()
        if len(cycles) == 0:
            return '()'
        return ''.join([str(c) for c in cycles]).replace(', ',',')

    def _latex_(self):
        return str(self)

    def __getitem__(self, i):
        """
        Return the ith permutation cycle in the disjoint cycle
        representation of self.

        INPUT:
            i -- integer

        OUTPUT:
            a permutation group element

        EXAMPLE:
            sage: G = PermutationGroup([[(1,2,3),(4,5)]],5)
            sage: g = G.gen(0)
            sage: g[0]
            (1,2,3)
            sage: g[1]
            (4,5)
        """
        return self.cycles()[i]

    def __cmp__(PermutationGroupElement self, PermutationGroupElement right):
        """
        Compare group elements self and right.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G.gen(0) < G.gen(1)
            False
            sage: G.gen(0) > G.gen(1)
            True
        """
        cdef int i
        cdef bint equal = 1
        for i from 0 <= i < self.n:
            if self.perm[i] != right.perm[i]:
                equal = 0
                break
        if equal:
            return 0
        r = right._gap_()
        G = r.parent()
        l = self._gap_(G)
        return cmp(l,r)

    def __call__(self, i):
        """
        Returns the image of the integer i under this permutation.
        Alternately, if i is a list, tuple or string, returns the
        result of self acting on i.

        EXAMPLE:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G
            Permutation Group with generators [(1,2,3)(4,5)]
            sage: g = G.gen(0)
            sage: g(5)
            4
            sage: g('abcde')
            'bcaed'
            sage: g([0,1,2,3,4])
            [1, 2, 0, 4, 3]
            sage: g(('who','what','when','where','why'))
            ('what', 'when', 'who', 'why', 'where')

            sage: g(x)
            Traceback (most recent call last):
            ...
            ValueError: Must be an integer, list, tuple or string.
            sage: g(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Must be an integer, list, tuple or string.

        """
        cdef int j
        if isinstance(i,(list,tuple,str)):
            permuted = [i[self.perm[j]] for j from 0 <= j < self.n]
            if PY_TYPE_CHECK(i, tuple):
                permuted = tuple(permuted)
            elif PY_TYPE_CHECK(i, str):
                permuted = ''.join(permuted)
            permuted += i[self.n:]
            return permuted
        else:
            if not isinstance(i, (int, Integer)):
                raise ValueError("Must be an integer, list, tuple or string.")
            j = i
            if 1 <= j <= self.n:
                return self.perm[j-1]+1
            else:
                return i

    def _r_action(self, left):
        """
        Return the right action of self on left.

        For example, if f=left is a polynomial, then this function
        returns f(sigma*x), which is image of f under the right action
        of sigma on the indeterminates.  This is a right action since
        the image of f(sigma*x) under tau is f(sigma*tau*x).

        INPUT:
            left -- element of space on which permutations act from the right

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: R.<x,y,z,u,v> = MPolynomialRing(QQ,5)
            sage: f = x^2 + y^2 - z^2 + 2*u^2
            sage: sigma, tau = G.gens()
            sage: f*sigma
            -x^2 + y^2 + z^2 + 2*v^2
            sage: f*tau
            y^2 + z^2 - u^2 + 2*v^2
            sage: f*(sigma*tau)
            2*x^2 - y^2 + z^2 + u^2
            sage: (f*sigma)*tau
            2*x^2 - y^2 + z^2 + u^2
        """
        if is_Polynomial(left):
            if self != 1:
                raise ValueError, "%s does not act on %s"%(self, left.parent())
            return left
        elif isinstance(left, PermutationGroupElement):
            return PermutationGroupElement(self._gap_()*left._gap_(),
                                           parent = None, check = True)
        elif is_MPolynomial(left):
            F = left.base_ring()
            R = left.parent()
            x = R.gens()
            vars = list(x)
            try:
                sigma_x  = [vars[int(self(i+1)-1)] for i in range(len(x))]
            except IndexError:
                raise TypeError, "%s does not act on %s"%(self, left.parent())
            return left(tuple(sigma_x))
        else:
            raise TypeError, "left (=%s) must be a polynomial."%left


    cdef MonoidElement _mul_c_impl(left, MonoidElement _right):
        cdef PermutationGroupElement prod = left._new_c()
        cdef PermutationGroupElement right = <PermutationGroupElement>_right
        cdef int i
        for i from 0 <= i < left.n:
            prod.perm[i] = right.perm[left.perm[i]]
        return prod

    def __invert__(self):
        """
        Return the inverse of this permutation.

        EXAMPLES:
            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: ~g
            (1,3,2)(4,5)
            sage: (~g) * g
            ()
        """
        cdef PermutationGroupElement inv = self._new_c()
        cdef int i
        for i from 0 <= i < self.n:
            inv.perm[self.perm[i]] = i
        return inv

    cpdef list(self):
        """
        Returns list of the images of the integers from 1 to n under
        this permutation as a list of Python ints.

        EXAMPLES:
            sage: G = SymmetricGroup(4)
            sage: x = G([2,1,4,3]); x
            (1,2)(3,4)
            sage: v = x.list(); v
            [2, 1, 4, 3]
            sage: type(v[0])
            <type 'int'>
            sage: x = G([2,1]); x
            (1,2)
            sage: x.list()
            [2, 1, 3, 4]
        """
        cdef int i
        return [self.perm[i]+1 for i from 0 <= i < self.n]

    def __hash__(self):
        """
        Return hash of this permutation.

        EXAMPLES:
            sage: G = SymmetricGroup(5)
            sage: s = G([2,1,5,3,4])
            sage: s.tuple()
            (2, 1, 5, 3, 4)
            sage: hash(s)
            1592966088          # 32-bit
            2865702456085625800 # 64-bit
            sage: hash(s.tuple())
            1592966088          # 32-bit
            2865702456085625800 # 64-bit
        """
        return hash(self.tuple())

    def tuple(self):
        """
        Return tuple of images of integers under self.

        EXAMPLES:
            sage: G = SymmetricGroup(5)
            sage: s = G([2,1,5,3,4])
            sage: s.tuple()
            (2, 1, 5, 3, 4)
        """
        if self.__tuple is None:
            self.__tuple = tuple(self.list())
        return self.__tuple

    def dict(self):
        """
        Returns list of the images of the integers from 1 to n under
        this permutation as a list of Python ints.

        NOTE:
            self.list() returns a zero-indexed list.  self.dict() return
            the actual mapping of the permutation, which will be indexed
            starting with 1.

        EXAMPLES:
            sage: G = SymmetricGroup(4)
            sage: g = G((1,2,3,4)); g
            (1,2,3,4)
            sage: v = g.dict(); v
            {1: 2, 2: 3, 3: 4, 4: 1}
            sage: type(v[1])
            <type 'int'>
            sage: x = G([2,1]); x
            (1,2)
            sage: x.dict()
            {1: 2, 2: 1, 3: 3, 4: 4}
        """
        cdef int i
        u = {}
        for i from 0 <= i < self.n:
            u[i+1] = self.perm[i]+1
        return u

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer $n$ for which $g^n = 1$.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.order()
            6

        TESTS:
            sage: prod(primes(150))
            1492182350939279320058875736615841068547583863326864530410
            sage: L = [tuple(range(sum(primes(p))+1, sum(primes(p))+1+p)) for p in primes(150)]
            sage: PermutationGroupElement(L).order()
            1492182350939279320058875736615841068547583863326864530410
        """
        order = None
        cdef long long order_c = 1
        cdef int cycle_len
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            k = self.perm[i]
            cycle_len = 1
            while k != i:
                seen[k] = 1
                k = self.perm[k]
                cycle_len += 1
            if order is not None:
                order = order.lcm(cycle_len)
            else:
                order_c = (order_c * cycle_len) / arith.c_gcd_longlong(order_c, cycle_len)
                if order_c > LONG_LONG_MAX / (self.n - i):
                    order = Integer(order_c)
        sage_free(seen)
        return int(order_c) if order is None else order

    def sign(self):
        """
        Returns the sign of self, which is $(-1)^{s}$, where $s$ is
        the number of swaps.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.sign()
            -1

        ALGORITHM:
            Only even cycles contribute to the sign, thus
            $$sign(sigma) = (-1)^{\sum_c len(c)-1}$$
            where the sum is over cycles in self.
        """
        cdef int cycle_len_sum = 0
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            k = self.perm[i]
            while k != i:
                seen[k] = 1
                k = self.perm[k]
                cycle_len_sum += 1
        sage_free(seen)
        return 1 - 2*(cycle_len_sum % 2) # == (-1)^cycle_len


    def orbit(self, n, bint sorted=True):
        """
        Returns the orbit of the integer $n$ under this group element,
        as a sorted list of integers.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.orbit(4)
            [4, 5]
            sage: g.orbit(3)
            [1, 2, 3]
            sage: g.orbit(10)
            [10]
        """
        cdef int i = n
        cdef int start = i
        if 1 <= i <= self.n:
            L = [i]
            i = self.perm[i-1]+1
            while i != start:
                PyList_Append(L,i)
                i = self.perm[i-1]+1
            if sorted:
                L.sort()
            return L
        else:
            return [n]

    def cycles(self):
        """
        Return self as a list of disjoint cycles.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5,6,7)'])
            sage: g = G.0
            sage: g.cycles()
            [(1,2,3), (4,5,6,7)]
            sage: a, b = g.cycles()
            sage: a(1), b(1)
            (2, 1)
        """
        L = []
        cdef PermutationGroupElement cycle
        cdef int i, j, k, next_k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            cycle = self._new_c()
            for j from 0 <= j < self.n: cycle.perm[j] = j
            k = cycle.perm[i] = self.perm[i]
            while k != i:
                seen[k] = 1
                next_k = cycle.perm[k] = self.perm[k]
                k = next_k
            PyList_Append(L, cycle)
        sage_free(seen)
        return L

    def cycle_tuples(self):
        """
        Return self as a list of disjoint cycles, represented
        as tuples rather than permutation group elements.
        """
        L = []
        cdef int i, k
        cdef bint* seen = <bint *>sage_malloc(sizeof(bint) * self.n)
        for i from 0 <= i < self.n: seen[i] = 0
        for i from 0 <= i < self.n:
            if seen[i] or self.perm[i] == i:
                continue
            cycle = [i+1]
            k = self.perm[i]
            while k != i:
                PyList_Append(cycle, k+1)
                seen[k] = 1
                k = self.perm[k]
            PyList_Append(L, tuple(cycle))
        sage_free(seen)
        return L

    def matrix(self):
        """
        Returns deg x deg permutation matrix associated
        to the permutation self

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.matrix()
            [0 1 0 0 0]
            [0 0 1 0 0]
            [1 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 1 0]
        """
        M = MatrixSpace(ZZ, self.n, self.n, sparse=True)
        cdef int i
        entries = {}
        for i from 0 <= i < self.n:
            entries[i, self.perm[i]] = 1
        return M(entries)

    def word_problem(g, words, display=True):
        """
        G and H are permutation groups, g in G, H is a subgroup of G
        generated by a list (words) of elements of G. If g is in H,
        return the expression for g as a word in the elements of
        (words).

        This function does not solve the word problem in SAGE. Rather
        it pushes it over to GAP, which has optimized algorithms for
        the word problem. Essentially, this function is a wrapper for the GAP
        functions "EpimorphismFromFreeGroup" and "PreImagesRepresentative".

        EXAMPLE:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: g1 = G.gens()[0]
            sage: g2 = G.gens()[1]
            sage: h = g1^2*g2*g1
            sage: h.word_problem([g1,g2], False)
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')
            sage: h.word_problem([g1,g2])
               x1^2*x2^-1*x1
               [['(1,2,3)(4,5)', 2], ['(3,4)', -1], ['(1,2,3)(4,5)', 1]]
            ('x1^2*x2^-1*x1', '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)')
        """
        import copy
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap

        G = gap(words[0].parent())
        g = words[0].parent()(g)
        gensH = gap(words)
        H = gensH.Group()
        hom = G.EpimorphismFromFreeGroup()
        ans = hom.PreImagesRepresentative(gap(g))

        l1 = str(ans)
        l2 = copy.copy(l1)
        l4 = []
        l3 = l1.split("*")
        for i in range(1,len(words)+1):
            l2 = l2.replace("x"+str(i),str(words[i-1]))

        if display:
            for i in range(len(l3)):    ## parsing the word for display
                if len(l3[i].split("^"))==2:
                    l4.append([l3[i].split("^")[0],int(l3[i].split("^")[1])])
                if len(l3[i].split("^"))==1:
                    l4.append([l3[i].split("^")[0],1])
            l5 = copy.copy(l4)
            for i in range(len(l4)):
                for j in range(1,len(words)+1):
                    l5[i][0] = l5[i][0].replace("x"+str(j),str(words[j-1]))
            print "         ",l1
            print "         ",l5
        return l1,l2
