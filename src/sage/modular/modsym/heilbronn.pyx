"""
Heilbronn matrix computation
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

import sage.rings.arith

import sage.misc.misc

include '../../ext/cdefs.pxi'
include '../../ext/interrupt.pxi'
include '../../ext/stdsage.pxi'

cimport p1list
import  p1list
cdef p1list.export export
export = p1list.export()

ctypedef long long llong

cdef int llong_prod_mod(int a, int b, int N):
    cdef int c
    c = <int> ( ((<llong> a) * (<llong> b)) % (<llong> N) )
    if c < 0:
        c = c + N
    return c

def HeilbronnCremonaList(int p):
    """
    The Heilbronn matrices of determinant p, as defined by Cremona.
    """
    cdef int r, x1, x2, y1, y2, a, b, c, x3, y3, q
    cdef float f
    cdef object ans

    if p == 2:
        return [[1,0,0,2], [2,0,0,1], [2,1,0,1], [1,0,1,2]]

    #if cremona.has_key(p):
    #    return cremona[p]
    assert sage.rings.arith.is_prime(p), "Input must be a pseudoprime."

    ans = [[1,0,0,p]]
    # WARNING: In C (which is what we're writing in below!), -p/2 means
    # "round toward 0", but in Python it means "round down" (right now),
    # and in Python 3.0 it'll mean "make a float".
    _sig_on
    for r from -p/2 <= r < p/2+1:
        x1=p; x2=-r; y1=0; y2=1; a=-p; b=r; c=0; x3=0; y3=0; q=0
        ans.append([x1,x2,y1,y2])
        while b:
            f = roundf(<float>a / <float> b)
            q = <int> f
            c = a - b*q
            a = -b
            b = c
            x3 = q*x2 - x1
            x1 = x2; x2 = x3; y3 = q*y2 - y1; y1 = y2; y2 = y3
            ans.append([x1,x2, y1,y2])
    _sig_off
    return ans

def HeilbronnMerelList(int n):
    """
    Set of Heilbronn matrices of determinant n, as defined by Merel.
    """
    cdef int a, q, d, b, c, bc
    cdef object H

    H = []
    _sig_on
    for a from 1 <= a <= n:
        ## We have ad-bc=n so c=0 and ad=n, or b=(ad-n)/c
        ## Must have ad - n >= 0, so d must be >= Ceiling(n/a).
        q = n/a
        if q*a == n:
          d = q
          for b from 0 <= b < a:
              H.append([a,b,0,d])
          for c from 1 <= c < d:
              H.append([a,0,c,d])
        for d from q+1 <= d <= n:
            bc = a*d-n
            ## Divisor c of bc must satisfy Floor(bc/c) lt a and c lt d.
            ## c ge (bc div a + 1)  <=>  Floor(bc/c) lt a  (for integers)
            ## c le d - 1           <=>  c lt d
            for c from bc/a + 1 <= c < d:
                if bc % c == 0:
                    H.append([a,bc/c,c,d])
    _sig_off
    return H

cdef struct list:
    int *v
    int i   # how many positions of list are filled
    int n   # how much memory has been allocated

cdef int* expand(int *v, int n, int new_length) except NULL:
    cdef int *w, i
    w = <int*>  sage_malloc(new_length*sizeof(int))
    if w == <int*> 0:
        return NULL
    if v:
        for i from 0 <= i < n:
            w[i] = v[i]
        sage_free(v)
    return w

cdef int list_append(list* L, int a) except -1:
    cdef int j
    if L.i >= L.n:
        j = 10 + 2*L.n
        L.v = expand(L.v, L.n, j)
        L.n = j
    L.v[L.i] = a
    L.i = L.i + 1

cdef int list_append4(list* L, int a, int b, int c, int d) except -1:
    list_append(L, a)
    list_append(L, b)
    list_append(L, c)
    list_append(L, d)

cdef void list_clear(list L):
    sage_free(L.v)

cdef void list_init(list* L):
    L.n = 16
    L.i = 0
    L.v = expand(<int*>0, 0, L.n)


cdef class Heilbronn:
    cdef int length
    cdef list list

    def __dealloc__(self):
        list_clear(self.list)

    def _initialize_list(self):
        raise NotImplementedError

    def __getitem__(self, int n):
        if n < 0 or n >= self.length:
            raise IndexError
        return [self.list.v[4*n], self.list.v[4*n+1], \
                self.list.v[4*n+2], self.list.v[4*n+3]]

    def __len__(self):
        return self.length

    def to_list(self):
        cdef int i
        L = []
        for i from 0 <= i < self.length:
            L.append([self.list.v[4*i], self.list.v[4*i+1], \
                      self.list.v[4*i+2], self.list.v[4*i+3]])
        return L

    def apply(self, int u, int v, int N):
        """
        Return a list of pairs ((c,d),m), which is obtained as follows:
          1) Compute the images (a,b) of the vector (u,v) (mod N) acted on by
             each of the HeilbronnCremona matrices in self.
          2) Reduce each (a,b) to canonical form (c,d) using p1normalize
          3) Sort.
          4) Create the list ((c,d),m), where m is the number of times
             that (c,d) appears in the list created in steps 1--3 above.
        Note that the pairs ((c,d),m) are sorted lexicographically by (c,d).
        """
        cdef int i, a, b, c, d, s
        cdef object X
        M = {}
        t = sage.misc.misc.verbose("start making list M.",level=5)
        _sig_on
        if N < 32768:   # use ints with no reduction modulo N
            for i from 0 <= i < self.length:
                a = u*self.list.v[4*i] + v*self.list.v[4*i+2]
                b = u*self.list.v[4*i+1] + v*self.list.v[4*i+3]
                export.c_p1_normalize_int(N, a, b, &c, &d, &s, 0)
                X = (c,d)
                if M.has_key(X):
                    M[X] = M[X] + 1
                else:
                    M[X] = 1
        elif N < 46340:    # use ints but reduce mod N so can add two
            for i from 0 <= i < self.length:
                a = (u * self.list.v[4*i])%N + (v * self.list.v[4*i+2])%N
                b = (u * self.list.v[4*i+1])%N + (v * self.list.v[4*i+3])%N
                export.c_p1_normalize_int(N, a, b, &c, &d, &s, 0)
                X = (c,d)
                if M.has_key(X):
                    M[X] = M[X] + 1
                else:
                    M[X] = 1
        else:
            for i from 0 <= i < self.length:
                a = llong_prod_mod(u,self.list.v[4*i],N) + llong_prod_mod(v,self.list.v[4*i+2], N)
                b = llong_prod_mod(u,self.list.v[4*i+1],N) + llong_prod_mod(v,self.list.v[4*i+3], N)
                export.c_p1_normalize_llong(N, a, b, &c, &d, &s, 0)
                X = (c,d)
                if M.has_key(X):
                    M[X] = M[X] + 1
                else:
                    M[X] = 1
        t = sage.misc.misc.verbose("finished making list M.",t, level=5)
        mul = []
        for x,y in M.items():
            mul.append((x,y))
        t = sage.misc.misc.verbose("finished making mul list.",t, level=5)
        _sig_off
        return mul

cdef class HeilbronnCremona(Heilbronn):
    cdef public int p

    def __init__(self, int p):
        if p <= 1 or not sage.rings.arith.is_prime(p):
            raise ValueError, "p must be >= 2 and prime"
        self.p = p
        self._initialize_list()

    def __repr__(self):
        return "The Cremona-Heilbronn matrices of determinant %s"%self.p

    def _initialize_list(self):
        cdef int r, x1, x2, y1, y2, a, b, c, x3, y3, q, n, p
        cdef float f
        cdef list *L
        list_init(&self.list)
        L = &self.list
        p = self.p

        list_append4(L, 1,0,0,p)

        # When p==2, then Heilbronn matrices are
        #    [[1,0,0,2], [2,0,0,1], [2,1,0,1], [1,0,1,2]]
        # which are not given by the algorithm below.
        if p == 2:
            list_append4(L, 2,0,0,1)
            list_append4(L, 2,1,0,1)
            list_append4(L, 1,0,1,2)
            self.length = 4
            return

        # NOTE: In C, -p/2 means "round toward 0", but in Python it
        # means "round down."
        _sig_on
        for r from -p/2 <= r < p/2+1:
            x1=p; x2=-r; y1=0; y2=1; a=-p; b=r; c=0; x3=0; y3=0; q=0
            list_append4(L, x1,x2,y1,y2)
            while b:
                f = roundf(<float>a / <float> b)
                q = <int> f
                c = a - b*q
                a = -b
                b = c
                x3 = q*x2 - x1
                x1 = x2; x2 = x3; y3 = q*y2 - y1; y1 = y2; y2 = y3
                list_append4(L, x1,x2, y1,y2)
        self.length = L.i/4
        _sig_off

    def __getitem__(self, int n):
        if n < 0 or n >= self.length:
            raise IndexError
        return [self.list.v[4*n], self.list.v[4*n+1], \
                self.list.v[4*n+2], self.list.v[4*n+3]]



cdef class HeilbronnMerel(Heilbronn):
    cdef public int n

    def __init__(self, int n):
        if n <= 0:
            raise ValueError, "n (=%s) must be >= 1"%n
        self.n = n
        self._initialize_list()

    def __repr__(self):
        return "The Merel-Heilbronn matrices of determinant %s"%self.n

    def _initialize_list(self):
        cdef int a, q, d, b, c, bc, n
        cdef list *L
        list_init(&self.list)
        L = &self.list
        n = self.n

        _sig_on
        for a from 1 <= a <= n:
            ## We have ad-bc=n so c=0 and ad=n, or b=(ad-n)/c
            ## Must have ad - n >= 0, so d must be >= Ceiling(n/a).
            q = n/a
            if q*a == n:
                d = q
                for b from 0 <= b < a:
                    list_append4(L, a,b,0,d)
                for c from 1 <= c < d:
                    list_append4(L, a,0,c,d)
            for d from q+1 <= d <= n:
                bc = a*d-n
                ## Divisor c of bc must satisfy Floor(bc/c) lt a and c lt d.
                ## c ge (bc div a + 1)  <=>  Floor(bc/c) lt a  (for integers)
                ## c le d - 1           <=>  c lt d
                for c from bc/a + 1 <= c < d:
                    if bc % c == 0:
                        list_append4(L,a,bc/c,c,d)
        self.length = L.i/4
        _sig_off
