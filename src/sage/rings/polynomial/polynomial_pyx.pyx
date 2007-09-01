"""
Polynomial over QQ

NOTE -- this is no longer used ?!

AUTHOR:
    -- William Stein (2004): first version
    -- Didier Deshom <dfdeshom@gmail.com> (2006-03-01): optimization to _pow.

todo: if this file is used, __pow__ should be modified to call generic_power.
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

START_PRIME = 40009  # used for modular algorithms

import sage.rings.polynomial_element as polynomial

include "../../ext/cdefs.pxi"
include "../../ext/gmp.pxi"
include "../../ext/stdsage.pxi"

cimport sage.structure.element
import  sage.structure.element

cimport sage.ext.arith_gmp
import  sage.ext.arith_gmp
cdef sage.ext.arith_gmp.functions ag
ag = sage.ext.arith_gmp.functions()

ctypedef unsigned long int ulong
cdef int py_mpq_set(mpq_t x, y) except -1:
    s = str(y)
    mpq_set_str(x, s, 0)
    mpq_canonicalize(x)


###########################################################################
# Polynomials over the integers modulo an int.
###########################################################################
# Representation: We represent polynomials using a C array of deg+1
# unsigned int's and a degree

cdef struct Pmodint_data:
    int degree
    ulong p
    ulong *v

ctypedef Pmodint_data Pmodint

cdef int Pmodint_init(Pmodint* f, ulong p, int degree) except -1:
    cdef int i
    if degree < -1:
        print "degree = ", degree
        raise ArithmeticError, "The degree must be >= -1."
    f.degree = degree
    f.p = p
    if degree == -1:
        f.v = <ulong*> 0
        return 0
    f.v = <ulong*>sage_malloc(sizeof(ulong)*(degree+1))
    if f.v == <ulong*>0:
        raise MemoryError, "Error allocating memory for polynomial."
    return 0

cdef int Pmodint_clear(Pmodint* f) except -1:
    cdef int i
    if f.v:
        sage_free(f.v)
        f.v = <ulong*> 0
    return 0

cdef int Pmodint_set(Pmodint* g, Pmodint f) except -1:
    """
    Set g equal to f.  We assume g has already been initialized.
    """
    cdef int i

    Pmodint_clear(g)
    Pmodint_init(g, f.p, f.degree)
    for i from 0 <= i <= f.degree:
        g.v[i] = f.v[i]
    return 0

cdef int Pmodint_cmp(Pmodint f, Pmodint g):
    cdef int i, j

    if f.degree < g.degree or f.p != g.p:
        return -1
    elif f.degree > g.degree:
        return 1
    for i from 0 <= i <= f.degree:
        if f.v[i] < g.v[i]:
            return -1
        elif f.v[i] > g.v[i]:
            return 1
    return 0

cdef int Pmodint_add_(Pmodint* sum, Pmodint f, Pmodint g) except -1:
    """
    Computes the sum of f and g and places the result in sum.`
    Note that sum must have been initialized.
    """
    if f.p != g.p:
        raise ArithmeticError, "Incompatible moduli."

    cdef int i, j

    Pmodint_clear(sum)
    if f.degree > g.degree:
        Pmodint_init(sum, f.p, f.degree)
        for j from 0 <= j <= g.degree:
            sum.v[j] = (f.v[j] + g.v[j]) % f.p
        for j from g.degree + 1 <= j <= f.degree:
            sum.v[j] = f.v[j]
    elif g.degree > f.degree:
        Pmodint_init(sum, f.p, g.degree)
        for j from 0 <= j <= f.degree:
            sum.v[j] = (g.v[j] + f.v[j]) % f.p
        for j from f.degree + 1 <= j <= g.degree:
            sum.v[j] = g.v[j]
    else:
        i = f.degree
        while i >= 0 and (f.v[i] + g.v[i]) % f.p == 0:
            i = i - 1
        Pmodint_init(sum, f.p, i)
        for j from 0 <= j <= i:
            sum.v[j] = (f.v[j] + g.v[j]) % f.p
    return 0


cdef int Pmodint_sub_(Pmodint* sum, Pmodint f, Pmodint g) except -1:
    """
    Computes the difference f - g and places the result in sum.`
    Note that sum must have been initialized.
    """
    if f.p != g.p:
        raise ArithmeticError, "Incompatible moduli."
    cdef ulong p
    p = f.p
    cdef int i, j

    Pmodint_clear(sum)
    if f.degree > g.degree:
        Pmodint_init(sum, p, f.degree)
        for j from 0 <= j <= g.degree:
            sum.v[j] = (f.v[j] +  (p - g.v[j])) % p
        for j from g.degree + 1 <= j <= f.degree:
            sum.v[j] = f.v[j]
    elif g.degree > f.degree:
        Pmodint_init(sum, g.p, g.degree)
        for j from 0 <= j <= f.degree:
            sum.v[j] = (f.v[j] + (p -g.v[j])) % p
        for j from f.degree + 1 <= j <= g.degree:
            sum.v[j] = f.p - g.v[j]
    else:
        i = f.degree
        while i >= 0 and f.v[i] == g.v[i]:
            i = i - 1
        Pmodint_init(sum, f.p, i)
        for j from 0 <= j <= i:
            sum.v[j] = (f.v[j] + (p - g.v[j])) % p
    return 0

cdef int Pmodint_scale(Pmodint* f, ulong a) except -1:
    """
    Replaces f by a*f.
    """
    cdef int i
    a = a % f.p
    for i from 0 <= i <= f.degree:
        f.v[i] = (f.v[i] * a) % f.p
    return 0

cdef int Pmodint_scalar_mul_(Pmodint* g, Pmodint f, ulong a) except -1:
    """
    Sets g fo a*f.
    """
    cdef int i
    a = a % f.p
    Pmodint_init(g, f.p, f.degree)
    for i from 0 <= i <= f.degree:
        g.v[i] = (f.v[i] * a) % f.p
    return 0


cdef int Pmodint_mul_(Pmodint* prod, Pmodint f, Pmodint g) except -1:
    """
    Computes the product of f and g and places the result in prod.`
    Note that prod must have been initialized.
    """
    if f.p != g.p:
        raise ArithmeticError, "Incompatible moduli."
    cdef int i, m, n, deg
    cdef ulong s, t
    Pmodint_clear(prod)
    if f.degree == -1 or g.degree == -1:
        Pmodint_init(prod, f.p, -1)
        return 0
    deg = f.degree + g.degree
    Pmodint_init(prod, f.p, deg)
    for n from 0 <= n <= deg:
        if f.degree < n:  # m = min(f.degree, n)
            m = f.degree
        else:
            m = n
        s = 0
        for i from 0 <= i <= m:
            if n - i <= g.degree and f.v[i] != 0 and g.v[n-i] != 0:
                t = (f.v[i] * g.v[n-i]) % f.p
                s = (s + t)%f.p
        prod.v[n] = s
    return 0

cdef int Pmodint_karatsuba_split(Pmodint *A, Pmodint *B, Pmodint f, int e) except -1:
    """
    Used by Pmodint_mul_karatsuba:
       A = f[e:]
       B = f[:e]
    """
    cdef int i, dA, dB
    if e > f.degree or e <= 0:
        s = "e = %s"%e
        raise IndexError, s

    Pmodint_clear(A)
    dA = f.degree - e
    while dA >= 0 and f.v[e+dA] == 0:
        dA = dA - 1
    Pmodint_init(A, f.p, dA)
    for i from e <= i <= e+dA:
        A.v[i-e] = f.v[i]

    Pmodint_clear(B)
    dB = e-1
    while dB >= 0 and f.v[dB] == 0:
        dB = dB - 1
    Pmodint_init(B, f.p, dB)
    for i from 0 <= i <= dB:
        B.v[i] = f.v[i]

    #print "f=%s, e=%s, A = %s, B=%s"%(Pmodint_repr(f), e, Pmodint_repr(A[0]), Pmodint_repr(B[0]))
    return 0

cdef int Pmodint_mul_karatsuba(Pmodint* prod, Pmodint f, Pmodint g) except -1:
    """
    Computes the product of f and g using Karatsuba multiplication.
    Note that prod must have been initialized.

    The basic idea is to use that

       (a*X+b)*(c*X+d) = ac*X^2+((a+b)*(c+d)-ac-bd)*X+bd

    where ac=a*c and bd=b*d, which requires three
    multiplications instead of the naive four.  (In my examples,
    strangely just doing the above with four multiplications
    does tend to speed things up noticeably.)
    Given f and g of arbitrary degree bigger than one, let e
    be min(deg(f),deg(g))/2.  Write

          f = a*X^e + b    and    g = c*X^e + d

    and use the identity

         (a*X^e+b)*(c*X^e+d) = ac*X^(2e) +((a+b)*(c+d)-ac-bd)*X^e + bd

    to recursively compute f*g.

    TIMINGS:
    On a Pentium M 1.8Ghz laptop:

    NOTES:
        * Karatsuba multiplication of polynomials is also implemented in PARI in
          src/basemath/polarit3.c.

        * A much more readable pure python version is in python/rings.py

    """

    if f.degree <= 25 or g.degree <= 25:
       Pmodint_mul_(prod, f, g)
       return 0

    cdef int e, i
    cdef ulong p
    # tmp variables for karatsuba.
    cdef Pmodint k_A, k_B, k_C, k_D, k_AC, k_BD, k_A_plus_B, k_C_plus_D, \
                 k_AC_plus_BD, k_t1, k_t2, k_t3, k_X, k_Y
    if f.degree < g.degree:
        e = f.degree + 1
    else:
        e = g.degree + 1
    e = e/2

    p = f.p
    Pmodint_init(&k_A,p,-1); Pmodint_init(&k_B,p,-1)
    Pmodint_init(&k_C,p,-1); Pmodint_init(&k_D,p,-1)
    Pmodint_init(&k_AC,p,-1); Pmodint_init(&k_BD,p,-1)
    Pmodint_init(&k_A_plus_B,p,-1); Pmodint_init(&k_C_plus_D,p,-1)
    Pmodint_init(&k_AC_plus_BD,p,-1)
    Pmodint_init(&k_t1,p,-1); Pmodint_init(&k_t2,p,-1); Pmodint_init(&k_t3,p,-1)
    Pmodint_init(&k_X,p,-1); Pmodint_init(&k_Y,p,-1)

    Pmodint_karatsuba_split(&k_A, &k_B, f, e)
    Pmodint_karatsuba_split(&k_C, &k_D, g, e)
    Pmodint_mul_karatsuba(&k_AC, k_A, k_C)
    Pmodint_mul_karatsuba(&k_BD, k_B, k_D)
    Pmodint_add_(&k_A_plus_B, k_A, k_B)
    Pmodint_add_(&k_C_plus_D, k_C, k_D)
    Pmodint_mul_karatsuba(&k_X, k_A_plus_B, k_C_plus_D)         # X = (A+B)*(C+D)
    Pmodint_add_(&k_AC_plus_BD, k_AC, k_BD)
    Pmodint_sub_(&k_Y, k_X, k_AC_plus_BD)

    Pmodint_init(&k_t2, p, 2*e + k_AC.degree)
    for i from 0 <= i < 2*e:
        k_t2.v[i] = 0
    for i from 0 <= i <= k_AC.degree:
        k_t2.v[2*e+i] = k_AC.v[i]

    Pmodint_init(&k_t1, p, e + k_Y.degree)
    for i from 0 <= i < e:
        k_t1.v[i] = 0
    for i from 0 <= i <= k_Y.degree:
        k_t1.v[e+i] = k_Y.v[i]

    Pmodint_add_(&k_t3, k_t1, k_t2)
    Pmodint_add_(prod, k_t3, k_BD)                  # prod = t1 + t2 + t3

    Pmodint_clear(&k_A); Pmodint_clear(&k_B); Pmodint_clear(&k_C); Pmodint_clear(&k_D)
    Pmodint_clear(&k_AC); Pmodint_clear(&k_BD)
    Pmodint_clear(&k_A_plus_B); Pmodint_clear(&k_C_plus_D); Pmodint_clear(&k_AC_plus_BD)
    Pmodint_clear(&k_t1); Pmodint_clear(&k_t2); Pmodint_clear(&k_t3)
    Pmodint_clear(&k_X); Pmodint_clear(&k_Y)
    return 0

cdef object Pmodint_list(Pmodint f):
    cdef int i
    v = []
    for i from 0 <= i <= f.degree:
        v.append(f.v[i])
    return v

cdef object Pmodint_repr(Pmodint f):
    s = " "
    X = "x"
    for n from f.degree >= n >= 0:
        x = f.v[n]
        if x != 0:
            if s != " ":
                s = s + " + "
            s = s + "%s*%s^%s"%(x,X,n)
    s=s.replace(" 1*"," ")
    s = s + "  (degree %s)"%f.degree
    if s==" ":
        return "0"
    return s


cdef int Pmodint_setitem(Pmodint* f, int n, ulong c) except -1:
    """
    Sets the coefficient of x^n in the polynomial to c.  Does not
    steal a reference to c.
    """
    cdef Pmodint g
    cdef int i, j
    c = c % f.p

    if n < 0:
        raise IndexError, "Negative polynomial coefficients not defined."
    if n <= f.degree:
        f.v[n] = c
        if n == f.degree and c == 0:
            # Setting the leading coefficient to 0 reduces the degree:
            i = f.degree-1
            while i >= 0 and f.v[i] == 0:
                i = i - 1
            Pmodint_init(&g, f.p, i)
            for j from 0 <= j <= i:
                g.v[j] = f.v[j]
            Pmodint_clear(f)   # delete f.v
            f.v = g.v
            f.degree = i
        #endif
        return 0
    if c == 0:
        return 0
    Pmodint_init(&g, f.p, n)
    for i from 0 <= i <= f.degree:
        g.v[i] = f.v[i]
    for i from f.degree+1 <= i < n:
        g.v[i] = 0
    g.v[n] = c
    Pmodint_clear(f)
    f.degree = n
    f.v = g.v
    return 0


cdef class Polynomial_modint:
    """
    Polynomial_modint():

    Polynomials over the integers modulo p.
    """
    cdef Pmodint poly    # polynomial over Z/pZ

    def __new__(self, int p):
        Pmodint_init(&self.poly, p, -1)

    def __init__(self, int p):
        pass

    def __dealloc__(self):
        #Pmodint_clear(&self.poly)
        pass

    def __add__(Polynomial_modint self, Polynomial_modint other):
        cdef Polynomial_modint f
        f = Polynomial_modint(self.poly.p)
        Pmodint_add_(&f.poly, self.poly, other.poly)
        return f

    def __sub__(Polynomial_modint self, Polynomial_modint other):
        cdef Polynomial_modint f
        f = Polynomial_modint(self.poly.p)
        Pmodint_sub_(&f.poly, self.poly, other.poly)
        return f

    def __mul__(Polynomial_modint self, Polynomial_modint other):
        cdef Polynomial_modint f
        f = Polynomial_modint(self.poly.p)
        Pmodint_mul_karatsuba(&f.poly, self.poly, other.poly)
        #Pmodint_mul_(&f.poly, self.poly, other.poly)
        return f

    # trick so that we can do type coercision.
    def __cmp(Polynomial_modint self, Polynomial_modint other):
            return Pmodint_cmp(self.poly, other.poly)
    def __cmp__(Polynomial_modint self, other):
        if not isinstance(other, Polynomial_modint):
            return -1
        return self.__cmp(other)

    def __repr__(self):
        return Pmodint_repr(self.poly)

    def __setitem__(self, n, object x):
        cdef int c
        c = long(x) % self.poly.p
        Pmodint_setitem(&self.poly, n, c)

    def __getitem__(self, int n):
        if n > self.poly.degree or n < 0:
            raise IndexError
        return self.poly.v[n]

    def copy(self):
        cdef Polynomial_modint f
        f = Polynomial_modint(f.poly.p)
        Pmodint_set(&f.poly, self.poly)
        return f

    def degree(self):
        return self.poly.degree

    def list(self):
        return Pmodint_list(self.poly)

    def prime(self):
        return self.poly.p

    def set_from_list(self, v):
        cdef int i
        Pmodint_clear(&self.poly)
        Pmodint_init(&self.poly, self.poly.p, len(v)+1)
        for i from 0 <= i < len(v):
            self[i] = v[i]

    def __pow__(self, n, m):
        _n = int(n)
        if _n != n:
            raise ValueError, "exponent must be an integer"
        ans = Polynomial_modint(self.prime())
        ans[0] = 1
        apow = self
        while _n != 0:
            if _n%2 != 0:
                ans = ans * apow
            apow = apow * apow
            _n = _n/2
        return ans


###########################################################################
# Polynomials over the rational numbers.
###########################################################################

# Representation: We represent polynomials using a C array of deg+1 GMP
# rationals.  Changing the degree entails re-allocating memory.

cdef struct PQ_data:
    int degree
    mpq_t* v

ctypedef PQ_data PQ

cdef int PQ_init(PQ* f, int degree) except -1:
    cdef int i
    if degree < -1:
        raise ArithmeticError, "The degree must be >= -1."
    f.degree = degree
    if degree == -1:
        f.v = <mpq_t*> 0
        return 0
    f.v = <mpq_t*>sage_malloc(sizeof(mpq_t)*(degree+1))
    if f.v == <mpq_t*>0:
        raise MemoryError, "Error allocating memory for polynomial."
    for i from 0 <= i <= degree:
        mpq_init(f.v[i])
    return 0

cdef int PQ_clear(PQ* f) except -1:
    cdef int i
    for i from 0 <= i <= f.degree:
        mpq_clear(f.v[i])
    if f.v:
        sage_free(f.v)
        f.v = <mpq_t*> 0
    return 0

cdef int PQ_set(PQ* g, PQ f) except -1:
    """
    Set g equal to f.  We assume g has already been initialized.
    """
    cdef int i

    PQ_clear(g)
    PQ_init(g, f.degree)
    for i from 0 <= i <= f.degree:
        mpq_set(g.v[i], f.v[i])
    return 0

cdef int PQ_cmp(PQ f, PQ g):
    cdef int i, j

    if f.degree < g.degree:
        return -1
    elif f.degree > g.degree:
        return 1
    for i from 0 <= i <= f.degree:
        j = mpq_cmp(f.v[i], g.v[i])
        if j != 0:
            return j
    return 0

cdef mpq_t PQ_add_tmp
mpq_init(PQ_add_tmp)
cdef int PQ_add_(PQ* sum, PQ f, PQ g) except -1:
    """
    Computes the sum of f and g and places the result in sum.`
    Note that sum must have been initialized.
    """
    cdef int i, j

    PQ_clear(sum)
    if f.degree > g.degree:
        PQ_init(sum, f.degree)
        for j from 0 <= j <= g.degree:
            mpq_add(sum.v[j], f.v[j], g.v[j])
        for j from g.degree + 1 <= j <= f.degree:
            mpq_set(sum.v[j], f.v[j])
    elif g.degree > f.degree:
        PQ_init(sum, g.degree)
        for j from 0 <= j <= f.degree:
            mpq_add(sum.v[j], g.v[j], f.v[j])
        for j from f.degree + 1 <= j <= g.degree:
            mpq_set(sum.v[j], g.v[j])
    else:
        i = f.degree
        if i == -1:
            PQ_init(sum, i)
            return 0
        mpq_neg(PQ_add_tmp, g.v[i])
        while i >= 0 and mpq_equal(f.v[i], PQ_add_tmp):
            i = i - 1
            if i >= 0:
                mpq_neg(PQ_add_tmp, g.v[i])
        PQ_init(sum, i)
        for j from 0 <= j <= i:
            mpq_add(sum.v[j], f.v[j], g.v[j])
    return 0

cdef mpq_t PQ_sub_tmp
mpq_init(PQ_sub_tmp)
cdef int PQ_sub_(PQ* sum, PQ f, PQ g) except -1:
    """
    Computes the sum of f and g and places the result in sum.`
    Note that sum must have been initialized.
    """
    cdef int i, j

    PQ_clear(sum)
    if f.degree > g.degree:
        PQ_init(sum, f.degree)
        for j from 0 <= j <= g.degree:
            mpq_sub(sum.v[j], f.v[j], g.v[j])
        for j from g.degree + 1 <= j <= f.degree:
            mpq_set(sum.v[j], f.v[j])
    elif g.degree > f.degree:
        PQ_init(sum, g.degree)
        for j from 0 <= j <= f.degree:
            mpq_sub(sum.v[j], f.v[j], g.v[j])
        for j from f.degree + 1 <= j <= g.degree:
            mpq_neg(PQ_sub_tmp, g.v[j])
            mpq_set(sum.v[j], PQ_sub_tmp)
    else:
        i = f.degree
        while i >= 0 and mpq_equal(f.v[i], g.v[i]):
            i = i - 1
        PQ_init(sum, i)
        for j from 0 <= j <= i:
            mpq_sub(sum.v[j], f.v[j], g.v[j])
    return 0


cdef int PQ_scale(PQ* f, mpq_t a) except -1:
    """
    Replaces f by a*f.
    """
    cdef int i
    for i from 0 <= i <= f.degree:
        mpq_mul(f.v[i], f.v[i], a)
    return 0

cdef int PQ_scalar_mul_(PQ* g, PQ f, mpq_t a) except -1:
    """
    Sets g fo a*f.
    """
    cdef int i
    PQ_init(g, f.degree)
    for i from 0 <= i <= f.degree:
        mpq_mul(g.v[i], a, f.v[i])
    return 0


cdef mpq_t PQ_mul_tmp, PQ_mul_tmp2
mpq_init(PQ_mul_tmp); mpq_init(PQ_mul_tmp2)

cdef int PQ_mul_(PQ* prod, PQ f, PQ g) except -1:
    """
    Computes the product of f and g and places the result in prod.`
    Note that prod must have been initialized.
    """
    cdef int i, m, n, deg
    PQ_clear(prod)
    if f.degree == -1 or g.degree == -1:
        PQ_init(prod, -1)
        return 0
    deg = f.degree + g.degree
    PQ_init(prod, deg)
    for n from 0 <= n <= deg:
        mpq_set_si(PQ_mul_tmp, 0, 1)
        if f.degree < n:  # m = min(f.degree, n)
            m = f.degree
        else:
            m = n
        for i from 0 <= i <= m:
            if n - i <= g.degree and mpq_sgn(f.v[i]) and mpq_sgn(g.v[n-i]):
                mpq_mul(PQ_mul_tmp2, f.v[i], g.v[n-i])
                mpq_add(PQ_mul_tmp, PQ_mul_tmp, PQ_mul_tmp2)
        mpq_set(prod.v[n], PQ_mul_tmp)
    return 0

cdef int PQ_karatsuba_split(PQ *A, PQ *B, PQ f, int e) except -1:
    """
    Used by PQ_mul_karatsuba:
       A = f[e:]
       B = f[:e]
    """
    cdef int i, dA, dB
    if e > f.degree or e <= 0:
        s = "e = %s"%e
        raise IndexError, s

    PQ_clear(A)
    dA = f.degree - e
    while dA >= 0 and mpq_sgn(f.v[e+dA]) == 0:
        dA = dA - 1
    PQ_init(A, dA)
    for i from e <= i <= e+dA:
        mpq_set(A.v[i-e], f.v[i])

    PQ_clear(B)
    dB = e-1
    while dB >= 0 and mpq_sgn(f.v[dB]) == 0:
        dB = dB - 1
    PQ_init(B, dB)
    for i from 0 <= i <= dB:
        mpq_set(B.v[i], f.v[i])

    #print "f=%s, e=%s, A = %s, B=%s"%(PQ_repr(f), e, PQ_repr(A[0]), PQ_repr(B[0]))
    return 0

cdef int PQ_mul_karatsuba(PQ* prod, PQ f, PQ g) except -1:
    """
    Computes the product of f and g using Karatsuba multiplication.
    Note that prod must have been initialized.

    The basic idea is to use that

       (a*X+b)*(c*X+d) = ac*X^2+((a+b)*(c+d)-ac-bd)*X+bd

    where ac=a*c and bd=b*d, which requires three
    multiplications instead of the naive four.  (In my examples,
    strangely just doing the above with four multiplications
    does tend to speed things up noticeably.)
    Given f and g of arbitrary degree bigger than one, let e
    be min(deg(f),deg(g))/2.  Write

          f = a*X^e + b    and    g = c*X^e + d

    and use the identity

         (a*X^e+b)*(c*X^e+d) = ac*X^(2e) +((a+b)*(c+d)-ac-bd)*X^e + bd

    to recursively compute f*g.

    TIMINGS:
    On a Pentium M 1.8Ghz laptop:

    NOTES:
        * Karatsuba multiplication of polynomials is also implemented in PARI in
          src/basemath/polarit3.c.

        * A much more readable pure python version is in python/rings.py

    """

    if f.degree <= 10 or g.degree <= 10:
        PQ_mul_(prod, f, g)
        return 0

    cdef int e, i
    # tmp variables for karatsuba.
    cdef PQ k_A, k_B, k_C, k_D, k_AC, k_BD, k_A_plus_B, k_C_plus_D,  k_AC_plus_BD, k_t1, k_t2, k_t3, k_X, k_Y
    if f.degree < g.degree:
        e = f.degree + 1
    else:
        e = g.degree + 1
    e = e/2

    PQ_init(&k_A,-1); PQ_init(&k_B,-1); PQ_init(&k_C,-1); PQ_init(&k_D,-1)
    PQ_init(&k_AC,-1); PQ_init(&k_BD,-1)
    PQ_init(&k_A_plus_B,-1); PQ_init(&k_C_plus_D,-1); PQ_init(&k_AC_plus_BD,-1)
    PQ_init(&k_t1,-1); PQ_init(&k_t2,-1); PQ_init(&k_t3,-1)
    PQ_init(&k_X,-1); PQ_init(&k_Y,-1)

    PQ_karatsuba_split(&k_A, &k_B, f, e)
    PQ_karatsuba_split(&k_C, &k_D, g, e)
    PQ_mul_karatsuba(&k_AC, k_A, k_C)
    PQ_mul_karatsuba(&k_BD, k_B, k_D)
    PQ_add_(&k_A_plus_B, k_A, k_B)
    PQ_add_(&k_C_plus_D, k_C, k_D)
    PQ_mul_karatsuba(&k_X, k_A_plus_B, k_C_plus_D)         # X = (A+B)*(C+D)
    PQ_add_(&k_AC_plus_BD, k_AC, k_BD)
    PQ_sub_(&k_Y, k_X, k_AC_plus_BD)

    PQ_init(&k_t2, 2*e + k_AC.degree)
    for i from 0 <= i < 2*e:
        mpq_set_si(k_t2.v[i], 0, 1)
    for i from 0 <= i <= k_AC.degree:
        mpq_set(k_t2.v[2*e+i], k_AC.v[i])

    PQ_init(&k_t1, e + k_Y.degree)
    for i from 0 <= i < e:
        mpq_set_si(k_t1.v[i], 0, 1)
    for i from 0 <= i <= k_Y.degree:
        mpq_set(k_t1.v[e+i], k_Y.v[i])

    PQ_add_(&k_t3, k_t1, k_t2)
    PQ_add_(prod, k_t3, k_BD)                  # prod = t1 + t2 + t3

    PQ_clear(&k_A); PQ_clear(&k_B); PQ_clear(&k_C); PQ_clear(&k_D)
    PQ_clear(&k_AC); PQ_clear(&k_BD)
    PQ_clear(&k_A_plus_B); PQ_clear(&k_C_plus_D); PQ_clear(&k_AC_plus_BD)
    PQ_clear(&k_t1); PQ_clear(&k_t2); PQ_clear(&k_t3)
    PQ_clear(&k_X); PQ_clear(&k_Y)
    return 0

cdef int PQ_num_den(mpz_t **vec, mpz_t den, PQ f) except -1:
    """
    Sets den equal to the least common multiple of the denominators
    of the coefficients of f, and vec equal to the coefficients
    of den*f.  This function assumes that vec has *not*
    been initialized in any way, but den better have been initialized.

    Raises a MemoryError exception, if this function fails to allocate
    the memory needed for the numerators.
    """
    cdef int i
    cdef mpz_t tmp

    # Allocate memory for answer, and raise exception on failure.
    vec[0] = <mpz_t *> sage_malloc(sizeof(mpz_t) * (f.degree+1))
    if vec[0] == <mpz_t *> 0:
        raise MemoryError

    # Go through and find the least commmon multiple of the denominators
    mpz_set_si(den, 1)
    for i from 0 <= i <= f.degree:
        mpz_lcm(den, den, mpq_denref(f.v[i]))

    # Fill vec with the products c*d, where c are the coefficients
    if mpz_sgn(den) == 0:   # denominator is 1
        for i from 0 <= i <= f.degree:
            mpz_init(vec[0][i])
            mpz_set(vec[0][i], mpq_numref(f.v[i]))
    else:
        mpz_init(tmp)
        for i from 0 <= i <= f.degree:
            mpz_init(vec[0][i])
            mpz_divexact(tmp, den, mpq_denref(f.v[i]))
            mpz_mul(vec[0][i], mpq_numref(f.v[i]), tmp)
        mpz_clear(tmp)

    return 0


cdef object to_str(mpz_t x):
    """
    Convert a GMP integer to a Python string.
    """
    cdef char *s
    s = mpz_get_str(NULL, 10, x)
    t = s
    free(s)
    return t

cdef int PQ_mul_modular_alg(PQ* prod, PQ f, PQ g) except -1:
    """
    Computes the product of f and g and places the result in prod.
    Uses a multi-modular method.
    Note that prod must have been initialized already.

    PERFORMANCE:
       On a 1.8Ghz Pentium-M laptop.  It's comparable to PARI, and
       when the degree is large it is much better.  MAGMA wins in
       all cases, but that's maybe due to it using an FFT algorithm?
       The times suggest this code is at least usable in practice.
       It's much faster than pure python.

       Multiplying sum_{n=0}^{10000} (n/19)*x^n by times itself.
       THIS Function     PARI       MAGMA    (time in seconds)
       1.56              7.78       0.12

       Multiplying x^20 + 19/5*x^7 + 1 by x^35 + 2*x + 1  10000 times.
       THIS              PARI       MAGMA      (PURE PYTHON)
       0.63              0.15       0.09        27 seconds

       Multiplying sum_{n=0}^100 n*x^n by sum_{n=0}^100 (100-n)*x^n 1000 times
       THIS              PARI       MAGMA
       0.64              0.748      0.21

       Multiplying sum_{n=1}^100 1/n * x^n by itself 100 times
       THIS              PARI       MAGMA
       1.14              0.716      0.333
    """
    PQ_clear(prod)
    if f.degree == -1 or g.degree == -1:
        PQ_init(prod, -1)
        return 0

    cdef Pmodint fmod, gmod, fgmod
    cdef mpz_t *f_vec, *g_vec, *fg_vec, *lift, *lift2, f_den, g_den, den, p, pr, pr2, TWO, B, H_f, H_g
    cdef int i, degree
    cdef ulong p_int

    mpz_init_set_si(TWO, 2)

    # Compute vectors over mpz_t elements from f and g, along with denominators;
    # input is assumed not-initialized
    mpz_init(f_den)
    mpz_init(g_den)
    PQ_num_den(&f_vec, f_den, f)
    PQ_num_den(&g_vec, g_den, g)

    # multiply denominators, to get denominator of product
    mpz_init(den)
    mpz_mul(den, f_den, g_den)

    # we do not need the denominators of each individual polynomial any more
    mpz_clear(f_den)
    mpz_clear(g_den)

    # Initialize prime and product of primes and bound
    mpz_init(B)
    # Compute bound B, so if we know the coefficients of f_den*f and g_den*g
    # mod B, then we know them.
    mpz_init(H_f)
    mpz_init(H_g)
    ag.mpz_height_vec(H_f, f_vec, f.degree+1)
    ag.mpz_height_vec(H_g, g_vec, g.degree+1)
    mpz_mul(B, H_f, H_g)
    mpz_mul(B, B, TWO)
    if f.degree >= g.degree:
        mpz_mul_si(B, B, (f.degree+1))
    else:
        mpz_mul_si(B, B, (g.degree+1))

    # p is the prime we are working modulo
    mpz_init(p)
    # pr is the product of the primes before p.
    # after the while loop below, it will be the product of all primes used in the computation.
    mpz_init(pr)
    # We compute the answer modulo enough primes so that the product is bigger than B.
    mpz_set_si(p, START_PRIME)
    mpz_set_si(pr, 1)
    degree = f.degree + g.degree

    first = True
    while first or mpz_cmp(pr, B) < 0:
        #print "p = ", to_str(p)
        p_int = mpz_get_ui(p)
        # Create the reductions of f and g mod p
        Pmodint_init(&fmod, p_int, -1)
        Pmodint_init(&gmod, p_int, -1)

        ag.mpzvec_to_intmod(&fmod.v, f_vec, f.degree+1, p_int)
        fmod.degree = f.degree
        ag.mpzvec_to_intmod(&gmod.v, g_vec, g.degree+1, p_int)
        gmod.degree = g.degree

        # Multiply the reductions
        Pmodint_init(&fgmod, p_int, -1)
        Pmodint_mul_karatsuba(&fgmod, fmod, gmod)

        # Now free up memory used by fmod and gmod
        Pmodint_clear(&fmod)
        Pmodint_clear(&gmod)

        # CRT with answer so far
        if first:
            ag.intmodvec_to_mpz(&lift, fgmod.v, degree+1)
            Pmodint_clear(&fgmod)
            first = False
        else:
            # Lift the mod p array to an array of GMP integers
            ag.intmodvec_to_mpz(&fg_vec, fgmod.v, degree+1)
            mpz_crt_vec(&lift2, lift, fg_vec, degree+1, pr, p)
            # De-allocate memory used in lift
            ag.mpzvec_clear(lift, degree+1)
            sage_free(lift)
            lift = lift2
            ag.mpzvec_clear(fg_vec, degree+1)
            sage_free(fg_vec)
            Pmodint_clear(&fgmod)
        #endif

        # Use another prime.
        mpz_mul(pr, pr, p)
        mpz_set_si(p,next_probab_prime_int(p_int))
    #end while

    PQ_init(prod, degree)

    mpz_init(pr2)
    mpz_fdiv_q(pr2, pr, TWO)
    for i from 0 <= i <= degree:
        if mpz_cmp(lift[i],pr2) > 0:
            mpz_sub(lift[i], lift[i], pr)
        mpz_set(mpq_numref(prod.v[i]), lift[i])
        mpz_set(mpq_denref(prod.v[i]), den)
    if mpz_cmp_si(den, 1) != 0:
        for i from 0 <= i <= degree:
            mpq_canonicalize(prod.v[i])
    mpz_clear(den)
    mpz_clear(p)
    mpz_clear(pr)
    mpz_clear(pr2)
    mpz_clear(TWO)
    mpz_clear(B)
    mpz_clear(H_f)
    mpz_clear(H_g)
    ag.mpzvec_clear(lift, degree+1)
    sage_free(lift)
    ag.mpzvec_clear(f_vec, f.degree+1)
    sage_free(f_vec)
    ag.mpzvec_clear(g_vec, g.degree+1)
    sage_free(g_vec)
    return 0

cdef object mpq_to_str(mpq_t x):
    cdef char *s
    s = mpq_get_str(NULL, 10, x)
    t = str(s)
    free(s)
    return t

cdef object PQ_list(PQ f):
    cdef int i
    v = []
    for i from 0 <= i <= f.degree:
        v.append(mpq_to_str(f.v[i]))
    return v

cdef object PQ_repr(PQ f):
    s = " "
    X = "x"
    for n from f.degree >= n >= 0:
        x = mpq_to_str(f.v[n])
        if x != "0":
            if s != " ":
                s = s + " + "
            s = s + "%s*%s^%s"%(x,X,n)
    s=s.replace(" 1*"," ")
    s = s + "  (degree %s)"%f.degree
    if s==" ":
        return "0"
    return s


cdef int PQ_setitem(PQ* f, int n, mpq_t c) except -1:
    """
    Sets the coefficient of x^n in the polynomial to c.  Does not
    steal a reference to c.
    """
    cdef PQ g
    cdef int i, j

    if n < 0:
        raise IndexError, "Negative polynomial coefficients not defined."
    if n <= f.degree:
        mpq_set(f.v[n], c)
        if n == f.degree and mpq_sgn(c) == 0:
            # Setting the leading coefficient to 0 reduces the degree:
            i = f.degree-1
            while i >= 0 and mpq_sgn(f.v[i]) == 0:
                i = i - 1
            PQ_init(&g, i)
            for j from 0 <= j <= i:
                mpq_set(g.v[j], f.v[j])
            PQ_clear(f)   # delete f.v
            f.v = g.v
            f.degree = i
        #endif
        return 0
    if mpq_sgn(c) == 0:
        return 0
    PQ_init(&g, n)
    for i from 0 <= i <= f.degree:
        mpq_set(g.v[i],f.v[i])
    for i from f.degree+1 <= i < n:
        mpq_set_si(g.v[i], 0, 1)
    mpq_set(g.v[n], c)
    PQ_clear(f)
    f.degree = n
    f.v = g.v
    return 0


cdef mpq_t class_tmp
mpq_init(class_tmp)
cdef class Polynomial_rational(sage.structure.element.RingElement):
    """
    Polynomial_rational():

    Create the zero polynomial over the rational numbers.
    """
    cdef PQ pq    # polynomial over Q

    def __new__(self):
        PQ_init(&self.pq, -1)

    def __dealloc__(self):
        PQ_clear(&self.pq)

    def __add__(Polynomial_rational self, Polynomial_rational other):
        cdef Polynomial_rational f
        f = Polynomial_rational()
        PQ_add_(&f.pq, self.pq, other.pq)
        return f

    def __sub__(Polynomial_rational self, Polynomial_rational other):
        cdef Polynomial_rational f
        f = Polynomial_rational()
        PQ_sub_(&f.pq, self.pq, other.pq)
        return f

    def __mul__(Polynomial_rational self, Polynomial_rational other):
        cdef Polynomial_rational f
        f = Polynomial_rational()
        # To switch to classical or karatsuba multiplication uncomment one of the
        # following and comment out the others.
        ##PQ_mul_(&f.pq, self.pq, other.pq)
        ##PQ_mul_karatsuba(&f.pq, self.pq, other.pq)
        PQ_mul_modular_alg(&f.pq, self.pq, other.pq)
        return f

    # trick so that we can do type coercion
    def __cmp(Polynomial_rational self, Polynomial_rational other):
            return PQ_cmp(self.pq, other.pq)

    def __cmp__(Polynomial_rational self, other):
        if not isinstance(other, Polynomial_rational):
            return -1
        return self.__cmp(other)

    def __repr__(self):
        return PQ_repr(self.pq)

    def __setitem__(self, n, object x):
        py_mpq_set(class_tmp, x)
        PQ_setitem(&self.pq, n, class_tmp)

    def __getitem__(self, int n):
        if n > self.pq.degree or n < 0:
            return '0'
        return mpq_to_str(self.pq.v[n])

    def __getslice__(self, Py_ssize_t i, Py_ssize_t j):
        cdef Py_ssize_t k
        z = []
        for k in i <= k < j:
            z.append(mpq_to_str(self.pq.v[k]))

    def copy(self):
        cdef Polynomial_rational f
        f = Polynomial_rational()
        PQ_set(&f.pq, self.pq)
        return f

    def degree(self):
        return self.pq.degree

    def __nonzero__(self):
        return self.pq.degree != -1

    def list(self):
        return PQ_list(self.pq)

    def set_from_list(self, v):
        cdef int i
        PQ_clear(&self.pq)
        PQ_init(&self.pq, len(v)-1)
        for i from 0 <= i < len(v):
            self[i] = v[i]

    def _pow(self, n):
        """
        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = x^3 + 2*x^2 + 323
            sage: f^3
            x^9 + 6*x^8 + 12*x^7 + 977*x^6 + 3876*x^5 + 3876*x^4 + 312987*x^3 + 625974*x^2 + 33698267

        AUTHORS:
            -- William Stein and Didier Deshom
        """
        # We may assume the input n is an integer (not a rational or something),
        # since the
        # function _pow_ is called by __pow__ after doing the conversion.
        ans = Polynomial_rational()
        ans[0] = 1
        apow = self
        if n < 0:
            raise ValueError, "power n must be nonnegative"
        while True:
            if n&1 > 0: ans = ans*apow
            n = n >> 1
            if n != 0:
                apow = apow*apow
            else:
                break
        return ans

#     def __div__(Polynomial_rational self, Polynomial_rational other):
#         """
#         Division with remainder.  Returns a tuple (quotient, remainder).
#         """
#         cdef Polynomial_rational A, B, Q, R, S, X
#         cdef mpq_t q, *quo_coeffs
#         cdef int n

#         if other.is_zero():
#             raise ZeroDivisionError
#         # This is algorithm 3.1.1 in Cohen GTM 138
#         A = self
#         B = other
#         R = A
#         mpq_set_si(class_tmp, 1, 1)
#         PQ_setitem(&X.pq, 1, class_tmp)

#         # The quotient will have degree deg(self) - deg(other)
#         n = self.pq.degree - other.pq.degree
#         quo_coeffs = <mpq_t> sage_malloc(sizeof(mpq_t)*n)

#         while R.pq.degree >= B.pq.degree:
#             # Here's what we do below:
#             # S =  (R.leading()/B.leading()) * X**(R.degree()-B.degree())
#             # Q = Q + S
#             # R = R - S*B
#             # Of course, those three steps can be done more quickly directly.
#             mpq_div(q, R.pq[R.pq.degree], B.pq[B.pq.degree])
#             n = R.pq.degree - B.pq.degree


#         Q = Polynomial_rational()
#         return (Q, R)


