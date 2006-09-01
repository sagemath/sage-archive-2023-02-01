"""
Sparse polynomial rings in pyrex
"""

# NOTE -- I implemented the algorithms in 2.2.4 of Knuth using
# circular linked lists.  I then tried squaring the polynomial sum
# n*x^n, for n < 1000, and it took 1.24 seconds, which is long.  So
# probably for many applications a dense polynomial implementation
# would be much better.  For multivariable poly rings, this would
# probably be useful.

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


include "cdefs.pxi"
include "gmp.pxi"
include 'interrupt.pxi'

cdef int py_mpq_set(mpq_t x, y) except -1:
    s = str(y)
    mpq_set_str(x, s, 0)
    mpq_canonicalize(x)

cdef mpq_t mpq_tmp       # pre-allocate an mpq variable for efficiency -- this could be a problem for multi-threading!!
mpq_init(mpq_tmp)

###########################################################################
# Abstract base class for Polynomials.
###########################################################################

cdef class Polynomial:
    def __new__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        pass

    def __getitem__(self, n):
        raise NotImplementedError

    def __setitem__(self, n, x):
        raise NotImplementedError

    def __getslice__(self, i, j):
        raise NotImplementedError

    def __call__(self, x):
        raise NotImplementedError

    def __mul__( self,  other):
        raise NotImplementedError

    def __add__( self,  other):
        raise NotImplementedError

    def __sub__( self,  other):
        raise NotImplementedError

    def __cmp__(self, other):
        raise NotImplementedError

    def __repr__(self):
        return "A Polynomial"

    def __hash__(self):
        raise NotImplementedError

    def __int__(self):
        raise NotImplementedError

    def __long__(self):
        raise NotImplementedError

    def __pow__(self, int n, int m):
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError

    def __pos__(self):
        raise NotImplementedError

    def __div__(self, other):
        """Division with remainder.  Returns a tuple (quotient, remainder)."""
        raise NotImplementedError

    def __floordiv__(self, right):
        """Quotient of division of self by other.  This is denoted //."""
        raise NotImplementedError

    def __mod__(self, other):
        """Remainder of division of self by other."""
        raise NotImplementedError

    def base_ring(self):
        return self.parent().base_ring()

    def copy(self):
        raise NotImplementedError

    def degree(self):
        raise NotImplementedError

    def derivative(self):
        raise NotImplementedError

    def factor(self):
        raise NotImplementedError

    def gcd( self,  other):
        raise NotImplementedError

    def is_irreducible(self):
        raise NotImplementedError

    def is_zero(self):
        raise NotImplementedError

    def leading_coefficient(self):
        raise NotImplementedError

    def randomize(self, degree, bound=10):
        raise NotImplementedError

    def resultant( self,  other):
        raise NotImplementedError

    def reverse(self):
        raise NotImplementedError

    def valuation(self):
        raise NotImplementedError

    def variable(self):
        return self.parent().variable()

    def xgcd( self,  other):
        raise NotImplementedError


###########################################################################
# Polynomials over the rational numbers.
###########################################################################

# Representation: We represent polynomials using a circular linked
# list, as described in Section 2.2.4 of Knuth, Art of Computer
# Programming, Vol 1.  Each coefficient of a polynomial is a GMP
# rational number.

# A term is a single term of a polynomial or a special node.
# From Knuth: "There is a special node at the end of every polynomial that
# has n = -1 and coef = 0.  This special node is a great convenience ...
# because it provides a convenient sentinel and it avoids the problem of
# an empty list (corresponding to the polynomial 0).

cdef struct polyrat_term:
    mpq_t coef
    int n           # the exponent, so this term is x**n
    polyrat_term* next

cdef int term_init(polyrat_term t) except -1:
    if t.n != -1:
        mpq_init(t.coef)
    return 0

cdef int term_clear(polyrat_term t) except -1:
    if t.n != -1:
        mpq_clear(t.coef)
    return 0

cdef polyrat_term *term_new(int n) except NULL:
    cdef polyrat_term *t
    t = <polyrat_term*> PyMem_Malloc(sizeof(polyrat_term))
    if t == <polyrat_term*>0:
        raise MemoryError, "Error allocating memory for polynomial."
    t.n = n
    if n != -1: mpq_init(t.coef)
    return t

cdef int term_delete(polyrat_term* t) except -1:
    if t.n != -1:
        mpq_clear(t.coef)
    PyMem_Free(t)
    return 0


# A polynomial a pointer to a polyrat_term, which is a pointer to
# the special node of the circular linked list that represents the polynomial.

ctypedef polyrat_term* polyrat

cdef polyrat polyrat_new() except NULL:
    cdef polyrat t
    t = <polyrat>PyMem_Malloc(sizeof(polyrat_term))
    if t == <polyrat>0:
        raise MemoryError, "Error allocating memory for polynomial."
    t = term_new(-1)
    t[0].next = t
    return t

cdef polyrat polyrat_copy(polyrat t) except NULL:
    cdef polyrat s, s2
    cdef polyrat_term *z
    s = polyrat_new()
    if t != t.next:
        s2 = s
        t = t.next
        while t.n != -1:
            z = term_new(t.n)
            mpq_set(z.coef, t.coef)
            s2.next = z
            s2 = z
            t = t.next
        s2.next = s
    return s

cdef int polyrat_free(polyrat t) except -1:
    cdef polyrat s
    t = t.next
    while t.n != -1:
        s = t.next
        term_clear(t[0])
        PyMem_Free(t)
        t = s
    term_clear(t[0])
    PyMem_Free(t)
    return 0

cdef int polyrat_add_into(polyrat P, polyrat Q) except -1:
    """
    This algorithm adds P to Q and replaces Q by the sum.
    """
    cdef polyrat_term *Q1, *Q2
    if P == Q:
        raise RuntimeError, "Input to polyrat_mul_into must be distinct."

    # The comments below are adapted from Section 2.2.4 of Knuth, Vol 1.

    # A1 [Initialize]  Now both P and Q point to the leading terms of
    # their polynomials.  Throughout most of this algorithm the
    # variable Q1 will be one step behind Q, in the sense that Q = Q1.next.
    P = P.next
    Q1 = Q
    Q = Q.next

    while 1:
        if P.n < Q.n:
            # A2. [P.n:Q.n]
            Q1 = Q
            Q = Q.next
        elif P.n == Q.n:
            # A3 [Add coefficients] We've found terms with equal exponents.
            if P.n < 0:
                return 0
            mpq_add(Q.coef, Q.coef, P.coef)
            if mpq_sgn(Q.coef) == 0:  # Q.coef == 0
                # A4 [Delete zero term.]
                Q2 = Q
                Q = Q.next
                Q1.next = Q
                term_delete(Q2)
                P = P.next
            else:
                P = P.next
                Q1 = Q
                Q = Q.next
        else: # P.n > Q.n
            # A5 [Insert new term.]  P contains a term that is not
            # present in Q, so we insert it in Q.
            Q2 = term_new(P.n)
            mpq_set(Q2.coef, P.coef)
            Q2.next = Q
            Q1.next = Q2
            Q1 = Q2
            P = P.next

cdef polyrat polyrat_add_(polyrat P, polyrat Q) except NULL:
    cdef polyrat R
    R = polyrat_copy(Q)
    polyrat_add_into(P, R)
    return R


cdef int polyrat_mul_into(polyrat M, polyrat P, polyrat Q) except -1:
    """
    This function replaces Q by Q + M*P.

    The pointer Q must not equal the pointer M or P, since otherwise
    as the output gets stored in Q the values of M or P might change,
    hence the result would be wrong.
    """
    cdef polyrat_term *Q1, *Q2
    cdef n

    if M == Q or P == Q:
        raise RuntimeError, "Input to polyrat_mul_into must be distinct."

    # The comments below are adapted from Section 2.2.4 of Knuth, Vol 1.

    # A1 [Initialize]  Now both P and Q point to the leading terms of
    # their polynomials.  Throughout most of this algorithm the
    # variable Q1 will be one step behind Q, in the sense that Q = Q1.next.

    while 1:
        M = M.next
        if M.n < 0:
            return 0
        # M2. [Multiply Cycle] Perform the algorithm for polyrat_add_into,
        # except that wherever the notation P.n appears in that algorithm,
        # replace it by "(if P.n < 0: -1, otherwise P.n + M.n)";
        # wherever P.coef appears in that algorithm, replace it
        # by "P.coef * M.coef".   Then go back to step M1.
        P = P.next
        Q1 = Q
        Q = Q.next
        while 1:
            if P.n == -1:
                n = -1
            else:
                n = P.n + M.n
            if n < Q.n:
                # A2. [P.n:Q.n]
                Q1 = Q
                Q = Q.next
            elif n == Q.n:
                # A3 [Add coefficients] We've found terms with equal exponents.
                if n < 0:
                    break
                mpq_mul(mpq_tmp, P.coef, M.coef)
                mpq_add(Q.coef, Q.coef, mpq_tmp)   # Q.coef += P.coef*M.coef

                if mpq_sgn(Q.coef) == 0:  # Q.coef == 0
                    # A4 [Delete zero term.]
                    Q2 = Q
                    Q = Q.next
                    Q1.next = Q
                    term_delete(Q2)
                    P = P.next
                else:
                    P = P.next
                    Q1 = Q
                    Q = Q.next
            else: # n > Q.n
                # A5 [Insert new term.]  P contains a term that is not
                # present in Q, so we insert it in Q.
                Q2 = term_new(n)
                mpq_mul(mpq_tmp, P.coef, M.coef)
                mpq_set(Q2.coef, mpq_tmp)
                Q2.next = Q
                Q1.next = Q2
                Q1 = Q2
                P = P.next


cdef polyrat polyrat_mul_(polyrat P, polyrat Q) except NULL:
    cdef polyrat R
    R = polyrat_new()
    polyrat_mul_into(P, Q, R)
    return R

cdef object polyrat_list(polyrat t):
    cdef polyrat s
    v = []
    t = t.next
    while t.n != -1:
        v.append((mpq_get_str(NULL, 10, t.coef), t.n))
        t = t.next
    return v

cdef int polyrat_setitem(polyrat t, int n, mpq_t c) except -1:
    """
    Sets the coefficient of x^n in the polynomial to c.  Does not
    steal a reference to c.
    """
    cdef polyrat_term *z, *prev
    if n < 0:
        raise IndexError, "Can only set coefficients of x^n for n>=0."
    if n > 1073741823:  # this is 2^30-1.  We can safely add two of these without overflow
        raise IndexError, "Exponent can be at most 1073741823."
    # Find position to insert x or change entry to x
    prev = t
    while t.next.n  >= n:
        prev = t
        t = t.next

    if t.n == n:
        # There is already a nonzero coefficient of x^n, so we change it, or
        # delete it if c = 0.
        if mpq_sgn(c) == 0:
            # delete t
            prev.next = t.next
            term_delete(t)
        else:
            mpq_set(t.coef, c)
    else:
        # There is no coefficient of x^n yet, so we create a new term
        # and insert it into the circular linked list, unless c=0.
        if mpq_sgn(c) == 0:
            return 0
        z = term_new(n)
        mpq_set(z.coef, c)
        z.next = t.next
        t.next = z
    return 0

cdef class Polynomial_rational(Polynomial):
    """
    Polynomial_rational():

    Create the zero polynomial over the rational numbers.
    """
    cdef polyrat t     # pointer into circular linked list of terms

    def __add__(Polynomial_rational self, Polynomial_rational other):
        cdef polyrat sum
        cdef Polynomial_rational f
        f = Polynomial_rational()
        polyrat_free(f.t)
        _sig_on
        f.t = polyrat_add_(self.t, other.t)
        _sig_off
        return f

    def __dealloc__(self):
        polyrat_free(self.t)

    def __init__(self):
        pass

    def __mul__(Polynomial_rational self, Polynomial_rational other):
        cdef polyrat sum
        cdef Polynomial_rational f
        f = Polynomial_rational()
        polyrat_free(f.t)
        _sig_on
        f.t = polyrat_mul_(self.t, other.t)
        _sig_off
        return f

    def __new__(self):
        self.t = polyrat_new()

    def __repr__(self):
        return str(self.list())

    def __setitem__(self, n, object x):
        py_mpq_set(mpq_tmp, x)
        polyrat_setitem(self.t, n, mpq_tmp)

    def copy(self):
        cdef Polynomial_rational f
        f = Polynomial_rational()
        polyrat_free(f.t)
        f.t = polyrat_copy(self.t)
        return f

    def list(self):
        return polyrat_list(self.t)


