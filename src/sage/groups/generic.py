r"""
Miscellaneous generic functions

A collection of functions implementing generic algorithms in arbitrary
groups, including additive and multiplicative groups.

In all cases the group operation is specified by a parameter
'operation', which is a string either one of the set of
multiplication_names or addition_names specified below, or 'other'.
In the latter case, the caller must provide an identity, inverse() and
op() functions.

Also included are a generic function for computing multiples (or
powers), and an iterator for general multiples and powers.

EXAMPLES:

Some examples in the multiplicative group of a finite field:

Discrete logs:
    sage: K = GF(3^6,'b')
    sage: b = K.gen()
    sage: a = b^210
    sage: discrete_log(a, b, K.order()-1)
    210

Linear relation finder:
    sage: F.<a>=GF(3^6,'a')
    sage: a.multiplicative_order().factor()
    2^3 * 7 * 13
    sage: b=a^7
    sage: c=a^13
    sage: linear_relation(b,c,'*')
    (13, 7)
    sage: b^13==c^7
    True

Orders of elements:
    sage: k.<a> = GF(5^5)
    sage: b = a^4
    sage: order_from_multiple(b,5^5-1,operation='*')
    781
    sage: order_from_bounds(b,(5^4,5^5),operation='*')
    781

Some examples in the group of points of an elliptic curve over a finite field:

Discrete logs:
    sage: F=GF(37^2,'a')
    sage: E=EllipticCurve(F,[1,1])
    sage: F.<a>=GF(37^2,'a')
    sage: E=EllipticCurve(F,[1,1])
    sage: P=E(25*a + 16 , 15*a + 7 )
    sage: P.order()
    672
    sage: Q=39*P; Q
    (36*a + 32 : 5*a + 12 : 1)
    sage: discrete_log(Q,P,P.order(),'+')
    39

Linear relation finder:
    sage: F.<a>=GF(3^6,'a')
    sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
    sage: P=E(a^5 + a^4 + a^3 + a^2 + a + 2 , 0)
    sage: Q=E(2*a^3 + 2*a^2 + 2*a , a^3 + 2*a^2 + 1)
    sage: linear_relation(P,Q,'+')
    (1, 2)
    sage: P == 2*Q
    True

Orders of elements:
    sage: k.<a> = GF(5^5)
    sage: E = EllipticCurve(k,[2,4])
    sage: P = E(3*a^4 + 3*a , 2*a + 1 )
    sage: M = E.cardinality(); M
    3227
    sage: plist = M.prime_factors()
    sage: order_from_multiple(P, M, plist, operation='+')
    3227
    sage: Q = E(0,2)
    sage: order_from_multiple(Q, M, plist, operation='+')
    7
    sage: order_from_bounds(Q, Hasse_bounds(5^5), operation='+')
    7

"""

###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                          John Cremona  <john.cremona@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from copy import copy
import math

import sage.misc.misc as misc
import sage.misc.search
import sage.rings.integer_ring as integer_ring
import sage.rings.integer

#
# Lists of names (as strings) which the user may use to identify one
# of the standard operations:
#
multiplication_names = ( 'multiplication', 'times', 'product', '*')
addition_names       = ( 'addition', 'plus', 'sum', '+')


#
# This function was moved from sage/rings/arith.py
#
def power(a, m, one=1):
    """
    The m-th power of a, where m is a non-negative
    integer and a is a Python object on which
    multiplication is defined.  The exponentiation
    is done using the standard binary powering algorithm.

    EXAMPLES:
        sage: power(2,5)
        32
        sage: power(RealField()('2.5'),4)
        39.0625000000000
        sage: power(0,0)
        Traceback (most recent call last):
        ...
        ArithmeticError: 0^0 is undefined.
        sage: power(2,-3)
        1/8
    """
    return sage.structure.element.generic_power(a,m,one)

def multiple(a, n, operation='*', identity=None, inverse=None, op=None):
    """
    Returns n*a [or a**n], where n is any integer and a is a Python
    object on which a group operaton such as addition or
    multiplication is defined.  Uses the standard binary algorithm.

    EXAMPLES:
        sage: multiple(2,5)
        32
        sage: multiple(RealField()('2.5'),4)
        39.0625000000000
        sage: multiple(2,-3)
        1/8
        sage: multiple(2,100,'+') == 100*2
        True
        sage: multiple(2,100) == 2**100
        True
        sage: multiple(2,-100,) == 2**-100
        True
        sage: R.<x>=ZZ[]
        sage: multiple(x,100)
        x^100
        sage: multiple(x,100,'+')
        100*x
        sage: multiple(x,-10)
        1/x^10

    Idempotence is detected making this fast:
        sage: multiple(1,10^1000)
        1

        sage: E=EllipticCurve('389a1')
        sage: P=E(-1,1)
        sage: multiple(P,10,'+')
        (645656132358737542773209599489/22817025904944891235367494656 : 525532176124281192881231818644174845702936831/3446581505217248068297884384990762467229696 : 1)
        sage: multiple(P,-10,'+')
        (645656132358737542773209599489/22817025904944891235367494656 : -528978757629498440949529703029165608170166527/3446581505217248068297884384990762467229696 : 1)


    """
    from operator import inv, mul, neg, add

    if operation in multiplication_names:
        identity = a.parent()(1)
        inverse  = inv
        op = mul
    elif operation in addition_names:
        identity = a.parent()(0)
        inverse  = neg
        op = add
    else:
        if identity==None or inverse==None or op==None:
            raise ValueError, "identity, inverse and operation must all be specified"

    if n == 0:
        return identity

    if n < 0:
        n = -n
        a = inverse(a)

    if n == 1:
        return a

    # check for idempotence, and store the result otherwise
    aa = op(a,a)
    if aa == a:
        return a

    if n == 2:
        return aa

    if n == 3:
        return op(aa,a)

    if n == 4:
        return op(aa,aa)

    # since we've computed a^2, let's start squaring there
    # so, let's keep the least-significant bit around, just
    # in case.
    m = n & 1
    n = n >> 1

    # One multiplication can be saved by starting with
    # the second-smallest power needed rather than with 1
    # we've already squared a, so let's start there.
    apow = aa
    while n&1 == 0:
        apow = op(apow,apow)
        n = n >> 1
    power = apow
    n = n >> 1

    # now multiply that least-significant bit in...
    if m:
        power = op(power,a)

    # and this is straight from the book.
    while n != 0:
        apow = op(apow,apow)
        if n&1 != 0:
            power = op(power,apow)
        n = n >> 1

    return power


#
# Generic iterator for looping through multiples or powers
#

class multiples:
    """
    Return an iterator which runs through P0+i*P for i in range(n)

    P and P0 must be Sage objects in some group; if the operation in
    multiplcation then the returned values are P0*P**i.

    EXAMPLES
        sage: list(multiples(1,10))
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list(multiples(1,10,100))
        [100, 101, 102, 103, 104, 105, 106, 107, 108, 109]

        sage: E=EllipticCurve('389a1')
        sage: P=E(-1,1)
        sage: for Q in multiples(P,5): print Q, Q.height()/P.height()
        (0 : 1 : 0) 0.000000000000000
        (-1 : 1 : 1) 1.00000000000000
        (10/9 : -35/27 : 1) 4.00000000000000
        (26/361 : -5720/6859 : 1) 9.00000000000000
        (47503/16641 : 9862190/2146689 : 1) 16.0000000000000

        sage: R.<x>=ZZ[]
        sage: list(multiples(x,5))
        [0, x, 2*x, 3*x, 4*x]
        sage: list(multiples(x,5,operation='*'))
        [1, x, x^2, x^3, x^4]
        sage: list(multiples(x,5,indexed=True))
        [(0, 0), (1, x), (2, 2*x), (3, 3*x), (4, 4*x)]
        sage: list(multiples(x,5,indexed=True,operation='*'))
        [(0, 1), (1, x), (2, x^2), (3, x^3), (4, x^4)]
        sage: for i,y in multiples(x,5,indexed=True): print "%s  times %s = %s"%(i,x,y)
        0  times x = 0
        1  times x = x
        2  times x = 2*x
        3  times x = 3*x
        4  times x = 4*x

        sage: for i,n in multiples(3,5,indexed=True,operation='*'):  print "3 to the power %s = %s"%(i,n)
        3 to the power 0 = 1
        3 to the power 1 = 3
        3 to the power 2 = 9
        3 to the power 3 = 27
        3 to the power 4 = 81
    """
    def __init__(self,P,n,P0=None,indexed=False, operation='+', op=None):
        """
        Create a multiples iterator

        INPUT:

            P -- step value: any Sage object on which a binary
                             operation is defined

            n -- number of multiples: non-negative integer
            P0 -- offset (default 0): Sage object which can be 'added' to P
            indexed -- boolean (default False)

            operation -- string: '+' (default ) or
                                '*' or other.

                         If other, a function op() must be supplied (a
                         function of 2 arguments) defining the
                         group binary operation; also either P0 must be supplied.

            If indexed==False the iterator delivers P0+i*P (if
            operation=='+') or P0*P**i (if
            operation=='*'), for i in range(n) .

            If indexed==True the iterator delivers tuples (i,P0+i*P)
            or (i,P0*P**i)
        """
        if n<0:
            raise ValueError, 'n cannot be negative in multiples'

        from operator import mul, add

        if operation in multiplication_names:
            if P0==None: P0 = P.parent()(1)
            self.op = mul
        elif operation in addition_names:
            if P0==None: P0 = P.parent()(0)
            self.op = add
        else:
            self.op = op
            if P0==None:
                raise ValueError, "P0 must be supplied when operation is neither addition nor multiplication"
            if op==None:
                raise ValueError, "op() must both be supplied when operation is neither addition nor multiplication"

        self.P=copy(P)
        self.Q=copy(P0)
        assert self.P != None and self.Q!=None
        self.i = 0
        self.bound = n
        self.indexed = indexed

    def next(self):
        if self.i >= self.bound:
            raise StopIteration
        i = self.i
        val = self.Q
        self.i +=1
        self.Q=self.op(self.Q,self.P)
        if self.indexed:
            return (i,val)
        else:
            return val
    def __iter__(self):
        return self


def bsgs(a, b, bounds, operation='*', identity=None, inverse=None, op=None):
    r"""
    Totally generic discrete baby-step giant-step function.

    Solves n*a=b (or a^n=b) with lb<=n<=ub where bounds==(lb,ub),
    raising an error if no such n exists.

    a and b must be elements of some group with given identity,
    inverse of x given by inverse(x), and group operation on x,y by
    op(x,y).

    If operation is '*' or '+' then the other
    arguments are provided automatically; otherwise they must be
    provided by the caller.

    INPUT:
        a    -- group element
        b    -- group element
        bounds -- a 2-tuple of integers (lower,upper) with 0<=lower<=upper
        operation -- string: '*', '+', 'other'
        identity -- the identity element of the group
        inverse()  -- function of 1 argument x returning inverse of x
        op() -- function of 2 arguments x,y returning x*y in group

    OUTPUT:
        Returns an integer $n$ such that $a^n = b$ (or $n*a = b$).
        If no such $n$ exists, this function raises a ValueError exception.

    NOTE: This is a generalization of discrete logarithm.  One
        situation where this version is useful is to find the order of
        an element in a group where we only have bounds on the group
        order (see the elliptic curve example below).

    ALGORITHM: Baby step giant step.  Time and space are soft
        O(sqrt(n)) where n is the difference between upper and lower
        bounds.

    EXAMPLES:
        sage: b = Mod(2,37);  a = b^20
        sage: bsgs(b, a, (0,36))
        20

        sage: p=next_prime(10^20)
        sage: a=Mod(2,p); b=a^(10^25)
        sage: bsgs(a, b, (10^25-10^6,10^25+10^6)) == 10^25
        True

        sage: K = GF(3^6,'b')
        sage: a = K.gen()
        sage: b = a^210
        sage: bsgs(a, b, (0,K.order()-1))
        210

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^500
        sage: bsgs(z,w,(0,229))
        40

        An additive example in an elliptic curve group:
        sage: F.<a>=GF(37^5,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: P=E.lift_x(a); P
        (a : 9*a^4 + 22*a^3 + 23*a^2 + 30 : 1)

        This will return a multiple of the order of P:
        sage: bsgs(P,P.parent()(0),Hasse_bounds(F.order()),operation='+')
        69327408

    AUTHOR:
        -- John Cremona (2008-03-15)
    """
    Z = integer_ring.ZZ

    from operator import inv, mul, neg, add

    if operation in multiplication_names:
        identity = a.parent()(1)
        inverse  = inv
        op = mul
    elif operation in addition_names:
        identity = a.parent()(0)
        inverse  = neg
        op = add
    else:
        if identity==None or inverse==None or op==None:
            raise ValueError, "identity, inverse and operation must be given"

    lb, ub = bounds
    if lb<0 or ub<lb:
        raise ValueError, "bsgs() requires 0<=lb<=ub"

    if a.is_zero() and not b.is_zero():
        raise ValueError, "No solution in bsgs()"

    ran = 1 + ub - lb   # the length of the interval

    c = op(inverse(b),multiple(a,lb,operation=operation))

    if ran < 30:    # use simple search for small ranges
        i = lb
        d = c
#        for i,d in multiples(a,ran,c,indexed=True,operation=operation):
        for i0 in range(ran):
            i = lb + i0
            if identity == d:        # identity == b^(-1)*a^i, so return i
                return Z(i)
            d = op(a,d)
        raise ValueError, "No solution in bsgs()"

    m = ran.isqrt()+1  # we need sqrt(ran) rounded up
    table = dict()     # will hold pairs (a^(lb+i),lb+i) for i in range(m)

    d=c
    for i0 in misc.srange(m):
        i = lb + i0
        if identity==d:        # identity == b^(-1)*a^i, so return i
            return Z(i)
        table[d] = i
        d=op(d,a)

    c = op(c,inverse(d))     # this is now a**(-m)
    d=identity
    for i in misc.srange(m):
        j = table.get(d)
        if not j==None:  # then d == b*a**(-i*m) == a**j
            return Z(i*m + j)
        d=op(c,d)

    raise ValueError, "Log of %s to the base %s does not exist in %s."%(b,a,bounds)

def discrete_log(a, base, ord=None, operation='*', identity=None, inverse=None, op=None):
    """
    Totally generic discrete log function.

    a and base must be elements of some group with identity given by
    identity, inverse of x by inverse(x), and group operation on x,y
    by op(x,y).

    If operation is '*' or '+' then the other
    arguments are provided automatically; otherwise they must be
    provided by the caller.

    WARNING: If x has a log method, it is likely to be vastly faster
    than using this function.  E.g., if x is an integer modulo n, use
    its log method instead!

    INPUT:
        a    -- group element
        base -- group element (the base)
        ord  -- integer (multiple of order of base, or None)
        operation -- string: '*', '+', 'other'
        identity -- the group's identity
        inverse()  -- function of 1 argument x returning inverse of x
        op() -- function of 2 arguments x,y returning x*y in group

    OUTPUT:
        Returns an integer $n$ such that $b^n = a$ (or $n*b = a$),
        assuming that ord is a multiple of the order of the base $b$.
        If ord is not specified an attempt is made to compute it.

        If no such $n$ exists, this function raises a ValueError exception.

    ALGORITHM: Baby step giant step.

    EXAMPLES:
        sage: b = Mod(2,37);  a = b^20
        sage: discrete_log(a, b)
        20
        sage: b = Mod(2,997);  a = b^20
        sage: discrete_log(a, b)
        20

        sage: K = GF(3^6,'b')
        sage: b = K.gen()
        sage: a = b^210
        sage: discrete_log(a, b, K.order()-1)
        210

        sage: b = Mod(1,37);  x = Mod(2,37)
        sage: discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of 2 found to base 1
        sage: b = Mod(1,997);  x = Mod(2,997)
        sage: discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of 2 found to base 1

        See trac\#2356:
        sage: F.<w> = GF(121)
        sage: v = w^120
        sage: v.log(w)
        0

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^50
        sage: discrete_log(w,z)
        50

        An example where the order is infinite: note that we must give
        an upper bound here:
        sage: K.<a> = QuadraticField(23)
        sage: eps = 5*a-24        # a fundamental unit
        sage: eps.multiplicative_order()
        +Infinity
        sage: eta = eps^100
        sage: discrete_log(eta,eps,1000)
        100

        In this case we cannot detect negative powers:
        sage: eta = eps^(-3)
        sage: discrete_log(eta,eps,100)
        Traceback (most recent call last):
        ...
        ValueError: No discrete log of -11515*a - 55224 found to base 5*a - 24

        But we can invert the base (and negate the result) instead:
        sage: - discrete_log(eta^-1,eps,100)
        -3

        An additive example: elliptic curve DLOG:
        sage: F=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: F.<a>=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: P=E(25*a + 16 , 15*a + 7 )
        sage: P.order()
        672
        sage: Q=39*P; Q
        (36*a + 32 : 5*a + 12 : 1)
        sage: discrete_log(Q,P,P.order(),'+')
        39

    AUTHOR:
        -- William Stein and David Joyner (2005-01-05)
        -- John Cremona (2008-02-29) rewrite using dict() and make generic

    """
    if ord == None:
        if operation in multiplication_names:
            try:
                ord = base.multiplicative_order()
            except:
                ord = base.order()
        elif operation in addition_names:
            try:
                ord = base.additive_order()
            except:
                ord = base.order()
        else:
            try:
                ord = base.order()
            except:
                raise ValueError, "ord must be specified"
    try:
        return bsgs(base,a,(0,ord),operation=operation)
    except ValueError:
        raise ValueError, "No discrete log of %s found to base %s"%(a,base)

################################################################
#
# Older version of discrete_log.  Works fine but has been
# superceded by the version which simply calls the more general
# bsgs() function.
#
################################################################

def old_discrete_log(a, base, ord=None, operation='*',
                         identity=None, inverse=None, op=None):
    r"""
    Totally generic discrete log function.

    a and base must be elements of some group with identity given by
    identity, inverse of x by inverse(x), and group operation on x,y
    by op(x,y).

    If operation is '*' or '+' then the other
    arguments are provided automatically; otherwise they must be
    provided by the caller.

    WARNING: If x has a log method, it is likely to be vastly faster
    than using this function.  E.g., if x is an integer modulo n, use
    its log method instead!

    INPUT:
        a    -- group element
        base -- group element (the base)
        ord  -- integer (multiple of order of base, or None)
        operation -- string: '*', '+', 'other'
        identity -- the group's identity
        inverse()  -- function of 1 argument x returning inverse of x
        op() -- function of 2 arguments x,y returning x*y in group

    OUTPUT:
        Returns an integer $n$ such that $b^n = a$ (or $n*b = a$),
        assuming that ord is a multiple of the order of the base $b$.
        If ord is not specified an attempt is made to compute it.

        If no such $n$ exists, this function raises a ValueError exception.

    ALGORITHM: Baby step giant step.

    EXAMPLES:
        sage: b = Mod(2,37);  a = b^20
        sage: old_discrete_log(a, b)
        20
        sage: b = Mod(2,997);  a = b^20
        sage: old_discrete_log(a, b)
        20

        sage: K = GF(3^6,'b')
        sage: b = K.gen()
        sage: a = b^210
        sage: old_discrete_log(a, b, K.order()-1)
        210

        sage: b = Mod(1,37);  x = Mod(2,37)
        sage: old_discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: Log of 2 to the base 1 does not exist.
        sage: b = Mod(1,997);  x = Mod(2,997)
        sage: old_discrete_log(x, b)
        Traceback (most recent call last):
        ...
        ValueError: Log of 2 to the base 1 does not exist.

        See trac\#2356:
        sage: F.<w> = GF(121)
        sage: v = w^120
        sage: v.log(w)
        0

        sage: K.<z>=CyclotomicField(230)
        sage: w=z^50
        sage: old_discrete_log(w,z)
        50

        An example where the order is infinite: note that we must give
        an upper bound here:
        sage: K.<a> = QuadraticField(23)
        sage: eps = 5*a-24        # a fundamental unit
        sage: eps.multiplicative_order()
        +Infinity
        sage: eta = eps^100
        sage: old_discrete_log(eta,eps,1000)
        100

        In this case we cannot detect negative powers:
        sage: eta = eps^(-3)
        sage: old_discrete_log(eta,eps,100)
        Traceback (most recent call last):
        ...
        ValueError: Log of -11515*a - 55224 to the base 5*a - 24 does not exist.

        But we can invert the base (and negate the result) instead:
        sage: - old_discrete_log(eta^-1,eps,100)
        -3

        An additive example: elliptic curve DLOG:
        sage: F=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: F.<a>=GF(37^2,'a')
        sage: E=EllipticCurve(F,[1,1])
        sage: P=E(25*a + 16 , 15*a + 7 )
        sage: P.order()
        672
        sage: Q=39*P; Q
        (36*a + 32 : 5*a + 12 : 1)
        sage: old_discrete_log(Q,P,P.order(),'+')
        39

    AUTHOR:
        -- William Stein and David Joyner (2005-01-05)
        -- John Cremona (2008-02-29) rewrite using dict() and make generic
    """
    Z = integer_ring.ZZ
    b = base

    from operator import inv, mul, neg, add

    if operation in multiplication_names:
        identity = b.parent()(1)
        inverse  = inv
        op = mul
        if ord==None:
            ord = b.multiplicative_order()
    elif operation in addition_names:
        identity = b.parent()(0)
        inverse  = neg
        op = add
        if ord==None:
            ord = b.order()
    else:
        if ord==None or identity==None or inverse==None or op==None:
            raise ValueError, "order, identity, inverse and operation must all be specified"

    ord = Z(ord)
    if ord < 100:
        c = identity
        for i in range(ord):
            if c == a:        # is b^i
                return Z(i)
            c = op(c,b)
        raise ValueError, "Log of %s to the base %s does not exist."%(a,b)

    m = ord.isqrt()+1  # we need sqrt(ord) rounded up
    table = dict()     # will hold pairs (b^j,j) for j in range(m)
    g = identity       # will run through b**j
    for j in range(m):
        if a==g:
            return Z(j)
        table[g] = j
        g = op(g,b)

    g = inverse(g)     # this is now b**(-m)
    h = op(a,g)        # will run through a*g**i = a*b**(-i*m)
    for i in range(1,m):
        j = table.get(h)
        if not j==None:  # then a*b**(-i*m) == b**j
            return Z(i*m + j)
        if i < m-1:
            h = op(h,g)

    raise ValueError, "Log of %s to the base %s does not exist."%(a,b)

discrete_log_generic = discrete_log


################################################################
#
# Generic linear relation finder
#
################################################################

def linear_relation(P, Q, operation='+', identity=None, inverse=None, op=None):
    """
    Function which solves the equation a*P=m*Q or P^a=Q^m.

    Additive version: returns (a,m) with minimal m>0 such that
    a*P==m*Q.  Special case: if <P> and <Q> intersect only in {0} then
    (a,m)=(0,Q.additive_order()).

    Multiplicative version: returns (a,m) with minimal m>0 such that
    P^a==Q^m.  Special case: if <P> and <Q> intersect only in {1} then
    (a,m)=(0,Q.multiplicative_order()).

    Works in general finite abelian groups:  uses bsgs()

    EXAMPLE:

    An additive example (in an elliptic curve group):
        sage: F.<a>=GF(3^6,'a')
        sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
        sage: P=E(a^5 + a^4 + a^3 + a^2 + a + 2 , 0)
        sage: Q=E(2*a^3 + 2*a^2 + 2*a , a^3 + 2*a^2 + 1)
        sage: linear_relation(P,Q,'+')
        (1, 2)
        sage: P == 2*Q
        True

    An multiplcative example (in a finite field's multiplicative group):
        sage: F.<a>=GF(3^6,'a')
        sage: a.multiplicative_order().factor()
        2^3 * 7 * 13
        sage: b=a^7
        sage: c=a^13
        sage: linear_relation(b,c,'*')
        (13, 7)
        sage: b^13==c^7
        True
    """

    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        op = mul
        try:
            n = P.multiplicative_order()
            m = Q.multiplicative_order()
        except:
            n = P.order()
            m = Q.order()
    elif operation in addition_names:
        op = add
        try:
            n = P.additive_order()
            m = Q.additive_order()
        except:
            n = P.order()
            m = Q.order()
    else:
        if op==None:
            raise ValueError, "operation must be specified"
        n = P.order()
        m = Q.order()

    g = sage.rings.arith.gcd(n,m)
    if g==1: return (m,Z(0))
    n1 = n//g
    m1 = m//g
    P1 = multiple(P,n1,operation=operation)  # has exact order g
    Q1 = multiple(Q,m1,operation=operation)  # has exact order g

    # now see if Q1 is a multiple of P1; the only multiples we
    # need check are h*Q1 where h divides g
    for h in g.divisors(): # positive divisors!
        try:
            Q2 = multiple(Q1,h,operation=operation)
            return (n1 * bsgs(P1,Q2,(0,g-1),operation=operation),
                    m1 * h)
        except ValueError:
            pass # to next h
    raise ValueError, "No solution found in linear_relation!"

################################################################
#
# Generic functions to find orders of elements
#
# 1. order_from_multiple: finds the order given a multiple of the order
#
# 2. order_from_bounds: finds the order given an interval containing a
# multiple of the order
#
################################################################

def order_from_multiple(P, m, plist=None, operation='+',
                         identity=None, inverse=None, op=None):
    """
    Generic function to find order of P given a multiple of the order

        INPUT:

            P -- a Sage object which is a group element

            m -- a Sage integer which is a multiple of the order of P,
            i.e. we require that m*P=0 (or P^m=1).

            plist -- a list of the prime factors of m, or None in
            which case this function will need to factor m.

            operation -- string: '+' (default ) or
                                 '*' or other.

                         If other, the following must be supplied:

                         identity: the identity element for the group;

                         inverse(): a function of one argument giving
                         the inverse of a group element;

                         op(): a function of 2 arguments defining
                         the group binary operation.

        NOTE: It is more efficient for the caller to factor m and cache the
        factors for subsequent calls.

        EXAMPLES:
            sage: k.<a> = GF(5^5)
            sage: b = a^4
            sage: order_from_multiple(b,5^5-1,operation='*')
            781
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E(3*a^4 + 3*a , 2*a + 1 )
            sage: M = E.cardinality(); M
            3227
            sage: plist = M.prime_factors()
            sage: order_from_multiple(P, M, plist, operation='+')
            3227
            sage: Q = E(0,2)
            sage: order_from_multiple(Q, M, plist, operation='+')
            7

            sage: K.<z>=CyclotomicField(230)
            sage: w=z^50
            sage: order_from_multiple(w,230,operation='*')
            23

    """
    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        op = mul
        identity = P.parent()(1)
    elif operation in addition_names:
        op = add
        identity = P.parent()(0)
    else:
        if op==None:
            raise ValueError, "operation and identity must be specified"

    M=Z(m)
    assert multiple(P,M,operation=operation) == identity

    if plist==None:
        plist = M.prime_factors()

    # For each p in plist we determine the power of p dividing
    # the order, accumulating the order in N

    N=Z(1)
    for p in plist:
        Q = multiple(P,M.prime_to_m_part(p),operation=operation)
        # so Q has p-power order
        while Q != identity:
            Q = multiple(Q,p,operation=operation)
            N *= p

    # now N is the exact order of self
    return N

def order_from_bounds(P, bounds, d=None, operation='+',
                         identity=None, inverse=None, op=None):
    """
    Generic function to find order of P given bounds on group order.

    upper and lower bounds for a multiple of the order (e.g. bounds on the
    order of the group of which P is an element)

        INPUT:

            P      -- a Sage object which is a group element

            bounds -- a 2-tuple (lb,ub) such that m*P=0 (or P^m=1) for
                      some m with lb<=m<=ub

            d      -- (optional) a positive integer; only m which are
            multiples of this will be considered.

            operation -- string: '+' (default ) or
                                 '*' or other.

                         If other, the following must be supplied:

                         identity: the identity element for the group;

                         inverse(): a function of one argument giving
                         the inverse of a group element;

                         op(): a function of 2 arguments defining
                         the group binary operation.


        NOTE: Typically lb and ub will be bounds on the group order,
        and from previous calculation we know that the group order is
        divisible by d.

        EXAMPLES:
            sage: k.<a> = GF(5^5)
            sage: b = a^4
            sage: order_from_bounds(b,(5^4,5^5),operation='*')
            781
            sage: E = EllipticCurve(k,[2,4])
            sage: P = E(3*a^4 + 3*a , 2*a + 1 )
            sage: bounds = Hasse_bounds(5^5)
            sage: Q = E(0,2)
            sage: order_from_bounds(Q, bounds, operation='+')
            7
            sage: order_from_bounds(P, bounds, 7, operation='+')
            3227

            sage: K.<z>=CyclotomicField(230)
            sage: w=z^50
            sage: order_from_bounds(w,(200,250),operation='*')
            23

    """
    from operator import mul, add
    Z = integer_ring.ZZ

    if operation in multiplication_names:
        op = mul
        identity = P.parent()(1)
    elif operation in addition_names:
        op = add
        identity = P.parent()(0)
    else:
        if op==None:
            raise ValueError, "operation and identity must be specified"

    Q = P
    if d == None: d = 1
    if d > 1:
        Q = multiple(P,d,operation=operation)
        lb, ub = bounds
        bounds = ( sage.rings.arith.integer_ceil(lb/d),
                   sage.rings.arith.integer_floor(ub/d) )

    # Use generic bsgs to find  n=d*m with lb<=n<=ub and n*P=0

    m = d * bsgs(Q, identity, bounds, operation=operation)

    # Now use the order_from_multiple() function to finish the job:

    return order_from_multiple(P, m, operation=operation)

def merge_points(P1,P2, operation='+',
                         identity=None, inverse=None, op=None):
    """
    Returns a group element whose order is the lcm of the given elements

    INPUT:
        P1 -- a pair (g1,n1) where g1 is a group element of order n1
        P2 -- a pair (g2,n2) where g2 is a group element of order n2
        operation -- string: '+' (default ) or
                             '*' or other.

                         If other, the following must be supplied:

                         identity: the identity element for the group;

                         inverse(): a function of one argument giving
                         the inverse of a group element;

                         op(): a function of 2 arguments defining
                         the group binary operation.


    OUTPUT:
         A pair (g3,n3) where g3 has order n3=lcm(n1,n2)

    EXAMPLES:
        sage: F.<a>=GF(3^6,'a')
        sage: b = a^7
        sage: c = a^13
        sage: ob = (3^6-1)//7
        sage: oc = (3^6-1)//13
        sage: merge_points((b,ob),(c,oc),operation='*')
        (a^4 + 2*a^3 + 2*a^2, 728)
        sage: d,od = merge_points((b,ob),(c,oc),operation='*')
        sage: od == d.multiplicative_order()
        True
        sage: od == lcm(ob,oc)
        True

        sage: E=EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1])
        sage: P=E(2*a^5 + 2*a^4 + a^3 + 2 , a^4 + a^3 + a^2 + 2*a + 2)
        sage: P.order()
        7
        sage: Q=E(2*a^5 + 2*a^4 + 1 , a^5 + 2*a^3 + 2*a + 2 )
        sage: Q.order()
        4
        sage: R,m = merge_points((P,7),(Q,4), operation='+')
        sage: R.order() == m
        True
        sage: m == lcm(7,4)
        True
    """
    from operator import mul, add
    Z = integer_ring.ZZ

    g1, n1 = P1
    g2, n2 = P2

    if operation in multiplication_names:
        op = mul
        identity = g1.parent()(1)
    elif operation in addition_names:
        op = add
        identity = g1.parent()(0)
    else:
        if op==None:
            raise ValueError, "operation and identity must be specified"

    assert multiple(g1,n1,operation=operation) == identity
    assert multiple(g2,n2,operation=operation) == identity

    # trivial cases
    if n1.divides(n2):
        return (g2,n2)
    if n2.divides(n1):
        return (g1,n1)

    m,k1,k2 = sage.rings.arith.xlcm(n1,n2);
    m1 = n1//k1
    m2 = n2//k2
    g1 = multiple(g1,m1,operation=operation)
    g2 = multiple(g2,m2,operation=operation)
    return (op(g1,g2), m)

