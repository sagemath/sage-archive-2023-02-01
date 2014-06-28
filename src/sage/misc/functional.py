# -*- coding: utf-8 -*-
"""
Functional notation

These are functions so that you can write foo(x) instead of x.foo()
in certain common cases.

AUTHORS:

- William Stein: Initial version

- David Joyner (2005-12-20): More Examples
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


import sage.misc.latex
import sage.interfaces.expect
import sage.interfaces.mathematica


from sage.rings.complex_double import CDF
from sage.rings.real_double import RDF, RealDoubleElement

import sage.rings.real_mpfr
import sage.rings.complex_field
import sage.rings.integer

import __builtin__

LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363 # a small overestimate of log(10,2)

##############################################################################
# There are many functions on elements of a ring, which mathematicians
# usually write f(x), e.g., it is weird to write x.log() and natural
# to write log(x).  The functions below allow for the more familiar syntax.
##############################################################################
def additive_order(x):
    """
    Returns the additive order of `x`.

    EXAMPLES::

        sage: additive_order(5)
        +Infinity
        sage: additive_order(Mod(5,11))
        11
        sage: additive_order(Mod(4,12))
        3
    """
    return x.additive_order()

def base_ring(x):
    """
    Returns the base ring over which x is defined.

    EXAMPLES::

        sage: R = PolynomialRing(GF(7), 'x')
        sage: base_ring(R)
        Finite Field of size 7
    """
    return x.base_ring()

def base_field(x):
    """
    Returns the base field over which x is defined.

    EXAMPLES::

        sage: R = PolynomialRing(GF(7), 'x')
        sage: base_ring(R)
        Finite Field of size 7
        sage: base_field(R)
        Finite Field of size 7

    This catches base rings which are fields as well, but does
    not implement a ``base_field`` method for objects which do
    not have one::

        sage: R.base_field()
        Traceback (most recent call last):
        ...
        AttributeError: 'PolynomialRing_dense_mod_p_with_category' object has no attribute 'base_field'
    """
    try:
        return x.base_field()
    except AttributeError:
        y = x.base_ring()
        if is_field(y):
            return y
        else:
            raise AttributeError("The base ring of %s is not a field"%x)

def basis(x):
    """
    Returns the fixed basis of x.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: S = V.subspace([[1,2,0],[2,2,-1]])
        sage: basis(S)
        [
        (1, 0, -1),
        (0, 1, 1/2)
        ]
    """
    return x.basis()

def category(x):
    """
    Returns the category of x.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: category(V)
        Category of vector spaces over Rational Field
    """
    try:
        return x.category()
    except AttributeError:
        import sage.categories.all
        return sage.categories.all.Objects()

def ceil(x):
    """
    Returns the ceiling (least integer) function of x.

    EXAMPLES::

        sage: ceil(3.5)
        4
        sage: ceil(7/2)
        4
        sage: ceil(-3.5)
        -3
        sage: ceil(RIF(1.3,2.3))
        3.?
    """
    try:
        return x.ceil()
    except AttributeError:
        return sage.rings.all.ceil(x)

def characteristic_polynomial(x, var='x'):
    """
    Returns the characteristic polynomial of x in the given variable.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: charpoly(A)
        x^3 - 15*x^2 - 18*x
        sage: charpoly(A, 't')
        t^3 - 15*t^2 - 18*t

    ::

        sage: k.<alpha> = GF(7^10); k
        Finite Field in alpha of size 7^10
        sage: alpha.charpoly('T')
        T^10 + T^6 + T^5 + 4*T^4 + T^3 + 2*T^2 + 3*T + 3
        sage: characteristic_polynomial(alpha, 'T')
        T^10 + T^6 + T^5 + 4*T^4 + T^3 + 2*T^2 + 3*T + 3

    Ensure the variable name of the polynomial does not conflict with
    variables used within the matrix, and that non-integral powers of
    variables don't confuse the computation (:trac:`14403`)::

        sage: y = var('y')
        sage: a = matrix([[x,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        sage: characteristic_polynomial(a).list()
        [x, -3*x - 1, 3*x + 3, -x - 3, 1]
        sage: b = matrix([[y^(1/2),0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        sage: charpoly(b).list()
        [sqrt(y), -3*sqrt(y) - 1, 3*sqrt(y) + 3, -sqrt(y) - 3, 1]
    """
    try:
        return x.charpoly(var)
    except AttributeError:
        raise NotImplementedError("computation of charpoly of x (={}) not implemented".format(x))

charpoly = characteristic_polynomial


def coerce(P, x):
    """
    Attempts to coerce x to type P if possible.

    EXAMPLES::

        sage: type(5)
        <type 'sage.rings.integer.Integer'>
        sage: type(coerce(QQ,5))
        <type 'sage.rings.rational.Rational'>
    """
    try:
        return P._coerce_(x)
    except AttributeError:
        return P(x)


def acos(x):
    """
    Returns the arc cosine of x.

    EXAMPLES::

        sage: acos(.5)
        1.04719755119660
        sage: acos(sin(pi/3))
        arccos(1/2*sqrt(3))
        sage: acos(sin(pi/3)).simplify_full()
        1/6*pi
    """
    try: return x.acos()
    except AttributeError: return RDF(x).acos()

def asin(x):
    """
    Returns the arc sine of x.

    EXAMPLES::

        sage: asin(.5)
        0.523598775598299
        sage: asin(sin(pi/3))
        arcsin(1/2*sqrt(3))
        sage: asin(sin(pi/3)).simplify_full()
        1/3*pi
    """
    try: return x.asin()
    except AttributeError: return RDF(x).asin()

def atan(x):
    """
    Returns the arc tangent of x.

    EXAMPLES::

        sage: z = atan(3);z
        arctan(3)
        sage: n(z)
        1.24904577239825
        sage: atan(tan(pi/4))
        1/4*pi
    """
    try: return x.atan()
    except AttributeError: return RDF(x).atan()

## def cuspidal_submodule(x):
##     return x.cuspidal_submodule()

## def cuspidal_subspace(x):
##     return x.cuspidal_subspace()

def cyclotomic_polynomial(n, var='x'):
    """
    Returns the `n^{th}` cyclotomic polynomial.

    EXAMPLES::

        sage: cyclotomic_polynomial(3)
        x^2 + x + 1
        sage: cyclotomic_polynomial(4)
        x^2 + 1
        sage: cyclotomic_polynomial(9)
        x^6 + x^3 + 1
        sage: cyclotomic_polynomial(10)
        x^4 - x^3 + x^2 - x + 1
        sage: cyclotomic_polynomial(11)
        x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
    """
    return sage.rings.all.ZZ[var].cyclotomic_polynomial(n)

def decomposition(x):
    """
    Returns the decomposition of x.

    EXAMPLES::

        sage: M = matrix([[2, 3], [3, 4]])
        sage: M.decomposition()
        [
        (Ambient free module of rank 2 over the principal ideal domain Integer Ring, True)
        ]

    ::

        sage: G.<a,b> = DirichletGroup(20)
        sage: c = a*b
        sage: d = c.decomposition(); d
        [Dirichlet character modulo 4 of conductor 4 mapping 3 |--> -1,
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4]
        sage: d[0].parent()
        Group of Dirichlet characters of modulus 4 over Cyclotomic Field of order 4 and degree 2
    """
    return x.decomposition()

def denominator(x):
    """
    Returns the denominator of x.

    EXAMPLES::

        sage: denominator(17/11111)
        11111
        sage: R.<x> = PolynomialRing(QQ)
        sage: F = FractionField(R)
        sage: r = (x+1)/(x-1)
        sage: denominator(r)
        x - 1
    """
    if isinstance(x, (int, long)):
        return 1
    return x.denominator()

def det(x):
    """
    Returns the determinant of x.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: det(A)
        0
    """
    return x.det()

def dimension(x):
    """
    Returns the dimension of x.

    EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: S = V.subspace([[1,2,0],[2,2,-1]])
        sage: dimension(S)
        2
    """
    return x.dimension()

dim = dimension

def discriminant(x):
    """
    Returns the discriminant of x.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: S = R.quotient(x^29 - 17*x - 1, 'alpha')
        sage: K = S.number_field()
        sage: discriminant(K)
        -15975100446626038280218213241591829458737190477345113376757479850566957249523
    """
    return x.discriminant()

disc = discriminant

# This is dangerous since it gets the scoping all wrong ??
#import __builtin__
#def eval(x):
#    try:
#        return x._eval_()
#    except AttributeError:
#        return __builtin__.eval(x)

def eta(x):
    r"""
    Returns the value of the eta function at `x`, which must be
    in the upper half plane.

    The `\eta` function is

    .. math::

                    \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty}(1-e^{2\pi inz})



    EXAMPLES::

        sage: eta(1+I)
        0.742048775837 + 0.19883137023*I
    """
    try: return x.eta()
    except AttributeError: return CDF(x).eta()

def exp(x):
    """
    Returns the value of the exponentiation function at x.

    EXAMPLES::

        sage: exp(3)
        e^3
        sage: exp(0)
        1
        sage: exp(2.5)
        12.1824939607035
        sage: exp(pi*i)
        -1
    """
    try: return x.exp()
    except AttributeError: return RDF(x).exp()

def factor(x, *args, **kwds):
    """
    Returns the (prime) factorization of x.

    EXAMPLES::

        sage: factor(factorial(10))
        2^8 * 3^4 * 5^2 * 7
        sage: n = next_prime(10^6); n
        1000003
        sage: factor(n)
        1000003

        Note that this depends on the type of x::

        sage: factor(55)
        5 * 11
        sage: factor(x^2+2*x+1)
        (x + 1)^2
        sage: factor(55*x^2+110*x+55)
        55*(x + 1)^2

    """
    try: return x.factor(*args, **kwds)
    except AttributeError: return sage.rings.all.factor(x, *args, **kwds)

factorization = factor
factorisation = factor

def fcp(x, var='x'):
    """
    Returns the factorization of the characteristic polynomial of x.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: fcp(A, 'x')
        x * (x^2 - 15*x - 18)
    """
    try: return x.fcp(var)
    except AttributeError: return factor(charpoly(x, var))

## def floor(x):
##     try:
##         return x.floor()
##     except AttributeError:
##         return sage.rings.all.floor(x)

def gen(x):
    """
    Returns the generator of x.

    EXAMPLES::

        sage: R.<x> = QQ[]; R
        Univariate Polynomial Ring in x over Rational Field
        sage: gen(R)
        x
        sage: gen(GF(7))
        1
        sage: A = AbelianGroup(1, [23])
        sage: gen(A)
        f
    """
    return x.gen()

def gens(x):
    """
    Returns the generators of x.

    EXAMPLES::

        sage: R.<x,y> = SR[]
        sage: R
        Multivariate Polynomial Ring in x, y over Symbolic Ring
        sage: gens(R)
        (x, y)
        sage: A = AbelianGroup(5, [5,5,7,8,9])
        sage: gens(A)
        (f0, f1, f2, f3, f4)
    """
    return x.gens()

def hecke_operator(x,n):
    """
    Returns the n-th Hecke operator T_n acting on x.

    EXAMPLES::

        sage: M = ModularSymbols(1,12)
        sage: hecke_operator(M,5)
        Hecke operator T_5 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
    """
    return x.hecke_operator(n)

ideal = sage.rings.ideal.Ideal

def image(x):
    """
    Returns the image of x.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: image(A)
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]
    """
    return x.image()

def symbolic_sum(expression, *args, **kwds):
    r"""
    Returns the symbolic sum `\sum_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` - a symbolic expression

    - ``v`` - a variable or variable name

    - ``a`` - lower endpoint of the sum

    - ``b`` - upper endpoint of the sum

    - ``algorithm`` - (default: 'maxima')  one of
      - 'maxima' - use Maxima (the default)
      - 'maple' - (optional) use Maple
      - 'mathematica' - (optional) use Mathematica

    EXAMPLES::

        sage: k, n = var('k,n')
        sage: sum(k, k, 1, n).factor()
        1/2*(n + 1)*n

    ::

        sage: sum(1/k^4, k, 1, oo)
        1/90*pi^4

    ::

        sage: sum(1/k^5, k, 1, oo)
        zeta(5)

    A well known binomial identity::

        sage: sum(binomial(n,k), k, 0, n)
        2^n

    The binomial theorem::

        sage: x, y = var('x, y')
        sage: sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
        (x + y)^n

    ::

        sage: sum(k * binomial(n, k), k, 1, n)
        2^(n - 1)*n

    ::

        sage: sum((-1)^k*binomial(n,k), k, 0, n)
        0

    ::

        sage: sum(2^(-k)/(k*(k+1)), k, 1, oo)
        -log(2) + 1

    Another binomial identity (trac #7952)::

        sage: t,k,i = var('t,k,i')
        sage: sum(binomial(i+t,t),i,0,k)
        binomial(k + t + 1, t + 1)

    Summing a hypergeometric term::

        sage: sum(binomial(n, k) * factorial(k) / factorial(n+1+k), k, 0, n)
        1/2*sqrt(pi)/factorial(n + 1/2)

    We check a well known identity::

        sage: bool(sum(k^3, k, 1, n) == sum(k, k, 1, n)^2)
        True

    A geometric sum::

        sage: a, q = var('a, q')
        sage: sum(a*q^k, k, 0, n)
        (a*q^(n + 1) - a)/(q - 1)

    The geometric series::

        sage: assume(abs(q) < 1)
        sage: sum(a*q^k, k, 0, oo)
        -a/(q - 1)

    A divergent geometric series.  Don't forget
    to forget your assumptions::

        sage: forget()
        sage: assume(q > 1)
        sage: sum(a*q^k, k, 0, oo)
        Traceback (most recent call last):
        ...
        ValueError: Sum is divergent.

    This summation only Mathematica can perform::

        sage: sum(1/(1+k^2), k, -oo, oo, algorithm = 'mathematica')     # optional - mathematica
        pi*coth(pi)

    Use Maple as a backend for summation::

        sage: sum(binomial(n,k)*x^k, k, 0, n, algorithm = 'maple')      # optional - maple
        (x + 1)^n

    Python ints should work as limits of summation (trac #9393)::

        sage: sum(x, x, 1r, 5r)
        15

    .. note::

       #. Sage can currently only understand a subset of the output of Maxima, Maple and
          Mathematica, so even if the chosen backend can perform the summation the
          result might not be convertable into a Sage expression.

    """
    if hasattr(expression, 'sum'):
        return expression.sum(*args, **kwds)
    elif len(args) <= 1:
        return sum(expression, *args)
    else:
        from sage.symbolic.ring import SR
        return SR(expression).sum(*args, **kwds)


def integral(x, *args, **kwds):
    """
    Returns an indefinite or definite integral of an object x.

    First call x.integral() and if that fails make an object and
    integrate it using Maxima, maple, etc, as specified by algorithm.

    For symbolic expression calls
    :func:`sage.calculus.calculus.integral` - see this function for
    available options.

    EXAMPLES::

        sage: f = cyclotomic_polynomial(10)
        sage: integral(f)
        1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x

    ::

        sage: integral(sin(x),x)
        -cos(x)

    ::

        sage: y = var('y')
        sage: integral(sin(x),y)
        y*sin(x)

    ::

        sage: integral(sin(x), x, 0, pi/2)
        1
        sage: sin(x).integral(x, 0,pi/2)
        1
        sage: integral(exp(-x), (x, 1, oo))
        e^(-1)

    Numerical approximation::

        sage: h = integral(tan(x)/x, (x, 1, pi/3)); h
        integrate(tan(x)/x, x, 1, 1/3*pi)
        sage: h.n()
        0.07571599101...

    Specific algorithm can be used for integration::

        sage: integral(sin(x)^2, x, algorithm='maxima')
        1/2*x - 1/4*sin(2*x)
        sage: integral(sin(x)^2, x, algorithm='sympy')
        -1/2*cos(x)*sin(x) + 1/2*x

    TESTS:

    A symbolic integral from :trac:`11445` that was incorrect in
    earlier versions of Maxima::

        sage: f = abs(x - 1) + abs(x + 1) - 2*abs(x)
        sage: integrate(f, (x, -Infinity, Infinity))
        2

    Another symbolic integral, from :trac:`11238`, that used to return
    zero incorrectly; with Maxima 5.26.0 one gets
    ``1/2*sqrt(pi)*e^(1/4)``, whereas with 5.29.1, and even more so
    with 5.33.0, the expression is less pleasant, but still has the
    same value.  Unfortunately, the computation takes a very long time
    with the default settings, so we temporarily use the Maxima
    setting ``domain: real``::

        sage: sage.calculus.calculus.maxima('domain: real')
        real
        sage: f = exp(-x) * sinh(sqrt(x))
        sage: t = integrate(f, x, 0, Infinity); t            # long time
        1/4*sqrt(pi)*(erf(1) - 1)*e^(1/4) - 1/4*(sqrt(pi)*(erf(1) - 1) - sqrt(pi) + 2*e^(-1) - 2)*e^(1/4) + 1/4*sqrt(pi)*e^(1/4) - 1/2*e^(1/4) + 1/2*e^(-3/4)
        sage: t.simplify_exp()  # long time
        1/2*sqrt(pi)*e^(1/4)
        sage: sage.calculus.calculus.maxima('domain: complex')
        complex

    An integral which used to return -1 before maxima 5.28. See :trac:`12842`::

        sage: f = e^(-2*x)/sqrt(1-e^(-2*x))
        sage: integrate(f, x, 0, infinity)
        1

    This integral would cause a stack overflow in earlier versions of
    Maxima, crashing sage. See :trac:`12377`. We don't care about the
    result here, just that the computation completes successfully::

        sage: y = (x^2)*exp(x) / (1 + exp(x))^2
        sage: _ = integrate(y, x, -1000, 1000)

    """
    if hasattr(x, 'integral'):
        return x.integral(*args, **kwds)
    else:
        from sage.symbolic.ring import SR
        return SR(x).integral(*args, **kwds)

integrate = integral


def integral_closure(x):
    """
    Returns the integral closure of x.

    EXAMPLES::

        sage: integral_closure(QQ)
        Rational Field
        sage: K.<a> = QuadraticField(5)
        sage: O2 = K.order(2*a); O2
        Order in Number Field in a with defining polynomial x^2 - 5
        sage: integral_closure(O2)
        Maximal Order in Number Field in a with defining polynomial x^2 - 5
    """
    return x.integral_closure()

def interval(a, b):
    r"""
    Integers between a and b *inclusive* (a and b integers).

    EXAMPLES::

        sage: I = interval(1,3)
        sage: 2 in I
        True
        sage: 1 in I
        True
        sage: 4 in I
        False
    """
    return range(a,b+1)

def xinterval(a, b):
    r"""
    Iterator over the integers between a and b, *inclusive*.

    EXAMPLES::

        sage: I = xinterval(2,5); I
        xrange(2, 6)
        sage: 5 in I
        True
        sage: 6 in I
        False
    """
    return xrange(a, b+1)

def is_commutative(x):
    """
    Returns whether or not x is commutative.

    EXAMPLES::

        sage: R = PolynomialRing(QQ, 'x')
        sage: is_commutative(R)
        True
    """
    return x.is_commutative()

def is_even(x):
    """
    Returns whether or not an integer x is even, e.g., divisible by 2.

    EXAMPLES::

        sage: is_even(-1)
        False
        sage: is_even(4)
        True
        sage: is_even(-2)
        True
    """
    try: return x.is_even()
    except AttributeError: return x%2==0

def is_integrally_closed(x):
    """
    Returns whether x is integrally closed.

    EXAMPLES::

        sage: is_integrally_closed(QQ)
        True
        sage: K.<a> = NumberField(x^2 + 189*x + 394)
        sage: R = K.order(2*a)
        sage: is_integrally_closed(R)
        False
    """
    return x.is_integrally_closed()

def is_field(x):
    """
    Returns whether or not x is a field.

    EXAMPLES::

        sage: R = PolynomialRing(QQ, 'x')
        sage: F = FractionField(R)
        sage: is_field(F)
        True
    """
    return x.is_field()

def is_noetherian(x):
    """
    Returns whether or not x is a Noetherian
    object (has ascending chain condition).

    EXAMPLES::

        sage: from sage.misc.functional import is_noetherian
        sage: is_noetherian(ZZ)
        True
        sage: is_noetherian(QQ)
        True
        sage: A = SteenrodAlgebra(3)
        sage: is_noetherian(A)
        False
    """

    return x.is_noetherian()

def is_odd(x):
    """
    Returns whether or not x is odd. This is by definition the
    complement of is_even.

    EXAMPLES::

        sage: is_odd(-2)
        False
        sage: is_odd(-3)
        True
        sage: is_odd(0)
        False
        sage: is_odd(1)
        True
    """
    return not is_even(x)

## def j_invariant(x):
##     """
##     Return the j_invariant of x.

##     EXAMPLES:
##         sage: E = EllipticCurve([0, -1, 1, -10, -20])
##         sage: E
##         Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
##         sage: j_invariant(E)
##         -122023936/161051
##     """
##     return x.j_invariant()

def kernel(x):
    """
    Returns the left kernel of x.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,2)
        sage: A = M([1,2,3,4,5,6])
        sage: kernel(A)
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1 -2  1]
        sage: kernel(A.transpose())
        Vector space of degree 2 and dimension 0 over Rational Field
        Basis matrix:
        []

    Here are two corner cases:
        sage: M=MatrixSpace(QQ,0,3)
        sage: A=M([])
        sage: kernel(A)
        Vector space of degree 0 and dimension 0 over Rational Field
        Basis matrix:
        []
        sage: kernel(A.transpose()).basis()
        [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
    """
    return x.kernel()

def krull_dimension(x):
    """
    Returns the Krull dimension of x.

    EXAMPLES::

        sage: krull_dimension(QQ)
        0
        sage: krull_dimension(ZZ)
        1
        sage: krull_dimension(ZZ[sqrt(5)])
        1
        sage: U.<x,y,z> = PolynomialRing(ZZ,3); U
        Multivariate Polynomial Ring in x, y, z over Integer Ring
        sage: U.krull_dimension()
        4
    """
    return x.krull_dimension()

def lift(x):
    """
    Lift an object of a quotient ring `R/I` to `R`.

    EXAMPLES: We lift an integer modulo `3`.

    ::

        sage: Mod(2,3).lift()
        2

    We lift an element of a quotient polynomial ring.

    ::

        sage: R.<x> = QQ['x']
        sage: S.<xmod> = R.quo(x^2 + 1)
        sage: lift(xmod-7)
        x - 7
    """
    try:
        return x.lift()
    except AttributeError:
        raise ArithmeticError("no lift defined.")

def log(x,b=None):
    r"""
    Returns the log of x to the base b. The default base is e.

    INPUT:


    -  ``x`` - number

    -  ``b`` - base (default: None, which means natural
       log)


    OUTPUT: number

    .. note::

       In Magma, the order of arguments is reversed from in Sage,
       i.e., the base is given first. We use the opposite ordering, so
       the base can be viewed as an optional second argument.

    EXAMPLES::

        sage: log(e^2)
        2
        sage: log(16,2)
        4
        sage: log(3.)
        1.09861228866811
    """
    if b is None:
        if hasattr(x, 'log'):
            return x.log()
        return RDF(x)._log_base(1)
    else:
        if hasattr(x, 'log'):
            return x.log(b)
        return RDF(x).log(b)

def minimal_polynomial(x, var='x'):
    """
    Returns the minimal polynomial of x.

    EXAMPLES::

        sage: a = matrix(ZZ, 2, [1..4])
        sage: minpoly(a)
        x^2 - 5*x - 2
        sage: minpoly(a,'t')
        t^2 - 5*t - 2
        sage: minimal_polynomial(a)
        x^2 - 5*x - 2
        sage: minimal_polynomial(a,'theta')
        theta^2 - 5*theta - 2
    """
    try:
        return x.minpoly(var=var)
    except AttributeError:
        return x.minimal_polynomial(var=var)

minpoly = minimal_polynomial


def multiplicative_order(x):
    r"""
    Returns the multiplicative order of self, if self is a unit, or
    raise ``ArithmeticError`` otherwise.

    EXAMPLES::

        sage: a = mod(5,11)
        sage: multiplicative_order(a)
        5
        sage: multiplicative_order(mod(2,11))
        10
        sage: multiplicative_order(mod(2,12))
        Traceback (most recent call last):
        ...
        ArithmeticError: multiplicative order of 2 not defined since it is not a unit modulo 12
    """
    return x.multiplicative_order()

## def new_submodule(x):
##     return x.new_submodule()

## def new_subspace(x):
##     return x.new_subspace()

def ngens(x):
    """
    Returns the number of generators of x.

    EXAMPLES::

        sage: R.<x,y> = SR[]; R
        Multivariate Polynomial Ring in x, y over Symbolic Ring
        sage: ngens(R)
        2
        sage: A = AbelianGroup(5, [5,5,7,8,9])
        sage: ngens(A)
        5
        sage: ngens(ZZ)
        1
    """
    return x.ngens()

def norm(x):
    r"""
    Returns the norm of ``x``.

    For matrices and vectors, this returns the L2-norm. The L2-norm of a
    vector `\textbf{v} = (v_1, v_2, \dots, v_n)`, also called the Euclidean
    norm, is defined as

    .. MATH::

        |\textbf{v}|
        =
        \sqrt{\sum_{i=1}^n |v_i|^2}

    where `|v_i|` is the complex modulus of `v_i`. The Euclidean norm is often
    used for determining the distance between two points in two- or
    three-dimensional space.

    For complex numbers, the function returns the field norm. If
    `c = a + bi` is a complex number, then the norm of `c` is defined as the
    product of `c` and its complex conjugate

    .. MATH::

        \text{norm}(c)
        =
        \text{norm}(a + bi)
        =
        c \cdot \overline{c}
        =
        a^2 + b^2.

    The norm of a complex number is different from its absolute value.
    The absolute value of a complex number is defined to be the square
    root of its norm. A typical use of the complex norm is in the
    integral domain `\ZZ[i]` of Gaussian integers, where the norm of
    each Gaussian integer `c = a + bi` is defined as its complex norm.

    .. SEEALSO::

        - :meth:`sage.matrix.matrix2.Matrix.norm`

        - :meth:`sage.modules.free_module_element.FreeModuleElement.norm`

        - :meth:`sage.rings.complex_double.ComplexDoubleElement.norm`

        - :meth:`sage.rings.complex_number.ComplexNumber.norm`

        - :meth:`sage.symbolic.expression.Expression.norm`

    EXAMPLES:

    The norm of vectors::

        sage: z = 1 + 2*I
        sage: norm(vector([z]))
        sqrt(5)
        sage: v = vector([-1,2,3])
        sage: norm(v)
        sqrt(14)
        sage: _ = var("a b c d")
        sage: v = vector([a, b, c, d])
        sage: norm(v)
        sqrt(abs(a)^2 + abs(b)^2 + abs(c)^2 + abs(d)^2)

    The norm of matrices::

        sage: z = 1 + 2*I
        sage: norm(matrix([[z]]))
        2.2360679775
        sage: M = matrix(ZZ, [[1,2,4,3], [-1,0,3,-10]])
        sage: norm(M)
        10.6903311292
        sage: norm(CDF(z))
        5.0
        sage: norm(CC(z))
        5.00000000000000

    The norm of complex numbers::

        sage: z = 2 - 3*I
        sage: norm(z)
        13
        sage: a = randint(-10^10, 100^10)
        sage: b = randint(-10^10, 100^10)
        sage: z = a + b*I
        sage: bool(norm(z) == a^2 + b^2)
        True

    The complex norm of symbolic expressions::

        sage: a, b, c = var("a, b, c")
        sage: assume((a, 'real'), (b, 'real'), (c, 'real'))
        sage: z = a + b*I
        sage: bool(norm(z).simplify() == a^2 + b^2)
        True
        sage: norm(a + b).simplify()
        a^2 + 2*a*b + b^2
        sage: v = vector([a, b, c])
        sage: bool(norm(v).simplify() == sqrt(a^2 + b^2 + c^2))
        True
        sage: forget()
    """
    return x.norm()

def numerator(x):
    """
    Returns the numerator of x.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: F = FractionField(R)
        sage: r = (x+1)/(x-1)
        sage: numerator(r)
        x + 1
        sage: numerator(17/11111)
        17
    """
    if isinstance(x, (int, long)):
        return x
    return x.numerator()

# Following is the top-level numerical_approx function.
# Implement a ._numerical_approx(prec, digits, algorithm) method for your
# objects to enable the three top-level functions and three methods

def numerical_approx(x, prec=None, digits=None, algorithm=None):
    r"""
    Returns a numerical approximation of an object ``x`` with at
    least ``prec`` bits (or decimal ``digits``) of precision.

    .. note::

       Both upper case ``N`` and lower case ``n`` are aliases for
       :func:`numerical_approx`, and all three may be used as
       methods.

    INPUT:

    -  ``x`` - an object that has a numerical_approx
       method, or can be coerced into a real or complex field
    -  ``prec`` (optional) - an integer (bits of
       precision)
    -  ``digits`` (optional) - an integer (digits of
       precision)
    -  ``algorithm`` (optional) - a string specifying
       the algorithm to use for functions that implement
       more than one

    If neither the ``prec`` or ``digits`` are specified,
    the default is 53 bits of precision.  If both are
    specified, then ``prec`` is used.

    EXAMPLES::

        sage: numerical_approx(pi, 10)
        3.1
        sage: numerical_approx(pi, digits=10)
        3.141592654
        sage: numerical_approx(pi^2 + e, digits=20)
        12.587886229548403854
        sage: n(pi^2 + e)
        12.5878862295484
        sage: N(pi^2 + e)
        12.5878862295484
        sage: n(pi^2 + e, digits=50)
        12.587886229548403854194778471228813633070946500941
        sage: a = CC(-5).n(prec=100)
        sage: b = ComplexField(100)(-5)
        sage: a == b
        True
        sage: type(a) == type(b)
        True
        sage: numerical_approx(9)
        9.00000000000000

    You can also usually use method notation.  ::

        sage: (pi^2 + e).n()
        12.5878862295484
        sage: (pi^2 + e).N()
        12.5878862295484
        sage: (pi^2 + e).numerical_approx()
        12.5878862295484

    Vectors and matrices may also have their entries approximated.  ::

        sage: v = vector(RDF, [1,2,3])
        sage: v.n()
        (1.00000000000000, 2.00000000000000, 3.00000000000000)

        sage: v = vector(CDF, [1,2,3])
        sage: v.n()
        (1.00000000000000, 2.00000000000000, 3.00000000000000)
        sage: _.parent()
        Vector space of dimension 3 over Complex Field with 53 bits of precision
        sage: v.n(prec=75)
        (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)

        sage: u = vector(QQ, [1/2, 1/3, 1/4])
        sage: n(u, prec=15)
        (0.5000, 0.3333, 0.2500)
        sage: n(u, digits=5)
        (0.50000, 0.33333, 0.25000)

        sage: v = vector(QQ, [1/2, 0, 0, 1/3, 0, 0, 0, 1/4], sparse=True)
        sage: u = v.numerical_approx(digits=4)
        sage: u.is_sparse()
        True
        sage: u
        (0.5000, 0.0000, 0.0000, 0.3333, 0.0000, 0.0000, 0.0000, 0.2500)

        sage: A = matrix(QQ, 2, 3, range(6))
        sage: A.n()
        [0.000000000000000  1.00000000000000  2.00000000000000]
        [ 3.00000000000000  4.00000000000000  5.00000000000000]

        sage: B = matrix(Integers(12), 3, 8, srange(24))
        sage: N(B, digits=2)
        [0.00  1.0  2.0  3.0  4.0  5.0  6.0  7.0]
        [ 8.0  9.0  10.  11. 0.00  1.0  2.0  3.0]
        [ 4.0  5.0  6.0  7.0  8.0  9.0  10.  11.]

    Internally, numerical approximations of real numbers are stored in base-2.
    Therefore, numbers which look the same in their decimal expansion might be
    different::

        sage: x=N(pi, digits=3); x
        3.14
        sage: y=N(3.14, digits=3); y
        3.14
        sage: x==y
        False
        sage: x.str(base=2)
        '11.001001000100'
        sage: y.str(base=2)
        '11.001000111101'

    As an exceptional case, ``digits=1`` usually leads to 2 digits (one
    significant) in the decimal output (see :trac:`11647`)::

        sage: N(pi, digits=1)
        3.2
        sage: N(pi, digits=2)
        3.1
        sage: N(100*pi, digits=1)
        320.
        sage: N(100*pi, digits=2)
        310.


    In the following example, ``pi`` and ``3`` are both approximated to two
    bits of precision and then subtracted, which kills two bits of precision::

        sage: N(pi, prec=2)
        3.0
        sage: N(3, prec=2)
        3.0
        sage: N(pi - 3, prec=2)
        0.00

    TESTS::

        sage: numerical_approx(I)
        1.00000000000000*I
        sage: x = QQ['x'].gen()
        sage: F.<k> = NumberField(x^2+2, embedding=sqrt(CC(2))*CC.0)
        sage: numerical_approx(k)
        1.41421356237309*I

        sage: type(numerical_approx(CC(1/2)))
        <type 'sage.rings.complex_number.ComplexNumber'>

    The following tests :trac:`10761`, in which ``n()`` would break when
    called on complex-valued algebraic numbers.  ::

        sage: E = matrix(3, [3,1,6,5,2,9,7,3,13]).eigenvalues(); E
        [18.16815365088822?, -0.08407682544410650? - 0.2190261484802906?*I, -0.08407682544410650? + 0.2190261484802906?*I]
        sage: E[1].parent()
        Algebraic Field
        sage: [a.n() for a in E]
        [18.1681536508882, -0.0840768254441065 - 0.219026148480291*I, -0.0840768254441065 + 0.219026148480291*I]

    Make sure we've rounded up log(10,2) enough to guarantee
    sufficient precision (trac #10164)::

        sage: ks = 4*10**5, 10**6
        sage: check_str_length = lambda k: len(str(numerical_approx(1+10**-k,digits=k+1)))-1 >= k+1
        sage: check_precision = lambda k: numerical_approx(1+10**-k,digits=k+1)-1 > 0
        sage: all(check_str_length(k) and check_precision(k) for k in ks)
        True

    Testing we have sufficient precision for the golden ratio (:trac:`12163`), note
    that the decimal point adds 1 to the string length::

        sage: len(str(n(golden_ratio, digits=5000)))
        5001
        sage: len(str(n(golden_ratio, digits=5000000)))  # long time (4s on sage.math, 2012)
        5000001

    Check that :trac:`14778` is fixed::

        sage: n(0, algorithm='foo')
        0.000000000000000
    """
    if prec is None:
        if digits is None:
            prec = 53
        else:
            prec = int((digits+1) * LOG_TEN_TWO_PLUS_EPSILON) + 1
    try:
        return x._numerical_approx(prec, algorithm=algorithm)
    except AttributeError:
        from sage.rings.complex_double import is_ComplexDoubleElement
        from sage.rings.complex_number import is_ComplexNumber
        if not (is_ComplexNumber(x) or is_ComplexDoubleElement(x)):
            try:
                return sage.rings.real_mpfr.RealField(prec)(x)
            # Trac 10761: now catches ValueErrors as well as TypeErrors
            except (TypeError, ValueError):
                pass
        return sage.rings.complex_field.ComplexField(prec)(x)

n = numerical_approx
N = numerical_approx

def objgens(x):
    """
    EXAMPLES::

        sage: R, x = objgens(PolynomialRing(QQ,3, 'x'))
        sage: R
        Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        sage: x
        (x0, x1, x2)
    """
    return x.objgens()

def objgen(x):
    """
    EXAMPLES::

        sage: R, x = objgen(FractionField(QQ['x']))
        sage: R
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: x
        x
    """
    return x.objgen()

def one(R):
    """
    Returns the one element of the ring R.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: one(R)*x == x
        True
        sage: one(R) in R
        True
    """
    return R(1)

def order(x):
    """
    Returns the order of x. If x is a ring or module element, this is
    the additive order of x.

    EXAMPLES::

        sage: C = CyclicPermutationGroup(10)
        sage: order(C)
        10
        sage: F = GF(7)
        sage: order(F)
        7
    """
    return x.order()

def rank(x):
    """
    Returns the rank of x.

    EXAMPLES: We compute the rank of a matrix::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: rank(A)
        2

    We compute the rank of an elliptic curve::

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: rank(E)
        1
    """
    return x.rank()

def regulator(x):
    """
    Returns the regulator of x.

    EXAMPLES::

        sage: regulator(NumberField(x^2-2, 'a'))
        0.881373587019543
        sage: regulator(EllipticCurve('11a'))
        1.00000000000000
    """
    return x.regulator()

def round(x, ndigits=0):
    """
    round(number[, ndigits]) - double-precision real number

    Round a number to a given precision in decimal digits (default 0
    digits). If no precision is specified this just calls the element's
    .round() method.

    EXAMPLES::

        sage: round(sqrt(2),2)
        1.41
        sage: q = round(sqrt(2),5); q
        1.41421
        sage: type(q)
        <type 'sage.rings.real_double.RealDoubleElement'>
        sage: q = round(sqrt(2)); q
        1
        sage: type(q)
        <type 'sage.rings.integer.Integer'>
        sage: round(pi)
        3
        sage: b = 5.4999999999999999
        sage: round(b)
        5


    Since we use floating-point with a limited range, some roundings can't
    be performed::

        sage: round(sqrt(Integer('1'*1000)),2)
        +infinity

    IMPLEMENTATION: If ndigits is specified, it calls Python's builtin
    round function, and converts the result to a real double field
    element. Otherwise, it tries the argument's .round() method; if
    that fails, it reverts to the builtin round function, converted to
    a real double field element.

    .. note::

       This is currently slower than the builtin round function, since
       it does more work - i.e., allocating an RDF element and
       initializing it. To access the builtin version do
       ``import __builtin__; __builtin__.round``.
    """
    try:
        if ndigits:
            return RealDoubleElement(__builtin__.round(x, ndigits))
        else:
            try:
                return x.round()
            except AttributeError:
                return RealDoubleElement(__builtin__.round(x, 0))
    except ArithmeticError:
        if not isinstance(x, RealDoubleElement):
            return round(RDF(x), ndigits)
        else:
            raise

def quotient(x, y, *args, **kwds):
    """
    Returns the quotient object x/y, e.g., a quotient of numbers or of a
    polynomial ring x by the ideal generated by y, etc.

    EXAMPLES::

        sage: quotient(5,6)
        5/6
        sage: quotient(5.,6.)
        0.833333333333333
        sage: R.<x> = ZZ[]; R
        Univariate Polynomial Ring in x over Integer Ring
        sage: I = Ideal(R, x^2+1)
        sage: quotient(R, I)
        Univariate Quotient Polynomial Ring in xbar over Integer Ring with modulus x^2 + 1
    """
    try:
        return x.quotient(y, *args, **kwds)
    except AttributeError:
        return x/y

quo = quotient

def show(x, *args, **kwds):
    r"""
    Show a graphics object x.

    For additional ways to show objects in the notebook, look
    at the methods on the html object.  For example,
    html.table will produce an HTML table from a nested
    list.


    OPTIONAL INPUT:


    -  ``filename`` - (default: None) string


    SOME OF THESE MAY APPLY:

    - ``dpi`` - dots per inch

    - ``figsize``- [width, height] (same for square aspect)

    - ``axes`` - (default: True)

    - ``fontsize`` - positive integer

    - ``frame`` - (default: False) draw a MATLAB-like frame around the
      image

    EXAMPLES::

        sage: show(graphs(3))
        sage: show(list(graphs(3)))
    """
    if not isinstance(x, (sage.interfaces.expect.Expect, sage.interfaces.expect.ExpectElement)):
        try:
            return x.show(*args, **kwds)
        except AttributeError:
            pass
    if isinstance(x, sage.interfaces.mathematica.MathematicaElement):
        return x.show(*args, **kwds)

    import types
    if isinstance(x, types.GeneratorType):
        x = list(x)
    if isinstance(x, list):
        if len(x) > 0:
            from sage.graphs.graph import GenericGraph
            if isinstance(x[0], GenericGraph):
                import sage.graphs.graph_list as graphs_list
                graphs_list.show_graphs(x)
                return
    _do_show(x)

def _do_show(x):
    if sage.doctest.DOCTEST_MODE:
        return sage.misc.latex.latex(x)
    from latex import view
    view(x, mode='display')
    #raise AttributeError, "object %s does not support show."%(x, )

def sqrt(x):
    """
    Returns a square root of x.

    This function (``numerical_sqrt``) is deprecated.  Use ``sqrt(x,
    prec=n)`` instead.

    EXAMPLES::

        sage: numerical_sqrt(10.1)
        doctest:1: DeprecationWarning: numerical_sqrt is deprecated, use sqrt(x, prec=n) instead
        See http://trac.sagemath.org/5404 for details.
        3.17804971641414
        sage: numerical_sqrt(9)
        3
    """
    from sage.misc.superseded import deprecation
    deprecation(5404, "numerical_sqrt is deprecated, use sqrt(x, prec=n) instead")
    try: return x.sqrt()
    except (AttributeError, ValueError):
        try:
            return RDF(x).sqrt()
        except TypeError:
            return CDF(x).sqrt()

def isqrt(x):
    """
    Returns an integer square root, i.e., the floor of a square root.

    EXAMPLES::

        sage: isqrt(10)
        3
        sage: isqrt(10r)
        3
    """
    try:
        return x.isqrt()
    except AttributeError:
        from sage.functions.all import floor
        n = sage.rings.integer.Integer(floor(x))
        return n.isqrt()

def squarefree_part(x):
    """
    Returns the square free part of `x`, i.e., a divisor
    `z` such that `x = z y^2`, for a perfect square
    `y^2`.

    EXAMPLES::

        sage: squarefree_part(100)
        1
        sage: squarefree_part(12)
        3
        sage: squarefree_part(10)
        10
        sage: squarefree_part(216r) # see #8976
        6

    ::

        sage: x = QQ['x'].0
        sage: S = squarefree_part(-9*x*(x-6)^7*(x-3)^2); S
        -9*x^2 + 54*x
        sage: S.factor()
        (-9) * (x - 6) * x

    ::

        sage: f = (x^3 + x + 1)^3*(x-1); f
        x^10 - x^9 + 3*x^8 + 3*x^5 - 2*x^4 - x^3 - 2*x - 1
        sage: g = squarefree_part(f); g
        x^4 - x^3 + x^2 - 1
        sage: g.factor()
        (x - 1) * (x^3 + x + 1)
    """
    try:
        return x.squarefree_part()
    except AttributeError:
        pass
    F = factor(x)
    n = parent(x)(1)
    for p, e in F:
        if e%2 != 0:
            n *= p
    return n * F.unit()

## def square_root(x):
##     """
##     Return a square root of x with the same parent as x, if possible,
##     otherwise raise a ValueError.
##     EXAMPLES:
##         sage: square_root(9)
##         3
##         sage: square_root(100)
##         10
##     """
##     try:
##         return x.square_root()
##     except AttributeError:
##         raise NotImplementedError

def transpose(x):
    """
    Returns the transpose of x.

    EXAMPLES::

        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: transpose(A)
        [1 4 7]
        [2 5 8]
        [3 6 9]
    """
    return x.transpose()

## def vector(x, R):
##     r"""
##     Return the \sage vector over $R$ obtained from x, if possible.
##     """
##     try:
##         return x._vector_(R)
##     except AttributeError:
##         import sage.modules.free_module_element
##         return sage.modules.free_module_element.Vector(x, R)

def zero(R):
    """
    Returns the zero element of the ring R.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: zero(R) in R
        True
        sage: zero(R)*x == zero(R)
        True
    """
    return R(0)



#################################################################
# Generic parent
#################################################################
def parent(x):
    """
    Returns x.parent() if defined, or type(x) if not.

    EXAMPLE::

        sage: Z = parent(int(5))
        sage: Z(17)
        17
        sage: Z
        <type 'int'>
    """
    try:
        return x.parent()
    except AttributeError:
        return type(x)

