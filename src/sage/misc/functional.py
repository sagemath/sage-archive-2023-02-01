"""
Functional notation

These are function so that you can write foo(x) instead of x.foo() in
certain common cases.

AUTHORS: Initial version -- William Stein
         More Examples -- David Joyner, 2005-12-20
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

import math

import sage.rings.all
import sage.categories.all
import sage.misc.latex
import sage.server.support
import sage.interfaces.expect

##############################################################################
# There are many functions on elements of a ring, which mathematicians
# usually write f(x), e.g., it is weird to write x.log() and natural
# to write log(x).  The functions below allow for the more familiar syntax.
##############################################################################
def additive_order(x):
    """
    Return the additive order of $x$.
    """
    return x.additive_order()

def arg(x):
    """
    Return the argument of a complex number $x$.

    EXAMPLES:
        sage: z = 1+2*I
        sage: theta = arg(z)
        sage: cos(theta)*abs(z)
        1.0000000000000002
        sage: sin(theta)*abs(z)
        2.0000000000000000
    """
    R = sage.rings.all.PolynomialRing(sage.rings.all.QQ,'z',inject_variables=True)
    print z
    try: return x.arg()
    except AttributeError: return sage.rings.all.CC(x).arg()

def base_ring(x):
    """
    Return the base ring over which x is defined.

    EXAMPLES:
        sage: R = PolynomialRing(GF(7))
        sage: base_ring(R)
        Finite Field of size 7
    """
    print z
    return x.base_ring()

def base_field(x):
    """
    Return the base field over which x is defined.
    """
    return x.base_field()

def basis(x):
    """
    Return the fixed basis of x.

    EXAMPLES:
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
    Return the category of x.

    EXAMPLES:
        sage: V = VectorSpace(QQ,3)
        sage: category(V)
        Category of vector spaces over Rational Field
    """
    try:
        return x.category()
    except AttributeError:
        return sage.categories.all.Objects()

def ceil(x):
    try:
        return x.ceil()
    except AttributeError:
        return sage.rings.all.ceil(x)

def charpoly(x, var):
    """
    Return the characteristic polynomial of x in the given variable.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: charpoly(A)
        x^3 - 15*x^2 - 18*x

        sage: GF(7^10, 'alpha')
        Finite Field in alpha of size 7^10
        sage: alpha.charpoly('T')
        T^10 + T^6 + T^5 + 4*T^4 + T^3 + 2*T^2 + 3*T + 3
    """
    try:
        return x.charpoly(var)
    except AttributeError:
        raise NotImplementedError, "computation of charpoly of x (=%s) not implemented"%x

def acos(x):
    """
    Return the arc cosine of x.

    EXAMPLES:
        sage: acos(0.5)
        1.0471975511965979
        sage: acos(1 + I*1.0)
        0.90455689430238140 - 1.0612750619050357*I
    """
    try: return x.acos()
    except AttributeError: return sage.rings.all.RR(x).acos()

def asin(x):
    """
    Return the arc sine of x.

    EXAMPLES:
        sage: asin(0.5)
        0.52359877559829893
        sage: asin(1 + I*1.0)
        0.66623943249251527 + 1.0612750619050357*I
    """
    try: return x.asin()
    except AttributeError: return sage.rings.all.RR(x).asin()

def atan(x):
    """
    Return the arc tangent of x.

    EXAMPLES:
        sage: atan(1/2)
        0.46364760900080609
        sage: atan(1 + I)
        1.0172219678978514 + 0.40235947810852507*I
    """
    try: return x.atan()
    except AttributeError: return sage.rings.all.RR(x).atan()

## def cuspidal_submodule(x):
##     return x.cuspidal_submodule()

## def cuspidal_subspace(x):
##     return x.cuspidal_subspace()

def cyclotomic_polynomial(n):
    """
    EXAMPLES:
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
    return sage.rings.all.PolynomialRing(\
                  sage.rings.all.QQ).cyclotomic_polynomial(n)

def decomposition(x):
    """
    Return the decomposition of x.
    """
    return x.decomposition()

def denominator(x):
    """
    Return the denominator of x.

    EXAMPLES:
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

def derivative(x):
    """
    Return the derivative of a polynomial x.

    EXAMPLES:
        sage: f = cyclotomic_polynomial(10)
        sage: derivative(f)
        4*x^3 - 3*x^2 + 2*x - 1
        sage: R.<x> = PolynomialRing(GF(7))
        sage: f = x^7 + x
        sage: derivative(f)
        1
    """
    return x.derivative()

def det(x):
    """
    Return the determinant of x.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: det(A)
        0
    """
    return x.det()

def dimension(x):
    """
    Return the dimension of x.

    EXAMPLES:
        sage: V = VectorSpace(QQ,3)
        sage: S = V.subspace([[1,2,0],[2,2,-1]])
        sage: dimension(S)
        2
    """
    return x.dimension()

dim = dimension

def discriminant(x):
    """
    EXAMPLES:
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
    Return the value of the eta function at $x$, which must
    be in the upper half plane.

    The $\eta$ function is
        $$
           \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty}(1-e^{2\pi inz})
        $$

    EXAMPLES:
        sage: eta(1+I)
        0.74204877583656470 + 0.19883137022991071*I
    """
    try: return x.eta()
    except AttributeError: return sage.rings.all.CC(x).eta()

def exp(x):
    """
    Return the value of the exponentation function at x.
    """
    try: return x.exp()
    except AttributeError: return sage.rings.all.RR(x).exp()

def factor(x, *args, **kwds):
    """
    Return the prime factorization of x.

    EXAMPLES:
        sage: factor(factorial(10))
        2^8 * 3^4 * 5^2 * 7
        sage: n = next_prime(10^6); n
        1000003
        sage: factor(n)
        1000003
    """
    try: return x.factor(*args, **kwds)
    except AttributeError: return sage.rings.all.factor(x, *args, **kwds)

factorization = factor
factorisation = factor

def fcp(x, var):
    """
    Return the factorization of the characteristic polynomial
    of x.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: fcp(A, 'x')
        x * (x^2 - 15*x - 18)
    """
    try: return x.fcp(var)
    except AttributeError: return factor(charpoly(x, var))

gcd = sage.rings.arith.gcd


def floor(x):
    try:
        return x.floor()
    except AttributeError:
        return sage.rings.all.floor(x)

def gen(x):
    """
    Return the generator of x.
    """
    return x.gen()

def gens(x):
    """
    Return the generators of x.
    """
    return x.gens()

def hecke_operator(x,n):
    """
    Return the n-th Hecke operator T_n acting on x.

    EXAMPLES:
        sage: M = ModularSymbols(1,12)
        sage: hecke_operator(M,5)
        Hecke operator T_5 on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
    """
    return x.hecke_operator(n)

def ideal(*x):
    """
    Return the ideal generated by x where x is an element or list.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(QQ)
        sage: ideal(x^2-2*x+1, x^2-1)
        Principal ideal (x - 1) of Univariate Polynomial Ring in x over Rational Field
        sage: ideal([x^2-2*x+1, x^2-1])
        Principal ideal (x - 1) of Univariate Polynomial Ring in x over Rational Field
    """
    if isinstance(x[0], (list, tuple)):
        return sage.rings.all.Ideal(x[0])
    if len(x) == 1:
        try:
            return x[0].ideal()
        except AttributeError:
            pass
    return sage.rings.all.Ideal(x)

def image(x):
    """
    Return the image of x.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: image(A)
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1]
        [ 0  1  2]
    """
    return x.image()

def imag(x):
    """
    Return the imaginary part of x.
    """
    try: return x.imag()
    except AttributeError: return sage.rings.all.CC(x).imag()

def imaginary(x):
    """
    Return the imaginary part of a complex number.

    EXAMPLES:
        sage: z = 1+2*I
        sage: imaginary(z)
        2.0000000000000000
        sage: imag(z)
        2.0000000000000000
    """
    return imag(x)

def integral(x, var=None, algorithm='maxima'):
    """
    Return an indefinite integral of an object x.

    First call x.integrate() and if that fails make an
    object and integrate it using maxima, maple, etc, as
    specified by algorithm.

    EXAMPLES:
        sage: f = cyclotomic_polynomial(10)
        sage: integral(f)
        1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x
    """
    if var is None:
        try:
            return x.integral()
        except AttributeError:
            pass
    import sage.interfaces.all as I
    if var is None:
        var = 'x'
    if algorithm == 'maxima':
        return I.maxima(x).integrate(var)
    elif algorithm == 'mathematica':
        return I.mathematica(x).Integrate(var)
    else:
        raise ValueError, 'no algorithm %s'%algorithm

def integral_closure(x):
    return x.integral_closure()

def interval(a, b):
    r"""
    Integers between a and b \emph{inclusive} (a and b integers).

    EXAMPLES:
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
    Iterator over the integers between a and b, \emph{inclusive}.
    """
    return xrange(a, b+1)

def is_commutative(x):
    """
    EXAMPLES:
        sage: R = PolynomialRing(QQ)
        sage: is_commutative(R)
        True
    """
    return x.is_commutative()

def is_even(x):
    """
    Return whether or not an integer x is even, e.g., divisible by 2.

    EXAMPLES:
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
    return x.is_integrally_closed()

def is_field(x):
    """
    EXAMPLES:
        sage: R = PolynomialRing(QQ)
        sage: F = FractionField(R)
        sage: is_field(F)
        True
    """
    return x.is_field()

def is_noetherian(x):
    return x.is_noetherian()

def is_odd(x):
    """
    Return whether or not x is odd.  This is by definition the
    complement of is_even.

    EXAMPLES:
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
    Return the kernel of x.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: kernel(A)
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1 -2  1]
    """
    return x.kernel()

def krull_dimension(x):
    return x.krull_dimension()

def lift(x):
    """
    Lift an object of a quotient ring $R/I$ to $R$.

    EXAMPLES:
    We lift an integer modulo $3$.
        sage: Mod(2,3).lift()
        2

    We lift an element of a quotient polynomial ring.
        sage: R.<x> = QQ['x']
        sage: S.<xmod> = R/(x^2 + 1)
        sage: lift(xmod-7)
        x - 7
    """
    try:
        return x.lift()
    except AttributeError:
        raise ArithmeticError, "no lift defined."

lcm = sage.rings.arith.lcm

def log(x,b=None):
    r"""
    Return the log of x to the base b.  The default base is e.

    INPUT:
        x -- number
        b -- base (default: None, which means natural log)
    OUTPUT:
        number

    \note{In Magma, the order of arguments is reversed from in
    \sage, i.e., the base is given first.  We use the opposite
    ordering, so the base can be viewed as an optional second
    argument.}

    EXAMPLES:
        sage: log(10,2)
        3.3219280948873626
        sage: log(8,2)
        3.0000000000000000
        sage: log(10)
        2.3025850929940459
        sage: log(2.718)
        0.99989631572895199
    """
    if b is None:
        try: return x.log()
        except AttributeError:
            return sage.rings.all.RR(x).log()
    else:
        try: return x.log(b)
        except AttributeError:
            return log(x) / log(b)

def minimal_polynomial(x):
    """
    Return the minimal polynomial of x.
    """
    return x.minimal_polynomial()


def multiplicative_order(x):
    r"""
    Return the multiplicative order of self, if self is a unit, or raise
    \code{ArithmeticError} otherwise.
    """
    return x.multiplicative_order()

## def new_submodule(x):
##     return x.new_submodule()

## def new_subspace(x):
##     return x.new_subspace()

def ngens(x):
    """
    Return the number of generators of x.
    """
    return x.ngens()

def norm(x):
    """
    Return the norm of x.

    EXAMPLES:
        sage: z = 1+2*I
        sage: norm(z)
        5.0000000000000000
    """
    return x.norm()

def numerator(x):
    """
    Return the numerator of x.

    EXAMPLES:
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

def objgens(x, names=None):
    """
    EXAMPLES:
        sage: R, x = objgens(MPolynomialRing(QQ,3))
        sage: R
        Polynomial Ring in x0, x1, x2 over Rational Field
        sage: x
        (x0, x1, x2)
    """
    return x.objgens(names)

def objgen(x, names=None):
    """
    EXAMPLES:
        sage: R, x = objgen(FractionField(QQ['x']))
        sage: R
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: x
        x
    """
    return x.objgen(names)

def one(R):
    """
    Return the one element of the ring R.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(QQ)
        sage: one(R)*x == x
        True
        sage: one(R) in R
        True
    """
    return R(1)

def order(x):
    """
    Return the order of x.  If x is a ring or module element, this is
    the additive order of x.

    EXAMPLES:
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
    Return the rank of x.

    EXAMPLES:
    We compute the rank of a matrix:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: rank(A)
        2

    We compute the rank of an elliptic curve:
        sage: E=EllipticCurve([0,0,1,-1,0])
        sage: rank(E)
        1
    """
    return x.rank()

def real(x):
    """
    Return the real part of x.

    EXAMPLES:
        sage: z = 1+2*I
        sage: real(z)
        1.0000000000000000
    """
    try: return x.real()
    except AttributeError: return sage.rings.all.CC(x).real()

def regulator(x):
    """
    Return the regulator of x.
    """
    return x.regulator()

def quotient(x, y, *args, **kwds):
    """
    Return the quotient object x/y, e.g., a quotient of numbers or of
    a polynomial ring x by the ideal generated by y, etc.
    """
    try:
        return x.quotient(y, *args, **kwds)
    except AttributeError:
        return x/y

quo = quotient

def show(x, *args, **kwds):
    """
    Show a graphics object x.
    """
    if not isinstance(x, (sage.interfaces.expect.Expect, sage.interfaces.expect.ExpectElement)):
        try:
            return x.show(*args, **kwds)
        except AttributeError:
            pass
    if sage.server.support.EMBEDDED_MODE:
        print '<html><div class="math">%s</div></html>'%sage.misc.latex.latex(x)
        return sage.misc.latex.LatexExpr('') # so not visible output
    raise AttributeError, "object %s does not support show."%x

def sqrt(x):
    """
    Return a square root of x.

    EXAMPLES:
        sage: sqrt(10.1)
        3.1780497164141406
        sage: sqrt(9)
        3
    """
    try: return x.sqrt()
    except (AttributeError, ValueError): return sage.rings.all.CC(x).sqrt()

def isqrt(x):
    """
    Return an integer square root, i.e., the floor of a
    square root.

    EXAMPLES:
        sage: isqrt(10)
        3
    """
    try: return x.isqrt()
    except AttributeError:
        raise NotImplementedError

def square_free_part(x):
    """
    Return the square free part of $x$, i.e., a divisor $z$ such that $x = z y^2$,
    for a perfect square $y^2$.

    EXAMPLES:
        sage: square_free_part(100)
        1
        sage: square_free_part(12)
        3
        sage: square_free_part(10)
        10

        sage: x = QQ['x'].0
        sage: S = square_free_part(-9*x*(x-6)^7*(x-3)^2); S
        -9*x^2 + 54*x
        sage: S.factor()
        (-9) * (x - 6) * x

        sage: f = (x^3 + x + 1)^3*(x-1); f
        x^10 - x^9 + 3*x^8 + 3*x^5 - 2*x^4 - x^3 - 2*x - 1
        sage: g = square_free_part(f); g
        x^4 - x^3 + x^2 - 1
        sage: g.factor()
        (x - 1) * (x^3 + x + 1)
    """
    try:
        return x.square_free_part()
    except AttributeError:
        pass
    F = factor(x)
    n = x.parent()(1)
    for p, e in F:
        if e%2 != 0:
            n *= p
    return n * F.unit()

squarefree_part = square_free_part

def square_root(x):
    """
    Return a square root of x with the same parent as x, if possible,
    otherwise raise a ValueError.

    EXAMPLES:
        sage: square_root(9)
        3
        sage: square_root(100)
        10
    """
    try:
        return x.square_root()
    except AttributeError:
        raise NotImplementedError

def tan(x):
    """
    Return the tangent of x.

    EXAMPLES:
        sage: tan(3.1415)
        -0.000092653590058191322
        sage: tan(3.1415/4)
        0.99995367427815629
    """
    try: return x.tan()
    except AttributeError: return sage.rings.all.RR(x).tan()

def transpose(x):
    """
    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,2,3,4,5,6,7,8,9])
        sage: transpose(A)
        [1 4 7]
        [2 5 8]
        [3 6 9]
    """
    return x.transpose()

xgcd = sage.rings.arith.xgcd

def vector(x, R):
    """
    Return the \sage vector over $R$ obtained from x, if possible.
    """
    try:
        return x._vector_(R)
    except AttributeError:
        raise TypeError, "No known way to create a vector from %s"%x

def zero(R):
    """
    Return the zero element of the ring R.

    EXAMPLES:
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
    Return x.parent() if defined, or type(x) if not.

    EXAMPLE:
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

