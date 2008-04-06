r"""
Dimensions of spaces of modular forms

AUTHORS:
    -- William Stein
    -- Jordi Quer

ACKNOWLEDGEMENT:
    The dimension formulas and implementations in this module grew out
    of a program that Bruce Kaskel wrote (around 1996) in PARI, which
    Kevin Buzzard subsequently extended.  I (William Stein) then
    implemented it in C++ for Hecke.  I also implemented it in Magma.
    Also, the functions for dimensions of spaces with nontrivial
    character are based on a paper (that has no proofs) by Cohen and
    Oesterle (Springer Lecture notes in math, volume 627, pages
    69--78).  The formulas for $\Gamma_H(N)$ were found and implemented by
    Jordi Quer.

The formulas here are more complete than in Hecke or Magma.

Currently the input to each function below is an integer and either a
Dirichlet character $\eps$ or a congruence subgroup, which must be
either $\Gamma_0(N)$ or $\Gamma_1(N)$.  If the input is a Dirichlet
character $\eps$, the dimensions are for subspaces of
$M_k(\Gamma_1(N), \eps)$, where $N$ is the modulus of $\eps$.
"""

##########################################################################
#       Copyright (C) 2004,2005,2006,2007,2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##########################################################################

import math

from sage.rings.arith import (factor, euler_phi as phi, divisors, is_prime,
                              valuation, kronecker_symbol, gcd, euler_phi, lcm)
import sage.modular.congroup as congroup
from sage.misc.misc import mul
from sage.rings.all import Mod, Integer, IntegerRing, IntegerModRing, ZZ, moebius
from sage.rings.rational_field import frac
import dirichlet
Z = ZZ  # useful abbreviation.

def mu0(n):
    r"""
    Return value of the combinatorial function
    $$
       \mu_0(n) = \prod_{p^r|n} (p+1)p^{r-1},
    $$
    where the product is over the maximal prime power
    divisors of $n$.

    INPUT:
        n -- an integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu0(n) for n in [1..19]]
        [1, 3, 4, 6, 6, 12, 8, 12, 12, 18, 12, 24, 14, 24, 24, 24, 18, 36, 20]
    """
    return mul([(p+1)*(p**(r-1)) for p, r in factor(n)])

def mu20(n):
    r"""
    Return value of the arithmetic function $\mu_{2,0}(n)$,
    where $\mu_{2,0}(n) = 0$ if $n$ is divisible by 4, and
    $$
     \mu_{2,0}(n) = \prod_{p|n} \left(1 + \left(\frac{-4}{p}\right)\right)
    $$
    otherwise.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu20(n) for n in [1..19]]
        [1, 1, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 2, 0, 0]
    """
    if n%4 == 0:
        return Integer(0)
    return mul([1 + kronecker_symbol(-4,p) for p, _ in factor(n)])

def mu30(n):
    r"""
    Return value of the arithmetic function $\mu_{3,0}(n)$,
    where $\mu_{3,0}(n) = 0$ if $n$ is divisible by 2 or 9,
    and
    $$
     \mu_{3,0}(n) = \prod_{p|n} \left(1 + \left(\frac{-3}{p}\right)\right)
    $$
    otherwise.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu30(n) for n in [1..19]]
        [1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2]
    """
    if n%2==0 or n%9==0:
        return Integer(0)
    return mul([1 + kronecker_symbol(-3,p) for p, _ in factor(n)])

def c0(n):
    """
    Return the number of cusps of the modular curve $X_0(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.c0(n) for n in [1..19]]
        [1, 2, 2, 3, 2, 4, 2, 4, 4, 4, 2, 6, 2, 4, 4, 6, 2, 8, 2]
        sage: [sage.modular.dims.c0(n) for n in prime_range(2,100)]
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    """
    return sum([phi(gcd(d,n//d)) for d in divisors(n)])

def g0(n):
    """
    Return the genus of the modular curve $X_0(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.g0(n) for n in [1..23]]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 2]
        sage: [n for n in [1..200] if sage.modular.dims.g0(n) == 1]
        [11, 14, 15, 17, 19, 20, 21, 24, 27, 32, 36, 49]
    """
    return Integer(1 + frac(mu0(n),12) - frac(mu20(n),4) - \
               frac(mu30(n),3) - frac(c0(n),2))

def mu1(n):
    r"""
    Return value of the combinatorial function $\mu_1(n)$,
    where for $n = 1,2$, $\mu_1(n) = \mu_0(n)$, and for
    $n > 2$,
    $$
     \mu_{1}(n) = \frac{\phi(n) \cdot \mu_0(n)}{2}.
    $$

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu1(n) for n in [1..16]]
        [1, 3, 4, 6, 12, 12, 24, 24, 36, 36, 60, 48, 84, 72, 96, 96]
    """
    if n <= 2:
        return mu0(n)
    return (phi(n)*mu0(n))/2

def mu21(n):
    r"""
    Return value of $\mu_{2,1}(n)$. By definition, this equals
    $\mu_{2,0}(n)$ for $n<4$ and equals 0 otherwise.

    EXAMPLES:
        sage: [sage.modular.dims.mu21(n) for n in [1..16]]
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    """
    if n<4:
        return mu20(n)
    return Integer(0)

def mu31(n):
    r"""
    Return value of $\mu_{3,1}(n)$. By definition, this equals
    $\mu_{3,0}(n)$ for $n<4$ and equals 0 otherwise.

    EXAMPLES:
        sage: [sage.modular.dims.mu31(n) for n in [1..10]]
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    """
    if n<4:
        return mu30(n)
    return Integer(0)

def c1(n):
    """
    Return value of $c_1(n)$, which is the number of cusps on the
    modular curve $X_1(n)$.

    EXAMPLES:
        sage: [sage.modular.dims.c1(n) for n in prime_range(2,100)]
        [2, 2, 4, 6, 10, 12, 16, 18, 22, 28, 30, 36, 40, 42, 46, 52, 58, 60, 66, 70, 72, 78, 82, 88, 96]
    """
    if n<3:
        return c0(n)
    if n==4:
        return Integer(3)
    return Integer(sum([frac(phi(d)*phi(n/d),2) \
                    for d in divisors(n)]))

def g1(n):
    """
    Return the genus of the modular curve $X_1(n)$.

    EXAMPLES:
        sage: [sage.modular.dims.g1(n) for n in prime_range(2,100)]
        [0, 0, 0, 0, 1, 2, 5, 7, 12, 22, 26, 40, 51, 57, 70, 92, 117, 126, 155, 176, 187, 222, 247, 287, 345]
    """
    return Integer(1+frac(mu1(n),12)-frac(mu21(n),4)-frac(mu31(n),3) - frac(c1(n),2))


def eisen(p):
    """
    Return the Eisenstein number $n$ which is the numerator of $(p-1)/12$.

    INPUT:
        p -- a prime

    OUTPUT:
        Integer

    EXAMPLES:
        sage: [(p,sage.modular.dims.eisen(p)) for p in prime_range(24)]
        [(2, 1), (3, 1), (5, 1), (7, 1), (11, 5), (13, 1), (17, 4), (19, 3), (23, 11)]
    """
    if not is_prime(p):
        raise ValueError, "p must be prime"
    return frac(p-1,12).numerator()

def S0(n,k):
    r"""
    Return the dimension of the space of level $n$ weight $k$
    cuspforms on $\Gamma_0(n)$.

    INPUT:
        n -- an integer
        k -- an integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.S0(11,2)
        1
        sage: a = sage.modular.dims.S0(1,2); a
        0
        sage: type(a)
        <type 'sage.rings.integer.Integer'>
        sage: sage.modular.dims.S0(5,0)
        0
        sage: sage.modular.dims.S0(20,1)
        0
        sage: sage.modular.dims.S0(20,4)
        6
    """
    n = Integer(n); k = Integer(k)
    if n <= 0:
        raise ValueError, "n must be positive"
    if k<=0 or k%2!=0:
        return Integer(0)
    if k==2:
        return g0(n)
    return Integer((k-1)*(g0(n)-1) + \
               (frac(k,2)-1)*c0(n)+mu20(n)*(k//4)+mu30(n)*(k//3))  # // is floor div

def S1(n,k):
    """
    Return the dimension of the space of level $n$ weight $k$
    cuspforms on $\Gamma_1(n)$.

    INPUT:
        n -- an integer
        k -- an integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.S1(11,2)
        1
        sage: a = sage.modular.dims.S1(13,2); a
        2
        sage: type(a)
        <type 'sage.rings.integer.Integer'>
        sage: sage.modular.dims.S1(20,4)
        26
    """
    n = Integer(n); k = Integer(k)
    if n <= 0:
        raise ValueError, "n must be positive"
    if k<=0 or (n<=2 and k%2!=0):
        return Integer(0)
    if k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."
    if k==2:
        return g1(n)
    if n<=2:
        return S0(n,k)
    a = (k-1)*(g1(n)-1)+(frac(k,2)-1)*c1(n)
    if n == 4 and k%2!=0:
        a += frac(1,2)
    elif n == 3:
        a += k//3  # // is floor div
    return Integer(a)

def idxG0(N):
    r"""
    Return the index $[\SL_2(\Z):\Gamma_0(N)]$ of $\Gamma_0(N)$ in $\SL_2(\Z)$.

    INPUT:
        N -- a positive integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.idxG0(N) for N in [1..10]]
        [1, 3, 4, 6, 6, 12, 8, 12, 12, 18]
        sage: type(sage.modular.dims.idxG0(15))
        <type 'sage.rings.integer.Integer'>
    """
    N = Integer(N)
    if N <= 0:
        raise ValueError, "N must be positive"
    return mul([(p+1)*p**(e-1) for p, e in N.factor()])

def idxG1(N):
    r"""
    The index $[\SL_2(\Z):\Gamma_1(N)]$ of $\Gamma_1(N)$ in $\SL_2(\Z)$.

    INPUT:
        N -- a positive integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.idxG1(N) for N in [1..10]]
        [1, 3, 8, 12, 24, 24, 48, 48, 72, 72]
    """
    return phi(N)*idxG0(N)

#    Formula of Cohen-Oesterle for dim S_k(Gamma_1(N),eps).  REF:
#    Springer Lecture notes in math, volume 627, pages 69--78.  The
#    functions CO_delta and CO_nu, which were first written by Kevin
#    Buzzard, are used only by the function CohenOesterle.

def CO_delta(r,p,N,eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterle.

    INPUT:
        r -- positive integer
        p -- a prime
        N -- positive integer
        $\epsilon$ -- character

    OUTPUT:
        element of the base ring of the character

    EXAMPLES:
        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_delta(1,5,7,eps^3)
        2
    """
    if not is_prime(p):
        raise ValueError, "p must be prime"
    K = eps.base_ring()
    if p%4 == 3:
        return K(0)
    if p==2:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 4).
    # omega is a primitive 4th root of unity mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//4)
    # this n is within a p-power root of a "local" 4th root of 1 modulo p.
    n = Mod(int(omega.crt(Mod(1,N//(p**r)))),N)
    n = n**(p**(r-1))   # this is correct now
    t = eps(n)
    if t==K(1):
        return K(2)
    if t==K(-1):
        return K(-2)
    return K(0)

def CO_nu(r, p, N, eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterle.

    INPUT:
        r -- positive integer
        p -- a prime
        N -- positive integer
        eps -- character

    OUTPUT:
        element of the base ring of the character

    EXAMPLES:
        sage: G.<eps> = DirichletGroup(7)
        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_nu(1,7,7,eps)
        -1
    """
    K = eps.base_ring()
    if p%3==2:
        return K(0)
    if p==3:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 3)
    # omega is a cube root of 1 mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//3)
    n = Mod(omega.crt(Mod(1,N//(p**r))), N)  # within a p-power root of a "local" cube root of 1 mod p.
    n = n**(p**(r-1))  # this is right now
    t = eps(n)
    if t==K(1):
        return K(2)
    return K(-1)

def CohenOesterle(eps, k):
    r"""
    Compute the Cohen-Oesterle function associate to eps, $k$.  This is
    a summand in the formula for the dimension of the space of cusp
    forms of weight $2$ with character $\epsilon$.

    INPUT:
        eps -- Dirichlet character
        k -- integer

    OUTPUT:
        element of the base ring of eps.

    EXAMPLES:
        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CohenOesterle(eps, 2)
        -2/3
        sage: sage.modular.dims.CohenOesterle(eps, 4)
        -1
    """
    N    = eps.modulus()
    facN = factor(N)
    f    = eps.conductor()
    gamma_k = 0
    if k%4==2:
        gamma_k = frac(-1,4)
    elif k%4==0:
        gamma_k = frac(1,4)
    mu_k = 0
    if k%3==2:
        mu_k = frac(-1,3)
    elif k%3==0:
        mu_k = frac(1,3)
    def _lambda(r,s,p):
        """
        Used internally by the CohenOesterle function.

        INPUT:
            r, s, p -- integers
        OUTPUT:
            Integer

        EXAMPLES:   (indirect doctest)
            sage: K = CyclotomicField(3)
            sage: eps = DirichletGroup(7*43,K).0^2
            sage: sage.modular.dims.CohenOesterle(eps,2)
            -4/3
        """
        if 2*s<=r:
            if r%2==0:
                return p**(r//2) + p**((r//2)-1)
            return 2*p**((r-1)//2)
        return 2*(p**(r-s))
    #end def of lambda
    K = eps.base_ring()
    return K(frac(-1,2) * mul([_lambda(r,valuation(f,p),p) for p, r in facN]) + \
               gamma_k * mul([CO_delta(r,p,N,eps)         for p, r in facN]) + \
                mu_k    * mul([CO_nu(r,p,N,eps)            for p, r in facN]))

def dimension_cusp_forms_eps(eps, k=2):
    r"""
    The dimension of the space of cusp forms of weight $k$ and character
    $\epsilon$.

    INPUT:
        $\epsilon$ -- a Dirichlet character
        k -- integer, a weight $\geq 2$.

    OUTPUT:
        integer -- the dimension

    EXAMPLES:
        sage: G.<eps> = DirichletGroup(9)
        sage: [sage.modular.dims.dimension_cusp_forms_eps(eps, k) for k in [2..10]]
        [0, 1, 0, 3, 0, 5, 0, 7, 0]
        sage: [sage.modular.dims.dimension_cusp_forms_eps(eps^2, k) for k in [2..10]]
        [0, 0, 2, 0, 4, 0, 6, 0, 8]
    """
    if isinstance(eps, (int,long) ):
        return dimension_cusp_forms_gamma0(eps,k)

    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        return Z(0)
    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."

    N = eps.modulus()
    if eps.is_trivial():
        return Z(S0(N,k))
    if (eps.is_odd() and k%2==0) or (eps.is_even() and k%2==1):
        return Z(0)
    K = eps.base_ring()
    return Z ( K(idxG0(N)*frac(k-1,12)) + CohenOesterle(eps,k) )


def dimension_eis_eps(eps, k=2):
    r"""
    The dimension of the space of eisenstein series of weight $k$ and
    character $\varepsilon$.

    INPUT:
        eps -- a Dirichlet character
        k -- integer, a weight >= 2.

    OUTPUT:
        integer -- the dimension

    EXAMPLES:
        sage: G.<eps> = DirichletGroup(9)
        sage: [sage.modular.dims.dimension_eis_eps(eps, k) for k in [2..10]]
        [0, 2, 0, 2, 0, 2, 0, 2, 0]
        sage: [sage.modular.dims.dimension_eis_eps(eps^2, k) for k in [2..10]]
        [2, 0, 2, 0, 2, 0, 2, 0, 2]
    """
    if isinstance(eps, (int,long) ):
        return dimension_eis_gamma0(eps,k)

    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        if eps.is_trivial():
            return Z(1)    # the constants
        else:
            return Z(0)

    N = eps.modulus()
    if (eps.is_odd() and k%2==0) or (eps.is_even() and k%2==1):
        return Z(0)
    if eps.is_trivial():
        return dimension_eis(congroup.Gamma0(N),k)
    K = eps.base_ring()
    j = 2-k
    # We use the Cohen-Oesterle formula in a subtle way to
    # compute dim M_k(N,eps) (see Ch. 6 of my book on
    # computing with modular forms).
    alpha = -Z ( K(idxG0(N)*frac(j-1,12)) + CohenOesterle(eps,j) )
    if k == 1:
        return alpha
    else:
        return alpha - dimension_cusp_forms_eps(eps, k)


def dimension_cusp_forms_gamma0(N,k=2):
    r"""
    The dimension of the space $S_k(\Gamma_0(N))$ of cusp forms.

    INPUT:
        N -- integer (the level)
        k -- integer (the weight)

    OUTPUT:
        Integer -- the dimension

    EXAMPLES:
        sage: sage.modular.dims.dimension_cusp_forms_gamma0(23,2)
        2
        sage: sage.modular.dims.dimension_cusp_forms_gamma0(1,24)
        2
        sage: sage.modular.dims.dimension_cusp_forms_gamma0(11,3)
        0
        sage: sage.modular.dims.dimension_cusp_forms_gamma0(11,-1)
        0
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k%2 == 1:
        return Z(0)
    return Z(S0(N,k))

def dimension_cusp_forms_gamma1(N,k=2):
    r"""
    The dimension of the space $S_k(\Gamma_1(N))$ of cusp forms.

    INPUT:
        N -- integer
        k -- integer

    OUTPUT:
        integer -- the dimension

    EXAMPLES:
        sage: sage.modular.dims.dimension_cusp_forms_gamma1(23,2)
        12
        sage: sage.modular.dims.dimension_cusp_forms_gamma1(1,24)
        2
        sage: sage.modular.dims.dimension_cusp_forms_gamma1(11,-1)
        0
        sage: sage.modular.dims.dimension_cusp_forms_gamma1(11,0)
        0
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    if k == 1:
        raise NotImplementedError, "computation of dimensions of spaces of weight 1 modular forms not implemented in general."
    return Z(S1(N,k))


def mumu(N):
    """
    Return 0 if any cube divides $N$.  Otherwise return $(-2)^v$ where
    $v$ is the number of primes that exactly divide $N$.

    This is similar to the Moebius function.

    INPUT:
        N -- an integer at least 1
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.mumu(27)
        0
        sage: sage.modular.dims.mumu(6*25)
        4
        sage: sage.modular.dims.mumu(7*9*25)
        -2
        sage: sage.modular.dims.mumu(9*25)
        1
    """
    if N < 1:
        raise ValueError, "N must be at least 1"
    p = 1
    for _,r in factor(N):
        if r > 2:
            return Z(0)
        elif r == 1:
            p *= -2
    return Z(p)

def dimension_new_cusp_forms_gamma0(N, k=2, p=0):
    r"""
    Dimension of the $p$-new subspace of $S_k(\Gamma_0(N))$.
    If $p=0$, dimension of the new subspace.

    INPUT:
        N -- a positive integer
        k -- an integer
        p -- a prime number

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_new_cusp_forms_gamma0(100, 2, 5)
        5

    Independently compute the dimension 5 above:
        sage: m = ModularSymbols(100, 2,sign=1).cuspidal_subspace()
        sage: m.new_subspace(5)
        Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 18 for Gamma_0(100) of weight 2 with sign 1 over Rational Field
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k%2 == 1:
        return Z(0)
    if p==0 or N%p!=0:
        return sum([dimension_cusp_forms_gamma0(M,k)*mumu(N//M) \
                    for M in divisors(N)])
    return dimension_cusp_forms_gamma0(N,k) - \
           2*dimension_new_cusp_forms_gamma0(N//p,k)

def dimension_new_cusp_forms_gamma1(N,k=2,p=0):
    r"""
    Return the dimension of the $p$-new subspace of
    $S_k(\Gamma_1(N))$.  If $p=0$, return the dimension of the new
    subspace.

    INPUT:
        N -- a positive integer
        k -- an integer
        p -- a prime number

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_new_cusp_forms_gamma1(4*25, 2, 5)
        225
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."

    if p==0 or N%p!=0:
        return Integer(sum([dimension_cusp_forms_gamma1(M,k)*mumu(N//M) \
                    for M in divisors(N)]))
    return dimension_cusp_forms_gamma1(N,k) - \
           2*dimension_new_cusp_forms_gamma1(N//p,k)

def dimension_new_cusp_forms_group(group, k=2, p=0):
    """
    Return the dimension of the new space of cusp forms for the
    congruence subgroup group.  If $p$ is given, return the
    $p$-new subspace.

    INPUT:
        group -- a congruence subgroup
        k -- an integer (default: 2)
        p -- a prime (default: 0); just the $p$-new subspace if given

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_new_cusp_forms_group(Gamma0(33),2)
        1
        sage: sage.modular.dims.dimension_new_cusp_forms_group(Gamma0(33),2,3)
        1
        sage: sage.modular.dims.dimension_new_cusp_forms_group(Gamma0(33),2,11)
        3
        sage: sage.modular.dims.dimension_new_cusp_forms_group(Gamma1(33),2)
        19
        sage: sage.modular.dims.dimension_new_cusp_forms_group(Gamma1(33),2,11)
        21
        sage: sage.modular.dims.dimension_new_cusp_forms_group(GammaH(33,[1,2]),2)
        3
    """
    assert isinstance(group, congroup.CongruenceSubgroup), \
           "Argument 1 must be a congruence subgroup."
    if isinstance(group, congroup.Gamma0):
        return dimension_new_cusp_forms_gamma0(group.level(), k, p)
    elif isinstance(group, congroup.Gamma1):
        return dimension_new_cusp_forms_gamma1(group.level(), k, p)
    elif congroup.is_GammaH(group):
        return dimension_new_cusp_forms_H(group, k, p)
    else:
        raise NotImplementedError, "Computing of dimensions for congruence subgroups besides \
        Gamma0 and Gamma1 is not yet implemented."


def dimension_new_cusp_forms_eps(eps, k=2, p=0):
    r"""
    Dimension of the new subspace (or $p$-new subspace) of cusp forms of
    weight $k$ and character $\epsilon$.

    INPUT:
        eps -- a Dirichlet character
        k -- an integer (ddefault: 2)
        p -- a prime (default: 0); just the $p$-new subspace if given

    OUTPUT:
        Integer

    EXAMPLES:
        sage: G = DirichletGroup(9)
        sage: eps = G.0^3
        sage: eps.conductor()
        3
        sage: [sage.modular.dims.dimension_new_cusp_forms_eps(eps, k) for k in [2..10]]
        [0, 0, 0, 2, 0, 2, 0, 2, 0]
        sage: [sage.modular.dims.dimension_cusp_forms_eps(eps, k) for k in [2..10]]
        [0, 0, 0, 2, 0, 4, 0, 6, 0]
        sage: [sage.modular.dims.dimension_new_cusp_forms_eps(eps, k, 3) for k in [2..10]]
        [0, 0, 0, 2, 0, 2, 0, 2, 0]

    Double check using modular symbols (independent calculation):
        sage: [ModularSymbols(eps,k,sign=1).cuspidal_subspace().new_subspace().dimension()  for k in [2..10]]
        [0, 0, 0, 2, 0, 2, 0, 2, 0]
        sage: [ModularSymbols(eps,k,sign=1).cuspidal_subspace().new_subspace(3).dimension()  for k in [2..10]]
        [0, 0, 0, 2, 0, 2, 0, 2, 0]

    Another example at level 33:
        sage: G = DirichletGroup(33)
        sage: eps = G.1
        sage: eps.conductor()
        11
        sage: [sage.modular.dims.dimension_new_cusp_forms_eps(G.1, k) for k in [2..4]]
        [0, 4, 0]
        sage: [sage.modular.dims.dimension_new_cusp_forms_eps(G.1^2, k) for k in [2..4]]
        [2, 0, 6]
        sage: [sage.modular.dims.dimension_new_cusp_forms_eps(G.1^2, k, 3) for k in [2..4]]
        [2, 0, 6]
    """
    if not isinstance(eps, dirichlet.DirichletCharacter):
        raise TypeError, "eps = (%s) must be a DirichletCharacter"%eps
    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        return Z(0)

    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."


    N = eps.modulus()
    if p == 0 or N%p != 0 or valuation(eps.conductor(),p) == valuation(N,p):
        D = [eps.conductor()*d for d in divisors(N//eps.conductor())]
        return sum([dimension_cusp_forms_eps(eps.restrict(M), k)*mumu(N//M) for M in D])
    eps_p = eps.restrict(N//p)
    old = dimension_cusp_forms(eps_p, k)
    return dimension_cusp_forms(eps, k) - 2*old




######################################################################
# Computing dimensions of modular forms spaces for Gamma_H.
# Algorithms found and implemented by Jordi Quer.
######################################################################

def muH(N,H):
    r"""
    Return the degree of the covering map $X_H(N) \rightarrow X_0(1)$.

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.muH(33, GammaH(33,[2])._list_of_elements_in_H())
        96

    AUTHOR: Jordi Quer
    """
    lenHpm = len(H)
    if N-Integer(1) not in H: lenHpm*=Integer(2)
    return mul([(p**Integer(2)-Integer(1))*(p**(Integer(2)*r-Integer(2))) for p, r in factor(N)])//lenHpm

def nu2H(N,H):
    r"""
    Number of elliptic points of order 2 for the group $\Gamma_H(N)$

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.nu2H(33, GammaH(33,[2])._list_of_elements_in_H())
        0
        sage: sage.modular.dims.nu2H(5, GammaH(5,[2])._list_of_elements_in_H())
        2

    AUTHOR: Jordi Quer
    """
    if N%Integer(4) == Integer(0): return Integer(0)
    for p, r in factor(N):
        if p%Integer(4) ==Integer(3): return Integer(0)
    return (euler_phi(N)//len(H))*len([x for x in H if (x**Integer(2)+Integer(1))%N == Integer(0)])

def nu3H(N,H):
    r"""
    Number of elliptic points of order 3 for the group $\Gamma_H(N)$

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.nu3H(33, GammaH(33,[2])._list_of_elements_in_H())
        0
        sage: sage.modular.dims.nu3H(7, GammaH(7,[2])._list_of_elements_in_H())
        2

    AUTHOR: Jordi Quer
    """
    if N%Integer(9) == Integer(0): return Integer(0)
    for p, r in factor(N):
        if p%Integer(3) == Integer(2): return Integer(0)
    lenHpm = len(H)
    if N-Integer(1) not in H: lenHpm*=Integer(2)
    return (euler_phi(N)//lenHpm)*len([x for x in H if (x**Integer(2)+x+Integer(1))%N == Integer(0)])

def nuinfH(N,H):
    r"""
    Number of cusps for the group $\Gamma_H(N)$

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.nuinfH(33, GammaH(33,[2])._list_of_elements_in_H())
        8

    AUTHOR: Jordi Quer
    """
    c = Integer(0)
    for d in [d for d in divisors(N) if d**Integer(2)<=N]:
        Nd = lcm(d,N//d)
        Hd = set([x%Nd for x in H])
        lenHd = len(Hd)
        if Nd-Integer(1) not in Hd: lenHd*=Integer(2)
        sumand = euler_phi(d)*euler_phi(N//d)//lenHd
        if d**Integer(2)==N:
            c = c + sumand
        else:
            c = c + Integer(2)*sumand
    return c

def nuinfHreg(N,H):
    r"""
    Number of regular cusps for the group $\Gamma_H(N)$

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.nuinfHreg(33, GammaH(33,[2])._list_of_elements_in_H())
        0
        sage: sage.modular.dims.nuinfHreg(7, GammaH(7,[2])._list_of_elements_in_H())
        2

    AUTHOR: Jordi Quer
    """
    c = Integer(0)
    for d in [d for d in divisors(N) if d**Integer(2)<=N]:
        Nd = lcm(d,N//d)
        Hd = set([x%Nd for x in H])
        if Nd-Integer(1) not in Hd:
            sumand = euler_phi(d)*euler_phi(N//d)//(Integer(2)*len(Hd))
            if d**Integer(2)==N:
                c = c + sumand
            else:
                c = c + Integer(2)*sumand
    return c


def gH(N,H):
    r"""
    Return the genus of the curve $X_H(N)$.

    INPUT:
        N -- an integer
        H -- list of elements in $\ZZ$ of a subgroup of $(\ZZ/N\ZZ)^*$
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.gH(33, GammaH(33,[2])._list_of_elements_in_H())
        5
        sage: sage.modular.dims.gH(7, GammaH(7,[2])._list_of_elements_in_H())
        0
        sage: sage.modular.dims.gH(23, [1..22])
        2
        sage: sage.modular.dims.g0(23)
        2
        sage: sage.modular.dims.gH(23, [1])
        12
        sage: sage.modular.dims.g1(23)
        12

    AUTHOR: Jordi Quer
    """
    # (formulas are multiplied by 12 so that everything is an integer)
    return (Integer(12) + muH(N,H) - Integer(3)*nu2H(N,H) - Integer(4)*nu3H(N,H) - Integer(6)*nuinfH(N,H))//Integer(12)

def genus_H(H):
    r"""
    Return the genus of the modular curve $X_H(N)$, where $H$ is a
    congruence subgroup of the form $\Gamma_H(N)$.

    INPUT:
        H -- a congruence subgroup $\Gamma_H(N)$
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.genus_H(GammaH(33,[2]))
        5

    AUTHOR: Jordi Quer
    """
    if not congroup.is_GammaH(H):
        raise TypeError, "H must be a congruence subgroup GammaH(N)"
    return gH(H.level(),H._list_of_elements_in_H())

def lambda4(k):
    """
    Return 0 if $k$ is odd, 3 if $k$ is divisible by 4, and -3 if $k$
    is even but not divisible by 4.

    This is used for computations of the Cohen-Oesterle formulas.

    INPUT:
        k -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.lambda4(k) for k in [0..10]]
        [3, 0, -3, 0, 3, 0, -3, 0, 3, 0, -3]
    """
    if k%Integer(2) == Integer(1):
        return Integer(0)
    elif k%Integer(4) == Integer(2):
        return -Integer(3)
    else:
        return Integer(3)

def lambda3(k):
    """
    Return 0 if $k$ is odd, 4 if $k$ is divisible by 4, and -4 if $k$
    is even but not divisible by 4.

    This is used for computations of the Cohen-Oesterle formulas.

    INPUT:
        k -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.lambda3(k) for k in [0..10]]
        [4, 0, -4, 4, 0, -4, 4, 0, -4, 4, 0]
    """
    if k%Integer(3) == Integer(1):
        return Integer(0)
    elif k%Integer(3) == Integer(2):
        return -Integer(4)
    else:
        return Integer(4)

###############################################################
# dimensions of spaces of cusp and modular forms for Gamma_H
###############################################################
def dimension_cusp_forms_H(G,k):
    r"""
    Return the dimension of the space of weight $k$ cusp forms for the
    group $G = \Gamma_H$.

    INPUT:
        G -- a congruence subgroup GammaH
        k -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_cusp_forms_H(GammaH(33,[2]),2)
        5
        sage: sage.modular.dims.dimension_cusp_forms_H(GammaH(33,[2]),3)
        0
        sage: sage.modular.dims.dimension_cusp_forms_H(GammaH(33,[2,5]),2)
        3
    """
    if not congroup.is_GammaH(G):
        raise TypeError, "H must be a congruence subgroup GammaH(N)"

    N = G.level()
    H = G._list_of_elements_in_H()
    if k%Integer(2) == Integer(1) and N-Integer(1) in H:
        return Integer(0)
    dim = Integer(0)
    if k == Integer(2):
        dim = Integer(12)
    dim += (k-Integer(1))*muH(N,H)
    if k%Integer(2) == Integer(0):
        dim += lambda4(k)*nu2H(N,H)
    dim += lambda3(k)*nu3H(N,H)
    if k%Integer(2) == Integer(0):
        dim += -Integer(6)*nuinfH(N,H)
    else:
        dim += -Integer(6)*nuinfHreg(N,H)
    return dim//Integer(12)

def dimension_eis_H(G,k):
    r"""
    Return the dimension of the space of weight $k$ Eisenstein series
    for the group $G = \Gamma_H$.

    INPUT:
        G -- a congruence subgroup GammaH
        k -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_eis_H(GammaH(33,[2]),2)
        7
        sage: sage.modular.dims.dimension_eis_H(GammaH(33,[2]),3)
        0
        sage: sage.modular.dims.dimension_eis_H(GammaH(33,[2,5]),2)
        3
    """
    N = G.level()
    H = G._list_of_elements_in_H()
    if k%Integer(2) == Integer(1) and N-Integer(1) in H:
        return Integer(0)
    if k%Integer(2) == Integer(0):
        dim = nuinfH(N,H)
    else:
        dim = nuinfHreg(N,H)
    if k == Integer(2):
        dim -= Integer(1)
    return dim

def dimension_new_cusp_forms_H(G,k, p=0):
    r"""
    Return the dimension of the space of new (or $p$-new) weight $k$
    cusp forms for $G$.

    INPUT:
        G -- group of the form Gamma_H(N)
        k -- an integer at least 2 (the weight)
        p -- integer (default: 0); if nonzero, compute the $p$-new subspace.

    OUTPUT:
        Integer

    EXAMPLES:
        sage: from sage.modular.dims import *
        sage: dimension_new_cusp_forms_H(GammaH(33,[2]), 2)
        3
    """
    N = G.level()
    if p==0 or N%p!=0:
        return sum([dimension_cusp_forms_H(H,k) * mumu(N // H.level()) \
                for H in G.divisor_subgroups()])
    else:
        return dimension_cusp_forms_H(G,k) - \
               2*dimension_new_cusp_forms_H(G.restrict(N//p),k)


def dimension_cusp_forms_fromH(chi,k):
    r"""
    Compute the dimension of the space of cusp forms of weight $k$ and
    character $\chi$ using formulas for dimensions for congruence
    subgroups $\Gamma_H(N)$ instead of the formulas of Cohen-Oesterle.

    INPUT:
        chi -- Dirichlet character
        k -- integer (weight)

    OUTPUT:
        Integer

    EXAMPLES:
        sage: K = CyclotomicField(3)
        sage: eps = DirichletGroup(7*43,K).0^2

    The following two calculations use different algorithms.
        sage: dimension_cusp_forms(eps,2)
        28
        sage: sage.modular.dims.dimension_cusp_forms_fromH(eps,2)
        28
    """
    N = chi.modulus()
    n = chi.order()
    dim = Integer(0)
    for d in divisors(n):
        G = congroup.GammaH(N,(chi**d).kernel())
        dim = dim + moebius(d)*dimension_cusp_forms_H(G,k)
    return dim//euler_phi(n)

def dimension_modular_forms_fromH(chi,k):
    r"""
    Compute the dimension of the space of modular forms of weight $k$
    and character $\chi$ using formulas for dimensions for congruence
    subgroups $\Gamma_H(N)$ instead of the formulas of Cohen-Oesterle.

    INPUT:
        chi -- Dirichlet character
        k -- integer (weight)

    OUTPUT:
        Integer

    EXAMPLES:
        sage: K = CyclotomicField(3)
        sage: eps = DirichletGroup(7*43,K).0^2
        sage: dimension_modular_forms(eps,2)
        32
        sage: sage.modular.dims.dimension_modular_forms_fromH(eps,2)
        32
    """
    N = chi.modulus()
    n = chi.order()
    dim = Integer(0)
    for d in divisors(n):
        G = congroup.GammaH(N,(chi**d).kernel())
        dim = dim + moebius(d)*(dimension_eis_H(G,k) + dimension_cusp_forms_H(G,k))
    return dim//euler_phi(n)

####################################################################
# Functions exported to the global namespace.
# These have very flexible inputs.
####################################################################

def dimension_new_cusp_forms(X, k=2, p=0):
    """
    Return the dimension of the new (or $p$-new) subspace of cusp
    forms for the character or group $X$.

    INPUT:
        X -- integer, congruence subgroup or Dirichlet character
        k -- weight (integer)
        p -- 0 or a prime

    EXAMPLES:
        sage: dimension_new_cusp_forms(100,2)
        1

        sage: dimension_new_cusp_forms(Gamma0(100),2)
        1
        sage: dimension_new_cusp_forms(Gamma0(100),4)
        5

        sage: dimension_new_cusp_forms(Gamma1(100),2)
        141
        sage: dimension_new_cusp_forms(Gamma1(100),4)
        463

        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,2)
        2
        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,4)
        8

        sage: sum(dimension_new_cusp_forms(e,3) for e in DirichletGroup(30))
        12
        sage: dimension_new_cusp_forms(Gamma1(30),3)
        12
    """
    if isinstance(X, congroup.CongruenceSubgroup):
        return dimension_new_cusp_forms_group(X,k,p)
    elif isinstance(X, dirichlet.DirichletCharacter):
        return dimension_new_cusp_forms_eps(X,k,p)
    elif isinstance(X, (int,long,Integer)):
        return dimension_new_cusp_forms_gamma0(X,k,p)
    else:
        raise TypeError, "X (=%s) must be an integer, congruence subgroup or Diirichlet character"%X

def dimension_cusp_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup or Dirichlet character.

    INPUT:
        X -- congruence subgroup or Dirichlet character or integer
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_cusp_forms(5,4)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(13),2)
        2

        sage: dimension_cusp_forms(DirichletGroup(13).0^2,2)
        1
        sage: dimension_cusp_forms(DirichletGroup(13).0,3)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma0(11),0)
        0
        sage: dimension_cusp_forms(Gamma0(1),12)
        1
        sage: dimension_cusp_forms(Gamma0(1),2)
        0
        sage: dimension_cusp_forms(Gamma0(1),4)
        0

        sage: dimension_cusp_forms(Gamma0(389),2)
        32
        sage: dimension_cusp_forms(Gamma0(389),4)
        97
        sage: dimension_cusp_forms(Gamma0(2005),2)
        199
        sage: dimension_cusp_forms(Gamma0(11),1)
        0

        sage: dimension_cusp_forms(Gamma1(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(1),12)
        1
        sage: dimension_cusp_forms(Gamma1(1),2)
        0
        sage: dimension_cusp_forms(Gamma1(1),4)
        0

        sage: dimension_cusp_forms(Gamma1(389),2)
        6112
        sage: dimension_cusp_forms(Gamma1(389),4)
        18721
        sage: dimension_cusp_forms(Gamma1(2005),2)
        159201

        sage: dimension_cusp_forms(Gamma1(11),1)
        Traceback (most recent call last):
        ...
        NotImplementedError: computation of dimensions of spaces of weight 1 modular forms not implemented in general.

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_cusp_forms(e,2)
        0
        sage: dimension_cusp_forms(e^2,2)
        1
    """
    if isinstance(X, dirichlet.DirichletCharacter):
        return dimension_cusp_forms_eps(X, k)
    elif isinstance(X, congroup.Gamma0):
        return dimension_cusp_forms_gamma0(X.level(),k)
    elif isinstance(X, (Integer,int,long)):
        return dimension_cusp_forms_gamma0(Integer(X),k)
    elif isinstance(X, congroup.Gamma1):
        return dimension_cusp_forms_gamma1(X.level(),k)
    elif congroup.is_GammaH(X):
        return dimension_cusp_forms_H(X,k)
    elif not isinstance(X, congroup.CongruenceSubgroup):
        raise TypeError, "Argument 1 must be a congruence subgroup or Dirichlet character"
    else:
        raise NotImplementedError, "Computing of dimensions for congruence subgroups besides \
        Gamma0, Gamma1, and GammaH is not yet implemented."

def dimension_eis(X, k=2):
    """
    The dimension of the space of eisenstein series for the given
    congruence subgroup.

    INPUT:
        X -- congruence subgroup or Dirichlet character or integer
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_eis(5,4)
        2

        sage: dimension_eis(Gamma0(11),2)
        1
        sage: dimension_eis(Gamma1(13),2)
        11
        sage: dimension_eis(Gamma1(2006),2)
        3711

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2
        sage: dimension_eis(e,13)
        2

        sage: G = DirichletGroup(20)
        sage: dimension_eis(G.0,3)
        4
        sage: dimension_eis(G.1,3)
        6
        sage: dimension_eis(G.1^2,2)
        6

        sage: G = DirichletGroup(200)
        sage: e = prod(G.gens(), G(1))
        sage: e.conductor()
        200
        sage: dimension_eis(e,2)
        4

        sage: dimension_modular_forms(Gamma1(4), 11)
        6
    """
    if k <= 1:
        # TODO
        raise NotImplementedError, "Dimension of weight <= 1 Eisenstein series not yet implemented."
    if isinstance(X, (int,long,Integer)):
        if k%2 == 1: return 0
        d = c0(Integer(X))
        if k==2: d -= 1
        return Z(d)
    if isinstance(X, dirichlet.DirichletCharacter):
        return dimension_eis_eps(X, k)
    if isinstance(X, congroup.Gamma0):
        if k%2 == 1: return 0
        d = c0(X.level())
        if k==2: d -= 1
        return Z(d)
    elif isinstance(X, congroup.Gamma1):
        N = X.level()
        if N == 2 and k%2 == 1:
            d = 0  # level Gamma1(2) and odd weight is a special case.
        elif N == 4 and k%2 == 1:
            d = 2  # level Gamma1(4) and odd weight is a special case.
        else:
            d = c1(N)
        if k==2: d -= 1
        return Z(d)
    elif congroup.is_GammaH(X):
        return dimension_eis_H(X, k)
    elif isinstance(X, congroup.CongruenceSubgroup):
        raise NotImplementedError, "Computation of dimensions for congruence subgroups besides " + \
              "Gamma0 and Gamma1 is not yet implemented."
    else:
        raise TypeError


def dimension_modular_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup (either $\Gamma_0(N)$, $\Gamma_1(N)$, or $\Gamma_H(N)$)
    or Dirichlet character.

    INPUT:
        X -- congruence subgroup or Dirichlet character
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_modular_forms(Gamma0(11),2)
        2
        sage: dimension_modular_forms(Gamma0(11),0)
        1
        sage: dimension_modular_forms(Gamma1(13),2)
        13

        sage: e = DirichletGroup(20).1
        sage: dimension_modular_forms(e,3)
        9
        sage: dimension_cusp_forms(e,3)
        3
        sage: dimension_eis(e,3)
        6
        sage: dimension_modular_forms(11,2)
        2
    """
    if isinstance(X, (int, long, Integer)):
        X = congroup.Gamma0(X)
    elif not isinstance(X, congroup.CongruenceSubgroup) and \
         not isinstance(X, dirichlet.DirichletCharacter):
        raise TypeError, "Argument 1 must be a congruence subgroup or Dirichlet character."
    if k == 0:
        return 1
    if congroup.is_GammaH(X):
        return dimension_modular_forms_H(X, k)
    return dimension_cusp_forms(X, k) + dimension_eis(X, k)

def sturm_bound(level, weight=2):
    r"""
    Returns the Sturm bound for modular forms with given level and
    weight.

    INPUT:
        level -- an integer or a congruence subgroup
        weight -- an integer $\geq 2$ (default: 2)

    EXAMPLES:
        sage: sturm_bound(11,2)
        2
        sage: sturm_bound(389,2)
        65
        sage: sturm_bound(1,12)
        1
        sage: sturm_bound(100,2)
        30
        sage: sturm_bound(1,36)
        3
        sage: sturm_bound(11)
        2
        sage: sturm_bound(Gamma0(11))
        2
        sage: sturm_bound(Gamma0(13))
        3
        sage: sturm_bound(Gamma0(16))
        4
        sage: sturm_bound(GammaH(16,[13]))
        8
        sage: sturm_bound(GammaH(16,[15]))
        16
        sage: sturm_bound(Gamma1(16))
        32
        sage: sturm_bound(Gamma1(13))
        36
        sage: sturm_bound(Gamma1(13),5)
        72

    FURTHER DETAILS: This function returns a positive integer~$n$ such
    that the Hecke operators $T_1,\ldots, T_n$ acting on \emph{cusp
    forms} generate the Hecke algebra as a $\Z$-module when the
    character is trivial or quadratic.  Otherwise, $T_1,\ldots,T_n$
    generate the Hecke algebra at least as a $\Z[\eps]$-module, where
    $\Z[\eps]$ is the ring generated by the values of the Dirichlet
    character $\eps$.  Alternatively, this is a bound such that if two
    cusp forms associated to this space of modular symbols are
    congruent modulo $(\lambda, q^n)$, then they are congruent modulo
    $\lambda$.

    REFERENCES:
    See the Agashe-Stein appendix to Lario and Schoof, \emph{Some
    computations with Hecke rings and deformation rings},
    Experimental Math., 11 (2002), no. 2, 303-311.  This result
    originated in the paper Sturm, \emph{On the congruence of
    modular forms}, Springer LNM 1240, 275--280, 1987.

    REMARK:
    Kevin Buzzard pointed out to me (William Stein) in Fall 2002
    that the above bound is fine for $\Gamma_1(N)$ with character,
    as one sees by taking a power of $f$.  More precisely, if $f
    \con 0 \pmod{p}$ for first $s$ coefficients, then $f^r \con 0
    \pmod{p}$ for first $sr$ coefficents.  Since the weight of
    $f^r$ is $r\cdot k(f)$, it follows that if $s \geq b$, where
    $b$ is the Sturm bound for $\Gamma_0(N)$ at weight $k(f)$, then
    $f^r$ has valuation large enough to be forced to be $0$ at
    $r*k(f)$ by Sturm bound (which is valid if we choose $r$
    correctly).  Thus $f \con 0 \pmod{p}$.  Conclusion: For
    $\Gamma_1(N)$ with fixed character, the Sturm bound is
    \emph{exactly} the same as for $\Gamma_0(N)$.

    A key point is that we are finding $\Z[\eps]$ generators for the
    Hecke algebra here, not $\Z$-generators.  So if one wants
    generators for the Hecke algebra over $\Z$, this bound must be
    suitably modified (and I'm not sure what the modification is).

    AUTHOR:
        -- William Stein
    """
    if congroup.is_Gamma0(level):
        level = level.level()
    elif congroup.is_GammaH(level):
        N = level.level()
        index = euler_phi(N) // len(level._list_of_elements_in_H())
        return sturm_bound(N, weight=weight) * index
    elif congroup.is_Gamma1(level):
        N = level.level()
        return sturm_bound(N, weight=weight) * euler_phi(N)
    return Integer(int(math.ceil((weight * idxG0(level) / ZZ(12)))))

